/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck,
 * 2020 - today
 *
 * festigkeitslehre@uibk.ac.at
 *
 * Alexandros Stathas alexandros.stathas@boku.ac.at
 * Matthias Neuner matthias.neuner@uibk.ac.at
 *
 * This file is part of the MAteRialMOdellingToolbox (marmot).
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of marmot.
 * ---------------------------------------------------------------------
 */

#include "Marmot/MarmotEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain.h"

#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrainFactory.h"

#include <stdexcept>

using namespace Marmot::FastorStandardTensors;
using namespace Marmot::FastorIndices;

MarmotEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain::
  MarmotEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain( const std::string& materialName,
                                                                   const double*      matProperties_,
                                                                   int                nMaterialProperties_,
                                                                   int                materialNumber_ )
  : MarmotEquilibratedInterfaceMaterialFiniteStrainBase( matProperties_,
                                                         nMaterialProperties_,
                                                         materialNumber_,
                                                         "MarmotEquilibratedGradientEnhancedInterfaceMaterialFiniteStra"
                                                         "in" )
{
  // Both faces are driven by the same gradient-enhanced law, each with its own instance and its own history.
  const auto createBulkMaterial = [&]() {
    return std::unique_ptr< MarmotMaterialGradientEnhancedFiniteStrain >(
      MarmotLibrary::MarmotMaterialGradientEnhancedFiniteStrainFactory::createMaterial( materialName,
                                                                                        baseMaterialProperties.data(),
                                                                                        static_cast< int >(
                                                                                          baseMaterialProperties
                                                                                            .size() ),
                                                                                        materialNumber ) );
  };

  topMaterial    = createBulkMaterial();
  bottomMaterial = createBulkMaterial();

  finalizeStateLayout( topMaterial->getNumberOfRequiredStateVars(), bottomMaterial->getNumberOfRequiredStateVars() );
}

void MarmotEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain::setCharacteristicElementLength( double length )
{
  topMaterial->setCharacteristicElementLength( length );
  bottomMaterial->setCharacteristicElementLength( length );
}

MarmotEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain::Solution MarmotEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain::
  solveWithFields( double*                                                                 stateVars,
                   const MarmotEquilibratedInterfaceMaterialFiniteStrainBase::Deformation& deformation,
                   const TimeIncrement&                                                    timeIncrement,
                   double                                                                  nonlocalFieldTop,
                   double                                                                  nonlocalFieldBottom ) const
{
  const MarmotMaterialGradientEnhancedFiniteStrain::TimeIncrement materialTime{ timeIncrement.timeOld,
                                                                                timeIncrement.dT };

  // One face at a total displacement gradient and the nonlocal field, from a scratch copy of the committed state, see
  // MarmotEquilibratedInterfaceMaterialFiniteStrain.
  const auto evaluateSide = [&]( Side side, const Matrix3dRowMajor& G ) {
    const MarmotMaterialGradientEnhancedFiniteStrain& material = side == Side::Top ? *topMaterial : *bottomMaterial;

    const double*         committedStateVars = bulkStateVars( side, stateVars );
    const int             nStateVars         = material.getNumberOfRequiredStateVars();
    std::vector< double > trialStateVars( committedStateVars, committedStateVars + nStateVars );

    Tensor33d F( 0.0 );
    F.eye();
    for ( int i = 0; i < 3; ++i )
      for ( int j = 0; j < 3; ++j )
        F( i, j ) += G( i, j );

    MarmotMaterialGradientEnhancedFiniteStrain::ConstitutiveResponse< 3 > response;
    response.stateVars = trialStateVars.data();
    MarmotMaterialGradientEnhancedFiniteStrain::AlgorithmicModuli< 3 > moduli;
    const MarmotMaterialGradientEnhancedFiniteStrain::Deformation< 3 > materialDeformation{ F,
                                                                                            side == Side::Top
                                                                                              ? nonlocalFieldTop
                                                                                              : nonlocalFieldBottom };

    material.computeStress( response, moduli, materialDeformation, materialTime );

    SideTrial trial = firstPiolaFromKirchhoff( F, response.tau, moduli.dTau_dF );

    // F is held fixed while abar varies, so that only the Kirchhoff stress moves
    const Tensor33d FInvT        = Fastor::transpose( Fastor::inverse( F ) );
    const Tensor33d dP_dNonlocal = Fastor::einsum< ij, jk, to_ik >( moduli.dTau_dN, FInvT );
    for ( int i = 0; i < 3; ++i )
      for ( int j = 0; j < 3; ++j ) {
        trial.dP_dNonlocal( i, j ) = dP_dNonlocal( i, j );
        trial.dL_dF( 3 * i + j )   = moduli.dL_dF( i, j );
      }

    trial.dL_dN          = moduli.dL_dN;
    trial.L              = response.L;
    trial.nonlocalRadius = response.nonLocalRadius;
    trial.trialStateVars = std::move( trialStateVars );
    return trial;
  };

  return solve( stateVars, deformation, evaluateSide, true );
}

void MarmotEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain::computeStress(
  State&               state,
  Tangents&            tangents,
  const Deformation&   deformation,
  const TimeIncrement& timeIncrement )
{
  const Solution solution = solveWithFields( state.stateVars,
                                             deformation,
                                             timeIncrement,
                                             deformation.nonlocalField,
                                             deformation.nonlocalField );

  writeResponse( solution, state, tangents );

  state.localField     = solution.localField;
  state.nonlocalRadius = solution.nonlocalRadius;

  Eigen::Map< Vector21d >( tangents.dGeneralizedStress_dNonlocal.data() ) = solution.dStress_dNonlocal;
  Eigen::Map< Eigen::Matrix< double, 1, 21 > >(
    tangents.dLocalField_dGeneralizedStrain.data() ) = solution.dLocalField_dStrain;
  tangents.dLocalField_dNonlocal                     = solution.dLocalField_dNonlocal;
}

void MarmotEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain::computeStress(
  FaceFieldsState&             state,
  FaceFieldsTangents&          tangents,
  const FaceFieldsDeformation& deformation,
  const TimeIncrement&         timeIncrement )
{
  const Solution solution = solveWithFields( state.stateVars,
                                             deformation,
                                             timeIncrement,
                                             deformation.nonlocalFieldTop,
                                             deformation.nonlocalFieldBottom );

  writeResponse( solution, state, tangents );

  state.localFieldTop        = solution.plus.L;
  state.localFieldBottom     = solution.minus.L;
  state.nonlocalRadiusTop    = solution.plus.nonlocalRadius;
  state.nonlocalRadiusBottom = solution.minus.nonlocalRadius;

  Eigen::Map< Vector21d >( tangents.dGeneralizedStress_dNonlocalTop.data() )    = solution.dStress_dNonlocalPlus;
  Eigen::Map< Vector21d >( tangents.dGeneralizedStress_dNonlocalBottom.data() ) = solution.dStress_dNonlocalMinus;
  Eigen::Map< Eigen::Matrix< double, 1, 21 > >(
    tangents.dLocalFieldTop_dGeneralizedStrain.data() ) = solution.dLocalFieldPlus_dStrain;
  Eigen::Map< Eigen::Matrix< double, 1, 21 > >(
    tangents.dLocalFieldBottom_dGeneralizedStrain.data() ) = solution.dLocalFieldMinus_dStrain;
  tangents.dLocalFieldTop_dNonlocalTop                     = solution.dLocalField_dNonlocalFaces( 0, 0 );
  tangents.dLocalFieldTop_dNonlocalBottom                  = solution.dLocalField_dNonlocalFaces( 0, 1 );
  tangents.dLocalFieldBottom_dNonlocalTop                  = solution.dLocalField_dNonlocalFaces( 1, 0 );
  tangents.dLocalFieldBottom_dNonlocalBottom               = solution.dLocalField_dNonlocalFaces( 1, 1 );
}

void MarmotEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain::initializeYourself( double* stateVars,
                                                                                          int     nStateVars )
{
  for ( int i = 0; i < nStateVars; ++i )
    stateVars[i] = 0.0;

  topMaterial->initializeYourself( bulkStateVars( Side::Top, stateVars ), topMaterial->getNumberOfRequiredStateVars() );
  bottomMaterial->initializeYourself( bulkStateVars( Side::Bottom, stateVars ),
                                      bottomMaterial->getNumberOfRequiredStateVars() );
}

double MarmotEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain::getDensity( const double* stateVars ) const
{
  return topMaterial->getDensity( bulkStateVars( Side::Top, stateVars ) );
}
