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

#include "Marmot/MarmotEquilibratedInterfaceMaterialFiniteStrain.h"

#include "Marmot/MarmotMaterialFiniteStrainFromGradientEnhanced.h"

#include <stdexcept>

using namespace Marmot::FastorStandardTensors;

MarmotEquilibratedInterfaceMaterialFiniteStrain::MarmotEquilibratedInterfaceMaterialFiniteStrain(
  const std::string& materialName,
  const double*      matProperties_,
  int                nMaterialProperties_,
  int                materialNumber_ )
  : MarmotEquilibratedInterfaceMaterialFiniteStrainBase( matProperties_,
                                                         nMaterialProperties_,
                                                         materialNumber_,
                                                         "MarmotEquilibratedInterfaceMaterialFiniteStrain" )
{
  // A gradient-enhanced material is accepted too, with its nonlocal field held fixed. The adapter refuses a material
  // whose stress depends on the field, so the substitution is never silent.
  const auto createBulkMaterial = [&]() {
    return MarmotLibrary::createFiniteStrainMaterialAllowingGradientEnhanced( materialName,
                                                                              baseMaterialProperties.data(),
                                                                              static_cast< int >(
                                                                                baseMaterialProperties.size() ),
                                                                              materialNumber );
  };

  topMaterial    = createBulkMaterial();
  bottomMaterial = createBulkMaterial();

  finalizeStateLayout( topMaterial->getNumberOfRequiredStateVars(), bottomMaterial->getNumberOfRequiredStateVars() );
}

void MarmotEquilibratedInterfaceMaterialFiniteStrain::setCharacteristicElementLength( double length )
{
  topMaterial->setCharacteristicElementLength( length );
  bottomMaterial->setCharacteristicElementLength( length );
}

void MarmotEquilibratedInterfaceMaterialFiniteStrain::computeStress( State&               state,
                                                                     Tangents&            tangents,
                                                                     const Deformation&   deformation,
                                                                     const TimeIncrement& timeIncrement )
{
  const MarmotMaterialFiniteStrain::TimeIncrement materialTime{ timeIncrement.timeOld, timeIncrement.dT };

  // One face at a total displacement gradient, from a scratch copy of the committed state. The persistent state
  // variables are never touched: a local Newton iteration re-evaluates the same increment many times, and a material
  // that mutated its state in place would integrate it once per trial.
  const auto evaluateSide = [&]( Side side, const Matrix3dRowMajor& G ) {
    MarmotMaterialFiniteStrain& material = side == Side::Top ? *topMaterial : *bottomMaterial;

    const double*         committedStateVars = bulkStateVars( side, state.stateVars );
    const int             nStateVars         = material.getNumberOfRequiredStateVars();
    std::vector< double > trialStateVars( committedStateVars, committedStateVars + nStateVars );

    Tensor33d F( 0.0 );
    F.eye();
    for ( int i = 0; i < 3; ++i )
      for ( int j = 0; j < 3; ++j )
        F( i, j ) += G( i, j );

    MarmotMaterialFiniteStrain::ConstitutiveResponse< 3 > response( Tensor33d( 0.0 ), 0.0, 0.0, trialStateVars.data() );
    MarmotMaterialFiniteStrain::AlgorithmicModuli< 3 >    moduli;
    const MarmotMaterialFiniteStrain::Deformation< 3 >    materialDeformation{ F };

    material.computeStress( response, moduli, materialDeformation, materialTime );

    SideTrial trial      = firstPiolaFromKirchhoff( F, response.tau, moduli.dTau_dF );
    trial.trialStateVars = std::move( trialStateVars );
    return trial;
  };

  const Solution solution = solve( state.stateVars, deformation, evaluateSide, false );

  writeResponse( solution, state, tangents );
}

void MarmotEquilibratedInterfaceMaterialFiniteStrain::initializeYourself( double* stateVars, int nStateVars )
{
  for ( int i = 0; i < nStateVars; ++i )
    stateVars[i] = 0.0;

  topMaterial->initializeYourself( bulkStateVars( Side::Top, stateVars ), topMaterial->getNumberOfRequiredStateVars() );
  bottomMaterial->initializeYourself( bulkStateVars( Side::Bottom, stateVars ),
                                      bottomMaterial->getNumberOfRequiredStateVars() );
}

double MarmotEquilibratedInterfaceMaterialFiniteStrain::getDensity() const
{
  return topMaterial->getDensity( nullptr );
}
