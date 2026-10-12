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

#include "Marmot/MarmotInterfaceMaterialFiniteStrain.h"

#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotInterfaceFiniteStrainHelperFunctions.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"
#include "Marmot/MarmotUtils.h"

#include "Fastor/Fastor.h"

#include <stdexcept>
#include <string>

using namespace Marmot;
using namespace Marmot::FastorIndices;
using namespace Marmot::FastorStandardTensors;

MarmotInterfaceMaterialFiniteStrain::MarmotInterfaceMaterialFiniteStrain( const std::string& materialName,
                                                                          const double*      matProperties_,
                                                                          int                nMaterialProperties_,
                                                                          int                materialNumber_ )
  : materialProperties( matProperties_ ),
    nMaterialProperties( nMaterialProperties_ ),
    h( checkedMaterialProperty( matProperties_, nMaterialProperties_, 0 ) ),
    materialNumber( materialNumber_ )
{
  if ( h <= 0.0 )
    throw std::invalid_argument( "MarmotInterfaceMaterialFiniteStrain requires a positive interface thickness h." );

  baseMaterialProperties.assign( materialProperties + 1, materialProperties + nMaterialProperties );

  baseMaterial = std::unique_ptr< MarmotMaterialFiniteStrain >(
    MarmotLibrary::MarmotMaterialFiniteStrainFactory::createMaterial( materialName,
                                                                      baseMaterialProperties.data(),
                                                                      static_cast< int >(
                                                                        baseMaterialProperties.size() ),
                                                                      materialNumber ) );

  stateLayout.add( "baseMaterialStateVars", baseMaterial->getNumberOfRequiredStateVars() );
  stateLayout.finalize();
}

void MarmotInterfaceMaterialFiniteStrain::setCharacteristicElementLength( double length )
{
  baseMaterial->setCharacteristicElementLength( length );
}

void MarmotInterfaceMaterialFiniteStrain::computeStress( State&               state,
                                                         Tangents&            tangents,
                                                         const Deformation&   deformation,
                                                         const TimeIncrement& timeIncrement )
{
  using namespace Marmot::Materials::InterfaceFiniteStrainHelperFunctions;

  const Tensor3d& normal = deformation.normal;

  const NormalFibre fibre = normalFibre( deformation.U, deformation.surfaceGradient, normal, deformation.connector );

  double* baseMaterialStateVars = stateLayout.getPtr( state.stateVars, "baseMaterialStateVars" );

  MarmotMaterialFiniteStrain::ConstitutiveResponse< 3 > response( Tensor33d( 0.0 ), 0.0, 0.0, baseMaterialStateVars );
  MarmotMaterialFiniteStrain::AlgorithmicModuli< 3 >    materialTangents;
  MarmotMaterialFiniteStrain::Deformation< 3 >          materialDeformation{ fibre.F };
  const MarmotMaterialFiniteStrain::TimeIncrement materialTimeIncrement{ timeIncrement.timeOld, timeIncrement.dT };

  baseMaterial->computeStress( response, materialTangents, materialDeformation, materialTimeIncrement );

  Tensor33d   P( 0.0 );
  Tensor3333d A( 0.0 );
  firstPiolaKirchhoff( response.tau, materialTangents.dTau_dF, fibre.F, P, A );

  const ReducedResponse reduced = reduceToInterface( P, A, normal, fibre, h );

  state.force         = reduced.force;
  state.surfaceStress = reduced.surfaceStress;

  tangents.Q_ik   = reduced.Q_ik;
  tangents.H_ikR  = reduced.H_ikR;
  tangents.H_iRk  = reduced.H_iRk;
  tangents.A_iRkP = reduced.A_iRkP;
}

void MarmotInterfaceMaterialFiniteStrain::initializeYourself( double* stateVars, int )
{
  baseMaterial->initializeYourself( stateLayout.getPtr( stateVars, "baseMaterialStateVars" ),
                                    baseMaterial->getNumberOfRequiredStateVars() );
}

double MarmotInterfaceMaterialFiniteStrain::getDensity() const
{
  return baseMaterial->getDensity( nullptr );
}
