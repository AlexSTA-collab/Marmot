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
  if ( !baseMaterial )
    throw std::invalid_argument( "Unknown base material for MarmotInterfaceMaterialFiniteStrain: " + materialName );

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
  using namespace Fastor;

  const Tensor3d& normal    = deformation.normal;
  const Tensor3d& connector = deformation.connector;

  const Tensor3d  jumpU                      = deformation.U( seq( 0, 3 ) ) - deformation.U( seq( 3, last ) );
  const Tensor9d  averageSurfaceGradientFlat = 0.5 * ( deformation.surfaceGradient( seq( 0, 9 ) ) +
                                                      deformation.surfaceGradient( seq( 9, last ) ) );
  const Tensor33d averageSurfaceGradient( reshape< 3, 3 >( averageSurfaceGradientFlat ) );

  // The reduced normal derivative is read off the connector between the paired faces, and not off the jump divided
  // by the constitutive thickness. With dX = ell N + dTau, the tangential offset dTau is carried by the surface map
  // rather than by the jump, so that the image of the unit normal fibre is
  //
  //     gbar = ( [u] - <grad_s u> dTau ) / ell .
  //
  // Leaving the correction out costs the frame indifference of an inclined interface: a rigid rotation of the body
  // would cause stresses. For a purely normal connector it reduces to [u] / h.
  const double ell = inner( connector, normal );

  if ( ell <= 0.0 )
    throw std::invalid_argument(
      "MarmotInterfaceMaterialFiniteStrain: the paired faces of the interface element must be separated along the "
      "interface normal, and the normal must point from the first side to the second. A coincident-face mesh "
      "carries no normal fibre to differentiate along." );

  const Tensor3d dTau = connector - ell * normal;

  const Tensor3d normalGradient = ( 1.0 / ell ) * ( jumpU - einsum< ij, j, to_i >( averageSurfaceGradient, dTau ) );

  Tensor33d F( 0.0 );
  F.eye();
  F += averageSurfaceGradient + einsum< i, j, to_ij >( normalGradient, normal );

  double* baseMaterialStateVars = stateLayout.getPtr( state.stateVars, "baseMaterialStateVars" );

  MarmotMaterialFiniteStrain::ConstitutiveResponse< 3 > response( Tensor33d( 0.0 ), 0.0, 0.0, baseMaterialStateVars );
  MarmotMaterialFiniteStrain::AlgorithmicModuli< 3 >    materialTangents;
  MarmotMaterialFiniteStrain::Deformation< 3 >          materialDeformation{ F };
  const MarmotMaterialFiniteStrain::TimeIncrement materialTimeIncrement{ timeIncrement.timeOld, timeIncrement.dT };

  baseMaterial->computeStress( response, materialTangents, materialDeformation, materialTimeIncrement );

  const Tensor33d FInvT = transpose( inverse( F ) );
  const Tensor33d P     = einsum< ij, jk, to_ik >( response.tau, FInvT );

  Tensor33d T( 0.0 );
  T.eye();
  T -= einsum< i, j, to_ij >( normal, normal );

  const Tensor3333d A = einsum< imkl, mj, to_ijkl >( materialTangents.dTau_dF, FInvT ) -
                        einsum< im, kj, ml, to_ijkl >( response.tau, FInvT, FInvT );

  // Work conjugates of the reduced kinematics. With
  //
  //     F = I + <grad_s u> + ( [u] - <grad_s u> dTau ) (x) N / ell ,
  //
  // the internal work per unit reference midsurface area is h P : dF, so that the jump is conjugate to
  //
  //     t      = ( h / ell ) P N                                            and
  //     Sigma  = h P T - ( h / ell ) ( P N ) (x) dTau
  //
  // is conjugate to the average surface gradient. The thicknesses h and ell coincide whenever the faces are
  // separated by the layer; ell is mesh geometry, h is the constitutive layer.
  const double hOverEll = h / ell;

  const Tensor3d    traction = einsum< ij, j, to_i >( P, normal );
  const Tensor33d   ANN      = einsum< ijkl, j, l, to_ik >( A, normal, normal );
  const Tensor333d  ANT      = einsum< ijkl, j, ml, to_ikm >( A, normal, T );
  const Tensor333d  ATN      = einsum< ijkl, mj, l, to_imk >( A, T, normal );
  const Tensor3333d ATT      = einsum< ijkl, mj, nl, to_imkn >( A, T, T );

  // ANN with the tangential offset appended, in both index orders the tangent blocks need
  const Tensor333d ANNdTauJumpAverage = einsum< ik, m, to_ikm >( ANN, dTau );
  const Tensor333d ANNdTauAverageJump = einsum< ik, m, to_imk >( ANN, dTau );

  state.force         = hOverEll * traction;
  state.surfaceStress = h * einsum< ij, kj, to_ik >( P, T ) - hOverEll * einsum< i, k, to_ik >( traction, dTau );

  tangents.Q_ik   = ( hOverEll / ell ) * ANN;
  tangents.H_ikR  = hOverEll * ( ANT - ( 1.0 / ell ) * ANNdTauJumpAverage );
  tangents.H_iRk  = hOverEll * ( ATN - ( 1.0 / ell ) * ANNdTauAverageJump );
  tangents.A_iRkP = h * ATT -
                    hOverEll * ( einsum< imk, n, to_imkn >( ATN, dTau ) + einsum< m, ikn, to_imkn >( dTau, ANT ) ) +
                    ( hOverEll / ell ) * einsum< imk, n, to_imkn >( ANNdTauAverageJump, dTau );
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
