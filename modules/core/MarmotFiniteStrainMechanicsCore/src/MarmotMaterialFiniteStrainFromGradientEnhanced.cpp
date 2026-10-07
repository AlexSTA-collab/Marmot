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

#include "Marmot/MarmotMaterialFiniteStrainFromGradientEnhanced.h"

#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrainFactory.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

MarmotMaterialFiniteStrainFromGradientEnhanced::MarmotMaterialFiniteStrainFromGradientEnhanced(
  std::unique_ptr< MarmotMaterialGradientEnhancedFiniteStrain > wrappedMaterial,
  const double*                                                 matProperties_,
  int                                                           nMaterialProperties_,
  int                                                           materialNumber_,
  double                                                        nonlocalFieldValue )
  : MarmotMaterialFiniteStrain( matProperties_, nMaterialProperties_, materialNumber_ ),
    wrapped( std::move( wrappedMaterial ) ),
    nonlocalField( nonlocalFieldValue )
{
  if ( !wrapped )
    throw std::invalid_argument(
      "MarmotMaterialFiniteStrainFromGradientEnhanced: no gradient-enhanced material to wrap." );

  stateLayout.add( "gradientEnhancedMaterialStateVars", wrapped->getNumberOfRequiredStateVars() );
  stateLayout.finalize();
}

void MarmotMaterialFiniteStrainFromGradientEnhanced::computeStress( ConstitutiveResponse< 3 >& response,
                                                                    AlgorithmicModuli< 3 >&    tangents,
                                                                    const Deformation< 3 >&    deformation,
                                                                    const TimeIncrement&       timeIncrement ) const
{
  MarmotMaterialGradientEnhancedFiniteStrain::ConstitutiveResponse< 3 > wrappedResponse;
  wrappedResponse.dissipation = response.dissipation;
  wrappedResponse.stateVars   = response.stateVars;

  MarmotMaterialGradientEnhancedFiniteStrain::AlgorithmicModuli< 3 > wrappedTangents;

  const MarmotMaterialGradientEnhancedFiniteStrain::Deformation< 3 > wrappedDeformation{ deformation.F, nonlocalField };
  const MarmotMaterialGradientEnhancedFiniteStrain::TimeIncrement wrappedTime{ timeIncrement.time, timeIncrement.dT };

  wrapped->computeStress( wrappedResponse, wrappedTangents, wrappedDeformation, wrappedTime );

  // The coupling to the nonlocal field has nowhere to go here, so it must be absent rather than discarded. A
  // material whose stress responds to the nonlocal field would get a wrong tangent and converge to a wrong state,
  // silently.
  double maximumCoupling = 0.0;
  double maximumTangent  = 0.0;

  for ( int i = 0; i < 3; ++i )
    for ( int j = 0; j < 3; ++j ) {
      maximumCoupling = std::max( maximumCoupling, std::abs( wrappedTangents.dTau_dN( i, j ) ) );

      for ( int k = 0; k < 3; ++k )
        for ( int l = 0; l < 3; ++l )
          maximumTangent = std::max( maximumTangent, std::abs( wrappedTangents.dTau_dF( i, j, k, l ) ) );
    }

  if ( maximumCoupling > 1.0e-12 * std::max( 1.0, maximumTangent ) )
    throw std::invalid_argument(
      MakeString() << "MarmotMaterialFiniteStrainFromGradientEnhanced: the stress of the wrapped material depends on "
                      "the nonlocal field (max |dTau_dN| = "
                   << maximumCoupling
                   << "), so it cannot be driven by an element that carries only displacements. Make the material "
                      "local, or use a gradient-enhanced element." );

  response.tau                  = wrappedResponse.tau;
  response.elasticEnergyDensity = wrappedResponse.elasticEnergyDensity;
  response.dissipation          = wrappedResponse.dissipation;

  tangents.dTau_dF = wrappedTangents.dTau_dF;
}

void MarmotMaterialFiniteStrainFromGradientEnhanced::initializeYourself( double* stateVars, int nStateVars )
{
  wrapped->initializeYourself( stateVars, nStateVars );
}

double MarmotMaterialFiniteStrainFromGradientEnhanced::getDensity( const double* stateVars ) const
{
  return wrapped->getDensity( stateVars );
}

void MarmotMaterialFiniteStrainFromGradientEnhanced::setCharacteristicElementLength( double length )
{
  wrapped->setCharacteristicElementLength( length );
}

std::unique_ptr< MarmotMaterialFiniteStrain > MarmotLibrary::createFiniteStrainMaterialAllowingGradientEnhanced(
  const std::string& materialName,
  const double*      materialProperties,
  int                nMaterialProperties,
  int                materialNumber )
{
  std::string notFiniteStrain;

  try {
    return std::unique_ptr< MarmotMaterialFiniteStrain >(
      MarmotMaterialFiniteStrainFactory::createMaterial( materialName,
                                                         materialProperties,
                                                         nMaterialProperties,
                                                         materialNumber ) );
  }
  catch ( const std::invalid_argument& error ) {
    // not a plain finite-strain material; it may still be a gradient-enhanced one
    notFiniteStrain = error.what();
  }

  std::unique_ptr< MarmotMaterialGradientEnhancedFiniteStrain > gradientEnhanced;

  try {
    gradientEnhanced.reset( MarmotMaterialGradientEnhancedFiniteStrainFactory::createMaterial( materialName,
                                                                                               materialProperties,
                                                                                               nMaterialProperties,
                                                                                               materialNumber ) );
  }
  catch ( const std::invalid_argument& error ) {
    throw std::invalid_argument( MakeString() << notFiniteStrain << "\n" << error.what() );
  }

  return std::make_unique< MarmotMaterialFiniteStrainFromGradientEnhanced >( std::move( gradientEnhanced ),
                                                                             materialProperties,
                                                                             nMaterialProperties,
                                                                             materialNumber );
}
