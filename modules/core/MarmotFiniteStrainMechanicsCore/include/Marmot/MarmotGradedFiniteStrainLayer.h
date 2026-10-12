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

#pragma once

#include "Marmot/MarmotStateHelpers.h"
#include "Marmot/MarmotThicknessProfile.h"
#include "Marmot/MarmotUtils.h"

#include <Eigen/Dense>

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

/**
 * @brief Integration points through the thickness of a layer with finite-strain bulk materials of a given kind.
 *
 * This is the part that all graded finite-strain interface materials share: the integration points through the
 * thickness read from the material properties (see Marmot::Materials::ThicknessProfile), one bulk material and its
 * state at every point, and the layout of the state variables. The bulk materials are finite-strain materials or
 * gradient-enhanced finite-strain materials, which differ in the evaluation of a point, so that this is done by the
 * derived classes.
 *
 * The state variables are, for every integration point through the thickness, the correction @f$w@f$ of its last
 * converged state, which starts the next iteration, followed by the state variables of its bulk material.
 *
 * @tparam BulkMaterial type of the bulk material
 * @tparam Factory factory of the bulk materials, with `createMaterial( name, properties, n, materialNumber )`
 */
template < class BulkMaterial, class Factory >
class MarmotGradedFiniteStrainLayer {

public:
  /// number of the material
  const int materialNumber;

  /// Layout of the state variables
  MarmotStateLayoutDynamic stateLayout;

  /**
   * @brief Construct the integration points through the thickness from a registered bulk material.
   *
   * @param[in] materialName name of the bulk material, registered in the factory
   * @param[in] matProperties material properties, see Marmot::Materials::ThicknessProfile
   * @param[in] nMaterialProperties number of material properties
   * @param[in] materialNumber_ number of the material
   *
   * @throws std::invalid_argument if the properties are inconsistent or the bulk material is not registered
   */
  MarmotGradedFiniteStrainLayer( const std::string& materialName,
                                 const double*      matProperties,
                                 int                nMaterialProperties,
                                 int                materialNumber_ )
    : materialNumber( materialNumber_ ), h( Marmot::checkedMaterialProperty( matProperties, nMaterialProperties, 0 ) )
  {
    auto profile = Marmot::Materials::ThicknessProfile::fromMaterialProperties( matProperties, nMaterialProperties );

    tangentialGradient = profile.tangentialGradient;

    cumulativeWeights.resize( profile.points.size(), profile.points.size() );
    for ( std::size_t q = 0; q < profile.points.size(); q++ )
      for ( std::size_t p = 0; p < profile.points.size(); p++ )
        cumulativeWeights( q, p ) = profile.cumulativeWeights[q][p];

    std::size_t offset = 0;
    for ( auto& profilePoint : profile.points ) {
      ThicknessPoint point;
      point.weight     = profilePoint.weight;
      point.position   = profilePoint.position;
      point.properties = std::move( profilePoint.properties );

      point.material = std::unique_ptr< BulkMaterial >(
        Factory::createMaterial( materialName,
                                 point.properties.data(),
                                 static_cast< int >( point.properties.size() ),
                                 materialNumber ) );
      if ( !point.material )
        throw std::invalid_argument( "Unknown base material for a graded finite-strain interface material: " +
                                     materialName );

      point.nBaseStateVars = point.material->getNumberOfRequiredStateVars();
      point.offset         = offset;
      offset += 3 + point.nBaseStateVars;
      points.push_back( std::move( point ) );
    }

    for ( std::size_t q = 0; q < points.size(); q++ ) {
      stateLayout.add( "correction" + std::to_string( q ), 3 );
      stateLayout.add( "baseMaterialStateVars" + std::to_string( q ), points[q].nBaseStateVars );
    }
    stateLayout.finalize();
  }

  /// Default destructor
  virtual ~MarmotGradedFiniteStrainLayer() = default;

  /**
   * @brief Set the characteristic element length at the considered quadrature point.
   *
   * @param[in] length characteristic length, handed to all bulk materials
   */
  void setCharacteristicElementLength( double length )
  {
    for ( auto& point : points )
      point.material->setCharacteristicElementLength( length );
  }

  /**
   * @brief Get a view to the state variables.
   *
   * @param[in] stateName name of the state variable
   * @param[in] stateVars pointer to the state variable array
   * @return view to the requested state variable
   */
  StateView getStateView( const std::string& stateName, double* stateVars ) const
  {
    return stateLayout.getStateView( stateVars, stateName );
  }

  /**
   * @brief Get the total number of required state variables.
   *
   * @return number of state variables of all integration points through the thickness
   */
  virtual int getNumberOfRequiredStateVars() const { return stateLayout.totalSize(); }

  /**
   * @brief Initialize the state variables: zero corrections and initial state of every bulk material.
   *
   * @param[in,out] stateVars state variables of the material point
   * @param[in] nStateVars number of state variables, which this material does not need
   */
  virtual void initializeYourself( double* stateVars, int nStateVars )
  {
    (void)nStateVars;
    for ( auto& point : points ) {
      std::fill( stateVars + point.offset, stateVars + point.offset + 3, 0.0 );
      point.material->initializeYourself( stateVars + point.offset + 3, static_cast< int >( point.nBaseStateVars ) );
    }
  }

  /**
   * @brief Get the density of the material.
   *
   * @param[in] stateVars state variables of the material point, or nullptr
   * @return through-thickness mean of the density of the bulk materials
   */
  double getDensity( const double* stateVars ) const
  {
    double density = 0.;
    for ( const auto& point : points )
      density += point.weight * point.material->getDensity( stateVars ? stateVars + point.offset + 3 : nullptr );
    return density;
  }

  /**
   * @brief Number of integration points through the thickness.
   *
   * @return the number of integration points
   */
  int numberOfThicknessPoints() const { return static_cast< int >( points.size() ); }

  /**
   * @brief Thickness of the layer.
   *
   * @return the constitutive thickness h
   */
  double thickness() const { return h; }

  /**
   * @brief Weight of an integration point in the integration through the thickness.
   *
   * @param[in] index index of the point, from the bottom to the top
   * @return the weight, the weights of all points sum up to one
   */
  double stationWeight( int index ) const { return points[index].weight; }

  /**
   * @brief Position of an integration point in the thickness.
   *
   * @param[in] index index of the point, from the bottom to the top
   * @return the coordinate y in [0, 1]
   */
  double stationPosition( int index ) const { return points[index].position; }

protected:
  /// Integration point through the thickness.
  struct ThicknessPoint {
    double                          weight;         ///< weight, the weights sum up to one
    double                          position;       ///< coordinate y in [0, 1]
    std::vector< double >           properties;     ///< properties of the bulk material
    std::unique_ptr< BulkMaterial > material;       ///< bulk material at the position
    std::size_t                     offset;         ///< offset of the state of this point
    std::size_t                     nBaseStateVars; ///< number of state variables of the bulk material
  };

  const double                  h;                  ///< thickness of the layer
  bool                          tangentialGradient; ///< the tangential gradient of the correction is kept
  std::vector< ThicknessPoint > points;             ///< integration points through the thickness, bottom to top

  /// weights of the cumulative integral from the bottom to every point, see Marmot::Materials::ThicknessProfile
  Eigen::MatrixXd cumulativeWeights;
};
