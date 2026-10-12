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

#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotStateHelpers.h"

#include <Eigen/Dense>

#include <memory>
#include <string>
#include <vector>

/**
 * @brief Integration points through the thickness of a layer with finite-strain bulk materials.
 *
 * This is the part that all graded finite-strain interface materials share (see
 * MarmotGradedInterfaceMaterialFiniteStrain and MarmotGradedEquilibratedInterfaceMaterialFiniteStrain): the
 * integration points through the thickness read from the material properties (see
 * Marmot::Materials::ThicknessProfile), one bulk material and its state at every point, the layout of the state
 * variables and the update of a single point.
 *
 * The state variables are, for every integration point through the thickness, the correction @f$w@f$ of its last
 * converged state, which starts the next iteration, followed by the state variables of its bulk material.
 */
class MarmotGradedInterfaceMaterialFiniteStrainBase {

public:
  /// number of the material
  const int materialNumber;

  /// Layout of the state variables
  MarmotStateLayoutDynamic stateLayout;

  /**
   * @brief Construct the integration points through the thickness from a registered finite-strain bulk material.
   *
   * @param[in] materialName name of the bulk material, registered in the MarmotMaterialFiniteStrainFactory
   * @param[in] matProperties material properties, see Marmot::Materials::ThicknessProfile
   * @param[in] nMaterialProperties number of material properties
   * @param[in] materialNumber_ number of the material
   *
   * @throws std::invalid_argument if the properties are inconsistent or the bulk material is not registered
   */
  MarmotGradedInterfaceMaterialFiniteStrainBase( const std::string& materialName,
                                                 const double*      matProperties,
                                                 int                nMaterialProperties,
                                                 int                materialNumber_ );

  /// Default destructor
  virtual ~MarmotGradedInterfaceMaterialFiniteStrainBase() = default;

  /**
   * @brief Set the characteristic element length at the considered quadrature point.
   *
   * @param[in] length characteristic length, handed to all bulk materials
   */
  void setCharacteristicElementLength( double length );

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
  virtual void initializeYourself( double* stateVars, int nStateVars );

  /**
   * @brief Get the density of the material.
   *
   * @return through-thickness mean of the density of the bulk materials
   */
  virtual double getDensity() const;

  /**
   * @brief Number of integration points through the thickness.
   *
   * @return the number of integration points
   */
  int numberOfThicknessPoints() const { return static_cast< int >( points.size() ); }

protected:
  /// Integration point through the thickness.
  struct ThicknessPoint {
    double                                        weight;         ///< weight, the weights sum up to one
    double                                        position;       ///< coordinate y in [0, 1]
    std::vector< double >                         properties;     ///< properties of the bulk material
    std::unique_ptr< MarmotMaterialFiniteStrain > material;       ///< bulk material at the position
    std::size_t                                   offset;         ///< offset of the state of this point
    std::size_t                                   nBaseStateVars; ///< number of state variables of the bulk material
  };

  /// Response of a bulk material at an integration point to a deformation gradient.
  struct BulkResponse {
    Marmot::FastorStandardTensors::Tensor33d   P; ///< first Piola-Kirchhoff stress
    Marmot::FastorStandardTensors::Tensor3333d A; ///< derivative of the stress with respect to F
    std::vector< double > state; ///< correction and state variables of the bulk material after the update
  };

  /**
   * @brief Update the bulk material of a point for a deformation gradient, starting from the committed state.
   *
   * The state variables passed are not modified, the updated ones are returned in the response, where the entries
   * of the correction are the committed ones.
   *
   * @param[in] q index of the integration point
   * @param[in] stateVars state variables of the material point at the beginning of the increment
   * @param[in] F total deformation gradient of the point
   * @param[in] timeIncrement time at the beginning of the increment and time increment
   * @return stress, tangent and updated state of the point
   */
  BulkResponse evaluateBulk( int                                              q,
                             const double*                                    stateVars,
                             const Marmot::FastorStandardTensors::Tensor33d&  F,
                             const MarmotMaterialFiniteStrain::TimeIncrement& timeIncrement ) const;

  const double                  h;                  ///< thickness of the layer
  bool                          tangentialGradient; ///< the tangential gradient of the correction is kept
  std::vector< ThicknessPoint > points;             ///< integration points through the thickness, bottom to top

  /// weights of the cumulative integral from the bottom to every point, see Marmot::Materials::ThicknessProfile
  Eigen::MatrixXd cumulativeWeights;
};
