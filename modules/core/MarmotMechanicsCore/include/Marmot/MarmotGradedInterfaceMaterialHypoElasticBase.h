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

#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotStateHelpers.h"
#include "Marmot/MarmotTypedefs.h"

#include <Eigen/Dense>

#include <memory>
#include <string>
#include <vector>

/**
 * @brief Integration points through the thickness of a layer with hypoelastic bulk materials.
 *
 * This is the part that all graded hypoelastic interface materials share (see
 * MarmotGradedInterfaceMaterialHypoElastic and MarmotGradedXInterfaceMaterialHypoElastic): the integration points
 * through the thickness read from the material properties (see Marmot::Materials::ThicknessProfile), one bulk
 * material and its state at every point, the layout of the state variables and the update of a single point.
 *
 * The state variables are, for every integration point through the thickness, the stress in Voigt notation followed
 * by the state variables of its bulk material.
 */
class MarmotGradedInterfaceMaterialHypoElasticBase {

public:
  /// number of the material
  const int materialNumber;

  /// Layout of the state variables
  MarmotStateLayoutDynamic stateLayout;

  /// Characteristic element length
  double characteristicElementLength = 0.0;

  /**
   * @brief Construct the integration points through the thickness from a registered hypoelastic bulk material.
   *
   * @param[in] materialName name of the bulk material, registered in the MarmotMaterialHypoElasticFactory
   * @param[in] matProperties material properties, see Marmot::Materials::ThicknessProfile
   * @param[in] nMaterialProperties number of material properties
   * @param[in] materialNumber_ number of the material
   *
   * @throws std::invalid_argument if the properties are inconsistent, for example a non-positive thickness or an
   * unknown mode, if the positions of the stations are not increasing from 0 to 1, or if the bulk material is not
   * registered
   */
  MarmotGradedInterfaceMaterialHypoElasticBase( const std::string& materialName,
                                                const double*      matProperties,
                                                int                nMaterialProperties,
                                                int                materialNumber_ );

  /// Default destructor
  virtual ~MarmotGradedInterfaceMaterialHypoElasticBase() = default;

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
   * @brief Initialize the state variables: zero stress and initial state of every bulk material.
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
  virtual double getDensity();

  /**
   * @brief Number of integration points through the thickness.
   *
   * @return the number of integration points
   */
  int numberOfThicknessPoints() const { return static_cast< int >( points.size() ); }

protected:
  /**
   * @brief Operator S that maps a vector w to the strain sym( w (x) n ) in Voigt notation.
   *
   * With the engineering shear strains of the Voigt notation, the transpose of S maps a stress in Voigt notation
   * to the traction sigma n.
   *
   * @param[in] n the normal
   * @return the 6x3 operator
   */
  static Eigen::Matrix< double, 6, 3 > normalStrainOperator( const Eigen::Vector3d& n );

  /**
   * @brief Operator that maps a displacement gradient, stored row by row in a 9-vector, to its strain in Voigt
   * notation.
   *
   * @return the 6x9 operator
   */
  static Eigen::Matrix< double, 6, 9 > gradientToVoigtStrain();

  /**
   * @brief Operator that maps a stress in Voigt notation to the row by row 9-vector of the stress tensor.
   *
   * @return the 9x6 operator
   */
  static Eigen::Matrix< double, 9, 6 > voigtStressToVector();

  /// Integration point through the thickness.
  struct ThicknessPoint {
    double                                       weight;     ///< weight, the weights of all points sum up to one
    double                                       position;   ///< coordinate y in [0, 1]
    std::vector< double >                        properties; ///< properties of the bulk material, which it refers to
    std::unique_ptr< MarmotMaterialHypoElastic > material;   ///< bulk material at the position
    std::size_t                                  offset; ///< offset of the state of this point in the state variables
    std::size_t                                  nBaseStateVars; ///< number of state variables of the bulk material
  };

  /// Response of an integration point through the thickness to the strain of a trial value of the correction w.
  struct PointEvaluation {
    Eigen::Vector3d       w = Eigen::Vector3d::Zero(); ///< first-order strain correction w = du1/dy
    Marmot::Vector6d      sigma;                       ///< stress
    Marmot::Matrix6d      D;                           ///< tangent
    Eigen::Matrix3d       Qinv;                        ///< inverse of the acoustic tensor
    Eigen::Vector3d       a;                           ///< traction of the point, which is equal at all points
    std::vector< double > state;                       ///< stress and state variables after the update
  };

  /// Response of a bulk material at an integration point to a strain increment.
  struct BulkResponse {
    Marmot::Vector6d      sigma; ///< stress
    Marmot::Matrix6d      D;     ///< tangent
    std::vector< double > state; ///< stress and state variables after the update
  };

  /**
   * @brief Update the bulk material of a point for a strain increment, starting from the committed state.
   *
   * The state variables passed are not modified, the updated ones are returned in the response.
   *
   * @param[in] q index of the integration point
   * @param[in] stateVars state variables of the material point at the beginning of the increment
   * @param[in] strainIncrement strain increment of the point in Voigt notation
   * @param[in] timeInfo time at the end of the increment and time increment
   * @return stress, tangent and updated state of the point
   */
  BulkResponse evaluateBulk( int                                        q,
                             const double*                              stateVars,
                             const Marmot::Vector6d&                    strainIncrement,
                             const MarmotMaterialHypoElastic::timeInfo& timeInfo ) const;

  const double                  h;                  ///< thickness of the layer
  bool                          tangentialGradient; ///< the tangential gradient of the correction is kept
  std::vector< ThicknessPoint > points; ///< integration points through the thickness, from the bottom to the top

  /// weights of the cumulative integral from the bottom to every point, see Marmot::Materials::ThicknessProfile
  Eigen::MatrixXd cumulativeWeights;
};
