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
#include "Marmot/MarmotInterfaceMaterialFiniteStrain.h"
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrain.h"
#include "Marmot/MarmotStateHelpers.h"

#include <Eigen/Dense>

#include <memory>
#include <string>
#include <vector>

/**
 * @brief Gradient-enhanced finite-strain interface material with a nonlocal field at every station through the
 * thickness of the layer.
 *
 * This is the gradient-enhanced counterpart of MarmotGradedInterfaceMaterialFiniteStrain. The layer has the
 * through-thickness coordinate @f$y \in [0,1]@f$, and its bulk material is a
 * MarmotMaterialGradientEnhancedFiniteStrain at every station, with a nonlocal field @f$ \bar N_q @f$ of its own,
 * so that the damage, which the nonlocal field drives, varies through the thickness. The stations are the
 * integration points of the thickness, with the weights of the trapezoidal rule (the mode 2 of
 * Marmot::Materials::ThicknessProfile), so that every point sees the field of its station. For two stations at the
 * faces the material is the one of MarmotEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain with a nonlocal
 * field for each face.
 *
 * The kinematics, the microproblem and the reduction to the interface are those of
 * MarmotGradedInterfaceMaterialFiniteStrain, with the stress of the point @f$ P_q(F_q, \bar N_q) @f$. The
 * microproblem is solved for given nonlocal fields. Since the fields enter the stress of their point, the
 * corrections @f$ w_q @f$ depend on them, and so do the force and the surface stress. The material returns, in
 * addition to the response of the interface,
 *
 * - the local driving field @f$ L_q @f$ and the nonlocal radius @f$ l_q @f$ of every station, and
 * - the derivatives of the force, of the surface stress and of the local fields with respect to the jump, to the
 *   average surface gradient and to the nonlocal fields, all condensed over the corrections.
 *
 * The material properties are the ones of Marmot::Materials::ThicknessProfile with the mode 2: @f$[h, 2, 1, K, m,
 * (p_k, \text{properties}) \times K]@f$, where the properties are the ones of the bulk material.
 *
 * The state variables are, for every station, the correction @f$w@f$ of its last converged state, which starts the
 * next iteration, followed by the state variables of its bulk material.
 */
class MarmotGradedGradientEnhancedInterfaceMaterialFiniteStrain {

public:
  using Tensor3d    = MarmotInterfaceMaterialFiniteStrain::Tensor3d;    ///< first order tensor in 3D
  using Tensor33d   = MarmotInterfaceMaterialFiniteStrain::Tensor33d;   ///< second order tensor in 3D
  using Tensor333d  = MarmotInterfaceMaterialFiniteStrain::Tensor333d;  ///< third order tensor in 3D
  using Tensor3333d = MarmotInterfaceMaterialFiniteStrain::Tensor3333d; ///< fourth order tensor in 3D
  using Tensor6d    = MarmotInterfaceMaterialFiniteStrain::Tensor6d;    ///< vector of six components, one per side
  using Tensor18d   = MarmotInterfaceMaterialFiniteStrain::Tensor18d;   ///< vector of 18 components, one per side

  using TimeIncrement = MarmotInterfaceMaterialFiniteStrain::TimeIncrement; ///< time data, as in the base

  /// number of the material
  const int materialNumber;

  /// Layout of the state variables
  MarmotStateLayoutDynamic stateLayout;

  /// Interface response: the updated values on exit.
  struct State {
    Tensor3d        force;         ///< force, conjugate to the jump of the displacements
    Tensor33d       surfaceStress; ///< surface stress, conjugate to the average surface gradient
    double*         stateVars;     ///< pointer to the state variables
    Eigen::VectorXd localFields;   ///< local driving field of every station
    Eigen::VectorXd nonlocalRadii; ///< nonlocal radius of every station
  };

  /// Algorithmic tangents, set by computeStress().
  struct Tangents {
    Tensor33d   Q_ik;   ///< derivative of the force with respect to the jump
    Tensor333d  H_ikR;  ///< derivative of the force with respect to the average surface gradient
    Tensor333d  H_iRk;  ///< derivative of the surface stress with respect to the jump
    Tensor3333d A_iRkP; ///< derivative of the surface stress with respect to the average surface gradient

    /// derivative of the force with respect to the nonlocal field of every station, 3 x K
    Eigen::MatrixXd dForce_dNonlocal;
    /// derivative of the surface stress with respect to the nonlocal field of every station, 9 x K, row major
    Eigen::MatrixXd dSurfaceStress_dNonlocal;
    /// derivative of the local field of every station with respect to the jump, K x 3
    Eigen::MatrixXd dLocal_dJump;
    /// derivative of the local field of every station with respect to the average surface gradient, K x 9
    Eigen::MatrixXd dLocal_dAverage;
    /// derivative of the local field of every station with respect to the nonlocal field of every station, K x K
    Eigen::MatrixXd dLocal_dNonlocal;
  };

  /// Kinematics of the interface at a quadrature point, and the nonlocal field of every station.
  struct Deformation {
    Tensor6d        U;               ///< total displacements of the two interface sides, top first
    Tensor18d       surfaceGradient; ///< total surface gradients of the two interface sides, top first
    Tensor3d        normal;          ///< interface normal, pointing from the first to the second side
    Tensor3d        connector;       ///< reference connector between the paired faces, X+ - X-
    Eigen::VectorXd nonlocalFields;  ///< nonlocal field of every station, from the bottom to the top
  };

  /**
   * @brief Construct the material from a registered gradient-enhanced finite-strain bulk material.
   *
   * @param[in] materialName name of the bulk material, registered in the
   * MarmotMaterialGradientEnhancedFiniteStrainFactory
   * @param[in] matProperties material properties, as described for this class
   * @param[in] nMaterialProperties number of material properties
   * @param[in] materialNumber_ number of the material
   *
   * @throws std::invalid_argument if the properties are inconsistent, if the mode is not 2, or if the bulk material is
   * not registered
   */
  MarmotGradedGradientEnhancedInterfaceMaterialFiniteStrain( const std::string& materialName,
                                                             const double*      matProperties,
                                                             int                nMaterialProperties,
                                                             int                materialNumber_ );

  /// Default destructor
  virtual ~MarmotGradedGradientEnhancedInterfaceMaterialFiniteStrain() = default;

  /**
   * @brief Set the characteristic element length at the considered quadrature point.
   *
   * @param[in] length characteristic length, handed to all bulk materials
   */
  void setCharacteristicElementLength( double length );

  /**
   * @brief Solve the microproblem of the layer and compute the interface response and the nonlocal couplings.
   *
   * @param[in,out] state interface response; the state variables are updated
   * @param[out] tangents tangents of the interface and of the coupling to the nonlocal fields
   * @param[in] deformation total displacements and surface gradients of both sides, the normal, the connector and
   * the nonlocal fields of the stations
   * @param[in] timeIncrement time at the beginning of the increment and time increment
   *
   * @throws std::invalid_argument if the number of nonlocal fields is not the number of stations
   * @throws Marmot::StressUpdateFailed if the iteration for the microproblem does not converge, or if the acoustic
   * tensor of a station is singular
   */
  virtual void computeStress( State&               state,
                              Tangents&            tangents,
                              const Deformation&   deformation,
                              const TimeIncrement& timeIncrement );

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
   * @return number of state variables of all stations
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
   * @param[in] stateVars state variables of the material point
   * @return through-thickness mean of the density of the bulk materials
   */
  virtual double getDensity( const double* stateVars ) const;

  /**
   * @brief Number of stations through the thickness, which is the number of nonlocal fields.
   *
   * @return the number of stations
   */
  int numberOfStations() const { return static_cast< int >( points.size() ); }

  /**
   * @brief Thickness of the layer.
   *
   * @return the constitutive thickness h
   */
  double thickness() const { return h; }

  /**
   * @brief Weight of a station in the integration through the thickness.
   *
   * @param[in] station index of the station, from the bottom to the top
   * @return the weight, the weights of all stations sum up to one
   */
  double stationWeight( int station ) const { return points[station].weight; }

  /**
   * @brief Position of a station in the thickness.
   *
   * @param[in] station index of the station, from the bottom to the top
   * @return the coordinate y in [0, 1]
   */
  double stationPosition( int station ) const { return points[station].position; }

private:
  /// Station through the thickness.
  struct Station {
    double                                                        weight;     ///< weight, the weights sum up to one
    double                                                        position;   ///< coordinate y in [0, 1]
    std::vector< double >                                         properties; ///< properties of the bulk material
    std::unique_ptr< MarmotMaterialGradientEnhancedFiniteStrain > material;   ///< bulk material of the station
    std::size_t                                                   offset;     ///< offset of the state of the station
    std::size_t nBaseStateVars; ///< number of state variables of the bulk material
  };

  const double           h;      ///< thickness of the layer
  std::vector< Station > points; ///< stations through the thickness, from the bottom to the top
};
