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
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotStateHelpers.h"

#include <memory>
#include <string>
#include <vector>

/**
 * @brief Finite-strain interface material with a stiffness that varies through the thickness of the layer.
 *
 * The layer has the thickness @f$h@f$ and the through-thickness coordinate @f$y \in [0,1]@f$, which starts at the
 * second side of the interface and points along the normal @f$N@f$ to the first side. The bulk material at the
 * height @f$y@f$ is a finite-strain material with its own properties and its own history. It is the finite-strain
 * counterpart of MarmotGradedInterfaceMaterialHypoElastic.
 *
 * The kinematics of the interface are those of MarmotInterfaceMaterialFiniteStrain, which give the deformation
 * gradient of the layer
 * @f[
 *   F^0 = I + \langle \nabla_s u \rangle + \bar{g} \otimes N ,
 * @f]
 * with the reduced normal derivative @f$\bar{g}@f$ of the displacement read off the connector between the paired
 * faces. The deformation gradient at the height @f$y@f$ is this one plus the correction of the microproblem,
 * @f[
 *   F(y) = F^0 + w(y) \otimes N , \qquad w = \partial_Y u^1 ,
 * @f]
 * where @f$u^1@f$ is the first-order displacement of the layer. The microproblem follows from the variational form
 * with a vanishing macroscopic test function. It requires that the traction
 * @f[
 *   t = P(y) \, N
 * @f]
 * is constant through the thickness and that @f$u^1@f$ vanishes at both sides, hence that
 * @f$\int_0^1 w \, \mathrm{d}y = 0@f$. It is solved at every stress update by the Newton iteration of
 * Marmot::Materials::ThicknessMicroproblem::solve() for the values of @f$w@f$ at the integration points through the
 * thickness, with the acoustic tensor @f$Q_q = A_q[N, N]@f$ of the first Piola-Kirchhoff tangent @f$A_q@f$ of the
 * point. The terms of @f$u^1@f$ with the tangential gradient, which are of higher order in @f$h@f$, are neglected.
 *
 * The force and the surface stress of the interface are the ones of MarmotInterfaceMaterialFiniteStrain for the
 * through-thickness mean @f$\langle P \rangle@f$ of the stress. The tangent is the one of the mean, condensed over
 * the corrections of the points,
 * @f[
 *   \bar{A} = \sum_q \omega_q A_q \left( I + \frac{\partial ( w_q \otimes N )}{\partial F^0} \right) .
 * @f]
 * For a single homogeneous layer the response is identical to the one of MarmotInterfaceMaterialFiniteStrain, and
 * for two layers of equal thickness with one integration point each, to the one of
 * MarmotEquilibratedInterfaceMaterialFiniteStrain with equal surface gradients of both faces.
 *
 * The material properties are the ones of Marmot::Materials::ThicknessProfile: @f$[h, \mathrm{mode}, n_G, K, m,
 * (p_k, \text{properties}) \times K]@f$, where the properties are the ones of the bulk material.
 *
 * The state variables are, for every integration point through the thickness, the correction @f$w@f$ of its
 * last converged state, which starts the next iteration, followed by the state variables of its bulk material.
 */
class MarmotGradedInterfaceMaterialFiniteStrain {

public:
  using Tensor3d    = MarmotInterfaceMaterialFiniteStrain::Tensor3d;    ///< first order tensor in 3D
  using Tensor33d   = MarmotInterfaceMaterialFiniteStrain::Tensor33d;   ///< second order tensor in 3D
  using Tensor333d  = MarmotInterfaceMaterialFiniteStrain::Tensor333d;  ///< third order tensor in 3D
  using Tensor3333d = MarmotInterfaceMaterialFiniteStrain::Tensor3333d; ///< fourth order tensor in 3D
  using Tensor6d    = MarmotInterfaceMaterialFiniteStrain::Tensor6d;    ///< vector of six components, one per side
  using Tensor18d   = MarmotInterfaceMaterialFiniteStrain::Tensor18d;   ///< vector of 18 components, one per side

  using State         = MarmotInterfaceMaterialFiniteStrain::State;         ///< interface response, as in the base
  using Tangents      = MarmotInterfaceMaterialFiniteStrain::Tangents;      ///< tangents, as in the base
  using Deformation   = MarmotInterfaceMaterialFiniteStrain::Deformation;   ///< interface kinematics, as in the base
  using TimeIncrement = MarmotInterfaceMaterialFiniteStrain::TimeIncrement; ///< time data, as in the base

  /// number of the material
  const int materialNumber;

  /// Layout of the state variables
  MarmotStateLayoutDynamic stateLayout;

  /**
   * @brief Construct the graded interface material from a registered finite-strain bulk material.
   *
   * @param[in] materialName name of the bulk material, registered in the MarmotMaterialFiniteStrainFactory
   * @param[in] matProperties material properties, as described for this class
   * @param[in] nMaterialProperties number of material properties
   * @param[in] materialNumber_ number of the material
   *
   * @throws std::invalid_argument if the properties are inconsistent or the bulk material is not registered
   */
  MarmotGradedInterfaceMaterialFiniteStrain( const std::string& materialName,
                                             const double*      matProperties,
                                             int                nMaterialProperties,
                                             int                materialNumber_ );

  /// Default destructor
  virtual ~MarmotGradedInterfaceMaterialFiniteStrain() = default;

  /**
   * @brief Set the characteristic element length at the considered quadrature point.
   *
   * @param[in] length characteristic length, handed to all bulk materials
   */
  void setCharacteristicElementLength( double length );

  /**
   * @brief Solve the microproblem of the layer and compute the interface response.
   *
   * @param[in,out] state interface response; the state variables are updated
   * @param[out] tangents tangents of the interface
   * @param[in] deformation total displacements and surface gradients of both sides, the normal and the connector
   * @param[in] timeIncrement time at the beginning of the increment and time increment
   *
   * @throws Marmot::StressUpdateFailed if the iteration for the microproblem does not converge, or if the acoustic
   * tensor of an integration point is singular
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

private:
  /// Integration point through the thickness.
  struct ThicknessPoint {
    double                                        weight;         ///< weight, the weights sum up to one
    double                                        position;       ///< coordinate y in [0, 1]
    std::vector< double >                         properties;     ///< properties of the bulk material
    std::unique_ptr< MarmotMaterialFiniteStrain > material;       ///< bulk material at the position
    std::size_t                                   offset;         ///< offset of the state of this point
    std::size_t                                   nBaseStateVars; ///< number of state variables of the bulk material
  };

  const double                  h;      ///< thickness of the layer
  std::vector< ThicknessPoint > points; ///< integration points through the thickness, from the bottom to the top
};
