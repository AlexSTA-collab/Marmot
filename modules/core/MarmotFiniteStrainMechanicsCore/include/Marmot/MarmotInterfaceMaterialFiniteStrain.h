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

#include <memory>
#include <string>
#include <vector>

/**
 * Finite-strain interface material.
 *
 * The interface is a thin layer of the thickness h between two faces. The layer is described by a bulk material of
 * finite strain, which is evaluated for the deformation gradient
 *
 *   F = I + <grad_s u> + ( [u] - <grad_s u> dTau ) (x) N / ell,
 *
 * where [u] is the jump of the total displacement, <grad_s u> the average of the total surface gradients of the two
 * faces, N the interface normal, ell the separation of the faces along the normal, and dTau the tangential part of
 * the connector between the faces. The connector dX = X+ - X- = ell N + dTau is taken from the mesh, so that a
 * rigid rotation of an inclined interface does not cause stresses.
 *
 * The material returns the force t = (h / ell) P N, which is conjugate to the jump, and the surface stress
 * Sigma = h P T - (h / ell) (P N) (x) dTau, which is conjugate to the average surface gradient, where P is the first
 * Piola-Kirchhoff stress of the bulk material and T = I - N (x) N the tangential projection. Their derivatives are
 * returned as tangents.
 *
 * In contrast to the small-strain interface materials, which receive increments, this material receives the TOTAL
 * displacements and surface gradients. The state variables are the ones of the bulk material.
 */
class MarmotInterfaceMaterialFiniteStrain {

protected:
  const double*         materialProperties;     ///< material properties [h, properties of the base material...]
  const int             nMaterialProperties;    ///< number of material properties
  const double          h;                      ///< thickness of the interface, the first material property
  std::vector< double > baseMaterialProperties; ///< properties handed to the base material
  std::unique_ptr< MarmotMaterialFiniteStrain > baseMaterial; ///< finite-strain material of the bulk

public:
  using Tensor3d    = Marmot::FastorStandardTensors::Tensor3d;    ///< first order tensor in 3D
  using Tensor33d   = Marmot::FastorStandardTensors::Tensor33d;   ///< second order tensor in 3D
  using Tensor333d  = Marmot::FastorStandardTensors::Tensor333d;  ///< third order tensor in 3D
  using Tensor3333d = Marmot::FastorStandardTensors::Tensor3333d; ///< fourth order tensor in 3D
  using Tensor6d    = Marmot::FastorStandardTensors::Tensor6d;    ///< vector of six components, one block per side
  using Tensor18d   = Marmot::FastorStandardTensors::Tensor18d;   ///< vector of 18 components, one block per side

  /// number of the material
  const int materialNumber;

  /**
   * @brief Construct an interface material from a registered finite-strain material.
   *
   * The material properties are [h, properties of the base material...], where h is the thickness of the
   * interface.
   *
   * @param[in] materialName name of the base material, registered in the MarmotMaterialFiniteStrainFactory
   * @param[in] matProperties_ material properties
   * @param[in] nMaterialProperties_ number of material properties, at least two
   * @param[in] materialNumber_ number of the material
   *
   * @throws std::invalid_argument if no property is given, if the thickness is not positive, or if the base material
   * is not registered
   */
  MarmotInterfaceMaterialFiniteStrain( const std::string& materialName,
                                       const double*      matProperties_,
                                       int                nMaterialProperties_,
                                       int                materialNumber_ );

  /// Default destructor
  virtual ~MarmotInterfaceMaterialFiniteStrain() = default;

  /// Layout of the state variables
  MarmotStateLayoutDynamic stateLayout;

  /**
   * Set the characteristic element length at the considered quadrature point.
   * It is needed for the regularization of materials with softening behavior
   * based on the mesh-adjusted softening modulus, and is handed on to the base material.
   *
   * @param[in] length characteristic length
   */
  void setCharacteristicElementLength( double length );

  /// Interface response: the updated values on exit.
  struct State {
    Tensor3d  force;         ///< force, conjugate to the jump of the displacements
    Tensor33d surfaceStress; ///< surface stress, conjugate to the average surface gradient
    double*   stateVars;     ///< pointer to the state variables
  };

  /// Algorithmic tangents, set by computeStress().
  struct Tangents {
    Tensor33d   Q_ik;   ///< derivative of the force with respect to the jump
    Tensor333d  H_ikR;  ///< derivative of the force with respect to the average surface gradient
    Tensor333d  H_iRk;  ///< derivative of the surface stress with respect to the jump
    Tensor3333d A_iRkP; ///< derivative of the surface stress with respect to the average surface gradient
  };

  /// Kinematics of the interface at a quadrature point.
  struct Deformation {
    Tensor6d  U;               ///< total displacements of the two interface sides, top first
    Tensor18d surfaceGradient; ///< total surface gradients of the two interface sides, top first
    Tensor3d  normal;          ///< interface normal, pointing from the first to the second side
    Tensor3d  connector;       ///< reference connector between the paired faces, X+ - X-
  };

  /// Time at the beginning of the increment and time increment.
  struct TimeIncrement {
    double timeOld; ///< time at the beginning of the increment
    double dT;      ///< time increment
  };

  /**
   * @brief Compute the force, the surface stress and their tangents from the total kinematics of the interface.
   *
   * @param[in,out] state the response, and the state variables, which are updated
   * @param[out] tangents the algorithmic tangents
   * @param[in] deformation the kinematics of the interface
   * @param[in] timeIncrement the time at the beginning of the increment and the time increment
   *
   * @throws std::invalid_argument if the faces are not separated along the normal
   */
  virtual void computeStress( State&               state,
                              Tangents&            tangents,
                              const Deformation&   deformation,
                              const TimeIncrement& timeIncrement );

  /**
   * @brief Get a view to the state variables.
   *
   * @param stateName Name of the state variable.
   * @param stateVars Pointer to the state variable array.
   * @return StateView to access the requested state variable.
   */
  StateView getStateView( const std::string& stateName, double* stateVars ) const
  {
    return stateLayout.getStateView( stateVars, stateName );
  }

  /**
   * @brief Get the total number of required state variables.
   *
   * @return Total number of required state variables.
   */
  virtual int getNumberOfRequiredStateVars() const { return stateLayout.totalSize(); }

  /**
   * @brief Initialize the state variables at a material point.
   *
   * @param[in,out] stateVars the state variables
   * @param[in] nStateVars number of state variables
   */
  virtual void initializeYourself( double* stateVars, int nStateVars );

  /**
   * @brief Get the density of the material.
   *
   * @return the density of the base material
   */
  virtual double getDensity() const;
};
