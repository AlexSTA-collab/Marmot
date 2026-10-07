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
 * Equilibrated two-sided interface material at finite strain.
 *
 * The total-Lagrangian counterpart of MarmotEquilibratedXInterfaceMaterialHypoElastic. The two faces keep their own
 * surface gradients A+ and A-, each driving its own bulk material of finite strain, and the normal-gradient jump
 *
 *   z = g+ - g-
 *
 * is an internal variable, condensed at every stress update by a local Newton iteration that enforces traction
 * equilibrium. At finite strain the condition is the continuity of the nominal traction on the reference normal,
 *
 *   P+ N = P- N ,
 *
 * which balances across the layer in the reference configuration. The shared average normal gradient comes from the
 * connector between the faces,
 *
 *   gbar = ( w - Abar d_tau ) / ell ,   w = [u] ,   Abar = (A+ + A-)/2 ,
 *
 * and g+ = gbar + z/2, g- = gbar - z/2. Each face then carries
 *
 *   F+- = I + A+- + g+- (x) N ,
 *
 * so that F is affine in (w, A+, A-, z), and the operators of the condensation are the same constant matrices as in
 * the small-strain material. Only the constitutive pair changes: (sigma, C) becomes (P, dP/dF).
 *
 * The Jacobian of the local residual r(z) = P+ N - P- N is the average one-sided acoustic tensor
 *
 *   <Q> = 1/2 ( N A+ N + N A- N ) ,   A = dP/dF ,
 *
 * so only the average is ever inverted, never one side alone.
 *
 * The material properties are [ p0, p1, h, ... ], as in the small-strain material: the first two properties of the
 * bulk law, the thickness of the layer normal to the interface, then the rest of the properties of the bulk law.
 *
 * The material does not have the gauge of the small-strain material, which fixes the undetermined component of z once
 * both faces reach a perfectly plastic state and <Q> loses its transverse stiffness. Elastic, viscoelastic and
 * hardening responses keep <Q> well conditioned. For near-perfect plasticity the local residual stalls, and the best
 * iterate is accepted.
 *
 * In contrast to the small-strain materials, this material receives the TOTAL displacements and surface gradients.
 */
class MarmotEquilibratedInterfaceMaterialFiniteStrain {

protected:
  const double*         materialProperties;     ///< material properties [p0, p1, h, properties of the bulk law...]
  const int             nMaterialProperties;    ///< number of material properties
  const double          h;                      ///< thickness of the layer, the third material property
  std::vector< double > baseMaterialProperties; ///< properties handed to the bulk materials
  std::unique_ptr< MarmotMaterialFiniteStrain > topMaterial;    ///< finite-strain material of the top face
  std::unique_ptr< MarmotMaterialFiniteStrain > bottomMaterial; ///< finite-strain material of the bottom face

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
   * @brief Construct the material from a registered finite-strain material.
   *
   * @param[in] materialName name of the bulk material, registered in the MarmotMaterialFiniteStrainFactory
   * @param[in] matProperties_ material properties [p0, p1, h, ...]
   * @param[in] nMaterialProperties_ number of material properties, at least three
   * @param[in] materialNumber_ number of the material
   *
   * @throws std::invalid_argument if fewer than three properties are given, if the thickness is not positive, or if
   * the bulk material is not registered
   */
  MarmotEquilibratedInterfaceMaterialFiniteStrain( const std::string& materialName,
                                                   const double*      matProperties_,
                                                   int                nMaterialProperties_,
                                                   int                materialNumber_ );

  /// Default destructor
  virtual ~MarmotEquilibratedInterfaceMaterialFiniteStrain() = default;

  /// Layout of the state variables
  MarmotStateLayoutDynamic stateLayout;

  /**
   * Set the characteristic element length at the considered quadrature point.
   * It is needed for the regularization of materials with softening behavior, and is handed on to the bulk
   * materials of both faces.
   *
   * @param[in] length characteristic length
   */
  void setCharacteristicElementLength( double length );

  /// Generalized stress conjugate to (w, A+, A-), based on the first Piola-Kirchhoff stress.
  struct State {
    Tensor3d  generalizedForce;   ///< conjugate to the jump w
    Tensor33d surfaceStressPlus;  ///< conjugate to A+
    Tensor33d surfaceStressMinus; ///< conjugate to A-
    double*   stateVars;          ///< pointer to the state variables
  };

  /**
   * @brief Full 3x3 block tangent over x = (w, A+, A-) after static condensation of z.
   *
   * No major symmetry is assumed and, in general, the cross blocks are nonzero. A block of a vector and a tensor is
   * stored row major as a tensor of the next higher rank, for instance Q_wAp(i, j, k) = d(f_i)/d(A+_jk), which is a
   * 3x9 matrix, and Q_Apw(i, j, k) = d(S+_ij)/d(w_k), which is a 9x3 matrix.
   */
  struct Tangents {
    Tensor33d   Q_ww;   ///< d(f)/d(w), 3x3
    Tensor333d  Q_wAp;  ///< d(f)/d(A+), 3x9
    Tensor333d  Q_wAm;  ///< d(f)/d(A-), 3x9
    Tensor333d  Q_Apw;  ///< d(S+)/d(w), 9x3
    Tensor3333d Q_ApAp; ///< d(S+)/d(A+), 9x9
    Tensor3333d Q_ApAm; ///< d(S+)/d(A-), 9x9
    Tensor333d  Q_Amw;  ///< d(S-)/d(w), 9x3
    Tensor3333d Q_AmAp; ///< d(S-)/d(A+), 9x9
    Tensor3333d Q_AmAm; ///< d(S-)/d(A-), 9x9
  };

  /// Kinematics of the interface at a quadrature point.
  struct Deformation {
    Tensor6d  U;               ///< total displacements of the two interface sides, top first
    Tensor18d surfaceGradient; ///< total surface gradients A+ and A- of the two interface sides
    Tensor3d  normal;          ///< interface normal, pointing from the first to the second side
    Tensor3d  connector;       ///< reference connector between the paired faces, X+ - X-
  };

  /// Time at the beginning of the increment and time increment.
  struct TimeIncrement {
    double timeOld; ///< time at the beginning of the increment
    double dT;      ///< time increment
  };

  /**
   * @brief Compute the generalized stress and the condensed tangent from the total kinematics of the interface.
   *
   * @param[in,out] state the response, and the state variables, which are updated
   * @param[out] tangents the condensed tangent
   * @param[in] deformation the kinematics of the interface
   * @param[in] timeIncrement the time at the beginning of the increment and the time increment
   *
   * @throws std::invalid_argument if the normal vanishes or the faces are not separated along the normal
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
   * @return the density of the bulk material
   */
  virtual double getDensity() const;
};
