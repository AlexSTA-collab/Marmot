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
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotStateHelpers.h"

#include <limits>
#include <memory>
#include <string>
#include <vector>

/**
 * @brief Equilibrated two-sided hypoelastic interface material.
 *
 * The top and bottom surface gradients @f$A^+@f$ and @f$A^-@f$ are kept separate and drive two independent
 * embedded bulk materials. The normal-gradient jump
 * @f[
 *   z = g^+ - g^-
 * @f]
 * is an internal variable, which is condensed at every stress update by a local Newton iteration that enforces
 * the traction equilibrium
 * @f[
 *   \sigma^+ n = \sigma^- n .
 * @f]
 * The average surface gradient @f$\bar{A} = (A^+ + A^-)/2@f$ together with the geometry of the connector between
 * the two faces provides the shared normal-gradient average
 * @f[
 *   \bar{g} = \frac{1}{\ell} \left( w - \bar{A} \, d_\tau \right), \qquad w = [\![u]\!],
 * @f]
 * and @f$g^\pm = \bar{g} \pm z/2@f$.
 *
 * The Jacobian of the equilibrium residual @f$r(z) = t^+(z) - t^-(z)@f$ is the average one-sided acoustic tensor
 * @f$\langle Q \rangle = \tfrac{1}{2} ( n C^+ n + n C^- n )@f$. It stays invertible where @f$Q^+@f$ or @f$Q^-@f$
 * degenerates individually, for instance at yield, so only @f$\langle Q \rangle@f$ is ever inverted.
 *
 * Where both faces flow, even @f$\langle Q \rangle@f$ loses the stiffness along one direction and the local problem
 * has no unique solution. The material then fixes the undetermined component of @f$z@f$ (the "gauge"), so that
 * the local solve and the condensed tangent use the same two-dimensional subspace.
 *
 * After the local convergence the generalized stress and the fully condensed tangent are returned. The tangent has
 * nonzero cross blocks between @f$A^+@f$ and @f$A^-@f$ that the equilibrium constraint introduces.
 */
class MarmotEquilibratedXInterfaceMaterialHypoElastic {

protected:
  const double*         materialProperties;     ///< material properties [E, nu, h, properties of the base material...]
  const int             nMaterialProperties;    ///< number of material properties
  const double          h;                      ///< thickness of the interface layer, the third material property
  std::vector< double > baseMaterialProperties; ///< properties handed to the base materials of both faces
  std::unique_ptr< MarmotMaterialHypoElastic > topMaterial;    ///< hypoelastic material of the top face
  std::unique_ptr< MarmotMaterialHypoElastic > bottomMaterial; ///< hypoelastic material of the bottom face

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
   * @brief Construct an equilibrated interface material from a registered hypoelastic material.
   *
   * The material properties are [E, nu, h, properties of the base material...], where h is the thickness of the
   * interface layer. Both faces receive E, nu and the properties after h.
   *
   * @param[in] materialName name of the base material, registered in the MarmotMaterialHypoElasticFactory
   * @param[in] matProperties_ material properties
   * @param[in] nMaterialProperties_ number of material properties, at least three
   * @param[in] materialNumber_ number of the material
   *
   * @throws std::invalid_argument if fewer than three properties are given, if h is not positive, or if the base
   * material is not registered
   */
  MarmotEquilibratedXInterfaceMaterialHypoElastic( const std::string& materialName,
                                                   const double*      matProperties_,
                                                   int                nMaterialProperties_,
                                                   int                materialNumber_ );

  /// Default destructor
  virtual ~MarmotEquilibratedXInterfaceMaterialHypoElastic() = default;

  /// Layout of the state variables
  MarmotStateLayoutDynamic stateLayout;

  /// Characteristic element length
  double characteristicElementLength = 0.0;

  /**
   * @brief Set the characteristic element length at the considered quadrature point.
   *
   * It is needed for the regularization of base materials with softening behavior and is handed to the base
   * materials of both faces.
   *
   * @param[in] length characteristic length, assigned to characteristicElementLength
   */
  void setCharacteristicElementLength( double length );

  /**
   * @brief Memory of the local solve that has to survive the iterations of an increment.
   *
   * The state variables are restored at the beginning of every iteration of the global solution procedure, so
   * they cannot carry information from one iteration to the next. The decision whether the gauge is active is
   * such information: if it were taken anew in every iteration, a few points whose decision keeps changing can
   * hold the global Newton iteration in a limit cycle. The caller therefore stores one IncrementMemory per
   * quadrature point outside of the state variables and passes it to computeStress().
   */
  struct IncrementMemory {
    /// time at the beginning of the increment in which the decision was taken, NaN if none was taken yet
    double timeOld = std::numeric_limits< double >::quiet_NaN();

    /// decision taken in the increment
    bool gaugeActive = false;

    /**
     * @brief Freeze the decision of an increment: the first decision taken in an increment is kept for the
     * whole increment, including all its iterations and cutbacks.
     *
     * @param[in] timeOldOfIncrement time at the beginning of the current increment
     * @param[in] proposedDecision decision proposed by the current evaluation
     * @return the decision that holds for the increment
     */
    bool decide( double timeOldOfIncrement, bool proposedDecision )
    {
      if ( timeOldOfIncrement != timeOld ) {
        timeOld     = timeOldOfIncrement;
        gaugeActive = proposedDecision;
      }
      return gaugeActive;
    }
  };

  /**
   * @brief State of the interface at a quadrature point: previous values on entry, updated values on exit.
   */
  struct State {
    /// equilibrated generalized force conjugate to the raw mesh jump w = [u]
    Tensor3d generalizedForce;

    /// generalized surface resultant conjugate to A+
    Tensor33d surfaceStressPlus;

    /// generalized surface resultant conjugate to A-
    Tensor33d surfaceStressMinus;

    /// pointer to the state variables
    double* stateVars;

    /// memory that survives the iterations of an increment, or nullptr to take the decision anew at every call
    IncrementMemory* incrementMemory = nullptr;
  };

  /**
   * @brief Full 3x3 block tangent after static condensation of the internal normal-gradient jump z.
   *
   * The tangent is given over the external generalized strain x = (w, A+, A-) and the generalized stress
   * p = (f, S+, S-). No major symmetry is assumed and, in general, the cross blocks are nonzero:
   * d(S+)/d(A-) != 0 and d(S-)/d(A+) != 0. A block of a vector and a tensor is stored row major as a tensor of
   * the next higher rank, for instance Q_wAp(i, j, k) = d(f_i)/d(A+_jk), which is a 3x9 matrix, and
   * Q_Apw(i, j, k) = d(S+_ij)/d(w_k), which is a 9x3 matrix.
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

  /**
   * @brief Increment of the interface kinematics.
   */
  struct Deformation {
    Tensor6d  dU;             ///< displacement increment of the two interface sides
    Tensor18d dSurfaceStrain; ///< surface displacement gradient increment of the two interface sides
    Tensor3d  normal;         ///< interface normal
    /// connector between the two faces, from the bottom to the top; zero for coincident faces, for which the
    /// material uses d = 0, the normal separation h and the tangential separation 0
    Tensor3d separationVector;
  };

  /// Time at the beginning of the increment and time increment.
  struct TimeIncrement {
    double timeOld; ///< time at the beginning of the increment
    double dT;      ///< time increment
  };

  /**
   * @brief Compute the generalized stress and the condensed tangent for an increment.
   *
   * @param[in,out] state generalized stress and state variables
   * @param[out] tangents condensed tangent
   * @param[in] deformation increment of the interface kinematics
   * @param[in] timeIncrement time and time increment
   *
   * @throws std::invalid_argument if the interface normal is zero, or if the connector between the faces has no
   * positive normal component
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
   * @return StateView to access the requested state variable
   */
  StateView getStateView( const std::string& stateName, double* stateVars ) const
  {
    return stateLayout.getStateView( stateVars, stateName );
  }

  /**
   * @brief Get the total number of required state variables.
   *
   * @return the number of state variables of the base materials of both faces and of the interface itself
   */
  virtual int getNumberOfRequiredStateVars() const { return stateLayout.totalSize(); }

  /**
   * @brief Initialize the state variables at a material point.
   *
   * Initializes the state variables of both base materials and sets the internal normal-gradient jump to zero.
   *
   * @param[in,out] stateVars state variables of the material point
   * @param[in] nStateVars number of state variables
   */
  virtual void initializeYourself( double* stateVars, int nStateVars );

  /**
   * @brief Get the density of the material.
   *
   * @return the density of the base material
   */
  virtual double getDensity();
};
