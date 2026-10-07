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
#include "Marmot/MarmotStateHelpers.h"

#include <Eigen/Dense>

#include <functional>
#include <string>
#include <vector>

/**
 * Common core of the equilibrated two-sided interface materials at finite strain.
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
 *
 * The class holds what MarmotEquilibratedInterfaceMaterialFiniteStrain and
 * MarmotEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain have in common: the parsing of the material card,
 * the layout of the state variables, the geometry of the connector, and the local solution for z with the static
 * condensation. The derived materials supply the bulk materials of the two faces through a function that evaluates
 * one face at a total displacement gradient.
 */
class MarmotEquilibratedInterfaceMaterialFiniteStrainBase {

protected:
  const double*         materialProperties;     ///< material properties [p0, p1, h, properties of the bulk law...]
  const int             nMaterialProperties;    ///< number of material properties
  const double          h;                      ///< thickness of the layer, the third material property
  std::vector< double > baseMaterialProperties; ///< properties handed to the bulk materials

public:
  using Tensor3d    = Marmot::FastorStandardTensors::Tensor3d;    ///< first order tensor in 3D
  using Tensor33d   = Marmot::FastorStandardTensors::Tensor33d;   ///< second order tensor in 3D
  using Tensor333d  = Marmot::FastorStandardTensors::Tensor333d;  ///< third order tensor in 3D
  using Tensor3333d = Marmot::FastorStandardTensors::Tensor3333d; ///< fourth order tensor in 3D
  using Tensor6d    = Marmot::FastorStandardTensors::Tensor6d;    ///< vector of six components, one block per side
  using Tensor18d   = Marmot::FastorStandardTensors::Tensor18d;   ///< vector of 18 components, one block per side

  using Vector3d         = Eigen::Matrix< double, 3, 1 >;                    ///< vector of three components
  using Vector9d         = Eigen::Matrix< double, 9, 1 >;                    ///< vector of nine components
  using Vector21d        = Eigen::Matrix< double, 21, 1 >;                   ///< generalized strain x = (w, A+, A-)
  using Matrix3dRowMajor = Eigen::Matrix< double, 3, 3, Eigen::RowMajor >;   ///< row-major 3x3 matrix
  using Matrix9dRowMajor = Eigen::Matrix< double, 9, 9, Eigen::RowMajor >;   ///< row-major 9x9 matrix
  using Matrix21d        = Eigen::Matrix< double, 21, 21, Eigen::RowMajor >; ///< tangent over (w, A+, A-)

  /// number of the material
  const int materialNumber;

  /// Layout of the state variables
  MarmotStateLayoutDynamic stateLayout;

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
   * @brief Thickness of the layer normal to the interface.
   *
   * @return the thickness h of the material card
   */
  double thickness() const { return h; }

  /// Default destructor
  virtual ~MarmotEquilibratedInterfaceMaterialFiniteStrainBase() = default;

protected:
  /// The two faces of the layer.
  enum class Side {
    Top,   ///< the face on the side of the normal, carrying A+
    Bottom ///< the face against the normal, carrying A-
  };

  /// Trial response of one face at a total displacement gradient.
  struct SideTrial {
    Matrix3dRowMajor P; ///< first Piola-Kirchhoff stress
    Matrix9dRowMajor A; ///< derivative of the first Piola-Kirchhoff stress with respect to F
    Matrix3dRowMajor dP_dNonlocal   = Matrix3dRowMajor::Zero(); ///< derivative of P with respect to the nonlocal field
    Vector9d         dL_dF          = Vector9d::Zero();         ///< derivative of the local field with respect to F
    double           dL_dN          = 0.0; ///< derivative of the local field with respect to the nonlocal field
    double           L              = 0.0; ///< local driving field
    double           nonlocalRadius = 0.0; ///< nonlocal radius
    std::vector< double > trialStateVars;  ///< state variables after the trial
  };

  /// Result of the local solution and the condensation.
  struct Solution {
    Vector3d  z;                                     ///< normal-gradient jump
    SideTrial plus;                                  ///< response of the top face
    SideTrial minus;                                 ///< response of the bottom face
    Vector21d p;                                     ///< generalized stress (f, S+, S-)
    Matrix21d tangent;                               ///< condensed tangent over (w, A+, A-)
    Vector21d dStress_dNonlocal = Vector21d::Zero(); ///< derivative of p with respect to the nonlocal field
    Eigen::Matrix< double, 1, 21 >
      dLocalField_dStrain = Eigen::Matrix< double, 1, 21 >::Zero(); ///< derivative of Lbar with respect to (w, A+, A-)
    double dLocalField_dNonlocal = 0.0; ///< derivative of Lbar with respect to the nonlocal field
    double localField            = 0.0; ///< average Lbar of the local driving field of the faces
    double nonlocalRadius        = 0.0; ///< average of the nonlocal radii of the faces
  };

  /**
   * @brief Construct the part of the material that is independent of the bulk law.
   *
   * @param[in] matProperties_ material properties [p0, p1, h, ...]
   * @param[in] nMaterialProperties_ number of material properties, at least three
   * @param[in] materialNumber_ number of the material
   * @param[in] className name of the derived class, for the error messages
   *
   * @throws std::invalid_argument if fewer than three properties are given or the thickness is not positive
   */
  MarmotEquilibratedInterfaceMaterialFiniteStrainBase( const double* matProperties_,
                                                       int           nMaterialProperties_,
                                                       int           materialNumber_,
                                                       const char*   className );

  /**
   * @brief Finish the layout of the state variables, once the bulk materials are known.
   *
   * @param[in] nTopStateVars number of state variables of the bulk material of the top face
   * @param[in] nBottomStateVars number of state variables of the bulk material of the bottom face
   */
  void finalizeStateLayout( int nTopStateVars, int nBottomStateVars );

  /**
   * @brief State variables of the bulk material of one face.
   *
   * @param[in] side the face
   * @param[in] stateVars state variables of the interface material
   * @return pointer to the state variables of the bulk material
   */
  double* bulkStateVars( Side side, double* stateVars ) const
  {
    return stateLayout.getPtr( stateVars, side == Side::Top ? "topMaterialStateVars" : "bottomMaterialStateVars" );
  }

  /**
   * @brief State variables of the bulk material of one face, for reading.
   *
   * @param[in] side the face
   * @param[in] stateVars state variables of the interface material
   * @return pointer to the state variables of the bulk material
   */
  const double* bulkStateVars( Side side, const double* stateVars ) const
  {
    return bulkStateVars( side, const_cast< double* >( stateVars ) );
  }

  /**
   * @brief First Piola-Kirchhoff stress and its tangent from the Kirchhoff stress of a bulk material.
   *
   * @param[in] F deformation gradient
   * @param[in] tau Kirchhoff stress
   * @param[in] dTau_dF derivative of the Kirchhoff stress with respect to F
   * @return a trial with P and A = dP/dF, the other members are zero
   */
  static SideTrial firstPiolaFromKirchhoff( const Tensor33d& F, const Tensor33d& tau, const Tensor3333d& dTau_dF );

  /**
   * @brief Solve for the normal-gradient jump and condense it.
   *
   * The local Newton iteration for z enforces P+ N = P- N, and every trial restarts from the committed state. The
   * converged z, the first Piola-Kirchhoff stresses and the state variables of the bulk materials are committed to
   * @p stateVars.
   *
   * @param[in,out] stateVars state variables of the interface material
   * @param[in] deformation kinematics of the interface
   * @param[in] evaluateSide function that evaluates the bulk material of a face at a total displacement gradient,
   * from a scratch copy of the committed state
   * @param[in] withNonlocal whether the derivatives with respect to the nonlocal field are condensed too
   * @return the solution
   *
   * @throws std::invalid_argument if the normal vanishes or the faces are not separated along the normal
   */
  Solution solve( double*                                                            stateVars,
                  const Deformation&                                                 deformation,
                  const std::function< SideTrial( Side, const Matrix3dRowMajor& ) >& evaluateSide,
                  bool                                                               withNonlocal ) const;

  /**
   * @brief Write the generalized stress and the condensed tangent of a solution to the response of the material.
   *
   * @param[in] solution the solution
   * @param[out] state the response
   * @param[out] tangents the condensed tangent
   */
  static void writeResponse( const Solution& solution, State& state, Tangents& tangents );
};
