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

/**
 * @file MarmotStabPressureXInterfaceMaterialHypoElastic.h
 * @brief Equilibrated two-sided hypoelastic interface material with a mixed, independent pressure.
 */
#pragma once

#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotStateHelpers.h"

#include <memory>
#include <string>
#include <vector>

/**
 * @brief Equilibrated two-sided hypoelastic interface material whose hydrostatic pressure is an independent field.
 *
 * The material has the kinematics of MarmotEquilibratedXInterfaceMaterialHypoElastic: the surface gradients
 * @f$A^\pm@f$ of the two faces drive two embedded bulk materials, and the normal-gradient jump @f$z@f$ is condensed
 * by a local Newton iteration. In the incompressible limit, which a von Mises layer with a vanishing hardening
 * modulus reaches (@f$K/\mu_\mathrm{p} \sim 10^7@f$), the pressure @f$K \operatorname{tr}\varepsilon@f$ of a
 * displacement formulation multiplies a trace that is driven to zero, and amplifies the discretisation error into
 * an element-to-element pressure checkerboard. Here the pressure @f$p@f$ is supplied by the element as a nodal field
 * instead, and the stress used for equilibrium is
 * @f[
 *   \tilde\sigma = \operatorname{dev}\sigma(\varepsilon) - p \, I .
 * @f]
 * The volumetric law is returned as the residual
 * @f[
 *   r_p = \langle \operatorname{tr}\Delta\varepsilon \rangle + \Delta p / K ,
 * @f]
 * where @f$\langle\cdot\rangle@f$ is the mean of the two faces, to be enforced weakly by the element, together with
 * the stabilisation of the pressure gradient, which belongs to the element.
 *
 * The law is incremental, as the hypoelastic base materials are: the strain increment and the pressure increment
 * accumulate from zero, so that enforcing the increment at every step is equivalent to the total relation
 * @f$p = -K \operatorname{tr}\varepsilon@f$.
 *
 * @par Restriction
 * Replacing the volumetric part of the stress is admissible only because the yield surface of the base material
 * is independent of the pressure, as the one of von Mises is: the deviatoric return map is then unaffected by
 * @f$p@f$. For a pressure-dependent law the pressure would have to enter the return map.
 *
 * @par Structure of the tangent
 * @f$\tilde\sigma@f$ is linear in the stress and in @f$p@f$, so @f$\partial\tilde\sigma/\partial G =
 * P_\mathrm{dev} C@f$ and @f$\partial\tilde\sigma/\partial p = -I@f$. Two terms vanish identically, which keeps the
 * coupling to the pressure purely geometric. The traction mismatch @f$R_t(\tilde\sigma^+ - \tilde\sigma^-)@f$ does
 * not depend on @f$p@f$, because the term @f$-pI@f$ cancels, so that @f$\partial z/\partial p = 0@f$. The mean
 * volumetric strain does not depend on @f$z@f$, because the volumetric parts of @f$g^+@f$ and @f$g^-@f$ cancel in
 * the mean, so that @f$\partial r_p/\partial z = 0@f$. The condensation of @f$z@f$ of the block over the generalized
 * strain @f$x = (w, A^+, A^-)@f$ is the one of the equilibrated material, built from @f$P_\mathrm{dev} C@f$.
 *
 * @par Local solve
 * The local solve is the one the results of the interface were produced with: a Newton iteration with a
 * backtracking line search on the deviatoric traction mismatch, which throws if it does not converge. The gauge of
 * MarmotEquilibratedXInterfaceMaterialHypoElastic is not applied.
 */
class MarmotStabPressureXInterfaceMaterialHypoElastic {

protected:
  const double*         materialProperties;     ///< material properties [E, nu, h, properties of the base material...]
  const int             nMaterialProperties;    ///< number of material properties
  double                h            = 0.0;     ///< thickness of the interface layer, the third material property
  double                bulkModulus  = 0.0;     ///< bulk modulus K of the layer, from E and nu
  double                shearModulus = 0.0;     ///< shear modulus of the layer, from E and nu
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
   * @brief Construct the material from a registered hypoelastic material.
   *
   * The material properties are [E, nu, h, properties of the base material...], where h is the thickness of the
   * interface layer. Both faces receive E, nu and the properties after h. E and nu also define the bulk modulus of
   * the volumetric law.
   *
   * @param[in] materialName name of the base material, registered in the MarmotMaterialHypoElasticFactory
   * @param[in] matProperties_ material properties
   * @param[in] nMaterialProperties_ number of material properties, at least three
   * @param[in] materialNumber_ number of the material
   *
   * @throws std::invalid_argument if fewer than three properties are given, if h is not positive, if nu is 0.5, or
   * if the base material is not registered
   */
  MarmotStabPressureXInterfaceMaterialHypoElastic( const std::string& materialName,
                                                   const double*      matProperties_,
                                                   int                nMaterialProperties_,
                                                   int                materialNumber_ );

  /// Default destructor
  virtual ~MarmotStabPressureXInterfaceMaterialHypoElastic() = default;

  /// Layout of the state variables
  MarmotStateLayoutDynamic stateLayout;

  /// Characteristic element length
  double characteristicElementLength = 0.0;

  /**
   * @brief Set the characteristic element length at the considered quadrature point.
   *
   * It is handed to the base materials of both faces.
   *
   * @param[in] length characteristic length, assigned to characteristicElementLength
   */
  void setCharacteristicElementLength( double length );

  /// @return the thickness h of the interface layer
  double getInterfaceThickness() const { return h; }

  /// @return the bulk modulus K of the layer
  double getBulkModulus() const { return bulkModulus; }

  /// @return the shear modulus of the layer
  double getShearModulus() const { return shearModulus; }

  /**
   * @brief Result of a stress update: the generalized stress, computed from @f$\tilde\sigma@f$, and the volumetric
   * residual.
   */
  struct State {
    Tensor3d  generalizedForce;   ///< generalized force conjugate to the raw mesh jump w = [u]
    Tensor33d surfaceStressPlus;  ///< generalized surface resultant conjugate to A+
    Tensor33d surfaceStressMinus; ///< generalized surface resultant conjugate to A-
    double    volumetricResidual; ///< r_p = <tr dEps> + dp / K
    double*   stateVars;          ///< pointer to the state variables
  };

  /**
   * @brief Tangent after static condensation of the internal normal-gradient jump z, over the generalized strain
   * x = (w, A+, A-) and the pressure.
   *
   * The 3x3 block structure over x is the one of MarmotEquilibratedXInterfaceMaterialHypoElastic::Tangents, with
   * the blocks of @f$P_\mathrm{dev} C@f$. The blocks of the pressure are purely geometric.
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
    Tensor3d    Q_wp;   ///< d(f)/d(p), 3x1
    Tensor33d   Q_App;  ///< d(S+)/d(p), 9x1
    Tensor33d   Q_Amp;  ///< d(S-)/d(p), 9x1
    Tensor3d    Q_pw;   ///< d(r_p)/d(w), 1x3
    Tensor33d   Q_pAp;  ///< d(r_p)/d(A+), 1x9
    Tensor33d   Q_pAm;  ///< d(r_p)/d(A-), 1x9
    double      Q_pp;   ///< d(r_p)/d(p) = 1 / K
  };

  /**
   * @brief Increment of the interface kinematics and of the pressure.
   */
  struct Deformation {
    Tensor6d  dU;             ///< displacement increment of the two interface sides
    Tensor18d dSurfaceStrain; ///< surface displacement gradient increment of the two interface sides
    Tensor3d  normal;         ///< interface normal
    /// connector between the two faces, from the bottom to the top; zero for coincident faces, for which the
    /// material uses d = 0, the normal separation h and the tangential separation 0
    Tensor3d separationVector;
    double   dPressure; ///< increment of the independent pressure p
  };

  /// Time at the beginning of the increment and time increment.
  struct TimeIncrement {
    double timeOld; ///< time at the beginning of the increment
    double dT;      ///< time increment
  };

  /**
   * @brief Compute the generalized stress, the volumetric residual and the condensed tangent for an increment.
   *
   * @param[in,out] state generalized stress, volumetric residual and state variables
   * @param[out] tangents condensed tangent
   * @param[in] deformation increment of the interface kinematics and of the pressure
   * @param[in] timeIncrement time and time increment
   *
   * @throws std::invalid_argument if the interface normal is zero, or if the connector between the faces has no
   * positive normal component
   * @throws Marmot::StressUpdateFailed if the local traction-equilibrium solve does not converge
   */
  virtual void computeStress( State&               state,
                              Tangents&            tangents,
                              const Deformation&   deformation,
                              const TimeIncrement& timeIncrement );

  /**
   * @brief Get a view to the state variables.
   *
   * A name that ends in Top or Bottom is handed to the base material of that face.
   *
   * @param[in] stateName name of the state variable
   * @param[in] stateVars pointer to the state variable array
   * @return StateView to access the requested state variable
   */
  StateView getStateView( const std::string& stateName, double* stateVars ) const
  {
    const std::size_t n = stateName.size();
    if ( n > 3 && stateName.compare( n - 3, 3, "Top" ) == 0 && topMaterial )
      return topMaterial->getStateView( stateName.substr( 0, n - 3 ),
                                        stateLayout.getPtr( stateVars, "topMaterialStateVars" ) );
    if ( n > 6 && stateName.compare( n - 6, 6, "Bottom" ) == 0 && bottomMaterial )
      return bottomMaterial->getStateView( stateName.substr( 0, n - 6 ),
                                           stateLayout.getPtr( stateVars, "bottomMaterialStateVars" ) );
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
   * Zeroes the state variables and initializes those of both base materials.
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
