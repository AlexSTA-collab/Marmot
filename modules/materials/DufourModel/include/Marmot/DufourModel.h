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
 * Alexander Dummer alexander.dummer@uibk.ac.at
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
#include "Marmot/MarmotFiniteStrainViscoelasticity.h"
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrain.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <string>
#include <tuple>
#include <vector>

namespace Marmot::Materials {

  /**
   * @class DufourModel
   * @brief Viscoelastic-viscoplastic, nonlocal damage model for adhesive layers at finite strain.
   *
   * The model combines a Pence-Gou hyperelastic law, an optional generalized Maxwell chain for the relaxation of the
   * elastic stress, an elasto-viscoplastic return map with a pressure-sensitive yield function and a non-associative
   * flow rule, and a damage variable that is driven by the plastic strain and by a nonlocal field. The Kirchhoff
   * stress is @f$ \boldsymbol{\tau} = (1 - \omega)\,\boldsymbol{\tau}_\mathrm{eff} @f$.
   *
   * @par Damage
   * The damage variable @f$ \omega = 1 - (1 - \omega_s)(1 - \omega_f) @f$ has two parts. The softening part
   * @f$ \omega_s = 1 - \exp(-\bar\alpha_P / \varepsilon_F) @f$ is always active. It depends on the weighted
   * plastic strain @f$ \bar\alpha_P = m\,\bar N + (1 - m)\,\alpha_P @f$, where @f$ \bar N @f$ is the nonlocal field and
   * @f$ m @f$ the weight of the nonlocal part. The failure part @f$ \omega_f @f$ is zero until a driver @f$ D @f$
   * reaches one, and then grows as @f$ \omega_f = 1 - \exp(-(D - 1)/d_F) @f$. The driver is a stress-weighted
   * accumulation of the plastic strain,
   * @f[
   *   D = c_{SW} \int \left[ \mathrm{e}^{E(T)} - \tfrac{1}{b_{SW}}\,\mathrm{e}^{-E(T)} \right] w(\dot\kappa)\,
   *       \mathrm{d}\bar\alpha_P ,
   * @f]
   * after Rice and Tracey and Smith et al., see Abrari Vajari et al., CMAME 400 (2022) 115467, eq. (65)-(66), with
   * the stress triaxiality @f$ T @f$ and without a Lode term. For @f$ c_{SW} = 0 @f$ the failure part is absent and
   * the model is the original softening law. The two-variable structure follows Nguyen et al., Int. J. Solids
   * Struct. 96 (2016).
   *
   * @par Material properties
   * | index | symbol | meaning |
   * |-------|--------|---------|
   * | 0, 1 | @f$ K, G @f$ | bulk and shear modulus |
   * | 2, 3 | @f$ f_t, f_c @f$ | tensile and compressive strength |
   * | 4..10 | @f$ Q_1, Q_2, b_1 \dots b_5 @f$ | hardening law of the yield strength |
   * | 11, 12 | @f$ \eta_{VP}, n @f$ | viscosity and rate exponent |
   * | 13, 14 | @f$ \nu_P^+, \nu_P^- @f$ | plastic Poisson ratios in tension and compression |
   * | 15, 16 | @f$ \varepsilon_F, \omega_\mathrm{max} @f$ | softening strain and largest damage |
   * | 17, 18 | @f$ l_d, m @f$ | nonlocal radius and weight of the nonlocal part |
   * | 19 | @f$ \rho @f$ | density (optional, zero if absent) |
   * | 20 | @f$ c_{SW} @f$ | scale of the driver, zero switches the failure part off |
   * | 21 | @f$ b_{SW} @f$ | divisor of the mirrored compression term, @f$ 10^6 @f$ switches it off |
   * | 22 | @f$ d_F @f$ | rate of the failure part after initiation; it does not move the point of initiation |
   * | 23 | @f$ n_M @f$ | number of Maxwell branches, an integer from 0 to nMaxwellMax |
   * | 24 + 3i | @f$ \gamma_{G,i} @f$ | relative deviatoric modulus of branch i |
   * | 25 + 3i | @f$ \gamma_{K,i} @f$ | relative volumetric modulus of branch i |
   * | 26 + 3i | @f$ \tau_i @f$ | relaxation time of branch i |
   * | E + 0..2 | @f$ b_0, b_1, b_2 @f$ | shape of the exponent, mandatory for @f$ c_{SW} \neq 0 @f$ |
   * | E + 3 | @f$ \dot\kappa_\mathrm{ref} @f$ | reference rate of the driver, zero switches the rate term off |
   * | E + 4 | @f$ s @f$ | rate exponent; exactly -1 means that the exponent n is used |
   * | E + 5 | @f$ T_\mathrm{cap} @f$ | cap on the triaxiality before @f$ E @f$ is evaluated, zero for no cap |
   * | E + 6..8 | @f$ A, T_c, w @f$ | amplitude, centre and width of an optional local bump on g |
   *
   * with @f$ E = 24 + 3 n_M @f$. The card has no entries beyond @f$ E + 8 @f$. The card is positional, and a card of
   * an older layout is refused when it is recognized, because a shifted entry would otherwise be read silently as
   * something else.
   *
   * @par State variables
   * The plastic deformation gradient Fp, the plastic strain alphaP, the damage omega, the driver damageDriver, the
   * weighted plastic strain alphaPBar of the previous increment, the reference stress PK2Ref and the states of the
   * Maxwell branches veDev and veVol, and the output variables triax, omegaS, omegaF and lode.
   *
   * The nonlocal radius is the value of the card and constant.
   *
   * @par Energy and dissipation
   * The elastic energy density of the response is @f$ (1 - \omega)\,\psi_\mathrm{eff} @f$ with the potential of the
   * elastic part. The dissipation is cumulative: the incoming value is incremented by
   * @f$ (1 - \omega)\,\boldsymbol{M} : \Delta\boldsymbol{\varepsilon}^p + \psi_\mathrm{eff}\,\Delta\omega @f$, with the
   * effective Mandel stress @f$ \boldsymbol{M} @f$, the plastic increment @f$ \Delta\boldsymbol{\varepsilon}^p @f$ of
   * the return map and the energy that the increase of the damage releases. A decrease of the damage, which the
   * nonlocal field can cause, does not dissipate. The viscous dissipation of the Maxwell branches is not included,
   * because the energy stored in the branches is not tracked.
   */
  class DufourModel : public MarmotMaterialGradientEnhancedFiniteStrain {

    using Tensor33d   = Marmot::FastorStandardTensors::Tensor33d;
    using Tensor3333d = Marmot::FastorStandardTensors::Tensor3333d;

    /// Largest number of Maxwell branches, which fixes the size of the state variables.
    inline const static int nMaxwellMax = 7;

    /// Index of the number of Maxwell branches on the card.
    inline const static int idxNMaxwell = 23;

    /// Index of the first Maxwell triplet on the card.
    inline const static int idxPronyBase = 24;

    /// Number of entries of the card after the Maxwell triplets.
    inline const static int nSwdfmExtra = 9;

  public:
    /// Sentinel value of the rate exponent that stands for the exponent n of the rate law.
    inline const static double swdfmSTieToN = -1.0;

    /// Rice-Tracey / Smith exponent of the driver for negative triaxialities, a micromechanical constant.
    inline const static double swdfmExponent = 1.3;

    /// Lower bound of the triaxiality used by the driver.
    inline const static double etaMin = -0.5;

    /// Upper bound of the triaxiality used by the driver, which keeps the weight bounded.
    inline const static double etaMax = 1.2;

    /// Lower bound of the rate of the driver, which keeps the logarithm of the rate weight finite.
    inline const static double swdfmKdotMin = 1e-8;

    /**
     * @brief Largest increase of the driver accepted in one global increment.
     *
     * A larger increase makes the material throw Marmot::StressUpdateFailed, so that the solver retakes the
     * increment with a smaller time step. This is a numerical tolerance, not a property of the material: it changes
     * how finely the crossing of the threshold D = 1 is resolved, never how much damage accumulates.
     */
    inline const static double dDMaxPerIncrement = 0.10;

    // elastic constants
    const double K; ///< bulk modulus
    const double G; ///< shear modulus

    // plasticity
    const double ft; ///< tensile strength
    const double fc; ///< compressive strength
    const double Q1; ///< hardening parameter
    const double Q2; ///< hardening parameter
    const double b1; ///< hardening parameter
    const double b2; ///< hardening parameter
    const double b3; ///< hardening parameter
    const double b4; ///< hardening parameter
    const double b5; ///< hardening parameter

    // rate dependence
    const double eta_VP; ///< viscosity
    const double n;      ///< rate exponent

    // viscoplastic flow
    const double nuP_plus;  ///< plastic Poisson ratio in tension
    const double nuP_minus; ///< plastic Poisson ratio in compression

    // damage
    const double epsF;     ///< softening strain
    const double omegaMax; ///< largest damage
    const double ld;       ///< nonlocal radius
    const double m;        ///< weight of the nonlocal field in the weighted plastic strain

    const double density; ///< density

    // driver of the failure part of the damage
    const double cSW; ///< scale of the driver
    const double bSW; ///< divisor of the mirrored compression term
    const double dF;  ///< rate of the failure part after initiation

    const double swdfmB0;      ///< shape of the exponent of the driver: E'(T) = (b0 + b1 T + b2 T^2)^2
    const double swdfmB1;      ///< shape of the exponent of the driver
    const double swdfmB2;      ///< shape of the exponent of the driver
    const double swdfmKdotRef; ///< reference rate of the driver
    const double swdfmS;       ///< rate exponent of the driver
    const double swdfmTCap;    ///< cap on the triaxiality before the exponent is evaluated
    const double swdfmBumpA;   ///< amplitude of the local bump on the weight g
    const double swdfmBumpTc;  ///< centre of the local bump on the weight g
    const double swdfmBumpW;   ///< width of the local bump on the weight g

    const int nMaxwell; ///< number of Maxwell branches

    /// properties of the deviatoric Maxwell branches
    const ContinuumMechanics::FiniteStrain::Viscoelasticity::MaxwellProperties maxwellDev;

    /// properties of the volumetric Maxwell branches
    const ContinuumMechanics::FiniteStrain::Viscoelasticity::MaxwellProperties maxwellVol;

    /**
     * @brief Construct the material from the material card.
     *
     * @param[in] materialProperties material properties, see the class description
     * @param[in] nMaterialProperties number of material properties
     * @param[in] materialNumber number of the material
     *
     * @throws std::invalid_argument if a required property is missing, if the card is of an older layout, if the
     * sum of the relative moduli of the Maxwell branches is not below one, if the driver is active without a shape of
     * the exponent, or if a rate exponent is given without a reference rate
     */
    DufourModel( const double* materialProperties, int nMaterialProperties, int materialNumber );

    /**
     * @brief Compute the Kirchhoff stress and the tangents for a deformation gradient and a nonlocal field.
     *
     * @param[in,out] response the response; the state variables are updated
     * @param[out] tangents the derivatives of the stress and of the local field
     * @param[in] deformation the deformation gradient and the nonlocal field
     * @param[in] timeIncrement the time at the beginning of the increment and the time increment
     *
     * @throws std::runtime_error if the return map does not converge
     * @throws Marmot::StressUpdateFailed if the driver of the damage grows by more than dDMaxPerIncrement
     */
    void computeStress( ConstitutiveResponse< 3 >& response,
                        AlgorithmicModuli< 3 >&    tangents,
                        const Deformation< 3 >&    deformation,
                        const TimeIncrement&       timeIncrement ) const override;

    /// @copydoc MarmotMaterialGradientEnhancedFiniteStrain::initializeYourself
    void initializeYourself( double* stateVars, int nStateVars ) override;

    /// @copydoc MarmotMaterialGradientEnhancedFiniteStrain::getDensity
    double getDensity( const double* stateVars ) const override;

    /**
     * @brief Failure part of the damage as a function of the driver.
     *
     * @f$ \omega_f = 1 - \exp(-(D - 1)/d_F) @f$ for @f$ D > 1 @f$ and zero otherwise.
     *
     * @param[in] D the driver
     * @return the failure part
     */
    double omegaFOfDriver( const double D ) const { return D > 1.0 ? 1.0 - std::exp( -( D - 1.0 ) / dF ) : 0.0; }

    /**
     * @brief Derivative of the failure part of the damage with respect to the driver.
     *
     * @param[in] D the driver
     * @return the derivative
     */
    double dOmegaFOfDriver( const double D ) const { return D > 1.0 ? std::exp( -( D - 1.0 ) / dF ) / dF : 0.0; }

    /**
     * @brief Softening part of the damage and its derivative with respect to the weighted plastic strain.
     *
     * One function for the value and the derivative, so that they cannot differ. At zero strain the derivative is
     * 1 / epsF, which fixes the initial softening slope.
     *
     * @param[in] kbar weighted plastic strain
     * @param[out] omega_s softening part of the damage, @f$ 1 - \exp(-\bar k / \varepsilon_F) @f$
     * @param[out] dOmega_s derivative of the softening part with respect to the weighted plastic strain
     */
    void softening( const double kbar, double& omega_s, double& dOmega_s ) const
    {
      const double x = kbar / epsF;
      const double E = std::exp( -x );
      omega_s        = 1.0 - E;
      dOmega_s       = E / epsF;
    }

    /**
     * @brief Whether the card gives a shape of the exponent of the driver.
     *
     * @return true unless b0, b1 and b2 are all zero
     */
    bool swdfmShapeGiven() const { return swdfmB0 != 0.0 || swdfmB1 != 0.0 || swdfmB2 != 0.0; }

    /**
     * @brief Logarithm of the weight g of the driver, with the factor of the rate.
     *
     * One function for the exponent, so that the weight, its tests and its documentation cannot differ. For a
     * negative triaxiality the exponent is linear with the slope b0^2, or swdfmExponent if the card gives no shape;
     * for a positive one it is the integral of (b0 + b1 T + b2 T^2)^2, so that it cannot decrease with the
     * triaxiality. The bump and the rate factor w = (kdot / kdotRef)^(-s) are added to the exponent.
     *
     * @param[in] T stress triaxiality, already clamped by the caller
     * @param[in] kdot rate of the weighted plastic strain of this increment, bounded below internally
     * @param[out] kdotDLogGDKdot if not null, the derivative of the result with respect to the logarithm of the rate
     * @return the logarithm of the weight
     */
    double swdfmLogG( const double T, const double kdot = 0.0, double* kdotDLogGDKdot = nullptr ) const;

    /**
     * @brief Stress triaxiality p / q of a stress.
     *
     * The triaxiality is a ratio of invariants, so it is the same for the Kirchhoff, the Cauchy and the effective
     * stress.
     *
     * @param[in] tau the stress
     * @return the triaxiality, zero for a vanishing deviator
     */
    static double triaxiality( const Tensor33d& tau );

    /**
     * @brief Lode angle parameter @f$ \zeta = \cos(3\theta) = 3\sqrt{3} J_3 / (2 J_2^{3/2}) @f$ of a stress.
     *
     * It is +1 for axisymmetric tension, 0 for pure shear and -1 for axisymmetric compression, as in Abrari Vajari et
     * al., CMAME 400 (2022) 115467, eq. (64). It is a ratio of invariants like the triaxiality.
     *
     * @param[in] tau the stress
     * @return the Lode parameter, zero for a vanishing deviator
     */
    static double lodeParameter( const Tensor33d& tau );

  private:
    /// The state of the material point and the time step that the helper functions of a stress update need.
    struct Context {
      double* stateVars; ///< state variables of the material point
      double  dT;        ///< time increment
    };

    /**
     * @brief Build the properties of the deviatoric or the volumetric Maxwell branches from the card.
     *
     * @param[in] materialProperties material properties
     * @param[in] nMaterialProperties number of material properties
     * @param[in] which 0 for the deviatoric weights, 1 for the volumetric weights
     * @return the properties of the branches
     *
     * @throws std::invalid_argument if there are too many branches or if their triplets are incomplete
     */
    static ContinuumMechanics::FiniteStrain::Viscoelasticity::MaxwellProperties makeMaxwellProperties(
      const double* materialProperties,
      int           nMaterialProperties,
      int           which );

    /**
     * @brief Entry of the block of the card after the Maxwell triplets, zero if absent.
     *
     * @param[in] materialProperties material properties
     * @param[in] nMaterialProperties number of material properties
     * @param[in] which index inside the block
     * @return the entry
     */
    static double swdfmExtra( const double* materialProperties, int nMaterialProperties, int which );

    /**
     * @brief Damage variable and its derivatives with respect to the local and the nonlocal plastic strain.
     *
     * Also records the split into the softening and the failure part and the Lode parameter in the state variables.
     *
     * @param[in] ctx state variables and time increment
     * @param[in] alphaP_local local plastic strain
     * @param[in] alphaP_nonlocal nonlocal field
     * @param[in] tau_eff effective stress, which gives the triaxiality
     * @param[in] D_old driver at the beginning of the increment
     * @param[in] alphaPBar_old weighted plastic strain at the beginning of the increment
     * @return omega, its derivatives with respect to the local and the nonlocal plastic strain, the new driver and the
     * new weighted plastic strain
     */
    std::tuple< double, double, double, double, double > computeOmega( const Context&   ctx,
                                                                       const double     alphaP_local,
                                                                       const double     alphaP_nonlocal,
                                                                       const Tensor33d& tau_eff,
                                                                       const double     D_old,
                                                                       const double     alphaPBar_old ) const;

    /**
     * @brief Relax the instantaneous elastic second Piola-Kirchhoff stress through the generalized Maxwell chain.
     *
     * The deviatoric and the volumetric part are relaxed independently with their own relative moduli and common
     * relaxation times. The recurrence of the branches is the one of
     * ContinuumMechanics::FiniteStrain::Viscoelasticity::evaluateGeneralizedMaxwellModel.
     *
     * @param[in] ctx state variables and time increment
     * @param[in] PK2_0 instantaneous second Piola-Kirchhoff stress
     * @param[in] dPK2_0_dCe its derivative with respect to the elastic right Cauchy-Green tensor
     * @param[in] commit false advances the branches on a scratch copy and leaves the history untouched, so that the
     * function may be called as often as the return map needs; true writes the advanced branches and the new
     * reference stress to the state variables, which is to be done once per increment with the converged state
     * @return the relaxed stress and its derivative with respect to the elastic right Cauchy-Green tensor
     */
    std::tuple< Tensor33d, Tensor3333d > applyViscoelasticity( const Context&     ctx,
                                                               const Tensor33d&   PK2_0,
                                                               const Tensor3333d& dPK2_0_dCe,
                                                               const bool         commit ) const;

    /**
     * @brief Mandel stress of the elastic deformation, with the viscoelastic relaxation, and its derivative.
     *
     * @param[in] ctx state variables and time increment
     * @param[in] Fe elastic deformation gradient
     * @return the Mandel stress and its derivative with respect to Fe
     */
    std::tuple< Tensor33d, Tensor3333d > computeMandelStress( const Context& ctx, const Tensor33d& Fe ) const;

    /**
     * @brief Yield strength and its derivative with respect to the plastic strain.
     *
     * @param[in] alphaP plastic strain
     * @return the yield strength and its derivative
     */
    std::tuple< double, double > computeBetaP( const double alphaP ) const;

    /**
     * @brief Yield function, flow potential and the plastic multiplier scale for a Mandel stress.
     *
     * @param[in] mandelStress Mandel stress
     * @param[in] betaP yield strength
     * @param[in] J determinant of the deformation gradient
     * @return the yield function, its derivative with respect to the Mandel stress and to the yield strength, the
     * derivative of the flow potential, its second derivative, the norm h of the flow direction and the derivative of h
     */
    std::tuple< double, Tensor33d, double, Tensor33d, Tensor3333d, double, Tensor33d > yieldFunctionFromStress(
      const Tensor33d& mandelStress,
      const double     betaP,
      const double     J ) const;

    /**
     * @brief Yield function of an elastic deformation and its derivatives with respect to Fe.
     *
     * @param[in] ctx state variables and time increment
     * @param[in] Fe elastic deformation gradient
     * @param[in] betaP yield strength
     * @param[in] J determinant of the deformation gradient
     * @return the yield function, its derivative with respect to Fe and to the yield strength, the derivative of the
     * flow potential and the derivative of h with respect to Fe
     */
    std::tuple< double, Tensor33d, double, Tensor33d, Tensor33d > yieldFunction( const Context&   ctx,
                                                                                 const Tensor33d& Fe,
                                                                                 const double     betaP,
                                                                                 const double     J ) const;

    /**
     * @brief Whether a trial elastic deformation violates the yield function.
     *
     * @param[in] ctx state variables and time increment
     * @param[in] Fe elastic deformation gradient
     * @param[in] betaP yield strength
     * @param[in] J determinant of the deformation gradient
     * @return true if the yield function is positive
     */
    bool isYielding( const Context& ctx, const Tensor33d& Fe, const double betaP, const double J ) const;

    /**
     * @brief Residual and Jacobian of the return map.
     *
     * The unknowns are the elastic deformation gradient, the plastic strain and the plastic multiplier.
     *
     * @param[in] ctx state variables and time increment
     * @param[in] X the unknowns
     * @param[in] FeTrial trial elastic deformation gradient
     * @param[in] alphaPTrial plastic strain at the beginning of the increment
     * @param[in] dt time increment
     * @param[in] J determinant of the deformation gradient
     * @return the residual and its Jacobian
     */
    std::tuple< Eigen::VectorXd, Eigen::MatrixXd > computeResidualVectorAndTangent( const Context&         ctx,
                                                                                    const Eigen::VectorXd& X,
                                                                                    const Tensor33d&       FeTrial,
                                                                                    const double           alphaPTrial,
                                                                                    const double           dt,
                                                                                    const double           J ) const;
  };

} // namespace Marmot::Materials
