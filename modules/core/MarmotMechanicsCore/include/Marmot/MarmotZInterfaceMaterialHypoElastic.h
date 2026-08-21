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

#include <memory>
#include <string>
#include <vector>

/**
 * Two-sided hypoelastic interface kernel with the normal-gradient jump
 * carried as an EXTERNAL (element) unknown.
 *
 * Relation to the rest of the family
 * ----------------------------------
 * MarmotEquilibratedXInterfaceMaterialHypoElastic parameterizes the same
 * layer, but treats
 *
 *   z_k = [[ u_{k,r} ]] n_r = g^+_k - g^-_k
 *
 * as an INTERNAL variable, eliminated at every stress update by a local
 * Newton iteration on the traction-continuity condition
 *
 *   F(z) = ( sigma^+(z) - sigma^-(z) ) n = 0 ,     F_{,z} = <Q> = (Q^+ + Q^-)/2 .
 *
 * That elimination requires <Q> to be invertible. Under (near-)perfect
 * plasticity it is not: once both faces flow, <Q> loses its transverse
 * stiffness along m_i = <N_ij> n_j (N = dev sigma / ||dev sigma||), the local
 * problem becomes bistable, and the global Newton stalls in a limit cycle.
 * The gauge/freeze machinery in the equilibrated kernel is a workaround for
 * exactly that, and it is a binary decision taken per quadrature point.
 *
 * This kernel removes the elimination instead of patching it. z is an INPUT,
 * supplied by the element from a globally interpolated nodal field, and the
 * equilibrium condition is returned as the residual conjugate to z rather
 * than solved for. There is no local Newton, no condensation, no gauge, and
 * no binary state decision anywhere in this file.
 *
 * Generalized strain / stress pair
 * --------------------------------
 * External generalized strain, size 24:
 *
 *   x = ( w , A^+ , A^- )   (3 + 9 + 9 = 21)  and  z  (3),
 *
 * with w = [[u]] the raw mesh jump and A^pm the two faces' surface gradients.
 * The connector geometry (normal n, separation vector d, ell = d.n,
 * d_tau = d - ell n) gives the shared average normal gradient
 *
 *   gbar = ( w - Abar d_tau ) / ell ,   Abar = (A^+ + A^-)/2 ,
 *
 * and each face's full gradient
 *
 *   G^pm = A^pm + ( gbar +- z/2 ) (x) n .
 *
 * The interface free energy per unit midsurface area is
 *
 *   phi = (h/2) ( psi^+(G^+) + psi^-(G^-) ) ,
 *
 * so the returned generalized stresses are exactly its gradient,
 *
 *   p_x = d phi / d x   ->  ( generalizedForce , surfaceStressPlus , surfaceStressMinus )
 *   p_z = d phi / d z   =  (h/4) ( sigma^+ - sigma^- ) n  =  (h/4) ( t^+ - t^- ) ,
 *
 * i.e. p_z IS the traction imbalance, up to the fixed factor h/4. The element
 * drives it to zero WEAKLY, as the variational equation belonging to the
 * nodal z field, instead of pointwise inside a local solve.
 *
 * Why that fixes perfect plasticity
 * ---------------------------------
 * The condensed tangent of the equilibrated kernel is the Schur complement
 *
 *   K_cond = K_xx - K_xz K_zz^{-1} K_zx ,     K_zz = (h/4) <Q> ,
 *
 * and forming it needs <Q>^{-1}. Leaving z in the global system needs only
 * the FULL block matrix [[K_xx, K_xz],[K_zx, K_zz]] to be regular, which is a
 * strictly weaker requirement: a direction m that <Q> cannot resist is still
 * controlled through the coupling K_zx as long as the surrounding mesh
 * resists it. Because z is interpolated with the midsurface shape functions,
 * neighbouring stations must ALSO be degenerate along the same m before a
 * genuine null mode appears.
 *
 * For the residual case where they are (a fully developed, uniformly oriented
 * shear band), a variationally consistent regularization is added: an elastic
 * spring of relative stiffness zeta acting on the z INCREMENT of the step,
 *
 *   phi_reg = (h/4) (zeta/2) dz . Q^e dz ,   Q^e_ik = mu delta_ik + (lambda+mu) n_i n_k ,
 *
 * contributing (h/4) zeta Q^e dz to p_z and (h/4) zeta Q^e to K_zz. It is
 * always on, smooth, and needs no degeneracy detector -- deliberately, since
 * a per-point binary criterion is what produced the limit cycle in the first
 * place. Where <Q> is healthy it is an O(zeta) perturbation; where <Q> is
 * degenerate it selects the minimum-increment solution and bounds the drift
 * of z along the undetermined direction (an unbounded drift would keep
 * accumulating spurious plastic strain on both faces).
 *
 * zeta defaults to 1e-6, is overridable through the environment variable
 * MARMOT_ZIFACE_REG, and can be set per instance with
 * setGradientJumpRegularization().
 *
 * Consistency with the equilibrated kernel
 * ----------------------------------------
 * With zeta = 0 and z solved to p_z = 0, this kernel reproduces
 * MarmotEquilibratedXInterfaceMaterialHypoElastic exactly: same p_x, and
 * K_xx - K_xz K_zz^{-1} K_zx equals its condensed tangent. That equivalence
 * is asserted in TestMarmotZInterfaceMaterialHypoElastic.
 *
 * Material properties: E, nu, h, followed by the base material's own
 * properties (identical layout to the rest of the interface family).
 */
class MarmotZInterfaceMaterialHypoElastic {

protected:
  const double*                                materialProperties;
  const int                                    nMaterialProperties;
  double                                       h = 0.0;
  std::vector< double >                        baseMaterialProperties;
  std::unique_ptr< MarmotMaterialHypoElastic > topMaterial;
  std::unique_ptr< MarmotMaterialHypoElastic > bottomMaterial;

  double shearModulus  = 0.0;
  double lameParameter = 0.0;
  double regularization;

public:
  /** Sizes of the generalized-strain partition x = (w, A+, A-) and of z. */
  static constexpr int nX = 21;
  static constexpr int nZ = 3;

  /** Offsets of the sub-blocks inside x. */
  static constexpr int offW  = 0;
  static constexpr int offAp = 3;
  static constexpr int offAm = 12;

  const int materialNumber;

  MarmotZInterfaceMaterialHypoElastic( const std::string& materialName,
                                       const double*      matProperties_,
                                       int                nMaterialProperties_,
                                       int                materialNumber_ );

  virtual ~MarmotZInterfaceMaterialHypoElastic() = default;

  MarmotStateLayoutDynamic stateLayout;

  double characteristicElementLength = 0.0;

  void setCharacteristicElementLength( double length );

  double getInterfaceThickness() const { return h; }

  /** Relative stiffness zeta of the elastic spring on the z increment. */
  double getGradientJumpRegularization() const { return regularization; }
  void   setGradientJumpRegularization( double zeta );

  /** Current (total) generalized stresses. All row-major. */
  struct State {
    double* generalizedForce;   /**< 3, conjugate to w.  */
    double* surfaceStressPlus;  /**< 9, conjugate to A+. */
    double* surfaceStressMinus; /**< 9, conjugate to A-. */
    double* tractionImbalance;  /**< 3, conjugate to z: (h/4)(t+ - t-) + regularization. */
    double* stateVars;
  };

  /**
   * Hessian of phi over ( x , z ), in four row-major blocks. No block is
   * condensed away and no inverse is ever formed here. With a symmetric
   * bulk tangent the assembled 24x24 is symmetric and K_zx = K_xz^T.
   */
  struct Tangents {
    double* K_xx; /**< 21 x 21 */
    double* K_xz; /**< 21 x  3 */
    double* K_zx; /**<  3 x 21 */
    double* K_zz; /**<  3 x  3 */
  };

  /** Increment inputs over the current step. */
  struct Deformation {
    const double* dU;                  /**< 6:  [ du+ (3) ; du- (3) ] */
    const double* dSurfaceStrain;      /**< 18: [ dA+ (9) ; dA- (9) ], row-major 3x3 each */
    const double* dNormalGradientJump; /**< 3:  dz */
    const double* normal;              /**< 3 */
    const double* separationVector;    /**< 3:  d, top minus bottom */

    Deformation( const double* dU_,
                 const double* dSurfaceStrain_,
                 const double* dNormalGradientJump_,
                 const double* normal_,
                 const double* separationVector_ )
      : dU( dU_ ),
        dSurfaceStrain( dSurfaceStrain_ ),
        dNormalGradientJump( dNormalGradientJump_ ),
        normal( normal_ ),
        separationVector( separationVector_ )
    {
    }

    /** Coincident interface faces: d = 0, so ell = h and d_tau = 0. */
    Deformation( const double* dU_,
                 const double* dSurfaceStrain_,
                 const double* dNormalGradientJump_,
                 const double* normal_ )
      : Deformation( dU_, dSurfaceStrain_, dNormalGradientJump_, normal_, zeroSeparationVector() )
    {
    }

  private:
    static const double* zeroSeparationVector()
    {
      static const double zero[3] = { 0.0, 0.0, 0.0 };
      return zero;
    }
  };

  struct TimeIncrement {
    double timeOld;
    double dT;
  };

  /**
   * Evaluate both faces at the supplied ( x , z ) increment and return the
   * generalized stresses and the four tangent blocks. Mutates the two bulk
   * state blocks and the committed one-sided Cauchy stresses in
   * state.stateVars -- call on a scratch copy for finite differences.
   */
  virtual void computeStress( State&               state,
                              Tangents&            tangents,
                              const Deformation&   deformation,
                              const TimeIncrement& timeIncrement );

  StateView getStateView( const std::string& stateName, double* stateVars ) const
  {
    // Expose the embedded bulk materials' internal variables for
    // postprocessing: a "Top"/"Bottom" suffix resolves into the respective
    // side's state block (e.g. "kappaTop" -> topMaterial's "kappa").
    const auto delegate =
      [&]( const std::string& base, MarmotMaterialHypoElastic* mat, const char* block ) -> StateView {
      return mat->getStateView( base, stateLayout.getPtr( stateVars, block ) );
    };
    const std::size_t n = stateName.size();
    if ( n > 3 && stateName.compare( n - 3, 3, "Top" ) == 0 && topMaterial )
      return delegate( stateName.substr( 0, n - 3 ), topMaterial.get(), "topMaterialStateVars" );
    if ( n > 6 && stateName.compare( n - 6, 6, "Bottom" ) == 0 && bottomMaterial )
      return delegate( stateName.substr( 0, n - 6 ), bottomMaterial.get(), "bottomMaterialStateVars" );

    return stateLayout.getStateView( stateVars, stateName );
  }

  virtual int getNumberOfRequiredStateVars() const { return stateLayout.totalSize(); }

  virtual void initializeYourself( double* stateVars, int nStateVars );

  virtual double getDensity();
};
