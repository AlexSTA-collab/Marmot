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

/**
 * @file YStabPressureMiniInterfaceFiniteElement.h
 * @brief Equilibrated two-sided interface element with a mixed, stabilised nodal pressure and a MINI bubble.
 *
 * This file defines the formulation `Marmot::Elements::StabPressureMiniFormulation` and the element
 * `Marmot::Elements::YStabPressureMiniInterfaceFiniteElement` (YIQUAD4_STABP_MINI), which is built on
 * InterfaceFiniteElementBase like the other interface elements.
 *
 * @par Motivation
 * Near perfect plasticity the von Mises layer reaches the incompressible limit (@f$K/\mu_\mathrm{p} \sim 10^7@f$).
 * A displacement-only formulation evaluates the pressure @f$p = K \operatorname{tr}\varepsilon@f$ from a trace that
 * is driven to zero, so that it multiplies the discretisation error by @f$K@f$ and produces a pressure checkerboard,
 * while the deviatoric stress is uniform. Neither the bending stiffness of the element nor the stability of the
 * interface functional is the cause.
 *
 * @par Formulation
 * The pressure is an independent, @f$C^0@f$-continuous nodal field with one scalar per node of the midsurface,
 * @f$p_h = N \hat p@f$. The equilibrium of the displacements uses @f$\tilde\sigma = \operatorname{dev}\sigma - p_h
 * I@f$, and the volumetric law is enforced weakly together with a Brezzi-Pitkaranta stabilisation of the pressure
 * gradient,
 * @f[
 *   R_u = W^T f + B^{+T} S^+ + B^{-T} S^- ,\qquad
 *   R_p = h \left[ N^T r_p + \tau \, \nabla N^T \nabla N \hat p \right] ,\qquad
 *   \tau = \frac{\gamma h_e^2}{2\mu} .
 * @f]
 * The stabilisation scale is the mesh size @f$h_e@f$ and not a physical length, and the term vanishes with @f$h_e@f$,
 * so it only removes the spurious pressure mode of the equal-order pairing. A mean-dilatation element-local fix
 * cannot do this, because the checkerboard is an element-to-element mode. The unknowns of the element are
 * @f$q_e = [u (24); \hat p (4)]@f$.
 *
 * @par Bubble
 * The displacement of the midsurface is enriched with the bubble
 * @f$b = (1 - \xi^2)(1 - \eta^2)@f$, which vanishes on the boundary of the element, with a vector amplitude
 * @f$\beta@f$. The amplitude is internal to the element. It enters the average surface gradient of both faces
 * equally, and it is condensed at every call: an inner Newton iteration solves @f$R_\beta(\beta) = 0@f$ with a
 * backtracking line search, because @f$\beta@f$ enters the return map and @f$R_\beta@f$ is nonlinear in it.
 *
 * @par State of the trial evaluations
 * The material commits its state in every evaluation, so every trial of the bubble iteration restarts from the
 * state at the beginning of the increment, and a rejected trial of the line search is undone before the next one is
 * evaluated. Otherwise the retry starts from the plastic state of the rejected trial, which makes the local
 * traction-equilibrium solve fail at a hardening modulus that is small against the yield stress.
 *
 * @par Restriction
 * Overriding the pressure of the material is admissible only for a yield surface that is independent of the
 * pressure, see MarmotStabPressureXInterfaceMaterialHypoElastic.
 *
 * @par Environment
 * The Brezzi-Pitkaranta coefficient @f$\gamma@f$ is 0.2. It can be changed with the second element property, or,
 * for sweeps, with the environment variable `MARMOT_STABP_GAMMA`, which the property overrules. @f$\gamma = 0@f$
 * switches the term off and leaves the bubble as the only stabilisation. The variable `MARMOT_MINI_BUBBLE=OFF`
 * holds @f$\beta@f$ at zero, so that the element is the plain equal-order displacement-pressure pairing with the
 * same pressure equation, stabilisation, quadrature and material. It is a diagnostic for ablation studies.
 */
// PLUG-IN VARIANT YIQUAD4_STABP_MINI_SQ (refinement/plugins/psiesqueeze), copied from YStabPressureMiniInterfaceFiniteElement.h
// of the paper build. The pressure equation carries the squeeze-flow term of the layer,
//   R_p += h c_sq grad N^T grad N p ,  c_sq = (h^2 / 12) (1/mu_t+ + 1/mu_t-) / 2 ,
// the resistance of a layer of thickness h and tangent shear modulus mu_t to an in-plane pressure gradient, which the
// linear kinematics through the thickness otherwise drop. mu_t of a face is mu H / (3 mu + H) where its von Mises state
// yielded in the last converged increment and mu where it did not: the yield state of each face is stored per quadrature
// point (squeezeYielded) when an increment is committed, so mu_t is fixed within an increment and the Newton iteration
// sees one pressure equation (deciding it from the current iterate lets it flip between iterations at a reversal). The Brezzi-Pitkaranta term is off by default
// (gamma = 0), so the bubble is the only stabilisation. MARMOT_SQUEEZE=OFF removes the squeeze-flow term; with
// MARMOT_STABP_GAMMA=0.2 the element is then the paper element.
#pragma once

#include "Marmot/InterfaceFiniteElementBase.h"
#include "Marmot/MarmotStabPressureXInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/YInterfaceFiniteElement.h"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <string>
#include <vector>

namespace Marmot::Elements {

  /**
   * @struct StabPressureMiniFormulation
   * @tparam nDim Spatial embedding dimension, 3.
   * @tparam nNodes Number of element nodes, 8.
   * @brief Formulation of the equilibrated two-sided interface element with a mixed, stabilised pressure and a MINI
   * bubble.
   *
   * It supplies the material, the state of a quadrature point, the nodal pressure as additional degree of freedom
   * and the evaluation of the whole element, see the file documentation.
   */
  template < int nDim, int nNodes >
  struct StabPressureMiniSqueezeFormulation {

    static_assert( nDim == 3 && nNodes == 8, "YIQUAD4_STABP_MINI_SQ is a three-dimensional element with 8 nodes" );

    static constexpr int nTensor         = nDim * nDim;
    static constexpr int nInterfaceNodes = nNodes / 2;
    static constexpr int nDofU           = nNodes * nDim;
    static constexpr int halfU           = nDofU / 2;

    /// the pressure is the additional field, one scalar per node of the midsurface
    static constexpr int nDofP      = nInterfaceNodes;
    static constexpr int nExtraDofs = nDofP;

    static constexpr int sizeLoadVector = nDofU + nExtraDofs;

    /// offsets of the displacements and of the pressure in the vectors of the element
    static constexpr int offU = 0;
    static constexpr int offP = nDofU;

    using VectorDim = Eigen::Matrix< double, nDim, 1 >;
    using TensorDim = Eigen::Matrix< double, nDim, nDim >;

    using InterfaceDisplSized       = Eigen::Matrix< double, 2 * nDim, 1 >;
    using InterfaceSurfaceGradSized = Eigen::Matrix< double, 2 * nTensor, 1 >;

    using Material = MarmotStabPressureXInterfaceMaterialHypoElastic; ///< interface material of the element

    /**
     * @brief Named state-variable manager for interface quadrature points.
     */
    class QPStateVarManager : public MarmotStateVarVectorManager {

      /*
       * Persistent state layout of a quadrature point.
       *
       * The amplitude of the bubble is a quantity of the element, parked in the state block of the first
       * quadrature point, because that is the only persistent storage an element is given.
       */
      inline const static auto layout = makeLayout( {
        { .name = "generalizedForce", .length = nDim },
        { .name = "surfaceStressPlus", .length = nTensor },
        { .name = "surfaceStressMinus", .length = nTensor },
        { .name = "volumetricResidual", .length = 1 },
        { .name = "bubbleAlpha", .length = nDim },
        { .name = "squeezeYielded", .length = 2 },
        // The material state starts at a multiple of four entries.
        { .name = "state block alignment padding", .length = ( 4 - ( ( nDim + 2 * nTensor + 1 + nDim + 2 ) % 4 ) ) % 4 },
        { .name = "begin of material state", .length = 0 },
      } );

    public:
      Eigen::Map< Eigen::Matrix< double, nDim, 1 > >    generalizedForce;
      Eigen::Map< Eigen::Matrix< double, nTensor, 1 > > surfaceStressPlus;
      Eigen::Map< Eigen::Matrix< double, nTensor, 1 > > surfaceStressMinus;
      Eigen::Map< Eigen::Matrix< double, 1, 1 > >       volumetricResidual;
      Eigen::Map< Eigen::Matrix< double, nDim, 1 > >    bubbleAlpha;
      Eigen::Map< Eigen::Matrix< double, 2, 1 > >       squeezeYielded;
      Eigen::Map< Eigen::VectorXd >                     materialStateVars;

      /**
       * @brief Number of state variables of the quadrature point itself, without the material state variables.
       *
       * @return the number of state variables up to the begin of the material state
       */
      static int getNumberOfRequiredStateVarsQuadraturePointOnly() { return layout.nRequiredStateVars; }

      /**
       * @brief Map the state variables of one quadrature point.
       *
       * @param[in] theStateVarVector first state variable of the quadrature point
       * @param[in] nStateVars number of state variables of the quadrature point, including the material state
       */
      QPStateVarManager( double* theStateVarVector, int nStateVars )
        : MarmotStateVarVectorManager( theStateVarVector, layout ),
          generalizedForce( &find( "generalizedForce" ) ),
          surfaceStressPlus( &find( "surfaceStressPlus" ) ),
          surfaceStressMinus( &find( "surfaceStressMinus" ) ),
          volumetricResidual( &find( "volumetricResidual" ) ),
          bubbleAlpha( &find( "bubbleAlpha" ) ),
          squeezeYielded( &find( "squeezeYielded" ) ),
          materialStateVars( &find( "begin of material state" ),
                             nStateVars - getNumberOfRequiredStateVarsQuadraturePointOnly() )
      {
      }
    };

    /**
     * @brief Data of the formulation at a quadrature point, besides the state variables.
     */
    struct QuadraturePointData : TwoSidedGeometryData< nDim > {
      /// surface gradient of the bubble b = (1 - xi^2)(1 - eta^2)
      VectorDim gradBubble = VectorDim::Zero();

      /// in-plane mesh size, the scale of the stabilisation
      double meshSize = 0.0;
    };

    /**
     * @brief Nodal fields: the displacement on every node, and the pressure on the nodes of the bottom side.
     *
     * @return the names of the fields of every node
     */
    static std::vector< std::vector< std::string > > nodeFields()
    {
      std::vector< std::vector< std::string > > fields;

      for ( int node = 0; node < nNodes; ++node ) {
        fields.push_back( { "displacement" } );
        if ( node < nInterfaceNodes ) {
          fields[node].push_back( "interfacePressure" );
        }
      }

      return fields;
    }

    /**
     * @brief Place of the degrees of freedom of the element in the degrees of freedom of the nodes.
     *
     * A bottom node holds [u(3), p(1)] in four slots, a top node holds u(3) in three, after the four bottom nodes.
     *
     * @return the position of every degree of freedom of the element
     */
    static std::vector< int > dofIndicesPermutationPattern()
    {
      std::vector< int > permutation( sizeLoadVector );

      for ( int node = 0; node < nInterfaceNodes; ++node ) {
        for ( int component = 0; component < nDim; ++component ) {
          permutation[offU + nDim * node + component]         = ( nDim + 1 ) * node + component;
          permutation[offU + halfU + nDim * node + component] = ( nDim + 1 ) * nInterfaceNodes + nDim * node +
                                                                component;
        }
        permutation[offP + node] = ( nDim + 1 ) * node + nDim;
      }

      return permutation;
    }

    /**
     * @brief Complete the geometry of a quadrature point with the connector between the two faces, the mesh size
     * and the surface gradient of the bubble.
     *
     * The surface gradient of the bubble follows from its derivative with respect to the parent coordinates with
     * the metric of the midsurface, @f$\nabla_s b = J G^{-1} \partial b / \partial \xi@f$.
     *
     * @tparam Element type of the element
     * @tparam QuadraturePoint type of the quadrature point
     * @param[in] element the element, which holds the coordinates of the nodes
     * @param[in,out] qp quadrature point, whose normal and projections may be reoriented
     *
     * @throws std::invalid_argument if the connector has no positive component along the normal
     */
    template < class Element, class QuadraturePoint >
    static void initializeQuadraturePoint( const Element& element, QuadraturePoint& qp )
    {
      initializeTwoSidedQuadraturePoint( element, qp, "YIQUAD4_STABP_MINI_SQ" );

      auto& data    = qp.formulationData;
      data.meshSize = std::sqrt( qp.sqrtDetG );

      const double xi  = qp.xi( 0 );
      const double eta = qp.xi( 1 );

      Eigen::Matrix< double, 2, 1 > dbdxi;
      dbdxi( 0 ) = -2.0 * xi * ( 1.0 - eta * eta );
      dbdxi( 1 ) = -2.0 * eta * ( 1.0 - xi * xi );

      data.gradBubble = qp.J * qp.G.inverse() * dbdxi;
    }

    /**
     * @brief Brezzi-Pitkaranta coefficient of the stabilisation.
     *
     * In order of precedence: the second element property, the environment variable `MARMOT_STABP_GAMMA`, 0 (this
     * variant: the squeeze-flow term replaces it).
     * Zero is admissible and removes the term.
     *
     * @tparam Element type of the element
     * @param[in] element the element, which holds its properties
     * @return the coefficient
     */
    template < class Element >
    static double stabilisationCoefficient( const Element& element )
    {
      if ( element.elementProperties.size() > 1 && element.elementProperties[1] > 0.0 ) {
        return element.elementProperties[1];
      }

      static const double environmentCoefficient = [] {
        if ( const char* value = std::getenv( "MARMOT_STABP_GAMMA" ) ) {
          char*        end    = nullptr;
          const double parsed = std::strtod( value, &end );
          if ( end != value && *end == '\0' && parsed >= 0.0 ) {
            return parsed;
          }
        }
        return 0.0;
      }();

      return environmentCoefficient;
    }

    /**
     * @brief Whether the bubble is held at zero, see the file documentation.
     *
     * @return true if the environment variable `MARMOT_MINI_BUBBLE` is `OFF`
     */
    static bool bubbleDisabled()
    {
      static const bool disabled = [] {
        const char* value = std::getenv( "MARMOT_MINI_BUBBLE" );
        return value && std::string( value ) == "OFF";
      }();

      return disabled;
    }

    /**
     * @brief Whether the squeeze-flow term is removed.
     *
     * @return true if the environment variable `MARMOT_SQUEEZE` is `OFF`
     */
    static bool squeezeDisabled()
    {
      static const bool disabled = [] {
        const char* value = std::getenv( "MARMOT_SQUEEZE" );
        return value && std::string( value ) == "OFF";
      }();

      return disabled;
    }

    /**
     * @brief Hardening modulus H of the von Mises base material, the fifth property [E, nu, h, sigma_y, H, ...] of the
     * interface material. The properties are protected; a pointer to member reaches them without changing Marmot.
     *
     * @param[in] material the interface material of the element
     * @return H, or a negative value if the material carries fewer than five properties
     */
    static double hardeningModulus( const Material& material )
    {
      struct Access : Material {
        static const double* Material::* properties() { return &Access::materialProperties; }
        static const int Material::*     count() { return &Access::nMaterialProperties; }
      };
      if ( material.*Access::count() < 5 ) {
        return -1.0;
      }
      return ( material.*Access::properties() )[4];
    }

    /**
     * @brief Evaluate the whole element: the inner iteration on the bubble, its condensation, the residual and the
     * tangent.
     *
     * The residual is the internal force, so that the tangent is its derivative. The residual of the bubble is
     * @f$R_\beta = \sum_q B_\beta^T (S^+ + S^-) J_0 W@f$. It is solved to zero, and then the displacements and
     * the pressure are condensed, @f$K \leftarrow K_{xx} - K_{x\beta} K_{\beta\beta}^{-1} K_{\beta x}@f$.
     *
     * @tparam Element type of the element
     * @tparam QType type of the nodal values
     * @tparam PeType type of the element residual
     * @tparam KeType type of the element tangent
     * @param[in,out] element the element, whose quadrature points hold the state
     * @param[in] QTotal total nodal values of the element, displacements first and then the pressure
     * @param[in] dQ increment of the nodal values of the element
     * @param[in] time time at the beginning of the increment
     * @param[in] dT time increment
     * @param[in,out] Pe element residual, to which the contribution is added
     * @param[in,out] Ke element tangent, to which the contribution is added
     */
    template < class Element, class QType, class PeType, class KeType >
    static void computeElementKernels( Element&     element,
                                       const QType& QTotal,
                                       const QType& dQ,
                                       double       time,
                                       double       dT,
                                       PeType&      Pe,
                                       KeType&      Ke )
    {
      using V3       = Eigen::Matrix< double, 3, 1 >;
      using M3       = Eigen::Matrix< double, 3, 3 >;
      using RhsSized = Eigen::Matrix< double, sizeLoadVector, 1 >;
      using KeSized  = Eigen::Matrix< double, sizeLoadVector, sizeLoadVector >;

      auto&        qps      = element.qps;
      const size_t nQp      = qps.size();
      const double gamma    = stabilisationCoefficient( element );
      const auto   pressure = QTotal.template segment< nDofP >( offP );
      const auto   dP       = dQ.template segment< nDofP >( offP );

      // The material commits its state in every evaluation, so every pass starts from the state at the beginning
      // of the increment.
      std::vector< std::vector< double > > materialStateIn( nQp );
      for ( size_t iq = 0; iq < nQp; ++iq ) {
        const auto& state   = qps[iq].managedStateVars->materialStateVars;
        materialStateIn[iq] = std::vector< double >( state.data(), state.data() + state.size() );
      }
      // yield state of the two faces at the end of the last converged increment: fixes mu_t within the increment
      std::vector< std::array< double, 2 > > yieldedIn( nQp );
      for ( size_t iq = 0; iq < nQp; ++iq ) {
        yieldedIn[iq] = { qps[iq].managedStateVars->squeezeYielded( 0 ), qps[iq].managedStateVars->squeezeYielded( 1 ) };
      }

      // Everything one pass over the quadrature points produces.
      struct Pass {
        RhsSized                                       Pe;
        KeSized                                        Ke;
        V3                                             Rb;
        M3                                             Kbb;
        Eigen::Matrix< double, 3, nDofU >              Kbu;
        Eigen::Matrix< double, nDofU, 3 >              Kub;
        Eigen::Matrix< double, 3, nDofP >              Kbp;
        Eigen::Matrix< double, nDofP, 3 >              Kpb;
        std::array< V3, 8 >                            force;
        std::array< Eigen::Matrix< double, 9, 1 >, 8 > surfacePlus;
        std::array< Eigen::Matrix< double, 9, 1 >, 8 > surfaceMinus;
        std::array< double, 8 >                        volumetricResidual;
        std::array< std::array< double, 2 >, 8 >       yielded{};
      };

      auto assemble = [&]( const V3& beta ) {
        Pass pass;
        pass.Pe.setZero();
        pass.Ke.setZero();
        pass.Rb.setZero();
        pass.Kbb.setZero();
        pass.Kbu.setZero();
        pass.Kub.setZero();
        pass.Kbp.setZero();
        pass.Kpb.setZero();

        for ( size_t iq = 0; iq < nQp; ++iq ) {
          auto& qp   = qps[iq];
          auto& data = qp.formulationData;
          auto& msv  = qp.managedStateVars->materialStateVars;
          std::copy( materialStateIn[iq].begin(), materialStateIn[iq].end(), msv.data() );

          const auto& Njump = qp.NmatJump;
          const auto& Bside = qp.BmatSide;
          const auto& Np    = qp.N;
          const auto& Gp    = qp.gradN;

          Eigen::Matrix< double, nTensor, nDofU > BPlus      = Eigen::Matrix< double, nTensor, nDofU >::Zero();
          Eigen::Matrix< double, nTensor, nDofU > BMinus     = Eigen::Matrix< double, nTensor, nDofU >::Zero();
          BMinus.template block< nTensor, halfU >( 0, 0 )    = Bside;
          BPlus.template block< nTensor, halfU >( 0, halfU ) = Bside;

          // the bubble perturbs the mean in-plane displacement: Abar_ij += beta_i (grad_s b)_j, equally on both faces
          Eigen::Matrix< double, nTensor, nDim > Bbubble = Eigen::Matrix< double, nTensor, nDim >::Zero();
          for ( int i = 0; i < nDim; ++i ) {
            for ( int j = 0; j < nDim; ++j ) {
              Bbubble( nDim * i + j, i ) = data.gradBubble( j );
            }
          }

          const InterfaceDisplSized dUSides = qp.interpolateToSides( dQ );

          InterfaceSurfaceGradSized dSurface = qp.surfaceGradientsAtSides( dQ );
          dSurface.template segment< nTensor >( 0 ) += Bbubble * beta;
          dSurface.template segment< nTensor >( nTensor ) += Bbubble * beta;

          element.material->setCharacteristicElementLength( qp.characteristicLength );

          typename Material::State state{ Material::Tensor3d( 0.0 ),
                                          Material::Tensor33d( 0.0 ),
                                          Material::Tensor33d( 0.0 ),
                                          0.0,
                                          msv.data() };

          typename Material::Tangents tangents{ Material::Tensor33d( 0.0 ),
                                                Material::Tensor333d( 0.0 ),
                                                Material::Tensor333d( 0.0 ),
                                                Material::Tensor333d( 0.0 ),
                                                Material::Tensor3333d( 0.0 ),
                                                Material::Tensor3333d( 0.0 ),
                                                Material::Tensor333d( 0.0 ),
                                                Material::Tensor3333d( 0.0 ),
                                                Material::Tensor3333d( 0.0 ),
                                                Material::Tensor3d( 0.0 ),
                                                Material::Tensor33d( 0.0 ),
                                                Material::Tensor33d( 0.0 ),
                                                Material::Tensor3d( 0.0 ),
                                                Material::Tensor33d( 0.0 ),
                                                Material::Tensor33d( 0.0 ),
                                                0.0 };

          typename Material::Deformation deformation{ Material::Tensor6d( dUSides.data() ),
                                                      Material::Tensor18d( dSurface.data() ),
                                                      Material::Tensor3d( qp.normal.data() ),
                                                      Material::Tensor3d( data.separationVector.data() ),
                                                      ( Np * dP )( 0, 0 ) };

          element.material->computeStress( state, tangents, deformation, typename Material::TimeIncrement{ time, dT } );

          // the response, as Eigen objects
          const V3                            f  = Eigen::Map< const V3 >( state.generalizedForce.data() );
          const Eigen::Matrix< double, 9, 1 > Sp = Eigen::Map< const Eigen::Matrix< double, 9, 1 > >(
            state.surfaceStressPlus.data() );
          const Eigen::Matrix< double, 9, 1 > Sm = Eigen::Map< const Eigen::Matrix< double, 9, 1 > >(
            state.surfaceStressMinus.data() );
          const double rp = state.volumetricResidual;

          const Eigen::Map< const Eigen::Matrix< double, 3, 3, Eigen::RowMajor > > Qww( tangents.Q_ww.data() );
          const Eigen::Map< const Eigen::Matrix< double, 3, 9, Eigen::RowMajor > > QwAp( tangents.Q_wAp.data() );
          const Eigen::Map< const Eigen::Matrix< double, 3, 9, Eigen::RowMajor > > QwAm( tangents.Q_wAm.data() );
          const Eigen::Map< const Eigen::Matrix< double, 9, 3, Eigen::RowMajor > > QApw( tangents.Q_Apw.data() );
          const Eigen::Map< const Eigen::Matrix< double, 9, 9, Eigen::RowMajor > > QApAp( tangents.Q_ApAp.data() );
          const Eigen::Map< const Eigen::Matrix< double, 9, 9, Eigen::RowMajor > > QApAm( tangents.Q_ApAm.data() );
          const Eigen::Map< const Eigen::Matrix< double, 9, 3, Eigen::RowMajor > > QAmw( tangents.Q_Amw.data() );
          const Eigen::Map< const Eigen::Matrix< double, 9, 9, Eigen::RowMajor > > QAmAp( tangents.Q_AmAp.data() );
          const Eigen::Map< const Eigen::Matrix< double, 9, 9, Eigen::RowMajor > > QAmAm( tangents.Q_AmAm.data() );
          const Eigen::Map< const V3 >                                             Qwp( tangents.Q_wp.data() );
          const Eigen::Map< const Eigen::Matrix< double, 9, 1 > >                  QApp( tangents.Q_App.data() );
          const Eigen::Map< const Eigen::Matrix< double, 9, 1 > >                  QAmp( tangents.Q_Amp.data() );
          const Eigen::Map< const Eigen::Matrix< double, 1, 3 > >                  Qpw( tangents.Q_pw.data() );
          const Eigen::Map< const Eigen::Matrix< double, 1, 9 > >                  QpAp( tangents.Q_pAp.data() );
          const Eigen::Map< const Eigen::Matrix< double, 1, 9 > >                  QpAm( tangents.Q_pAm.data() );
          const double                                                             Qpp = tangents.Q_pp;

          pass.force[iq]              = f;
          pass.surfacePlus[iq]        = Sp;
          pass.surfaceMinus[iq]       = Sm;
          pass.volumetricResidual[iq] = rp;

          const double thickness = element.material->getInterfaceThickness();
          const double mu        = element.material->getShearModulus();
          double       tau       = gamma * data.meshSize * data.meshSize / ( 2.0 * mu );
          if ( !squeezeDisabled() ) {
            // tangent shear modulus of each face from its yield state in the last converged increment; the yield state
            // of this pass (kappa grows) is stored for the next increment when the pass is committed
            const double H        = hardeningModulus( *element.material );
            const double muPlastic = H > 0.0 ? mu * H / ( 3.0 * mu + H ) : mu;
            double       compliance = 0.0;
            int          side_      = 0;
            for ( const char* side : { "kappaTop", "kappaBottom" } ) {
              const StateView view   = element.material->getStateView( side, msv.data() );
              const std::ptrdiff_t k = view.stateLocation - msv.data();
              const bool yielding    = view.stateSize == 1 && msv( k ) > materialStateIn[iq][k];
              pass.yielded[iq][side_] = yielding ? 1.0 : 0.0;
              compliance += 0.5 / ( yieldedIn[iq][side_] > 0.5 ? muPlastic : mu );
              ++side_;
            }
            // R_p = h [ N^T r_p + tau grad N^T grad N p ]: tau = h^2 / 12 (1/mu_t)_mean gives h^3 / 12 (1/mu_t)_mean
            tau += thickness * thickness / 12.0 * compliance;
          }

          const Eigen::Matrix< double, nDofU, 1 > Ru = Njump.transpose() * f + BPlus.transpose() * Sp +
                                                       BMinus.transpose() * Sm;
          const Eigen::Matrix< double, nDofP, 1 > Rp = thickness * ( Np.transpose() * rp +
                                                                     tau * ( Gp.transpose() * ( Gp * pressure ) ) );

          pass.Pe.template segment< nDofU >( offU ) += Ru * qp.J0xW;
          pass.Pe.template segment< nDofP >( offP ) += Rp * qp.J0xW;
          pass.Rb += ( Bbubble.transpose() * ( Sp + Sm ) ) * qp.J0xW;

          const Eigen::Matrix< double, nDofU, nDofU > Kuu = Njump.transpose() * Qww * Njump +
                                                            Njump.transpose() * QwAp * BPlus +
                                                            Njump.transpose() * QwAm * BMinus +
                                                            BPlus.transpose() * QApw * Njump +
                                                            BPlus.transpose() * QApAp * BPlus +
                                                            BPlus.transpose() * QApAm * BMinus +
                                                            BMinus.transpose() * QAmw * Njump +
                                                            BMinus.transpose() * QAmAp * BPlus +
                                                            BMinus.transpose() * QAmAm * BMinus;
          const Eigen::Matrix< double, nDofU, nDofP > Kup = ( Njump.transpose() * Qwp + BPlus.transpose() * QApp +
                                                              BMinus.transpose() * QAmp ) *
                                                            Np;
          const Eigen::Matrix< double, nDofP, nDofU > Kpu = thickness *
                                                            ( Np.transpose() *
                                                              ( Qpw * Njump + QpAp * BPlus + QpAm * BMinus ) );
          const Eigen::Matrix< double, nDofP, nDofP > Kpp = thickness * ( Np.transpose() * Qpp * Np +
                                                                          tau * ( Gp.transpose() * Gp ) );

          pass.Ke.template block< nDofU, nDofU >( offU, offU ) += Kuu * qp.J0xW;
          pass.Ke.template block< nDofU, nDofP >( offU, offP ) += Kup * qp.J0xW;
          pass.Ke.template block< nDofP, nDofU >( offP, offU ) += Kpu * qp.J0xW;
          pass.Ke.template block< nDofP, nDofP >( offP, offP ) += Kpp * qp.J0xW;

          // the coupling of the bubble with the displacements and the pressure
          const Eigen::Matrix< double, nTensor, nDim > dSpdb = ( QApAp + QApAm ) * Bbubble;
          const Eigen::Matrix< double, nTensor, nDim > dSmdb = ( QAmAp + QAmAm ) * Bbubble;
          const Eigen::Matrix< double, nDim, nDim >    dfdb  = ( QwAp + QwAm ) * Bbubble;

          pass.Kbb += ( Bbubble.transpose() * ( dSpdb + dSmdb ) ) * qp.J0xW;
          pass.Kbu += ( Bbubble.transpose() *
                        ( ( QApw + QAmw ) * Njump + ( QApAp + QAmAp ) * BPlus + ( QApAm + QAmAm ) * BMinus ) ) *
                      qp.J0xW;
          pass.Kub += ( Njump.transpose() * dfdb + BPlus.transpose() * dSpdb + BMinus.transpose() * dSmdb ) * qp.J0xW;
          pass.Kbp += ( Bbubble.transpose() * ( QApp + QAmp ) * Np ) * qp.J0xW;
          pass.Kpb += ( thickness * ( Np.transpose() * ( QpAp + QpAm ) * Bbubble ) ) * qp.J0xW;
        }

        return pass;
      };

      // The material state left in the state vector after a pass belongs to the amplitude of that pass. It is
      // saved whenever a pass is accepted, so that a rejected trial is undone without evaluating the material again.
      std::vector< std::vector< double > > acceptedMaterialState( nQp );

      auto saveAcceptedState = [&]() {
        for ( size_t iq = 0; iq < nQp; ++iq ) {
          const auto& state         = qps[iq].managedStateVars->materialStateVars;
          acceptedMaterialState[iq] = std::vector< double >( state.data(), state.data() + state.size() );
        }
      };

      auto restoreAcceptedState = [&]() {
        for ( size_t iq = 0; iq < nQp; ++iq ) {
          auto& state = qps[iq].managedStateVars->materialStateVars;
          std::copy( acceptedMaterialState[iq].begin(), acceptedMaterialState[iq].end(), state.data() );
        }
      };

      auto commit = [&]( const Pass& pass ) {
        for ( size_t iq = 0; iq < nQp; ++iq ) {
          auto& managed = *qps[iq].managedStateVars;
          std::copy( acceptedMaterialState[iq].begin(),
                     acceptedMaterialState[iq].end(),
                     managed.materialStateVars.data() );
          managed.generalizedForce        = pass.force[iq];
          managed.surfaceStressPlus       = pass.surfacePlus[iq];
          managed.surfaceStressMinus      = pass.surfaceMinus[iq];
          managed.volumetricResidual( 0 ) = pass.volumetricResidual[iq];
          managed.squeezeYielded( 0 )     = pass.yielded[iq][0];
          managed.squeezeYielded( 1 )     = pass.yielded[iq][1];
        }
      };

      if ( bubbleDisabled() ) {
        const Pass pass = assemble( V3::Zero() );
        saveAcceptedState();
        commit( pass );
        qps[0].managedStateVars->bubbleAlpha = V3::Zero();
        Pe += pass.Pe;
        Ke += pass.Ke;
        return;
      }

      // Inner Newton iteration on the amplitude of the bubble, R_beta(beta) = 0. The amplitude enters the strain and
      // therefore the return map, so that R_beta is nonlinear in it, and a single step from zero would leave the
      // residual inconsistent with the condensed tangent. The backtracking line search is required for robustness: a
      // full step can increase the residual at the onset of yield, drive the strain to a state in which the local solve
      // of the material cannot converge, and the failure propagates as a cutback of the increment. A trial that throws
      // is treated like a trial that does not improve.
      V3   beta = qps[0].managedStateVars->bubbleAlpha;
      Pass pass = assemble( beta );
      saveAcceptedState();

      for ( int iteration = 0; iteration < 20; ++iteration ) {
        const double residualNorm = pass.Rb.norm();
        const double tolerance    = 1.0e-10 * std::max( 1.0, pass.Pe.norm() );
        if ( residualNorm <= tolerance ) {
          break;
        }

        const V3 dBeta = -pass.Kbb.fullPivLu().solve( pass.Rb );

        double alpha    = 1.0;
        bool   accepted = false;
        for ( int lineSearchIteration = 0; lineSearchIteration < 10; ++lineSearchIteration ) {
          restoreAcceptedState();
          try {
            const Pass trial = assemble( beta + alpha * dBeta );
            if ( trial.Rb.norm() <= tolerance || trial.Rb.norm() < residualNorm ) {
              beta += alpha * dBeta;
              pass = trial;
              saveAcceptedState();
              accepted = true;
              break;
            }
          }
          catch ( const std::exception& ) {
            // the trial drove the local solve of the material past convergence: treat it like a step that does not
            // improve, and halve it
          }
          alpha *= 0.5;
        }

        if ( !accepted ) {
          break; // keep the last good amplitude; the global iteration continues from there
        }
      }

      // The accepted pass evaluated everything at this amplitude, so the state is written from it. The material
      // state is restored from the saved copy, which undoes a rejected trial that ran afterwards.
      commit( pass );

      // one factorisation for the three right-hand sides
      const auto                              KbbFactorisation = pass.Kbb.fullPivLu();
      const V3                                KbbInvRb         = KbbFactorisation.solve( pass.Rb );
      const Eigen::Matrix< double, 3, nDofU > KbbInvKbu        = KbbFactorisation.solve( pass.Kbu );
      const Eigen::Matrix< double, 3, nDofP > KbbInvKbp        = KbbFactorisation.solve( pass.Kbp );

      // R_beta is zero after the iteration, so the condensed residual is the residual of the displacements and the
      // pressure; the term of the remainder keeps the residual consistent with the tangent if the iteration stopped
      Pe.template segment< nDofU >( offU ) += pass.Pe.template segment< nDofU >( offU ) - pass.Kub * KbbInvRb;
      Pe.template segment< nDofP >( offP ) += pass.Pe.template segment< nDofP >( offP ) - pass.Kpb * KbbInvRb;

      Ke += pass.Ke;
      Ke.template block< nDofU, nDofU >( offU, offU ) -= pass.Kub * KbbInvKbu;
      Ke.template block< nDofU, nDofP >( offU, offP ) -= pass.Kub * KbbInvKbp;
      Ke.template block< nDofP, nDofU >( offP, offU ) -= pass.Kpb * KbbInvKbu;
      Ke.template block< nDofP, nDofP >( offP, offP ) -= pass.Kpb * KbbInvKbp;

      qps[0].managedStateVars->bubbleAlpha = beta;
    }
  };

  /**
   * @brief Equilibrated two-sided interface element with a mixed, stabilised nodal pressure and a MINI bubble.
   *
   * @tparam nDim Spatial embedding dimension, 3.
   * @tparam nNodes Number of element nodes, 8.
   */
  template < int nDim, int nNodes >
  using YStabPressureMiniSqueezeInterfaceFiniteElement = InterfaceFiniteElementBase<
    nDim,
    nNodes,
    StabPressureMiniSqueezeFormulation< nDim, nNodes > >;

} // namespace Marmot::Elements
