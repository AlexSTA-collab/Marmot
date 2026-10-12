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
 * @file GYSInterfaceFiniteElementFiniteStrain.h
 * @brief Gradient-enhanced interface finite element at finite strain with a nonlocal field at every station through
 * the thickness.
 *
 * This file defines the formulation `Marmot::Elements::GradientEnhancedStationFieldsFormulation` and the element
 * `Marmot::Elements::GYSInterfaceFiniteElementFiniteStrain`, which is the common base of the interface elements
 * (InterfaceFiniteElementBase.h) with this formulation. The layer is divided by stations through its thickness, and
 * every station carries its own bulk material and its own nonlocal field, so that the damage varies through the
 * thickness. It generalizes the element with a nonlocal field for each face
 * (`GradientEnhancedFaceFieldsFormulation`), which has two stations.
 */
#pragma once

#include "Marmot/InterfaceFiniteElementBase.h"
#include "Marmot/MarmotGradedGradientEnhancedInterfaceMaterialFiniteStrain.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/YInterfaceFiniteElement.h"

#include <Eigen/Dense>

#include <string>
#include <vector>

namespace Marmot::Elements {

  /**
   * @struct GradientEnhancedStationFieldsFormulation
   * @tparam nDim Spatial embedding dimension.
   * @tparam nNodes Number of element nodes.
   * @tparam nStations Number of stations through the thickness, which is the number of nonlocal fields; at least two.
   * @brief Formulation of the interface element with the kinematics of the displacement jump at finite strain and a
   * nonlocal field at every station through the thickness.
   *
   * The interface material (`MarmotGradedGradientEnhancedInterfaceMaterialFiniteStrain`) receives the total
   * displacements and surface gradients of both interface sides, the connector between the faces, and the nonlocal
   * fields of all stations at the quadrature point. It returns the force, the surface stress, the local driving
   * field and the nonlocal radius of every station, together with the condensed tangents.
   *
   * The stations @f$ y_s @f$ run from the bottom face to the top face. The nonlocal fields interpolate linearly
   * between them, and the station @f$ s @f$ carries the weight @f$ \omega_s @f$ of the trapezoidal rule, so that its
   * nonlocal equation, reduced to the midsurface, is
   *
   * @f[
   *   h \omega_s \left[ \mathbf{N}^\mathsf{T} ( \bar N_s - L_s ) + l_s^2 \nabla\mathbf{N}^\mathsf{T}
   *   \nabla\bar N_s \right] + \frac{l_s^2}{h} \left[ \frac{\mathbf{N}^\mathsf{T} ( \bar N_s - \bar N_{s+1} )}
   *   {y_{s+1} - y_s} + \frac{\mathbf{N}^\mathsf{T} ( \bar N_s - \bar N_{s-1} )}{y_s - y_{s-1}} \right] = 0 ,
   * @f]
   *
   * where @f$ L_s @f$ is the local driving field of the station. The last term is the derivative of the field through
   * the thickness, which couples the neighbouring stations, with the nonlocal radius of the station. For two
   * stations it is the element with a nonlocal field for each face.
   *
   * The degrees of freedom of the nonlocal fields are placed after the displacements in the vectors of the element,
   * the bottom station first, each with one degree of freedom on every node of a side. On the nodes of the element,
   * the fields of the bottom station and of the inner stations belong to the nodes of the bottom side, and the field
   * of the top station to the nodes of the top side.
   */
  template < int nDim, int nNodes, int nStations >
  struct GradientEnhancedStationFieldsFormulation {

    static_assert( nStations >= 2, "A layer with a nonlocal field at its stations needs at least two stations" );

    static constexpr int nTensor         = nDim * nDim;
    static constexpr int nInterfaceNodes = nNodes / 2;
    static constexpr int nDofD           = nNodes * nDim;   ///< degrees of freedom of the displacements
    static constexpr int nDofA           = nInterfaceNodes; ///< degrees of freedom of the field of a station

    /// number of degrees of freedom in addition to the displacements: a nonlocal field per station on every node
    static constexpr int nExtraDofs = nStations * nDofA;

    using ForceSized                = Eigen::Matrix< double, nDim, 1 >;
    using SurfaceStressSized        = Eigen::Matrix< double, nTensor, 1 >;
    using InterfaceDisplSized       = Eigen::Matrix< double, 2 * nDim, 1 >;
    using InterfaceSurfaceGradSized = Eigen::Matrix< double, 2 * nTensor, 1 >;

    using QMatrixSized      = Eigen::Matrix< double, nDim, nDim, Eigen::RowMajor >;
    using HMatrixSized      = Eigen::Matrix< double, nDim, nTensor, Eigen::RowMajor >;
    using HAverageJumpSized = Eigen::Matrix< double, nTensor, nDim, Eigen::RowMajor >;
    using ZMatrixSized      = Eigen::Matrix< double, nTensor, nTensor, Eigen::RowMajor >;

    using Material = MarmotGradedGradientEnhancedInterfaceMaterialFiniteStrain; ///< interface material

    /**
     * @brief Named state-variable manager for interface quadrature points.
     */
    class QPStateVarManager : public MarmotStateVarVectorManager {

      /*
       * Persistent state layout for the force, the surface stress, the total displacement, the total surface strain,
       * the local and the nonlocal field of every station, and the material state variables.
       */
      inline const static auto layout = makeLayout( {
        { .name = "force", .length = nDim },
        { .name = "surface stress", .length = nTensor },
        { .name = "displacement", .length = 2 * nDim },
        { .name = "surface strain", .length = 2 * nTensor },
        { .name = "local fields", .length = nStations },
        { .name = "nonlocal fields", .length = nStations },
        // The material state starts at a multiple of four entries.
        { .name   = "state block alignment padding",
          .length = ( 4 - ( ( nDim + nTensor + 2 * nDim + 2 * nTensor + 2 * nStations ) % 4 ) ) % 4 },
        { .name = "begin of material state", .length = 0 },
      } );

    public:
      Eigen::Map< ForceSized >                            force;
      Eigen::Map< SurfaceStressSized >                    surfaceStress;
      Eigen::Map< InterfaceDisplSized >                   displacement;
      Eigen::Map< InterfaceSurfaceGradSized >             surfaceStrain;
      Eigen::Map< Eigen::Matrix< double, nStations, 1 > > localFields;
      Eigen::Map< Eigen::Matrix< double, nStations, 1 > > nonlocalFields;

      Eigen::Map< Eigen::VectorXd > materialStateVars;

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
          force( &find( "force" ) ),
          surfaceStress( &find( "surface stress" ) ),
          displacement( &find( "displacement" ) ),
          surfaceStrain( &find( "surface strain" ) ),
          localFields( &find( "local fields" ) ),
          nonlocalFields( &find( "nonlocal fields" ) ),
          materialStateVars( &find( "begin of material state" ),
                             nStateVars - getNumberOfRequiredStateVarsQuadraturePointOnly() )
      {
      }
    };

    /// This formulation keeps the connector between the faces at the quadrature points.
    struct QuadraturePointData : TwoSidedGeometryData< nDim > {};

    /**
     * @brief Complete the geometry of a quadrature point with the connector between the two faces.
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
      initializeTwoSidedQuadraturePoint( element, qp, "GYSInterfaceFiniteElementFiniteStrain" );
    }

    /**
     * @brief Fields on the nodes of the element.
     *
     * The nodes of the bottom side carry the displacement and the fields of the bottom and of the inner stations,
     * the nodes of the top side the displacement and the field of the top station.
     *
     * @return the names of the fields of every node
     */
    static std::vector< std::vector< std::string > > nodeFields()
    {
      std::vector< std::string > bottom = { "displacement", "nonlocal damage" };
      for ( int s = 1; s < nStations - 1; ++s )
        bottom.push_back( "nonlocal damage " + std::to_string( s ) );

      std::vector< std::vector< std::string > > fields;
      for ( int A = 0; A < nInterfaceNodes; ++A )
        fields.push_back( bottom );
      for ( int A = 0; A < nInterfaceNodes; ++A )
        fields.push_back( { "displacement", "nonlocal damage" } );
      return fields;
    }

    /**
     * @brief Permutation from the order of the degrees of freedom of the element to the order node by node.
     *
     * @return for every degree of freedom of the element its index in the order node by node
     */
    static std::vector< int > dofIndicesPermutationPattern()
    {
      std::vector< int > permutation( nDofD + nExtraDofs );

      // a bottom node carries [ u(nDim), N_0, ..., N_(nStations-2) ], a top node [ u(nDim), N_(nStations-1) ]
      constexpr int bottomSize = nDim + nStations - 1;
      constexpr int topSize    = nDim + 1;

      for ( int A = 0; A < nInterfaceNodes; ++A ) {
        const int bottomStart = A * bottomSize;
        const int topStart    = nInterfaceNodes * bottomSize + A * topSize;

        for ( int c = 0; c < nDim; ++c ) {
          permutation[nDim * A + c]                          = bottomStart + c;
          permutation[nDim * nInterfaceNodes + nDim * A + c] = topStart + c;
        }
        for ( int s = 0; s < nStations - 1; ++s )
          permutation[nDofD + s * nDofA + A] = bottomStart + nDim + s;
        permutation[nDofD + ( nStations - 1 ) * nDofA + A] = topStart + nDim;
      }
      return permutation;
    }

    /**
     * @brief Evaluate the material at a quadrature point and add its contributions to the element.
     *
     * @tparam Element type of the element
     * @tparam QuadraturePoint type of the quadrature point
     * @tparam QType type of the total nodal values
     * @tparam PeType type of the element residual
     * @tparam KeType type of the element tangent
     * @param[in] element the element, which holds the material
     * @param[in,out] qp quadrature point, whose state variables are updated
     * @param[in] QTotal total nodal displacements and nonlocal fields of the element at the end of the increment
     * @param[in] dU_GPs increment of the displacements of both interface sides, which this formulation does not need
     * @param[in] dSurface_strain_GPs increment of the surface gradients, which this formulation does not need
     * @param[in] time time at the beginning of the increment
     * @param[in] dT time increment
     * @param[in,out] Pe element residual, to which the contribution is added
     * @param[in,out] Ke element tangent, to which the contribution is added
     */
    template < class Element, class QuadraturePoint, class QType, class PeType, class KeType >
    static void integrate( Element&                         element,
                           QuadraturePoint&                 qp,
                           const QType&                     QTotal,
                           const InterfaceDisplSized&       dU_GPs,
                           const InterfaceSurfaceGradSized& dSurface_strain_GPs,
                           double                           time,
                           double                           dT,
                           PeType&                          Pe,
                           KeType&                          Ke )
    {
      (void)dU_GPs;
      (void)dSurface_strain_GPs;

      using namespace Marmot;

      auto& material = element.material;
      auto& data     = qp.formulationData;

      if ( material->numberOfStations() != nStations )
        throw std::invalid_argument(
          "GYSInterfaceFiniteElementFiniteStrain: the material has " + std::to_string( material->numberOfStations() ) +
          " stations, but the element has " + std::to_string( nStations ) + " nonlocal fields" );

      const auto& Njump = qp.NmatJump;
      const auto& Bavg  = qp.BmatAverage;

      const InterfaceDisplSized       totalU_GPs               = qp.interpolateToSides( QTotal );
      const InterfaceSurfaceGradSized totalSurfaceGradient_GPs = qp.surfaceGradientsAtSides( QTotal );

      // nonlocal fields of the stations at this point
      std::array< Eigen::Matrix< double, nDofA, 1 >, nStations > nodalNonlocal;
      Eigen::VectorXd                                            nonlocal( nStations );
      for ( int s = 0; s < nStations; ++s ) {
        nodalNonlocal[s] = QTotal.template segment< nDofA >( nDofD + s * nDofA );
        nonlocal( s )    = ( qp.N * nodalNonlocal[s] )( 0 );
      }

      // the materials are fixed to three dimensions, a two-dimensional element embeds its kinematics
      const Eigen::Matrix< double, 6, 1 >  U3d               = expandSidesTo3D< nDim, 1 >( totalU_GPs );
      const Eigen::Matrix< double, 18, 1 > surfaceGradient3d = expandSidesTo3D< nDim, 2 >( totalSurfaceGradient_GPs );

      Material::State state{ Material::Tensor3d( 0.0 ),
                             Material::Tensor33d( 0.0 ),
                             qp.managedStateVars->materialStateVars.data(),
                             Eigen::VectorXd(),
                             Eigen::VectorXd() };

      Material::Tangents tangents{ Material::Tensor33d( 0.0 ),
                                   Material::Tensor333d( 0.0 ),
                                   Material::Tensor333d( 0.0 ),
                                   Material::Tensor3333d( 0.0 ),
                                   {},
                                   {},
                                   {},
                                   {},
                                   {} };

      Material::Deformation deformation{ Material::Tensor6d( U3d.data() ),
                                         Material::Tensor18d( surfaceGradient3d.data() ),
                                         embedVectorIn3D< nDim >( qp.normal.data() ),
                                         embedVectorIn3D< nDim >( data.separationVector.data() ),
                                         nonlocal };

      material->computeStress( state, tangents, deformation, { time, dT } );

      // pull the response back to the dimension of the element
      ForceSized         force;
      SurfaceStressSized surfaceStress;
      QMatrixSized       Q_ik;
      HMatrixSized       H_ikR;
      HAverageJumpSized  H_iRk;
      ZMatrixSized       A_iRkP;

      Eigen::Matrix< double, nDim, nStations >    dForce_dNonlocal;
      Eigen::Matrix< double, nTensor, nStations > dSurface_dNonlocal;
      Eigen::Matrix< double, nStations, nDim >    dLocal_dJump;
      Eigen::Matrix< double, nStations, nTensor > dLocal_dAverage;

      if constexpr ( nDim == 3 ) {
        force         = mapEigenToFastor( state.force );
        surfaceStress = Eigen::Map< const SurfaceStressSized >( state.surfaceStress.data() );
        Q_ik          = mapEigenToFastor( tangents.Q_ik );
        H_ikR         = mapEigenToFastor( tangents.H_ikR );
        H_iRk         = Eigen::Map< const HAverageJumpSized >( tangents.H_iRk.data() );
        A_iRkP        = mapEigenToFastor( tangents.A_iRkP );

        dForce_dNonlocal   = tangents.dForce_dNonlocal;
        dSurface_dNonlocal = tangents.dSurfaceStress_dNonlocal;
        dLocal_dJump       = tangents.dLocal_dJump;
        dLocal_dAverage    = tangents.dLocal_dAverage;
      }
      else {
        const Fastor::Tensor< double, 2 >          force2d         = reduceTo2D< U >( state.force );
        const Fastor::Tensor< double, 2, 2 >       surfaceStress2d = reduceTo2D< U, U >( state.surfaceStress );
        const Fastor::Tensor< double, 2, 2 >       Q2d             = reduceTo2D< U, U >( tangents.Q_ik );
        const Fastor::Tensor< double, 2, 2, 2 >    H2d             = reduceTo2D< U, U, U >( tangents.H_ikR );
        const Fastor::Tensor< double, 2, 2, 2 >    HAverageJump2d  = reduceTo2D< U, U, U >( tangents.H_iRk );
        const Fastor::Tensor< double, 2, 2, 2, 2 > A2d             = reduceTo2D< U, U, U, U >( tangents.A_iRkP );

        force         = mapEigenToFastor( force2d );
        surfaceStress = Eigen::Map< const SurfaceStressSized >( surfaceStress2d.data() );
        Q_ik          = mapEigenToFastor( Q2d );
        H_ikR         = mapEigenToFastor( H2d );
        H_iRk         = Eigen::Map< const HAverageJumpSized >( HAverageJump2d.data() );
        A_iRkP        = mapEigenToFastor( A2d );

        // the 3D derivatives are indexed by (i, R) with i the first index of the tensor
        for ( int s = 0; s < nStations; ++s ) {
          for ( int i = 0; i < nDim; ++i ) {
            dForce_dNonlocal( i, s ) = tangents.dForce_dNonlocal( i, s );
            dLocal_dJump( s, i )     = tangents.dLocal_dJump( s, i );
            for ( int R = 0; R < nDim; ++R ) {
              dSurface_dNonlocal( nDim * i + R, s ) = tangents.dSurfaceStress_dNonlocal( 3 * i + R, s );
              dLocal_dAverage( s, nDim * i + R )    = tangents.dLocal_dAverage( s, 3 * i + R );
            }
          }
        }
      }

      qp.managedStateVars->force          = force;
      qp.managedStateVars->surfaceStress  = surfaceStress;
      qp.managedStateVars->displacement   = totalU_GPs;
      qp.managedStateVars->surfaceStrain  = totalSurfaceGradient_GPs;
      qp.managedStateVars->localFields    = state.localFields;
      qp.managedStateVars->nonlocalFields = nonlocal;

      // momentum
      Pe.template head< nDofD >() += ( Njump.transpose() * force + Bavg.transpose() * surfaceStress ) * qp.J0xW;

      Ke.template topLeftCorner< nDofD, nDofD >() += ( Njump.transpose() * Q_ik * Njump +
                                                       Njump.transpose() * H_ikR * Bavg +
                                                       Bavg.transpose() * H_iRk * Njump +
                                                       Bavg.transpose() * A_iRkP * Bavg ) *
                                                     qp.J0xW;

      // the stress depends on the nonlocal field of every station
      for ( int f = 0; f < nStations; ++f )
        Ke.template block< nDofD, nDofA >( 0,
                                           nDofD + f * nDofA ) += ( Njump.transpose() * dForce_dNonlocal.col( f ) +
                                                                    Bavg.transpose() * dSurface_dNonlocal.col( f ) ) *
                                                                  qp.N * qp.J0xW;

      // the nonlocal equation of every station, reduced to the midsurface, and the derivative through the thickness
      const double                                thickness = material->thickness();
      const Eigen::Matrix< double, nDofA, nDofA > NN        = qp.N.transpose() * qp.N;
      const Eigen::Matrix< double, nDofA, nDofA > gradGrad  = qp.gradN.transpose() * qp.gradN;

      for ( int s = 0; s < nStations; ++s ) {
        const double weight = thickness * material->stationWeight( s );
        const double l      = state.nonlocalRadii( s );
        const int    offS   = nDofD + s * nDofA;

        Pe.template segment< nDofA >( offS ) += weight *
                                                ( qp.N.transpose() * ( nonlocal( s ) - state.localFields( s ) ) +
                                                  l * l * qp.gradN.transpose() * ( qp.gradN * nodalNonlocal[s] ) ) *
                                                qp.J0xW;

        // the local field depends on the displacements and on the fields of all stations
        Ke.template block< nDofA, nDofD >( offS,
                                           0 ) -= weight * qp.N.transpose() *
                                                  ( dLocal_dJump.row( s ) * Njump + dLocal_dAverage.row( s ) * Bavg ) *
                                                  qp.J0xW;

        for ( int f = 0; f < nStations; ++f ) {
          const int offF = nDofD + f * nDofA;
          Ke.template block< nDofA, nDofA >( offS, offF ) -= weight * NN * tangents.dLocal_dNonlocal( s, f ) * qp.J0xW;
        }

        Ke.template block< nDofA, nDofA >( offS, offS ) += weight * ( NN + l * l * gradGrad ) * qp.J0xW;

        // the derivative through the thickness: the intervals on both sides of the station
        for ( const int neighbour : { s - 1, s + 1 } ) {
          if ( neighbour < 0 || neighbour >= nStations )
            continue;

          const double distance = std::abs( material->stationPosition( neighbour ) - material->stationPosition( s ) );
          const double coupling = l * l / ( thickness * distance );
          const int    offN     = nDofD + neighbour * nDofA;

          Pe.template segment< nDofA >( offS ) += coupling * qp.N.transpose() *
                                                  ( nonlocal( s ) - nonlocal( neighbour ) ) * qp.J0xW;
          Ke.template block< nDofA, nDofA >( offS, offS ) += coupling * NN * qp.J0xW;
          Ke.template block< nDofA, nDofA >( offS, offN ) -= coupling * NN * qp.J0xW;
        }
      }
    }
  };

  /**
   * @brief Gradient-enhanced interface element at finite strain with a nonlocal field at every station through the
   * thickness.
   *
   * @tparam nDim Spatial embedding dimension.
   * @tparam nNodes Number of element nodes.
   * @tparam nStations Number of stations, which is the number of nonlocal fields.
   */
  template < int nDim, int nNodes, int nStations >
  using GYSInterfaceFiniteElementFiniteStrain = InterfaceFiniteElementBase<
    nDim,
    nNodes,
    GradientEnhancedStationFieldsFormulation< nDim, nNodes, nStations > >;

} // namespace Marmot::Elements
