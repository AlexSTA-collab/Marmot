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
 * @file GYHInterfaceFiniteElementFiniteStrain.h
 * @brief Two-sided gradient-enhanced interface finite element at finite strain with a nonlocal field at every
 * station through the thickness, for a layer with a stiffness varying through the thickness.
 *
 * This file defines the formulation `Marmot::Elements::GradientEnhancedTwoSidedStationFieldsFormulation` and the
 * element `Marmot::Elements::GYHInterfaceFiniteElementFiniteStrain`. It is the two-sided counterpart of
 * `GYSInterfaceFiniteElementFiniteStrain`: the surface gradients of the top and of the bottom face are given to the
 * material separately, and a surface stress is returned for every face, as for the equilibrated elements. The
 * nonlocal fields and their equations are those of the element with one-sided kinematics.
 */
#pragma once

#include "Marmot/GYInterfaceFiniteElementFiniteStrain.h"
#include "Marmot/GYSInterfaceFiniteElementFiniteStrain.h"
#include "Marmot/MarmotGradedEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain.h"

#include <Eigen/Dense>

namespace Marmot::Elements {

  /**
   * @struct GradientEnhancedTwoSidedStationFieldsFormulation
   * @tparam nDim Spatial embedding dimension.
   * @tparam nNodes Number of element nodes.
   * @tparam nStations Number of stations through the thickness, which is the number of nonlocal fields; at least two.
   * @brief Formulation of the two-sided interface element at finite strain with a nonlocal field at every station.
   *
   * The interface material (`MarmotGradedEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain`) receives the
   * total displacements and surface gradients of both interface sides, the connector between the faces, and the
   * nonlocal fields of all stations at the quadrature point. It returns the generalized force, conjugate to the
   * jump, the surface stresses of both faces, conjugate to their surface gradients, the local driving field and the
   * nonlocal radius of every station, and the condensed tangents. The momentum residual is assembled with the
   * operator of the generalized strain, and the nonlocal equations, the coupling through the thickness, and the
   * degrees of freedom are those of GradientEnhancedStationFieldsFormulation.
   */
  template < int nDim, int nNodes, int nStations >
  struct GradientEnhancedTwoSidedStationFieldsFormulation
    : GradientEnhancedStationFieldsFormulation< nDim, nNodes, nStations > {

    using Base = GradientEnhancedStationFieldsFormulation< nDim, nNodes, nStations >; ///< one-sided formulation

    static constexpr int nTensor      = Base::nTensor;      ///< components of a second order tensor
    static constexpr int nDofD        = Base::nDofD;        ///< degrees of freedom of the displacements
    static constexpr int nDofA        = Base::nDofA;        ///< degrees of freedom of the field of a station
    static constexpr int nGeneralized = nDim + 2 * nTensor; ///< components of ( w, A+, A- )

    using ForceSized                = typename Base::ForceSized;                ///< force
    using SurfaceStressSized        = typename Base::SurfaceStressSized;        ///< surface stress of a face
    using InterfaceDisplSized       = typename Base::InterfaceDisplSized;       ///< displacements of both sides
    using InterfaceSurfaceGradSized = typename Base::InterfaceSurfaceGradSized; ///< surface gradients of both sides

    using Material = MarmotGradedEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain; ///< interface material

    /**
     * @brief Named state-variable manager for interface quadrature points.
     */
    class QPStateVarManager : public MarmotStateVarVectorManager {

      /*
       * Persistent state layout for the generalized force, the two surface stresses, the total displacement, the
       * total surface strain, the local and the nonlocal field of every station, and the material state variables.
       */
      /// layout of the state variables of the quadrature point
      inline const static auto layout = makeLayout( {
        { .name = "generalizedForce", .length = nDim },
        { .name = "surfaceStressPlus", .length = nTensor },
        { .name = "surfaceStressMinus", .length = nTensor },
        { .name = "displacement", .length = 2 * nDim },
        { .name = "surface strain", .length = 2 * nTensor },
        { .name = "local fields", .length = nStations },
        { .name = "nonlocal fields", .length = nStations },
        // The material state starts at a multiple of four entries.
        { .name   = "state block alignment padding",
          .length = ( 4 - ( ( nDim + 2 * nTensor + 2 * nDim + 2 * nTensor + 2 * nStations ) % 4 ) ) % 4 },
        { .name = "begin of material state", .length = 0 },
      } );

    public:
      Eigen::Map< ForceSized >                            generalizedForce;   ///< force conjugate to the jump
      Eigen::Map< SurfaceStressSized >                    surfaceStressPlus;  ///< surface stress of the top face
      Eigen::Map< SurfaceStressSized >                    surfaceStressMinus; ///< surface stress of the bottom face
      Eigen::Map< InterfaceDisplSized >                   displacement;       ///< total displacements of the sides
      Eigen::Map< InterfaceSurfaceGradSized >             surfaceStrain;      ///< total surface gradients of the sides
      Eigen::Map< Eigen::Matrix< double, nStations, 1 > > localFields;        ///< local field of every station
      Eigen::Map< Eigen::Matrix< double, nStations, 1 > > nonlocalFields;     ///< local field of every station

      Eigen::Map< Eigen::VectorXd > materialStateVars; ///< state variables of the material

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
          displacement( &find( "displacement" ) ),
          surfaceStrain( &find( "surface strain" ) ),
          localFields( &find( "local fields" ) ),
          nonlocalFields( &find( "nonlocal fields" ) ),
          materialStateVars( &find( "begin of material state" ),
                             nStateVars - getNumberOfRequiredStateVarsQuadraturePointOnly() )
      {
      }
    };

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
          "GYHInterfaceFiniteElementFiniteStrain: the material has " + std::to_string( material->numberOfStations() ) +
          " stations, but the element has " + std::to_string( nStations ) + " nonlocal fields" );

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

      typename Material::State state{ typename Material::Tensor3d( 0.0 ),
                                      typename Material::Tensor33d( 0.0 ),
                                      typename Material::Tensor33d( 0.0 ),
                                      qp.managedStateVars->materialStateVars.data(),
                                      Eigen::VectorXd(),
                                      Eigen::VectorXd() };

      typename Material::Tangents tangents;

      typename Material::Deformation deformation{ typename Material::Tensor6d( U3d.data() ),
                                                  typename Material::Tensor18d( surfaceGradient3d.data() ),
                                                  embedVectorIn3D< nDim >( qp.normal.data() ),
                                                  embedVectorIn3D< nDim >( data.separationVector.data() ),
                                                  nonlocal };

      material->computeStress( state, tangents, deformation, { time, dT } );

      qp.managedStateVars->displacement   = totalU_GPs;
      qp.managedStateVars->surfaceStrain  = totalSurfaceGradient_GPs;
      qp.managedStateVars->localFields    = state.localFields;
      qp.managedStateVars->nonlocalFields = nonlocal;

      // momentum: the generalized stress and the condensed tangent of the displacements
      assembleTwoSidedResponse< nDim, nNodes >( qp,
                                                state,
                                                tangents,
                                                Pe.template head< nDofD >(),
                                                Ke.template topLeftCorner< nDofD, nDofD >() );

      // the couplings to the fields at the dimension of the element
      const auto reduce = [&]( const Eigen::Matrix< double, 21, 1 >& v ) {
        return reduceGeneralizedStrain< nDim >( Fastor::Tensor< double, 21 >( v.data() ) );
      };

      Eigen::Matrix< double, nGeneralized, nStations > dStress_dNonlocal;
      Eigen::Matrix< double, nStations, nGeneralized > dLocal_dStrain;
      for ( int s = 0; s < nStations; ++s ) {
        dStress_dNonlocal.col( s ) = reduce( tangents.dStress_dNonlocal.col( s ) );
        dLocal_dStrain.row( s )    = reduce( tangents.dLocal_dStrain.row( s ).transpose() ).transpose();
      }

      const auto Bx = generalizedStrainOperator< nDim, nNodes >( qp );

      // the stress depends on the nonlocal field of every station
      for ( int f = 0; f < nStations; ++f )
        Ke.template block< nDofD, nDofA >( 0, nDofD + f * nDofA ) += Bx.transpose() * dStress_dNonlocal.col( f ) *
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
        Ke.template block< nDofA, nDofD >( offS, 0 ) -= weight * qp.N.transpose() * dLocal_dStrain.row( s ) * Bx *
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
   * @brief Two-sided gradient-enhanced interface element at finite strain with a nonlocal field at every station
   * through the thickness.
   *
   * @tparam nDim Spatial embedding dimension.
   * @tparam nNodes Number of element nodes.
   * @tparam nStations Number of stations, which is the number of nonlocal fields.
   */
  template < int nDim, int nNodes, int nStations >
  using GYHInterfaceFiniteElementFiniteStrain = InterfaceFiniteElementBase<
    nDim,
    nNodes,
    GradientEnhancedTwoSidedStationFieldsFormulation< nDim, nNodes, nStations > >;

} // namespace Marmot::Elements
