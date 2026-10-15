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
 * @file GradientEnhancedStationFieldsFormulation.h
 * @brief Layout of a finite-strain interface element with a nonlocal field at every station through the thickness.
 *
 * This file defines the formulation `Marmot::Elements::GradientEnhancedStationFieldsFormulation`, which holds the
 * state-variable layout, the geometry of the quadrature points and the degrees of freedom of the nonlocal fields of
 * the elements GYH (GYHInterfaceFiniteElementFiniteStrain.h). The layer is divided by stations through its
 * thickness, and every station carries its own bulk material and its own nonlocal field, so that the damage varies
 * through the thickness. The residual and the tangent are assembled by the derived formulation.
 */
#pragma once

#include "Marmot/InterfaceFiniteElementBase.h"
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
   * @brief Base of the formulation of the interface element at finite strain with a nonlocal field at every station
   * through the thickness: layout of the state variables, geometry and degrees of freedom.
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
      initializeTwoSidedQuadraturePoint( element, qp, "GYHInterfaceFiniteElementFiniteStrain" );
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
  };

} // namespace Marmot::Elements
