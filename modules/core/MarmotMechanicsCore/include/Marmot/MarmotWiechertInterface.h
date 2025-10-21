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
#include "Marmot/MarmotKelvinChain.h"
#include "Marmot/MarmotNumericalIntegration.h"
#include "Marmot/MarmotTypedefs.h"
#include "autodiff/forward/real.hpp"
#include <fstream>  // **needed for std::ofstream**
#include <functional>
#include <iostream> // for std::cout / std::cerr
//
namespace Marmot::Materials {

  namespace WiechertInterface {

    typedef Eigen::VectorXd          Properties;
    typedef Eigen::Map< Properties > mapProperties;

    typedef Eigen::Matrix< double, 3, Eigen::Dynamic > StateVarMatrixRu;
    typedef Eigen::Matrix< double, 9, Eigen::Dynamic > StateVarMatrixRs;

    typedef Eigen::Map< StateVarMatrixRu > mapStateVarMatrixRu;
    typedef Eigen::Map< StateVarMatrixRs > mapStateVarMatrixRs;

    // template < int k >
    // Properties computeElasticModuli_Ru( std::function< autodiff::Real< k, double >( autodiff::Real< k, double > ) >
    // phi,
    //                                  Properties retardationTimes_Ru,
    //                                  bool       gaussQuadrature = false )
    //{
    //   Properties elasticModuli_Ru( retardationTimes_Ru.size() );
    //   double     spacing = retardationTimes_Ru( 1 ) / retardationTimes_Ru( 0 );
    //   for ( int i = 0; i < retardationTimes_Ru.size(); i++ ) {
    //     double tau = retardationTimes_Ru( i );
    //     if ( !gaussQuadrature ) {
    //
    //       elasticModuli_Ru( i ) = 1. / ( log( spacing ) * KelvinChain::evaluatePostWidderFormula< k >( phi, tau ) );
    //
    //    }
    //    else {
    //      elasticModuli_Ru( i ) = 1. /
    //                           ( log( spacing ) / 2. *
    //                             ( KelvinChain::evaluatePostWidderFormula< k >( phi, tau * pow( spacing, -sqrt( 3. )
    //                             / 6. ) ) +
    //                               KelvinChain::evaluatePostWidderFormula< k >( phi, tau * pow( spacing, sqrt( 3. )
    //                               / 6. ) ) ) );
    //    }
    //
    //      return elasticModuli_Ru;
    //    }

    // template < int k >
    // Properties computeElasticModuli_Rs( std::function< autodiff::Real< k, double >( autodiff::Real< k, double > ) >
    // phi,
    //                                  Properties retardationTimes_Rs,
    //                                  bool       gaussQuadrature = false )
    //{
    //   Properties elasticModuli_Rs( retardationTimes_Rs.size() );
    //   double     spacing = retardationTimes_Rs( 1 ) / retardationTimes_Rs( 0 );
    //
    //      for ( int i = 0; i < retardationTimes_Rs.size(); i++ ) {
    //        double tau = retardationTimes_Rs( i );
    //        if ( !gaussQuadrature ) {
    //          elasticModuli_Rs( i ) = 1. / ( log( spacing ) * KelvinChain::evaluatePostWidderFormula< k >( phi, tau )
    //          );
    //       }
    //       else {
    //         elasticModuli_Rs( i ) = 1. /
    //                              ( log( spacing ) / 2. *
    //                               ( KelvinChain::evaluatePostWidderFormula< k >( phi, tau * pow( spacing, -sqrt( 3. )
    //                               / 6. ) ) +
    //                                 KelvinChain::evaluatePostWidderFormula< k >( phi, tau * pow( spacing, sqrt( 3. )
    //                                 / 6. ) ) ) );
    //       }
    //     }
    //
    //     return elasticModuli_Rs;
    //  }

    // Properties generateRelaxationTimes( int n, double min, double spacing );

    Properties initializeElasticModuliRu( int nMaxwellRu, double nRu );

    Properties initializeElasticModuliRs( int nMaxwellRs, double nRs );

    Properties initializeRelaxationTimesRu( int nMaxwellRu, double mRu );

    Properties initializeRelaxationTimesRs( int nMaxwellRs, double mRs );

    void updateStateVarMatrixRu( const double                   dT,
                                 Properties                     elasticModuli_Ru,
                                 Properties                     relaxationTimes_Ru,
                                 Eigen::Ref< StateVarMatrixRu > stateVarsRu,
                                 const Marmot::Vector3d&        dforce,
                                 const Marmot::Matrix3d&        unitH_inv_ij );

    void updateStateVarMatrixRs( const double                   dT,
                                 Properties                     elasticModuliRs,
                                 Properties                     relaxationTimesRs,
                                 Eigen::Ref< StateVarMatrixRs > stateVarsRs,
                                 const Marmot::Vector9d&        dsurface_stress,
                                 const Marmot::Matrix9d&        unitZ_ijkl );

    void evaluateWiechertRu( const double      dT,
                             Properties        elasticModuliRu,
                             Properties        relaxationTimesRu,
                             StateVarMatrixRu  stateVarsRu,
                             double&           uniaxialStiffnessRu,
                             Marmot::Vector3d& dforce_v,
                             const double      factor );

    void evaluateWiechertRs( const double      dT,
                             Properties        elasticModuliRs,
                             Properties        retardationTimesRs,
                             StateVarMatrixRs  stateVarsRs,
                             double&           uniaxialStiffnessRs,
                             Marmot::Vector9d& dsurface_stress,
                             const double      factor );

    void computeLambdaAndBeta( double dT, double tau, double& lambda, double& beta );

  } // namespace WiechertInterface
} // namespace Marmot::Materials
