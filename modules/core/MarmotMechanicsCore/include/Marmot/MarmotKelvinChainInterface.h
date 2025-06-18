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
#include "MarmotKelvinChain.h"
#include "Marmot/MarmotNumericalIntegration.h"
#include "Marmot/MarmotTypedefs.h"
#include "autodiff/forward/real.hpp"
#include <functional>

namespace Marmot::Materials {

  namespace KelvinChainInterface {

    typedef Eigen::VectorXd          Properties;
    typedef Eigen::Map< Properties > mapProperties;

    typedef Eigen::Matrix< double, 3, Eigen::Dynamic > StateVarMatrix_Ju;
    typedef Eigen::Matrix< double, 9, Eigen::Dynamic > StateVarMatrix_Js;
    
    typedef Eigen::Map< StateVarMatrix_Ju >               mapStateVarMatrix_Ju;
    typedef Eigen::Map< StateVarMatrix_Js >               mapStateVarMatrix_Js;
    
    using namespace KelvinChain;
    
    template < int k >
    Properties computeElasticModuli_Ju( std::function< autodiff::Real< k, double >( autodiff::Real< k, double > ) > phi,
                                     Properties retardationTimes_Ju,
                                     bool       gaussQuadrature = false )
    {
      Properties elasticModuli_Ju( retardationTimes_Ju.size() );
      double     spacing = retardationTimes_Ju( 1 ) / retardationTimes_Ju( 0 );

      for ( int i = 0; i < retardationTimes_Ju.size(); i++ ) {
        double tau = retardationTimes_Ju( i );
        if ( !gaussQuadrature ) {
          elasticModuli_Ju( i ) = 1. / ( log( spacing ) * evaluatePostWidderFormula< k >( phi, tau ) );
        }
        else {
          elasticModuli_Ju( i ) = 1. /
                               ( log( spacing ) / 2. *
                                 ( evaluatePostWidderFormula< k >( phi, tau * pow( spacing, -sqrt( 3. ) / 6. ) ) +
                                   evaluatePostWidderFormula< k >( phi, tau * pow( spacing, sqrt( 3. ) / 6. ) ) ) );
        }
      }

      return elasticModuli_Ju;
    }

    template < int k >
    Properties computeElasticModuli_Js( std::function< autodiff::Real< k, double >( autodiff::Real< k, double > ) > phi,
                                     Properties retardationTimes_Js,
                                     bool       gaussQuadrature = false )
    {
      Properties elasticModuli_Js( retardationTimes_Js.size() );
      double     spacing = retardationTimes_Js( 1 ) / retardationTimes_Js( 0 );

      for ( int i = 0; i < retardationTimes_Js.size(); i++ ) {
        double tau = retardationTimes_Js( i );
        if ( !gaussQuadrature ) {
          elasticModuli_Js( i ) = 1. / ( log( spacing ) * evaluatePostWidderFormula< k >( phi, tau ) );
        }
        else {
          elasticModuli_Js( i ) = 1. /
                               ( log( spacing ) / 2. *
                                 ( evaluatePostWidderFormula< k >( phi, tau * pow( spacing, -sqrt( 3. ) / 6. ) ) +
                                   evaluatePostWidderFormula< k >( phi, tau * pow( spacing, sqrt( 3. ) / 6. ) ) ) );
        }
      }

      return elasticModuli_Js;
    }

    Properties generateRetardationTimes( int n, double min, double spacing );

    void updateStateVarMatrix_Ju(    const double                    dT,
                                     Properties                      elasticModuli_Ju,
                                     Properties                      retardationTimes_Ju,
                                     Eigen::Ref< StateVarMatrix_Ju > stateVars_Ju,
                                     const Marmot::Vector3d&         dforce,
                                     const Marmot::Matrix3d&         unitH_ij );

    void updateStateVarMatrix_Js(    const double                 dT,
                                     Properties                   elasticModuli_Js,
                                     Properties                   retardationTimes_Js,
                                     Eigen::Ref< StateVarMatrix_Js > stateVars_Js,
                                     const Marmot::Vector9d&      dsurface_stress,
                                     const Marmot::Matrix9d&      unitZ_inv_ijkl );

    void evaluateKelvinChain_Ju(    const double         dT,
                                    Properties           elasticModuli_Ju,
                                    Properties           retardationTimes_Ju,
                                    StateVarMatrix_Ju    stateVars_Ju,
                                    double&              uniaxialCompliance_Ju,
                                    Marmot::Vector3d&    dJumpu,
                                    const double         factor );

    void evaluateKelvinChain_Js(    const double         dT,
                                    Properties           elasticModuli_Js,
                                    Properties           retardationTimes_Js,
                                    StateVarMatrix_Js    stateVars_Js,
                                    double&              uniaxialCompliance_Js,
                                    Marmot::Vector9d&    dsurface_strain,
                                    const double         factor );

    void computeLambdaAndBeta( double dT, double tau, double& lambda, double& beta );

  } // namespace KelvinChainInterface
} // namespace Marmot::Materials
