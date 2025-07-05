#include "Marmot/MarmotKelvinChainInterface.h"
#include <iostream>
#include <ostream>

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Eigen;

  namespace KelvinChainInterface {

    Properties generateRetardationTimes( int n, double min, double spacing )
    { //std::cout<<"Inside generateRetardationTimes"<<std::endl;
      Properties retardationTimes( n );
      for ( int i = 0; i < n; i++ )
        retardationTimes( i ) = min * std::pow( spacing, i );
      return retardationTimes;
    }

    void evaluateKelvinChain_Ju( double            dT,
                                 Properties        elasticModuli_Ju,
                                 Properties        retardationTimes_Ju,
                                 StateVarMatrix_Ju stateVars_Ju,
                                 double&           uniaxialCompliance_Ju,
                                 Vector3d&         dJumpU,
                                 const double      factor )
    {
      for ( int i = 0; i < retardationTimes_Ju.size(); i++ ) {
        const double& tau = retardationTimes_Ju( i );
        const double& D   = elasticModuli_Ju( i );

        double lambda, beta;
        computeLambdaAndBeta( dT, tau, lambda, beta );
        uniaxialCompliance_Ju += ( 1. - lambda ) / D * factor;
        dJumpU += ( 1. - beta ) * stateVars_Ju.col( i ) * factor;
      }
    }
    
    void evaluateKelvinChain_Js( double            dT,
                                 Properties        elasticModuli_Js,
                                 Properties        retardationTimes_Js,
                                 StateVarMatrix_Js stateVars_Js,
                                 double&           uniaxialCompliance_Js,
                                 Vector9d&         dsurface_strain,
                                 const double      factor )
    {
      for ( int i = 0; i < retardationTimes_Js.size(); i++ ) {
        const double& tau = retardationTimes_Js( i );
        const double& D   = elasticModuli_Js( i );

        double lambda, beta;
        computeLambdaAndBeta( dT, tau, lambda, beta );
        uniaxialCompliance_Js += ( 1. - lambda ) / D * factor;
        dsurface_strain += ( 1. - beta ) * stateVars_Js.col( i ) * factor;
      }
    }

    void updateStateVarMatrix_Ju( double                    dT,
                                  Properties                elasticModuli_Ju,
                                  Properties                retardationTimes_Ju,
                                  Ref< StateVarMatrix_Ju >  stateVars_Ju,
                                  const Vector3d&           dforce,
                                  const Matrix3d&           unitH_ij )
    {

      if ( dT <= 1e-14 )
        return;
      for ( int i = 0; i < retardationTimes_Ju.size(); i++ ) {
        const double& tau = retardationTimes_Ju( i );
        const double& D   = elasticModuli_Ju( i );
        double        lambda, beta;
        computeLambdaAndBeta( dT, tau, lambda, beta );
        stateVars_Ju.col( i ) = ( lambda / D ) * unitH_ij * dforce + beta * stateVars_Ju.col( i );
      }
    }

    void updateStateVarMatrix_Js( double                   dT,
                                  Properties               elasticModuli_Js,
                                  Properties               retardationTimes_Js,
                                  Ref< StateVarMatrix_Js > stateVars_Js,
                                  const Vector9d&          dsurface_stress,
                                  const Matrix9d&          unitZ_inv_ijkl )
    {

      if ( dT <= 1e-14 )
        return;
      for ( int i = 0; i < retardationTimes_Js.size(); i++ ) {
        const double& tau = retardationTimes_Js( i );
        const double& D   = elasticModuli_Js( i );
        double        lambda, beta;
        computeLambdaAndBeta( dT, tau, lambda, beta );
        stateVars_Js.col( i ) = ( lambda / D ) * unitZ_inv_ijkl * dsurface_stress + beta * stateVars_Js.col( i );
      }
    }

    void computeLambdaAndBeta( double dT, double tau, double& lambda, double& beta )
    {
      const double dT_tau = dT / tau;
      // respect extreme values according to Jirasek Bazant
      if ( dT_tau >= 30.0 ) {
        beta   = 0.;
        lambda = 1. / dT_tau;
      }
      else if ( dT_tau < 1e-6 ) {
        beta   = 1.0;
        lambda = 1 - 0.5 * dT_tau + 1. / 6 * dT_tau * dT_tau;
      }
      else {
        beta   = std::exp( -dT_tau );
        lambda = ( 1 - beta ) / dT_tau;
      }
    }

  } // namespace KelvinChainInterface
} // namespace Marmot::Materials
