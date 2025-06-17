#include "Marmot/MarmotKelvinChainInterface.h"

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Eigen;

  namespace KelvinChainInterface {

    Properties generateRetardationTimes( int n, double min, double spacing )
    {
      Properties retardationTimes( n );
      for ( int i = 0; i < n; i++ )
        retardationTimes( i ) = min * std::pow( spacing, i );
      return retardationTimes;
    }

    void evaluateKelvinChain( double         dT,
                              Properties     elasticModuli_Ju,
                              Properties     retardationTimes_Ju,
                              StateVarMatrix stateVars_Ju,
                              double&        uniaxialCompliance_Ju,
                              const Vector3d&     dforce,
                              const double   factor )
    {

      for ( int i = 0; i < retardationTimes_Ju.size(); i++ ) {
        const double& tau = retardationTimes_Ju( i );
        const double& D   = elasticModuli_Ju( i );

        double lambda, beta;
        computeLambdaAndBeta( dT, tau, lambda, beta );
        uniaxialCompliance_Ju += ( 1. - lambda ) / D * factor;
        dforce += ( 1. - beta ) * stateVars_Ju.col( i ) * factor;
      }
    }
    
    void evaluateKelvinChain( double         dT,
                              Properties     elasticModuli_Js,
                              Properties     retardationTimes_Js,
                              StateVarMatrix stateVars_Js,
                              double&        uniaxialCompliance_Js,
                              const Vector9d&     dsurface_stress,
                              const double   factor )
    {

      for ( int i = 0; i < retardationTimes_Js.size(); i++ ) {
        const double& tau = retardationTimes_Js( i );
        const double& D   = elasticModuli_Js( i );

        double lambda, beta;
        computeLambdaAndBeta( dT, tau, lambda, beta );
        uniaxialCompliance_Js += ( 1. - lambda ) / D * factor;
        dsurface_stress += ( 1. - beta ) * stateVars_Js.col( i ) * factor;
      }
    }

    void updateStateVarMatrix( double                dT,
                               Properties            elasticModuli_Ju,
                               Properties            retardationTimes_Ju,
                               Ref< StateVarMatrix > stateVars_Ju,
                               const Vector3d&       dforce,
                               const Matrix3d&       unitH_ij )
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

    void updateStateVarMatrix( double                dT,
                               Properties            elasticModuli_Js,
                               Properties            retardationTimes_Js,
                               Ref< StateVarMatrix > stateVars_Js,
                               const Vector9d&       dsurface_stress,
                               const Matrix9d&       unitZ_inv_ijkl )
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
