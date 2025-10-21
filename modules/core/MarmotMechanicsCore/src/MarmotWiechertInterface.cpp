#include "Marmot/MarmotWiechertInterface.h"
#include <iostream>
#include <ostream>

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Eigen;

  namespace WiechertInterface {

    // Properties generateRetardationTimes( int n, double min, double spacing )
    //{ //std::cout<<"Inside generateRetardationTimes"<<std::endl;
    //   Properties retardationTimes( n );
    //   for ( int i = 0; i < n; i++ )
    //     retardationTimes( i ) = min * std::pow( spacing, i );
    //   return retardationTimes;
    // }
    Properties initializeElasticModuliRu( int nMaxwellRu, double nRu )
    {
      Properties elasticModuliRu( nMaxwellRu );
      elasticModuliRu << nRu;
      return elasticModuliRu;
    }

    Properties initializeElasticModuliRs( int nMaxwellRs, double nRs )
    {
      Properties elasticModuliRs( nMaxwellRs );
      elasticModuliRs << nRs;
      return elasticModuliRs;
    }

    Properties initializeRelaxationTimesRu( int nMaxwellRu, double mRu )
    {
      Properties relaxationTimesRu( nMaxwellRu );
      relaxationTimesRu << mRu;
      return relaxationTimesRu;
    }

    Properties initializeRelaxationTimesRs( int nMaxwellRs, double mRs )
    {
      Properties relaxationTimesRs( nMaxwellRs );
      relaxationTimesRs << mRs;
      return relaxationTimesRs;
    }

    void evaluateWiechertRu( double           dT,
                             Properties       elasticModuliRu,
                             Properties       relaxationTimesRu,
                             StateVarMatrixRu stateVarsRu,
                             double&          uniaxialStiffnessRu,
                             Vector3d&        dforce,
                             const double     factor )
    {
      for ( int i = 0; i < relaxationTimesRu.size(); i++ ) {
        const double& tau = relaxationTimesRu( i );
        const double& D   = elasticModuliRu( i );
        double        lambda, beta;
        computeLambdaAndBeta( dT, tau, lambda, beta );
        uniaxialStiffnessRu += lambda * D * factor;
        dforce += ( 1. - beta ) * stateVarsRu.col( i ) * factor;
      }
    }

    void evaluateWiechertRs( double           dT,
                             Properties       elasticModuliRs,
                             Properties       relaxationTimesRs,
                             StateVarMatrixRs stateVarsRs,
                             double&          uniaxialStiffnessRs,
                             Vector9d&        dsurfaceStress,
                             const double     factor )
    {
      for ( int i = 0; i < relaxationTimesRs.size(); i++ ) {
        const double& tau = relaxationTimesRs( i );
        const double& D   = elasticModuliRs( i );

        double lambda, beta;
        computeLambdaAndBeta( dT, tau, lambda, beta );

        uniaxialStiffnessRs += lambda * D * factor;
        dsurfaceStress += ( 1. - beta ) * stateVarsRs.col( i ) * factor;
      }
    }

    void updateStateVarMatrixRu( double                  dT,
                                 Properties              elasticModuliRu,
                                 Properties              relaxationTimesRu,
                                 Ref< StateVarMatrixRu > stateVarsRu,
                                 const Vector3d&         dforce,
                                 const Matrix3d&         unitH_inv_ij )
    {

      if ( dT <= 1e-14 )
        return;
      for ( int i = 0; i < relaxationTimesRu.size(); i++ ) {
        const double& tau = relaxationTimesRu( i );
        const double& D   = elasticModuliRu( i );
        double        lambda, beta;
        computeLambdaAndBeta( dT, tau, lambda, beta );
        stateVarsRu.col( i ) = ( lambda * D ) * unitH_inv_ij * dforce + beta * stateVarsRu.col( i );
      }
    }

    void updateStateVarMatrixRs( double                  dT,
                                 Properties              elasticModuliRs,
                                 Properties              relaxationTimesRs,
                                 Ref< StateVarMatrixRs > stateVarsRs,
                                 const Vector9d&         dsurfaceStress,
                                 const Matrix9d&         unitZ_ijkl )
    {

      if ( dT <= 1e-14 )
        return;
      for ( int i = 0; i < relaxationTimesRs.size(); i++ ) {
        const double& tau = relaxationTimesRs( i );
        const double& D   = elasticModuliRs( i );
        double        lambda, beta;
        computeLambdaAndBeta( dT, tau, lambda, beta );
        stateVarsRs.col( i ) = ( lambda * D ) * unitZ_ijkl * dsurfaceStress + beta * stateVarsRs.col( i );
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

  } // namespace WiechertInterface
} // namespace Marmot::Materials
