#include "Marmot/FiniteStrainJ2Plasticity.h"
#include "Marmot/MarmotAutomaticDifferentiationForFastor.h"
#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotEnergyDensityFunctions.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotMicromorphicTensorBasics.h"
#include "Marmot/MarmotNumericalDifferentiationForFastor.h"
#include "Marmot/MarmotStressMeasures.h"
#include "Marmot/MarmotVoigt.h"
#include <Fastor/expressions/linalg_ops/unary_trans_op.h>
#include <Fastor/tensor/Tensor.h>
#include <Fastor/tensor_algebra/einsum_explicit.h>
#include <Fastor/tensor_algebra/indicial.h>
#include <autodiff/forward/dual/dual.hpp>
#include <map>

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Fastor;
  using namespace FastorIndices;
  using namespace FastorStandardTensors;

  FiniteStrainJ2Plasticity::FiniteStrainJ2Plasticity( const double* materialProperties,
                                                      int           nMaterialProperties,
                                                      int           materialLabel )
    : MarmotMaterialFiniteStrain( materialProperties, nMaterialProperties, materialLabel ),
      K( materialProperties[0] ),
      G( materialProperties[1] ),
      fy( materialProperties[2] ),
      fyInf( materialProperties[3] ),
      eta( materialProperties[4] ),
      H( materialProperties[5] ),
      implementationType( materialProperties[6] )
  {
  }

  void FiniteStrainJ2Plasticity::computeStress( ConstitutiveResponse< 3 >& response,
                                                AlgorithmicModuli< 3 >&    tangents,
                                                const Deformation< 3 >&    deformation,
                                                const TimeIncrement&       timeIncrement )
  {
    switch ( implementationType ) {
    case 0: computeStressWithScalarReturnMapping( response, tangents, deformation, timeIncrement ); break;
    case 1: computeStressWithFullReturnMapping( response, tangents, deformation, timeIncrement ); break;
    default: throw std::invalid_argument( "implementation type not supported" );
    };
  }
  void FiniteStrainJ2Plasticity::computeStressWithScalarReturnMapping( ConstitutiveResponse< 3 >& response,
                                                                       AlgorithmicModuli< 3 >&    tangents,
                                                                       const Deformation< 3 >&    deformation,
                                                                       const TimeIncrement&       timeIncrement )
  {
    throw std::invalid_argument( "not implemented yet" );
  }
  void FiniteStrainJ2Plasticity::computeStressWithFullReturnMapping( ConstitutiveResponse< 3 >& response,
                                                                     AlgorithmicModuli< 3 >&    tangents,
                                                                     const Deformation< 3 >&    deformation,
                                                                     const TimeIncrement&       timeIncrement )
  {

    auto&           Fp = stateVars->Fp;
    const Tensor33d FpOld( Fp );
    double&         alphaP    = stateVars->alphaP;
    const double    alphaPOld = alphaP;

    using namespace Marmot;
    using namespace Fastor;
    using namespace Eigen;
    using namespace autodiff;
    using namespace FastorIndices;
    using namespace FastorStandardTensors;

    Tensor33d FeTrial = deformation.F % Fastor::inverse( FpOld );
    double    betaP, dBetaP_dAlphaP;
    std::tie( betaP, dBetaP_dAlphaP ) = computeBetaP( alphaPOld );

    Tensor33d dFp = Spatial3D::I;
    Tensor33d Fe  = FeTrial;

    if ( isYielding( FeTrial, betaP ) ) {

      size_t counter = 0;

      using mV9d = Eigen::Map< Eigen::Matrix< double, 9, 1 > >;
      VectorXd X( 11 );
      X.segment( 0, 9 )  = mV9d( FeTrial.data() );
      X( 9 )             = alphaP;
      X( 10 )            = 0;
      VectorXd        dX = VectorXd::Zero( 11 );
      VectorXd        R  = VectorXd::Zero( 11 );
      Eigen::MatrixXd dR_dX( 11, 11 );

      std::tie( R, dR_dX ) = computeResidualVector( X, FeTrial, alphaP );

      while ( R.norm() > 1e-12 || dX.norm() > 1e-12 ) {
        /* std::cout << "residual: " << R.norm() << std::endl; */
        /* std::cout << "R: " << std::endl << R.norm() << std::endl; */
        /* std::cout << "dR_dX: " << std::endl << dR_dX << std::endl; */
        if ( counter > 10 )
          throw std::runtime_error( "inner newton not converged" );

        dX = -dR_dX.colPivHouseholderQr().solve( R );
        X += dX;
        std::tie( R, dR_dX ) = computeResidualVector( X, FeTrial, alphaP );
        counter += 1;
      }

      // update plastic deformation increment
      Fe              = X.segment( 0, 9 ).data();
      dFp             = Fastor::inverse( Fe ) % FeTrial;
      alphaP          = X( 9 );
      Tensor33d FpNew = dFp % Fp;
      memcpy( Fp.data(), FpNew.data(), 9 * sizeof( double ) );

      using namespace ContinuumMechanics;
      double      psi_;
      Tensor33d   Ce, dPsi_dCe;
      Tensor3333d dCe_dFe, d2Psi_dCedCe;
      std::tie( Ce, dCe_dFe ) = DeformationMeasures::FirstOrderDerived::CauchyGreen( Fe );

      // compute energy density, first and second partial derivatives wrt Cauchy Green deformation
      std::tie( psi_, dPsi_dCe, d2Psi_dCedCe ) = EnergyDensityFunctions::SecondOrderDerived::PenceGouPotentialB( Ce,
                                                                                                                 K,
                                                                                                                 G );

      // compute Kirchhoff stress
      Tensor33d   PK2 = 2. * dPsi_dCe;
      Tensor3333d dTau_dPK2, dTau_dFe_;
      std::tie( response.tau, dTau_dPK2, dTau_dFe_ ) = StressMeasures::FirstOrderDerived::KirchhoffStressFromPK2( PK2,
                                                                                                                  Fe );
      response.rho                                   = 1.0;
      response.elasticEnergyDensity                  = psi_;

      // compute tangent operator
      Tensor3333d dTau_dFe = einsum< ijKL, KLMN >( einsum< ijKL, IJKL >( dTau_dPK2, 2.0 * d2Psi_dCedCe ), dCe_dFe ) +
                             dTau_dFe_;
      using mM9d = Eigen::Map< Eigen::Matrix< double, 9, 9 > >;

      MatrixXd dYdDeformation              = MatrixXd::Zero( 11, 11 );
      dYdDeformation.block< 9, 9 >( 0, 0 ) = mM9d( Tensor3333d( einsum< KI, JL, to_IJKL >( Fastor::inverse( FpOld ),
                                                                                           Spatial3D::I ) )
                                                     .data() )
                                               .transpose();
      MatrixXd dXdDeformation = dR_dX.colPivHouseholderQr().solve( dYdDeformation );

      Tensor3333d dFe_dF = Tensor3333d( Matrix9d( dXdDeformation.block< 9, 9 >( 0, 0 ).transpose() ).data() );

      tangents.dTau_dF = einsum< IJKL, KLMN >( dTau_dFe, dFe_dF );
    }
    else {
      using namespace ContinuumMechanics;
      double      psi_;
      Tensor33d   C, dPsi_dC;
      Tensor3333d dC_dF, d2Psi_dCdC;
      std::tie( C, dC_dF ) = DeformationMeasures::FirstOrderDerived::CauchyGreen( Fe );

      // compute energy density, first and second partial derivatives wrt Cauchy Green deformation
      std::tie( psi_, dPsi_dC, d2Psi_dCdC ) = EnergyDensityFunctions::SecondOrderDerived::PenceGouPotentialB( C, K, G );

      // compute Kirchhoff stress
      Tensor33d   PK2 = 2. * dPsi_dC;
      Tensor3333d dTau_dPK2, dTau_dFe_;
      std::tie( response.tau, dTau_dPK2, dTau_dFe_ ) = StressMeasures::FirstOrderDerived::KirchhoffStressFromPK2( PK2,
                                                                                                                  Fe );
      response.rho                                   = 1.0;
      response.elasticEnergyDensity                  = psi_;

      // compute tangent operator
      Tensor3333d dTau_dFe = einsum< ijKL, KLMN >( einsum< ijKL, KLMN >( dTau_dPK2, 2.0 * d2Psi_dCdC ), dC_dF ) +
                             dTau_dFe_;

      tangents.dTau_dF = dTau_dFe;
    }
  }

  StateView FiniteStrainJ2Plasticity::getStateView( const std::string& stateName )
  {
    return stateVars->getStateView( stateName );
  }

  void FiniteStrainJ2Plasticity::assignStateVars( double* stateVars_, int nStateVars )
  {
    if ( nStateVars < getNumberOfRequiredStateVars() )
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": Not sufficient stateVars!" );

    this->stateVars = std::make_unique< FiniteStrainJ2PlasticityStateVarManager >( stateVars_ );
  }

  void FiniteStrainJ2Plasticity::initializeYourself()
  {
    stateVars->Fp.eye();
  }
} // namespace Marmot::Materials
