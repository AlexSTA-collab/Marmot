#include "Marmot/LinearViscoElasticInterface.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialHypoElasticInterface.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotUtility.h"
#include "Marmot/MarmotViscoelasticity.h"
#include "Marmot/MarmotVoigt.h"
#include "Marmot/MarmotWiechertInterface.h"

#include "Fastor/Fastor.h"
#include "Marmot/MarmotInterfaceMaterialHelperFunctions.h"
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <Fastor/expressions/linalg_ops/unary_norm_op.h>
#include <Fastor/tensor/TensorMap.h>

#include "autodiff/forward/real.hpp"
#include <iostream>
#include <map>
#include <string>

using namespace Marmot;
using namespace Eigen;
using namespace Fastor;

using Tensor1D = Fastor::Tensor< double, 3 >;
using Tensor2D = Fastor::Tensor< double, 3, 3 >;
using Tensor3D = Fastor::Tensor< double, 3, 3, 3 >;
using Tensor4D = Fastor::Tensor< double, 3, 3, 3, 3 >;

namespace Marmot::Materials {

  LinearViscoElasticInterface::LinearViscoElasticInterface( const double* materialProperties,
                                                            int           nMaterialProperties,
                                                            int           materialNumber )
    : MarmotMaterialHypoElasticInterface( materialProperties, nMaterialProperties, materialNumber ),
      // clang-format off
      // elasticity parameters
      E_M( materialProperties[0] ),
      nu_M( materialProperties[1] ),
      E_I( materialProperties[2] ),
      nu_I( materialProperties[3] ),
      E_0( materialProperties[4] ),
      nu_0( materialProperties[5] ),
      h( materialProperties[6] ),

      // equivalent viscoelastic parameters common for Jump u
      mRu( materialProperties[7] ),
      nRu( materialProperties[8] ),
      nMaxwellRu( static_cast< size_t >( materialProperties[9] ) ),
      minTauRu( materialProperties[10] ),

      // equivalent viscoelastic parameters common for surface strain
      mRs( materialProperties[11] ),
      nRs( materialProperties[12] ),
      nMaxwellRs( static_cast< size_t >( materialProperties[13] ) ),
      minTauRs( materialProperties[14] ),

      timeToDays( materialProperties[15] )
  // clang-format on
  {

    // retardationTimes_Ru = Marmot::Materials::WiechertInterface::generateRetardationTimes(nMaxwell_Ru, minTau_Ru,
    // sqrt(10)); retardationTimes_Rs = Marmot::Materials::WiechertInterface::generateRetardationTimes(nMaxwell_Rs,
    // minTau_Rs, sqrt(10));
    relaxationTimesRu = Marmot::Materials::WiechertInterface::initializeRelaxationTimesRu( nMaxwellRu, mRu );
    relaxationTimesRs = Marmot::Materials::WiechertInterface::initializeRelaxationTimesRs( nMaxwellRs, mRs );
    elasticModuliRu   = Marmot::Materials::WiechertInterface::initializeElasticModuliRu( nMaxwellRu, nRu );
    elasticModuliRs   = Marmot::Materials::WiechertInterface::initializeElasticModuliRs( nMaxwellRs, nRs );

    using namespace Marmot::ContinuumMechanics::Viscoelasticity;
    auto phiRu_ = [&]( autodiff::Real< powerLawApproximationOrder, double > tau ) {
      return ComplianceFunctions::powerLaw( tau, mRu, nRu );
    };

    auto phiRs_ = [&]( autodiff::Real< powerLawApproximationOrder, double > tau ) {
      return ComplianceFunctions::powerLaw( tau, mRs, nRs );
    };

    // elasticModuli_Ru =
    // Marmot::Materials::WiechertInterface::computeElasticModuli_Ru<powerLawApproximationOrder>(phiRu_,
    // retardationTimes_Ru); elasticModuli_Rs =
    // Marmot::Materials::WiechertInterface::computeElasticModuli_Rs<powerLawApproximationOrder>(phiRs_,
    // retardationTimes_Rs);

    zerothWiechertStiffnessRu = 1.0; // m_Ru*(1. - n_Ru )*pow( 2., n_Ru )*pow(minTau_Ru/sqrt(10.), n_Ru);
    zerothWiechertStiffnessRs = 1.0; // m_Rs*(1. - n_Rs )*pow( 2., n_Rs )*pow(minTau_Rs/sqrt(10.), n_Rs);
  }

  void LinearViscoElasticInterface::computeStress( double*       force,
                                                   double*       surfaceStress,
                                                   double*       H_inv_ij,
                                                   double*       Z_ijkl,
                                                   double*       H_inv_nF_ijk,
                                                   double*       Yn_H_inv_Fn_ijkl,
                                                   const double* dU,
                                                   const double* dSurfaceStrain,
                                                   const double* normal,
                                                   const double* timeOld,
                                                   const double  dT,
                                                   double&       pNewDT )
  {
    using namespace Marmot::Materials::InterfaceMaterialHelperFunctions;

    // map to force, surface stress, displacement, surface strain, normal and tangent stiffness
    // use Fastor because we really need to use the einsum
    // std::cout<<"Inside proper file\n";

    Fastor::Tensor< double, 3 >          forceFtensor( force );
    Fastor::Tensor< double, 3, 3 >       surfaceStressFtensor( surfaceStress );
    Fastor::Tensor< double, 3, 3 >       H_inv_ij_Ftensor( H_inv_ij );
    Fastor::Tensor< double, 3, 3, 3, 3 > Z_ijkl_Ftensor( Z_ijkl );
    Fastor::Tensor< double, 3, 3, 3 >    H_inv_nF_ijk_Ftensor( H_inv_nF_ijk );
    Fastor::Tensor< double, 3, 3, 3, 3 > Yn_H_inv_Fn_ijkl_Ftensor( Yn_H_inv_Fn_ijkl );
    auto                                 dUFtensorConst = Fastor::TensorMap< const double, 6, 1 >( dU );
    auto dSurfaceStrainFtensorConst                     = Fastor::TensorMap< const double, 18, 1 >( dSurfaceStrain );
    auto normalFtensorConst                             = Fastor::TensorMap< const double, 3 >( normal );

    Fastor::Tensor< double, 6, 1 >  dUFtensor( dUFtensorConst.data() );
    Fastor::Tensor< double, 18, 1 > dSurfaceStrainFtensor( dSurfaceStrainFtensorConst.data() );
    Fastor::Tensor< double, 3 >     normalFtensor( normalFtensorConst.data() );

    auto [Z_ijkl_mat,
          H_inv_ij_mat,
          H_inv_nF_ijk_mat,
          Yn_H_inv_Fn_ijkl_mat] = calculateInterfaceMaterialParameters( normalFtensor,
                                                                        E_M,
                                                                        nu_M,
                                                                        E_I,
                                                                        nu_I,
                                                                        E_0,
                                                                        nu_0 );

    // Assign the material matrices to a larger structure. (Not necessary ...)
    Eigen::Matrix< double, 21, 21 > Cel = Eigen::Matrix< double, 21, 21 >::Zero();

    // handle zero strain increment
    if ( Fastor::norm( dUFtensor ) < 1e-14 && Fastor::norm( dSurfaceStrainFtensor ) < 1e-14 && dT == 0 ) {
      Z_ijkl_Ftensor           = h / 2. * Z_ijkl_mat;
      Yn_H_inv_Fn_ijkl_Ftensor = h / 2. * Yn_H_inv_Fn_ijkl_mat;
      H_inv_ij_Ftensor         = 2. / h * H_inv_ij_mat;
      H_inv_nF_ijk_Ftensor     = H_inv_nF_ijk_mat;
      return;
    }
    // visco elastic step

    Eigen::Ref< WiechertInterface::mapStateVarMatrixRu > creepStateVarsRu( stateVarManager->MaxwellStateVarsRu );
    Eigen::Ref< WiechertInterface::mapStateVarMatrixRs > creepStateVarsRs( stateVarManager->MaxwellStateVarsRs );

    const double dTimeDays = dT * timeToDays;

    Vector3d creepRuIncrement = Vector3d::Zero();
    double   creepRuStiffness = 0;

    WiechertInterface::evaluateWiechertRu( dTimeDays,
                                           elasticModuliRu,
                                           relaxationTimesRu,
                                           creepStateVarsRu,
                                           creepRuStiffness,
                                           creepRuIncrement,
                                           1.0 );

    Vector9d creepRsIncrement = Vector9d::Zero();
    double   creepRsStiffness = 0;

    WiechertInterface::evaluateWiechertRs( dTimeDays,
                                           elasticModuliRs,
                                           relaxationTimesRs,
                                           creepStateVarsRs,
                                           creepRsStiffness,
                                           creepRsIncrement,
                                           1.0 );

    using namespace Marmot::ContinuumMechanics::Viscoelasticity;

    // Evaluate effective compliances due to the displacement jump and the surface stress

    auto [H_inv_ij_effective,
          Z_ijkl_effective,
          unitH_inv_voigt_full,
          unitZ_voigt_full] = calculateEffectiveProperties( zerothWiechertStiffnessRu,
                                                            creepRuStiffness,
                                                            zerothWiechertStiffnessRs,
                                                            creepRsStiffness,
                                                            normalFtensor,
                                                            E_M,
                                                            nu_M,
                                                            E_I,
                                                            nu_I,
                                                            E_0,
                                                            nu_0 );

    enum { i, j, k, l };
    // Calculate jump increment
    Tensor1D jumpUFtensor = dUFtensor( Fastor::seq( 0, 3 ), 0 ) - dUFtensor( Fastor::seq( 3, Fastor::last ), 0 );

    // Calculate average surface strain increment
    Fastor::Tensor< double, 9, 1 > averageDsurfaceStrainFtensor = 1. / 2. *
                                                                  ( dSurfaceStrainFtensor( Fastor::seq( 0, 9 ), 0 ) +
                                                                    dSurfaceStrainFtensor( Fastor::seq( 9,
                                                                                                        Fastor::last ),
                                                                                           0 ) );
    auto averageDsurfaceStrainFtensorReshape = Fastor::reshape< 3, 3 >( averageDsurfaceStrainFtensor );

    Fastor::TensorMap< double, 3 > creepRuIncrementFastor( creepRuIncrement.data() );

    // std::cout<<"creep_Ru_increment_fastor:\n"<<creep_Ru_increment_fastor<<'\n';
    Tensor1D
      dForce_i = Fastor::einsum< Fastor::Index< i, j >, Fastor::Index< j >, Fastor::OIndex< i > >( H_inv_ij_effective,
                                                                                                   jumpUFtensor ) -
                 creepRuIncrementFastor;

    Fastor::TensorMap< double, 3, 3 > creepRsIncrementFastor( creepRsIncrement.data() );

    // std::cout<<"creep_Rs_increment_fastor:\n"<<creep_Rs_increment_fastor<<'\n';
    Tensor2D dSurfaceStress_ij = Fastor::einsum< Fastor::Index< i, j, k, l >,
                                                 Fastor::Index< k, l >,
                                                 Fastor::OIndex< i, j > >( Z_ijkl_effective,
                                                                           averageDsurfaceStrainFtensorReshape ) +
                                 creepRsIncrementFastor;

    forceFtensor += 2. / h * dForce_i;
    surfaceStressFtensor += h / 2. * dSurfaceStress_ij;

    Z_ijkl_Ftensor           = h / 2. * Z_ijkl_effective;
    Yn_H_inv_Fn_ijkl_Ftensor = h / 2. * Yn_H_inv_Fn_ijkl_mat;
    H_inv_ij_Ftensor         = 2. / h * H_inv_ij_effective;
    H_inv_nF_ijk_Ftensor     = H_inv_nF_ijk_mat;

    std::copy( forceFtensor.data(), forceFtensor.data() + 3, force );
    std::copy( surfaceStressFtensor.data(), surfaceStressFtensor.data() + 3 * 3, surfaceStress );

    // Use already available functionality convert Fastor tensors to Eigen matricfes/vectors
    // Tranform to Eigen matrices to work with the internal machinery of KelvinChainInterface ...
    Eigen::Map< Eigen::Matrix< double, 3, Eigen::RowMajor > > jumpUVoigtFull( jumpUFtensor.data() );
    Eigen::Map< Eigen::Matrix< double, 9, Eigen::RowMajor > > averageDsurfaceStrainFtensorReshapeVoigtFull(
      averageDsurfaceStrainFtensorReshape.data() );

    WiechertInterface::updateStateVarMatrixRu( dTimeDays,
                                               elasticModuliRu,
                                               relaxationTimesRu,
                                               creepStateVarsRu,
                                               jumpUVoigtFull,
                                               unitH_inv_voigt_full );

    WiechertInterface::updateStateVarMatrixRs( dTimeDays,
                                               elasticModuliRs,
                                               relaxationTimesRs,
                                               creepStateVarsRs,
                                               averageDsurfaceStrainFtensorReshapeVoigtFull,
                                               unitZ_voigt_full );

    return;
  };

  void LinearViscoElasticInterface::assignStateVars( double* stateVars_, int nStateVars )
  {
    if ( nStateVars < getNumberOfRequiredStateVars() )
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": Not sufficient stateVars!" );

    this->stateVarManager = std::make_unique< LinearViscoElasticInterfaceStateVarManager >( stateVars_,
                                                                                            nMaxwellRu,
                                                                                            nMaxwellRs );

    MarmotMaterial::assignStateVars( stateVars_, nStateVars );
  }

  StateView LinearViscoElasticInterface::getStateView( const std::string& stateName )
  {
    return stateVarManager->getStateView( stateName );
  }

  int LinearViscoElasticInterface::getNumberOfRequiredStateVars()
  {
    return LinearViscoElasticInterfaceStateVarManager::layout.nRequiredStateVars + 3 * nMaxwellRu + 9 * nMaxwellRs;
  }
} // namespace Marmot::Materials
