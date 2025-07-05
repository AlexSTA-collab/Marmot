#include "Marmot/LinearViscoElasticInterface.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotMaterialHypoElasticInterface.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotKelvinChainInterface.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotViscoelasticity.h"
#include "Marmot/MarmotUtility.h"
#include "Marmot/MarmotVoigt.h"

#include "Fastor/Fastor.h"
#include "Marmot/interface_material_helper_functions.h"
#include <Eigen/src/Core/util/Constants.h>
#include <Fastor/tensor/TensorMap.h>

#include "autodiff/forward/real.hpp"
#include <iostream>
#include <map>
#include <string>

using namespace Marmot;
using namespace Eigen;
using namespace Fastor;

using Tensor1D = Fastor::Tensor<double,3>;
using Tensor2D = Fastor::Tensor<double,3,3>;
using Tensor3D = Fastor::Tensor<double,3,3,3>;
using Tensor4D = Fastor::Tensor<double,3,3,3,3 >;

namespace Marmot::Materials {


  
  LinearViscoElasticInterface::LinearViscoElasticInterface( const double* materialProperties,
                                                            int           nMaterialProperties, 
                                                            int                materialNumber )
    : MarmotMaterialHypoElasticInterface( materialProperties, nMaterialProperties, materialNumber ),
    // clang-format off 
    // elasticity parameters
    E_M                     (materialProperties[0]),
    nu_M                    (materialProperties[1]),
    E_I                     (materialProperties[2]),
    nu_I                    (materialProperties[3]),
    E_0                     (materialProperties[4]),
    nu_0                    (materialProperties[5]),
    h                       (materialProperties[6]),

    // equivalent viscoelastic parameters common for Jump u
    m_Ju                      (materialProperties[7]),
    n_Ju                      (materialProperties[8]),
    nKelvin_Ju               (static_cast<size_t>(materialProperties[9])),
    minTau_Ju                 (materialProperties[10]),
    
    // equivalent viscoelastic parameters common for surface strain
    m_Js                      (materialProperties[11]),
    n_Js                      (materialProperties[12]),
    nKelvin_Js                (static_cast<size_t>(materialProperties[13])),
    minTau_Js                 (materialProperties[14]),
    

    timeToDays              (materialProperties[15])
  // clang-format on
  {

    retardationTimes_Ju = KelvinChainInterface::generateRetardationTimes(nKelvin_Ju, minTau_Ju, sqrt(10));
    retardationTimes_Js = KelvinChainInterface::generateRetardationTimes(nKelvin_Js, minTau_Js, sqrt(10));

    using namespace Marmot::ContinuumMechanics::Viscoelasticity ;
    auto phiJu_ = [&](autodiff::Real<powerLawApproximationOrder, double> tau){
        return ComplianceFunctions::powerLaw( tau, m_Ju, n_Ju);
    };

     auto phiJs_ = [&](autodiff::Real<powerLawApproximationOrder, double> tau){
        return ComplianceFunctions::powerLaw( tau, m_Js, n_Js);
    };


    elasticModuli_Ju = KelvinChainInterface::computeElasticModuli_Ju<powerLawApproximationOrder>(phiJu_, retardationTimes_Ju);
    elasticModuli_Js = KelvinChainInterface::computeElasticModuli_Js<powerLawApproximationOrder>(phiJs_, retardationTimes_Js);

    zerothKelvinChainCompliance_Ju = m_Ju*(1. - n_Ju )*pow( 2., n_Ju )*pow(minTau_Ju/sqrt(10.), n_Ju);  
    zerothKelvinChainCompliance_Js = m_Js*(1. - n_Js )*pow( 2., n_Js )*pow(minTau_Js/sqrt(10.), n_Js);
  }

  void LinearViscoElasticInterface::computeStress( double*  force,
                                                   double*  surface_stress,
                                                   double* dStress_dStrain,
                                                   const double* dU,
                                                   const double* dSurface_strain,
                                                   const double* normal,
                                                   const double* timeOld,
                                                   const double  dT,
                                                   double&       pNewDT)
  {
    // map to force, surface stress, displacement, surface strain, normal and tangent stiffness 
    // use Fastor because we really need to use the einsum  
    Fastor::Tensor<double, 3>  force_ftensor(force);
    Fastor::Tensor<double,3,3>  surface_stress_ftensor(surface_stress);
    Fastor::Tensor<double,21,21> dStress_dStrain_ftensor(dStress_dStrain);
    auto dU_ftensor_const = Fastor::TensorMap<const double,6,1>(dU);
    auto dSurface_strain_ftensor_const = Fastor::TensorMap<const double,18,1>(dSurface_strain);
    auto normal_ftensor_const = Fastor::TensorMap<const double, 3>(normal); 
    
    Fastor::Tensor<double,6,1> dU_ftensor(dU_ftensor_const.data());
    Fastor::Tensor<double,18,1> dSurface_strain_ftensor(dSurface_strain_ftensor_const.data());
    Fastor::Tensor<double,3> normal_ftensor(normal_ftensor_const.data());

    auto [Z_ijkl, H_inv_ij, H_inv_nF_ijk, Yn_H_inv_Fn_ijkl] = calculate_interface_material_parameters(normal_ftensor, E_M, nu_M, E_I, nu_I, E_0, nu_0); 
    
    //Assign the material matrices to a larger structure. (Not necessary ...)
    Eigen::Matrix<double, 21, 21> Cel = Eigen::Matrix<double, 21,21>::Zero();


    // handle zero strain increment
    if ( Fastor::norm(dU_ftensor) <1e-14 && Fastor::norm(dSurface_strain_ftensor) <1e-14 && dT==0  ) {
        Cel.block(0,0,9,9) = convert4thOrderTensorToMatrix(h/2.*Z_ijkl);
        Cel.block(12,12,9,9) = convert4thOrderTensorToMatrix(h/2.*Yn_H_inv_Fn_ijkl);
        Cel.block(0,9,9,3) = convert3rdOrderTensorToMatrix(H_inv_nF_ijk);
        Cel.block(9,9,3,3) = convert2ndOrderTensorToMatrix(2./h*H_inv_ij);
    
        Fastor::Tensor<double,21,21> Cel_fastor(Cel.data()); 
        dStress_dStrain_ftensor = Cel_fastor;
        std::copy(dStress_dStrain_ftensor.data(), dStress_dStrain_ftensor.data() + 21*21, dStress_dStrain);
      return;
    }
    //visco elastic step
  
    Eigen::Ref< KelvinChainInterface::mapStateVarMatrix_Ju> creepStateVars_Ju( stateVarManager->kelvinStateVars_Ju );
    Eigen::Ref< KelvinChainInterface::mapStateVarMatrix_Js> creepStateVars_Js( stateVarManager->kelvinStateVars_Js );
     
    const double dTimeDays = dT * timeToDays;
    
    Vector3d creep_Ju_Increment = Vector3d::Zero();
    double creep_Ju_compliance = 0;

    KelvinChainInterface::evaluateKelvinChain_Ju(dTimeDays,
                                                 elasticModuli_Ju,
                                                 retardationTimes_Ju,
                                                 creepStateVars_Ju,
                                                 creep_Ju_compliance,
                                                 creep_Ju_Increment,
                                                 1.0);
    
    Vector9d creep_Js_Increment = Vector9d::Zero();
    double creep_Js_compliance = 0;
 
    KelvinChainInterface::evaluateKelvinChain_Js(dTimeDays,
                                                 elasticModuli_Js,
                                                 retardationTimes_Js,
                                                 creepStateVars_Js,
                                                 creep_Js_compliance,
                                                 creep_Js_Increment,
                                                 1.0);
    
    using namespace Marmot::ContinuumMechanics::Viscoelasticity;

    //Evaluate effective compliances due to the displacement jump and the surface stress

    
    auto [unitH_ij, unitZ_inv_ijkl, E_eff, H_eff, nu_eff] = calculate_unitcompliance_interface(normal_ftensor, E_M, nu_M, E_I, nu_I, E_0, nu_0);

    double barC_Ju = std::pow(1./H_eff,-1) + zerothKelvinChainCompliance_Ju + creep_Ju_compliance;
    double barC_Js = 1./E_eff + zerothKelvinChainCompliance_Js + creep_Js_compliance;

    auto [H_inv_ij_effective, Z_ijkl_effective] = calculate_effective_properties(barC_Ju, barC_Js, normal_ftensor, nu_eff);

    enum {i,j,k,l};
    // Calculate jump increment
    Tensor1D jumpU_ftensor = dU_ftensor(Fastor::seq(0,3),0)-dU_ftensor(Fastor::seq(3,Fastor::last),0);    
    
    // Calculate average surface strain increment
    Fastor::Tensor<double, 9,1> average_dSurface_strain_ftensor = 1./2.*(dSurface_strain_ftensor(Fastor::seq(0,9),0)+dSurface_strain_ftensor(Fastor::seq(9,Fastor::last),0));
    auto average_dSurface_strain_ftensor_reshape = Fastor::reshape<3,3>(average_dSurface_strain_ftensor);

    Fastor::TensorMap<double,3> creep_Ju_increment_fastor(creep_Ju_Increment.data());
    Tensor1D creep_Ju_increment_fastor_tensor = creep_Ju_increment_fastor;
    Tensor1D jU_increment_difference_fastor = jumpU_ftensor- creep_Ju_increment_fastor_tensor; 
    Tensor1D  deltaforce_i = Fastor::einsum<Fastor::Index<i,j>, Fastor::Index<j>, Fastor::OIndex<i>>(H_inv_ij_effective, jU_increment_difference_fastor);

    Fastor::TensorMap<double,3,3> creep_Js_increment_fastor(creep_Js_Increment.data());    
    Tensor2D creep_Js_increment_fastor_tensor = creep_Js_increment_fastor;
    Tensor2D jdSurface_strain_increment_difference_fastor = average_dSurface_strain_ftensor_reshape - creep_Js_increment_fastor_tensor; 
    Tensor2D deltasurface_stress_ij = Fastor::einsum<Fastor::Index<i,j,k,l>,Fastor::Index<k,l>, Fastor::OIndex<i,j>>(Z_ijkl_effective, jdSurface_strain_increment_difference_fastor); 
    
    force_ftensor     += 2./h * deltaforce_i;
    surface_stress_ftensor += h/2.*deltasurface_stress_ij;
    
    // We write the variational form using the averaged incremental matrices
    Cel.block(0,0,9,9) = convert4thOrderTensorToMatrix(h/2.*Z_ijkl_effective);
    Cel.block(12,12,9,9) = convert4thOrderTensorToMatrix(h/2.*Yn_H_inv_Fn_ijkl);
    Cel.block(0,9,9,3) = convert3rdOrderTensorToMatrix(H_inv_nF_ijk);
    Cel.block(9,9,3,3) = convert2ndOrderTensorToMatrix(2./h*H_inv_ij_effective);

    Fastor::Tensor<double,21,21> Cel_fastor(Cel.data()); 
    dStress_dStrain_ftensor = Cel_fastor;
    
    std::copy(force_ftensor.data(), force_ftensor.data() + 3, force);
    std::copy(surface_stress_ftensor.data(), surface_stress_ftensor.data() + 3*3, surface_stress);
    std::copy(dStress_dStrain_ftensor.data(), dStress_dStrain_ftensor.data() + 21*21, dStress_dStrain);
    
    // Use already available functionality convert Fastor tensors to Eigen matricfes/vectors

    // Tranform to Eigen matrices to work with the internal machinery of KelvinChainInterface ...
    Eigen::Map<Eigen::Matrix<double,3,3, Eigen::RowMajor>>unitH_voigt_full(unitH_ij.data());
    Eigen::Map<Eigen::Matrix<double,9,9, Eigen::RowMajor>>unitZ_inv_voigt_full(unitZ_inv_ijkl.data());
    // Tranform to Eigen vectors to work with the internal machinery of KelvinChainInterface ...
    Eigen::Map<Eigen::Matrix<double,3, Eigen::RowMajor>> deltaforce_voigt_full(deltaforce_i.data());
    Eigen::Map<Eigen::Matrix<double,9, Eigen::RowMajor>> deltasurface_stress_voigt_full(deltasurface_stress_ij.data());
    
    KelvinChainInterface::updateStateVarMatrix_Ju( dTimeDays,
                                                   elasticModuli_Ju,
                                                   retardationTimes_Ju,
                                                   creepStateVars_Ju,
                                                   deltaforce_voigt_full,
                                                   unitH_voigt_full);
    
    KelvinChainInterface::updateStateVarMatrix_Js( dTimeDays,
                                                   elasticModuli_Js,
                                                   retardationTimes_Js,
                                                   creepStateVars_Js,
                                                   deltasurface_stress_voigt_full,
                                                   unitZ_inv_voigt_full);
  
    
    return;
      };
  
  void LinearViscoElasticInterface::assignStateVars( double* stateVars_, int nStateVars )
  {
    if ( nStateVars < getNumberOfRequiredStateVars() )
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": Not sufficient stateVars!" );

    this->stateVarManager = std::make_unique< LinearViscoElasticInterfaceStateVarManager >( stateVars_, nKelvin_Ju, nKelvin_Js );

    MarmotMaterial::assignStateVars( stateVars_, nStateVars );
  }

  StateView LinearViscoElasticInterface::getStateView( const std::string& stateName )
  {
    return stateVarManager->getStateView( stateName );
  }

  int LinearViscoElasticInterface::getNumberOfRequiredStateVars()
  {
    return LinearViscoElasticInterfaceStateVarManager::layout.nRequiredStateVars + nKelvin_Ju * 3 + nKelvin_Js * 9;
  }    
} // namespace Marmot::Materials
