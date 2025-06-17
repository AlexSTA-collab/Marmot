#include "Marmot/LinearViscoElasticInterface.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotMaterialHypoElasticInterface.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotKelvinChain.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTypedefs.h"
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
                                                            int materialNumber )
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
    

    timeToDays              (materialProperties[15]),
  {
    retardationTimes_Ju = KelvinChainInterface::generateRetardationTimes(nKelvin_Ju, minTau_Ju, sqrt(10));
    retardationTimes_Js = KelvinChainInterface::generateRetardationTimes(nKelvin_Js, minTau_Js, sqer(10));

    using namespace Marmot::ContinuumMechanics::Viscoelasticity;
    auto phiJu_ = [&](autodiff::Real<powerLawApproximationOrder, double> tau){
        return ComplianceFunctions::powerLaw( tau, m_Ju, n_Ju);
    };

     auto phiJs_ = [&](autodiff::Real<powerLawApproximationOrder, double> tau){
        return ComplianceFunctions::powerLaw( tau, m_Js, n_Js);
    };


    elasticModuli_H_inv_ij = KelvinChainInterface::computeElasticModuli<powerLawApproximationOrder>(phiJu_, retardationTimes_Ju);
    elasticModuli_Zijkl = KelvinChainInterface::computeElasticModuli<powerLawApproximationOrder>(phiJs_, retardationTimes_Js);

    zerothKelvinCompliance_Ju = m_Ju*(1. - n_Ju )*pow( 2., n_Ju )*pow(minTau_Ju/sqrt(10.), n_Ju);  
    zerothKelvinCompliance_Js = m_Js*(1. - n_Js )*pow( 2., n_Js )*pow(minTau_Js/sqrt(10.), n_Js);
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

    Cel.block(0,0,9,9) = convert4thOrderTensorToMatrix(h/2.*Z_ijkl);
    Cel.block(12,12,9,9) = convert4thOrderTensorToMatrix(h/2.*Yn_H_inv_Fn_ijkl);
    Cel.block(0,9,9,3) = convert3rdOrderTensorToMatrix(H_inv_nF_ijk);
    Cel.block(9,9,3,3) = convert2ndOrderTensorToMatrix(2./h*H_inv_ij);
    
    Fastor::Tensor<double,21,21> Cel_fastor(Cel.data()); 
    // handle zero strain increment
    if ( Fastor::norm(dU_ftensor) <1e-14 && Fastor::norm(dSurface_strain_ftensor) <1e-14  ) {
      dStress_dStrain_ftensor = Cel_fastor;
      std::copy(dStress_dStrain_ftensor.data(), dStress_dStrain_ftensor.data() + 21*21, dStress_dStrain);
      return;
    }
    //visco elastic step
  
    Eigen::Ref< KelvinChainInterface::mapStateVarMatrix_Ju> creepStateVars_Ju( stateVarManager->kelvinStateVars_Ju );
    Eigen::Ref< KelvinChainInterface::mapStateVarMatric_Js> creepStateVars_Js( stateVarManager->kelvinStateVars_Js );
  
    const double dTimeDays = dT * timeToDays;
    
    auto [unitH_ij, unitZ_inv_ijkl, E_eff, nu_eff] = calculate_unitcompliance_interface(normal_ftensor, E_M, nu_M, E_I, nu_I, E_0, nu_0);
    
    // Use already available functionality convert Fastor tensors to Eigen matrices
    Eigen::Map<Eigen::Matrix<double,3,3, Eigen::RowMajor>>unitH_ij_eigen(unitH_ij.data());
    Eigen::Map<Eigen::Matrix<double,9,9, Eigen::RowMajor>>unitZ_inv_ab_eigen(unitZ_inv_ijkl.data());

    Vector3d creep_Ju_Increment = Vector3d::Zero();
    double creep_Ju_compliance = 0;

    KelvinChainInterface::evaluateKelvinChain(dTimeDays,
                                              elasticModuli_H_inv_ij,
                                              retardationTimes_Ju,
                                              creepStateVars_Ju,
                                              creep_Ju_compliance,
                                              creep_Ju_increment,
                                              1.0);
    
    Vector9d creep_Js_Increment = Vector9d::Zero();
    double creep_Js_compliance = 0;
 
    KelvinChainInterface::evaluateKelvinChain(dTimeDays,
                                              elasticModuli_Zijkl,
                                              retardationTimes_Js,
                                              creepStateVars_Js,
                                              creep_Js_compliance,
                                              creep_Js_increment,
                                              1.0);
    
    using namespace Marmot::ContinuumMechanics::Viscoelasticity;

    //Evaluate effective compliances due to the displacement jump and the surface stress

    double barC_Ju = 1./Ebar + zerothKelvinCompliance_Ju + creep_Ju_compliance;
    double barC_Js = 1./Ebar + zerothKelvinCompliance_Js + creep_Js_compliance;

    auto [H_inv_ij_effective, Z_ijkl_effective] = calculate_effective_properties(C_Ju, C_Js, normal_ftensor, E_eff, nu_eff);

    Vector3d deltaforce = H_inv_ij_effective*();

    
    enum {i,j,k,l};
    
    Tensor1D jumpU_ftensor = dU_ftensor(Fastor::seq(0,3),0)-dU_ftensor(Fastor::seq(3,Fastor::last),0);
    
    //std:: cout << "jumpU:" << jumpU_ftensor << std::endl;    
    
    Fastor::Tensor<double, 9,1> average_dSurface_strain_ftensor = 1./2.*(dSurface_strain_ftensor(Fastor::seq(0,9),0)+dSurface_strain_ftensor(Fastor::seq(9,Fastor::last),0));
    auto average_dSurface_strain_ftensor_reshape = Fastor::reshape<3,3>(average_dSurface_strain_ftensor);
    force_ftensor     += 2./h*Fastor::einsum<Fastor::Index<i,j>,Fastor::Index<j>,Fastor::OIndex<i>>(H_inv_ij, jumpU_ftensor);
    force_ftensor     += Fastor::einsum<Fastor::Index<i,j,k>,Fastor::Index<j,k>,Fastor::OIndex<i>>(H_inv_nF_ijk, average_dSurface_strain_ftensor_reshape);
    
    surface_stress_ftensor += h/2.*Fastor::einsum<Fastor::Index<i,j,k,l>,Fastor::Index<k,l>,Fastor::OIndex<i,j>>(Z_ijkl, average_dSurface_strain_ftensor_reshape);
    surface_stress_ftensor += h/2.*Fastor::einsum<Fastor::Index<i,j,k,l>,Fastor::Index<k,l>,Fastor::OIndex<i,j>>(Yn_H_inv_Fn_ijkl, average_dSurface_strain_ftensor_reshape);
    surface_stress_ftensor += Fastor::einsum<Fastor::Index<i,j,k>,Fastor::Index<k>,Fastor::OIndex<i,j>>(H_inv_nF_ijk, jumpU_ftensor); 
    
    dStress_dStrain_ftensor = Cel_fastor;
    
    std::copy(force_ftensor.data(), force_ftensor.data() + 3, force);
    std::copy(surface_stress_ftensor.data(), surface_stress_ftensor.data() + 3*3, surface_stress);
    std::copy(dStress_dStrain_ftensor.data(), dStress_dStrain_ftensor.data() + 21*21, dStress_dStrain);

      };

} // namespace Marmot::Materials

int main() {
    const double E_M= 200.0;  // Young's modulus (GPa)
    const double nu_M = 0.3;   // Poisson's ratio
    const double E_I = 200.0;
    const double nu_I = 0.3;
    const double E_0 = 200.0*0.1;
    const double nu_0 = 0.3;
    const double h = 1e0;
    
    const double material_props[] = {E_M, nu_M, E_I, nu_I, E_0, nu_0, h};
    
    Tensor1D force(0.f);
    Tensor2D surface_stress(0.f);
    Fastor::Tensor<double, 21,21> dStress_dStrain(0.f);
    const Fastor::Tensor<double, 6,1> dU(0.0);
    const Fastor::Tensor<double, 18,1> dSurface_strain(0.1);
    const Tensor1D normal = {0.0,1.0,0.0};
   
    double zero = 0.0;  
    const double* timeOld = &zero;
    double dT = 0.1;
    double pNewDT = 1.e36;
    
    //Marmot::Materials::LinearElasticInterface(material_props, 6, 0).computeStress( force, surface_stress, dStress_dStrain, dU, dSurface_strain, normal, timeOld, dT, pNewDT);

    double* forcePtr = force.data();
    double* surface_stressPtr = surface_stress.data();
    double* dStress_dStrainPtr = dStress_dStrain.data();
    const double* dUPtr = dU.data();
    const double* dSurface_strainPtr = dSurface_strain.data();
    const double* normalPtr = normal.data();
    
    Marmot::Materials::LinearElasticInterface(material_props, 6, 0).computeStress( forcePtr, surface_stressPtr, dStress_dStrainPtr, dUPtr, dSurface_strainPtr, normalPtr, timeOld, dT, pNewDT);


    return 0;
}
