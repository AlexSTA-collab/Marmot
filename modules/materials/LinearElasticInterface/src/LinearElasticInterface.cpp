#include "Marmot/LinearElasticInterface.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotUtility.h"
#include "Marmot/MarmotVoigt.h"

#include "Fastor/Fastor.h"
#include "Marmot/interface_material_helper_functions.h"
#include <Fastor/tensor/TensorMap.h>

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Eigen;
  using namespace Fastor;

  using Tensor1D = Fastor::Tensor<double,3>;
  using Tensor2D = Fastor::Tensor<double,3,3>;
  using Tensor3D = Fastor::Tensor<double,3,3,3>;
  using Tensor4D = Fastor::Tensor<double,3,3,3,3 >;
  
  LinearElasticInterface::LinearElasticInterface( const double* materialProperties, int nMaterialProperties, int materialNumber )
    : MarmotMaterialHypoElasticInterface::MarmotMaterialHypoElasticInterface( materialProperties, nMaterialProperties, materialNumber )
  {
    assert( nMaterialProperties == 7 || nMaterialProperties == 9 || nMaterialProperties == 13 );
  }

  void LinearElasticInterface::computeStress( double*  force,
                                              double*  surface_stress,
                                              double* dStress_dStrain,
                                              const double* dU,
                                              const double* dSurface_strain,
                                              const double* normal,
                                              const double* timeOld,
                                              const double  dT,
                                              double&       pNewDT)
  {
    // elasticity parameters
    const double& E_M  = this->materialProperties[0];
    const double& nu_M = this->materialProperties[1];
    const double& E_I  = this->materialProperties[2];
    const double& nu_I = this->materialProperties[3];
    const double& E_0  = this->materialProperties[4];
    const double& nu_0 = this->materialProperties[5];
    const double& h = this->materialProperties[6];

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
 

    //std:: cout<<"Before Calculation inside material"<<std::endl;
    //std:: cout << "force:" << force[0]<< "," << force[3-1] << std::endl;
    //std:: cout << "surface_stress:" << surface_stress[0] << "," << surface_stress[9-1]<< std::endl;
    //std:: cout << "dS_dE:" << dStress_dStrain[0] << "," <<  dStress_dStrain[21*21-1]<< std::endl;
    //std:: cout << "dU:" << dU[0] << "," << dU[6-1] << std::endl;
    //std:: cout << "dSurface_strain:" << dSurface_strain[0] << "," << dSurface_strain[18-1]<< std::endl;
    //std:: cout << "normal:" << normal[0] << "," << normal[3-1] << std::endl;

    auto [Z_ijkl, H_inv_ij, H_inv_nF_ijk, Yn_H_inv_Fn_ijkl] = calculate_interface_material_parameters(normal_ftensor, E_M, nu_M, E_I, nu_I, E_0, nu_0); 
    
    //std:: cout << "Z_ijkl:" << Z_ijkl << std::endl;
    //std:: cout << "H_inv_ij:" << H_inv_ij << std::endl;
    //std:: cout << "H_inv_nF_ijk:" << H_inv_nF_ijk << std::endl;
    //std:: cout << "Yn_H_inv_nF_ijkl:" << Yn_H_inv_Fn_ijkl << std::endl;
    
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
    //elastic step
    enum {i,j,k,l};
    
    Tensor1D jumpU_ftensor = dU_ftensor(Fastor::seq(0,3),0)-dU_ftensor(Fastor::seq(3,Fastor::last),0);    
    
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
