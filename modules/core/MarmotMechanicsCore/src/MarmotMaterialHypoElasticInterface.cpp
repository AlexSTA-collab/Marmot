#include "Marmot/MarmotMaterialHypoElasticInterface.h"
#include "Marmot/HughesWinget.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotKinematics.h"
#include "Marmot/MarmotLowerDimensionalStress.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotVoigt.h"
#include <iostream>
#include "Fastor/Fastor.h"
#include "Marmot/interface_material_helper_functions.h"

using namespace Eigen;
using namespace Fastor;

using Tensor1D = Fastor::Tensor<double,3>;
using Tensor2D = Fastor::Tensor<double,3,3>;
using Tensor3D = Fastor::Tensor<double,3,3,3>;
using Tensor4D = Fastor::Tensor<double,3,3,3,3 >;

void MarmotMaterialHypoElasticInterface::setCharacteristicElementLength( double length )
{
  characteristicElementLength = length;
}

void MarmotMaterialHypoElasticInterface::computeStress( Tensor1D&  force,
                                                        Tensor2D&  surface_stress,
                                                        Fastor::Tensor<double, 21,21>& dStress_dStrain,
                                                        const Fastor::Tensor<double, 6,1>& dU,
                                                        const Fastor::Tensor<double, 18,1>& dSurface_strain,
                                                        const Tensor1D& normal,
                                                        const double* timeOld,
                                                        const double  dT,
                                                        double&       pNewDT  )
 {
    // elasticity parameters
    const double& E_M  = this->materialProperties[0];
    const double& nu_M = this->materialProperties[1];
    const double& E_I  = this->materialProperties[2];
    const double& nu_I = this->materialProperties[3];
    const double& E_0  = this->materialProperties[4];
    const double& nu_0 = this->materialProperties[5];
    
    // map to force, surface stress, displacement, surface strain, normal and tangent stiffness 
    // use Fastor because we really need to use the einsum  
    Tensor1D  F = force ;
    Tensor2D  S = surface_stress;
    auto  dS_dE = dStress_dStrain;
    
    auto [Z_ijkl, H_inv_ij, H_inv_nF_ijk, Yn_H_inv_Fn_ijkl] = calculate_interface_material_parameters(normal, E_M, nu_M, E_I, nu_I, E_0, nu_0); 
    std:: cout << "Z_ijkl:" << Z_ijkl << std::endl;
    std:: cout << "H_inv_ij:" << H_inv_ij << std::endl;
    std:: cout << "H_inv_nF_ijk:" << H_inv_nF_ijk << std::endl;
    std:: cout << "Yn_H_inv_nF_ijkl:" << Yn_H_inv_Fn_ijkl << std::endl;
    
    //Assign the material matrices to a larger structure. (Not necessary ...)
    Eigen::Matrix<double, 21, 21> Cel = Eigen::Matrix<double, 21,21>::Zero();

    Cel.block(0,0,9,9) = convert4thOrderTensorToMatrix(Z_ijkl);
    Cel.block(12,12,9,9) = convert4thOrderTensorToMatrix(Yn_H_inv_Fn_ijkl);
    Cel.block(0,9,9,3) = convert3rdOrderTensorToMatrix(H_inv_nF_ijk);
    Cel.block(9,9,3,3) = convert2ndOrderTensorToMatrix(H_inv_ij);
    
    Fastor::Tensor<double,21,21> Cel_fastor(Cel.data()); 
    // handle zero strain increment
    if ( Fastor::norm(dU) <1e-14 && Fastor::norm(dSurface_strain) <1e-14  ) {
      dS_dE = Cel_fastor;
      return;
    }
    //elastic step
    enum {i,j,k,l};
    
    Tensor1D jumpU = dU(Fastor::seq(0,3))-dU(Fastor::seq(3,Fastor::last));
    Fastor::Tensor<double, 9,1> average_dSurface_strain = 1./2.*(dSurface_strain(Fastor::seq(0,9))+dSurface_strain(Fastor::seq(9,Fastor::last)));
    auto average_dSurface_strain_reshape = Fastor::reshape<3,3>(average_dSurface_strain);
    force     += Fastor::einsum<Fastor::Index<i,j>,Fastor::Index<j>,Fastor::OIndex<i>>(H_inv_ij, jumpU);
    force     += Fastor::einsum<Fastor::Index<i,j,k>,Fastor::Index<j,k>,Fastor::OIndex<i>>(H_inv_nF_ijk, average_dSurface_strain_reshape);
    
    surface_stress += Fastor::einsum<Fastor::Index<i,j,k,l>,Fastor::Index<k,l>,Fastor::OIndex<i,j>>(Z_ijkl, average_dSurface_strain_reshape);
    surface_stress += Fastor::einsum<Fastor::Index<i,j,k,l>,Fastor::Index<k,l>,Fastor::OIndex<i,j>>(Yn_H_inv_Fn_ijkl, average_dSurface_strain_reshape);
    surface_stress += Fastor::einsum<Fastor::Index<i,j,k>,Fastor::Index<k>,Fastor::OIndex<i,j>>(H_inv_nF_ijk, jumpU);
    
    dS_dE = Cel_fastor;
    std:: cout << "dS_dE:" << dS_dE << std::endl;
    std:: cout << "force:" << force << std::endl;
    std:: cout << "surface_stress:" << surface_stress << std::endl;    
} // namespace Marmot::Materials



