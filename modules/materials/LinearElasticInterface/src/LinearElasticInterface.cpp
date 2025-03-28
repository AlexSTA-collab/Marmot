#include "../include/Marmot/LinearElasticInterface.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotUtility.h"
#include "Marmot/MarmotVoigt.h"

#include "Fastor/Fastor.h"
#include "../include/Marmot/interface_material_helper_functions.h"

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
    assert( nMaterialProperties == 6 || nMaterialProperties == 8 || nMaterialProperties == 12 );
  }

  void LinearElasticInterface::computeStress( Tensor1D&  force,
                                              Tensor2D&  surface_stress,
                                              Fastor::Tensor<double,21,21>& dStress_dStrain,
                                              const Fastor::Tensor<double,6,1>& dU,
                                              const Fastor::Tensor<double,18,1>& dSurface_strain,
                                              const Tensor1D& normal,
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
   
      };

} // namespace Marmot::Materials

int main() {
    const double E_M= 200.0;  // Young's modulus (GPa)
    const double nu_M = 0.3;   // Poisson's ratio
    const double E_I = 200.0;
    const double nu_I = 0.3;
    const double E_0 = 200.0*0.1;
    const double nu_0 = 0.3;
    
    const double material_props[] = {E_M, nu_M, E_I, nu_I, E_0, nu_0};
    
    Tensor1D force(0.f);
    Tensor2D surface_stress(0.f);
    Fastor::Tensor<double, 21,21> dStress_dStrain(0.f);
    const Fastor::Tensor<double, 6,1> dU(0.0);
    const Fastor::Tensor<double, 18,1> dSurface_strain(0.1);
    Tensor1D normal = {0.0,1.0,0.0};
   
    double zero = 0.0;  
    const double* timeOld = &zero;
    double dT = 0.1;
    double pNewDT = 1.e36;
    
  Marmot::Materials::LinearElasticInterface(material_props, 6, 0).computeStress( force, surface_stress, dStress_dStrain, dU, dSurface_strain, normal, timeOld, dT, pNewDT);  
    return 0;
}
