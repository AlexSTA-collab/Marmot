#include "Marmot/MarmotMaterialInterface.h"
#include "Marmot/interface_material_helper_functions.h"
#include "Marmot/HughesWinget.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotKinematics.h"
#include "Marmot/MarmotLowerDimensionalStress.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotVoigt.h"
#include <iostream>

using namespace Eigen;
using namespace Fastor;

using Tensor1D = Fastor::Tensor<double,3>;
using Tensor2D = Fastor::Tensor<double,3,3>;
using Tensor3D = Fastor::Tensor<double,3,3,3>;
using Tensor4D = Fastor::Tensor<double,3,3,3,3 >;

using namespace Marmot;
using namespace Marmot::ContinuumMechanics::TensorUtility;

void MarmotMaterialInterface::assignStateVars( double* stateVars, int nStateVars )
  {
    if ( nStateVars < getNumberOfRequiredStateVars() )
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": Not sufficient stateVars!" );

    managedStateVars = std::make_unique< InterfaceVonMisesModelStateVarManager >( stateVars );
  }

  StateView MarmotMaterialInterface::getStateView( const std::string& stateName )
  {
    return managedStateVars->getStateView( stateName );
  }

void MarmotMaterialInterface::setCharacteristicElementLength( double length )
{
  characteristicElementLength = length;
}

void MarmotMaterialInterface::computeStress(double Tensor1D&  force,
                                            double Tensor2D&  surface_stress,
                                            Fastor<double, 21,21>& dStress_dStrain,
                                            const Fastor<double, 6,1>& dU,
                                            const Fastor<double, 18,1>& dSurface_strain,
                                            const Tensor1D& normal,
                                            const double* timeOld,
                                            const double  dT,
                                            double&       pNewDT )

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
    
    auto [Z_ijkl, H_inv_ij, H_inv_nF_ijk, Yn_H_inv_Fn_ijkl] = calculate_interface_material_parameters(E_M, nu_M, E_I, nu_I, E_0, nu_0, normal); 
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
    
    // handle zero strain increment
    if ( dE.isZero( 1e-14 ) ) {
      dS_dE = Cel;
      return;
    }
    //elastic step
    enum {i,j,k,l};
    
    Tensor1D jumpU = dU(seq(0,3))-dU(seq(3,last));
    Fastor<double, 9,1> average_dSurface_strain = 1./2.*(dSurface_strain(seq(0,9))+dSurface_strain(seq(9,last)));
    auto average_dSurface_strain = reshape<3,3>(average_dSurface_strain_j);
    force     += Fastor::einsum<Fastor::Index<i,j>,Fastor::Index<j>,Fastor::OIndex<i>>(H_inv_ij, jumpU);
    force     += Fastor::einsum<Fastor::Index<i,j,k>,Fastor::Index<j,k>,Fastor::OIndex<i>>(H_inv_ijk, average_dSurface_strain);
    
    surface_stress += Fastor::einsum<Fastor::Index<i,j,k,l>,Fastor::Index<k,l>,Fastor::OIndex<i,j>>(Z_ijkl, average_dSurface_strain);
    surface_stress +=Fastor::einsum<Fastor::Index<i,j,k,l>,Fastor::Index<k,l>,Fastor::OIndex<i,j>>(Yn_H_inv_Fn_ijkl, average_dSurface_strain);
    surface_stress += Fastor::einsum<Fastor::Index<i,j,k>,Fastor::Index<k>,Fastor::OIndex<i,j>>(H_inv_nF_ijk, jumpU);
    
    dS_dE = Cel;
    std:: cout << "dS_dE:" << dS_dE << std::endl;
    std:: cout << "force:" << force << std::endl;
    std:: cout << "surface_stress:" << surface_stress << std::endl;    
} // namespace Marmot::Materials


int main() {
    double E_M= 200.0;  // Young's modulus (GPa)
    double nu_M = 0.3;   // Poisson's ratio
    double E_I = 200.0;
    double nu_I = 0.3;
    double E_0 = 200.0*0.1;
    double nu_0 = 0.3;
    
    double material_props[] = {E_M, nu_M, E_I, nu_I, E_0, nu_0};
    
    Tensor1D force(0.f);
    Tensor2D surface_stress(0.f);
    Fastor<double, 21,21> dStress_dStrain(0.f);
    Fastor<double, 6,1> dU(0.f);
    Fastor<double, 18,1> dSurface_strain(0.1)
    Tensor1D normal = {0.0,1.0,0.0};
    
    double timeOld = 0.0;
    double dT = 0.1;
    double pNewDT = 1.e36;
    
    MarmotMaterialInterface(material_props, 6, 0).computeStress( force, surface_stress, dStress_dStrain, dU, dSurface_strain, normal, timeOld, dT, pNewDT);  
    return 0;
}
