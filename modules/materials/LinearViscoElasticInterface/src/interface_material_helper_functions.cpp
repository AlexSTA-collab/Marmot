/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck,
 * 2020 - today
 *
 * festigkeitslehre@uibk.ac.at
 *
 * Alexandros Stathas alexandros.stathas@boku.ac.at
 *
 * This file is part of the MAteRialMOdellingToolbox (marmot).
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of marmot.
 * ---------------------------------------------------------------------
 */
//#`:pragma once
//#include "Marmot/MarmotTypedefs.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/src/Core/Matrix.h>
#include <Fastor/Fastor.h>
#include <Fastor/tensor/AbstractTensorFunctions.h>
#include <Fastor/tensor_algebra/einsum.h>
#include <Fastor/tensor_algebra/indicial.h>
#include <ostream>
#include <tuple>
#include <unsupported/Eigen/CXX11/Tensor>
#include <cassert>
#include <cmath>
#include <iostream>
#include "../include/Marmot/interface_material_helper_functions.h"

using namespace Eigen;
using namespace Fastor;

using Tensor1D = Fastor::Tensor<double,3>;
using Tensor2D = Fastor::Tensor<double,3,3>;
using Tensor3D = Fastor::Tensor<double,3,3,3>;
using Tensor4D = Fastor::Tensor<double,3,3,3,3 >;

void assert_consistent_arrays_indices(
   const Tensor4D& L , 
   Tensor4D& A , 
   Tensor4D& B , 
   const Tensor4D& M , 
   const Tensor2D& J , 
   const Tensor2D& I
    );


void assert_equivalent_F_Falt_Y(
    Tensor4D F ,
    Tensor4D Y ,
    Tensor4D A_0 ,
    Tensor4D L_0 ,
    Tensor4D A_M ,
    Tensor4D L_M ,
    Tensor4D A_I ,
    Tensor4D L_I 
    );

enum {a,i,b,j,k,l,m,n};

Tensor2D compute_inv( const Tensor2D& I, Tensor2D& Q )
{
  // Solve G in batch mode using LU decomposition
  Eigen::Matrix3d Im (I.data());
  Eigen::Matrix3d Qm (Q.data());
  Eigen::Matrix3d G_mat = Qm.fullPivLu().solve( Im );

  // Store result back into G (reinterpret as a 3D Fastor tensor)
  Tensor2D G(0);
  Eigen::Map< Eigen::Matrix< double, 3, 3 , Eigen::RowMajor > >( G.data() ) = G_mat;

  return G;
}
std::tuple<const Tensor4D, Tensor4D, const Tensor4D, Tensor4D, Tensor2D>interface_geometry_system_couplings( const Tensor2D& I,
                                          const Tensor2D& J,
                                          const Tensor2D& N,
                                          const Tensor2D& T,
                                          const Tensor4D& L,
                                          const Tensor4D& M
                                         )
{

  // **Compute Q = einsum("aibjq,ijq->abq", L_expanded, N)**
  Tensor2D Q = Fastor::einsum<Fastor::Index<a,i,b,j>,Fastor::Index<i,j>,Fastor::OIndex<a,b>>(L, N);
  Tensor2D G = compute_inv(I, Q );

  Tensor4D A = Fastor::einsum<Fastor::Index<a,b>,Fastor::Index<i,j>,Fastor::OIndex<a,i,b,j>>(G, N );
  Tensor4D LA = Fastor::einsum<Fastor::Index<a,i,m,n>,Fastor::Index<m,n,b,j>,Fastor::OIndex<a,i,b,j>>(L,A);
  Tensor4D LAL = Fastor::einsum<Fastor::Index<a,i,m,n>,Fastor::Index<m,n,b,j>,Fastor::OIndex<a,i,b,j>>(LA, L );
  Tensor4D B = L - LAL;

  // check that the expressions are consistent
  assert_consistent_arrays_indices(L, A, B, M, J, I);
  
  return std::make_tuple(M, B, L, A, G);
}

std::tuple<Tensor4D, Tensor2D, Tensor3D, Tensor4D> calculate_material_matrices(
                                 const Tensor1D& normal,
                                 const Tensor2D& I, 
                                 const Tensor2D& J, 
                                 const Tensor2D& N, 
                                 const Tensor2D& T, 
                                 const Tensor4D& C_0_aibj,
                                 const Tensor4D& C_M_aibj,
                                 const Tensor4D& C_I_aibj,
                                 const Tensor4D& S_0_aibj, 
                                 const Tensor4D& S_M_aibj, 
                                 const Tensor4D& S_I_aibj
                                )
  {
    Tensor2D G_0 ;
    Tensor4D A_0 ;
    Tensor4D B_0 ;
    Tensor4D M_0 ;
    Tensor4D L_0 ;

    Tensor2D G_M ;
    Tensor4D A_M ;
    Tensor4D B_M ;
    Tensor4D M_M ;
    Tensor4D L_M ;

    Tensor2D G_I ;
    Tensor4D A_I ;
    Tensor4D B_I ;
    Tensor4D M_I ;
    Tensor4D L_I ;
    
    //std:: cout << "C_0_aibj:" << C_0_aibj << std::endl;
    //std:: cout << "C_M_aibj:" << C_M_aibj << std::endl;
    //std:: cout << "C_I_aibj:" << C_I_aibj << std::endl;
    
    std::tie(M_0, B_0, L_0, A_0, G_0) = interface_geometry_system_couplings(
        I, J, N, T, C_0_aibj, S_0_aibj
    );

    std::tie(M_M, B_M, L_M, A_M, G_M) = interface_geometry_system_couplings(
        I, J, N, T, C_M_aibj, S_M_aibj
    );

    std::tie(M_I, B_I, L_I, A_I, G_I) = interface_geometry_system_couplings(
        I, J, N, T, C_I_aibj, S_I_aibj
    );

    //std:: cout << "M_M:" << M_M << std::endl;
    //std:: cout << "B_M:" << B_M << std::endl;
    
    //std:: cout << "M_I:" << M_I << std::endl;
    //std:: cout << "B_I:" << B_I << std::endl;

    //std:: cout << "M_0:" << M_0 << std::endl;
    //std:: cout << "B_0:" << B_0 << std::endl;
    
    //std:: cout<<"F_1st:" << Fastor::einsum<Fastor::Index<a,i,m,n>,Fastor::Index<m,n,b,j>,Fastor::OIndex<a,i,b,j>>(M_0, B_0) << std::endl;
    
    Tensor4D F = 2.0 * Fastor::einsum<Fastor::Index<a,i,m,n>,Fastor::Index<m,n,b,j>,Fastor::OIndex<a,i,b,j>>(M_0, B_0);
    F-= Fastor::einsum<Fastor::Index<a,i,m,n>,Fastor::Index<m,n,b,j>,Fastor::OIndex<a,i,b,j>>(M_M, B_M);
    F-= Fastor::einsum<Fastor::Index<a,i,m,n>,Fastor::Index<m,n,b,j>,Fastor::OIndex<a,i,b,j>>(M_I, B_I);

    Tensor4D Y = Fastor::einsum<Fastor::Index<a,i,m,n>,Fastor::Index<m,n,b,j>,Fastor::OIndex<a,i,b,j>>(L_M, A_M);
    Y+= Fastor::einsum<Fastor::Index<a,i,m,n>,Fastor::Index<m,n,b,j>,Fastor::OIndex<a,i,b,j>>(L_I, A_I);
    Y-= 2.0 * Fastor::einsum<Fastor::Index<a,i,m,n>,Fastor::Index<m,n,b,j>,Fastor::OIndex<a,i,b,j>>(L_0, A_0);  

    assert_equivalent_F_Falt_Y(F , Y , A_0 , L_0 , A_M , L_M , A_I , L_I );

    Tensor2D H = 2.0 * G_0 - G_M - G_I;
    
    //std:: cout << "from material B_M:" << B_M << std::endl;
    //std:: cout << "from material B_I:" << B_I << std::endl;
    //std:: cout << "from material B_0:" << B_0 << std::endl;
    
    Tensor4D Z = B_M + B_I - 2.0 * B_0;
    Tensor2D H_inv = compute_inv(I, H );
    //std:: cout << "from material H_inv:" << H_inv << std::endl;    
    //Tensor2D I_approx = H * H_inv;
    //std:: cout << "I_approx H H_inv:" << I_approx << std::endl;
    
    Tensor3D nF = Fastor::einsum<Fastor::Index<a>,Fastor::Index<a,i,b,j>,Fastor::OIndex<i,b,j>>(normal, F);
    Tensor3D Fn = Fastor::einsum<Fastor::Index<a,i,b,j>,Fastor::Index<j>,Fastor::OIndex<a,i,b>>(F, normal);
    Tensor3D Yn = Fastor::einsum<Fastor::Index<a,i,b,j>,Fastor::Index<j>,Fastor::OIndex<a,i,b>>(Y, normal);
    Tensor3D H_inv_nF = Fastor::einsum<Fastor::Index<a,b>,Fastor::Index<b,i,j>,Fastor::OIndex<a,i,j>>(H_inv, nF);
    Tensor4D Yn_H_inv_Fn = Fastor::einsum<Fastor::Index<a,i,m>,Fastor::Index<m,n>,Fastor::Index<n,b,j>,Fastor::OIndex<a,i,b,j>>(Yn, H_inv, Fn);
    //std:: cout << "F:" << F << std::endl;
    return std::make_tuple(Z, H_inv, H_inv_nF, Yn_H_inv_Fn);
  }

void assert_consistent_arrays_indices(
   const Tensor4D& L , 
   Tensor4D& A , 
   Tensor4D& B , 
   const Tensor4D& M , 
   const Tensor2D& J , 
   const Tensor2D& I
    )
{
    //P_norm = np.einsum( "abq,ijq->aibjq", J, N ) #aibj
    //P_tan  = np.einsum( "abq,ijq->aibjq", J, T ) #aibj

    Tensor4D JI = Fastor::einsum<Fastor::Index<a,b>,Fastor::Index<i,j>,Fastor::OIndex<a,i,b,j>>(J, I);
    Tensor4D LM = Fastor::einsum<Fastor::Index<a,i,m,n>,Fastor::Index<m,n,b,j>,Fastor::OIndex<a,i,b,j>>(L, M);

    //Fastor::print(L);
    //Fastor::print(M);
    //Fastor::print(LM);
    //Fastor::print(JI);
  
  //Tensor4D ML = Fastor::einsum<Fastor::Index<a,i,m,n>,Fastor::Index<m,n,b,j>,Fastor::OIndex<a,i,b,j>>(M, L);
   
    double atol = 1e-8;
    Eigen::Map<const Eigen::Matrix<double,Eigen::Dynamic, 1>> LM_flat(LM.data(),LM.size());
    Eigen::Map<const Eigen::Matrix<double,Eigen::Dynamic, 1>> JI_flat(JI.data(), JI.size());
    double max_diff1 = (LM_flat-JI_flat).cwiseAbs().maxCoeff();
    if (max_diff1 > atol){
      std::cerr << "Assertion failed: Arrays not  LM, JI, within the tolerance."<< std::endl;
      std::abort();
    }

    Tensor4D LA = Fastor::einsum<Fastor::Index<a,i,m,n>,Fastor::Index<m,n,b,j>,Fastor::OIndex<a,i,b,j>>(L, A);
    Tensor4D BM = Fastor::einsum<Fastor::Index<a,i,m,n>,Fastor::Index<m,n,b,j>,Fastor::OIndex<a,i,b,j>>(B, M);
    Tensor4D LA_plus_BM = LA + BM;
    Eigen::Map<const Eigen::Matrix<double,Eigen::Dynamic, 1>> LA_plus_BM_flat(LA_plus_BM.data(),LA_plus_BM.size());
  
    double max_diff2 = (LA_plus_BM_flat-JI_flat).cwiseAbs().maxCoeff();
    if (max_diff2 > atol){
      std::cerr << "Arrays LA+BM, JI are not equal within the tolerance."<< std::endl;
      std::abort(); 
  }

    //std::cout<<"Passed the assertion about proper definition of matrices A, B";
}

void assert_equivalent_F_Falt_Y(
    Tensor4D F ,
    Tensor4D Y ,
    Tensor4D A_0 ,
    Tensor4D L_0 ,
    Tensor4D A_M ,
    Tensor4D L_M ,
    Tensor4D A_I ,
    Tensor4D L_I 
    )
{
    Tensor4D F_alt = -2.0 * Fastor::einsum<Fastor::Index<a,i,m,n>,Fastor::Index<m,n,b,j>,Fastor::OIndex<a,i,b,j>>(A_0, L_0);
    F_alt+= Fastor::einsum<Fastor::Index<a,i,m,n>,Fastor::Index<m,n,b,j>,Fastor::OIndex<a,i,b,j>>(A_M, L_M);
    F_alt+= Fastor::einsum<Fastor::Index<a,i,m,n>,Fastor::Index<m,n,b,j>,Fastor::OIndex<a,i,b,j>>(A_I, L_I);
    
    double atol = 1e-8;
    Eigen::Map<const Eigen::Matrix<double,Eigen::Dynamic,1>> F_alt_flat(F_alt.data(), F_alt.size());
    Eigen::Map<const Eigen::Matrix<double,Eigen::Dynamic,1>> F_flat(F.data(), F.size());
   
    double max_diff = (F_alt_flat-F_flat).cwiseAbs().maxCoeff();
    if (max_diff > atol){
      std::cerr << "Arrays are not equal within the tolerance."<< std::endl;
      std::abort();
    //assert np.allclose(
    //    Y.T, F, atol=10.0 ** (-8)
    //), "Arrays are not equal within the tolerance."
    }
}

template <typename T>
Tensor4D voigt_full_to_tensor(const Matrix<T, 9, 9>& M_matrix) 
{
    std::vector<std::pair<int, int>> matrix_map = {
        {0, 0}, {0, 1}, {0, 2},
        {1, 0}, {1, 1}, {1, 2},
        {2, 0}, {2, 1}, {2, 2}   
    };

    // Initialize a 4th-order tensor with zeros (Fastor tensor)
    Tensor4D M_tensor(0.0);

    // Map Voigt indices to tensor indices and apply proper scaling
    for (int I = 0; I < 9; ++I) {
        for (int J = 0; J < 9; ++J) {
            int i = matrix_map[I].first;
            int j = matrix_map[I].second;
            int m = matrix_map[J].first;
            int n = matrix_map[J].second;

            T value = M_matrix(I, J);
      
            // Assign only necessary components (avoid redundant copies)
            M_tensor(i, j, m, n) += value;
        }
    }

    return M_tensor;
}

// Convert 4th-order Fastor tensor (3x3x3x3) to Eigen 9x9 matrix
Eigen::Matrix<double,9,9> convert4thOrderTensorToMatrix(const Tensor4D& tensor) {
    Eigen::Matrix<double,9,9> matrix(9, 9);
    
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                for (int l = 0; l < 3; ++l) {
                    int row = 3 * i + j;  // Convert (i, j) to single index
                    int col = 3 * k + l;  // Convert (k, l) to single index
                    matrix(row, col) = tensor(i, j, k, l);
                }
            }
        }
    }
    
    return matrix;
}

// Convert 3rd-order Fastor tensor (3x3x3) to Eigen 9x3 matrix
Eigen::Matrix<double,9,3> convert3rdOrderTensorToMatrix(const Tensor3D& tensor) {
    Eigen::Matrix<double,9,3> matrix(9, 3);
    
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                int row = 3 * i + j;  // Convert (i, j) to single index
                matrix(row, k) = tensor(i, j, k);
            }
        }
    }
    
    return matrix;
}

// Convert 2nd-order Fastor tensor (3x3) to Eigen 3x3 matrix
Eigen::Matrix<double,3,3>convert2ndOrderTensorToMatrix(const Tensor2D& tensor){
    Eigen::Matrix<double,3,3> matrix(3, 3);
    
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
                matrix(i, j) = tensor(i, j);
            }
        }   
    return matrix;
}

// Function to create the isotropic elasticity tensor in Voigt notation (6x6)
Eigen::Matrix<double, 9, 9> create_isotropic_elasticity_tensor(const double E, const double nu) {
    Eigen::Matrix<double, 9, 9> C_voigt_full;
    C_voigt_full.setZero();  // Initialize

    double lambda = (E * nu) / ((1 + nu) * (1 - 2 * nu));
    double mu = E / (2 * (1 + nu));

    C_voigt_full(0, 0) = C_voigt_full(1, 1) = C_voigt_full(2, 2) = lambda + 2 * mu;
    C_voigt_full(0, 1) = C_voigt_full(1, 0) = lambda;
    C_voigt_full(0, 2) = C_voigt_full(2, 0) = lambda;
    C_voigt_full(1, 2) = C_voigt_full(2, 1) = lambda;

    C_voigt_full(3, 3) = C_voigt_full(4, 4) = C_voigt_full(5, 5) = C_voigt_full(6, 6) = C_voigt_full(7, 7) = C_voigt_full(8,8) = mu;

    return C_voigt_full;
}

std::tuple<Tensor4D, Tensor2D, Tensor3D, Tensor4D> calculate_interface_material_parameters(
                                                                                           const Tensor1D& normal,
                                                                                           const double& E_M,
                                                                                           const double& nu_M,
                                                                                           const double& E_I,
                                                                                           const double& nu_I,
                                                                                           const double& E_0,
                                                                                           const double& nu_0)
{
    Eigen::Matrix<double, 9, 9> II = Eigen::Matrix<double, 9, 9>::Identity();
    Eigen::Matrix<double, 9, 9> C_M_voigt_full = create_isotropic_elasticity_tensor(E_M, nu_M);
    Eigen::Matrix<double, 9, 9> C_I_voigt_full = create_isotropic_elasticity_tensor(E_I, nu_I);
    Eigen::Matrix<double, 9, 9> C_0_voigt_full = create_isotropic_elasticity_tensor(E_0, nu_0);
    
    Eigen::Matrix<double, 9, 9> S_M_voigt_full = C_M_voigt_full.fullPivLu().solve( II );
    Eigen::Matrix<double, 9, 9> S_I_voigt_full = C_I_voigt_full.fullPivLu().solve( II );
    Eigen::Matrix<double, 9, 9> S_0_voigt_full = C_0_voigt_full.fullPivLu().solve( II );

    
    //Tensor1D normal = {0.0,1.0,0.0};
    
    Tensor2D I = {{1.0,0.0,0.0},
                  {0.0,1.0,0.0},
                  {0.0,0.0,1.0}
                 };

    Tensor2D J = {{1.0,0.0,0.0},
                  {0.0,1.0,0.0},
                  {0.0,0.0,1.0}
                 };
    Tensor2D N = Fastor::einsum<Fastor::Index<i>, Fastor::Index<j>, Fastor::OIndex<i,j>>(normal, normal);
    
    Tensor2D T = I - N;

    Tensor4D C_M_aibj = voigt_full_to_tensor(C_M_voigt_full);
    Tensor4D C_I_aibj = voigt_full_to_tensor(C_I_voigt_full);
    Tensor4D C_0_aibj = voigt_full_to_tensor(C_0_voigt_full);

    Tensor4D S_M_aibj = voigt_full_to_tensor(S_M_voigt_full);
    Tensor4D S_I_aibj = voigt_full_to_tensor(S_I_voigt_full);
    Tensor4D S_0_aibj = voigt_full_to_tensor(S_0_voigt_full);
    // Step 3: Print one sample value
      
    auto [Z, H_inv, H_inv_nF, Yn_H_inv_Fn] = calculate_material_matrices(normal, I, J, N, T, C_0_aibj, C_M_aibj, C_I_aibj, S_0_aibj, S_M_aibj, S_I_aibj);
    
    return {Z, H_inv, H_inv_nF, Yn_H_inv_Fn}; 
}


std::tuple<Tensor2D, Tensor4D, double, double>calculate_unitcompliance_interface(
                                                                                           const Tensor1D& normal,
                                                                                           const double& E_M,
                                                                                           const double& nu_M,
                                                                                           const double& E_I,
                                                                                           const double& nu_I,
                                                                                           const double& E_0,
                                                                                           const double& nu_0)
{
    Tensor2D I = {{1.0,0.0,0.0},
                  {0.0,1.0,0.0},
                  {0.0,0.0,1.0}
                 };

    Tensor2D N = Fastor::einsum<Fastor::Index<i>, Fastor::Index<j>, Fastor::OIndex<i,j>>(normal, normal);
    Tensor2D T = I-N;

    double mu_bar = (E_M)/(2.*(1.+nu_M))+ (E_I)/(2.*(1.+nu_I)) -  2*(E_0)/(2.*(1.+nu_0));
    double lambda_bar = (E_M*nu_M)/((1+nu_M)*(1-2*nu_M))+ (E_I*nu_I)/((1+nu_I)*(1-2*nu_I))- 2*(E_0*nu_0)/((1+nu_0)*(1-2*nu_0));

    double E_bar = mu_bar*(3.*lambda_bar+2.*mu_bar)/(lambda_bar+mu_bar);
    double nu_bar = lambda_bar/(2.*(lambda_bar+mu_bar));
  
    Eigen::Matrix<double,9,9> unitC_bar_voigt_full = create_isotropic_elasticity_tensor(1.0 , nu_bar);
    Tensor4D unitC_bar_aibj = voigt_full_to_tensor(unitC_bar_voigt_full);
    Tensor2D Q_ab = Fastor::einsum<Fastor::Index<a,i,b,j>, Fastor::Index<i>, Fastor::Index<j>, Fastor::OIndex<a,b>>(unitC_bar_aibj, normal, normal);
    Tensor2D unitH_bar_ab = compute_inv(I, Q_ab );

    Tensor4D unitCs_bar_klmn = Fastor::einsum<Fastor::Index<k,a>, Fastor::Index<l,i>, Fastor::Index<m,b>, Fastor::Index<n,j>, Fastor::Index<a,i,b,j>, Fastor::OIndex<k,l,m,n>>(T,T,T,T,unitC_bar_aibj); 

    return  {unitH_bar_ab, unitCs_bar_klmn, E_bar, nu_bar}; 
}

std::tuple<Eigen::Matrix<double,3,3>, Eigen::Matrix<double,9,9>>calculate_effective_properties(const double& barC_Ju, 
                                                             const double& barC_Js, 
                                                             const Tensor1D& normal,
                                                             const double& E_bar,
                                                             const double& nu_bar
                                                             )
{
    Tensor2D I = {{1.0,0.0,0.0},
                  {0.0,1.0,0.0},
                  {0.0,0.0,1.0}
                 };

    Tensor2D N = Fastor::einsum<Fastor::Index<i>, Fastor::Index<j>, Fastor::OIndex<i,j>>(normal, normal);
    Tensor2D T = I-N;

        // create  the 9x9 matrix representation of the constitutive modulus (can be done with Voigt 6x6 but needs extra space for 9x9) 
    Eigen::Matrix<double, 9, 9> C_Ju_full = create_isotropic_elasticity_tensor(1./barC_Ju, nu_bar);
    Eigen::Matrix<double, 9, 9> C_Js_full = create_isotropic_elasticity_tensor(1./barC_Js, nu_bar);
    
    Tensor4D C_Ju_aibj = voigt_full_to_tensor(C_Ju_full);
    Tensor4D C_Js_aibj = voigt_full_to_tensor(C_Js_full);

    Tensor2D H_inv_ab = Fastor::einsum<Fastor::Index<a,i,b,j>, Fastor::Index<i>, Fastor::Index<j>, Fastor::OIndex<a,b>>(C_Ju_aibj, normal, normal);

    Tensor4D Cs_klmn = Fastor::einsum<Fastor::Index<k,a>, Fastor::Index<l,i>, Fastor::Index<m,b>, Fastor::Index<n,j>, Fastor::Index<a,i,b,j>, Fastor::OIndex<k,l,m,n>>(T,T,T,T,C_Js_aibj);

    Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> H_inv_ij_eigen(H_inv_ab.data());
    Eigen::Map<Eigen::Matrix<double, 9, 9, Eigen::RowMajor>> Z_ijkl_eigen(Cs_klmn.data());
    return  {H_inv_ij_eigen, Z_ijkl_eigen}; 
}
//int main() {
//    double E_M= 200.0;  // Young's modulus (GPa)
//    double nu_M = 0.3;   // Poisson's ratio
//    double E_I = 200.0;
//    double nu_I = 0.3;
//    double E_0 = 200.0*0.1;
//    double nu_0 = 0.3;
      
//    auto Ce = calculate_interface_material_parameters(E_M, nu_M, E_I, nu_I, E_0, nu_0); 
//    return 0;
//}
