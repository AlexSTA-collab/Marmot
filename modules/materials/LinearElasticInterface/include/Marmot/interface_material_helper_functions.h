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
#pragma once
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/src/Core/Matrix.h>
#include <Fastor/Fastor.h>
#include <Fastor/tensor_algebra/indicial.h>
#include <unsupported/Eigen/CXX11/Tensor>
#include <tuple>
#include <cmath>

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

Tensor2D compute_inv( const Tensor2D& I, Tensor2D& Q );

std::tuple<const Tensor4D, Tensor4D, const Tensor4D, Tensor4D, Tensor2D>interface_geometry_system_couplings( const Tensor2D& I,
                                          const Tensor2D& J,
                                          const Tensor2D& N,
                                          const Tensor2D& T,
                                          const Tensor4D& L,
                                          const Tensor4D& M
                                         );

std::tuple<Tensor4D, Tensor2D, Tensor3D, Tensor4D> calculate_material_matrices(const Tensor1D& normal,
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
                                );

template <typename T>
Tensor4D voigt_full_to_tensor(const Matrix<T, 9, 9>& M_matrix);

Eigen::Matrix<double, 9, 9> create_isotropic_elasticity_tensor(double E, double nu);

Eigen::Matrix<double,9,9> convert4thOrderTensorToMatrix(const Tensor4D& tensor);
Eigen::Matrix<double,9,3> convert3rdOrderTensorToMatrix(const Tensor3D& tensor);
Eigen::Matrix<double,3,3>convert2ndOrderTensorToMatrix(const Tensor2D& tensor);

std::tuple<Tensor4D, Tensor2D, Tensor3D, Tensor4D>  calculate_interface_material_parameters(
                                                                    const Tensor1D& normal,
                                                                    const double& E_M,
                                                                    const double& nu_M,
                                                                    const double& E_I,
                                                                    const double& nu_I,
                                                                    const double& E_0,
                                                                    const double& nu_0);

std::tuple<Tensor2D, Tensor4D, double, double>calculate_unitcompliance_interface(
                                                                                 const Tensor1D& normal,
                                                                                 const double& E_M,
                                                                                 const double& nu_M,
                                                                                 const double& E_I,
                                                                                 const double& nu_I,
                                                                                 const double& E_0,
                                                                                 const double& nu_0);


std::tuple<Tensor2D, Tensor4D>calculate_effective_properties(const double& barC_Ju, 
                                                             const double& barC_Js, 
                                                             const Tensor1D& normal,
                                                             const double& E_bar,
                                                             const double& nu_bar
                                                             );

