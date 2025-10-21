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
// #`:pragma once
// #include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotInterfaceMaterialHelperFunctions.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotTypedefs.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/src/Core/Matrix.h>
#include <Fastor/Fastor.h>
#include <Fastor/expressions/linalg_ops/unary_norm_op.h>
#include <Fastor/expressions/linalg_ops/unary_trans_op.h>
#include <Fastor/tensor/AbstractTensorFunctions.h>
#include <Fastor/tensor_algebra/einsum.h>
#include <Fastor/tensor_algebra/indicial.h>
#include <cassert>
#include <cmath>
#include <iostream>
#include <ostream>
#include <tuple>
#include <unsupported/Eigen/CXX11/Tensor>

using namespace Eigen;
using namespace Fastor;

using Tensor1D = Fastor::Tensor< double, 3 >;
using Tensor2D = Fastor::Tensor< double, 3, 3 >;
using Tensor3D = Fastor::Tensor< double, 3, 3, 3 >;
using Tensor4D = Fastor::Tensor< double, 3, 3, 3, 3 >;

enum { a, i, b, j, k, l, m, n, I, J };
namespace Marmot::Materials {

  namespace InterfaceMaterialHelperFunctions {
    Tensor2D compute_inv( const Tensor2D& I, Tensor2D& Q )
    {
      // Solve G in batch mode using LU decomposition
      Eigen::Matrix3d Im( I.data() );
      Eigen::Matrix3d Qm( Q.data() );
      Eigen::Matrix3d G_mat = Qm.fullPivLu().solve( Im );

      // Store result back into G (reinterpret as a 3D Fastor tensor)
      Tensor2D G( 0 );
      Eigen::Map< Eigen::Matrix< double, 3, 3, Eigen::RowMajor > >( G.data() ) = G_mat;

      return G;
    }
    std::tuple< Tensor4D, const Tensor4D, Tensor4D, Tensor2D > interfaceGeometrySystemCouplings( const Tensor2D& I,
                                                                                                 const Tensor2D& N,
                                                                                                 const Tensor2D& T,
                                                                                                 const Tensor4D& L )
    {
      Tensor2D Q = Fastor::einsum< Fastor::Index< a, i, b, j >, Fastor::Index< i, j >, Fastor::OIndex< a, b > >( L, N );
      Tensor2D G = compute_inv( I, Q );

      Tensor4D A = Fastor::einsum< Fastor::Index< a, b >, Fastor::Index< i, j >, Fastor::OIndex< a, i, b, j > >( G, N );
      Tensor4D LA = Fastor::
        einsum< Fastor::Index< a, i, m, n >, Fastor::Index< m, n, b, j >, Fastor::OIndex< a, i, b, j > >( L, A );
      Tensor4D LAL = Fastor::
        einsum< Fastor::Index< a, i, m, n >, Fastor::Index< m, n, b, j >, Fastor::OIndex< a, i, b, j > >( LA, L );
      Tensor4D B = L - LAL;

      return std::make_tuple( B, L, A, G );
    }
    std::tuple< Tensor4D,
                Tensor4D,
                Tensor4D,
                Tensor4D,
                Tensor4D,
                Tensor4D,
                Tensor4D,
                Tensor4D,
                Tensor2D,
                Tensor2D,
                Tensor2D,
                Tensor4D,
                Tensor4D,
                Tensor4D >
    calculateFY( const Tensor2D& I,
                 const Tensor2D& N,
                 const Tensor2D& T,
                 const Tensor4D& C_0_aibj,
                 const Tensor4D& C_M_aibj,
                 const Tensor4D& C_I_aibj )
    {
      Tensor2D G_0;
      Tensor4D A_0;
      Tensor4D B_0;
      Tensor4D L_0;

      Tensor2D G_M;
      Tensor4D A_M;
      Tensor4D B_M;
      Tensor4D L_M;

      Tensor2D G_I;
      Tensor4D A_I;
      Tensor4D B_I;
      Tensor4D L_I;

      std::tie( B_0, L_0, A_0, G_0 ) = interfaceGeometrySystemCouplings( I, N, T, C_0_aibj );

      std::tie( B_M, L_M, A_M, G_M ) = interfaceGeometrySystemCouplings( I, N, T, C_M_aibj );

      std::tie( B_I, L_I, A_I, G_I ) = interfaceGeometrySystemCouplings( I, N, T, C_I_aibj );

      Tensor4D F = -2.0 * Fastor::einsum< Fastor::Index< a, i, m, n >,
                                          Fastor::Index< m, n, b, j >,
                                          Fastor::OIndex< a, i, b, j > >( A_0, L_0 );
      F += Fastor::
        einsum< Fastor::Index< a, i, m, n >, Fastor::Index< m, n, b, j >, Fastor::OIndex< a, i, b, j > >( A_M, L_M );
      F += Fastor::
        einsum< Fastor::Index< a, i, m, n >, Fastor::Index< m, n, b, j >, Fastor::OIndex< a, i, b, j > >( A_I, L_I );

      Tensor4D Y = Fastor::
        einsum< Fastor::Index< a, i, m, n >, Fastor::Index< m, n, b, j >, Fastor::OIndex< a, i, b, j > >( L_M, A_M );
      Y += Fastor::
        einsum< Fastor::Index< a, i, m, n >, Fastor::Index< m, n, b, j >, Fastor::OIndex< a, i, b, j > >( L_I, A_I );
      Y -= 2.0 * Fastor::einsum< Fastor::Index< a, i, m, n >,
                                 Fastor::Index< m, n, b, j >,
                                 Fastor::OIndex< a, i, b, j > >( L_0, A_0 );
      return std::make_tuple( F, Y, A_0, L_0, A_M, L_M, A_I, L_I, G_0, G_M, G_I, B_0, B_M, B_I );
    }

    std::tuple< Tensor4D, Tensor2D, Tensor3D, Tensor4D > calculateMaterialMatrices( const Tensor1D& normal,
                                                                                    const Tensor2D& I,
                                                                                    const Tensor2D& N,
                                                                                    const Tensor2D& T,
                                                                                    const Tensor4D& C_0_aibj,
                                                                                    const Tensor4D& C_M_aibj,
                                                                                    const Tensor4D& C_I_aibj )
    {

      auto [F, Y, A_0, L_0, A_M, L_M, A_I, L_I, G_0, G_M, G_I, B_0, B_M, B_I] = calculateFY( I,
                                                                                             N,
                                                                                             T,
                                                                                             C_0_aibj,
                                                                                             C_M_aibj,
                                                                                             C_I_aibj );

      assert_equivalent_F_Falt_Y( F, Y, A_0, L_0, A_M, L_M, A_I, L_I );

      Tensor2D H = 2.0 * G_0 - G_M - G_I;

      Tensor4D Z     = B_M + B_I - 2.0 * B_0;
      Tensor2D H_inv = compute_inv( I, H );

      Tensor3D
        nF = Fastor::einsum< Fastor::Index< a >, Fastor::Index< a, i, b, j >, Fastor::OIndex< i, b, j > >( normal, F );
      Tensor3D
        Fn = Fastor::einsum< Fastor::Index< a, i, b, j >, Fastor::Index< j >, Fastor::OIndex< a, i, b > >( F, normal );
      Tensor3D
        Yn = Fastor::einsum< Fastor::Index< a, i, b, j >, Fastor::Index< j >, Fastor::OIndex< a, i, b > >( Y, normal );
      Tensor3D
        H_inv_nF = Fastor::einsum< Fastor::Index< a, b >, Fastor::Index< b, i, j >, Fastor::OIndex< a, i, j > >( H_inv,
                                                                                                                 nF );
      Tensor4D Yn_H_inv_Fn = Fastor::einsum< Fastor::Index< a, i, m >,
                                             Fastor::Index< m, n >,
                                             Fastor::Index< n, b, j >,
                                             Fastor::OIndex< a, i, b, j > >( Yn, H_inv, Fn );
      return std::make_tuple( Z, H_inv, H_inv_nF, Yn_H_inv_Fn );
    }

    void assert_equivalent_F_Falt_Y( Tensor4D F,
                                     Tensor4D Y,
                                     Tensor4D A_0,
                                     Tensor4D L_0,
                                     Tensor4D A_M,
                                     Tensor4D L_M,
                                     Tensor4D A_I,
                                     Tensor4D L_I )
    {
      Tensor4D F_alt = -2.0 * Fastor::einsum< Fastor::Index< a, i, m, n >,
                                              Fastor::Index< m, n, b, j >,
                                              Fastor::OIndex< a, i, b, j > >( A_0, L_0 );
      F_alt += Fastor::
        einsum< Fastor::Index< a, i, m, n >, Fastor::Index< m, n, b, j >, Fastor::OIndex< a, i, b, j > >( A_M, L_M );
      F_alt += Fastor::
        einsum< Fastor::Index< a, i, m, n >, Fastor::Index< m, n, b, j >, Fastor::OIndex< a, i, b, j > >( A_I, L_I );

      double                                                         atol = 1e-8;
      Eigen::Map< const Eigen::Matrix< double, Eigen::Dynamic, 1 > > F_alt_flat( F_alt.data(), F_alt.size() );
      Eigen::Map< const Eigen::Matrix< double, Eigen::Dynamic, 1 > > F_flat( F.data(), F.size() );

      double max_diff = ( F_alt_flat - F_flat ).cwiseAbs().maxCoeff();
      if ( max_diff > atol ) {
        std::cerr << "Arrays are not equal within the tolerance." << std::endl;
        std::abort();
      }
    }

    // Convert 4th-order Fastor tensor (3x3x3x3) to Eigen 9x9 matrix
    Eigen::Matrix< double, 9, 9 > convert4thOrderTensorToMatrix( const Tensor4D& tensor )
    {
      Eigen::Matrix< double, 9, 9 > matrix( 9, 9 );

      for ( int i = 0; i < 3; ++i ) {
        for ( int j = 0; j < 3; ++j ) {
          for ( int k = 0; k < 3; ++k ) {
            for ( int l = 0; l < 3; ++l ) {
              int row            = 3 * i + j; // Convert (i, j) to single index
              int col            = 3 * k + l; // Convert (k, l) to single index
              matrix( row, col ) = tensor( i, j, k, l );
            }
          }
        }
      }

      return matrix;
    }

    // Convert 3rd-order Fastor tensor (3x3x3) to Eigen 9x3 matrix
    Eigen::Matrix< double, 9, 3 > convert3rdOrderTensorToMatrix( const Tensor3D& tensor )
    {
      Eigen::Matrix< double, 9, 3 > matrix( 9, 3 );

      for ( int i = 0; i < 3; ++i ) {
        for ( int j = 0; j < 3; ++j ) {
          for ( int k = 0; k < 3; ++k ) {
            int row          = 3 * i + j; // Convert (i, j) to single index
            matrix( row, k ) = tensor( i, j, k );
          }
        }
      }

      return matrix;
    }

    // Convert 2nd-order Fastor tensor (3x3) to Eigen 3x3 matrix
    Eigen::Matrix< double, 3, 3 > convert2ndOrderTensorToMatrix( const Tensor2D& tensor )
    {
      Eigen::Matrix< double, 3, 3 > matrix( 3, 3 );

      for ( int i = 0; i < 3; ++i ) {
        for ( int j = 0; j < 3; ++j ) {
          matrix( i, j ) = tensor( i, j );
        }
      }
      return matrix;
    }

    Tensor4D voigtToStiffness( const Eigen::Matrix< double, 6, 6 >& voigtStiffness )
    {
      using namespace Marmot::ContinuumMechanics::TensorUtility::IndexNotation;

      Fastor::Tensor< double, 3, 3, 3, 3 > stiffness;
      stiffness.zeros(); // Set to zero

      int row, col;
      for ( int i = 0; i < 3; ++i )
        for ( int j = 0; j < 3; ++j ) {
          row = toVoigt< 3 >( i, j );
          for ( int k = 0; k < 3; ++k )
            for ( int l = 0; l < 3; ++l ) {
              col = toVoigt< 3 >( k, l );
              stiffness( i, j, k, l ) += voigtStiffness( row, col );
            }
        }
      return stiffness;
    }

    std::tuple< Tensor4D, Tensor2D, Tensor3D, Tensor4D > calculateInterfaceMaterialParameters( const Tensor1D& normal,
                                                                                               const double&   E_M,
                                                                                               const double&   nu_M,
                                                                                               const double&   E_I,
                                                                                               const double&   nu_I,
                                                                                               const double&   E_0,
                                                                                               const double&   nu_0 )
    {
      using namespace Marmot::ContinuumMechanics::Elasticity::Isotropic;
      Eigen::Matrix< double, 6, 6 > C_M_voigt_full = stiffnessTensor( E_M, nu_M );
      Eigen::Matrix< double, 6, 6 > C_I_voigt_full = stiffnessTensor( E_I, nu_I );
      Eigen::Matrix< double, 6, 6 > C_0_voigt_full = stiffnessTensor( E_0, nu_0 );

      Tensor2D I = { { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } };

      Tensor2D N = Fastor::einsum< Fastor::Index< i >, Fastor::Index< j >, Fastor::OIndex< i, j > >( normal, normal );

      Tensor2D T = I - N;

      Tensor4D C_M_aibj = voigtToStiffness( C_M_voigt_full );
      Tensor4D C_I_aibj = voigtToStiffness( C_I_voigt_full );
      Tensor4D C_0_aibj = voigtToStiffness( C_0_voigt_full );

      auto [Z,
            H_inv,
            H_inv_nF,
            Yn_H_inv_Fn] = calculateMaterialMatrices( normal, I, N, T, C_0_aibj, C_M_aibj, C_I_aibj );

      return { Z, H_inv, H_inv_nF, Yn_H_inv_Fn };
    }

    std::tuple< Tensor2D, Tensor4D, Eigen::Matrix< double, 3, 3 >, Eigen::Matrix< double, 9, 9 > > calculateEffectiveProperties(
      double&         zerothWienertStiffness_Ru,
      double&         creep_Ru_stiffness,
      double&         zerothWienertStiffness_Rs,
      double&         creep_Rs_stiffness,
      const Tensor1D& normal,
      const double&   E_M,
      const double&   nu_M,
      const double&   E_I,
      const double&   nu_I,
      const double&   E_0,
      const double&   nu_0 )
    {
      using namespace Marmot::ContinuumMechanics::Elasticity::Isotropic;
      Tensor2D N = Fastor::einsum< Fastor::Index< i >, Fastor::Index< j >, Fastor::OIndex< i, j > >( normal, normal );

      double E_bar = E_M + E_I - 2. * E_0;
      double H_bar = 2. * 1. / E_0 - 1. / E_M - 1 / E_I;

      double                        nu_bar               = nu_0;
      Eigen::Matrix< double, 6, 6 > unitC_bar_voigt_full = stiffnessTensor( 1, nu_bar );
      Tensor4D                      unitC_bar_tensor     = voigtToStiffness( unitC_bar_voigt_full );

      Tensor2D unitQ_bar_tensor = Fastor::
        einsum< Fastor::Index< a, i, b, j >, Fastor::Index< i, j >, Fastor::OIndex< a, b > >( unitC_bar_tensor, N );
      Eigen::Map< Eigen::Matrix< double, 3, 3, Eigen::RowMajor > > unitQ_bar_voigt_full( unitQ_bar_tensor.data() );

      Eigen::Matrix< double, 3, 3 >     II                   = Eigen::Matrix< double, 3, 3 >::Identity();
      Eigen::Matrix< double, 3, 3 >     unitG_bar_voigt_full = unitQ_bar_voigt_full.fullPivLu().solve( II );
      Fastor::TensorMap< double, 3, 3 > unitG_bar_tensor( unitG_bar_voigt_full.data() );

      Tensor4D unitZ_bar_tensor = unitC_bar_tensor;
      Tensor4D unitA_bar_tensor = Fastor::
        einsum< Fastor::Index< a, b >, Fastor::Index< i, j >, Fastor::OIndex< a, i, b, j > >( unitG_bar_tensor, N );
      Tensor4D unitCA_tensor = Fastor::einsum< Fastor::Index< a, i, m, n >,
                                               Fastor::Index< m, n, k, l >,
                                               Fastor::OIndex< a, i, k, l > >( unitC_bar_tensor, unitA_bar_tensor );
      unitZ_bar_tensor -= Fastor::einsum< Fastor::Index< m, n, k, l >,
                                          Fastor::Index< k, l, b, j >,
                                          Fastor::OIndex< m, n, b, j > >( unitCA_tensor, unitC_bar_tensor );
      Eigen::Map< Eigen::Matrix< double, 9, 9, Eigen::RowMajor > > unitZ_bar_voigt_full( unitZ_bar_tensor.data() );

      double barE_Ru = 1. / H_bar + 0. * zerothWienertStiffness_Ru + creep_Ru_stiffness;
      double barE_Rs = E_bar + 0. * zerothWienertStiffness_Rs - creep_Rs_stiffness;

      Tensor2D H_inv_Ru_tensor = barE_Ru * unitQ_bar_tensor;

      Tensor4D Z_Rs_tensor = barE_Rs * unitZ_bar_tensor;

      return { H_inv_Ru_tensor, Z_Rs_tensor, unitQ_bar_voigt_full, unitZ_bar_voigt_full };
    }
  } // namespace InterfaceMaterialHelperFunctions
} // namespace Marmot::Materials