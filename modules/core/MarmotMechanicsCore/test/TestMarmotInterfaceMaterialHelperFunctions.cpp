#include "Fastor/Fastor.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotInterfaceMaterialHelperFunctions.h"
#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
#include <functional>

using namespace Marmot::Testing;
using namespace Marmot::ContinuumMechanics::Elasticity::Isotropic;
using namespace Marmot::Materials::InterfaceMaterialHelperFunctions;
using namespace Marmot;
// #include "Marmot/MarmotVoigt.h"

// Define Fastor index variables
// stexpr int i = 0, j = 1, k = 2, l = 3, m = 4, n = 5, a = 6, b = 7;
enum { a, i, b, j, k, l, m, n, I, J };

void interfaceGeometrySystemCouplingsTestFunction()
{
  // Test function implementation
  Tensor2D I = { { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } };

  Tensor1D normal = { 1.0, 0.0, 0.0 };

  Tensor2D N = Fastor::einsum< Fastor::Index< i >, Fastor::Index< j >, Fastor::OIndex< i, j > >( normal, normal );

  Tensor2D T = I - N;
  // Evaluate for a simple case (nu=0) without intercoupling terms that complicate the evaluation of B
  double E  = 1.0;                      // Young's modulus
  double nu = 0.;                       // Poisson's ratio
  double mu = E / ( 2. * ( 1. + nu ) ); // Shear modulus
  // Construct the stiffness matrix in Voigt notation using Marmot functionality
  const Eigen::Matrix< double, 6, 6 > C_voigt_full = stiffnessTensor( E, nu );
  // Construct the stiffness matrix in tensorial form
  Tensor4D C_aibj = Marmot::Materials::InterfaceMaterialHelperFunctions::voigtToStiffness( C_voigt_full );
  // Construct Q_ij = C_ijkl*N_kl
  Tensor2D Q_ik = Fastor::einsum< Fastor::Index< i, j, k, l >, Fastor::Index< j, l >, Fastor::OIndex< i, k > >( C_aibj,
                                                                                                                N );
  // Call the interface geometry function that constructs the principal B, G interface material matrices for each layer
  auto [B_ijkl,
        C_ijkl,
        A_ijkl,
        G_ik] = Marmot::Materials::InterfaceMaterialHelperFunctions::interfaceGeometrySystemCouplings( I,
                                                                                                       N,
                                                                                                       T,
                                                                                                       C_aibj );

  // Construct the analytical expression for Q,G,A,tensors
  Tensor4D I1_ijkl = Fastor::einsum< Fastor::Index< i, k >, Fastor::Index< j, l >, Fastor::OIndex< i, j, k, l > >( I,
                                                                                                                   I );
  Tensor4D I2_ijkl = Fastor::einsum< Fastor::Index< i, l >, Fastor::Index< j, k >, Fastor::OIndex< i, j, k, l > >( I,
                                                                                                                   I );
  Tensor4D IB_ijkl = I1_ijkl + I2_ijkl;
  Tensor2D
    Q_calc_ik = mu *
                Fastor::einsum< Fastor::Index< i, j, k, l >, Fastor::Index< j, l >, Fastor::OIndex< i, k > >( IB_ijkl,
                                                                                                              N );

  double diff_Q_Q_calc = Fastor::norm( ( Q_ik - Q_calc_ik ) );
  throwExceptionOnFailure( checkIfEqual( diff_Q_Q_calc, 0.0, 1e-10 ), "error in the Q tensor calculation" );

  Tensor2D G_calc_ik     = 1 / mu * I - 1 / ( 2 * mu ) * N;
  double   diff_G_G_calc = Fastor::norm( ( G_ik - G_calc_ik ) );
  throwExceptionOnFailure( checkIfEqual( diff_G_G_calc, 0.0, 1e-10 ), "error in the G tensor calculation" );

  // Tensor4D A_ijkl = Fastor::einsum<Fastor::Index<i,k>, Fastor::Index<j,l>, Fastor::OIndex<i,j,k,l>>(G_ik,N);
  Tensor4D
    A_calc_ijkl = ( 1 / mu ) *
                  Fastor::einsum< Fastor::Index< i, k >, Fastor::Index< j, l >, Fastor::OIndex< i, j, k, l > >( I, N );
  A_calc_ijkl -= ( 1 / ( 2 * mu ) ) *
                 Fastor::einsum< Fastor::Index< i, k >, Fastor::Index< j, l >, Fastor::OIndex< i, j, k, l > >( N, N );
  double diff_A_A_calc = Fastor::norm( ( A_ijkl - A_calc_ijkl ) );
  throwExceptionOnFailure( checkIfEqual( diff_A_A_calc, 0.0, 1e-10 ), "error in the A tensor calculation" );

  // Calculate analytical form of the B tensor
  Tensor4D Ia_ijkl = Fastor::einsum< Fastor::Index< i, k >, Fastor::Index< j, l >, Fastor::OIndex< i, j, k, l > >( I,
                                                                                                                   N );
  Tensor4D Ib_ijkl = -( 1. / 2. ) *
                     Fastor::einsum< Fastor::Index< i, k >, Fastor::Index< j, l >, Fastor::OIndex< i, j, k, l > >( N,
                                                                                                                   N );
  Tensor4D Ic_ijkl = Fastor::einsum< Fastor::Index< j, k >, Fastor::Index< i, l >, Fastor::OIndex< i, j, k, l > >( I,
                                                                                                                   N );
  Tensor4D Id_ijkl = -( 1. / 2. ) *
                     Fastor::einsum< Fastor::Index< j, k >, Fastor::Index< i, l >, Fastor::OIndex< i, j, k, l > >( N,
                                                                                                                   N );
  Tensor4D Ie_ijkl = Fastor::einsum< Fastor::Index< i, l >, Fastor::Index< j, k >, Fastor::OIndex< i, j, k, l > >( I,
                                                                                                                   N );
  Tensor4D If_ijkl = -( 1. / 2. ) *
                     Fastor::einsum< Fastor::Index< i, l >, Fastor::Index< j, k >, Fastor::OIndex< i, j, k, l > >( N,
                                                                                                                   N );
  Tensor4D Ig_ijkl = Fastor::einsum< Fastor::Index< j, l >, Fastor::Index< i, k >, Fastor::OIndex< i, j, k, l > >( I,
                                                                                                                   N );
  Tensor4D Ih_ijkl = -( 1. / 2. ) *
                     Fastor::einsum< Fastor::Index< j, l >, Fastor::Index< i, k >, Fastor::OIndex< i, j, k, l > >( N,
                                                                                                                   N );

  Tensor4D IA_ijkl = Ia_ijkl + Ib_ijkl + Ic_ijkl + Id_ijkl + Ie_ijkl + If_ijkl + Ig_ijkl + Ih_ijkl;

  Tensor4D B_calc_ijkl = mu * IB_ijkl - mu * IA_ijkl;

  double diff_B_B_calc = Fastor::norm( ( B_ijkl - B_calc_ijkl ) );
  throwExceptionOnFailure( checkIfEqual( diff_B_B_calc, 0.0, 1e-10 ), "error in the B tensor calculation" );
}

void calculateMaterialMatricesTestFunction()
{
  using namespace Marmot::ContinuumMechanics::Elasticity::Isotropic;
  // Define material properties for the interphase
  const double E_M  = 1.0; // Young's modulus top layer
  const double nu_M = 0.3; // Poisson's ratio top layer
  const double E_I  = 1.0; // Young's modulus bottom layer
  const double nu_I = 0.3; // Poisson's ratio bottom layer
  const double E_0  = 1e2; // Young's modulus interphase
  const double nu_0 = 0.3; // Poisson's ratio interphase

  Eigen::Matrix< double, 6, 6 > C_M_voigt_full = stiffnessTensor( E_M, nu_M );
  Eigen::Matrix< double, 6, 6 > C_I_voigt_full = stiffnessTensor( E_I, nu_I );
  Eigen::Matrix< double, 6, 6 > C_0_voigt_full = stiffnessTensor( E_0, nu_0 );

  Tensor2D I      = { { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } };
  Tensor1D normal = { 1.0, 0.0, 0.0 };

  Tensor2D N = Fastor::einsum< Fastor::Index< i >, Fastor::Index< j >, Fastor::OIndex< i, j > >( normal, normal );

  Tensor2D T = I - N;

  Tensor4D C_M_aibj = voigtToStiffness( C_M_voigt_full );
  Tensor4D C_I_aibj = voigtToStiffness( C_I_voigt_full );
  Tensor4D C_0_aibj = voigtToStiffness( C_0_voigt_full );

  auto [F, Y, A_0, L_0, A_M, L_M, A_I, L_I, G_0, G_M, G_I, B_0, B_M, B_I] = calculateFY( I,
                                                                                         N,
                                                                                         T,
                                                                                         C_0_aibj,
                                                                                         C_M_aibj,
                                                                                         C_I_aibj );

  Tensor4D F_alt = -2.0 * Fastor::einsum< Fastor::Index< a, i, m, n >,
                                          Fastor::Index< m, n, b, j >,
                                          Fastor::OIndex< a, i, b, j > >( A_0, L_0 );
  F_alt += Fastor::
    einsum< Fastor::Index< a, i, m, n >, Fastor::Index< m, n, b, j >, Fastor::OIndex< a, i, b, j > >( A_M, L_M );
  F_alt += Fastor::
    einsum< Fastor::Index< a, i, m, n >, Fastor::Index< m, n, b, j >, Fastor::OIndex< a, i, b, j > >( A_I, L_I );

  Eigen::Map< const Eigen::Matrix< double, Eigen::Dynamic, 1 > > F_alt_flat( F_alt.data(), F_alt.size() );
  Eigen::Map< const Eigen::Matrix< double, Eigen::Dynamic, 1 > > F_flat( F.data(), F.size() );

  double max_diff = ( F_alt_flat - F_flat ).cwiseAbs().maxCoeff();
  throwExceptionOnFailure( checkIfEqual( max_diff, 0.0, 1e-10 ), "error in the F tensor calculation" );
}

void voigtToStiffnessTestFunction()
{
  const double E_M  = 1.0; // Young's modulus top layer
  const double nu_M = 0.3; // Poisson's ratio top layer

  Eigen::Matrix< double, 6, 6 > C_M_voigt_full = stiffnessTensor( E_M, nu_M );

  Tensor4D C_M_ijkl = voigtToStiffness( C_M_voigt_full );

  // Compare with alternative description of the stiffness tensor
  // Conversion tensors V_{Jkl} and W_{ijI}
  Fastor::TensorMap< const double, 6, 6 > Cvoigt( C_M_voigt_full.data() );
  Fastor::Tensor< double, 6, 3, 3 >       V;
  V.zeros();
  Fastor::Tensor< double, 3, 3, 6 > W;
  W.zeros();

  // V: maps symmetric tensor -> Voigt strain
  V( 0, 0, 0 ) = 1.0;                // 11
  V( 1, 1, 1 ) = 1.0;                // 22
  V( 2, 2, 2 ) = 1.0;                // 33
  V( 3, 1, 2 ) = V( 3, 2, 1 ) = 1.0; // 23
  V( 4, 0, 2 ) = V( 4, 2, 0 ) = 1.0; // 13
  V( 5, 0, 1 ) = V( 5, 1, 0 ) = 1.0; // 12

  // W: maps Voigt stress -> symmetric tensor
  W( 0, 0, 0 ) = 1.0;                // 11
  W( 1, 1, 1 ) = 1.0;                // 22
  W( 2, 2, 2 ) = 1.0;                // 33
  W( 1, 2, 3 ) = W( 2, 1, 3 ) = 1.0; // 23
  W( 0, 2, 4 ) = W( 2, 0, 4 ) = 1.0; // 13
  W( 0, 1, 5 ) = W( 1, 0, 5 ) = 1.0; // 12

  Tensor4D stiffness = einsum< Fastor::Index< i, j, I >,
                               Fastor::Index< I, J >,
                               Fastor::Index< J, k, l >,
                               Fastor::OIndex< i, j, k, l > >( W, Cvoigt, V );
  // Compare the computed stiffness tensor with the expected one
  double normDiff = Fastor::norm( C_M_ijkl - stiffness );
  throwExceptionOnFailure( checkIfEqual( normDiff, 0.0, 1e-10 ),
                           "error in the stiffness tensor conversion from Voigt notation" );
}

int main()
{
  auto tests = std::vector< std::function< void() > >{ interfaceGeometrySystemCouplingsTestFunction,
                                                       calculateMaterialMatricesTestFunction,
                                                       voigtToStiffnessTestFunction };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
