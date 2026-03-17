#include "Marmot/VonMisesInterface.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialHypoElasticInterface.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotUtility.h"
#include "Marmot/MarmotViscoelasticity.h"
#include "Marmot/MarmotVoigt.h"
#include "Marmot/MarmotWiechertInterface.h"
#include "Marmot/VonMises.h"

#include "Fastor/Fastor.h"
#include "Marmot/MarmotInterfaceMaterialHelperFunctions.h"
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <Fastor/expressions/linalg_ops/unary_norm_op.h>

#include "autodiff/forward/real.hpp"
#include <iostream>
#include <map>
#include <string>

using namespace Marmot;
using namespace Eigen;

using Tensor1D = Fastor::Tensor< double, 3 >;
using Tensor2D = Fastor::Tensor< double, 3, 3 >;
using Tensor3D = Fastor::Tensor< double, 3, 3, 3 >;
using Tensor4D = Fastor::Tensor< double, 3, 3, 3, 3 >;

namespace Marmot::Materials {

  VonMisesInterface::VonMisesInterface( const double* materialProperties, int nMaterialProperties, int materialNumber )
    : MarmotMaterialHypoElasticInterface( materialProperties, nMaterialProperties, materialNumber ),
      // clang-format off
      // elasticity parameters
      E_0( materialProperties[0] ),
      nu_0( materialProperties[1] ),
      h( materialProperties[2] ),
      // plasticity parameters
      yieldStress( materialProperties[3] ),
      HLin( materialProperties[4] ),
      deltaYieldStress( materialProperties[5] ),
      delta( materialProperties[6] ),
      // Re-map properties for VonMisesModel: [E, nu, yieldStress, HLin, deltaYieldStress, delta]
      // (skip h at index 2)
      vonMisesProps{ materialProperties[0], materialProperties[1],
                     materialProperties[3], materialProperties[4],
                     materialProperties[5], materialProperties[6] },
      // Instantiate VonMisesModel once using the re-mapped properties stored in vonMisesProps
      vonMisesModel( vonMisesProps.data(), 6, materialNumber )
  // clang-format on
  {
  }
  void VonMisesInterface::computeStress( double*       force,
                                         double*       averageStress,
                                         double*       H_inv_ij,
                                         double*       Z_ijkl,
                                         double*       H_inv_nF_ijk,
                                         double*       Yn_H_inv_Fn_ijkl,
                                         const double* dU,
                                         const double* dSurfaceDispGradient,
                                         const double* normal,
                                         const double* timeOld,
                                         const double  dT,
                                         double&       pNewDT )
  {
    using namespace Marmot::Materials::InterfaceMaterialHelperFunctions;
    enum { i, j, k, l };

    // map to force, surface stress, displacement, surface strain, normal and tangent stiffness
    // use Fastor because we really need to use the einsum

    Fastor::Tensor< double, 3 >          forceFtensor( force );
    Fastor::Tensor< double, 3, 3 >       averageStressFtensor( averageStress );
    Fastor::Tensor< double, 3, 3 >       H_inv_ij_Ftensor( H_inv_ij );
    Fastor::Tensor< double, 3, 3, 3, 3 > Z_ijkl_Ftensor( Z_ijkl );
    Fastor::Tensor< double, 3, 3, 3 >    H_inv_nF_ijk_Ftensor( H_inv_nF_ijk );
    Fastor::Tensor< double, 3, 3, 3, 3 > Yn_H_inv_Fn_ijkl_Ftensor( Yn_H_inv_Fn_ijkl );
    Fastor::Tensor< double, 6, 1 >       dUFtensor( dU );
    Fastor::Tensor< double, 18, 1 >      dSurfaceDispGradientFtensor( dSurfaceDispGradient );
    Fastor::Tensor< double, 3 >          normalFtensor( normal );

    // handle zero strain increment: C_ep_voigt is already initialized to Cel in the constructor
    if ( Fastor::norm( dUFtensor ) < 1e-14 && Fastor::norm( dSurfaceDispGradientFtensor ) < 1e-14 && timeOld == 0 ) {

      pNewDT = 1.0;

      auto [unitZ_ijkl,
            unitH_inv_ij,
            unitH_inv_nF_ijk,
            unitYn_H_inv_Fn_ijkl] = calculateInterfaceMaterialParameters( normalFtensor, nu_0 );

      Z_ijkl_Ftensor           = h * E_0 * unitZ_ijkl;
      Yn_H_inv_Fn_ijkl_Ftensor = h * E_0 * unitYn_H_inv_Fn_ijkl;
      H_inv_ij_Ftensor         = 1. / h * E_0 * unitH_inv_ij;
      H_inv_nF_ijk_Ftensor     = E_0 * unitH_inv_nF_ijk;

      std::copy( H_inv_ij_Ftensor.data(), H_inv_ij_Ftensor.data() + 9, H_inv_ij );
      std::copy( Z_ijkl_Ftensor.data(), Z_ijkl_Ftensor.data() + 81, Z_ijkl );
      std::copy( H_inv_nF_ijk_Ftensor.data(), H_inv_nF_ijk_Ftensor.data() + 27, H_inv_nF_ijk );
      std::copy( Yn_H_inv_Fn_ijkl_Ftensor.data(), Yn_H_inv_Fn_ijkl_Ftensor.data() + 81, Yn_H_inv_Fn_ijkl );
      return;
    }

    // Evaluate average stress on the layer using Von Mises yield criterion.
    // vonMisesModel.computeStress updates averageStress and writes the new
    // elastoplastic tangent into C_ep (state var), which persists across increments.

    // Displacement jump: top(0:3) - bottom(3:6)
    Fastor::Tensor< double, 3 > dJumpU = dUFtensor( Fastor::seq( 0, 3 ), 0 ) -
                                         dUFtensor( Fastor::seq( 3, Fastor::last ), 0 );

    // Average surface strain: 0.5*(top(0:9) + bottom(9:18)), reshaped to 3x3
    Fastor::Tensor< double, 9, 1 >
         dSurfaceDispGradientAvgFlat = 0.5 * ( dSurfaceDispGradientFtensor( Fastor::seq( 0, 9 ), 0 ) +
                                            dSurfaceDispGradientFtensor( Fastor::seq( 9, Fastor::last ), 0 ) );
    auto dSurfaceDispGradientAvg     = Fastor::Tensor< double, 3, 3 >(
      Fastor::reshape< 3, 3 >( dSurfaceDispGradientAvgFlat ) );
    // Average displacement gradient: jump contribution (normal-to-layer) + surface strain
    Fastor::Tensor< double, 3, 3 >
      dU_kl_Jump = ( 1. / h ) *
                   Fastor::einsum< Fastor::Index< i >, Fastor::Index< j >, Fastor::OIndex< i, j > >( dJumpU,
                                                                                                     normalFtensor );

    // Symmetrize and include surface strain to obtain the full average strain increment
    Fastor::Tensor< double, 3, 3 > dDispGradAvg = ( dU_kl_Jump + dSurfaceDispGradientAvg );
    Fastor::Tensor< double, 3, 3 > dStrainAvg   = 0.5 * ( dDispGradAvg + Fastor::transpose( dDispGradAvg ) );

    // Convert 3x3 strain tensor to Voigt 6-vector (with factor 2 on shear components)
    Eigen::Map< const Eigen::Matrix< double, 3, 3, Eigen::RowMajor > > dStrainAvgEigen( dStrainAvg.data() );
    const Marmot::Vector6d dStrainAvgVoigt = Marmot::ContinuumMechanics::VoigtNotation::strainToVoigt(
      dStrainAvgEigen );

    // VonMisesModel updates stress in-place (incremental hypoelastic-plastic);
    // read the current 3x3 averageStress, symmetrize, convert to Voigt 6-vector
    Eigen::Map< const Eigen::Matrix< double, 3, 3, Eigen::RowMajor > > averageStressCurrent( averageStress );
    const Eigen::Matrix< double, 3, 3 >                                averageStressSym = 0.5 *
                                                           ( averageStressCurrent + averageStressCurrent.transpose() );
    Marmot::Vector6d averageStressVoigt = Marmot::ContinuumMechanics::VoigtNotation::stressToVoigt( averageStressSym );

    auto& C_ep = managedStateVars->C_ep_voigt;
    vonMisesModel.computeStress( averageStressVoigt.data(), C_ep.data(), dStrainAvgVoigt.data(), timeOld, dT, pNewDT );

    // Expand the updated Voigt 6-vector back to full 3x3 and write into averageStress for the FE model
    // voigtToStress returns column-major Eigen matrix; convert to row-major before copying
    // so that Fastor (row-major) reads the layout correctly
    const Eigen::Matrix< double, 3, 3, Eigen::RowMajor >
      averageStressFull = h * Marmot::ContinuumMechanics::VoigtNotation::voigtToStress( averageStressVoigt );
    std::copy( averageStressFull.data(), averageStressFull.data() + 9, averageStress );

    // Reload Fastor tensor from the updated 3x3 buffer
    averageStressFtensor = Fastor::Tensor< double, 3, 3 >( averageStress );

    auto [Z_ijkl_ep,
          H_inv_ij_ep,
          H_inv_nF_ijk_ep,
          Yn_H_inv_Fn_ijkl_ep] = calculateInterfaceMaterialParameters( normalFtensor, C_ep );

    Fastor::Tensor< double, 3, 3 >       H_inv_ij_Ftensor_scaled = ( 1.0 / h ) * H_inv_ij_ep;
    Fastor::Tensor< double, 3, 3, 3, 3 > Z_ijkl_Ftensor_scaled   = (h)*Z_ijkl_ep;
    Fastor::Tensor< double, 3, 3, 3 >    H_inv_nF_Ftensor_scaled = H_inv_nF_ijk_ep;
    Fastor::Tensor< double, 3, 3, 3, 3 > Yn_H_Fn_Ftensor_scaled  = h * Yn_H_inv_Fn_ijkl_ep;

    forceFtensor = 1. / h *
                   Fastor::einsum< Fastor::Index< i, j >,
                                   Fastor::Index< j >,
                                   Fastor::OIndex< i > >( averageStressFtensor, normalFtensor );

    std::copy( forceFtensor.data(), forceFtensor.data() + 3, force );

    std::copy( H_inv_ij_Ftensor_scaled.data(), H_inv_ij_Ftensor_scaled.data() + 9, H_inv_ij );
    std::copy( Z_ijkl_Ftensor_scaled.data(), Z_ijkl_Ftensor_scaled.data() + 81, Z_ijkl );
    std::copy( H_inv_nF_Ftensor_scaled.data(), H_inv_nF_Ftensor_scaled.data() + 27, H_inv_nF_ijk );
    std::copy( Yn_H_Fn_Ftensor_scaled.data(), Yn_H_Fn_Ftensor_scaled.data() + 81, Yn_H_inv_Fn_ijkl );

    return;
  };
  void VonMisesInterface::assignStateVars( double* stateVars, int nStateVars )
  {
    if ( nStateVars < getNumberOfRequiredStateVars() )
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": Not sufficient stateVars!" );

    managedStateVars = std::make_unique< VonMisesInterfaceStateVarManager >( stateVars );

    // Also assign the kappa state var pointer to vonMisesModel so it shares the same memory
    vonMisesModel.assignStateVars( &managedStateVars->kappa, 1 );

    // If C_ep_voigt state var is still zero (first ever assignment), initialize it to elastic stiffness
    if ( managedStateVars->C_ep_voigt.isZero() )
      managedStateVars->C_ep_voigt = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E_0, nu_0 );

    return MarmotMaterialHypoElasticInterface::assignStateVars( stateVars, nStateVars );
  }

  StateView VonMisesInterface::getStateView( const std::string& stateName )
  {
    return managedStateVars->getStateView( stateName );
  }

  double VonMisesInterface::getDensity()
  {
    if ( this->nMaterialProperties < 7 )
      throw std::runtime_error( MakeString() << __PRETTY_FUNCTION__ << ": No density given! nMaterialProperties < 7" );
    return this->materialProperties[6];
  }

} // namespace Marmot::Materials
