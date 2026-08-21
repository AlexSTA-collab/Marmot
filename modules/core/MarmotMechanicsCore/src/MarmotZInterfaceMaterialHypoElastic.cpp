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

#include "Marmot/MarmotZInterfaceMaterialHypoElastic.h"

#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Marmot;

namespace {

  using Vector3d           = Eigen::Matrix< double, 3, 1 >;
  using Vector9d           = Eigen::Matrix< double, 9, 1 >;
  using Vector21d          = Eigen::Matrix< double, 21, 1 >;
  using Matrix3dRowMajor   = Eigen::Matrix< double, 3, 3, Eigen::RowMajor >;
  using Matrix9dRowMajor   = Eigen::Matrix< double, 9, 9, Eigen::RowMajor >;
  using Matrix9x3RowMajor  = Eigen::Matrix< double, 9, 3, Eigen::RowMajor >;
  using Matrix9x21RowMajor = Eigen::Matrix< double, 9, 21, Eigen::RowMajor >;
  using Matrix21x3RowMajor = Eigen::Matrix< double, 21, 3, Eigen::RowMajor >;
  using Matrix3x21RowMajor = Eigen::Matrix< double, 3, 21, Eigen::RowMajor >;
  using Matrix21dRowMajor  = Eigen::Matrix< double, 21, 21, Eigen::RowMajor >;

  constexpr int flatIndex( int i, int j )
  {
    return 3 * i + j;
  }

  /** Identical to the corresponding helper in the sibling interface kernels:
   * lifts a Voigt tangent to the full 9x9 gradient tangent d vec(sigma) / d vec(G). */
  Matrix9dRowMajor fullGradientTangent( const Matrix6d& tangentVoigt )
  {
    Matrix9dRowMajor tangentFull = Matrix9dRowMajor::Zero();

    for ( int k = 0; k < 3; ++k ) {
      for ( int l = 0; l < 3; ++l ) {
        Matrix3dRowMajor dGradient = Matrix3dRowMajor::Zero();
        dGradient( k, l )          = 1.0;

        const Matrix3dRowMajor dStrain      = 0.5 * ( dGradient + dGradient.transpose() );
        const Vector6d         dStrainVoigt = ContinuumMechanics::VoigtNotation::strainToVoigt( dStrain );
        const Vector6d         dStressVoigt = tangentVoigt * dStrainVoigt;
        const Eigen::Matrix3d  dStress      = ContinuumMechanics::VoigtNotation::voigtToStress( dStressVoigt );

        const int column = flatIndex( k, l );
        for ( int i = 0; i < 3; ++i )
          for ( int j = 0; j < 3; ++j )
            tangentFull( flatIndex( i, j ), column ) = dStress( i, j );
      }
    }

    return tangentFull;
  }

  struct InterfaceGeometry {
    double   normalSeparation;
    Vector3d tangentialSeparation;
  };

  InterfaceGeometry evaluateInterfaceGeometry( const Vector3d& normal,
                                               const Vector3d& separationVector,
                                               double          constitutiveThickness )
  {
    constexpr double tolerance = 1.0e-12;

    if ( constitutiveThickness <= 0.0 )
      throw std::invalid_argument( "MarmotZInterfaceMaterialHypoElastic: interface thickness h must be positive." );

    if ( separationVector.norm() <= tolerance )
      return { constitutiveThickness, Vector3d::Zero() };

    const double normalSeparation = separationVector.dot( normal );
    if ( normalSeparation <= tolerance )
      throw std::invalid_argument(
        "MarmotZInterfaceMaterialHypoElastic: the top-bottom connector must have a positive normal component." );

    return { normalSeparation, separationVector - normalSeparation * normal };
  }

  /** Maps a 3-vector g to vec(g (x) n) / normalSeparation. With
   * normalSeparation = 1 this is B_n; with normalSeparation = ell it is B_w. */
  Matrix9x3RowMajor jumpToGradient( const Vector3d& normal, double normalSeparation )
  {
    Matrix9x3RowMajor map = Matrix9x3RowMajor::Zero();
    for ( int k = 0; k < 3; ++k )
      for ( int l = 0; l < 3; ++l )
        map( flatIndex( k, l ), k ) = normal( l ) / normalSeparation;
    return map;
  }

  /** Maps vec(A) to vec(A) - (1/normalSeparation) vec((A d_tau) (x) n), i.e. I_A + B_tau. */
  Matrix9dRowMajor surfaceToGradient( const Vector3d& normal,
                                      const Vector3d& tangentialSeparation,
                                      double          normalSeparation )
  {
    Matrix9dRowMajor map = Matrix9dRowMajor::Zero();

    for ( int k = 0; k < 3; ++k ) {
      for ( int l = 0; l < 3; ++l ) {
        const int row = flatIndex( k, l );
        for ( int b = 0; b < 3; ++b ) {
          const double identityPart     = l == b ? 1.0 : 0.0;
          map( row, flatIndex( k, b ) ) = identityPart - tangentialSeparation( b ) * normal( l ) / normalSeparation;
        }
      }
    }

    return map;
  }

  /** One side's constitutive response, evaluated from a scratch copy of the
   * OLD (committed) state -- the persistent state is written only on commit. */
  struct SideTrial {
    Matrix3dRowMajor      stress;
    Matrix9dRowMajor      CFull;
    std::vector< double > trialStateVars;
  };

  SideTrial evaluateSideTrial( MarmotMaterialHypoElastic&                 material,
                               const double*                              oldStateVars,
                               int                                        nStateVars,
                               const Matrix3dRowMajor&                    stressCurrent,
                               const Vector6d&                            strainIncrementVoigt,
                               const MarmotMaterialHypoElastic::timeInfo& timeInfo )
  {
    SideTrial trial;
    trial.trialStateVars.assign( oldStateVars, oldStateVars + nStateVars );

    const Eigen::Matrix3d stressCurrentSym( 0.5 * ( stressCurrent + stressCurrent.transpose() ) );
    const Vector6d        stressVoigt = ContinuumMechanics::VoigtNotation::stressToVoigt( stressCurrentSym );

    Matrix6d                           tangentVoigt = Matrix6d::Zero();
    MarmotMaterialHypoElastic::state3D baseState{ stressVoigt, 0.0, 0.0, trial.trialStateVars.data() };

    material.computeStress( baseState, tangentVoigt, strainIncrementVoigt, timeInfo );

    trial.stress = ContinuumMechanics::VoigtNotation::voigtToStress( baseState.stress );
    trial.CFull  = fullGradientTangent( tangentVoigt );
    return trial;
  }

  double defaultRegularization()
  {
    static const double value = []() {
      const char* e = std::getenv( "MARMOT_ZIFACE_REG" );
      return e ? std::atof( e ) : 1.0e-6;
    }();
    return value;
  }

} // namespace

MarmotZInterfaceMaterialHypoElastic::MarmotZInterfaceMaterialHypoElastic( const std::string& materialName,
                                                                          const double*      matProperties_,
                                                                          int                nMaterialProperties_,
                                                                          int                materialNumber_ )
  : materialProperties( matProperties_ ),
    nMaterialProperties( nMaterialProperties_ ),
    regularization( defaultRegularization() ),
    materialNumber( materialNumber_ )
{
  if ( nMaterialProperties < 3 )
    throw std::invalid_argument(
      "MarmotZInterfaceMaterialHypoElastic requires at least E, nu, and interface thickness h." );

  h = materialProperties[2];
  if ( h <= 0.0 )
    throw std::invalid_argument( "MarmotZInterfaceMaterialHypoElastic requires h > 0." );

  // Reference elastic acoustic tensor for the regularization, from the same
  // (E, nu) the whole interface family reads out of slots 0 and 1. nu is
  // clamped away from 1/2 so a near-incompressible layer cannot produce an
  // infinite reference stiffness; only the regularizer sees the clamp.
  const double youngsModulus = materialProperties[0];
  const double poissonRatio  = std::min( std::max( materialProperties[1], -0.999999 ), 0.499999 );
  shearModulus               = youngsModulus / ( 2.0 * ( 1.0 + poissonRatio ) );
  lameParameter              = youngsModulus * poissonRatio / ( ( 1.0 + poissonRatio ) * ( 1.0 - 2.0 * poissonRatio ) );

  baseMaterialProperties.reserve( nMaterialProperties - 1 );
  baseMaterialProperties.push_back( materialProperties[0] );
  baseMaterialProperties.push_back( materialProperties[1] );
  baseMaterialProperties.insert( baseMaterialProperties.end(),
                                 materialProperties + 3,
                                 materialProperties + nMaterialProperties );

  topMaterial = std::unique_ptr< MarmotMaterialHypoElastic >(
    MarmotLibrary::MarmotMaterialHypoElasticFactory::createMaterial( materialName,
                                                                     baseMaterialProperties.data(),
                                                                     static_cast< int >(
                                                                       baseMaterialProperties.size() ),
                                                                     materialNumber ) );
  bottomMaterial = std::unique_ptr< MarmotMaterialHypoElastic >(
    MarmotLibrary::MarmotMaterialHypoElasticFactory::createMaterial( materialName,
                                                                     baseMaterialProperties.data(),
                                                                     static_cast< int >(
                                                                       baseMaterialProperties.size() ),
                                                                     materialNumber ) );

  stateLayout.add( "topMaterialStateVars", topMaterial->getNumberOfRequiredStateVars() );
  stateLayout.add( "bottomMaterialStateVars", bottomMaterial->getNumberOfRequiredStateVars() );
  // Each face carries its own committed Cauchy stress: recovering them from a
  // single generalized force is exact only where t+ == t- holds identically,
  // and any residual imbalance would compound over increments.
  stateLayout.add( "stressPlus", 6 );
  stateLayout.add( "stressMinus", 6 );
  stateLayout.finalize();
}

void MarmotZInterfaceMaterialHypoElastic::setCharacteristicElementLength( double length )
{
  characteristicElementLength = length;
  if ( topMaterial )
    topMaterial->setCharacteristicElementLength( length );
  if ( bottomMaterial )
    bottomMaterial->setCharacteristicElementLength( length );
}

void MarmotZInterfaceMaterialHypoElastic::setGradientJumpRegularization( double zeta )
{
  if ( zeta < 0.0 )
    throw std::invalid_argument( "MarmotZInterfaceMaterialHypoElastic: the regularization zeta must be >= 0." );
  regularization = zeta;
}

void MarmotZInterfaceMaterialHypoElastic::computeStress( State&               state,
                                                         Tangents&            tangents,
                                                         const Deformation&   deformation,
                                                         const TimeIncrement& timeIncrement )
{
  if ( !topMaterial || !bottomMaterial )
    throw std::logic_error( "MarmotZInterfaceMaterialHypoElastic has no base material." );

  Vector3d     normal     = Eigen::Map< const Vector3d >( deformation.normal );
  const double normalNorm = normal.norm();
  if ( normalNorm <= 1.0e-12 )
    throw std::invalid_argument( "MarmotZInterfaceMaterialHypoElastic: interface normal is zero." );
  normal /= normalNorm;

  const Vector3d          separation    = Eigen::Map< const Vector3d >( deformation.separationVector );
  const InterfaceGeometry geometry      = evaluateInterfaceGeometry( normal, separation, h );
  const double            ell           = geometry.normalSeparation;
  const Vector3d&         dTangential   = geometry.tangentialSeparation;
  const double            sideThickness = 0.5 * h;

  const Eigen::Map< const Eigen::Matrix< double, 6, 1 > >  dU( deformation.dU );
  const Eigen::Map< const Eigen::Matrix< double, 18, 1 > > dSurfaceGradient( deformation.dSurfaceStrain );
  const Vector3d dz = Eigen::Map< const Vector3d >( deformation.dNormalGradientJump );

  // --- purely geometric operators: vec(G^pm) = Bx^pm x + Bz^pm z ---
  const Matrix9x3RowMajor Bw                = jumpToGradient( normal, ell );
  const Matrix9x3RowMajor Bn                = jumpToGradient( normal, 1.0 );
  const Matrix9dRowMajor  Asrf              = surfaceToGradient( normal, dTangential, ell );
  const Matrix9dRowMajor  I9                = Matrix9dRowMajor::Identity();
  const Matrix9dRowMajor  half_I_plus_Asrf  = 0.5 * ( I9 + Asrf );
  const Matrix9dRowMajor  half_Asrf_minus_I = 0.5 * ( Asrf - I9 );

  Matrix9x21RowMajor BxPlus         = Matrix9x21RowMajor::Zero();
  Matrix9x21RowMajor BxMinus        = Matrix9x21RowMajor::Zero();
  BxPlus.block< 9, 3 >( 0, offW )   = Bw;
  BxPlus.block< 9, 9 >( 0, offAp )  = half_I_plus_Asrf;
  BxPlus.block< 9, 9 >( 0, offAm )  = half_Asrf_minus_I;
  BxMinus.block< 9, 3 >( 0, offW )  = Bw;
  BxMinus.block< 9, 9 >( 0, offAp ) = half_Asrf_minus_I;
  BxMinus.block< 9, 9 >( 0, offAm ) = half_I_plus_Asrf;

  const Matrix9x3RowMajor BzPlus  = 0.5 * Bn;
  const Matrix9x3RowMajor BzMinus = -0.5 * Bn;

  // --- assemble the generalized strain increment x = (w, A+, A-) ---
  Vector21d x             = Vector21d::Zero();
  x.segment< 3 >( offW )  = dU.segment< 3 >( 0 ) - dU.segment< 3 >( 3 ); // w = [[u]] = u+ - u-
  x.segment< 9 >( offAp ) = dSurfaceGradient.segment< 9 >( 0 );
  x.segment< 9 >( offAm ) = dSurfaceGradient.segment< 9 >( 9 );

  const Vector9d dGPlus  = BxPlus * x + BzPlus * dz;
  const Vector9d dGMinus = BxMinus * x + BzMinus * dz;

  const Matrix3dRowMajor GPlus  = Eigen::Map< const Matrix3dRowMajor >( dGPlus.data() );
  const Matrix3dRowMajor GMinus = Eigen::Map< const Matrix3dRowMajor >( dGMinus.data() );

  const Vector6d strainIncrementPlusVoigt = ContinuumMechanics::VoigtNotation::strainToVoigt(
    Matrix3dRowMajor( 0.5 * ( GPlus + GPlus.transpose() ) ) );
  const Vector6d strainIncrementMinusVoigt = ContinuumMechanics::VoigtNotation::strainToVoigt(
    Matrix3dRowMajor( 0.5 * ( GMinus + GMinus.transpose() ) ) );

  // --- committed per-side state ---
  double* topStateVars     = stateLayout.getPtr( state.stateVars, "topMaterialStateVars" );
  double* bottomStateVars  = stateLayout.getPtr( state.stateVars, "bottomMaterialStateVars" );
  double* stressPlusStore  = stateLayout.getPtr( state.stateVars, "stressPlus" );
  double* stressMinusStore = stateLayout.getPtr( state.stateVars, "stressMinus" );

  const Matrix3dRowMajor stressCurrentPlus = ContinuumMechanics::VoigtNotation::voigtToStress(
    Vector6d( Eigen::Map< const Vector6d >( stressPlusStore ) ) );
  const Matrix3dRowMajor stressCurrentMinus = ContinuumMechanics::VoigtNotation::voigtToStress(
    Vector6d( Eigen::Map< const Vector6d >( stressMinusStore ) ) );

  const MarmotMaterialHypoElastic::timeInfo timeInfo{ timeIncrement.timeOld, timeIncrement.dT };

  // --- one constitutive evaluation per face. No local iteration: the
  //     traction imbalance is returned, not eliminated. ---
  const SideTrial plusTrial  = evaluateSideTrial( *topMaterial,
                                                 topStateVars,
                                                 topMaterial->getNumberOfRequiredStateVars(),
                                                 stressCurrentPlus,
                                                 strainIncrementPlusVoigt,
                                                 timeInfo );
  const SideTrial minusTrial = evaluateSideTrial( *bottomMaterial,
                                                  bottomStateVars,
                                                  bottomMaterial->getNumberOfRequiredStateVars(),
                                                  stressCurrentMinus,
                                                  strainIncrementMinusVoigt,
                                                  timeInfo );

  std::copy( plusTrial.trialStateVars.begin(), plusTrial.trialStateVars.end(), topStateVars );
  std::copy( minusTrial.trialStateVars.begin(), minusTrial.trialStateVars.end(), bottomStateVars );

  Eigen::Map< Vector6d > stressPlusStoreMap( stressPlusStore );
  Eigen::Map< Vector6d > stressMinusStoreMap( stressMinusStore );
  stressPlusStoreMap = ContinuumMechanics::VoigtNotation::stressToVoigt(
    Eigen::Matrix3d( 0.5 * ( plusTrial.stress + plusTrial.stress.transpose() ) ) );
  stressMinusStoreMap = ContinuumMechanics::VoigtNotation::stressToVoigt(
    Eigen::Matrix3d( 0.5 * ( minusTrial.stress + minusTrial.stress.transpose() ) ) );

  // --- regularization: an elastic spring of relative stiffness zeta on the
  //     z increment of this step, Q^e_ik = mu delta_ik + (lambda + mu) n_i n_k ---
  const Matrix3dRowMajor QElastic = shearModulus * Matrix3dRowMajor::Identity() +
                                    ( lameParameter + shearModulus ) * ( normal * normal.transpose() );
  const Matrix3dRowMajor QReg = ( 0.25 * h * regularization ) * QElastic;

  // --- generalized stresses p = d phi / d(x, z) ---
  const Eigen::Map< const Vector9d > vecSigmaPlus( plusTrial.stress.data() );
  const Eigen::Map< const Vector9d > vecSigmaMinus( minusTrial.stress.data() );

  const Vector21d pX = sideThickness * ( BxPlus.transpose() * vecSigmaPlus + BxMinus.transpose() * vecSigmaMinus );
  const Vector3d  pZ = sideThickness * ( BzPlus.transpose() * vecSigmaPlus + BzMinus.transpose() * vecSigmaMinus ) +
                      QReg * dz;

  Eigen::Map< Vector3d >( state.generalizedForce )   = pX.segment< 3 >( offW );
  Eigen::Map< Vector9d >( state.surfaceStressPlus )  = pX.segment< 9 >( offAp );
  Eigen::Map< Vector9d >( state.surfaceStressMinus ) = pX.segment< 9 >( offAm );
  Eigen::Map< Vector3d >( state.tractionImbalance )  = pZ;

  // --- Hessian blocks. Nothing is condensed and nothing is inverted. ---
  const Matrix9dRowMajor& CPlus  = plusTrial.CFull;
  const Matrix9dRowMajor& CMinus = minusTrial.CFull;

  Eigen::Map< Matrix21dRowMajor >( tangents.K_xx )  = sideThickness * ( BxPlus.transpose() * CPlus * BxPlus +
                                                                       BxMinus.transpose() * CMinus * BxMinus );
  Eigen::Map< Matrix21x3RowMajor >( tangents.K_xz ) = sideThickness * ( BxPlus.transpose() * CPlus * BzPlus +
                                                                        BxMinus.transpose() * CMinus * BzMinus );
  Eigen::Map< Matrix3x21RowMajor >( tangents.K_zx ) = sideThickness * ( BzPlus.transpose() * CPlus * BxPlus +
                                                                        BzMinus.transpose() * CMinus * BxMinus );
  // K_zz = (h/4) <Q> + (h/4) zeta Q^e, with <Q> = (Q+ + Q-)/2 the averaged
  // acoustic tensor. This is the block the equilibrated kernel had to invert.
  Eigen::Map< Matrix3dRowMajor >( tangents.K_zz ) = sideThickness * ( BzPlus.transpose() * CPlus * BzPlus +
                                                                      BzMinus.transpose() * CMinus * BzMinus ) +
                                                    QReg;
}

void MarmotZInterfaceMaterialHypoElastic::initializeYourself( double* stateVars, int nStateVars )
{
  if ( !topMaterial || !bottomMaterial ) {
    for ( int i = 0; i < nStateVars; ++i )
      stateVars[i] = 0.0;
    return;
  }

  for ( int i = 0; i < nStateVars; ++i )
    stateVars[i] = 0.0;

  topMaterial->initializeYourself( stateLayout.getPtr( stateVars, "topMaterialStateVars" ),
                                   topMaterial->getNumberOfRequiredStateVars() );
  bottomMaterial->initializeYourself( stateLayout.getPtr( stateVars, "bottomMaterialStateVars" ),
                                      bottomMaterial->getNumberOfRequiredStateVars() );
}

double MarmotZInterfaceMaterialHypoElastic::getDensity()
{
  if ( !topMaterial )
    return -1;

  return topMaterial->getDensity( nullptr );
}
