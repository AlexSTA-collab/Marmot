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
 * Alexander Dummer alexander.dummer@uibk.ac.at
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

#include "Marmot/DufourModel.h"

#include "Marmot/MarmotDeformationMeasures.h"
#include "Marmot/MarmotEnergyDensityFunctions.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotFiniteStrainPlasticity.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotStressMeasures.h"
#include "Marmot/MarmotUtils.h"

#include <array>
#include <cmath>
#include <cstring>
#include <stdexcept>

namespace Marmot::Materials {

  using namespace Marmot;
  using namespace Fastor;
  using namespace FastorIndices;
  using namespace FastorStandardTensors;

  DufourModel::DufourModel( const double* materialProperties, int nMaterialProperties, int materialNumber )
    : MarmotMaterialGradientEnhancedFiniteStrain( materialProperties, nMaterialProperties, materialNumber ),
      K( checkedMaterialProperty( materialProperties, nMaterialProperties, 0 ) ),
      G( checkedMaterialProperty( materialProperties, nMaterialProperties, 1 ) ),
      ft( checkedMaterialProperty( materialProperties, nMaterialProperties, 2 ) ),
      fc( checkedMaterialProperty( materialProperties, nMaterialProperties, 3 ) ),
      Q1( checkedMaterialProperty( materialProperties, nMaterialProperties, 4 ) ),
      Q2( checkedMaterialProperty( materialProperties, nMaterialProperties, 5 ) ),
      b1( checkedMaterialProperty( materialProperties, nMaterialProperties, 6 ) ),
      b2( checkedMaterialProperty( materialProperties, nMaterialProperties, 7 ) ),
      b3( checkedMaterialProperty( materialProperties, nMaterialProperties, 8 ) ),
      b4( checkedMaterialProperty( materialProperties, nMaterialProperties, 9 ) ),
      b5( checkedMaterialProperty( materialProperties, nMaterialProperties, 10 ) ),
      eta_VP( checkedMaterialProperty( materialProperties, nMaterialProperties, 11 ) ),
      n( checkedMaterialProperty( materialProperties, nMaterialProperties, 12 ) ),
      nuP_plus( checkedMaterialProperty( materialProperties, nMaterialProperties, 13 ) ),
      nuP_minus( checkedMaterialProperty( materialProperties, nMaterialProperties, 14 ) ),
      epsF( checkedMaterialProperty( materialProperties, nMaterialProperties, 15 ) ),
      omegaMax( checkedMaterialProperty( materialProperties, nMaterialProperties, 16 ) ),
      ld( checkedMaterialProperty( materialProperties, nMaterialProperties, 17 ) ),
      m( checkedMaterialProperty( materialProperties, nMaterialProperties, 18 ) ),
      density( nMaterialProperties > 19 ? materialProperties[19] : 0.0 ),
      // The driver of the failure part is optional. Without cSW the damage is the original softening law. Without
      // bSW the compression term is switched off, so that only cSW has to be given to activate the driver.
      cSW( nMaterialProperties > 20 ? materialProperties[20] : 0.0 ),
      bSW( nMaterialProperties > 21 ? materialProperties[21] : 1e6 ),
      dF( nMaterialProperties > 22 ? materialProperties[22] : 1.0 ),
      swdfmB0( swdfmExtra( materialProperties, nMaterialProperties, 0 ) ),
      swdfmB1( swdfmExtra( materialProperties, nMaterialProperties, 1 ) ),
      swdfmB2( swdfmExtra( materialProperties, nMaterialProperties, 2 ) ),
      swdfmKdotRef( swdfmExtra( materialProperties, nMaterialProperties, 3 ) ),
      swdfmS( swdfmExtra( materialProperties, nMaterialProperties, 4 ) ),
      swdfmTCap( swdfmExtra( materialProperties, nMaterialProperties, 5 ) ),
      swdfmBumpA( swdfmExtra( materialProperties, nMaterialProperties, 6 ) ),
      swdfmBumpTc( swdfmExtra( materialProperties, nMaterialProperties, 7 ) ),
      swdfmBumpW( swdfmExtra( materialProperties, nMaterialProperties, 8 ) ),
      nMaxwell( nMaterialProperties > idxNMaxwell ? static_cast< int >( materialProperties[idxNMaxwell] ) : 0 ),
      maxwellDev( makeMaxwellProperties( materialProperties, nMaterialProperties, 0 ) ),
      maxwellVol( makeMaxwellProperties( materialProperties, nMaterialProperties, 1 ) )
  {
    if ( maxwellDev.sumGamma >= 1.0 || maxwellVol.sumGamma >= 1.0 )
      throw std::invalid_argument( "DufourModel: the sum of the Maxwell relative moduli must stay below 1, "
                                   "otherwise the equilibrium branch has non-positive stiffness" );

    // The card is positional. A card of an older layout has the entry dF at the index of the number of Maxwell
    // branches and more entries than the layout ends with, so that it would be read as a model without Maxwell
    // branches, silently. Such a card is recognized and refused.
    const double atNMaxwell = nMaterialProperties > idxNMaxwell ? materialProperties[idxNMaxwell] : 0.0;
    const bool   plausible  = atNMaxwell >= 0.0 && atNMaxwell <= static_cast< double >( nMaxwellMax ) &&
                           atNMaxwell == std::floor( atNMaxwell );
    if ( !plausible )
      throw std::invalid_argument( "DufourModel: entry [" + std::to_string( idxNMaxwell ) +
                                   "] must be the number of Maxwell branches, an integer in 0.." +
                                   std::to_string( nMaxwellMax ) + ", but it is " + std::to_string( atNMaxwell ) +
                                   ". This is the signature of a card of an older layout." );

    const int expected = idxPronyBase + 3 * static_cast< int >( atNMaxwell ) + nSwdfmExtra;
    if ( nMaterialProperties > expected )
      throw std::invalid_argument( "DufourModel: the card has " + std::to_string( nMaterialProperties ) +
                                   " entries but the layout ends at " + std::to_string( expected ) +
                                   " for this number of Maxwell branches. Extra entries mean a card of an older "
                                   "layout or a parameter that no longer exists." );

    // The shape of the exponent is mandatory while the driver is active. Without it E vanishes, so that g = 1 at every
    // positive triaxiality and the dependence on the stress state would be absent silently.
    if ( cSW != 0.0 && !swdfmShapeGiven() )
      throw std::invalid_argument( "DufourModel: the driver is active (cSW != 0) but the card gives no shape of the "
                                   "exponent (b0, b1, b2). Without it g = 1 at every positive triaxiality." );

    // A rate exponent without a reference rate is inert, because the rate factor would never be applied. It is more
    // likely a card written wrongly than a deliberate choice.
    if ( swdfmS != 0.0 && swdfmKdotRef <= 0.0 )
      throw std::invalid_argument( "DufourModel: a rate exponent of the driver is given but the reference rate is zero "
                                   "or absent, which switches the rate term off. Give the reference rate or clear "
                                   "the exponent." );

    stateLayout.add( "Fp", 9 );
    stateLayout.add( "alphaP", 1 );
    stateLayout.add( "omega", 1 );
    stateLayout.add( "damageDriver", 1 );
    stateLayout.add( "alphaPBar", 1 );
    stateLayout.add( "PK2Ref", 9 );
    stateLayout.add( "veDev", nMaxwellMax * 9 );
    stateLayout.add( "veVol", nMaxwellMax );
    stateLayout.add( "triax", 1 );
    stateLayout.add( "omegaS", 1 );
    stateLayout.add( "omegaF", 1 );
    stateLayout.add( "lode", 1 );
    stateLayout.finalize();
  }

  ContinuumMechanics::FiniteStrain::Viscoelasticity::MaxwellProperties DufourModel::makeMaxwellProperties(
    const double* materialProperties,
    int           nMaterialProperties,
    int           which )
  {
    const int n = nMaterialProperties > idxNMaxwell ? static_cast< int >( materialProperties[idxNMaxwell] ) : 0;
    if ( n <= 0 )
      return ContinuumMechanics::FiniteStrain::Viscoelasticity::createMaxwellProperties( 0, nullptr );
    if ( n > nMaxwellMax )
      throw std::invalid_argument( "DufourModel: too many Maxwell branches for the state variables" );
    if ( nMaterialProperties < idxPronyBase + 3 * n )
      throw std::invalid_argument( "DufourModel: incomplete (gammaG, gammaK, tau) triplets on the material card" );

    // The helper of Marmot expects interleaved pairs (gamma, tau), so the triplets are repacked.
    std::vector< double > pairs( 2 * n );
    for ( int i = 0; i < n; i++ ) {
      pairs[2 * i]     = materialProperties[idxPronyBase + 3 * i + which];
      pairs[2 * i + 1] = materialProperties[idxPronyBase + 2 + 3 * i];
    }
    return ContinuumMechanics::FiniteStrain::Viscoelasticity::createMaxwellProperties( n, pairs.data() );
  }

  double DufourModel::swdfmExtra( const double* materialProperties, int nMaterialProperties, int which )
  {
    const int n = nMaterialProperties > idxNMaxwell ? static_cast< int >( materialProperties[idxNMaxwell] ) : 0;
    const int i = idxPronyBase + 3 * ( n > 0 ? n : 0 ) + which;
    return nMaterialProperties > i ? materialProperties[i] : 0.0;
  }

  double DufourModel::swdfmLogG( const double T, const double kdot, double* kdotDLogGDKdot ) const
  {
    if ( kdotDLogGDKdot )
      *kdotDLogGDKdot = 0.0;

    double E;
    if ( T <= 0.0 ) {
      const double slope0 = swdfmB0 != 0.0 ? swdfmB0 * swdfmB0 : swdfmExponent;
      E                   = slope0 * T;
    }
    else {
      const double Tc = swdfmTCap > 0.0 ? std::min( T, swdfmTCap ) : T;

      E = ( ( ( ( swdfmB2 * swdfmB2 / 5.0 ) * Tc + swdfmB1 * swdfmB2 / 2.0 ) * Tc +
              ( swdfmB1 * swdfmB1 + 2.0 * swdfmB0 * swdfmB2 ) / 3.0 ) *
              Tc +
            swdfmB0 * swdfmB1 ) *
            Tc * Tc +
          swdfmB0 * swdfmB0 * Tc;
    }

    // a Gaussian bump on g, local in the triaxiality and equal to one outside of its band
    if ( swdfmBumpA != 0.0 && swdfmBumpW > 0.0 ) {
      const double z = ( T - swdfmBumpTc ) / swdfmBumpW;
      E += std::log1p( swdfmBumpA * std::exp( -z * z ) );
    }

    if ( swdfmKdotRef > 0.0 ) {
      const double s = ( swdfmS == swdfmSTieToN ) ? n : swdfmS;
      E += -s * std::log( std::max( kdot, swdfmKdotMin ) / swdfmKdotRef );
      if ( kdotDLogGDKdot && kdot > swdfmKdotMin )
        *kdotDLogGDKdot = -s;
    }
    return E;
  }

  double DufourModel::triaxiality( const Tensor33d& tau )
  {
    const double p = ( tau( 0, 0 ) + tau( 1, 1 ) + tau( 2, 2 ) ) / 3.0;

    double J2 = 0.0;
    for ( int i = 0; i < 3; i++ ) {
      for ( int j = 0; j < 3; j++ ) {
        const double s_ij = tau( i, j ) - ( i == j ? p : 0.0 );
        J2 += s_ij * s_ij;
      }
    }
    J2 *= 0.5;

    const double q = std::sqrt( std::max( 3.0 * J2, 0.0 ) );

    if ( q < 1e-12 * std::max( 1.0, std::abs( p ) ) )
      return 0.0;

    return p / q;
  }

  double DufourModel::lodeParameter( const Tensor33d& tau )
  {
    const double p = ( tau( 0, 0 ) + tau( 1, 1 ) + tau( 2, 2 ) ) / 3.0;

    double s[3][3];
    for ( int i = 0; i < 3; i++ )
      for ( int j = 0; j < 3; j++ )
        s[i][j] = tau( i, j ) - ( i == j ? p : 0.0 );

    double J2 = 0.0;
    for ( int i = 0; i < 3; i++ )
      for ( int j = 0; j < 3; j++ )
        J2 += s[i][j] * s[i][j];
    J2 *= 0.5;

    const double J3 = s[0][0] * ( s[1][1] * s[2][2] - s[1][2] * s[2][1] ) -
                      s[0][1] * ( s[1][0] * s[2][2] - s[1][2] * s[2][0] ) +
                      s[0][2] * ( s[1][0] * s[2][1] - s[1][1] * s[2][0] );

    // near-hydrostatic or unstressed: the deviatoric state is undefined, and reported as shear-like
    const double den = std::pow( std::max( J2, 0.0 ), 1.5 );
    if ( den < 1e-24 * std::max( 1.0, std::abs( p ) * std::abs( p ) * std::abs( p ) ) )
      return 0.0;

    // clamped because the round-off can push the parameter a few ULP past one in near-axisymmetric states
    return std::min( std::max( 1.5 * std::sqrt( 3.0 ) * J3 / den, -1.0 ), 1.0 );
  }

  std::tuple< double, double, double, double, double > DufourModel::computeOmega( const Context&   ctx,
                                                                                  const double     alphaP_local,
                                                                                  const double     alphaP_nonlocal,
                                                                                  const Tensor33d& tau_eff,
                                                                                  const double     D_old,
                                                                                  const double     alphaPBar_old ) const
  {
    const double alphaP_weighted = alphaP_nonlocal * m + alphaP_local * ( 1 - m );

    const double dAlphaP_weighted_dAlphaP_local    = 1 - m;
    const double dAlphaP_weighted_dAlphaP_nonlocal = m;

    if ( alphaP_weighted < 0.0 )
      return { 0.0, 0.0, 0.0, D_old, alphaPBar_old };

    const double alphaPBar_new = alphaP_weighted;

    // (1) The softening variable is always active, and reproduces the response before the initiation of failure.
    double omega_s, dOmega_s;
    softening( alphaP_weighted, omega_s, dOmega_s );

    // (2) The failure variable is zero until the driver reaches one, and steep afterwards.
    double omega_f  = 0.0;
    double dOmega_f = 0.0;
    double D        = D_old;

    if ( cSW != 0.0 ) {
      const double T = std::min( std::max( triaxiality( tau_eff ), etaMin ), etaMax );

      // The increment of the driver is needed before g, because g depends on the rate at which the driver accumulates.
      const double dAlphaPBar = std::max( alphaP_weighted - alphaPBar_old, 0.0 );
      const double kdot       = ctx.dT > 0.0 ? dAlphaPBar / ctx.dT : 0.0;

      // Rice-Tracey weight with the exponent E(T) and the rate factor; the compression term keeps its mirrored form.
      double       kdotDLogGDKdot = 0.0;
      const double logG           = swdfmLogG( T, kdot, &kdotDLogGDKdot );

      double g = std::exp( logG ) - std::exp( -logG ) / bSW;
      g        = std::max( g, 0.0 ); // the damage is irreversible under monotonic loading

      D = D_old + cSW * g * dAlphaPBar;

      // Ask for a smaller step rather than let the driver leap past the threshold.
      if ( D - D_old > dDMaxPerIncrement )
        throw Marmot::StressUpdateFailed( "DufourModel: damage increment " + std::to_string( D - D_old ) +
                                          " exceeds the limit " + std::to_string( dDMaxPerIncrement ) +
                                          "; a smaller time step is requested" );

      if ( D > 1.0 ) {
        omega_f = omegaFOfDriver( D );
        // dD/dAlphaPBar has the direct term and the rate sensitivity of g:
        //   dD/dk = cSW ( g + dk * dg/dkdot / dt ),  dk / ( kdot dt ) == 1,
        // so that kdot cancels exactly, and there is no division by a small rate.
        double dD_dAlphaPBar = cSW * g;
        if ( kdotDLogGDKdot != 0.0 )
          dD_dAlphaPBar += cSW * kdotDLogGDKdot * ( std::exp( logG ) + std::exp( -logG ) / bSW );
        dOmega_f = dOmegaFOfDriver( D ) * dD_dAlphaPBar;
      }
    }

    const double omega                   = 1.0 - ( 1.0 - omega_s ) * ( 1.0 - omega_f );
    const double dOmega_dAlphaP_weigthed = ( 1.0 - omega_f ) * dOmega_s + ( 1.0 - omega_s ) * dOmega_f;

    // Output variables, written from the very values that formed omega, so that the split cannot differ from the
    // omega that the stress update used.
    stateLayout.getAs< double& >( ctx.stateVars, "omegaS" ) = omega_s;
    stateLayout.getAs< double& >( ctx.stateVars, "omegaF" ) = omega_f;
    stateLayout.getAs< double& >( ctx.stateVars, "lode" )   = lodeParameter( tau_eff );

    // The triaxiality, the Lode parameter and the rate factor are held fixed in the tangent: the contributions of
    // their derivatives are omitted, which leaves the structure of the tangent unchanged. For the rate factor this is
    // mild, because it varies by a factor of about 1.7 over four decades of the rate.
    const double dOmega_dAlphaP_local    = dOmega_dAlphaP_weigthed * dAlphaP_weighted_dAlphaP_local;
    const double dOmega_dAlphaP_nonlocal = dOmega_dAlphaP_weigthed * dAlphaP_weighted_dAlphaP_nonlocal;

    if ( omega > omegaMax )
      return { omegaMax, 0.0, 0.0, D, alphaPBar_new };

    return { omega, dOmega_dAlphaP_local, dOmega_dAlphaP_nonlocal, D, alphaPBar_new };
  }

  std::tuple< Tensor33d, Tensor3333d > DufourModel::applyViscoelasticity( const Context&     ctx,
                                                                          const Tensor33d&   PK2_0,
                                                                          const Tensor3333d& dPK2_0_dCe,
                                                                          const bool         commit ) const
  {
    // without Maxwell branches the model is the hyperelastic-viscoplastic model exactly
    if ( nMaxwell == 0 )
      return { PK2_0, dPK2_0_dCe };

    using namespace ContinuumMechanics::FiniteStrain::Viscoelasticity;

    // Split the instantaneous stress and its tangent into a volumetric and a deviatoric part. The contraction below
    // gives d(tr PK2)/dCe because 2 d2Psi/dCdC is major symmetric.
    const double    p0       = trace( PK2_0 ) / 3.0;
    const Tensor33d PK2Vol_0 = p0 * Spatial3D::I;
    const Tensor33d PK2Dev_0 = PK2_0 - PK2Vol_0;

    const Tensor3333d dPK2Vol_dCe = ( 1.0 / 3.0 ) *
                                    Fastor::outer( Spatial3D::I,
                                                   Tensor33d( einsum< ijkl, kl >( dPK2_0_dCe, Spatial3D::I ) ) );
    const Tensor3333d dPK2Dev_dCe = dPK2_0_dCe - dPK2Vol_dCe;

    double* PK2Ref = stateLayout.getPtr( ctx.stateVars, "PK2Ref" );
    double* veDev  = stateLayout.getPtr( ctx.stateVars, "veDev" );
    double* veVol  = stateLayout.getPtr( ctx.stateVars, "veVol" );

    // The recurrence is driven by the increment since the last converged increment.
    const Tensor33d PK2RefOld( PK2Ref );
    const double    pRefOld = trace( PK2RefOld ) / 3.0;
    const Tensor33d dPK2Vol = ( p0 - pRefOld ) * Spatial3D::I;
    const Tensor33d dPK2Dev = PK2Dev_0 - ( PK2RefOld - pRefOld * Spatial3D::I );

    // The helper writes through its pointer, so that an evaluation without commit gets a scratch buffer.
    std::array< double, nMaxwellMax * 9 > devState{};
    std::array< double, nMaxwellMax * 9 > volState{};
    std::memcpy( devState.data(), veDev, nMaxwell * 9 * sizeof( double ) );
    for ( int i = 0; i < nMaxwell; i++ )
      for ( int d = 0; d < 3; d++ )
        volState[i * 9 + d * 4] = veVol[i]; // rebuild q_i * I

    Tensor33d   PK2Dev = PK2Dev_0, PK2Vol = PK2Vol_0;
    Tensor3333d tanDev = dPK2Dev_dCe, tanVol = dPK2Vol_dCe;

    evaluateGeneralizedMaxwellModel( PK2Dev, tanDev, dPK2Dev, ctx.dT, maxwellDev, devState.data() );
    evaluateGeneralizedMaxwellModel( PK2Vol, tanVol, dPK2Vol, ctx.dT, maxwellVol, volState.data() );

    if ( commit ) {
      std::memcpy( veDev, devState.data(), nMaxwell * 9 * sizeof( double ) );
      for ( int i = 0; i < nMaxwell; i++ )
        veVol[i] = volState[i * 9]; // spherical, so that the magnitude is kept only
      std::memcpy( PK2Ref, PK2_0.data(), 9 * sizeof( double ) );
    }

    return { Tensor33d( PK2Dev + PK2Vol ), Tensor3333d( tanDev + tanVol ) };
  }

  std::tuple< double, Tensor33d, double, Tensor33d, Tensor3333d, double, Tensor33d > DufourModel::
    yieldFunctionFromStress( const Tensor33d& mandelStress, const double betaP, const double J ) const
  {
    const double eta = fc / ft;
    Tensor33d    dev = deviatoric( mandelStress );
    const double J2  = 0.5 * Fastor::inner( dev, dev );
    const double I1  = trace( mandelStress );
    const double p   = I1 / 3.0;

    const double thetaPlus  = 9 * ( ( 1 - 2 * nuP_plus ) / ( 1 + nuP_plus ) ) / 2.0;
    const double thetaMinus = 9 * ( ( 1 - 2 * nuP_minus ) / ( 1 + nuP_minus ) ) / 2.0;

    const double A = eta - 1.0;
    const double B = sqrt( std::max( A * A * I1 * I1 + 12.0 * eta * J2, 1e-15 ) );

    const double f = ( A * I1 + B ) / ( 2.0 * eta ) / J - betaP;

    const double phiI1      = A * ( 1.0 + A * I1 / B ) / ( 2.0 * eta );
    const double phiJ2      = 3.0 / B;
    Tensor33d    df_dMandel = ( phiI1 * Spatial3D::I + phiJ2 * dev ) / J;

    double df_dBetaP = -1.0;

    // The flow rule is non-associative, so that the second derivative of the yield function is not needed.
    const double g = sqrt( std::max( 3 * J2 + thetaPlus * Math::macauly( p ) * Math::macauly( p ) +
                                       thetaMinus * Math::macauly( -p ) * Math::macauly( -p ),
                                     1e-15 ) );

    Tensor33d Q = 3.0 * dev +
                  2.0 * ( thetaPlus * Math::macauly( p ) - thetaMinus * Math::macauly( -p ) ) * Spatial3D::I / 3.0;

    Tensor3333d dQ_dMandel = 3.0 * ( Spatial3D::ISymm - Spatial3D::IHyd / 3.0 ) +
                             2.0 *
                               ( thetaPlus * Math::heavisideExclude0( p ) +
                                 thetaMinus * Math::heavisideExclude0( -p ) ) *
                               Spatial3D::IHyd / 9.0;

    Tensor33d   dg_dMandel          = Q / ( 2.0 * g );
    Tensor3333d d2g_dMandel_dMandel = dQ_dMandel / ( 2.0 * g ) - Fastor::outer( Q, Q ) / ( 4.0 * g * g * g );

    double    h          = sqrt( 2.0 * Fastor::inner( dg_dMandel, dg_dMandel ) / 3.0 );
    Tensor33d dh_dMandel = 2.0 * einsum< ijkl, kl, to_ij >( d2g_dMandel_dMandel, dg_dMandel ) / ( 3.0 * h );

    return { f, df_dMandel, df_dBetaP, dg_dMandel, d2g_dMandel_dMandel, h, dh_dMandel };
  }

  std::tuple< double, Tensor33d, double, Tensor33d, Tensor33d > DufourModel::yieldFunction( const Context&   ctx,
                                                                                            const Tensor33d& Fe,
                                                                                            const double     betaP,
                                                                                            const double     J ) const
  {
    Tensor33d   mandelStress;
    Tensor3333d dMandel_dFe;
    std::tie( mandelStress, dMandel_dFe ) = computeMandelStress( ctx, Fe );

    double      f, h, df_dBetaP;
    Tensor33d   df_dMandel, dh_dMandel, df_dFe, dg_dMandel, dh_dFe;
    Tensor3333d d2g_dMandel_dMandel;

    std::tie( f,
              df_dMandel,
              df_dBetaP,
              dg_dMandel,
              d2g_dMandel_dMandel,
              h,
              dh_dMandel ) = yieldFunctionFromStress( mandelStress, betaP, J );
    dh_dFe                 = einsum< mn, mnij, to_ij >( dh_dMandel, dMandel_dFe );
    df_dFe                 = einsum< mn, mnij, to_ij >( df_dMandel, dMandel_dFe );

    return { f, df_dFe, df_dBetaP, dg_dMandel, dh_dFe };
  }

  std::tuple< double > DufourModel::yieldFunctionNominal( const Context&   ctx,
                                                          const Tensor33d& Fe,
                                                          const double     betaP,
                                                          const double     omega,
                                                          const double     J ) const
  {
    Tensor33d   mandelStressN;
    Tensor3333d dMandel_dFe;
    std::tie( mandelStressN, dMandel_dFe ) = computeNominalMandelStress( ctx, Fe, omega );

    const double eta = fc / ft;
    Tensor33d    dev = deviatoric( mandelStressN );
    const double J2  = 0.5 * Fastor::inner( dev, dev );
    const double I1  = trace( mandelStressN );

    const double A = eta - 1.0;
    const double B = sqrt( std::max( A * A * I1 * I1 + 12.0 * eta * J2, 1e-15 ) );

    const double f = ( A * I1 + B ) / ( 2.0 * eta ) / J - betaP;

    return { f };
  }

  bool DufourModel::isYielding( const Context& ctx, const Tensor33d& Fe, const double betaP, const double J ) const
  {
    return std::get< 0 >( yieldFunction( ctx, Fe, betaP, J ) ) > 0.0;
  }

  std::tuple< Tensor33d, Tensor3333d > DufourModel::computeMandelStress( const Context& ctx, const Tensor33d& Fe ) const
  {
    using namespace Marmot::ContinuumMechanics;
    Tensor33d   Ce;
    Tensor3333d dCe_dFe;
    std::tie( Ce, dCe_dFe ) = DeformationMeasures::FirstOrderDerived::rightCauchyGreen( Fe );

    double      psi_;
    Tensor33d   dPsi_dCe;
    Tensor3333d d2Psi_dCedCe, dMandel_dCe;

    std::tie( psi_, dPsi_dCe, d2Psi_dCedCe ) = EnergyDensityFunctions::SecondOrderDerived::PenceGouPotentialB( Ce,
                                                                                                               K,
                                                                                                               G );

    // The viscoelastic relaxation happens here, upstream of the return map, so that the yield function sees the
    // relaxed stress. It is called inside the Newton iteration, so that the history is not committed.
    auto [PK2, dPK2_dCe] = applyViscoelasticity( ctx,
                                                 Tensor33d( 2.0 * dPsi_dCe ),
                                                 Tensor3333d( 2.0 * d2Psi_dCedCe ),
                                                 false );

    const Tensor33d mandel = Ce % PK2;
    dMandel_dCe            = einsum< Ii, iJKL, to_IJKL >( Ce, dPK2_dCe ) +
                  einsum< IK, iL, iJ, to_IJKL >( Spatial3D::I, Spatial3D::I, PK2 );
    Tensor3333d dMandel_dFe = einsum< IJKL, KLMN >( dMandel_dCe, dCe_dFe );
    return { mandel, dMandel_dFe };
  }

  std::tuple< Tensor33d, Tensor3333d > DufourModel::computeNominalMandelStress( const Context&   ctx,
                                                                                const Tensor33d& Fe,
                                                                                const double&    omega ) const
  {
    using namespace Marmot::ContinuumMechanics;
    Tensor33d   Ce;
    Tensor3333d dCe_dFe;
    std::tie( Ce, dCe_dFe ) = DeformationMeasures::FirstOrderDerived::rightCauchyGreen( Fe );

    double      psi_;
    Tensor33d   dPsi_dCe;
    Tensor3333d d2Psi_dCedCe, dMandelN_dCe;

    std::tie( psi_, dPsi_dCe, d2Psi_dCedCe ) = EnergyDensityFunctions::SecondOrderDerived::PenceGouPotentialB( Ce,
                                                                                                               K,
                                                                                                               G );

    auto [PK2, dPK2_dCe] = applyViscoelasticity( ctx,
                                                 Tensor33d( 2.0 * dPsi_dCe ),
                                                 Tensor3333d( 2.0 * d2Psi_dCedCe ),
                                                 false );

    const Tensor33d mandelN = ( Ce % PK2 ) * ( 1.0 - omega );

    dMandelN_dCe = ( einsum< Ii, iJKL, to_IJKL >( Ce, dPK2_dCe ) +
                     einsum< IK, iL, iJ, to_IJKL >( Spatial3D::I, Spatial3D::I, PK2 ) ) *
                   ( 1.0 - omega );

    Tensor3333d dMandelN_dFe = einsum< IJKL, KLMN >( dMandelN_dCe, dCe_dFe );

    return { mandelN, dMandelN_dFe };
  }

  std::tuple< double, double > DufourModel::computeBetaP( const double alphaP ) const
  {
    const double beta = ft + Q1 * alphaP * exp( -b1 * alphaP ) + Q2 * ( 1.0 - exp( -b2 * alphaP ) ) +
                        b3 * alphaP * alphaP * alphaP + b4 * alphaP * alphaP + b5 * alphaP;
    const double dBetaP_dAlphaP = Q1 * exp( -b1 * alphaP ) * ( 1.0 - b1 * alphaP ) + Q2 * b2 * exp( -b2 * alphaP ) +
                                  3.0 * b3 * alphaP * alphaP + 2.0 * b4 * alphaP + b5;
    return { beta, dBetaP_dAlphaP };
  }

  std::tuple< Eigen::VectorXd, Eigen::MatrixXd > DufourModel::computeResidualVectorAndTangent( const Context& ctx,
                                                                                               const Eigen::VectorXd& X,
                                                                                               const Tensor33d& FeTrial,
                                                                                               const double alphaPTrial,
                                                                                               const double dt,
                                                                                               const double J ) const
  {
    const int idxA = 9;
    const int idxF = 10;
    using namespace Eigen;
    using mV9d = Eigen::Map< Eigen::Matrix< double, 9, 1 > >;
    using mM9d = Eigen::Map< Eigen::Matrix< double, 9, 9 > >;
    VectorXd R( 11 );
    MatrixXd dR_dX( 11, 11 );
    dR_dX.setZero();

    // initialize the residual
    R.segment< 9 >( 0 ) = -mV9d( const_cast< double* >( FeTrial.data() ) );
    R( 9 )              = -alphaPTrial;

    Tensor33d Fe( X.segment( 0, 9 ).data() );

    double       h;
    const double dLambda = X( 10 );
    const double alphaP  = X( 9 );

    double betaP, dBetaP_dAlphaP;
    std::tie( betaP, dBetaP_dAlphaP ) = computeBetaP( alphaP );

    // Mandel stress
    Tensor33d   mandelStress;
    Tensor3333d dMandel_dFe, d2g_dMandel_dMandel;
    std::tie( mandelStress, dMandel_dFe ) = computeMandelStress( ctx, Fe );

    double    f, df_dBetaP;
    Tensor33d df_dMandel, dg_dMandel, dh_dMandel, df_dFe, dh_dFe;
    std::tie( f,
              df_dMandel,
              df_dBetaP,
              dg_dMandel,
              d2g_dMandel_dMandel,
              h,
              dh_dMandel ) = yieldFunctionFromStress( mandelStress, betaP, J );

    Tensor33d   dGp = dLambda * dg_dMandel;
    Tensor33d   dFp;
    Tensor3333d ddFp_ddGp;
    std::tie( dFp, ddFp_ddGp ) = ContinuumMechanics::FiniteStrain::Plasticity::FlowIntegration::FirstOrderDerived::
      exponentialMap( dGp );

    Tensor3333d ddGp_dFe      = dLambda * einsum< ijmn, mnkL >( d2g_dMandel_dMandel, dMandel_dFe );
    Tensor33d   ddFp_ddLambda = einsum< IJKL, KL >( ddFp_ddGp, dg_dMandel );

    Tensor3333d ddFp_dFe = einsum< iImn, mnkL >( ddFp_ddGp, ddGp_dFe );

    Tensor3333d dFeTrial_dFe = einsum< iI, IJKL >( Fe, ddFp_dFe ) +
                               einsum< IK, JL, to_IJKL >( Spatial3D::I, transpose( Spatial3D::I % dFp ) );

    Tensor33d dFe_ddLambda = einsum< Ii, iJ >( Fe, ddFp_ddLambda );

    std::tie( f, df_dFe, df_dBetaP, dg_dMandel, dh_dFe ) = yieldFunction( ctx, Fe, betaP, J );

    double beta_min           = 1e-12;
    double sgn_beta           = ( betaP >= 0 ) ? 1.0 : -1.0;
    double betaP_cap          = sgn_beta * std::max( std::abs( betaP ), beta_min );
    double dBetaP_dAlphaP_cap = ( std::abs( betaP ) > beta_min ) ? sgn_beta * dBetaP_dAlphaP : 0.0;

    double ratio = Math::macauly( f + betaP_cap ) / betaP_cap;

    const double    hmin  = 1e-8;
    const double    hsafe = std::max( h, hmin );
    const Tensor33d hZero( 0.0 );
    const Tensor33d dh_dMandel_safe = ( h > hmin ) ? dh_dMandel : hZero;
    const Tensor33d dh_dFe_safe     = einsum< mn, mnij, to_ij >( dh_dMandel_safe, dMandel_dFe );

    double q  = dLambda * hsafe / ( dt * eta_VP );
    double qn = std::pow( std::max( q, 1e-30 ), n );

    // residual
    R.segment< 9 >( 0 ) += mV9d( Tensor33d( einsum< iJ, JK >( Fe, dFp ) ).data() );
    R( idxA ) += ( alphaP - dLambda * hsafe );
    R( idxF ) = ratio - qn;

    // Jacobian
    // dR_dFe
    dR_dX.block< 9, 9 >( 0, 0 )    = mM9d( dFeTrial_dFe.data() ).transpose();
    dR_dX.block< 9, 1 >( 0, idxF ) = mV9d( dFe_ddLambda.data() );

    // dR_dalphaP
    Tensor33d dRa_dFe              = -dLambda * einsum< mn, mnij, to_ij >( dh_dMandel_safe, dMandel_dFe );
    dR_dX.block< 1, 9 >( idxA, 0 ) = mV9d( dRa_dFe.data() ).transpose();
    dR_dX( idxA, idxA )            = 1.0;
    dR_dX( idxA, idxF )            = -hsafe;

    // dR_dLambda
    Tensor33d dRl_dFe = df_dFe / betaP_cap - ( n * qn / hsafe ) * dh_dFe_safe;

    dR_dX.block< 1, 9 >( idxF, 0 ) = mV9d( dRl_dFe.data() ).transpose();
    dR_dX( idxF, idxA )            = -dBetaP_dAlphaP_cap * ratio / betaP_cap;
    dR_dX( idxF, idxF )            = -n * qn / std::max( dLambda, 1e-30 );

    return { R, dR_dX };
  }

  void DufourModel::computeStress( ConstitutiveResponse< 3 >& response,
                                   AlgorithmicModuli< 3 >&    tangents,
                                   const Deformation< 3 >&    deformation,
                                   const TimeIncrement&       timeIncrement ) const
  {
    using namespace Eigen;

    const Context ctx{ response.stateVars, timeIncrement.dT };

    double* FpPtr = stateLayout.getPtr( ctx.stateVars, "Fp" );

    const Tensor33d FpOld( FpPtr );
    double&         alphaP       = stateLayout.getAs< double& >( ctx.stateVars, "alphaP" );
    const double    alphaPOld    = alphaP;
    double&         omega        = stateLayout.getAs< double& >( ctx.stateVars, "omega" );
    double&         driver       = stateLayout.getAs< double& >( ctx.stateVars, "damageDriver" );
    const double    driverOld    = driver;
    double&         alphaPBar    = stateLayout.getAs< double& >( ctx.stateVars, "alphaPBar" );
    const double    alphaPBarOld = alphaPBar;
    double&         triax        = stateLayout.getAs< double& >( ctx.stateVars, "triax" );

    // The interaction radius is the value of the card, constant in space and time.
    response.nonLocalRadius      = ld;
    const double alphaP_nonlocal = deformation.N;

    Tensor33d FeTrial = deformation.F % Fastor::inverse( FpOld );
    double    betaP, dBetaP_dAlphaP;
    std::tie( betaP, dBetaP_dAlphaP ) = computeBetaP( alphaPOld );

    Tensor33d dFp;
    dFp.eye();
    Tensor33d Fe = FeTrial;

    const double J = determinant( deformation.F );

    if ( isYielding( ctx, FeTrial, betaP, J ) ) {

      size_t counter = 0;

      using mV9d = Eigen::Map< Eigen::Matrix< double, 9, 1 > >;
      VectorXd X( 11 );
      X.segment( 0, 9 ) = mV9d( FeTrial.data() );
      X( 9 )            = alphaPOld;
      Tensor33d   mandelTrial;
      Tensor3333d dMandelTrial_dFe;
      std::tie( mandelTrial, dMandelTrial_dFe ) = computeMandelStress( ctx, FeTrial );
      const double hTrial                       = std::get< 5 >( yieldFunctionFromStress( mandelTrial, betaP, J ) );
      X( 10 )                                   = timeIncrement.dT * eta_VP / std::max( hTrial, 1e-8 );
      VectorXd        dX                        = VectorXd::Zero( 11 );
      VectorXd        R                         = VectorXd::Zero( 11 );
      Eigen::MatrixXd dR_dX( 11, 11 );

      std::tie( R, dR_dX ) = computeResidualVectorAndTangent( ctx, X, FeTrial, alphaPOld, timeIncrement.dT, J );

      while ( R.norm() > 1e-12 || dX.norm() > 1e-12 ) {

        if ( counter > 15 )
          throw std::runtime_error( "inner newton not converged" );

        dX = -dR_dX.colPivHouseholderQr().solve( R );
        X += dX;
        std::tie( R, dR_dX ) = computeResidualVectorAndTangent( ctx, X, FeTrial, alphaPOld, timeIncrement.dT, J );
        counter += 1;
      }
      // update the plastic deformation
      Fe              = X.segment( 0, 9 ).data();
      dFp             = Fastor::inverse( Fe ) % FeTrial;
      alphaP          = X( 9 );
      response.L      = alphaP;
      Tensor33d FpNew = dFp % FpOld;
      std::memcpy( FpPtr, FpNew.data(), 9 * sizeof( double ) );

      using namespace ContinuumMechanics;
      double      psi_, dOmega_dAlphaP_local, dOmega_dAlphaP_nonlocal;
      Tensor33d   Ce, dPsi_dCe, tau_eff;
      Tensor3333d dCe_dFe, d2Psi_dCedCe, dTau_dPK2_eff, dTau_dFe_partial_eff;
      std::tie( Ce, dCe_dFe ) = DeformationMeasures::FirstOrderDerived::rightCauchyGreen( Fe );

      // energy density and its first and second derivatives with respect to the right Cauchy-Green tensor
      std::tie( psi_, dPsi_dCe, d2Psi_dCedCe ) = EnergyDensityFunctions::SecondOrderDerived::PenceGouPotentialB( Ce,
                                                                                                                 K,
                                                                                                                 G );
      // The Kirchhoff stress comes first: the damage evolution depends on the triaxiality, so that the stress state has
      // to be known before omega can be evaluated. The viscoelastic relaxation is applied here as well, without commit.
      auto [PK2_eff, dPK2_dCe_eff] = applyViscoelasticity( ctx,
                                                           Tensor33d( 2.0 * dPsi_dCe ),
                                                           Tensor3333d( 2.0 * d2Psi_dCedCe ),
                                                           false );

      std::tie( tau_eff,
                dTau_dPK2_eff,
                dTau_dFe_partial_eff ) = StressMeasures::FirstOrderDerived::KirchhoffStressFromPK2( PK2_eff, Fe );

      // damage variable
      std::tie( omega,
                dOmega_dAlphaP_local,
                dOmega_dAlphaP_nonlocal,
                driver,
                alphaPBar ) = computeOmega( ctx, alphaP, alphaP_nonlocal, tau_eff, driverOld, alphaPBarOld );

      // the triaxiality that the driver saw at this point, clamped as in computeOmega
      triax = std::min( std::max( triaxiality( tau_eff ), etaMin ), etaMax );

      response.tau                  = tau_eff * ( 1.0 - omega );
      response.elasticEnergyDensity = psi_;

      // tangent operator
      using mM9d = Eigen::Map< Eigen::Matrix< double, 9, 9 > >;

      MatrixXd dYdDeformation              = MatrixXd::Zero( 11, 11 );
      dYdDeformation.block< 9, 9 >( 0, 0 ) = mM9d( Tensor3333d( einsum< IK, JL, to_IJKL >( Spatial3D::I,
                                                                                           transpose( Fastor::inverse(
                                                                                             FpOld ) ) ) )
                                                     .data() )
                                               .transpose();

      const int       idxF   = 10;
      const Tensor33d Finv_T = transpose( Fastor::inverse( deformation.F ) );

      double betaP_c, dbeta_c;
      std::tie( betaP_c, dbeta_c ) = computeBetaP( alphaP );
      Tensor33d   mandel_c;
      Tensor3333d dmandel_c;
      std::tie( mandel_c, dmandel_c ) = computeMandelStress( ctx, Fe );
      const double f_c                = std::get< 0 >( yieldFunctionFromStress( mandel_c, betaP_c, J ) );
      const double betaP_cap_c        = std::max( betaP_c, 1e-12 );
      const double ratio_conv         = Math::macauly( f_c + betaP_cap_c ) / betaP_cap_c;

      dYdDeformation.block< 1, 9 >( idxF, 0 ) = ratio_conv * mV9d( Finv_T.data() ).transpose();

      MatrixXd dXdDeformation = dR_dX.colPivHouseholderQr().solve( dYdDeformation );

      Tensor3333d dFe_dF = Tensor3333d( Matrix9d( dXdDeformation.block< 9, 9 >( 0, 0 ).transpose() ).data() );

      Tensor3333d dPK2_dFe_eff = einsum< ijKL, KLMN >( dPK2_dCe_eff, dCe_dFe );
      Tensor3333d dPK2_dF_eff  = einsum< ijKL, KLMN >( dPK2_dFe_eff, dFe_dF );

      Tensor3333d dTau_dF_eff = einsum< IJKL, KLMN >( dTau_dPK2_eff, dPK2_dF_eff ) +
                                einsum< ijKL, KLMN >( dTau_dFe_partial_eff, dFe_dF );

      Tensor33d dAlphaP_local_dF = Tensor33d( Vector9d( dXdDeformation.block< 1, 9 >( 9, 0 ).transpose() ).data() );

      // the local driving strain is the plastic strain
      const Tensor33d& dLocal_dF = dAlphaP_local_dF;

      tangents.dTau_dF = ( 1 - omega ) * dTau_dF_eff - dOmega_dAlphaP_local * Fastor::outer( tau_eff, dLocal_dF );
      tangents.dTau_dN = -tau_eff * dOmega_dAlphaP_nonlocal;
      tangents.dL_dF   = dLocal_dF;
    }
    else {
      using namespace Marmot::ContinuumMechanics;
      double      psi_, dOmega_dAlphaP_local, dOmega_dAlphaP_nonlocal;
      Tensor33d   Ce, dPsi_dCe, tau_eff;
      Tensor3333d dCe_dFe, d2Psi_dCedCe, dTau_dPK2_eff, dTau_dFe_partial_eff;
      std::tie( Ce, dCe_dFe ) = DeformationMeasures::FirstOrderDerived::rightCauchyGreen( Fe );

      std::tie( psi_, dPsi_dCe, d2Psi_dCedCe ) = EnergyDensityFunctions::SecondOrderDerived::PenceGouPotentialB( Ce,
                                                                                                                 K,
                                                                                                                 G );
      // the Kirchhoff stress comes first, the triaxiality enters the damage evolution
      auto [PK2_eff, dPK2_dCe_eff] = applyViscoelasticity( ctx,
                                                           Tensor33d( 2.0 * dPsi_dCe ),
                                                           Tensor3333d( 2.0 * d2Psi_dCedCe ),
                                                           false );

      std::tie( tau_eff,
                dTau_dPK2_eff,
                dTau_dFe_partial_eff ) = StressMeasures::FirstOrderDerived::KirchhoffStressFromPK2( PK2_eff, Fe );

      std::tie( omega,
                dOmega_dAlphaP_local,
                dOmega_dAlphaP_nonlocal,
                driver,
                alphaPBar ) = computeOmega( ctx, alphaPOld, alphaP_nonlocal, tau_eff, driverOld, alphaPBarOld );

      triax = std::min( std::max( triaxiality( tau_eff ), etaMin ), etaMax );

      response.tau                  = tau_eff * ( 1.0 - omega );
      response.elasticEnergyDensity = psi_;
      response.L                    = alphaP;

      // tangent operator
      Tensor3333d dPK2_dFe    = einsum< ijKL, KLMN >( dPK2_dCe_eff, dCe_dFe );
      Tensor3333d dFe_dF      = einsum< IK, JL, to_IJKL >( Spatial3D::I, transpose( Fastor::inverse( FpOld ) ) );
      Tensor3333d dPK2_dF_eff = einsum< ijKL, KLMN >( dPK2_dFe, dFe_dF );

      Tensor3333d dTau_dF_eff = einsum< IJKL, KLMN >( dTau_dPK2_eff, dPK2_dF_eff ) +
                                einsum< ijKL, KLMN >( dTau_dFe_partial_eff, dFe_dF );

      Tensor33d dAlphaP_local_dF = Tensor33d( 0.0 );

      tangents.dTau_dF = ( 1 - omega ) * dTau_dF_eff;
      tangents.dTau_dN = -tau_eff * dOmega_dAlphaP_nonlocal;
      tangents.dL_dF   = dAlphaP_local_dF;
    }

    tangents.dL_dN = 0.0;

    // Advance the Maxwell branches exactly once per increment, with the converged elastic deformation. Everything above
    // ran without commit, so that the return map and the evaluations after its convergence could evaluate the relaxed
    // stress as often as needed without corrupting the history.
    if ( nMaxwell > 0 ) {
      using namespace Marmot::ContinuumMechanics;
      Tensor33d   CeFinal;
      Tensor3333d dCeFinal_dFe;
      std::tie( CeFinal, dCeFinal_dFe ) = DeformationMeasures::FirstOrderDerived::rightCauchyGreen( Fe );

      double      psiFinal;
      Tensor33d   dPsiFinal_dCe;
      Tensor3333d d2PsiFinal_dCedCe;
      std::tie( psiFinal,
                dPsiFinal_dCe,
                d2PsiFinal_dCedCe ) = EnergyDensityFunctions::SecondOrderDerived::PenceGouPotentialB( CeFinal, K, G );

      applyViscoelasticity( ctx, Tensor33d( 2.0 * dPsiFinal_dCe ), Tensor3333d( 2.0 * d2PsiFinal_dCedCe ), true );
    }
  }

  void DufourModel::initializeYourself( double* stateVars, int nStateVars )
  {
    for ( int i = 0; i < nStateVars; ++i )
      stateVars[i] = 0.0;

    // all other variables start at zero: the plastic strain, the damage, the driver, the weighted plastic strain, the
    // reference stress, the branches of the Maxwell chain and the output variables
    double* Fp = stateLayout.getPtr( stateVars, "Fp" );
    for ( int i = 0; i < 3; ++i )
      Fp[3 * i + i] = 1.0;
  }

  double DufourModel::getDensity( const double* ) const
  {
    return density;
  }

} // namespace Marmot::Materials
