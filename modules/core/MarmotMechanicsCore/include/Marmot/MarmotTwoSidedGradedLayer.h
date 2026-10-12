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

#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotThicknessMicroproblem.h"
#include "Marmot/MarmotTwoSidedInterfaceGeometry.h"

#include <Eigen/Dense>

#include <algorithm>
#include <functional>
#include <vector>

namespace Marmot::Materials::TwoSidedGradedLayer {

  /// Maximum number of passes of the iteration for the derivative of the correction of the closure.
  constexpr int maxClosurePasses = 8;

  /// Largest change of the closure terms, relative to their norm, which ends the closure iteration.
  constexpr double closureTolerance = 1e-10;

  /// Generalized strain ( w, A+, A- ) of a two-sided interface layer.
  using Vector21d = Eigen::Matrix< double, 21, 1 >;

  /// Tangent over the generalized strain.
  using Matrix21d = Eigen::Matrix< double, 21, 21 >;

  /**
   * @brief Data of a layer with several integration points through the thickness and two surface gradients.
   *
   * @tparam nStrain number of components of the strain-like measure of a point: 6 for the strain in Voigt notation
   * (small strains), 9 for the deformation gradient (finite strains)
   */
  template < int nStrain >
  struct Layer {
    using StrainVector = Eigen::Matrix< double, nStrain, 1 >;

    Eigen::Vector3d                     normal;               ///< unit normal, from the bottom to the top face
    double                              normalSeparation;     ///< separation of the faces along the normal
    Eigen::Vector3d                     tangentialSeparation; ///< tangential part of the connector between the faces
    double                              thickness;         ///< thickness h of the layer, which scales the stress and u1
    std::vector< double >               weights;           ///< weights of the integration points, which sum up to one
    std::vector< double >               positions;         ///< coordinates y of the integration points
    const Eigen::MatrixXd*              cumulativeWeights; ///< weights of the cumulative integral from the bottom
    bool                                tangentialGradient; ///< the closure of the tangential gradient is kept
    Eigen::Matrix< double, nStrain, 9 > gradientToStrain;   ///< maps the 9-vector of a gradient to the strain-like
                                                            ///< measure
    Eigen::Matrix< double, 9, nStrain > stressToVector;     ///< maps the stress-like measure to the 9-vector of the
                                                            ///< stress tensor
    StrainVector strainOffset; ///< strain-like measure of the unloaded state: zero or the identity
  };

  /// Response of an integration point to a trial value of the correction.
  template < int nStrain >
  struct PointResponse {
    Eigen::Matrix< double, nStrain, 1 >       stress; ///< stress-like measure: Cauchy or first Piola-Kirchhoff
    Eigen::Matrix< double, nStrain, nStrain > D;      ///< derivative of the stress with respect to the strain
    std::vector< double >                     state;  ///< state of the point after the update
  };

  /// Response of an integration point, as the microproblem requires it.
  template < int nStrain >
  struct Point {
    Eigen::Vector3d                           w = Eigen::Vector3d::Zero(); ///< correction w = du1/dy
    Eigen::Matrix< double, nStrain, 1 >       stress;                      ///< stress-like measure
    Eigen::Matrix< double, nStrain, nStrain > D;                           ///< derivative of the stress
    Eigen::Matrix3d                           Qinv;                        ///< inverse of the acoustic tensor
    Eigen::Vector3d                           a;                           ///< traction plus the closure force
    std::vector< double >                     state;                       ///< state after the update
  };

  /// Result of the solution of the layer.
  template < int nStrain >
  struct Solution {
    std::vector< Point< nStrain > > points; ///< responses of the integration points at the solution
    Vector21d                       p;      ///< generalized stress ( force, S+, S- )
    Matrix21d                       K;      ///< tangent over the generalized strain ( w, A+, A- ), condensed
  };

  /**
   * @brief Solve the microproblem of a layer with two surface gradients and compute its generalized stress and tangent.
   *
   * With the layer coordinate y, the strain-like measure at the integration point q is
   * @f[
   *   e_q = e_0 + E \, B_q x + \hat{G}_q x + S w_q , \qquad
   *   B_q = \left[ \frac{1}{\ell} \, n \otimes \cdot \;\Big|\; y_q I + \tfrac{1}{2} B_\tau \;\Big|\;
   *   (1 - y_q) I + \tfrac{1}{2} B_\tau \right] ,
   * @f]
   * with the generalized strain @f$ x = ( [\![u]\!], A^+, A^- ) @f$, the operator @f$ S = E \, ( \cdot \otimes n ) @f$
   * of the correction, the tangential gradient @f$ \hat{G}_q x @f$ of the correction, and the traction of the point,
   * plus the closure force, is equal at all points. The generalized stress and the tangent are the weighted sums
   * over the points, with the tangent condensed over the corrections. See MarmotGradedXInterfaceMaterialHypoElastic
   * for the derivation of the closure of the tangential gradient.
   *
   * @tparam nStrain number of components of the strain-like measure of a point
   * @tparam EvaluateBulk callable `PointResponse< nStrain >( int q, const StrainVector& strain )` that updates the
   * bulk material of the point q for the strain-like measure from the committed state
   * @param[in] layer geometry and integration points of the layer
   * @param[in] x generalized strain: increments for small strains, totals for finite strains
   * @param[in] initialCorrections first estimate of the corrections of all points
   * @param[in] tolerance accepted component of the Newton step of the microproblem
   * @param[in] evaluateBulk update of the bulk material of a point
   * @return the response of the points, the generalized stress and the condensed tangent
   *
   * @throws StressUpdateFailed if the iteration does not converge or an acoustic tensor is singular
   */
  template < int nStrain, class EvaluateBulk >
  Solution< nStrain > solve( const Layer< nStrain >&               layer,
                             const Vector21d&                      x,
                             const std::vector< Eigen::Vector3d >& initialCorrections,
                             double                                tolerance,
                             EvaluateBulk&&                        evaluateBulk )
  {
    using namespace Marmot::Materials::TwoSidedInterface;

    using StrainVector = Eigen::Matrix< double, nStrain, 1 >;
    using Matrix3xN    = Eigen::Matrix< double, 3, nStrain >;
    using MatrixNxN    = Eigen::Matrix< double, nStrain, nStrain >;
    using MatrixNx21   = Eigen::Matrix< double, nStrain, 21 >;
    using Matrix3x21   = Eigen::Matrix< double, 3, 21 >;
    using Matrix9x21   = Eigen::Matrix< double, 9, 21 >;

    const int nq = static_cast< int >( layer.weights.size() );

    const Eigen::Vector3d& n   = layer.normal;
    const double           ell = layer.normalSeparation;
    const double           h   = layer.thickness;

    const Matrix9x9 I9   = Matrix9x9::Identity();
    const Matrix9x9 Hsrf = surfaceToGradient( n, layer.tangentialSeparation, ell ) - I9;
    const Matrix9x3 Bw   = jumpToGradient( n, ell );
    const Matrix9x3 N9   = jumpToGradient( n, 1.0 );

    // operator of the correction and of the traction, the transposed one of the stress is its 3 x nStrain image
    const Eigen::Matrix< double, nStrain, 3 > S  = layer.gradientToStrain * N9;
    Eigen::Matrix< double, 3, 9 >             Rt = Eigen::Matrix< double, 3, 9 >::Zero();
    for ( int i = 0; i < 3; i++ )
      for ( int b = 0; b < 3; b++ )
        Rt( i, 3 * i + b ) = n( b );
    const Matrix3xN St = Rt * layer.stressToVector;

    // generalized strain operator and strain of every point
    std::vector< Matrix9x21 >   Bx( nq );
    std::vector< MatrixNx21 >   EBx( nq );
    std::vector< StrainVector > strain( nq );
    for ( int q = 0; q < nq; q++ ) {
      const double y = layer.positions[q];
      Bx[q].setZero();
      Bx[q].template block< 9, 3 >( 0, 0 )  = Bw;
      Bx[q].template block< 9, 9 >( 0, 3 )  = y * I9 + 0.5 * Hsrf;
      Bx[q].template block< 9, 9 >( 0, 12 ) = ( 1. - y ) * I9 + 0.5 * Hsrf;
      EBx[q]                                = layer.gradientToStrain * Bx[q];
      strain[q]                             = layer.strainOffset + EBx[q] * x;
    }

    const Eigen::Matrix3d dA = Eigen::Map< const Matrix3dRowMajor >( x.data() + 3 ) -
                               Eigen::Map< const Matrix3dRowMajor >( x.data() + 12 );

    struct ClosureTerms {
      std::vector< StrainVector >    strain; ///< strain of the tangential gradient of the correction
      std::vector< Eigen::Vector3d > force;  ///< change of the traction from the bottom to the point
    };

    // Terms of the closure of the tangential gradient of the correction, linear in the difference of the surface
    // gradients, for the derivative of the correction and the tangent of the points.
    const auto closure = [&]( const std::vector< Matrix3xN >& dw,
                              const std::vector< MatrixNxN >& D,
                              const Eigen::Matrix3d&          surfaceGradientDifference ) {
      ClosureTerms terms;
      terms.strain.assign( nq, StrainVector::Zero() );
      terms.force.assign( nq, Eigen::Vector3d::Zero() );

      std::vector< Eigen::Matrix3d > derivativeOfCorrection( nq, Eigen::Matrix3d::Zero() ); // column r: d_r w_q
      std::vector< Eigen::Vector3d > source( nq, Eigen::Vector3d::Zero() );
      for ( int r = 0; r < 3; r++ ) {
        // tangential derivative of the strain-like measure along e_r
        const Eigen::Vector3d         column = surfaceGradientDifference.col( r );
        Eigen::Matrix< double, 9, 1 > gradient;
        for ( int k = 0; k < 3; k++ )
          for ( int l = 0; l < 3; l++ )
            gradient( 3 * k + l ) = column( k ) * n( l ) / ell;
        const StrainVector s = layer.gradientToStrain * gradient;

        for ( int q = 0; q < nq; q++ ) {
          const Eigen::Vector3d dwr          = dw[q] * s;
          derivativeOfCorrection[q].col( r ) = dwr;

          // derivative of the stress along e_r, of which the divergence drives the traction through the thickness
          const Eigen::Matrix< double, 9, 1 > dStress = layer.stressToVector * ( D[q] * ( s + S * dwr ) );
          for ( int k = 0; k < 3; k++ )
            source[q]( k ) += h * dStress( 3 * k + r );
        }
      }

      for ( int q = 0; q < nq; q++ ) {
        Eigen::Matrix3d gradient = Eigen::Matrix3d::Zero();
        for ( int p = 0; p < nq; p++ ) {
          terms.force[q] += ( *layer.cumulativeWeights )( q, p ) * source[p];
          gradient += ( *layer.cumulativeWeights )( q, p ) * h * derivativeOfCorrection[p];
        }
        // the tangential gradient of u1 as a gradient: entry (k, l) is the derivative of u1_k along e_l
        Eigen::Matrix< double, 9, 1 > flat;
        for ( int k = 0; k < 3; k++ )
          for ( int l = 0; l < 3; l++ )
            flat( 3 * k + l ) = gradient( k, l );
        terms.strain[q] = layer.gradientToStrain * flat;
      }
      return terms;
    };

    const auto makeEvaluate = [&]( const ClosureTerms& terms ) {
      return [&, terms]( int q, const Eigen::Vector3d& w ) {
        auto bulk = evaluateBulk( q, StrainVector( strain[q] + terms.strain[q] + S * w ) );

        Point< nStrain > point;
        point.w      = w;
        point.stress = bulk.stress;
        point.D      = bulk.D;
        point.state  = std::move( bulk.state );
        point.a      = St * point.stress + terms.force[q];

        point.Qinv = ( St * point.D * S ).inverse();
        if ( !point.Qinv.allFinite() )
          throw StressUpdateFailed( MakeString() << __PRETTY_FUNCTION__
                                                 << ": the acoustic tensor of an integration point through the "
                                                    "thickness is singular" );
        return point;
      };
    };

    const auto condenseOverStrain = [&]( const std::vector< Point< nStrain > >& points ) {
      std::vector< Matrix3xN > B( nq );
      for ( int q = 0; q < nq; q++ )
        B[q] = St * points[q].D;
      return ThicknessMicroproblem::condensedCorrectionTangents< nStrain >( layer.weights, points, B );
    };

    const bool closureActive = layer.tangentialGradient && dA.cwiseAbs().maxCoeff() > 0.;

    ClosureTerms             terms;
    std::vector< Matrix3xN > dwStrain( nq );
    std::vector< MatrixNxN > tangentPoint( nq );

    ClosureTerms none;
    none.strain.assign( nq, StrainVector::Zero() );
    none.force.assign( nq, Eigen::Vector3d::Zero() );

    if ( closureActive ) {
      // the derivative of the correction of the state at the start of the increment starts the iteration
      const auto                      evaluateNone = makeEvaluate( none );
      std::vector< Point< nStrain > > start;
      start.reserve( nq );
      for ( int q = 0; q < nq; q++ )
        start.push_back( evaluateNone( q, initialCorrections[q] ) );

      dwStrain = condenseOverStrain( start );
      for ( int q = 0; q < nq; q++ )
        tangentPoint[q] = start[q].D;
      terms = closure( dwStrain, tangentPoint, dA );
    }
    else {
      terms = none;
    }

    Solution< nStrain > solution;
    auto&               points = solution.points;

    for ( int pass = 0; pass < maxClosurePasses; pass++ ) {
      points   = ThicknessMicroproblem::solve< Point< nStrain > >( layer.weights,
                                                                 makeEvaluate( terms ),
                                                                 initialCorrections,
                                                                 tolerance );
      dwStrain = condenseOverStrain( points );
      for ( int q = 0; q < nq; q++ )
        tangentPoint[q] = points[q].D;

      if ( !closureActive )
        break;

      const ClosureTerms updated = closure( dwStrain, tangentPoint, dA );

      double change = 0.;
      double scale  = 0.;
      for ( int q = 0; q < nq; q++ ) {
        change = std::max( change, ( updated.strain[q] - terms.strain[q] ).cwiseAbs().maxCoeff() );
        change = std::max( change, ( updated.force[q] - terms.force[q] ).cwiseAbs().maxCoeff() );
        scale  = std::max( scale, updated.strain[q].cwiseAbs().maxCoeff() );
        scale  = std::max( scale, updated.force[q].cwiseAbs().maxCoeff() );
      }
      terms = updated;
      if ( change <= closureTolerance * scale + 1e-300 )
        break;
    }

    // linear dependence of the closure terms on the surface gradients A+ and A-, with the derivative of the
    // correction of the converged state frozen
    std::vector< MatrixNx21 > closureStrainTangent( nq, MatrixNx21::Zero() );
    std::vector< Matrix3x21 > closureForceTangent( nq, Matrix3x21::Zero() );
    for ( int k = 0; k < ( layer.tangentialGradient ? 9 : 0 ); k++ ) {
      Eigen::Matrix3d unit        = Eigen::Matrix3d::Zero();
      unit( k / 3, k % 3 )        = 1.0;
      const ClosureTerms response = closure( dwStrain, tangentPoint, unit );
      for ( int q = 0; q < nq; q++ ) {
        closureStrainTangent[q].col( 3 + k )  = response.strain[q];
        closureStrainTangent[q].col( 12 + k ) = -response.strain[q];
        closureForceTangent[q].col( 3 + k )   = response.force[q];
        closureForceTangent[q].col( 12 + k )  = -response.force[q];
      }
    }

    // condensation of the corrections over the generalized strain
    std::vector< Matrix3x21 > B( nq );
    for ( int q = 0; q < nq; q++ )
      B[q] = St * points[q].D * ( EBx[q] + closureStrainTangent[q] ) + closureForceTangent[q];
    const auto dw = ThicknessMicroproblem::condensedCorrectionTangents< 21 >( layer.weights, points, B );

    solution.p = Vector21d::Zero();
    solution.K = Matrix21d::Zero();
    for ( int q = 0; q < nq; q++ ) {
      const double                         scale   = h * layer.weights[q];
      const Eigen::Matrix< double, 9, 21 > dStress = layer.stressToVector * points[q].D *
                                                     ( EBx[q] + closureStrainTangent[q] + S * dw[q] );
      solution.p += scale * Bx[q].transpose() * ( layer.stressToVector * points[q].stress );
      solution.K += scale * Bx[q].transpose() * dStress;
    }
    return solution;
  }

} // namespace Marmot::Materials::TwoSidedGradedLayer
