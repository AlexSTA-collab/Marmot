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

    /// derivative of the stress with respect to the nonlocal field of the point, if the material has one
    Eigen::Matrix< double, nStrain, 1 > dStress_dNonlocal = Eigen::Matrix< double, nStrain, 1 >::Zero();
    /// derivative of the local driving field with respect to the strain
    Eigen::Matrix< double, 1, nStrain > dLocal_dStrain   = Eigen::Matrix< double, 1, nStrain >::Zero();
    double                              dLocal_dNonlocal = 0.; ///< derivative of the local field with respect to the
                                                               ///< nonlocal field of the point
    double local  = 0.;                                        ///< local driving field
    double radius = 0.;                                        ///< nonlocal radius
  };

  /// Response of an integration point, as the microproblem requires it.
  template < int nStrain >
  struct Point : PointResponse< nStrain > {
    Eigen::Vector3d w = Eigen::Vector3d::Zero(); ///< correction w = du1/dy
    Eigen::Matrix3d Qinv;                        ///< inverse of the acoustic tensor
    Eigen::Vector3d a;                           ///< traction plus the closure force
  };

  /// Result of the solution of the layer.
  template < int nStrain >
  struct Solution {
    std::vector< Point< nStrain > > points; ///< responses of the integration points at the solution
    Vector21d                       p;      ///< generalized stress ( force, S+, S- )
    Matrix21d                       K;      ///< tangent over the generalized strain ( w, A+, A- ), condensed
  };

  /**
   * @brief Solver for the microproblem of a layer with two surface gradients.
   *
   * With the layer coordinate y, the strain-like measure at the integration point q is
   * @f[
   *   e_q = e_0 + E \, B_q x + \hat{G}_q x + S w_q , \qquad
   *   B_q = \left[ \frac{1}{\ell} \, n \otimes \cdot \;\Big|\; y_q I + \tfrac{1}{2} B_\tau \;\Big|\;
   *   (1 - y_q) I + \tfrac{1}{2} B_\tau \right] ,
   * @f]
   * with the generalized strain @f$ x = ( [\![u]\!], A^+, A^- ) @f$, the operator @f$ S = E \, ( \cdot \otimes n ) @f$
   * of the correction, and the tangential gradient @f$ \hat{G}_q x @f$ of the correction. The traction of the point,
   * plus the closure force, is equal at all points. The generalized stress and the tangent are the weighted sums
   * over the points, with the tangent condensed over the corrections. See MarmotGradedXInterfaceMaterialHypoElastic
   * for the derivation of the closure of the tangential gradient.
   *
   * A material with further variables that enter the stress of the points, for example nonlocal fields, adds their
   * derivative of the traction to the condensation and assembles the derivatives of the generalized stress with
   * respect to them from the condensed corrections and the accessors of the points.
   *
   * @tparam nStrain number of components of the strain-like measure of a point
   */
  template < int nStrain >
  class LayerSolver {

  public:
    using StrainVector = Eigen::Matrix< double, nStrain, 1 >;
    using Matrix3xN    = Eigen::Matrix< double, 3, nStrain >;
    using MatrixNxN    = Eigen::Matrix< double, nStrain, nStrain >;
    using MatrixNx3    = Eigen::Matrix< double, nStrain, 3 >;
    using MatrixNx21   = Eigen::Matrix< double, nStrain, 21 >;
    using Matrix3x21   = Eigen::Matrix< double, 3, 21 >;
    using Matrix9x21   = Eigen::Matrix< double, 9, 21 >;
    using Matrix9x3    = Marmot::Materials::TwoSidedInterface::Matrix9x3;
    using Matrix9x9    = Marmot::Materials::TwoSidedInterface::Matrix9x9;

    /**
     * @brief Set up the operators of the layer for a generalized strain.
     *
     * @param[in] layer geometry and integration points of the layer, which has to outlive the solver
     * @param[in] x generalized strain: increments for small strains, totals for finite strains
     */
    LayerSolver( const Layer< nStrain >& layer, const Vector21d& x ) : layer_( layer ), x_( x )
    {
      using namespace Marmot::Materials::TwoSidedInterface;

      const int nq = static_cast< int >( layer.weights.size() );

      const Matrix9x9 I9   = Matrix9x9::Identity();
      const Matrix9x9 Hsrf = surfaceToGradient( layer.normal, layer.tangentialSeparation, layer.normalSeparation ) - I9;
      const Matrix9x3 Bw   = jumpToGradient( layer.normal, layer.normalSeparation );
      N9_                  = jumpToGradient( layer.normal, 1.0 );

      S_                               = layer.gradientToStrain * N9_;
      Eigen::Matrix< double, 3, 9 > Rt = Eigen::Matrix< double, 3, 9 >::Zero();
      for ( int i = 0; i < 3; i++ )
        for ( int b = 0; b < 3; b++ )
          Rt( i, 3 * i + b ) = layer.normal( b );
      St_ = Rt * layer.stressToVector;

      Bx_.resize( nq );
      EBx_.resize( nq );
      strain_.resize( nq );
      for ( int q = 0; q < nq; q++ ) {
        const double y = layer.positions[q];
        Bx_[q].setZero();
        Bx_[q].template block< 9, 3 >( 0, 0 )  = Bw;
        Bx_[q].template block< 9, 9 >( 0, 3 )  = y * I9 + 0.5 * Hsrf;
        Bx_[q].template block< 9, 9 >( 0, 12 ) = ( 1. - y ) * I9 + 0.5 * Hsrf;
        EBx_[q]                                = layer.gradientToStrain * Bx_[q];
        strain_[q]                             = layer.strainOffset + EBx_[q] * x;
      }

      dA_ = Eigen::Map< const Matrix3dRowMajor >( x.data() + 3 ) -
            Eigen::Map< const Matrix3dRowMajor >( x.data() + 12 );
    }

    /**
     * @brief Solve the microproblem and the closure of the tangential gradient.
     *
     * @tparam EvaluateBulk callable `PointResponse< nStrain >( int q, const StrainVector& strain )` that updates the
     * bulk material of the point q for the strain-like measure from the committed state
     * @param[in] initialCorrections first estimate of the corrections of all points
     * @param[in] tolerance accepted component of the Newton step of the microproblem
     * @param[in] evaluateBulk update of the bulk material of a point
     *
     * @throws StressUpdateFailed if the iteration does not converge or an acoustic tensor is singular
     */
    template < class EvaluateBulk >
    void solve( const std::vector< Eigen::Vector3d >& initialCorrections,
                double                                tolerance,
                EvaluateBulk&&                        evaluateBulk )
    {
      const int nq = static_cast< int >( layer_.weights.size() );

      const auto makeEvaluate = [&]( const ClosureTerms& terms ) {
        return [&, terms]( int q, const Eigen::Vector3d& w ) {
          auto bulk = evaluateBulk( q, StrainVector( strain_[q] + terms.strain[q] + S_ * w ) );

          Point< nStrain > point;
          static_cast< PointResponse< nStrain >& >( point ) = std::move( bulk );
          point.w                                           = w;
          point.a                                           = St_ * point.stress + terms.force[q];

          point.Qinv = ( St_ * point.D * S_ ).inverse();
          if ( !point.Qinv.allFinite() )
            throw StressUpdateFailed( MakeString() << __PRETTY_FUNCTION__
                                                   << ": the acoustic tensor of an integration point through the "
                                                      "thickness is singular" );
          return point;
        };
      };

      const bool closureActive = layer_.tangentialGradient && dA_.cwiseAbs().maxCoeff() > 0.;

      ClosureTerms none;
      none.strain.assign( nq, StrainVector::Zero() );
      none.force.assign( nq, Eigen::Vector3d::Zero() );

      std::vector< Matrix3xN > dw( nq );
      tangentPoint_.resize( nq );

      ClosureTerms terms = none;
      if ( closureActive ) {
        // the derivative of the correction of the state at the start of the increment starts the iteration
        const auto                      evaluateNone = makeEvaluate( none );
        std::vector< Point< nStrain > > start;
        start.reserve( nq );
        for ( int q = 0; q < nq; q++ )
          start.push_back( evaluateNone( q, initialCorrections[q] ) );

        dw = condenseOverStrain( start );
        for ( int q = 0; q < nq; q++ )
          tangentPoint_[q] = start[q].D;
        terms = closure( dw, tangentPoint_, dA_ );
      }

      for ( int pass = 0; pass < maxClosurePasses; pass++ ) {
        points_ = ThicknessMicroproblem::solve< Point< nStrain > >( layer_.weights,
                                                                    makeEvaluate( terms ),
                                                                    initialCorrections,
                                                                    tolerance );
        dw      = condenseOverStrain( points_ );
        for ( int q = 0; q < nq; q++ )
          tangentPoint_[q] = points_[q].D;

        if ( !closureActive )
          break;

        const ClosureTerms updated = closure( dw, tangentPoint_, dA_ );

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
      closureStrainTangent_.assign( nq, MatrixNx21::Zero() );
      closureForceTangent_.assign( nq, Matrix3x21::Zero() );
      for ( int k = 0; k < ( layer_.tangentialGradient ? 9 : 0 ); k++ ) {
        Eigen::Matrix3d unit        = Eigen::Matrix3d::Zero();
        unit( k / 3, k % 3 )        = 1.0;
        const ClosureTerms response = closure( dw, tangentPoint_, unit );
        for ( int q = 0; q < nq; q++ ) {
          closureStrainTangent_[q].col( 3 + k )  = response.strain[q];
          closureStrainTangent_[q].col( 12 + k ) = -response.strain[q];
          closureForceTangent_[q].col( 3 + k )   = response.force[q];
          closureForceTangent_[q].col( 12 + k )  = -response.force[q];
        }
      }
    }

    /// Responses of the integration points at the solution.
    const std::vector< Point< nStrain > >& points() const { return points_; }

    /// Derivative of the gradient of the point q with respect to the generalized strain, 9 x 21.
    const Matrix9x21& gradientOperator( int q ) const { return Bx_[q]; }

    /// Derivative of the strain-like measure of the point q with respect to the generalized strain, with the closure.
    MatrixNx21 strainOperator( int q ) const { return EBx_[q] + closureStrainTangent_[q]; }

    /// Operator S that maps a correction to the strain-like measure.
    const MatrixNx3& correctionOperator() const { return S_; }

    /// Operator that maps the correction to the 9-vector of the gradient, the normal fibre n (x) w.
    const Matrix9x3& correctionGradientOperator() const { return N9_; }

    /**
     * @brief Derivative of the corrections with respect to the generalized strain and further variables.
     *
     * @param[in] extraTraction for every point, the derivative of its traction with respect to the further
     * variables, 3 x nExtra, or an empty vector if there are none
     * @return for every point, the derivative of its correction with respect to ( x, further variables ),
     * 3 x ( 21 + nExtra )
     */
    std::vector< Eigen::MatrixXd > condensedCorrections( const std::vector< Eigen::MatrixXd >& extraTraction ) const
    {
      const int nq     = static_cast< int >( layer_.weights.size() );
      const int nExtra = extraTraction.empty() ? 0 : static_cast< int >( extraTraction[0].cols() );

      std::vector< Eigen::Matrix< double, 3, Eigen::Dynamic > > B( nq );
      for ( int q = 0; q < nq; q++ ) {
        B[q].resize( 3, 21 + nExtra );
        B[q].leftCols( 21 ) = St_ * points_[q].D * ( EBx_[q] + closureStrainTangent_[q] ) + closureForceTangent_[q];
        if ( nExtra > 0 )
          B[q].rightCols( nExtra ) = extraTraction[q];
      }
      const auto dw = ThicknessMicroproblem::condensedCorrectionTangents< Eigen::Dynamic >( layer_.weights,
                                                                                            points_,
                                                                                            B );

      std::vector< Eigen::MatrixXd > result( nq );
      for ( int q = 0; q < nq; q++ )
        result[q] = dw[q];
      return result;
    }

    /// Generalized stress ( force, S+, S- ) of the points at the solution.
    Vector21d stress() const
    {
      Vector21d p = Vector21d::Zero();
      for ( std::size_t q = 0; q < points_.size(); q++ )
        p += layer_.thickness * layer_.weights[q] * Bx_[q].transpose() * ( layer_.stressToVector * points_[q].stress );
      return p;
    }

    /**
     * @brief Tangent of the generalized stress over the generalized strain.
     *
     * @param[in] dw derivative of the corrections, as returned by condensedCorrections()
     * @return the condensed tangent, 21 x 21
     */
    Matrix21d tangent( const std::vector< Eigen::MatrixXd >& dw ) const
    {
      Matrix21d K = Matrix21d::Zero();
      for ( std::size_t q = 0; q < points_.size(); q++ )
        K += layer_.thickness * layer_.weights[q] * Bx_[q].transpose() * layer_.stressToVector * points_[q].D *
             ( EBx_[q] + closureStrainTangent_[q] + S_ * dw[q].leftCols( 21 ) );
      return K;
    }

  private:
    /// Terms of the closure of the tangential gradient of the correction at all integration points.
    struct ClosureTerms {
      std::vector< StrainVector >    strain; ///< strain of the tangential gradient of the correction
      std::vector< Eigen::Vector3d > force;  ///< change of the traction from the bottom to the point
    };

    /// Derivative of the corrections with respect to the strain of the layer, over the points.
    std::vector< Matrix3xN > condenseOverStrain( const std::vector< Point< nStrain > >& points ) const
    {
      const int nq = static_cast< int >( points.size() );

      std::vector< Matrix3xN > B( nq );
      for ( int q = 0; q < nq; q++ )
        B[q] = St_ * points[q].D;
      return ThicknessMicroproblem::condensedCorrectionTangents< nStrain >( layer_.weights, points, B );
    }

    /// Terms of the closure, linear in the difference of the surface gradients, for the derivative of the
    /// correction and the tangent of the points.
    ClosureTerms closure( const std::vector< Matrix3xN >& dw,
                          const std::vector< MatrixNxN >& D,
                          const Eigen::Matrix3d&          surfaceGradientDifference ) const
    {
      const int    nq  = static_cast< int >( layer_.weights.size() );
      const auto&  n   = layer_.normal;
      const double h   = layer_.thickness;
      const double ell = layer_.normalSeparation;

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
        const StrainVector s = layer_.gradientToStrain * gradient;

        for ( int q = 0; q < nq; q++ ) {
          const Eigen::Vector3d dwr          = dw[q] * s;
          derivativeOfCorrection[q].col( r ) = dwr;

          // derivative of the stress along e_r, of which the divergence drives the traction through the thickness
          const Eigen::Matrix< double, 9, 1 > dStress = layer_.stressToVector * ( D[q] * ( s + S_ * dwr ) );
          for ( int k = 0; k < 3; k++ )
            source[q]( k ) += h * dStress( 3 * k + r );
        }
      }

      for ( int q = 0; q < nq; q++ ) {
        Eigen::Matrix3d gradient = Eigen::Matrix3d::Zero();
        for ( int p = 0; p < nq; p++ ) {
          terms.force[q] += ( *layer_.cumulativeWeights )( q, p ) * source[p];
          gradient += ( *layer_.cumulativeWeights )( q, p ) * h * derivativeOfCorrection[p];
        }
        // the tangential gradient of u1 as a gradient: entry (k, l) is the derivative of u1_k along e_l
        Eigen::Matrix< double, 9, 1 > flat;
        for ( int k = 0; k < 3; k++ )
          for ( int l = 0; l < 3; l++ )
            flat( 3 * k + l ) = gradient( k, l );
        terms.strain[q] = layer_.gradientToStrain * flat;
      }
      return terms;
    }

    const Layer< nStrain >& layer_; ///< geometry and integration points of the layer
    Vector21d               x_;     ///< generalized strain
    Eigen::Matrix3d         dA_;    ///< difference A+ - A- of the surface gradients

    MatrixNx3 S_;  ///< operator of the correction
    Matrix9x3 N9_; ///< normal fibre operator
    Matrix3xN St_; ///< operator of the traction

    std::vector< Matrix9x21 >   Bx_;     ///< derivative of the gradient of a point with respect to x
    std::vector< MatrixNx21 >   EBx_;    ///< derivative of the strain-like measure of a point with respect to x
    std::vector< StrainVector > strain_; ///< strain-like measure of a point without correction and closure

    std::vector< Point< nStrain > > points_;       ///< responses of the points at the solution
    std::vector< MatrixNxN >        tangentPoint_; ///< tangent of the points at the solution

    std::vector< MatrixNx21 > closureStrainTangent_; ///< derivative of the closure strain with respect to x
    std::vector< Matrix3x21 > closureForceTangent_;  ///< derivative of the closure force with respect to x
  };

  /**
   * @brief Solve the microproblem of a layer and compute its generalized stress and tangent.
   *
   * @tparam nStrain number of components of the strain-like measure of a point
   * @tparam EvaluateBulk callable `PointResponse< nStrain >( int q, const StrainVector& strain )`
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
    LayerSolver< nStrain > solver( layer, x );
    solver.solve( initialCorrections, tolerance, std::forward< EvaluateBulk >( evaluateBulk ) );

    const auto dw = solver.condensedCorrections( {} );

    Solution< nStrain > solution;
    solution.points = solver.points();
    solution.p      = solver.stress();
    solution.K      = solver.tangent( dw );
    return solution;
  }

} // namespace Marmot::Materials::TwoSidedGradedLayer
