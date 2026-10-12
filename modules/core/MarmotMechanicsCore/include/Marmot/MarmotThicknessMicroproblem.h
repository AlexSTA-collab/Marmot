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

#include <Eigen/Dense>

#include <algorithm>
#include <vector>

namespace Marmot::Materials::ThicknessMicroproblem {

  /**
   * @brief Maximum number of iterations of the microproblem.
   */
  constexpr int maxIterations = 60;

  /**
   * @brief Maximum number of halvings of a step of the line search.
   */
  constexpr int maxBacktracks = 12;

  /**
   * @brief Solve the microproblem through the thickness of a layer by a damped Newton iteration.
   *
   * At every integration point q of the thickness the strain is the strain of the layer plus the correction
   * @f$ \mathrm{sym}( w_q \otimes n ) @f$ (small strains) or @f$ w_q \otimes N @f$ (finite strains), and the
   * traction @f$ a_q @f$ of the point is its stress contracted with the normal. The microproblem requires that
   *
   * - the traction is the same at all points, @f$ a_q = t @f$, and
   * - the corrections vanish on average, @f$ \sum_q \omega_q w_q = 0 @f$,
   *
   * for the unknown corrections @f$ w_q @f$ and the unknown common traction @f$ t @f$. With the inverse
   * @f$ Q_q^{-1} @f$ of the acoustic tensor of the point, a Newton step is
   * @f[
   *   \delta t = \Big( \sum_q \omega_q Q_q^{-1} \Big)^{-1}
   *     \Big( \sum_q \omega_q Q_q^{-1} ( a_q - t ) - \sum_q \omega_q w_q \Big) , \qquad
   *   \delta w_q = Q_q^{-1} ( \delta t - ( a_q - t ) ) .
   * @f]
   * The step is damped by a backtracking line search on the squared norm of the residual, in which the tractions
   * are scaled to strains by the stiffest acoustic tensor. A bulk material that fails to update at a trial step
   * makes the step shorter. A linear material converges in one step.
   *
   * @tparam Point response of an integration point, which has the members `w` (correction), `a` (traction) and
   * `Qinv` (inverse of the acoustic tensor), as Eigen::Vector3d, Eigen::Vector3d and Eigen::Matrix3d
   * @tparam Evaluate callable `Point( int q, const Eigen::Vector3d& w )` that evaluates the point q for the
   * correction w, starting from the state of the previous increment
   * @param[in] weights weights of the integration points
   * @param[in] evaluate evaluation of an integration point
   * @param[in] initialCorrections first estimate of the corrections of all points
   * @param[in] tolerance largest accepted component of the Newton step, in units of strain
   * @return the responses of all integration points at the solution
   *
   * @throws StressUpdateFailed if the iteration does not converge
   */
  template < class Point, class Evaluate >
  std::vector< Point > solve( const std::vector< double >&          weights,
                              Evaluate&&                            evaluate,
                              const std::vector< Eigen::Vector3d >& initialCorrections,
                              double                                tolerance )
  {
    const int nq = static_cast< int >( weights.size() );

    std::vector< Point > evaluations;
    evaluations.reserve( nq );
    for ( int q = 0; q < nq; q++ )
      evaluations.push_back( evaluate( q, initialCorrections[q] ) );

    // scale of the acoustic tensors, which turns a misfit of the tractions into a strain
    double stiffness = 0.;
    for ( const auto& evaluation : evaluations )
      stiffness = std::max( stiffness, evaluation.Qinv.inverse().cwiseAbs().maxCoeff() );

    Eigen::Matrix3d M   = Eigen::Matrix3d::Zero();
    Eigen::Vector3d rhs = Eigen::Vector3d::Zero();
    for ( int q = 0; q < nq; q++ ) {
      M += weights[q] * evaluations[q].Qinv;
      rhs += weights[q] * evaluations[q].Qinv * evaluations[q].a;
    }
    // the first traction is the one of the linearization at the initial corrections
    Eigen::Vector3d t = M.partialPivLu().solve( rhs );

    const auto merit = [&]( const std::vector< Point >& trial, const Eigen::Vector3d& tTrial ) {
      double          sum  = 0.;
      Eigen::Vector3d mean = Eigen::Vector3d::Zero();
      for ( int q = 0; q < nq; q++ ) {
        sum += weights[q] * ( ( trial[q].a - tTrial ) / stiffness ).squaredNorm();
        mean += weights[q] * trial[q].w;
      }
      return 0.5 * ( sum + mean.squaredNorm() );
    };

    for ( int iteration = 0; iteration < maxIterations; iteration++ ) {

      Eigen::Vector3d mean = Eigen::Vector3d::Zero();
      M.setZero();
      rhs.setZero();
      for ( int q = 0; q < nq; q++ ) {
        mean += weights[q] * evaluations[q].w;
        M += weights[q] * evaluations[q].Qinv;
        rhs += weights[q] * evaluations[q].Qinv * ( evaluations[q].a - t );
      }
      const Eigen::Vector3d dt = M.partialPivLu().solve( rhs - mean );

      std::vector< Eigen::Vector3d > steps( nq );
      double                         maxStep = 0.;
      for ( int q = 0; q < nq; q++ ) {
        steps[q] = evaluations[q].Qinv * ( dt - ( evaluations[q].a - t ) );
        maxStep  = std::max( maxStep, steps[q].cwiseAbs().maxCoeff() );
      }

      if ( maxStep <= tolerance )
        return evaluations;

      const double meritBefore = merit( evaluations, t );

      double alpha = 1.0;
      for ( int backtrack = 0;; backtrack++ ) {
        try {
          std::vector< Point > trial;
          trial.reserve( nq );
          for ( int q = 0; q < nq; q++ )
            trial.push_back( evaluate( q, evaluations[q].w + alpha * steps[q] ) );

          if ( merit( trial, t + alpha * dt ) <= ( 1. - 2e-4 * alpha ) * meritBefore || backtrack >= maxBacktracks ) {
            evaluations = std::move( trial );
            t += alpha * dt;
            break;
          }
        }
        catch ( const std::exception& ) {
          // a bulk material could not update the strain of this trial: take a shorter step
          if ( backtrack >= maxBacktracks )
            throw;
        }
        alpha *= 0.5;
      }
    }

    throw StressUpdateFailed( MakeString()
                              << __PRETTY_FUNCTION__ << ": the microproblem through the thickness did not converge" );
  }

  /**
   * @brief Derivative of the corrections with respect to the strain of the layer, with the corrections condensed.
   *
   * The traction equilibrium of a point is @f$ Q_q \, \delta w_q + B_q \, \delta E = \delta t @f$, where @f$ B_q @f$
   * is the derivative of the traction of the point with respect to the strain of the layer, and the corrections
   * vanish on average. The result is
   * @f[
   *   \delta t = \Big( \sum_q \omega_q Q_q^{-1} \Big)^{-1} \sum_q \omega_q Q_q^{-1} B_q \, \delta E , \qquad
   *   \frac{\partial w_q}{\partial E} = Q_q^{-1} ( \partial_E t - B_q ) .
   * @f]
   *
   * @tparam nStrain number of components of the strain of the layer
   * @tparam Point response of an integration point, with the member `Qinv`
   * @param[in] weights weights of the integration points
   * @param[in] evaluations responses of the integration points at the solution of the microproblem
   * @param[in] B derivative of the traction of every point with respect to the strain of the layer
   * @return the derivative of the correction of every point with respect to the strain of the layer
   */
  template < int nStrain, class Point >
  std::vector< Eigen::Matrix< double, 3, nStrain > > condensedCorrectionTangents(
    const std::vector< double >&                              weights,
    const std::vector< Point >&                               evaluations,
    const std::vector< Eigen::Matrix< double, 3, nStrain > >& B )
  {
    const int nq = static_cast< int >( weights.size() );

    Eigen::Matrix3d                     M = Eigen::Matrix3d::Zero();
    Eigen::Matrix< double, 3, nStrain > G = Eigen::Matrix< double, 3, nStrain >::Zero();
    for ( int q = 0; q < nq; q++ ) {
      M += weights[q] * evaluations[q].Qinv;
      G += weights[q] * evaluations[q].Qinv * B[q];
    }
    const Eigen::Matrix< double, 3, nStrain > dt = M.partialPivLu().solve( G );

    std::vector< Eigen::Matrix< double, 3, nStrain > > dw( nq );
    for ( int q = 0; q < nq; q++ )
      dw[q] = evaluations[q].Qinv * ( dt - B[q] );
    return dw;
  }

} // namespace Marmot::Materials::ThicknessMicroproblem
