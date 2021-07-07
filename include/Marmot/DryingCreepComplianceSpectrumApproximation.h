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
#pragma once
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/Solidification.h"

namespace Marmot::Materials {

  namespace SpectrumApproximation {

    namespace DryingCreepB3 {

      enum class ApproximationOrder : int { firstOrder = 1, secondOrder = 2, thirdOrder = 3 };

      void computeApproximation( Eigen::Ref< Solidification::KelvinProperties > kelvinElasticModuli,
                                 Eigen::Ref< Solidification::KelvinProperties > kelvinRetardationTimes,
                                 const double                                   hEnv,
                                 const double                                   xi,
                                 const double                                   xiZero,
                                 const enum ApproximationOrder                  approximationOrder );

      Solidification::KelvinProperties generateRetardationTimes( double minTime, int numUnits );

      template < int N >
      struct Factorial {
        enum { value = N * Factorial< N - 1 >::value };
      };

      template <>
      struct Factorial< 0 >{
        enum { value = 1 };
      };

      template < typename T_ >
      T_ phi( T_ xi, double b, double xiZero )
      {

        return sqrt( f< T_ >( xi - xiZero, b ) - f< T_ >( -xiZero, b ) );
      }

      template < typename T_ >
      T_ T( T_ eta )
      {
        return tanh( eta );
      }

      template < typename T_ >
      T_ S( T_ xi )
      {
        return T( T_( sqrt( xi ) ) );
      }

      template < typename T_ >
      T_ f( T_ xi, double b )
      {
        return exp( b * S( xi ) );
      }
      template < int k = 1 >
      double computeSpectrum( double tau, double b, double xiZero )
      {
        typedef autodiff::forward::HigherOrderDual< k > currType;
        currType                                        tau_ = k * tau;

        double val = -pow( -double( tau_.val ), k ) / Factorial< k - 1 >::value;
        val *= autodiff::forward::derivative( phi< currType >,
                                              autodiff::forward::wrt< k >( tau_ ),
                                              autodiff::forward::at( tau_, b, xiZero ) );
        return val;
      }

    } // namespace DryingCreepB3
  }   // namespace SpectrumApproximation
} // namespace Marmot::Materials
