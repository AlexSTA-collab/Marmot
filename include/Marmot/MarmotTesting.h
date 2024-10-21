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

#include "autodiff/forward/dual/dual.hpp"
#include <Eigen/Core>
#include <iostream>

namespace Marmot::Testing {

  bool checkIfEqual( const double a, const double b, const double tol = 1e-15 );

  bool checkIfEqual( const autodiff::dual a, const autodiff::dual b, const double tol = 1e-15 );

  std::string getString( const double a );
  std::string getString( const autodiff::dual a );

  template < typename T >
  bool checkIfEqual( const Eigen::Matrix< T, -1, -1 >& a,
                     const Eigen::Matrix< T, -1, -1 >& b,
                     const double                      tol = 1e-15 )
  {
    if ( a.rows() != b.rows() || a.cols() != b.cols() ) {
      return false;
    }
    for ( int i = 0; i < a.rows(); i++ ) {
      for ( int j = 0; j < a.cols(); j++ ) {
        auto cond = checkIfEqual( a( i, j ), b( i, j ), tol );
        if ( !cond ) {
          std::cout << "  -> HINT:  a(" << i << "," << j << ") = " << getString( a( i, j ) ) << " !=  b(" << i << ","
                    << j << ") =" << getString( b( i, j ) ) << std::endl;
          return false;
        }
      }
    }
    return true;
  }

  void throwExceptionOnFailure( const bool condition, const std::string& message = "" );

} // namespace Marmot::Testing
