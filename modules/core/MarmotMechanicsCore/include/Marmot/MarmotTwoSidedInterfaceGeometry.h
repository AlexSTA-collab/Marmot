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

#include <Eigen/Dense>

namespace Marmot::Materials::TwoSidedInterface {

  /// Row-major 3x3 matrix, whose flat index 3 i + j is the one of the 9-vectors of the two-sided interface materials.
  using Matrix3dRowMajor = Eigen::Matrix< double, 3, 3, Eigen::RowMajor >;

  /// 9x3 matrix.
  using Matrix9x3 = Eigen::Matrix< double, 9, 3 >;

  /// 9x9 matrix.
  using Matrix9x9 = Eigen::Matrix< double, 9, 9 >;

  /// Split of the connector between the two faces into its normal and its tangential part.
  struct InterfaceGeometry {
    double          normalSeparation;     ///< component of the connector along the normal
    Eigen::Vector3d tangentialSeparation; ///< component of the connector in the tangent plane
  };

  /**
   * @brief Flat row-major index of the entry (i, j) of a 3x3 tensor.
   *
   * @param[in] i row
   * @param[in] j column
   * @return the flat index
   */
  constexpr int flatIndex( int i, int j )
  {
    return 3 * i + j;
  }

  /**
   * @brief Split the connector between the two faces into its normal and its tangential part.
   *
   * @param[in] normal unit normal of the interface
   * @param[in] separationVector connector from the bottom to the top face, zero for coincident faces
   * @param[in] thickness thickness of the layer, used as the normal separation of coincident faces
   * @return the normal separation and the tangential separation
   *
   * @throws std::invalid_argument if the connector has no positive normal component
   */
  InterfaceGeometry evaluateInterfaceGeometry( const Eigen::Vector3d& normal,
                                               const Eigen::Vector3d& separationVector,
                                               double                 thickness );

  /**
   * @brief Map a vector g to vec( g (x) n ) / normalSeparation.
   *
   * @param[in] normal unit normal of the interface
   * @param[in] normalSeparation normal separation of the two faces
   * @return the 9x3 map
   */
  Matrix9x3 jumpToGradient( const Eigen::Vector3d& normal, double normalSeparation );

  /**
   * @brief Map vec( A ) to vec( A ) - vec( ( A d_tau ) (x) n ) / normalSeparation, which is the identity plus the
   * correction of the tangential separation of the faces.
   *
   * @param[in] normal unit normal of the interface
   * @param[in] tangentialSeparation tangential part of the connector
   * @param[in] normalSeparation normal separation of the two faces
   * @return the 9x9 map
   */
  Matrix9x9 surfaceToGradient( const Eigen::Vector3d& normal,
                               const Eigen::Vector3d& tangentialSeparation,
                               double                 normalSeparation );

} // namespace Marmot::Materials::TwoSidedInterface
