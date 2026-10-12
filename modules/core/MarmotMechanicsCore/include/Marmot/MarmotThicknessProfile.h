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

#include <vector>

namespace Marmot::Materials {

  /**
   * @brief Integration points through the thickness of a layer, whose bulk properties vary with the position.
   *
   * The layer has the through-thickness coordinate @f$y \in [0,1]@f$. The material properties of the layer are
   *
   * | Index | Name | Meaning |
   * |-------|------|---------|
   * | 0 | h | thickness of the layer |
   * | 1 | mode | 0: layers with constant properties, 1: properties interpolated linearly between stations, 2: points at
   * the stations | | 2 | nGauss | number of Gauss points per layer (mode 0) or per interval between two stations (mode
   * 1), not used in mode 2 | | 3 | K | number of layers (mode 0) or of stations (modes 1 and 2, at least two) | | 4 | m
   * | number of properties of the bulk material |
   *
   * followed by K groups of the entry @f$p_k@f$ and the m properties of the bulk material. In the mode 0, @f$p_k@f$
   * is the relative thickness of the layer k, which are ordered from the bottom to the top and normalized by their
   * sum. In the modes 1 and 2, @f$p_k@f$ is the position of the station k, with @f$p_0 = 0@f$, @f$p_{K-1} = 1@f$ and
   * increasing positions in between. In the mode 2, the integration points are the stations themselves, with the
   * weights of the trapezoidal rule, and the properties of a point are the ones given at its station. This is the
   * layout of the faces of a two-sided layer and of the nonlocal fields at the stations.
   */
  struct ThicknessProfile {

    /// Integration point through the thickness.
    struct Point {
      double                weight;     ///< weight, the weights of all points sum up to one
      double                position;   ///< coordinate y in [0, 1]
      std::vector< double > properties; ///< properties of the bulk material at the position
    };

    double               thickness = 0.; ///< thickness h of the layer
    std::vector< Point > points;         ///< integration points, from the bottom to the top

    /**
     * @brief Read the integration points from the material properties of the layer.
     *
     * @param[in] materialProperties material properties, as in the table of this class
     * @param[in] nMaterialProperties number of material properties
     * @return the thickness and the integration points of the layer
     *
     * @throws std::invalid_argument if the properties are inconsistent, for example a non-positive thickness or an
     * unknown mode, or if the positions of the stations are not increasing from 0 to 1
     */
    static ThicknessProfile fromMaterialProperties( const double* materialProperties, int nMaterialProperties );
  };

  /**
   * @brief Nodes and weights of the Gauss-Legendre rule on [-1, 1].
   *
   * @param[in] n number of points, at least one
   * @return the nodes and the weights, which sum up to 2
   */
  std::pair< std::vector< double >, std::vector< double > > gaussLegendre( int n );

} // namespace Marmot::Materials
