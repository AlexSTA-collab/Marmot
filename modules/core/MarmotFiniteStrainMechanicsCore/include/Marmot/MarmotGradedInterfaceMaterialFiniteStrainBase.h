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
 * Matthias Neuner matthias.neuner@uibk.ac.at
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

#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotGradedFiniteStrainLayer.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"

#include <string>

/**
 * @brief Integration points through the thickness of a layer with finite-strain bulk materials.
 *
 * This is the part that all graded finite-strain interface materials share (see
 * MarmotGradedInterfaceMaterialFiniteStrain and MarmotGradedEquilibratedInterfaceMaterialFiniteStrain), with the
 * update of a single point added to the points of MarmotGradedFiniteStrainLayer. The state variables are, for every
 * integration point through the thickness, the correction @f$w@f$ of its last converged state, followed by the
 * state variables of its bulk material.
 */
class MarmotGradedInterfaceMaterialFiniteStrainBase
  : public MarmotGradedFiniteStrainLayer< MarmotMaterialFiniteStrain,
                                          MarmotLibrary::MarmotMaterialFiniteStrainFactory > {

public:
  /**
   * @brief Construct the integration points through the thickness from a registered finite-strain bulk material.
   *
   * @param[in] materialName name of the bulk material, registered in the MarmotMaterialFiniteStrainFactory
   * @param[in] matProperties material properties, see Marmot::Materials::ThicknessProfile
   * @param[in] nMaterialProperties number of material properties
   * @param[in] materialNumber_ number of the material
   *
   * @throws std::invalid_argument if the properties are inconsistent or the bulk material is not registered
   */
  MarmotGradedInterfaceMaterialFiniteStrainBase( const std::string& materialName,
                                                 const double*      matProperties,
                                                 int                nMaterialProperties,
                                                 int                materialNumber_ )
    : MarmotGradedFiniteStrainLayer( materialName, matProperties, nMaterialProperties, materialNumber_ )
  {
  }

  using MarmotGradedFiniteStrainLayer::getDensity;

  /**
   * @brief Get the density of the material.
   *
   * @return through-thickness mean of the density of the bulk materials
   */
  double getDensity() const { return getDensity( nullptr ); }

protected:
  /// Response of a bulk material at an integration point to a deformation gradient.
  struct BulkResponse {
    Marmot::FastorStandardTensors::Tensor33d   P; ///< first Piola-Kirchhoff stress
    Marmot::FastorStandardTensors::Tensor3333d A; ///< derivative of the stress with respect to F
    std::vector< double > state; ///< correction and state variables of the bulk material after the update
  };

  /**
   * @brief Update the bulk material of a point for a deformation gradient, starting from the committed state.
   *
   * The state variables passed are not modified, the updated ones are returned in the response, where the entries
   * of the correction are the committed ones.
   *
   * @param[in] q index of the integration point
   * @param[in] stateVars state variables of the material point at the beginning of the increment
   * @param[in] F total deformation gradient of the point
   * @param[in] timeIncrement time at the beginning of the increment and time increment
   * @return stress, tangent and updated state of the point
   */
  BulkResponse evaluateBulk( int                                              q,
                             const double*                                    stateVars,
                             const Marmot::FastorStandardTensors::Tensor33d&  F,
                             const MarmotMaterialFiniteStrain::TimeIncrement& timeIncrement ) const;
};
