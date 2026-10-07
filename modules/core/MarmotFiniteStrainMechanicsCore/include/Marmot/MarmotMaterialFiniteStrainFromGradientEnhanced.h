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

#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrain.h"

#include <memory>
#include <string>

/**
 * Presents a gradient-enhanced finite-strain material as a plain MarmotMaterialFiniteStrain, which holds the
 * nonlocal field at a fixed value (zero by default).
 *
 * This lets an element that carries only displacements, for instance an interface element of finite strain,
 * drive a material written against MarmotMaterialGradientEnhancedFiniteStrain.
 *
 * The presentation is exact, not an approximation, whenever the stress of the material does not depend on the
 * nonlocal field at the fixed value. Instead of trusting the material card, the adapter checks the derivative of the
 * stress with respect to the nonlocal field on every call and refuses the material if it does not vanish: dropping
 * this coupling silently would give a wrong tangent and a wrong solution, not a small error.
 *
 * What is necessarily lost is the nonlocal regularization itself: the local driving force and the nonlocal radius
 * of the material are discarded, because there is no field on the interface to receive them.
 */
class MarmotMaterialFiniteStrainFromGradientEnhanced : public MarmotMaterialFiniteStrain {

  std::unique_ptr< MarmotMaterialGradientEnhancedFiniteStrain > wrapped; ///< the gradient-enhanced material

  const double nonlocalField; ///< value at which the nonlocal field is held

public:
  /**
   * @brief Construct the adapter.
   *
   * @param[in] wrappedMaterial the gradient-enhanced material, which the adapter owns
   * @param[in] materialProperties properties of the material
   * @param[in] nMaterialProperties number of properties
   * @param[in] materialNumber number of the material
   * @param[in] nonlocalFieldValue value at which the nonlocal field is held
   *
   * @throws std::invalid_argument if there is no material to wrap
   */
  MarmotMaterialFiniteStrainFromGradientEnhanced(
    std::unique_ptr< MarmotMaterialGradientEnhancedFiniteStrain > wrappedMaterial,
    const double*                                                 materialProperties,
    int                                                           nMaterialProperties,
    int                                                           materialNumber,
    double                                                        nonlocalFieldValue = 0.0 );

  /**
   * @brief Compute the stress of the wrapped material at the fixed nonlocal field.
   *
   * @param[in,out] response the response; the Kirchhoff stress and the energies are taken from the wrapped material
   * @param[out] tangents the derivative of the stress with respect to the deformation gradient
   * @param[in] deformation the deformation gradient
   * @param[in] timeIncrement the time at the beginning of the increment and the time increment
   *
   * @throws std::invalid_argument if the stress of the wrapped material depends on the nonlocal field
   */
  void computeStress( ConstitutiveResponse< 3 >& response,
                      AlgorithmicModuli< 3 >&    tangents,
                      const Deformation< 3 >&    deformation,
                      const TimeIncrement&       timeIncrement ) const override;

  /// @copydoc MarmotMaterialFiniteStrain::initializeYourself
  void initializeYourself( double* stateVars, int nStateVars ) override;

  /// @copydoc MarmotMaterialFiniteStrain::getDensity
  double getDensity( const double* stateVars ) const override;

  /// @copydoc MarmotMaterialFiniteStrain::setCharacteristicElementLength
  void setCharacteristicElementLength( double length ) override;
};

namespace MarmotLibrary {

  /**
   * @brief Create a finite-strain material by name, falling back to the gradient-enhanced materials.
   *
   * The name is looked up among the finite-strain materials first and then among the gradient-enhanced ones, which
   * are wrapped in a MarmotMaterialFiniteStrainFromGradientEnhanced. Elements that carry only displacements can
   * therefore be handed the name of a gradient-enhanced material, and either work or say precisely why not.
   *
   * @param[in] materialName name of the material
   * @param[in] materialProperties properties of the material
   * @param[in] nMaterialProperties number of properties
   * @param[in] materialNumber number of the material
   * @return the material
   *
   * @throws std::invalid_argument if the name is registered in neither registry
   */
  std::unique_ptr< MarmotMaterialFiniteStrain > createFiniteStrainMaterialAllowingGradientEnhanced(
    const std::string& materialName,
    const double*      materialProperties,
    int                nMaterialProperties,
    int                materialNumber );

} // namespace MarmotLibrary
