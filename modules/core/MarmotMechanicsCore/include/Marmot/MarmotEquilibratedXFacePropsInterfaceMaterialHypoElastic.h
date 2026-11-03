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

#include "Marmot/MarmotEquilibratedXInterfaceMaterialHypoElastic.h"

#include <string>

/**
 * @brief Equilibrated two-sided hypoelastic interface material whose two faces carry different bulk properties.
 *
 * The material of MarmotEquilibratedXInterfaceMaterialHypoElastic hands the same properties to the bulk materials
 * of the top and of the bottom face, so the tangents of the two faces differ only through their histories. This
 * material lets them differ from the start, which describes a layer whose two halves are different materials, for
 * instance a biphasic layer. It adds nothing else: the local solve, the gauge and the condensed tangent are those of
 * the base class, and the mismatch of the face tangents enters through them.
 *
 * The material properties are
 * [E+, nu+, h, properties of the base material..., E-, nu-]
 * with the first part as in the base class for the top face (the face of the surface gradient A+) and the last two
 * entries the Young's modulus and Poisson's ratio of the bottom face. The base material of the bottom face receives
 * E-, nu- and the properties after h.
 */
class MarmotEquilibratedXFacePropsInterfaceMaterialHypoElastic : public MarmotEquilibratedXInterfaceMaterialHypoElastic {

public:
  /**
   * @brief Construct the material from a registered hypoelastic base material.
   *
   * @param[in] materialName name of the base material, registered in the MarmotMaterialHypoElasticFactory
   * @param[in] matProperties material properties, see above
   * @param[in] nMaterialProperties number of material properties, at least five
   * @param[in] materialNumber_ number of the material
   *
   * @throws std::invalid_argument if fewer than five properties are given, if h is not positive, or if the base
   * material is not registered
   */
  MarmotEquilibratedXFacePropsInterfaceMaterialHypoElastic( const std::string& materialName,
                                                            const double*      matProperties,
                                                            int                nMaterialProperties,
                                                            int                materialNumber_ );
};
