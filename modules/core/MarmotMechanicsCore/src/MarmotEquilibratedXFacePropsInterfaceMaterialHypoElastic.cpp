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

#include "Marmot/MarmotEquilibratedXFacePropsInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

#include <stdexcept>
#include <vector>

MarmotEquilibratedXFacePropsInterfaceMaterialHypoElastic::MarmotEquilibratedXFacePropsInterfaceMaterialHypoElastic(
  const std::string& materialName,
  const double*      matProperties,
  int                nMaterialProperties,
  int                materialNumber_ )
  : MarmotEquilibratedXInterfaceMaterialHypoElastic( materialName,
                                                     matProperties,
                                                     nMaterialProperties >= 5 ? nMaterialProperties - 2
                                                                              : nMaterialProperties,
                                                     materialNumber_ )
{
  if ( nMaterialProperties < 5 ) {
    throw std::invalid_argument( "MarmotEquilibratedXFacePropsInterfaceMaterialHypoElastic requires at least five "
                                 "properties: E+, nu+, h, E-, nu-." );
  }

  std::vector< double > bottomProperties;
  bottomProperties.reserve( nMaterialProperties - 3 );
  bottomProperties.push_back( matProperties[nMaterialProperties - 2] );
  bottomProperties.push_back( matProperties[nMaterialProperties - 1] );
  bottomProperties.insert( bottomProperties.end(), matProperties + 3, matProperties + nMaterialProperties - 2 );

  bottomMaterial = std::unique_ptr< MarmotMaterialHypoElastic >(
    MarmotLibrary::MarmotMaterialHypoElasticFactory::createMaterial( materialName,
                                                                     bottomProperties.data(),
                                                                     static_cast< int >( bottomProperties.size() ),
                                                                     materialNumber ) );
}
