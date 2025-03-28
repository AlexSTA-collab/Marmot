/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck
 * 2020 - today
 *
 * festigkeitslehre@uibk.ac.at
 *
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
#include "Marmot/MarmotMaterialHypoElasticInterface.h"
#include "Marmot/MarmotTypedefs.h"
#include <iostream>
#include <string>
#include <vector>
#include "Fastor/Fastor.h"

namespace Marmot::Materials {
  /**
   * \brief Implementation of a linear elastic interface material
   * for 3D stress states.
   *
   * For further information see \ref linearelasticinterface.
   */
  class LinearElasticInterface : public MarmotMaterialHypoElasticInterface {
  public:
    using MarmotMaterialHypoElasticInterface::MarmotMaterialHypoElasticInterface;
    using Tensor1D = Fastor::Tensor<double,3>;
    using Tensor2D = Fastor::Tensor<double,3,3>;

    LinearElasticInterface( const double* materialProperties, int nMaterialProperties, int materialNumber );

    void computeStress( Tensor1D&  force,
                        Tensor2D&  surface_stress,
                        Fastor::Tensor<double, 21,21>& dStress_dStrain,
                        const Fastor::Tensor<double, 6,1>& dU,
                        const Fastor::Tensor<double, 18,1>& dSurface_strain,
                        const Tensor1D& normal,
                        const double* timeOld,
                        const double  dT,
                        double&       pNewDT ) ;

    StateView getStateView( const std::string& result ) { return { nullptr, 0 }; };

    int getNumberOfRequiredStateVars() { return 0; }

    double getDensity();
  };
} // namespace Marmot::Materials
