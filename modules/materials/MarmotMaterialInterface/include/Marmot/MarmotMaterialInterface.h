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
 * Matthias   Neuner  matthias.neuner@boku.ac.at
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
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotTypedefs.h"
#include "interface_material_helper_functions.h"
#include "Marmot/MarmotUtils.h"
#include <string>
#include "Eigen/Dense"
#include "Eigen/Core"
#include "Fastor/Fastor.h"

using namespace Eigen;
using namespace Fastor;

using Tensor1D = Fastor::Tensor<double,3>;
using Tensor2D = Fastor::Tensor<double,3,3>;
using Tensor3D = Fastor::Tensor<double,3,3,3>;
using Tensor4D = Fastor::Tensor<double,3,3,3,3 >;


class MarmotMaterialInterface {

protected:
  const double* materialProperties;
  const int     nMaterialProperties;

  double* stateVars;
  int     nStateVars;

public:
  const int materialNumber;

  MarmotMaterialInterface( const double* materialProperties, int nMaterialProperties, int materialNumber );

  MarmotMaterialInterface();

  double* getAssignedStateVars();

  int getNumberOfAssignedStateVars(){return 0;};

  void initializeYourself();
  
  void computeStress( Tensor1D&  force,
                      Tensor2D&  surface_stress,
                      Fastor::Tensor<double, 21,21>& dStress_dStrain,
                      const Fastor::Tensor<double, 6,1>& dU,
                      const Fastor::Tensor<double, 18,1>& dSurface_strain,
                      const Tensor1D& normal,
                      const double* timeOld,
                      const double  dT,
                      double&       pNewDT  );
  
 class MarmotMaterialInterfaceStateVarManager : public MarmotStateVarVectorManager {

    public:
      inline const static auto layout = makeLayout( {
        { .name = "kappa", .length = 1 },
      } );

      double& kappa;

      MarmotMaterialInterfaceStateVarManager( double* theStateVarVector )
        : MarmotStateVarVectorManager( theStateVarVector, layout ), kappa( find( "kappa" ) ){};
    };
    std::unique_ptr< MarmotMaterialInterfaceStateVarManager > managedStateVars;

    int getNumberOfRequiredStateVars(){ return MarmotMaterialInterfaceStateVarManager::layout.nRequiredStateVars; }

    void assignStateVars( double* stateVars, int nStateVars );

    StateView getStateView( const std::string& result );

};
