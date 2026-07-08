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
 * This file is part of the MAteRialMOdellingToolbox (marmot).
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of marmot.
 * --------------------------------------------------------------------- */

#pragma once

#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotStateHelpers.h"

#include <memory>
#include <string>
#include <vector>

class MarmotInterfaceMaterialFiniteStrain {
protected:
  const double*                                 materialProperties;
  const int                                     nMaterialProperties;
  double                                        h = 0.0;
  std::vector< double >                         baseMaterialProperties;
  std::unique_ptr< MarmotMaterialFiniteStrain > baseMaterial;

public:
  const int materialNumber;

  MarmotInterfaceMaterialFiniteStrain( const std::string& materialName,
                                       const double*      materialProperties,
                                       int                nMaterialProperties,
                                       int                materialNumber );

  virtual ~MarmotInterfaceMaterialFiniteStrain() = default;

  MarmotStateLayoutDynamic stateLayout;

  struct State {
    double* force;
    double* surfaceStress;
    double* stateVars;
  };

  struct Tangents {
    double* Q_ik;
    double* H_ikR;
    double* H_iRk;
    double* A_iRkP;
  };

  struct Deformation {
    const double* dU;
    const double* dSurfaceGradient;
    const double* normal;
  };

  struct TimeIncrement {
    double timeOld;
    double dT;
  };

  virtual void computeStress( State&               state,
                              Tangents&            tangents,
                              const Deformation&   deformation,
                              const TimeIncrement& timeIncrement );

  virtual int getNumberOfRequiredStateVars() const { return stateLayout.totalSize(); }

  virtual void initializeYourself( double* stateVars, int nStateVars );

  virtual double getDensity() const;
};
