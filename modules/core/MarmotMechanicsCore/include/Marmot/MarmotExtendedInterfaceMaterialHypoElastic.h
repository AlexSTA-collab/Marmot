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
 * ---------------------------------------------------------------------
 */

#pragma once

#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotStateHelpers.h"

#include <memory>
#include <string>
#include <vector>

/**
 * Extended hypoelastic interface material with independent top and bottom
 * interphase material points.
 *
 * Preferred property layout:
 *   [ h, nBottom, bottomProperties..., nTop, topProperties... ]
 *
 * A legacy single-material layout [ E, nu, h, ... ] is accepted as a fallback.
 * It assigns the same initial material properties to top and bottom, but still
 * creates two independent material instances with separate state-variable
 * storage.
 */
class MarmotExtendedInterfaceMaterialHypoElastic {

protected:
  const double* materialProperties;
  const int     nMaterialProperties;
  double        h = 0.0;

  std::vector< double > bottomMaterialProperties;
  std::vector< double > topMaterialProperties;

  std::unique_ptr< MarmotMaterialHypoElastic > bottomMaterial;
  std::unique_ptr< MarmotMaterialHypoElastic > topMaterial;

public:
  using TensorMap3d  = Marmot::FastorStandardTensors::TensorMap3d;
  using TensorMap33d = Marmot::FastorStandardTensors::TensorMap33d;
  using TensorMap6d  = Marmot::FastorStandardTensors::TensorMap6d;
  using TensorMap9d  = Marmot::FastorStandardTensors::TensorMap9d;
  using TensorMap18d = Marmot::FastorStandardTensors::TensorMap18d;

  const int materialNumber;

  MarmotExtendedInterfaceMaterialHypoElastic( const std::string& materialName,
                                              const double*      matProperties_,
                                              int                nMaterialProperties_,
                                              int                materialNumber_ );

  virtual ~MarmotExtendedInterfaceMaterialHypoElastic() = default;

  MarmotStateLayoutDynamic stateLayout;

  double characteristicElementLength;

  void setCharacteristicElementLength( double length );

  struct State {
    TensorMap3d  force;
    TensorMap33d averageSurfaceStress;
    TensorMap33d jumpSurfaceStress;
    double*      stateVars;
  };

  struct Tangents {
    double* forceJumpU;
    double* forceAverageSurfaceGradient;
    double* forceJumpSurfaceGradient;

    double* averageSurfaceStressJumpU;
    double* averageSurfaceStressAverageSurfaceGradient;
    double* averageSurfaceStressJumpSurfaceGradient;

    double* jumpSurfaceStressJumpU;
    double* jumpSurfaceStressAverageSurfaceGradient;
    double* jumpSurfaceStressJumpSurfaceGradient;
  };

  struct Deformation {
    TensorMap6d  dU;
    TensorMap18d dSurfaceStrain;
    TensorMap3d  normal;

    Deformation( const double* dU_, const double* dSurfaceStrain_, const double* normal_ )
      : dU( const_cast< double* >( dU_ ) ),
        dSurfaceStrain( const_cast< double* >( dSurfaceStrain_ ) ),
        normal( const_cast< double* >( normal_ ) )
    {
    }
  };

  struct TimeIncrement {
    double timeOld;
    double dT;
  };

  virtual void computeStress( State&               state,
                              Tangents&            tangents,
                              const Deformation&   deformation,
                              const TimeIncrement& timeIncrement );

  MarmotMaterialHypoElastic& getBottomMaterial() { return *bottomMaterial; }

  MarmotMaterialHypoElastic& getTopMaterial() { return *topMaterial; }

  StateView getStateView( const std::string& stateName, double* stateVars ) const
  {
    return stateLayout.getStateView( stateVars, stateName );
  }

  virtual int getNumberOfRequiredStateVars() const { return stateLayout.totalSize(); }

  virtual void initializeYourself( double* stateVars, int nStateVars );

  virtual double getDensity();
};
