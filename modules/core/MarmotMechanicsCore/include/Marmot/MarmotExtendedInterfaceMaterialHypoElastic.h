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
 * --------------------------------------------------------------------- */

#pragma once

#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotStateHelpers.h"

#include <memory>
#include <string>
#include <vector>

/**
 * Extended hypoelastic interface material with two independent material
 * points and a locally condensed moving bilinear kink.
 *
 * The top material state occupies alpha*h and the bottom material state
 * occupies (1-alpha)*h.  At every constitutive update the local variables
 *
 *   g     = [u_{,n}],
 *   alpha = h_top / h,
 *
 * are condensed from the incremental potential of the two sublayers.
 * Stationarity with respect to g enforces traction continuity, while
 * stationarity with respect to alpha determines the kink position.
 *
 * The implementation reconstructs the condensed incremental potential of
 * the underlying material by integrating the stress response along the
 * straight strain-increment path.  This is appropriate for potential-based
 * algorithmic updates such as associative Von Mises plasticity and linear
 * elasticity.  It must not be used unchanged for genuinely non-associated
 * models whose algorithmic stress map is not potential-derived.
 *
 * Preferred property layout:
 *   [ h, nBottom, bottomProperties..., nTop, topProperties... ]
 *
 * Legacy single-material layout:
 *   [ E, nu, h, remainingBaseMaterialProperties... ]
 *
 * The two layouts are distinguished structurally: the explicit layout is
 * selected iff the second property is an integer sublayer property count
 * >= 1, which no physically valid Poisson's ratio (nu < 1, including
 * nu == 0) can be.
 */
class MarmotExtendedInterfaceMaterialHypoElastic {

protected:
  const double* materialProperties;
  const int     nMaterialProperties;
  double        h = 0.0;

  std::string materialName;

  std::vector< double > bottomMaterialProperties;
  std::vector< double > topMaterialProperties;

  std::unique_ptr< MarmotMaterialHypoElastic > bottomMaterial;
  std::unique_ptr< MarmotMaterialHypoElastic > topMaterial;

public:
  using TensorMap3d  = Marmot::FastorStandardTensors::TensorMap3d;
  using TensorMap33d = Marmot::FastorStandardTensors::TensorMap33d;
  using TensorMap6d  = Marmot::FastorStandardTensors::TensorMap6d;
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