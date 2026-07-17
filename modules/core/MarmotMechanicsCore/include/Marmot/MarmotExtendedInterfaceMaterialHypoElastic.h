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
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotStateHelpers.h"

#include <memory>
#include <string>
#include <vector>

/**
 * Abstract base class for hypoelastic interface materials.
 *
 * The interface kinematics use the actual connector between paired mesh
 * points,
 *
 *   d = x_top - x_bottom = ell n + d_tau,
 *
 * and reconstruct the displacement-gradient increment as
 *
 *   dG = dG_s + (1/ell) ( d[u] - dG_s d_tau ) \otimes n.
 *
 * For coincident interface faces, ell falls back to the constitutive
 * interface thickness h and d_tau is set to zero.
 */
class MarmotExtendedInterfaceMaterialHypoElastic {

protected:
  const double*                                materialProperties;
  const int                                    nMaterialProperties;
  double                                       h = 0.0;
  std::vector< double >                        baseMaterialProperties;
  std::unique_ptr< MarmotMaterialHypoElastic > baseMaterial;

public:
  using TensorMap3d    = Marmot::FastorStandardTensors::TensorMap3d;
  using TensorMap33d   = Marmot::FastorStandardTensors::TensorMap33d;
  using TensorMap333d  = Marmot::FastorStandardTensors::TensorMap333d;
  using TensorMap3333d = Marmot::FastorStandardTensors::TensorMap3333d;
  using TensorMap6d    = Marmot::FastorStandardTensors::TensorMap6d;
  using TensorMap18d   = Marmot::FastorStandardTensors::TensorMap18d;

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
    /** Generalized force conjugate to the raw mesh displacement jump. */
    TensorMap3d force;

    /**
     * Generalized surface resultant conjugate to the average surface
     * gradient. For a nonzero tangential connector,
     *
     *   surfaceStress = h sigma - force \otimes d_tau.
     */
    TensorMap33d surfaceStress;

    double* stateVars;
  };

  struct Tangents {
    /** d(force_i) / d([u_k]). */
    TensorMap33d Q_ij;

    /** d(surfaceStress_ij) / d(<u_{k,l}>_s). */
    TensorMap3333d Z_ijkl;

    /** d(force_i) / d(<u_{k,l}>_s). */
    TensorMap333d H_ijk;

    /** d(surfaceStress_ij) / d([u_k]). */
    TensorMap333d K_ijk;
  };

  struct Deformation {
    TensorMap6d  dU;
    TensorMap18d dSurfaceStrain;
    TensorMap3d  normal;
    TensorMap3d  separationVector;

    // Fastor's const TensorMap cannot be used with slicing and norm operations.
    // These views are therefore mutable types but are exposed through const Deformation&.
    Deformation( const double* dU_,
                 const double* dSurfaceStrain_,
                 const double* normal_,
                 const double* separationVector_ )
      : dU( const_cast< double* >( dU_ ) ),
        dSurfaceStrain( const_cast< double* >( dSurfaceStrain_ ) ),
        normal( const_cast< double* >( normal_ ) ),
        separationVector( const_cast< double* >( separationVector_ ) )
    {
    }

    /**
     * Backward-compatible constructor for coincident interface faces.
     * The material then uses d = 0, ell = h, and d_tau = 0.
     */
    Deformation( const double* dU_, const double* dSurfaceStrain_, const double* normal_ )
      : Deformation( dU_, dSurfaceStrain_, normal_, zeroSeparationVector() )
    {
    }

  private:
    static double* zeroSeparationVector()
    {
      static double zero[3] = { 0.0, 0.0, 0.0 };
      return zero;
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

  StateView getStateView( const std::string& stateName, double* stateVars ) const
  {
    return stateLayout.getStateView( stateVars, stateName );
  }

  virtual int getNumberOfRequiredStateVars() const { return stateLayout.totalSize(); }

  virtual void initializeYourself( double* stateVars, int nStateVars );

  virtual double getDensity();
};