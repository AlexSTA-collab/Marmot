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
#include "Fastor/Fastor.h"
#include "Marmot/MarmotMaterialHypoElasticInterface.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotWiechertInterface.h"
#include <iostream>
#include <string>
#include <vector>

namespace Marmot::Materials {
  /**
   * \brief Implementation of a linear elastic interface material
   * for 3D stress states.
   *
   * For further information see \ref linearelasticinterface.
   * according to the LinearViscoelasticPowerLaw model by Bazant et al. (2015)

   * generalized for 3D stress states.

   *

   * For further information see \ref b4.

   */
  class LinearViscoElasticInterface : public MarmotMaterialHypoElasticInterface {

    /// \brief Young's modulus
    const double& E_M;

    /// \brief Poisson's ratio
    const double& nu_M;

    /// \brief Young's modulus
    const double& E_I;

    /// \brief Poisson's ratio
    const double& nu_I;

    /// \brief Young's modulus
    const double& E_0;

    /// \brief Poisson's ratio
    const double& nu_0;

    /// \brief height of the middle layer
    const double& h;

    /// \brief power law compliance parameter displacement jump
    const double& mRu;

    /// \brief power law exponent displacement jump
    const double& nRu;

    /// \brief number of Kelvin units to approximate the viscoelastic compliance for displacement jump
    const size_t nMaxwellRu;

    /// \brief minimal retardation time used in the viscoelastic Kelvin chain for displacement jump
    const double& minTauRu;

    /// \brief power law compliance parameter surface stress
    const double& mRs;

    /// \brief power law exponent surface stress
    const double& nRs;

    /// \brief number of Kelvin units to approximate the viscoelastic compliance for surface stress
    const size_t nMaxwellRs;

    /// \brief minimal retardation time used in the viscoelastic Kelvin chain for surface stress
    const double& minTauRs;

    /// \brief ratio of simulation time to days
    const double& timeToDays;

    class LinearViscoElasticInterfaceStateVarManager : public MarmotStateVarVectorManager {

    public:
      inline const static auto layout = makeLayout( {
        { .name = "MaxwellStateVarsRu", .length = 3 * 1 },
        { .name = "MaxwellStateVarsRs", .length = 9 * 1 },
      } );

      WiechertInterface::mapStateVarMatrixRu MaxwellStateVarsRu;
      WiechertInterface::mapStateVarMatrixRs MaxwellStateVarsRs;

      LinearViscoElasticInterfaceStateVarManager( double* theStateVarVector, int nMaxwellUnitsRu, int nMaxwellUnitsRs )
        : MarmotStateVarVectorManager( theStateVarVector, layout ),
          MaxwellStateVarsRu( &find( "MaxwellStateVarsRu" ), 3, nMaxwellUnitsRu ),
          MaxwellStateVarsRs( &find( "MaxwellStateVarsRs" ), 9, nMaxwellUnitsRs ){};
    };

    ::std::unique_ptr< LinearViscoElasticInterfaceStateVarManager > stateVarManager;

  public:
    using MarmotMaterialHypoElasticInterface::MarmotMaterialHypoElasticInterface;
    using Tensor1D = Fastor::Tensor< double, 3 >;
    using Tensor2D = Fastor::Tensor< double, 3, 3 >;

    LinearViscoElasticInterface( const double* materialProperties, int nMaterialProperties, int materialNumber );

    void computeStress( double*       force,
                        double*       surfaceStress,
                        double*       dStressDstrain,
                        const double* dU,
                        const double* dSurfaceStrain,
                        const double* normal,
                        const double* timeOld,
                        const double  dT,
                        double&       pNewDT );

    int getNumberOfRequiredStateVars();

    void assignStateVars( double* stateVars_, int nStateVars );

    StateView getStateView( const ::std::string& stateName );

  private:
    WiechertInterface::Properties elasticModuliRu;
    WiechertInterface::Properties elasticModuliRs;
    WiechertInterface::Properties relaxationTimesRu;
    WiechertInterface::Properties relaxationTimesRs;
    double                        zerothWiechertStiffnessRu;
    double                        zerothWiechertStiffnessRs;

    static constexpr int powerLawApproximationOrder = 2;
  };
} // namespace Marmot::Materials
