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
#include <iostream>
#include <string>
#include <vector>
#include "Fastor/Fastor.h"
#include "Marmot/MarmotKelvinChainInterface.h"
#include "Marmot/MarmotMaterialHypoElasticInterface.h"
#include "Marmot/MarmotStateVarVectorManager.h"

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
    const double& m_Ju;

    /// \brief power law exponent displacement jump
    const double& n_Ju;

    /// \brief number of Kelvin units to approximate the viscoelastic compliance for displacement jump
    const size_t nKelvin_Ju;

    /// \brief minimal retardation time used in the viscoelastic Kelvin chain for displacement jump
    const double& minTau_Ju;

    /// \brief power law compliance parameter surface stress
    const double& m_Js;

    /// \brief power law exponent surface stress
    const double& n_Js;

    /// \brief number of Kelvin units to approximate the viscoelastic compliance for surface stress
    const size_t nKelvin_Js;

    /// \brief minimal retardation time used in the viscoelastic Kelvin chain for surface stress
    const double& minTau_Js;

    /// \brief ratio of simulation time to days
    const double& timeToDays;

    class LinearViscoElasticInterfaceStateVarManager : public MarmotStateVarVectorManager {

    public:
      inline const static auto layout = makeLayout( {
        { .name = "kelvinStateVars_Ju", .length = 0 },
        { .name = "kelvinStateVars_Js", .length = 0 },
      } );

      KelvinChainInterface::mapStateVarMatrix_Ju kelvinStateVars_Ju;
      KelvinChainInterface::mapStateVarMatrix_Js kelvinStateVars_Js;
      
      LinearViscoElasticInterfaceStateVarManager( double* theStateVarVector, int nKelvinUnits_Ju, int nKelvinUnits_Js )
        : MarmotStateVarVectorManager( theStateVarVector, layout ),
          kelvinStateVars_Ju( &find( "kelvinStateVars_Ju" ), 3, nKelvinUnits_Ju ), 
          kelvinStateVars_Js( &find( "kelvinStateVars_Js" ), 9, nKelvinUnits_Js )
          {};
    };
    ::std::unique_ptr< LinearViscoElasticInterfaceStateVarManager > stateVarManager;

  public:
    using MarmotMaterialHypoElasticInterface::MarmotMaterialHypoElasticInterface;
    using Tensor1D = Fastor::Tensor<double,3>;
    using Tensor2D = Fastor::Tensor<double,3,3>;

    LinearViscoElasticInterface( const double* materialProperties, int nMaterialProperties, int materialNumber );
    
    void computeStress( double*  force,
                        double*  surface_stress,
                        double* dStress_dStrain,
                        const double* dU,
                        const double* dSurface_strain,
                        const double* normal,
                        const double* timeOld,
                        const double  dT,
                        double&       pNewDT);

    int getNumberOfRequiredStateVars();

    void assignStateVars( double* stateVars_, int nStateVars );

    StateView getStateView( const ::std::string& stateName );

  private:
    KelvinChainInterface::Properties elasticModuli_Ju;
    KelvinChainInterface::Properties elasticModuli_Js;
    KelvinChainInterface::Properties retardationTimes_Ju;
    KelvinChainInterface::Properties retardationTimes_Js;
    double                  zerothKelvinChainCompliance_Ju;
    double                  zerothKelvinChainCompliance_Js;

    static constexpr int powerLawApproximationOrder = 2;
  };
} // namespace Marmot::Materials
