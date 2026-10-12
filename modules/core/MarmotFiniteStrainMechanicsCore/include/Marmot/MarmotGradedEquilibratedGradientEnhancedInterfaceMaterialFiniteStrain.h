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

#include "Marmot/MarmotEquilibratedInterfaceMaterialFiniteStrainBase.h"
#include "Marmot/MarmotGradedFiniteStrainLayer.h"
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrain.h"
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrainFactory.h"

#include <Eigen/Dense>

#include <string>

/**
 * @brief Two-sided gradient-enhanced finite-strain interface material with a nonlocal field at every station through
 * the thickness of the layer.
 *
 * This material is the gradient-enhanced counterpart of MarmotGradedEquilibratedInterfaceMaterialFiniteStrain and
 * the two-sided counterpart of MarmotGradedGradientEnhancedInterfaceMaterialFiniteStrain. The surface gradients
 * @f$A^+@f$ and @f$A^-@f$ of the faces are kept separate and interpolated linearly through the thickness, the
 * microproblem and the closure of the tangential gradient of the correction are those of the first one, and the bulk
 * material of every station is a MarmotMaterialGradientEnhancedFiniteStrain with a nonlocal field @f$\bar N_q@f$ of
 * its own, so that the damage, which the nonlocal field drives, varies through the thickness. The stations are the
 * integration points of the thickness, with the weights of the trapezoidal rule (the mode 2 of
 * Marmot::Materials::ThicknessProfile). For two stations at the faces and without the closure of the tangential
 * gradient, the material is the one of MarmotEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain with a
 * nonlocal field for each face.
 *
 * The microproblem is solved for given nonlocal fields. Since the fields enter the stress of their point, the
 * corrections depend on them, and so do the generalized stress and the local fields. In addition to the response of
 * the interface the material returns the local driving field and the nonlocal radius of every station, and the
 * derivatives of the generalized stress and of the local fields with respect to the generalized strain
 * @f$([\![u]\!], A^+, A^-)@f$ and to the nonlocal fields, all condensed over the corrections.
 *
 * The material properties are the ones of Marmot::Materials::ThicknessProfile with the mode 2: @f$[h, 2, 1, K, m,
 * (p_k, \text{properties}) \times K]@f$, optionally followed by the flag of the tangential gradient, where the
 * properties are the ones of the bulk material.
 */
class MarmotGradedEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain
  : public MarmotGradedFiniteStrainLayer< MarmotMaterialGradientEnhancedFiniteStrain,
                                          MarmotLibrary::MarmotMaterialGradientEnhancedFiniteStrainFactory > {

public:
  using Tensor3d    = MarmotEquilibratedInterfaceMaterialFiniteStrainBase::Tensor3d;    ///< first order tensor in 3D
  using Tensor33d   = MarmotEquilibratedInterfaceMaterialFiniteStrainBase::Tensor33d;   ///< second order tensor in 3D
  using Tensor333d  = MarmotEquilibratedInterfaceMaterialFiniteStrainBase::Tensor333d;  ///< third order tensor in 3D
  using Tensor3333d = MarmotEquilibratedInterfaceMaterialFiniteStrainBase::Tensor3333d; ///< fourth order tensor in 3D
  using Tensor6d    = MarmotEquilibratedInterfaceMaterialFiniteStrainBase::Tensor6d;    ///< vector of six components
  using Tensor18d   = MarmotEquilibratedInterfaceMaterialFiniteStrainBase::Tensor18d;   ///< vector of 18 components

  using TimeIncrement = MarmotEquilibratedInterfaceMaterialFiniteStrainBase::TimeIncrement; ///< time data

  /// Response: the generalized stress, the state variables and the fields of the stations.
  struct State {
    Tensor3d        generalizedForce;   ///< conjugate to the jump w
    Tensor33d       surfaceStressPlus;  ///< conjugate to A+
    Tensor33d       surfaceStressMinus; ///< conjugate to A-
    double*         stateVars;          ///< pointer to the state variables
    Eigen::VectorXd localFields;        ///< local driving field of every station
    Eigen::VectorXd nonlocalRadii;      ///< nonlocal radius of every station
  };

  /// Condensed block tangent over ( w, A+, A- ) and the couplings to the nonlocal fields.
  struct Tangents {
    Tensor33d   Q_ww;   ///< d(f)/d(w), 3x3
    Tensor333d  Q_wAp;  ///< d(f)/d(A+), 3x9
    Tensor333d  Q_wAm;  ///< d(f)/d(A-), 3x9
    Tensor333d  Q_Apw;  ///< d(S+)/d(w), 9x3
    Tensor3333d Q_ApAp; ///< d(S+)/d(A+), 9x9
    Tensor3333d Q_ApAm; ///< d(S+)/d(A-), 9x9
    Tensor333d  Q_Amw;  ///< d(S-)/d(w), 9x3
    Tensor3333d Q_AmAp; ///< d(S-)/d(A+), 9x9
    Tensor3333d Q_AmAm; ///< d(S-)/d(A-), 9x9

    /// derivative of ( f, S+, S- ) with respect to the nonlocal field of every station, 21 x K
    Eigen::MatrixXd dStress_dNonlocal;
    /// derivative of the local field of every station with respect to ( w, A+, A- ), K x 21
    Eigen::MatrixXd dLocal_dStrain;
    /// derivative of the local field of every station with respect to the nonlocal field of every station, K x K
    Eigen::MatrixXd dLocal_dNonlocal;
  };

  /// Kinematics of the interface at a quadrature point, and the nonlocal field of every station.
  struct Deformation {
    Tensor6d        U;               ///< total displacements of the two interface sides, top first
    Tensor18d       surfaceGradient; ///< total surface gradients A+ and A- of the two interface sides
    Tensor3d        normal;          ///< interface normal, pointing from the first to the second side
    Tensor3d        connector;       ///< reference connector between the paired faces, X+ - X-
    Eigen::VectorXd nonlocalFields;  ///< nonlocal field of every station, from the bottom to the top
  };

  /**
   * @brief Construct the material from a registered gradient-enhanced finite-strain bulk material.
   *
   * @param[in] materialName name of the bulk material, registered in the
   * MarmotMaterialGradientEnhancedFiniteStrainFactory
   * @param[in] matProperties material properties, as described for this class
   * @param[in] nMaterialProperties number of material properties
   * @param[in] materialNumber_ number of the material
   *
   * @throws std::invalid_argument if the properties are inconsistent, if the mode is not 2, or if the bulk material is
   * not registered
   */
  MarmotGradedEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain( const std::string& materialName,
                                                                         const double*      matProperties,
                                                                         int                nMaterialProperties,
                                                                         int                materialNumber_ );

  /**
   * @brief Solve the microproblem of the layer and compute the generalized stress and the nonlocal couplings.
   *
   * @param[in,out] state generalized stress, fields of the stations and state variables
   * @param[out] tangents condensed tangent and couplings to the nonlocal fields
   * @param[in] deformation total displacements and surface gradients of both faces, the normal, the connector and
   * the nonlocal fields of the stations
   * @param[in] timeIncrement time at the beginning of the increment and time increment
   *
   * @throws std::invalid_argument if the number of nonlocal fields is not the number of stations, or the normal or
   * the connector is invalid
   * @throws Marmot::StressUpdateFailed if the iteration for the microproblem does not converge, or if the acoustic
   * tensor of a station is singular
   */
  virtual void computeStress( State&               state,
                              Tangents&            tangents,
                              const Deformation&   deformation,
                              const TimeIncrement& timeIncrement );

  /**
   * @brief Number of stations through the thickness, which is the number of nonlocal fields.
   *
   * @return the number of stations
   */
  int numberOfStations() const { return numberOfThicknessPoints(); }
};
