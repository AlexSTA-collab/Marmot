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
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrain.h"

#include <memory>
#include <string>

/**
 * Equilibrated two-sided interface material at finite strain, which carries a nonlocal field.
 *
 * The gradient-enhanced counterpart of MarmotEquilibratedInterfaceMaterialFiniteStrain. Both faces are driven by a
 * MarmotMaterialGradientEnhancedFiniteStrain and share ONE value of the nonlocal field abar.
 *
 * Sharing it is the assumption l > h of a nonlocal radius that exceeds the thickness of the layer, and it is the whole
 * simplification. The nonlocal equation carries l^2 on the gradient of abar, so a variation through the thickness
 * costs about (l/h)^2 times what a variation in the plane costs. The field is then flat across the layer, and the
 * faces of the layer pass no flux into the adjacent bodies. What survives is a surface Helmholtz equation on the
 * midsurface, with one value per node of the midsurface: no jump of abar, and so no connector correction for it.
 *
 * Besides the generalized stress, the material returns the average of the local fields of the two faces,
 * Lbar = ( L+ + L- ) / 2, and the average of their nonlocal radii, and three more tangents: d(p)/d(abar),
 * d(Lbar)/d(x) and d(Lbar)/d(abar). All three are condensed through z exactly as the generalized stress is, because
 * z solves P+ N = P- N for a given abar and therefore depends on it.
 *
 * The material properties and the local solution are the ones of the base class.
 */
class MarmotEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain
  : public MarmotEquilibratedInterfaceMaterialFiniteStrainBase {

protected:
  std::unique_ptr< MarmotMaterialGradientEnhancedFiniteStrain > topMaterial;    ///< material of the top face
  std::unique_ptr< MarmotMaterialGradientEnhancedFiniteStrain > bottomMaterial; ///< material of the bottom face

  /**
   * @brief Solve for the normal-gradient jump with a nonlocal field of each face.
   *
   * @param[in,out] stateVars state variables of the interface material
   * @param[in] deformation kinematics of the interface
   * @param[in] timeIncrement time at the beginning of the increment and time increment
   * @param[in] nonlocalFieldTop nonlocal field of the top face
   * @param[in] nonlocalFieldBottom nonlocal field of the bottom face
   * @return the solution, with the derivatives with respect to the field of each face
   */
  Solution solveWithFields( double*                                                                 stateVars,
                            const MarmotEquilibratedInterfaceMaterialFiniteStrainBase::Deformation& deformation,
                            const TimeIncrement&                                                    timeIncrement,
                            double                                                                  nonlocalFieldTop,
                            double nonlocalFieldBottom ) const;

public:
  using Tensor21d = Fastor::Tensor< double, 21 >; ///< vector over the generalized strain (w, A+, A-)

  /// Generalized stress, and the local field of the interface.
  struct State : MarmotEquilibratedInterfaceMaterialFiniteStrainBase::State {
    double localField;     ///< average Lbar of the local driving fields of the two faces
    double nonlocalRadius; ///< average of the nonlocal radii of the two faces
  };

  /// Condensed tangent, and the tangents of the coupling to the nonlocal field.
  struct Tangents : MarmotEquilibratedInterfaceMaterialFiniteStrainBase::Tangents {
    Tensor21d dGeneralizedStress_dNonlocal;   ///< d(f, S+, S-)/d(abar)
    Tensor21d dLocalField_dGeneralizedStrain; ///< d(Lbar)/d(w, A+, A-)
    double    dLocalField_dNonlocal;          ///< d(Lbar)/d(abar)
  };

  /// Kinematics of the interface at a quadrature point, and the nonlocal field.
  struct Deformation : MarmotEquilibratedInterfaceMaterialFiniteStrainBase::Deformation {
    double nonlocalField; ///< nonlocal field abar, shared by both faces
  };

  /// Generalized stress, and the local fields of the two faces.
  struct FaceFieldsState : MarmotEquilibratedInterfaceMaterialFiniteStrainBase::State {
    double localFieldTop;        ///< local driving field L+ of the top face
    double localFieldBottom;     ///< local driving field L- of the bottom face
    double nonlocalRadiusTop;    ///< nonlocal radius of the top face
    double nonlocalRadiusBottom; ///< nonlocal radius of the bottom face
  };

  /// Condensed tangent, and the tangents of the coupling to the nonlocal fields of the two faces.
  struct FaceFieldsTangents : MarmotEquilibratedInterfaceMaterialFiniteStrainBase::Tangents {
    Tensor21d dGeneralizedStress_dNonlocalTop;      ///< d(f, S+, S-)/d(N+)
    Tensor21d dGeneralizedStress_dNonlocalBottom;   ///< d(f, S+, S-)/d(N-)
    Tensor21d dLocalFieldTop_dGeneralizedStrain;    ///< d(L+)/d(w, A+, A-)
    Tensor21d dLocalFieldBottom_dGeneralizedStrain; ///< d(L-)/d(w, A+, A-)
    double    dLocalFieldTop_dNonlocalTop;          ///< d(L+)/d(N+)
    double    dLocalFieldTop_dNonlocalBottom;       ///< d(L+)/d(N-), through the normal-gradient jump
    double    dLocalFieldBottom_dNonlocalTop;       ///< d(L-)/d(N+), through the normal-gradient jump
    double    dLocalFieldBottom_dNonlocalBottom;    ///< d(L-)/d(N-)
  };

  /// Kinematics of the interface at a quadrature point, and the nonlocal field of each face.
  struct FaceFieldsDeformation : MarmotEquilibratedInterfaceMaterialFiniteStrainBase::Deformation {
    double nonlocalFieldTop;    ///< nonlocal field N+ of the top face
    double nonlocalFieldBottom; ///< nonlocal field N- of the bottom face
  };

  /**
   * @brief Construct the material from a registered gradient-enhanced finite-strain material.
   *
   * @param[in] materialName name of the bulk material, registered in the
   * MarmotMaterialGradientEnhancedFiniteStrainFactory
   * @param[in] matProperties_ material properties [p0, p1, h, ...]
   * @param[in] nMaterialProperties_ number of material properties, at least three
   * @param[in] materialNumber_ number of the material
   *
   * @throws std::invalid_argument if fewer than three properties are given, if the thickness is not positive, or if
   * the bulk material is not registered
   */
  MarmotEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain( const std::string& materialName,
                                                                   const double*      matProperties_,
                                                                   int                nMaterialProperties_,
                                                                   int                materialNumber_ );

  /**
   * Set the characteristic element length at the considered quadrature point.
   * It is handed on to the bulk materials of both faces.
   *
   * @param[in] length characteristic length
   */
  void setCharacteristicElementLength( double length );

  /**
   * @brief Compute the generalized stress, the local field and their tangents from the total kinematics.
   *
   * @param[in,out] state the response, and the state variables, which are updated
   * @param[out] tangents the condensed tangents
   * @param[in] deformation the kinematics of the interface and the nonlocal field
   * @param[in] timeIncrement the time at the beginning of the increment and the time increment
   *
   * @throws std::invalid_argument if the normal vanishes or the faces are not separated along the normal
   */
  virtual void computeStress( State&               state,
                              Tangents&            tangents,
                              const Deformation&   deformation,
                              const TimeIncrement& timeIncrement );

  /**
   * @brief Compute the generalized stress, the local fields and their tangents for a nonlocal field of each face.
   *
   * The faces have independent nonlocal fields, and so independent damage. They interact through the traction
   * equilibrium of the normal-gradient jump only, so that the derivatives of the stress and of the local field of a
   * face with respect to the field of the other face are the ones through that jump.
   *
   * @param[in,out] state the response, and the state variables, which are updated
   * @param[out] tangents the condensed tangents
   * @param[in] deformation the kinematics of the interface and the nonlocal fields
   * @param[in] timeIncrement the time at the beginning of the increment and the time increment
   *
   * @throws std::invalid_argument if the normal vanishes or the faces are not separated along the normal
   */
  virtual void computeStress( FaceFieldsState&             state,
                              FaceFieldsTangents&          tangents,
                              const FaceFieldsDeformation& deformation,
                              const TimeIncrement&         timeIncrement );

  /**
   * @brief Initialize the state variables at a material point.
   *
   * @param[in,out] stateVars the state variables
   * @param[in] nStateVars number of state variables
   */
  virtual void initializeYourself( double* stateVars, int nStateVars );

  /**
   * @brief Get the density of the material.
   *
   * @param[in] stateVars the state variables of the interface material
   * @return the density of the bulk material
   */
  virtual double getDensity( const double* stateVars ) const;
};
