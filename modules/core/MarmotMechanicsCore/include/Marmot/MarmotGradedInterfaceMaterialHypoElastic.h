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
#include "Marmot/MarmotInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotStateHelpers.h"

#include <memory>
#include <string>
#include <vector>

/**
 * @brief Hypoelastic interface material with a stiffness that varies through the thickness of the layer.
 *
 * The layer has the thickness @f$h@f$ and the through-thickness coordinate @f$y \in [0,1]@f$, which starts at the
 * second side of the interface (the bottom) and points along the normal @f$n@f$ to the first side (the top). The
 * bulk material at the height @f$y@f$ is a hypoelastic material with the stiffness @f$C(y)@f$.
 *
 * The interface kinematics are those of MarmotInterfaceMaterialHypoElastic, namely
 * @f[
 *   \varepsilon^0 = \mathrm{sym} \left( \frac{1}{h} [\![u]\!] \otimes n + \overline{\nabla_s u} \right) ,
 * @f]
 * with the displacement jump @f$[\![u]\!]@f$ and the average surface gradient @f$\overline{\nabla_s u}@f$ of the
 * two sides. The strain in the layer is this strain plus the correction of the microproblem,
 * @f[
 *   \varepsilon(y) = \varepsilon^0 + \mathrm{sym} \left( w(y) \otimes n \right) , \qquad w = \partial_y u^1 ,
 * @f]
 * where @f$u^1@f$ is the first-order displacement of the layer. The microproblem follows from the variational
 * form with a vanishing macroscopic test function. It requires that the traction
 * @f[
 *   t = \sigma(y) \, n
 * @f]
 * is constant through the thickness and that @f$u^1@f$ vanishes at both sides, hence that
 * @f$\int_0^1 w \, \mathrm{d}y = 0@f$. The terms of @f$u^1@f$ with the tangential gradient, which are of the
 * order @f$h@f$ higher, are neglected. For a stiffness that varies through the thickness the strain is therefore not
 * homogeneous.
 *
 * The microproblem is solved at every call of computeStress() by a Newton iteration for the values of @f$w@f$ at
 * the integration points through the thickness. With the acoustic tensor @f$Q_q = n \cdot C_q \cdot n@f$ of the
 * integration point @f$q@f$ with the weight @f$\omega_q@f$, a Newton step is given by
 * @f[
 *   t = \Big( \sum_q \omega_q Q_q^{-1} \Big)^{-1} \Big( \sum_q \omega_q Q_q^{-1} a_q - \sum_q \omega_q w_q \Big) ,
 *   \qquad
 *   \delta w_q = Q_q^{-1} ( t - a_q ) ,
 * @f]
 * where @f$a_q = \sigma_q n@f$ is the traction of the integration point. A linear material converges in one step.
 * After the convergence the tangent of the homogenized layer
 * @f[
 *   \bar{C} = \sum_q \omega_q C_q \left( I + \frac{\partial \, \mathrm{sym}(w_q \otimes n)}{\partial \varepsilon^0}
 *   \right)
 * @f]
 * is condensed and converted to the tangent terms of the interface in the same way as in
 * MarmotInterfaceMaterialHypoElastic, with the force @f$\langle \sigma \rangle n@f$ and the surface stress
 * @f$h \langle \sigma \rangle@f$ of the through-thickness mean of the stress. For a single homogeneous layer, the
 * response is identical to the one of MarmotInterfaceMaterialHypoElastic.
 *
 * The material properties are
 *
 * | Index | Name | Meaning |
 * |-------|------|---------|
 * | 0 | h | thickness of the layer |
 * | 1 | mode | 0: layers with constant properties, 1: properties interpolated linearly between stations |
 * | 2 | nGauss | number of Gauss points per layer (mode 0) or per interval between two stations (mode 1) |
 * | 3 | K | number of layers (mode 0) or of stations (mode 1, at least two) |
 * | 4 | m | number of properties of the bulk material |
 *
 * followed by K groups of the entry @f$p_k@f$ and the m properties of the bulk material. In the mode 0, @f$p_k@f$
 * is the relative thickness of the layer k, which are ordered from the bottom to the top and normalized by their
 * sum. In the mode 1, @f$p_k@f$ is the position @f$y@f$ of the station k, with @f$p_0 = 0@f$, @f$p_{K-1} = 1@f$ and
 * increasing positions in between.
 *
 * The state variables are, for every integration point through the thickness, the stress in Voigt notation followed
 * by the state variables of its bulk material.
 *
 * @note Like the other interface materials, this class is the material of the interface finite elements (see
 * Marmot::Elements::DisplacementJumpFormulation), and the name given to the element is the name of the bulk material.
 */
class MarmotGradedInterfaceMaterialHypoElastic {

public:
  using Tensor3d    = MarmotInterfaceMaterialHypoElastic::Tensor3d;    ///< first order tensor in 3D
  using Tensor33d   = MarmotInterfaceMaterialHypoElastic::Tensor33d;   ///< second order tensor in 3D
  using Tensor333d  = MarmotInterfaceMaterialHypoElastic::Tensor333d;  ///< third order tensor in 3D
  using Tensor3333d = MarmotInterfaceMaterialHypoElastic::Tensor3333d; ///< fourth order tensor in 3D
  using Tensor6d    = MarmotInterfaceMaterialHypoElastic::Tensor6d;    ///< vector of six components, one per side
  using Tensor18d   = MarmotInterfaceMaterialHypoElastic::Tensor18d;   ///< vector of 18 components, one per side

  using State         = MarmotInterfaceMaterialHypoElastic::State;         ///< interface response, as in the base
  using Tangents      = MarmotInterfaceMaterialHypoElastic::Tangents;      ///< tangent terms, as in the base
  using Deformation   = MarmotInterfaceMaterialHypoElastic::Deformation;   ///< interface kinematics, as in the base
  using TimeIncrement = MarmotInterfaceMaterialHypoElastic::TimeIncrement; ///< time data, as in the base

  /// number of the material
  const int materialNumber;

  /// Layout of the state variables
  MarmotStateLayoutDynamic stateLayout;

  /// Characteristic element length
  double characteristicElementLength = 0.0;

  /**
   * @brief Construct the graded interface material from a registered hypoelastic bulk material.
   *
   * @param[in] materialName name of the bulk material, registered in the MarmotMaterialHypoElasticFactory
   * @param[in] matProperties material properties, as in the table of this class
   * @param[in] nMaterialProperties number of material properties
   * @param[in] materialNumber_ number of the material
   *
   * @throws std::invalid_argument if the properties are inconsistent, for example a non-positive thickness or an
   * unknown mode, if the positions of the stations are not increasing from 0 to 1, or if the bulk material is not
   * registered
   */
  MarmotGradedInterfaceMaterialHypoElastic( const std::string& materialName,
                                            const double*      matProperties,
                                            int                nMaterialProperties,
                                            int                materialNumber_ );

  /// Default destructor
  virtual ~MarmotGradedInterfaceMaterialHypoElastic() = default;

  /**
   * @brief Set the characteristic element length at the considered quadrature point.
   *
   * @param[in] length characteristic length, handed to all bulk materials
   */
  void setCharacteristicElementLength( double length );

  /**
   * @brief Solve the microproblem of the layer and compute the interface response.
   *
   * @param[in,out] state previous interface response and state variables on entry, updated ones on exit
   * @param[out] tangents tangent terms of the interface
   * @param[in] deformation increments of the jump and of the surface gradients, and the interface normal
   * @param[in] timeIncrement time at the beginning of the increment and time increment
   *
   * @throws Marmot::StressUpdateFailed if the iteration for the microproblem does not converge, or if the acoustic
   * tensor of an integration point is singular
   */
  virtual void computeStress( State&               state,
                              Tangents&            tangents,
                              const Deformation&   deformation,
                              const TimeIncrement& timeIncrement );

  /**
   * @brief Get a view to the state variables.
   *
   * @param[in] stateName name of the state variable
   * @param[in] stateVars pointer to the state variable array
   * @return view to the requested state variable
   */
  StateView getStateView( const std::string& stateName, double* stateVars ) const
  {
    return stateLayout.getStateView( stateVars, stateName );
  }

  /**
   * @brief Get the total number of required state variables.
   *
   * @return number of state variables of all integration points through the thickness
   */
  virtual int getNumberOfRequiredStateVars() const { return stateLayout.totalSize(); }

  /**
   * @brief Initialize the state variables: zero stress and initial state of every bulk material.
   *
   * @param[in,out] stateVars state variables of the material point
   * @param[in] nStateVars number of state variables, which this material does not need
   */
  virtual void initializeYourself( double* stateVars, int nStateVars );

  /**
   * @brief Get the density of the material.
   *
   * @return through-thickness mean of the density of the bulk materials
   */
  virtual double getDensity();

  /**
   * @brief Number of integration points through the thickness.
   *
   * @return the number of integration points
   */
  int numberOfThicknessPoints() const { return static_cast< int >( points.size() ); }

private:
  /// Integration point through the thickness.
  struct ThicknessPoint {
    double                                       weight;     ///< weight, the weights of all points sum up to one
    double                                       position;   ///< coordinate y in [0, 1]
    std::vector< double >                        properties; ///< properties of the bulk material, which it refers to
    std::unique_ptr< MarmotMaterialHypoElastic > material;   ///< bulk material at the position
    std::size_t                                  offset; ///< offset of the state of this point in the state variables
    std::size_t                                  nBaseStateVars; ///< number of state variables of the bulk material
  };

  const double                  h;      ///< thickness of the layer
  std::vector< ThicknessPoint > points; ///< integration points through the thickness, from the bottom to the top
};
