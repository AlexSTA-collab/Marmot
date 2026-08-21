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

#include "Marmot/ZInterfaceFiniteElement.h"

namespace Marmot::Elements {

  ZInterfaceFiniteElement::ZInterfaceFiniteElement( int                                         elementID,
                                                    FiniteElement::Quadrature::IntegrationTypes integrationType,
                                                    SectionType                                 sectionType )
    : ParentGeometryElement(),
      elementProperties( Eigen::Map< const Eigen::VectorXd >( nullptr, 0 ) ),
      elLabel( elementID ),
      sectionType( sectionType )
  {
    const auto qpInfos = FiniteElement::Quadrature::getGaussPointInfo( this->shape, integrationType );
    for ( const auto& qpInfo : qpInfos ) {
      QuadraturePoint qp( qpInfo.xi, qpInfo.weight );
      qps.push_back( std::move( qp ) );
    }
  }

  int ZInterfaceFiniteElement::getNumberOfRequiredStateVars()
  {
    return qps[0].getNumberOfRequiredStateVars() * qps.size();
  }

  std::vector< std::vector< std::string > > ZInterfaceFiniteElement::getNodeFields()
  {
    using namespace std;

    static vector< vector< string > > nodeFields;

    if ( nodeFields.empty() ) {
      for ( int i = 0; i < nNodes; i++ ) {
        nodeFields.push_back( vector< string >() );
        nodeFields[i].push_back( "displacement" );
        if ( i < nInterfaceNodes ) {
          // bottom nodes represent the midsurface and carry the gradient jump
          nodeFields[i].push_back( "normalGradientJump" );
        }
      }
    }

    return nodeFields;
  }

  std::vector< int > ZInterfaceFiniteElement::getDofIndicesPermutationPattern()
  {
    static std::vector< int > permutationPattern;

    if ( permutationPattern.empty() ) {
      permutationPattern.resize( sizeLoadVector );

      // Canonical (node-major, field-minor) layout:
      //   bottom node A (0..3): 6 slots [disp(3), g(3)]  -> base 6*A
      //   top node   (4+A)    : 3 slots [disp(3)]        -> base 24 + 3*A
      for ( int A = 0; A < nInterfaceNodes; A++ )
        for ( int c = 0; c < 3; c++ )
          permutationPattern[offD + 3 * A + c] = 6 * A + c; // u- (bottom)

      for ( int A = 0; A < nInterfaceNodes; A++ )
        for ( int c = 0; c < 3; c++ )
          permutationPattern[offD + nSideDofU + 3 * A + c] = 24 + 3 * A + c; // u+ (top)

      for ( int A = 0; A < nInterfaceNodes; A++ )
        for ( int k = 0; k < 3; k++ )
          permutationPattern[offG + 3 * A + k] = 6 * A + 3 + k; // g on bottom node A
    }

    return permutationPattern;
  }

  void ZInterfaceFiniteElement::assignStateVars( double* stateVars, int nStateVars )
  {
    const int nQpStateVars = nStateVars / qps.size();

    for ( size_t i = 0; i < qps.size(); i++ ) {
      auto&   qp          = qps[i];
      double* qpStateVars = stateVars + ( i * nQpStateVars );
      qp.assignStateVars( qpStateVars, nQpStateVars );
    }
  }

  void ZInterfaceFiniteElement::assignProperty( const ElementProperties& elementPropertiesInfo )
  {
    new ( &elementProperties ) Eigen::Map< const Eigen::VectorXd >( elementPropertiesInfo.elementProperties,
                                                                    elementPropertiesInfo.nElementProperties );
  }

  /**
   * Push the element-level settings into a quadrature point's material.
   *
   * Called from BOTH initializeYourself and the two material assignments,
   * because the host framework does not fix their order: EdelweissFE's SOLID
   * section calls initializeElement() and only then setMaterial(), and never
   * calls setProperties() at all, while its PLANE section calls
   * setProperties() first. Applying this in one place only would silently
   * leave the characteristic element length -- and the regularization zeta --
   * at their defaults in production runs.
   */
  void ZInterfaceFiniteElement::applyMaterialSettings( QuadraturePoint& qp )
  {
    if ( !qp.material || qp.sqrtDetG <= 0.0 )
      return;

    qp.material->setCharacteristicElementLength( std::sqrt( qp.sqrtDetG ) );

    // Element-property slot 1 overrides zeta. Note that a SOLID section never
    // supplies element properties, so there zeta comes from MARMOT_ZIFACE_REG
    // or from the material default.
    if ( elementProperties.size() > 1 )
      qp.material->setGradientJumpRegularization( elementProperties[1] );
  }

  void ZInterfaceFiniteElement::assignProperty( const MarmotMaterialSection& section )
  {
    for ( auto& qp : qps ) {
      qp.material = std::make_unique< Material >( section.materialName,
                                                  section.materialProperties,
                                                  section.nMaterialProperties,
                                                  elLabel );
      applyMaterialSettings( qp );
    }
  }

  void ZInterfaceFiniteElement::assignMaterial( const std::string& materialName,
                                                const double*      materialProperties,
                                                int                nMaterialProperties )
  {
    for ( auto& qp : qps ) {
      qp.material = std::make_unique< Material >( materialName, materialProperties, nMaterialProperties, elLabel );
      applyMaterialSettings( qp );
    }
  }

  void ZInterfaceFiniteElement::assignNodeCoordinates( const double* coordinates )
  {
    ParentGeometryElement::assignNodeCoordinates( coordinates );
  }

  void ZInterfaceFiniteElement::initializeYourself()
  {
    // Slot 0 is the YIQUAD4 out-of-plane thickness (1.0 in 3D). Slot 1, when
    // present, overrides the gradient-jump regularization zeta of the material.
    const double thickness = elementProperties.size() > 0 ? elementProperties[0] : 1.0;

    for ( QuadraturePoint& qp : qps ) {
      const auto geom = this->evaluateAt( qp.xi, 0 );

      qp.N                 = geom.N;
      qp.dNdXi             = geom.dNdXi;
      qp.J                 = geom.J;
      qp.G                 = geom.G;
      qp.sqrtDetG          = geom.sqrtDetG;
      qp.detJ              = geom.sqrtDetG;
      qp.gradN             = geom.gradN;
      qp.normal            = geom.n;
      qp.normalProjection  = geom.normalProjection;
      qp.tangentProjection = geom.tangentProjection;

      qp.NmatSide    = geom.NmatSide;
      qp.BmatSide    = geom.BmatSide;
      qp.NmatJump    = geom.NmatJump;
      qp.BmatAverage = geom.BmatAverage;

      const VectorDim xBottom = qp.NmatSide * this->getSideCoordinates( 0 );
      const VectorDim xTop    = qp.NmatSide * this->getSideCoordinates( 1 );
      qp.separationVector     = xTop - xBottom;

      constexpr double geometryTolerance = 1.0e-12;
      if ( qp.separationVector.norm() > geometryTolerance && qp.separationVector.dot( qp.normal ) < 0.0 ) {
        qp.normal *= -1.0;
      }

      qp.normalProjection     = qp.normal * qp.normal.transpose();
      qp.tangentProjection    = TensorDim::Identity() - qp.normalProjection;
      qp.normalSeparation     = qp.separationVector.dot( qp.normal );
      qp.tangentialSeparation = qp.tangentProjection * qp.separationVector;

      if ( qp.separationVector.norm() > geometryTolerance && qp.normalSeparation <= geometryTolerance ) {
        throw std::invalid_argument(
          "ZInterfaceFiniteElement: paired faces have no positive separation in the interface-normal direction." );
      }

      qp.J0xW = qp.weight * qp.sqrtDetG * thickness;

      // B_X = [ J ; B+ ; B- ], the geometric map from the element displacement
      // DOFs to the generalised strain X_M of eq. (43). Bottom nodes carry u-,
      // top nodes u+, matching the YIQUAD4 ordering.
      qp.BX.setZero();
      qp.BX.block< 3, nDofD >( 0, 0 )                     = qp.NmatJump;
      qp.BX.block< nTensor, nSideDofU >( 3, nSideDofU )   = qp.BmatSide; // A+ from the top nodes
      qp.BX.block< nTensor, nSideDofU >( 3 + nTensor, 0 ) = qp.BmatSide; // A- from the bottom nodes

      // N_g: block-diagonal Q1 interpolation of a 3-vector over the 4 midsurface nodes.
      qp.Ng.setZero();
      for ( int A = 0; A < nInterfaceNodes; A++ )
        for ( int k = 0; k < 3; k++ )
          qp.Ng( k, 3 * A + k ) = qp.N( A );

      applyMaterialSettings( qp );
    }
  }

  void ZInterfaceFiniteElement::computeKernels( const double* QTotal_,
                                                const double* dQ_,
                                                double*       Pe_,
                                                double*       Ke_,
                                                double        time,
                                                double        dT )
  {
    Eigen::Map< const RhsSized > QTotal( QTotal_ );
    Eigen::Map< const RhsSized > dQ( dQ_ );
    Eigen::Map< KeSizedMatrix >  Ke( Ke_ );
    Eigen::Map< RhsSized >       Pe( Pe_ );

    Ke.setZero();
    Pe.setZero();

    for ( QuadraturePoint& qp : qps ) {
      const auto& BX = qp.BX; // 21x24
      const auto& Ng = qp.Ng; // 3x12

      const Eigen::Matrix< double, nDofD, 1 > dd = dQ.segment< nDofD >( offD );
      const Eigen::Matrix< double, nDofG, 1 > dg = dQ.segment< nDofG >( offG );

      // The material takes the two faces separately, as eq. (43) prescribes.
      // dU is ordered [ du+ ; du- ]; X_M is rebuilt inside the material from it.
      const Eigen::Matrix< double, nX, 1 > dX = BX * dd;
      const Eigen::Matrix< double, nZ, 1 > dz = Ng * dg;

      Eigen::Matrix< double, 6, 1 > dU;
      dU.segment< 3 >( 0 ) = qp.NmatSide * dd.segment< nSideDofU >( nSideDofU ); // u+ (top)
      dU.segment< 3 >( 3 ) = qp.NmatSide * dd.segment< nSideDofU >( 0 );         // u- (bottom)

      Eigen::Matrix< double, 18, 1 > dSurfaceStrain;
      dSurfaceStrain.segment< 9 >( 0 ) = dX.segment< 9 >( 3 );  // A+
      dSurfaceStrain.segment< 9 >( 9 ) = dX.segment< 9 >( 12 ); // A-

      Eigen::Matrix< double, 3, 1 > generalizedForce;
      Eigen::Matrix< double, 9, 1 > surfaceStressPlus, surfaceStressMinus;
      Eigen::Matrix< double, 3, 1 > tractionImbalance;

      Eigen::Matrix< double, nX, nX, Eigen::RowMajor > K_xx;
      Eigen::Matrix< double, nX, nZ, Eigen::RowMajor > K_xz;
      Eigen::Matrix< double, nZ, nX, Eigen::RowMajor > K_zx;
      Eigen::Matrix< double, nZ, nZ, Eigen::RowMajor > K_zz;

      Material::State         state{ generalizedForce.data(),
                             surfaceStressPlus.data(),
                             surfaceStressMinus.data(),
                             tractionImbalance.data(),
                             qp.managedStateVars->materialStateVars.data() };
      Material::Tangents      tangents{ K_xx.data(), K_xz.data(), K_zx.data(), K_zz.data() };
      Material::Deformation   deformation{ dU.data(),
                                         dSurfaceStrain.data(),
                                         dz.data(),
                                         qp.normal.data(),
                                         qp.separationVector.data() };
      Material::TimeIncrement timeIncrement{ time, dT };

      qp.material->computeStress( state, tangents, deformation, timeIncrement );

      qp.managedStateVars->generalizedForce   = generalizedForce;
      qp.managedStateVars->surfaceStressPlus  = surfaceStressPlus;
      qp.managedStateVars->surfaceStressMinus = surfaceStressMinus;
      qp.managedStateVars->tractionImbalance  = tractionImbalance;
      // Report the TOTAL gradient jump at the point, read straight off the
      // nodal field, rather than accumulating increments.
      qp.managedStateVars->normalGradientJump = Ng * QTotal.segment< nDofG >( offG );

      // p_X = ( f , S+ , S- ), work-conjugate to X_M.
      Eigen::Matrix< double, nX, 1 > pX;
      pX.segment< 3 >( 0 )  = generalizedForce;
      pX.segment< 9 >( 3 )  = surfaceStressPlus;
      pX.segment< 9 >( 12 ) = surfaceStressMinus;

      // Residuals. Pe holds MINUS the residual, as everywhere in this family.
      Pe.segment< nDofD >( offD ) -= BX.transpose() * pX * qp.J0xW;
      Pe.segment< nDofG >( offG ) -= Ng.transpose() * tractionImbalance * qp.J0xW;

      Ke.block< nDofD, nDofD >( offD, offD ) += BX.transpose() * K_xx * BX * qp.J0xW;
      Ke.block< nDofD, nDofG >( offD, offG ) += BX.transpose() * K_xz * Ng * qp.J0xW;
      Ke.block< nDofG, nDofD >( offG, offD ) += Ng.transpose() * K_zx * BX * qp.J0xW;
      Ke.block< nDofG, nDofG >( offG, offG ) += Ng.transpose() * K_zz * Ng * qp.J0xW;
    }
  }

  void ZInterfaceFiniteElement::setInitialConditions( StateTypes state, const double* values )
  {
    switch ( state ) {
    case MarmotElement::MarmotMaterialInitialization: {
      for ( QuadraturePoint& qp : qps ) {
        qp.material->initializeYourself( qp.managedStateVars->materialStateVars.data(),
                                         qp.managedStateVars->materialStateVars.size() );
      }
      break;
    }

    case MarmotElement::MarmotMaterialStateVars: {
      throw std::invalid_argument( "Please use initializeStateVars directly on material" );
    }

    default:
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__
                                                << ": invalid initial condition for ZInterfaceFiniteElement" );
    }
  }

  void ZInterfaceFiniteElement::computeDistributedLoad( MarmotElement::DistributedLoadTypes loadType,
                                                        double*                             P,
                                                        double*                             K,
                                                        const int                           elementFace,
                                                        const double*                       load,
                                                        const double*                       QTotal,
                                                        double                              time,
                                                        double                              dT )
  {
    throw std::invalid_argument(
      MakeString() << __PRETTY_FUNCTION__ << ": distributed loads are not implemented for ZInterfaceFiniteElement." );
  }

  void ZInterfaceFiniteElement::computeBodyForce( double*       P,
                                                  double*       K,
                                                  const double* load,
                                                  const double* QTotal,
                                                  double        time,
                                                  double        dT )
  {
    throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__
                                              << ": body forces are not implemented for ZInterfaceFiniteElement." );
  }

  void ZInterfaceFiniteElement::computeConsistentInertia( double* M )
  {
    throw std::runtime_error( MakeString()
                              << __PRETTY_FUNCTION__ << ": inertia is not implemented for ZInterfaceFiniteElement." );
  }

  void ZInterfaceFiniteElement::computeLumpedInertia( double* M )
  {
    throw std::runtime_error( MakeString()
                              << __PRETTY_FUNCTION__ << ": inertia is not implemented for ZInterfaceFiniteElement." );
  }

  std::vector< double > ZInterfaceFiniteElement::getCoordinatesAtCenter()
  {
    std::vector< double >   coords( nDim );
    Eigen::Map< VectorDim > coordsMap( coords.data() );
    const auto              centerXi = XiSized::Zero();
    const auto              Ncenter  = this->N( centerXi );
    const auto              Nmat     = this->NMatrix( Ncenter );
    coordsMap                        = Nmat * this->getSideCoordinates( 0 );
    return coords;
  }

  std::vector< std::vector< double > > ZInterfaceFiniteElement::getCoordinatesAtQuadraturePoints()
  {
    std::vector< std::vector< double > > listedCoords;
    for ( const auto& qp : qps ) {
      std::vector< double >   coords( nDim );
      Eigen::Map< VectorDim > coordsMap( coords.data() );
      coordsMap = qp.NmatSide * this->getSideCoordinates( 0 );
      listedCoords.push_back( coords );
    }
    return listedCoords;
  }

  int ZInterfaceFiniteElement::getNumberOfQuadraturePoints()
  {
    return static_cast< int >( qps.size() );
  }

} // namespace Marmot::Elements
