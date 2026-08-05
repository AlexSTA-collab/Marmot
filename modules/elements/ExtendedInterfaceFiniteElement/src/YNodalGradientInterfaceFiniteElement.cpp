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

#include "Marmot/YNodalGradientInterfaceFiniteElement.h"

namespace Marmot::Elements {

  YNodalGradientInterfaceFiniteElement::YNodalGradientInterfaceFiniteElement(
    int                                         elementID,
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

  int YNodalGradientInterfaceFiniteElement::getNumberOfRequiredStateVars()
  {
    return qps[0].getNumberOfRequiredStateVars() * qps.size();
  }

  std::vector< std::vector< std::string > > YNodalGradientInterfaceFiniteElement::getNodeFields()
  {
    using namespace std;

    static vector< vector< string > > nodeFields;

    if ( nodeFields.empty() ) {
      for ( int i = 0; i < nNodes; i++ ) {
        nodeFields.push_back( vector< string >() );
        nodeFields[i].push_back( "displacement" );
        if ( i < nInterfaceNodes ) {
          // bottom nodes represent the midsurface: carry g and t
          nodeFields[i].push_back( "normalGradientAverage" );
          nodeFields[i].push_back( "commonTraction" );
        }
      }
    }

    return nodeFields;
  }

  std::vector< int > YNodalGradientInterfaceFiniteElement::getDofIndicesPermutationPattern()
  {
    static std::vector< int > permutationPattern;

    if ( permutationPattern.empty() ) {
      permutationPattern.resize( sizeLoadVector );

      // Canonical (node-major, field-minor) layout:
      //   bottom node A (0..3): 9 slots [disp(3), g(3), t(3)]  -> base 9*A
      //   top node   (4+A)    : 3 slots [disp(3)]              -> base 36 + 3*A
      for ( int A = 0; A < nInterfaceNodes; A++ )
        for ( int c = 0; c < 3; c++ )
          permutationPattern[offD + 3 * A + c] = 9 * A + c; // u- (bottom)

      for ( int A = 0; A < nInterfaceNodes; A++ )
        for ( int c = 0; c < 3; c++ )
          permutationPattern[offD + nSideDofU + 3 * A + c] = 36 + 3 * A + c; // u+ (top)

      for ( int A = 0; A < nInterfaceNodes; A++ )
        for ( int k = 0; k < 3; k++ )
          permutationPattern[offG + 3 * A + k] = 9 * A + 3 + k; // g on bottom node A

      for ( int A = 0; A < nInterfaceNodes; A++ )
        for ( int k = 0; k < 3; k++ )
          permutationPattern[offT + 3 * A + k] = 9 * A + 6 + k; // t on bottom node A
    }

    return permutationPattern;
  }

  void YNodalGradientInterfaceFiniteElement::assignStateVars( double* stateVars, int nStateVars )
  {
    const int nQpStateVars = nStateVars / qps.size();

    for ( size_t i = 0; i < qps.size(); i++ ) {
      auto&   qp          = qps[i];
      double* qpStateVars = stateVars + ( i * nQpStateVars );
      qp.assignStateVars( qpStateVars, nQpStateVars );
    }
  }

  void YNodalGradientInterfaceFiniteElement::assignProperty( const ElementProperties& elementPropertiesInfo )
  {
    new ( &elementProperties ) Eigen::Map< const Eigen::VectorXd >( elementPropertiesInfo.elementProperties,
                                                                    elementPropertiesInfo.nElementProperties );
  }

  void YNodalGradientInterfaceFiniteElement::assignProperty( const MarmotMaterialSection& section )
  {
    for ( auto& qp : qps ) {
      qp.material = std::make_unique< Material >( section.materialName,
                                                  section.materialProperties,
                                                  section.nMaterialProperties,
                                                  elLabel );
    }
  }

  void YNodalGradientInterfaceFiniteElement::assignMaterial( const std::string& materialName,
                                                             const double*      materialProperties,
                                                             int                nMaterialProperties )
  {
    for ( auto& qp : qps ) {
      qp.material = std::make_unique< Material >( materialName, materialProperties, nMaterialProperties, elLabel );
    }
  }

  void YNodalGradientInterfaceFiniteElement::assignNodeCoordinates( const double* coordinates )
  {
    ParentGeometryElement::assignNodeCoordinates( coordinates );
  }

  void YNodalGradientInterfaceFiniteElement::initializeYourself()
  {
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
      qp.normalProjection  = qp.normal * qp.normal.transpose();
      qp.tangentProjection = TensorDim::Identity() - qp.normalProjection;

      const double normalSeparation = qp.separationVector.dot( qp.normal );
      if ( qp.separationVector.norm() > geometryTolerance && normalSeparation <= geometryTolerance ) {
        throw std::invalid_argument(
          "YNodalGradientInterfaceFiniteElement: paired faces have no positive normal separation." );
      }

      qp.dA = qp.weight * qp.sqrtDetG;

      // B_a = [ B_Abar ; B_dA ], with B_Abar = <>-tangential = BmatAverage,
      // B_dA = [u,r]T = A+ - A- = [ -Bside | Bside ].
      qp.Ba.setZero();
      qp.Ba.block< nTensor, nDofD >( 0, 0 )                   = qp.BmatAverage;
      qp.Ba.block< nTensor, nSideDofU >( nTensor, 0 )         = -qp.BmatSide;
      qp.Ba.block< nTensor, nSideDofU >( nTensor, nSideDofU ) = qp.BmatSide;

      // N_g = N_t : block-diagonal Q1 interpolation of a 3-vector over 4 nodes.
      qp.Ng.setZero();
      for ( int A = 0; A < nInterfaceNodes; A++ )
        for ( int k = 0; k < 3; k++ )
          qp.Ng( k, 3 * A + k ) = qp.N( A );

      if ( qp.material ) {
        qp.material->setCharacteristicElementLength( std::sqrt( qp.sqrtDetG ) );
        qp.hFac = qp.material->getInterfaceThickness();
      }
    }
  }

  void YNodalGradientInterfaceFiniteElement::computeKernels( const double* QTotal_,
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
      const auto& Ba = qp.Ba;       // 18x24
      const auto& J  = qp.NmatJump; // 3x24
      const auto& Ng = qp.Ng;       // 3x12 (= N_t)
      // Read the interface thickness DIRECTLY from the material at evaluation
      // time. It must NOT be cached in initializeYourself: the host framework
      // (EdelweissFE) calls initializeYourself BEFORE assigning the material
      // section, so a cached value would silently stay at its default (1.0)
      // -- which scales the compatibility [u] = h g and ALL h-weighted terms
      // wrongly and makes the interface orders of magnitude too compliant.
      const double h  = qp.material->getInterfaceThickness();
      const double dA = qp.dA;

      // ---- increment inputs to the constitutive kernel ----
      const Eigen::Matrix< double, nDofD, 1 > dd  = dQ.segment< nDofD >( offD );
      const Eigen::Matrix< double, nDofG, 1 > dqg = dQ.segment< nDofG >( offG );

      const Eigen::Matrix< double, nA, 1 >   aIncr = Ba * dd;  // [dAbar(9); dDeltaA(9)]  contiguous
      const Eigen::Matrix< double, nDim, 1 > dg    = Ng * dqg; // 3

      Eigen::Matrix< double, 9, 1 >                    sAbar, sDeltaA;
      Eigen::Matrix< double, 3, 1 >                    tMat;
      Eigen::Matrix< double, nA, nA, Eigen::RowMajor > H_aa;
      Eigen::Matrix< double, nA, 3, Eigen::RowMajor >  H_ag;
      Eigen::Matrix< double, 3, nA, Eigen::RowMajor >  H_ga;
      Eigen::Matrix< double, 3, 3, Eigen::RowMajor >   H_gg;

      Material::KernelInput in{ aIncr.data(), aIncr.data() + 9, dg.data(), qp.normal.data() };
      Material::KernelOutput
        out{ sAbar.data(), sDeltaA.data(), tMat.data(), H_aa.data(), H_ag.data(), H_ga.data(), H_gg.data() };
      Material::TimeIncrement ti{ time, dT };

      qp.material->computeMixedKernel( qp.managedStateVars->materialStateVars.data(), in, out, ti );

      // ---- totals (current iterate) for the compatibility / traction terms ----
      const Eigen::Matrix< double, nDofD, 1 > dTot = QTotal.segment< nDofD >( offD );
      const Eigen::Matrix< double, nDofG, 1 > gTot = QTotal.segment< nDofG >( offG );
      const Eigen::Matrix< double, nDofT, 1 > tTot = QTotal.segment< nDofT >( offT );

      const Eigen::Matrix< double, nDim, 1 > jump = J * dTot;      // [u]
      const Eigen::Matrix< double, nDim, 1 > gVal = Ng * gTot;     // g
      const Eigen::Matrix< double, nDim, 1 > tVal = Ng * tTot;     // N_t p

      const Eigen::Matrix< double, nDim, 1 > cP = jump - h * gVal; // compatibility defect [u] - h g

      // diagnostics
      qp.managedStateVars->commonTraction     = tMat;
      qp.managedStateVars->surfaceStressAvg   = sAbar;
      qp.managedStateVars->surfaceStressJump  = sDeltaA;
      qp.managedStateVars->normalGradientAtQp = gVal;

      Eigen::Matrix< double, nA, 1 > sA;
      sA.segment< 9 >( 0 ) = sAbar;
      sA.segment< 9 >( 9 ) = sDeltaA;

      // ---- residuals (Pe -= R * dA) ----
      const Eigen::Matrix< double, nDofD, 1 > R_d = h * ( Ba.transpose() * sA ) + J.transpose() * tVal;
      const Eigen::Matrix< double, nDofG, 1 > R_q = h * ( Ng.transpose() * ( tMat - tVal ) );
      const Eigen::Matrix< double, nDofT, 1 > R_p = Ng.transpose() * cP;

      Pe.segment< nDofD >( offD ) -= R_d * dA;
      Pe.segment< nDofG >( offG ) -= R_q * dA;
      Pe.segment< nDofT >( offT ) -= R_p * dA;

      // ---- tangent blocks ----
      const Eigen::Matrix< double, nDofD, nDofD > K_dd = h * ( Ba.transpose() * H_aa * Ba );
      const Eigen::Matrix< double, nDofD, nDofG > K_dq = h * ( Ba.transpose() * H_ag * Ng );
      const Eigen::Matrix< double, nDofD, nDofT > K_dp = J.transpose() * Ng;
      const Eigen::Matrix< double, nDofG, nDofD > K_qd = h * ( Ng.transpose() * H_ga * Ba );
      const Eigen::Matrix< double, nDofG, nDofG > K_qq = h * ( Ng.transpose() * H_gg * Ng );
      const Eigen::Matrix< double, nDofG, nDofT > K_qp = -h * ( Ng.transpose() * Ng );
      const Eigen::Matrix< double, nDofT, nDofD > K_pd = Ng.transpose() * J;
      const Eigen::Matrix< double, nDofT, nDofG > K_pq = -h * ( Ng.transpose() * Ng );

      Ke.block< nDofD, nDofD >( offD, offD ) += K_dd * dA;
      Ke.block< nDofD, nDofG >( offD, offG ) += K_dq * dA;
      Ke.block< nDofD, nDofT >( offD, offT ) += K_dp * dA;
      Ke.block< nDofG, nDofD >( offG, offD ) += K_qd * dA;
      Ke.block< nDofG, nDofG >( offG, offG ) += K_qq * dA;
      Ke.block< nDofG, nDofT >( offG, offT ) += K_qp * dA;
      Ke.block< nDofT, nDofD >( offT, offD ) += K_pd * dA;
      Ke.block< nDofT, nDofG >( offT, offG ) += K_pq * dA;
      // K_pp = 0
    }
  }

  void YNodalGradientInterfaceFiniteElement::setInitialConditions( StateTypes state, const double* values )
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
      throw std::invalid_argument(
        MakeString() << __PRETTY_FUNCTION__ << ": invalid initial condition for YNodalGradientInterfaceFiniteElement" );
    }
  }

  void YNodalGradientInterfaceFiniteElement::computeDistributedLoad( MarmotElement::DistributedLoadTypes loadType,
                                                                     double*                             P,
                                                                     double*                             K,
                                                                     const int                           elementFace,
                                                                     const double*                       load,
                                                                     const double*                       QTotal,
                                                                     double                              time,
                                                                     double                              dT )
  {
    throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__
                                              << ": distributed loads are not implemented for "
                                                 "YNodalGradientInterfaceFiniteElement." );
  }

  void YNodalGradientInterfaceFiniteElement::computeBodyForce( double*       P,
                                                               double*       K,
                                                               const double* load,
                                                               const double* QTotal,
                                                               double        time,
                                                               double        dT )
  {
    throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__
                                              << ": body forces are not implemented for "
                                                 "YNodalGradientInterfaceFiniteElement." );
  }

  void YNodalGradientInterfaceFiniteElement::computeConsistentInertia( double* M )
  {
    throw std::runtime_error(
      MakeString() << __PRETTY_FUNCTION__ << ": inertia is not implemented for YNodalGradientInterfaceFiniteElement." );
  }

  void YNodalGradientInterfaceFiniteElement::computeLumpedInertia( double* M )
  {
    throw std::runtime_error(
      MakeString() << __PRETTY_FUNCTION__ << ": inertia is not implemented for YNodalGradientInterfaceFiniteElement." );
  }

  std::vector< double > YNodalGradientInterfaceFiniteElement::getCoordinatesAtCenter()
  {
    std::vector< double >   coords( nDim );
    Eigen::Map< VectorDim > coordsMap( coords.data() );
    const auto              centerXi = XiSized::Zero();
    const auto              Ncenter  = this->N( centerXi );
    const auto              Nmat     = this->NMatrix( Ncenter );
    coordsMap                        = Nmat * this->getSideCoordinates( 0 );
    return coords;
  }

  std::vector< std::vector< double > > YNodalGradientInterfaceFiniteElement::getCoordinatesAtQuadraturePoints()
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

  int YNodalGradientInterfaceFiniteElement::getNumberOfQuadraturePoints()
  {
    return static_cast< int >( qps.size() );
  }

} // namespace Marmot::Elements
