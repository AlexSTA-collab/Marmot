/* ---------------------------------------------------------------------
 *  marmot - MAteRialMOdellingToolbox
 *  Alexandros Stathas alexandros.stathas@boku.ac.at
 *  LGPL 2.1 or later; see LICENSE.md at the top level directory of marmot.
 * ---------------------------------------------------------------------
 */

#include "Marmot/ZStabPressureInterfaceFiniteElement.h"

namespace Marmot::Elements {

  ZStabPressureInterfaceFiniteElement::ZStabPressureInterfaceFiniteElement(
    int                                         elementID,
    FiniteElement::Quadrature::IntegrationTypes integrationType,
    SectionType                                 sectionType_ )
    : ParentGeometryElement(),
      elementProperties( Eigen::Map< const Eigen::VectorXd >( nullptr, 0 ) ),
      elLabel( elementID ),
      sectionType( sectionType_ )
  {
    for ( const auto& qpInfo : FiniteElement::Quadrature::getGaussPointInfo( this->shape, integrationType ) )
      qps.push_back( QuadraturePoint( qpInfo.xi, qpInfo.weight ) );
  }

  int ZStabPressureInterfaceFiniteElement::getNumberOfRequiredStateVars()
  {
    return qps[0].getNumberOfRequiredStateVars() * qps.size();
  }

  std::vector< std::vector< std::string > > ZStabPressureInterfaceFiniteElement::getNodeFields()
  {
    static std::vector< std::vector< std::string > > nodeFields;
    if ( nodeFields.empty() )
      for ( int i = 0; i < nNodes; i++ ) {
        nodeFields.push_back( { "displacement" } );
        if ( i < nInterfaceNodes ) {
          nodeFields[i].push_back( "normalGradientJump" );
          nodeFields[i].push_back( "interfacePressure" );
          nodeFields[i].push_back( "interfacePressureJump" );
        }
      }
    return nodeFields;
  }

  std::vector< int > ZStabPressureInterfaceFiniteElement::getDofIndicesPermutationPattern()
  {
    static std::vector< int > perm;
    if ( perm.empty() ) {
      perm.resize( sizeLoadVector );
      // node-major: bottom node A carries [u(3), g(3), pbar(1), [p](1)] = 8 slots
      for ( int A = 0; A < nInterfaceNodes; A++ ) {
        for ( int c = 0; c < 3; c++ )
          perm[offU + 3 * A + c] = 8 * A + c;
        for ( int k = 0; k < 3; k++ )
          perm[offG + 3 * A + k] = 8 * A + 3 + k;
        perm[offPm + A] = 8 * A + 6;
        perm[offPj + A] = 8 * A + 7;
      }
      for ( int A = 0; A < nInterfaceNodes; A++ )
        for ( int c = 0; c < 3; c++ )
          perm[offU + nSideDofU + 3 * A + c] = 32 + 3 * A + c;
    }
    return perm;
  }

  void ZStabPressureInterfaceFiniteElement::assignStateVars( double* stateVars, int nStateVars )
  {
    const int n = nStateVars / qps.size();
    for ( size_t i = 0; i < qps.size(); i++ )
      qps[i].assignStateVars( stateVars + i * n, n );
  }

  void ZStabPressureInterfaceFiniteElement::assignProperty( const ElementProperties& info )
  {
    new ( &elementProperties ) Eigen::Map< const Eigen::VectorXd >( info.elementProperties, info.nElementProperties );
    if ( elementProperties.size() > 1 && !std::getenv( "MARMOT_STABP_GAMMA" ) )
      stabGamma = elementProperties[1];
  }

  void ZStabPressureInterfaceFiniteElement::applyMaterialSettings( QuadraturePoint& qp )
  {
    if ( !qp.material || qp.sqrtDetG <= 0.0 )
      return;
    qp.material->setCharacteristicElementLength( std::sqrt( qp.sqrtDetG ) );
    if ( elementProperties.size() > 2 )
      qp.material->setGradientJumpRegularization( elementProperties[2] );
  }

  void ZStabPressureInterfaceFiniteElement::assignProperty( const MarmotMaterialSection& s )
  {
    for ( auto& qp : qps ) {
      qp.material = std::make_unique< Material >( s.materialName,
                                                  s.materialProperties,
                                                  s.nMaterialProperties,
                                                  elLabel );
      applyMaterialSettings( qp );
    }
  }

  void ZStabPressureInterfaceFiniteElement::assignMaterial( const std::string& name, const double* props, int n )
  {
    for ( auto& qp : qps ) {
      qp.material = std::make_unique< Material >( name, props, n, elLabel );
      applyMaterialSettings( qp );
    }
  }

  void ZStabPressureInterfaceFiniteElement::assignNodeCoordinates( const double* c )
  {
    ParentGeometryElement::assignNodeCoordinates( c );
  }

  void ZStabPressureInterfaceFiniteElement::initializeYourself()
  {
    const double thickness = elementProperties.size() > 0 ? elementProperties[0] : 1.0;
    for ( QuadraturePoint& qp : qps ) {
      const auto geom = this->evaluateAt( qp.xi, 0 );
      qp.N            = geom.N;
      qp.dNdXi        = geom.dNdXi;
      qp.J            = geom.J;
      qp.G            = geom.G;
      qp.sqrtDetG     = geom.sqrtDetG;
      qp.detJ         = geom.sqrtDetG;
      qp.gradN        = geom.gradN;
      qp.normal       = geom.n;
      qp.NmatSide     = geom.NmatSide;
      qp.BmatSide     = geom.BmatSide;
      qp.NmatJump     = geom.NmatJump;
      qp.BmatAverage  = geom.BmatAverage;

      const VectorDim xB   = qp.NmatSide * this->getSideCoordinates( 0 );
      const VectorDim xT   = qp.NmatSide * this->getSideCoordinates( 1 );
      qp.separationVector  = xT - xB;
      constexpr double tol = 1.0e-12;
      if ( qp.separationVector.norm() > tol && qp.separationVector.dot( qp.normal ) < 0.0 )
        qp.normal *= -1.0;
      qp.normalProjection     = qp.normal * qp.normal.transpose();
      qp.tangentProjection    = TensorDim::Identity() - qp.normalProjection;
      qp.normalSeparation     = qp.separationVector.dot( qp.normal );
      qp.tangentialSeparation = qp.tangentProjection * qp.separationVector;
      if ( qp.separationVector.norm() > tol && qp.normalSeparation <= tol )
        throw std::invalid_argument( "ZStabPressureInterfaceFiniteElement: no positive normal separation." );

      qp.J0xW     = qp.weight * qp.sqrtDetG * thickness;
      qp.meshSize = std::sqrt( qp.sqrtDetG ); // mesh size h_e for the stabilisation scale

      qp.BX.setZero();
      qp.BX.block< 3, nDofU >( 0, 0 )                     = qp.NmatJump;
      qp.BX.block< nTensor, nSideDofU >( 3, nSideDofU )   = qp.BmatSide;
      qp.BX.block< nTensor, nSideDofU >( 3 + nTensor, 0 ) = qp.BmatSide;

      qp.Ng.setZero();
      for ( int A = 0; A < nInterfaceNodes; A++ )
        for ( int k = 0; k < 3; k++ )
          qp.Ng( k, 3 * A + k ) = qp.N( A );

      for ( int A = 0; A < nInterfaceNodes; A++ ) {
        qp.Np( 0, A ) = qp.N( A );
        for ( int i = 0; i < nDim; ++i )
          qp.gradNp( i, A ) = qp.gradN( i, A );
      }

      // MINI bubble b = (1-xi^2)(1-eta^2): vanishes on the element boundary, so
      // beta is element-local and can be condensed.
      {
        const double xi = qp.xi( 0 ), eta = qp.xi( 1 );
        qp.bubble = ( 1.0 - xi * xi ) * ( 1.0 - eta * eta );
        Eigen::Matrix< double, 2, 1 > dbdxi;
        dbdxi << -2.0 * xi * ( 1.0 - eta * eta ), -2.0 * eta * ( 1.0 - xi * xi );
        qp.gradBubble = qp.J * qp.G.inverse() * dbdxi;
      }

      applyMaterialSettings( qp );
    }
  }

  void ZStabPressureInterfaceFiniteElement::computeKernels( const double* QTotal_,
                                                            const double* dQ_,
                                                            double*       Pe_,
                                                            double*       Ke_,
                                                            double        time,
                                                            double        dT )
  {
    Eigen::Map< const RhsSized > QTotal( QTotal_ ), dQ( dQ_ );
    Eigen::Map< KeSizedMatrix >  Ke( Ke_ );
    Eigen::Map< RhsSized >       Pe( Pe_ );

    using V3  = Eigen::Matrix< double, 3, 1 >;
    using M33 = Eigen::Matrix< double, 3, 3 >;

    /** Everything the element needs at a given bubble amplitude beta. */
    struct Acc {
      RhsSized                                   R;   // residual over the 44 global dofs
      KeSizedMatrix                              K;
      V3                                         Rb;  // d(energy)/d(beta)
      M33                                        Kbb;
      Eigen::Matrix< double, 3, sizeLoadVector > Kbq; // d Rb / d q
      Eigen::Matrix< double, sizeLoadVector, 3 > Kqb; // d R  / d beta
    };

    // The material COMMITS its state on every computeStress call, so the inner
    // Newton on beta -- which evaluates at many trial betas -- would accumulate
    // plastic state across iterations. Snapshot the committed state once and
    // restore it before every trial, so each is measured from the same t_n.
    std::vector< std::vector< double > > committed;
    committed.reserve( qps.size() );
    for ( auto& qp : qps ) {
      const auto& m = qp.managedStateVars->materialStateVars;
      committed.emplace_back( m.data(), m.data() + m.size() );
    }

    auto assembleAll = [&]( const V3& beta ) -> Acc {
      for ( size_t i = 0; i < qps.size(); ++i )
        std::copy( committed[i].begin(), committed[i].end(), qps[i].managedStateVars->materialStateVars.data() );
      Acc A;
      A.R.setZero();
      A.K.setZero();
      A.Rb.setZero();
      A.Kbb.setZero();
      A.Kbq.setZero();
      A.Kqb.setZero();

      for ( QuadraturePoint& qp : qps ) {
        const auto& BX = qp.BX;
        const auto& Ng = qp.Ng;
        const auto& Np = qp.Np;

        const Eigen::Matrix< double, nDofU, 1 > du  = dQ.segment< nDofU >( offU );
        const Eigen::Matrix< double, nDofG, 1 > dg  = dQ.segment< nDofG >( offG );
        const double                            dpm = ( Np * dQ.segment< nDofP >( offPm ) )( 0 );
        const double                            dpj = ( Np * dQ.segment< nDofP >( offPj ) )( 0 );

        // MINI: ubar += b*beta  =>  Abar_ij += beta_i (grad_s b)_j, EQUALLY on both faces.
        Eigen::Matrix< double, 9, 3 > Bbub = Eigen::Matrix< double, 9, 3 >::Zero();
        for ( int i = 0; i < 3; ++i )
          for ( int j = 0; j < 3; ++j )
            Bbub( 3 * i + j, i ) = qp.gradBubble( j );
        // beta enters x only through the two surface-gradient blocks
        Eigen::Matrix< double, nX, 3 > Bxb = Eigen::Matrix< double, nX, 3 >::Zero();
        Bxb.block< 9, 3 >( 3, 0 )          = Bbub;
        Bxb.block< 9, 3 >( 12, 0 )         = Bbub;

        const Eigen::Matrix< double, nX, 1 > dX = BX * du + Bxb * beta;
        const Eigen::Matrix< double, nZ, 1 > dz = Ng * dg;

        Eigen::Matrix< double, 6, 1 > dU;
        dU.segment< 3 >( 0 ) = qp.NmatSide * du.segment< nSideDofU >( nSideDofU );
        dU.segment< 3 >( 3 ) = qp.NmatSide * du.segment< nSideDofU >( 0 );
        Eigen::Matrix< double, 18, 1 > dSurf;
        dSurf.segment< 9 >( 0 ) = dX.segment< 9 >( 3 );
        dSurf.segment< 9 >( 9 ) = dX.segment< 9 >( 12 );

        Eigen::Matrix< double, 3, 1 >                    f, rZ;
        Eigen::Matrix< double, 9, 1 >                    sP, sM;
        double                                           rPm = 0.0, rPj = 0.0;
        Eigen::Matrix< double, nX, nX, Eigen::RowMajor > K_xx;
        Eigen::Matrix< double, nX, nZ, Eigen::RowMajor > K_xz;
        Eigen::Matrix< double, nX, 1 >                   K_xpm, K_xpj;
        Eigen::Matrix< double, nZ, nX, Eigen::RowMajor > K_zx;
        Eigen::Matrix< double, nZ, nZ, Eigen::RowMajor > K_zz;
        Eigen::Matrix< double, nZ, 1 >                   K_zpm, K_zpj;
        Eigen::Matrix< double, 1, nX >                   K_pmx, K_pjx;
        Eigen::Matrix< double, 1, nZ >                   K_pmz, K_pjz;
        double                                           K_pmpm = 0.0, K_pjpj = 0.0;

        Material::Response resp{ f.data(), sP.data(), sM.data(), rZ.data(), &rPm, &rPj };
        Material::Tangents tan{ K_xx.data(),
                                K_xz.data(),
                                K_xpm.data(),
                                K_xpj.data(),
                                K_zx.data(),
                                K_zz.data(),
                                K_zpm.data(),
                                K_zpj.data(),
                                K_pmx.data(),
                                K_pmz.data(),
                                &K_pmpm,
                                K_pjx.data(),
                                K_pjz.data(),
                                &K_pjpj };
        Material::Deformation
          def{ dU.data(), dSurf.data(), dz.data(), qp.normal.data(), qp.separationVector.data(), dpm, dpj };
        Material::TimeIncrement ti{ time, dT };
        qp.material->computeStress( qp.managedStateVars->materialStateVars.data(), resp, tan, def, ti );

        qp.managedStateVars->generalizedForce   = f;
        qp.managedStateVars->surfaceStressPlus  = sP;
        qp.managedStateVars->surfaceStressMinus = sM;
        qp.managedStateVars->tractionImbalance  = rZ;
        qp.managedStateVars->normalGradientJump = Ng * QTotal.segment< nDofG >( offG );
        qp.managedStateVars->pressureMean( 0 )  = ( Np * QTotal.segment< nDofP >( offPm ) )( 0 );
        qp.managedStateVars->pressureJump( 0 )  = ( Np * QTotal.segment< nDofP >( offPj ) )( 0 );

        Eigen::Matrix< double, nX, 1 > pX;
        pX.segment< 3 >( 0 )  = f;
        pX.segment< 9 >( 3 )  = sP;
        pX.segment< 9 >( 12 ) = sM;

        const double                                h    = qp.material->getInterfaceThickness();
        const double                                mu   = qp.material->getShearModulus();
        const double                                stab = stabGamma * qp.meshSize * qp.meshSize / ( 2.0 * mu );
        const Eigen::Matrix< double, nDofP, nDofP > S    = stab * ( qp.gradNp.transpose() * qp.gradNp );

        const Eigen::Matrix< double, nDofP, 1 > pmTot = QTotal.segment< nDofP >( offPm );
        const Eigen::Matrix< double, nDofP, 1 > pjTot = QTotal.segment< nDofP >( offPj );

        A.R.segment< nDofU >( offU ) += BX.transpose() * pX * qp.J0xW;
        A.R.segment< nDofG >( offG ) += Ng.transpose() * rZ * qp.J0xW;
        A.R.segment< nDofP >( offPm ) += h * ( Np.transpose() * rPm + S * pmTot ) * qp.J0xW;
        A.R.segment< nDofP >( offPj ) += h * ( Np.transpose() * rPj + S * pjTot ) * qp.J0xW;

        A.K.block< nDofU, nDofU >( offU, offU ) += BX.transpose() * K_xx * BX * qp.J0xW;
        A.K.block< nDofU, nDofG >( offU, offG ) += BX.transpose() * K_xz * Ng * qp.J0xW;
        A.K.block< nDofU, nDofP >( offU, offPm ) += BX.transpose() * K_xpm * Np * qp.J0xW;
        A.K.block< nDofU, nDofP >( offU, offPj ) += BX.transpose() * K_xpj * Np * qp.J0xW;
        A.K.block< nDofG, nDofU >( offG, offU ) += Ng.transpose() * K_zx * BX * qp.J0xW;
        A.K.block< nDofG, nDofG >( offG, offG ) += Ng.transpose() * K_zz * Ng * qp.J0xW;
        A.K.block< nDofG, nDofP >( offG, offPm ) += Ng.transpose() * K_zpm * Np * qp.J0xW;
        A.K.block< nDofG, nDofP >( offG, offPj ) += Ng.transpose() * K_zpj * Np * qp.J0xW;
        A.K.block< nDofP, nDofU >( offPm, offU ) += h * ( Np.transpose() * K_pmx * BX ) * qp.J0xW;
        A.K.block< nDofP, nDofG >( offPm, offG ) += h * ( Np.transpose() * K_pmz * Ng ) * qp.J0xW;
        A.K.block< nDofP, nDofP >( offPm, offPm ) += h * ( Np.transpose() * K_pmpm * Np + S ) * qp.J0xW;
        A.K.block< nDofP, nDofU >( offPj, offU ) += h * ( Np.transpose() * K_pjx * BX ) * qp.J0xW;
        A.K.block< nDofP, nDofG >( offPj, offG ) += h * ( Np.transpose() * K_pjz * Ng ) * qp.J0xW;
        A.K.block< nDofP, nDofP >( offPj, offPj ) += h * ( Np.transpose() * K_pjpj * Np + S ) * qp.J0xW;

        // ---- bubble blocks: beta reaches everything through x ----
        A.Rb += Bxb.transpose() * pX * qp.J0xW;
        A.Kbb += Bxb.transpose() * K_xx * Bxb * qp.J0xW;
        A.Kbq.block< 3, nDofU >( 0, offU ) += Bxb.transpose() * K_xx * BX * qp.J0xW;
        A.Kbq.block< 3, nDofG >( 0, offG ) += Bxb.transpose() * K_xz * Ng * qp.J0xW;
        A.Kbq.block< 3, nDofP >( 0, offPm ) += Bxb.transpose() * K_xpm * Np * qp.J0xW;
        A.Kbq.block< 3, nDofP >( 0, offPj ) += Bxb.transpose() * K_xpj * Np * qp.J0xW;
        A.Kqb.block< nDofU, 3 >( offU, 0 ) += BX.transpose() * K_xx * Bxb * qp.J0xW;
        A.Kqb.block< nDofG, 3 >( offG, 0 ) += Ng.transpose() * K_zx * Bxb * qp.J0xW;
        A.Kqb.block< nDofP, 3 >( offPm, 0 ) += h * ( Np.transpose() * K_pmx * Bxb ) * qp.J0xW;
        A.Kqb.block< nDofP, 3 >( offPj, 0 ) += h * ( Np.transpose() * K_pjx * Bxb ) * qp.J0xW;
      }
      return A;
    };

    // ---- inner Newton on beta: R_beta(beta) = 0. Nonlinear, because beta enters
    //      the strain and therefore the plastic return map. Backtracking line
    //      search, as in YIQUAD4_STABP_MINI.
    V3  beta = bubbleDisabled() ? V3::Zero() : V3( qps[0].managedStateVars->bubbleAlpha );
    Acc A    = assembleAll( beta );

    if ( !bubbleDisabled() ) {
      const double scale = std::max( 1.0, A.R.segment< nDofU >( offU ).norm() );
      for ( int it = 0; it < 20 && A.Rb.norm() > 1.0e-12 * scale; ++it ) {
        const V3 dBeta = A.Kbb.fullPivLu().solve( V3( -A.Rb ) );
        if ( !dBeta.allFinite() )
          break;
        double alpha    = 1.0;
        bool   accepted = false;
        for ( int ls = 0; ls < 12; ++ls ) {
          try {
            Acc trial = assembleAll( V3( beta + alpha * dBeta ) );
            if ( trial.Rb.norm() < A.Rb.norm() ) {
              beta += alpha * dBeta;
              A        = trial;
              accepted = true;
              break;
            }
          }
          catch ( ... ) {
          }
          alpha *= 0.5;
        }
        if ( !accepted )
          break;
      }
      qps[0].managedStateVars->bubbleAlpha = beta;
    }
    // Leave the persistent state at the ACCEPTED beta: the loop above may have
    // ended on a rejected trial, whose state must not be what gets committed.
    A = assembleAll( beta );

    // ---- static condensation of beta ----
    Ke = A.K;
    Pe = -A.R;
    if ( !bubbleDisabled() ) {
      const Eigen::FullPivLU< M33 > lu( A.Kbb );
      if ( lu.isInvertible() ) {
        Ke -= A.Kqb * lu.solve( A.Kbq );
        Pe += A.Kqb * lu.solve( A.Rb );
      }
    }
  }

  void ZStabPressureInterfaceFiniteElement::setInitialConditions( StateTypes state, const double* )
  {
    switch ( state ) {
    case MarmotElement::MarmotMaterialInitialization:
      for ( QuadraturePoint& qp : qps )
        qp.material->initializeYourself( qp.managedStateVars->materialStateVars.data(),
                                         qp.managedStateVars->materialStateVars.size() );
      break;
    case MarmotElement::MarmotMaterialStateVars:
      throw std::invalid_argument( "Please use initializeStateVars directly on material" );
    default: throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid initial condition" );
    }
  }

  void ZStabPressureInterfaceFiniteElement::computeDistributedLoad( MarmotElement::DistributedLoadTypes,
                                                                    double*,
                                                                    double*,
                                                                    const int,
                                                                    const double*,
                                                                    const double*,
                                                                    double,
                                                                    double )
  {
    throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": distributed loads not implemented." );
  }
  void ZStabPressureInterfaceFiniteElement::computeBodyForce( double*,
                                                              double*,
                                                              const double*,
                                                              const double*,
                                                              double,
                                                              double )
  {
    throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": body forces not implemented." );
  }
  void ZStabPressureInterfaceFiniteElement::computeConsistentInertia( double* )
  {
    throw std::runtime_error( MakeString() << __PRETTY_FUNCTION__ << ": inertia not implemented." );
  }
  void ZStabPressureInterfaceFiniteElement::computeLumpedInertia( double* )
  {
    throw std::runtime_error( MakeString() << __PRETTY_FUNCTION__ << ": inertia not implemented." );
  }

  std::vector< double > ZStabPressureInterfaceFiniteElement::getCoordinatesAtCenter()
  {
    std::vector< double >   c( nDim );
    Eigen::Map< VectorDim > m( c.data() );
    m = this->NMatrix( this->N( XiSized::Zero() ) ) * this->getSideCoordinates( 0 );
    return c;
  }
  std::vector< std::vector< double > > ZStabPressureInterfaceFiniteElement::getCoordinatesAtQuadraturePoints()
  {
    std::vector< std::vector< double > > out;
    for ( const auto& qp : qps ) {
      std::vector< double >   c( nDim );
      Eigen::Map< VectorDim > m( c.data() );
      m = qp.NmatSide * this->getSideCoordinates( 0 );
      out.push_back( c );
    }
    return out;
  }
  int ZStabPressureInterfaceFiniteElement::getNumberOfQuadraturePoints()
  {
    return (int)qps.size();
  }

} // namespace Marmot::Elements
