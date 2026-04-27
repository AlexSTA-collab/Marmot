#include "Marmot/InterfaceFiniteElement.h"

#include <cstddef>
#include <iostream>
#include <stdexcept>

namespace Marmot::Elements {

  template < int nDim, int nNodes >
  InterfaceFiniteElement< nDim, nNodes >::InterfaceFiniteElement(
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

  template < int nDim, int nNodes >
  int InterfaceFiniteElement< nDim, nNodes >::getNumberOfRequiredStateVars()
  {
    return qps[0].getNumberOfRequiredStateVars() * qps.size();
  }

  template < int nDim, int nNodes >
  std::vector< std::vector< std::string > > InterfaceFiniteElement< nDim, nNodes >::getNodeFields()
  {
    using namespace std;

    static vector< vector< string > > nodeFields;

    if ( nodeFields.empty() ) {
      for ( int i = 0; i < nNodes; i++ ) {
        nodeFields.push_back( vector< string >() );
        nodeFields[i].push_back( "displacement" );
      }
    }

    return nodeFields;
  }

  template < int nDim, int nNodes >
  std::vector< int > InterfaceFiniteElement< nDim, nNodes >::getDofIndicesPermutationPattern()
  {
    static std::vector< int > permutationPattern;

    if ( permutationPattern.empty() ) {
      for ( int i = 0; i < nNodes * nDim; i++ )
        permutationPattern.push_back( i );
    }

    return permutationPattern;
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::assignStateVars( double* stateVars, int nStateVars )
  {
    const int nQpStateVars = nStateVars / qps.size();

    for ( size_t i = 0; i < qps.size(); i++ ) {
      auto&   qp          = qps[i];
      double* qpStateVars = stateVars + ( i * nQpStateVars );
      qp.assignStateVars( qpStateVars, nQpStateVars );
    }
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::assignProperty( const ElementProperties& elementPropertiesInfo )
  {
    new ( &elementProperties ) Eigen::Map< const Eigen::VectorXd >( elementPropertiesInfo.elementProperties,
                                                                    elementPropertiesInfo.nElementProperties );
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::assignProperty( const MarmotMaterialSection& section )
  {
    for ( auto& qp : qps ) {
      qp.material = std::unique_ptr< Material >( dynamic_cast< Material* >(
        MarmotLibrary::MarmotMaterialHypoElasticInterfaceFactory::createMaterial( std::to_string(
                                                                                    section.materialCode ),
                                                                                  section.materialProperties,
                                                                                  section.nMaterialProperties,
                                                                                  elLabel ) ) );

      if ( !qp.material ) {
        throw std::invalid_argument(
          MakeString() << __PRETTY_FUNCTION__
                       << ": invalid material assigned; cannot cast to MarmotMaterialHypoElasticInterface!" );
      }
    }
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::assignMaterial( const std::string& materialName,
                                                               const double*      materialProperties,
                                                               int                nMaterialProperties )
  {
    for ( auto& qp : qps ) {
      qp.material = std::unique_ptr< Material >( dynamic_cast< Material* >(
        MarmotLibrary::MarmotMaterialHypoElasticInterfaceFactory::createMaterial( materialName,
                                                                                  materialProperties,
                                                                                  nMaterialProperties,
                                                                                  elLabel ) ) );

      if ( !qp.material ) {
        throw std::invalid_argument(
          MakeString() << __PRETTY_FUNCTION__
                       << ": invalid material assigned; cannot cast to MarmotMaterialHypoElasticInterface!" );
      }
    }
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::assignNodeCoordinates( const double* coordinates )
  {
    ParentGeometryElement::assignNodeCoordinates( coordinates );
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::initializeYourself()
  {
    const double thickness = elementProperties.size() > 0 ? elementProperties[0] : 1.0;

    for ( QuadraturePoint& qp : qps ) {
      const bool fullyProjectedB = ( nDim == 3 );
      const auto geom            = this->evaluateAt( qp.xi, 0, fullyProjectedB );

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

      qp.J0xW = qp.weight * qp.sqrtDetG * thickness;

      if ( qp.material ) {
        if constexpr ( nDim == 3 ) {
          qp.material->setCharacteristicElementLength( std::sqrt( qp.sqrtDetG ) );
        }
        else if constexpr ( nDim == 2 ) {
          qp.material->setCharacteristicElementLength( qp.sqrtDetG );
        }
      }
    }
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::computeYourself( const double* QTotal_,
                                                                const double* dQ_,
                                                                double*       Pe_,
                                                                double*       Ke_,
                                                                const double* time,
                                                                double        dT,
                                                                double&       pNewDT )
  {
    Eigen::Map< const RhsSized > dQ( dQ_ );
    Eigen::Map< KeSizedMatrix >  Ke( Ke_ );
    Eigen::Map< RhsSized >       Pe( Pe_ );

    constexpr int halfSize = nNodes * nDim / 2;

    using BWholeSized = Eigen::Matrix< double, nTensor, sizeLoadVector >;

    static bool printAllInterfaceDebugOnce = true;

    if ( printAllInterfaceDebugOnce ) {
      std::cout << "\n===== CXX Interface ALL GP geometry/material/K debug =====\n";
      std::cout.precision( 16 );
    }

    for ( size_t qpIndex = 0; qpIndex < qps.size(); ++qpIndex ) {
      QuadraturePoint& qp = qps[qpIndex];

      const auto& Nside = qp.NmatSide;
      const auto& Bside = qp.BmatSide;
      const auto& Njump = qp.NmatJump;
      const auto& Bavg  = qp.BmatAverage;

      BWholeSized Bbottom;
      BWholeSized Btop;
      BWholeSized Bsum;

      Bbottom.setZero();
      Btop.setZero();

      Bbottom.template block< nTensor, halfSize >( 0, 0 )     = Bside;
      Btop.template block< nTensor, halfSize >( 0, halfSize ) = Bside;

      Bsum = Bbottom + Btop;

      const auto dQBottom = dQ.template segment< halfSize >( 0 );
      const auto dQTop    = dQ.template segment< halfSize >( halfSize );

      InterfaceDisplSized dU_GPs;
      dU_GPs.template segment< nDim >( 0 )    = Nside * dQTop;
      dU_GPs.template segment< nDim >( nDim ) = Nside * dQBottom;

      InterfaceSurfaceGradSized dSurface_strain_GPs;
      dSurface_strain_GPs.template segment< nTensor >( 0 )       = Bside * dQTop;
      dSurface_strain_GPs.template segment< nTensor >( nTensor ) = Bside * dQBottom;

      ForceSized         force          = qp.managedStateVars->force;
      SurfaceStressSized surface_stress = qp.managedStateVars->surfaceStress;

      QMatrixSized Q_ij;
      ZMatrixSized Z_ijkl;
      HMatrixSized H_ijk;
      YMatrixSized Y_ijkl;

      Q_ij.setZero();
      Z_ijkl.setZero();
      H_ijk.setZero();
      Y_ijkl.setZero();

      qp.material->computeStress( force.data(),
                                  surface_stress.data(),
                                  Q_ij.data(),
                                  Z_ijkl.data(),
                                  H_ijk.data(),
                                  Y_ijkl.data(),
                                  dU_GPs.data(),
                                  dSurface_strain_GPs.data(),
                                  qp.normal.data(),
                                  time,
                                  dT,
                                  pNewDT );

      KeSizedMatrix K_jumpu_jumpv;
      KeSizedMatrix K_grad_s_u_grad_s_v_Z;
      KeSizedMatrix K_grad_s_u_grad_s_v_Y;
      KeSizedMatrix K_grad_s_u_jump_v;
      KeSizedMatrix K_jump_u_grad_s_v;
      KeSizedMatrix K_total_qp_unweighted;
      KeSizedMatrix K_total_qp_weighted;

      K_jumpu_jumpv         = Njump.transpose() * Q_ij * Njump;
      K_grad_s_u_grad_s_v_Z = Bavg.transpose() * Z_ijkl * Bavg;
      K_grad_s_u_grad_s_v_Y = Bavg.transpose() * Y_ijkl * Bavg;
      K_grad_s_u_jump_v     = Njump.transpose() * H_ijk * Bavg;
      K_jump_u_grad_s_v     = Bavg.transpose() * H_ijk.transpose() * Njump;

      K_total_qp_unweighted = K_jumpu_jumpv + K_grad_s_u_grad_s_v_Z + K_grad_s_u_grad_s_v_Y + K_grad_s_u_jump_v +
                              K_jump_u_grad_s_v;

      K_total_qp_weighted = K_total_qp_unweighted * qp.J0xW;

      if ( printAllInterfaceDebugOnce ) {
        Eigen::Matrix< double, nDim, 1 > xiPythonEquivalent;

        if constexpr ( nDim == 3 ) {
          xiPythonEquivalent( 0 ) = 0.5 * ( qp.xi( 0 ) + 1.0 );
          xiPythonEquivalent( 1 ) = 0.5 * ( qp.xi( 1 ) + 1.0 );
          xiPythonEquivalent( 2 ) = 0.0;
        }
        else {
          xiPythonEquivalent( 0 ) = 0.5 * ( qp.xi( 0 ) + 1.0 );
          if constexpr ( nDim > 1 )
            xiPythonEquivalent( 1 ) = 0.0;
        }

        std::cout << "\n--- CXX GP " << qpIndex << " ---\n";
        std::cout << "qp.index = " << qpIndex << "\n";
        std::cout << "qp.xi = " << qp.xi.transpose() << "\n";

        if constexpr ( nDim == 3 ) {
          std::cout << "qp.xi_python_equivalent = " << xiPythonEquivalent.template segment< 2 >( 0 ).transpose()
                    << "\n";
        }
        else {
          std::cout << "qp.xi_python_equivalent = " << xiPythonEquivalent( 0 ) << "\n";
        }

        std::cout << "qp.weight = " << qp.weight << "\n";
        std::cout << "qp.normal = " << qp.normal.transpose() << "\n";
        std::cout << "qp.normal.norm = " << qp.normal.norm() << "\n";
        std::cout << "qp.sqrtDetG = " << qp.sqrtDetG << "\n";
        std::cout << "qp.J0xW = " << qp.J0xW << "\n";

        std::cout << "\ndNdXi shape = (" << qp.dNdXi.rows() << ", " << qp.dNdXi.cols() << ")\n";
        std::cout << "dNdXi =\n" << qp.dNdXi << "\n";

        std::cout << "\ngradN shape = (" << qp.gradN.rows() << ", " << qp.gradN.cols() << ")\n";
        std::cout << "gradN =\n" << qp.gradN << "\n";

        std::cout << "\nNside shape = (" << Nside.rows() << ", " << Nside.cols() << ")\n";
        std::cout << "Nside =\n" << Nside << "\n";

        std::cout << "\nNjump shape = (" << Njump.rows() << ", " << Njump.cols() << ")\n";
        std::cout << "Njump =\n" << Njump << "\n";

        std::cout << "\nBside shape = (" << Bside.rows() << ", " << Bside.cols() << ")\n";
        std::cout << "Bside =\n" << Bside << "\n";

        std::cout << "\nBavg shape = (" << Bavg.rows() << ", " << Bavg.cols() << ")\n";
        std::cout << "Bavg =\n" << Bavg << "\n";

        std::cout << "\nBbottom shape = (" << Bbottom.rows() << ", " << Bbottom.cols() << ")\n";
        std::cout << "Bbottom =\n" << Bbottom << "\n";

        std::cout << "\nBtop shape = (" << Btop.rows() << ", " << Btop.cols() << ")\n";
        std::cout << "Btop =\n" << Btop << "\n";

        std::cout << "\nBsum shape = (" << Bsum.rows() << ", " << Bsum.cols() << ")\n";
        std::cout << "Bsum =\n" << Bsum << "\n";

        std::cout << "\ndU_GPs shape = (" << dU_GPs.rows() << ", " << dU_GPs.cols() << ")\n";
        std::cout << "dU_GPs =\n" << dU_GPs << "\n";

        std::cout << "\ndSurface_strain_GPs shape = (" << dSurface_strain_GPs.rows() << ", "
                  << dSurface_strain_GPs.cols() << ")\n";
        std::cout << "dSurface_strain_GPs =\n" << dSurface_strain_GPs << "\n";

        std::cout << "\nforce shape = (" << force.rows() << ", " << force.cols() << ")\n";
        std::cout << "force =\n" << force << "\n";

        std::cout << "\nsurface_stress shape = (" << surface_stress.rows() << ", " << surface_stress.cols() << ")\n";
        std::cout << "surface_stress =\n" << surface_stress << "\n";

        std::cout << "\nQ_ij shape = (" << Q_ij.rows() << ", " << Q_ij.cols() << ")\n";
        std::cout << "Q_ij =\n" << Q_ij << "\n";

        std::cout << "\nZ_ijkl shape = (" << Z_ijkl.rows() << ", " << Z_ijkl.cols() << ")\n";
        std::cout << "Z_ijkl =\n" << Z_ijkl << "\n";

        std::cout << "\nY_ijkl shape = (" << Y_ijkl.rows() << ", " << Y_ijkl.cols() << ")\n";
        std::cout << "Y_ijkl =\n" << Y_ijkl << "\n";

        std::cout << "\nH_ijk shape = (" << H_ijk.rows() << ", " << H_ijk.cols() << ")\n";
        std::cout << "H_ijk =\n" << H_ijk << "\n";

        std::cout << "\nK_jumpu_jumpv shape = (" << K_jumpu_jumpv.rows() << ", " << K_jumpu_jumpv.cols() << ")\n";
        std::cout << "K_jumpu_jumpv =\n" << K_jumpu_jumpv << "\n";

        std::cout << "\nK_grad_s_u_grad_s_v_Z shape = (" << K_grad_s_u_grad_s_v_Z.rows() << ", "
                  << K_grad_s_u_grad_s_v_Z.cols() << ")\n";
        std::cout << "K_grad_s_u_grad_s_v_Z =\n" << K_grad_s_u_grad_s_v_Z << "\n";

        std::cout << "\nK_grad_s_u_grad_s_v_Y shape = (" << K_grad_s_u_grad_s_v_Y.rows() << ", "
                  << K_grad_s_u_grad_s_v_Y.cols() << ")\n";
        std::cout << "K_grad_s_u_grad_s_v_Y =\n" << K_grad_s_u_grad_s_v_Y << "\n";

        std::cout << "\nK_grad_s_u_jump_v shape = (" << K_grad_s_u_jump_v.rows() << ", " << K_grad_s_u_jump_v.cols()
                  << ")\n";
        std::cout << "K_grad_s_u_jump_v =\n" << K_grad_s_u_jump_v << "\n";

        std::cout << "\nK_jump_u_grad_s_v shape = (" << K_jump_u_grad_s_v.rows() << ", " << K_jump_u_grad_s_v.cols()
                  << ")\n";
        std::cout << "K_jump_u_grad_s_v =\n" << K_jump_u_grad_s_v << "\n";

        std::cout << "\nK_total_qp_unweighted shape = (" << K_total_qp_unweighted.rows() << ", "
                  << K_total_qp_unweighted.cols() << ")\n";
        std::cout << "K_total_qp_unweighted =\n" << K_total_qp_unweighted << "\n";

        std::cout << "\nK_total_qp_weighted shape = (" << K_total_qp_weighted.rows() << ", "
                  << K_total_qp_weighted.cols() << ")\n";
        std::cout << "K_total_qp_weighted =\n" << K_total_qp_weighted << "\n";
      }

      if ( pNewDT < 1.0 )
        return;

      qp.managedStateVars->force         = force;
      qp.managedStateVars->surfaceStress = surface_stress;
      qp.managedStateVars->displacement += dU_GPs;
      qp.managedStateVars->surfaceStrain += dSurface_strain_GPs;

      Pe -= Njump.transpose() * force * qp.J0xW;
      Pe -= Bavg.transpose() * surface_stress * qp.J0xW;

      Ke += K_total_qp_weighted;
    }

    if ( printAllInterfaceDebugOnce ) {
      std::cout << "\n===== END CXX Interface ALL GP geometry/material/K debug =====\n";
      std::cout.flush();
      printAllInterfaceDebugOnce = false;
    }
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::setInitialConditions( StateTypes state, const double* values )
  {
    switch ( state ) {
    case MarmotElement::MarmotMaterialInitialization: {
      for ( QuadraturePoint& qp : qps ) {
        qp.material->initializeYourself();
      }
      break;
    }

    case MarmotElement::MarmotMaterialStateVars: {
      throw std::invalid_argument( "Please use initializeStateVars directly on material" );
    }

    default:
      throw std::invalid_argument( MakeString()
                                   << __PRETTY_FUNCTION__ << ": invalid initial condition for InterfaceFiniteElement" );
    }
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::computeDistributedLoad( MarmotElement::DistributedLoadTypes loadType,
                                                                       double*                             P,
                                                                       double*                             K,
                                                                       const int                           elementFace,
                                                                       const double*                       load,
                                                                       const double*                       QTotal,
                                                                       const double*                       time,
                                                                       double                              dT )
  {
    throw std::invalid_argument(
      MakeString() << __PRETTY_FUNCTION__ << ": distributed loads are not implemented for InterfaceFiniteElement." );
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::computeBodyForce( double*       P,
                                                                 double*       K,
                                                                 const double* load,
                                                                 const double* QTotal,
                                                                 const double* time,
                                                                 double        dT )
  {
    throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__
                                              << ": body forces are not implemented for InterfaceFiniteElement." );
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::computeConsistentInertia( double* M )
  {
    Eigen::Map< KeSizedMatrix > Me( M );
    Me.setZero();
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::computeLumpedInertia( double* M )
  {
    Eigen::Map< RhsSized > Me( M );
    Me.setZero();
  }

  template < int nDim, int nNodes >
  std::vector< double > InterfaceFiniteElement< nDim, nNodes >::getCoordinatesAtCenter()
  {
    std::vector< double > coords( nDim );

    Eigen::Map< VectorDim > coordsMap( coords.data() );

    const auto centerXi = XiSized::Zero();
    const auto Ncenter  = this->N( centerXi );
    const auto Nmat     = this->NMatrix( Ncenter );

    const auto xSide = this->getSideCoordinates( 0 );

    coordsMap = Nmat * xSide;

    return coords;
  }

  template < int nDim, int nNodes >
  std::vector< std::vector< double > > InterfaceFiniteElement< nDim, nNodes >::getCoordinatesAtQuadraturePoints()
  {
    std::vector< std::vector< double > > listedCoords;

    for ( const auto& qp : qps ) {
      std::vector< double > coords( nDim );

      Eigen::Map< VectorDim > coordsMap( coords.data() );

      const auto xSide = this->getSideCoordinates( 0 );
      coordsMap        = qp.NmatSide * xSide;

      listedCoords.push_back( coords );
    }

    return listedCoords;
  }

  template < int nDim, int nNodes >
  int InterfaceFiniteElement< nDim, nNodes >::getNumberOfQuadraturePoints()
  {
    return static_cast< int >( qps.size() );
  }

  template class InterfaceFiniteElement< 2, 4 >;
  template class InterfaceFiniteElement< 3, 8 >;

} // namespace Marmot::Elements