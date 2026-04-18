#include "Marmot/InterfaceFiniteElement.h"
#include <cstddef>
#include <iostream>

namespace Marmot::Elements {

  template < int nDim, int nNodes >
  InterfaceFiniteElement< nDim, nNodes >::InterfaceFiniteElement(
    int                                         elementID,
    FiniteElement::Quadrature::IntegrationTypes integrationType,
    SectionType                                 sectionType )
    : coordinates( nullptr ), elementProperties( nullptr, 0 ), elLabel( elementID ), sectionType( sectionType )
  {
    auto qpInfos = FiniteElement::Quadrature::getGaussPointInfo( shape, integrationType );
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
    if ( nodeFields.empty() )
      for ( int i = 0; i < nNodes; i++ ) {
        nodeFields.push_back( vector< string >() );
        nodeFields[i].push_back( "displacement" );
      }

    return nodeFields;
  }

  template < int nDim, int nNodes >
  std::vector< int > InterfaceFiniteElement< nDim, nNodes >::getDofIndicesPermutationPattern()
  {
    static std::vector< int > permutationPattern;
    if ( permutationPattern.empty() )
      for ( int i = 0; i < nNodes * nDim; i++ )
        permutationPattern.push_back( i );

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
      qp.material = std::unique_ptr< MarmotMaterialHypoElasticInterface >(
        dynamic_cast< MarmotMaterialHypoElasticInterface* >(
          MarmotLibrary::MarmotMaterialHypoElasticInterfaceFactory::createMaterial( std::to_string(
                                                                                      section.materialCode ),
                                                                                    section.materialProperties,
                                                                                    section.nMaterialProperties,
                                                                                    elLabel ) ) );

      if ( !qp.material )
        throw std::invalid_argument(
          MakeString() << __PRETTY_FUNCTION__
                       << ": invalid material assigned; cannot cast to MarmotMaterialHypoElasticInterface!" );

      if constexpr ( nDim == 3 )
        qp.material->setCharacteristicElementLength( std::cbrt( 8 * qp.detJ ) );
      if constexpr ( nDim == 2 )
        qp.material->setCharacteristicElementLength( std::sqrt( 4 * qp.detJ ) );
      if constexpr ( nDim == 1 )
        qp.material->setCharacteristicElementLength( 2 * qp.detJ );
    }
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::assignMaterial( const std::string& materialName,
                                                               const double*      materialProperties,
                                                               int                nMaterialProperties )
  {
    for ( auto& qp : qps ) {
      qp.material = std::unique_ptr< MarmotMaterialHypoElasticInterface >(
        dynamic_cast< MarmotMaterialHypoElasticInterface* >(
          MarmotLibrary::MarmotMaterialHypoElasticInterfaceFactory::createMaterial( materialName,
                                                                                    materialProperties,
                                                                                    nMaterialProperties,
                                                                                    elLabel ) ) );

      if ( !qp.material )
        throw std::invalid_argument(
          MakeString() << __PRETTY_FUNCTION__
                       << ": invalid material assigned; cannot cast to MarmotMaterialHypoElasticInterface!" );

      if constexpr ( nDim == 3 )
        qp.material->setCharacteristicElementLength( std::cbrt( 8 * qp.detJ ) );
      if constexpr ( nDim == 2 )
        qp.material->setCharacteristicElementLength( std::sqrt( 4 * qp.detJ ) );
      if constexpr ( nDim == 1 )
        qp.material->setCharacteristicElementLength( 2 * qp.detJ );
    }
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::assignNodeCoordinates( const double* coordinates )
  {
    new ( &this->coordinates ) Eigen::Map< const CoordinateVector >( coordinates );
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::initializeYourself()
  {
    constexpr int nNodesLoc = nNodes / 2;
    constexpr int nDimLoc   = nDim - 1;

    for ( QuadraturePoint& qp : qps ) {

      Fastor::Tensor< double, nNodesLoc >          N_fastor;
      Fastor::Tensor< double, nDimLoc, nNodesLoc > dNdXi_fastor;

      if constexpr ( nDim == 2 && nNodesLoc == 2 ) {
        auto N_val  = Marmot::FiniteElement::Spatial1D::Bar2::N( qp.xi( 0 ) );
        auto dN_val = Marmot::FiniteElement::Spatial1D::Bar2::dNdXi( qp.xi( 0 ) );
        for ( int i = 0; i < nNodesLoc; ++i ) {
          N_fastor( i )        = N_val( i );
          dNdXi_fastor( 0, i ) = dN_val( i );
        }
      }
      else if constexpr ( nDim == 2 && nNodesLoc == 3 ) {
        auto N_val  = Marmot::FiniteElement::Spatial1D::Bar3::N( qp.xi( 0 ) );
        auto dN_val = Marmot::FiniteElement::Spatial1D::Bar3::dNdXi( qp.xi( 0 ) );
        for ( int i = 0; i < nNodesLoc; ++i ) {
          N_fastor( i )        = N_val( i );
          dNdXi_fastor( 0, i ) = dN_val( i );
        }
      }
      else if constexpr ( nDim == 3 && nNodesLoc == 4 ) {
        Eigen::Vector2d xi2d( qp.xi( 0 ), qp.xi( 1 ) );
        auto            N_val  = Marmot::FiniteElement::Spatial2D::Quad4::N( xi2d );
        auto            dN_val = Marmot::FiniteElement::Spatial2D::Quad4::dNdXi( xi2d );
        for ( int i = 0; i < nNodesLoc; ++i ) {
          N_fastor( i )        = N_val( i );
          dNdXi_fastor( 0, i ) = dN_val( 0, i );
          dNdXi_fastor( 1, i ) = dN_val( 1, i );
        }
      }
      else if constexpr ( nDim == 3 && nNodesLoc == 8 ) {
        Eigen::Vector2d xi2d( qp.xi( 0 ), qp.xi( 1 ) );
        auto            N_val  = Marmot::FiniteElement::Spatial2D::Quad8::N( xi2d );
        auto            dN_val = Marmot::FiniteElement::Spatial2D::Quad8::dNdXi( xi2d );
        for ( int i = 0; i < nNodesLoc; ++i ) {
          N_fastor( i )        = N_val( i );
          dNdXi_fastor( 0, i ) = dN_val( 0, i );
          dNdXi_fastor( 1, i ) = dN_val( 1, i );
        }
      }
      else {
        throw std::invalid_argument( "Unsupported element configuration in InterfaceFiniteElement" );
      }

      Fastor::Tensor< double, nDim, nNodesLoc > X_fastor;
      for ( int i = 0; i < nNodesLoc; ++i ) {
        for ( int d = 0; d < nDim; ++d ) {
          X_fastor( d, i ) = this->coordinates( i * nDim + d );
        }
      }

      // Fastor NumPy-like tensor mapping for Jacobian
      // J = X * dNdXi^T => J_{ik} = X_{ia} dNdXi_{ka}
      // where i=nDim, k=nDimLoc, a=nNodesLoc
      enum { i_fastor, j_fastor, k_fastor, l_fastor, m_fastor, n_fastor, a_fastor, b_fastor };
      Fastor::Tensor< double, nDim, nDimLoc >
        J_fastor = Fastor::einsum< Fastor::Index< i_fastor, a_fastor >,
                                   Fastor::Index< k_fastor, a_fastor > >( X_fastor, dNdXi_fastor );

      // Temporarily transfer J back to Eigen for cross products
      Eigen::Matrix< double, nDim, nDimLoc > J;
      for ( int idx = 0; idx < nDim; ++idx ) {
        for ( int jdx = 0; jdx < nDimLoc; ++jdx ) {
          J( idx, jdx ) = J_fastor( idx, jdx );
        }
      }

      if constexpr ( nDim == 2 && nDimLoc == 1 ) {
        double t_x = J( 0, 0 );
        double t_y = J( 1, 0 );
        qp.detJ    = std::sqrt( t_x * t_x + t_y * t_y );
        if ( qp.detJ > 1e-16 ) {
          qp.normal( 0 ) = -t_y / qp.detJ;
          qp.normal( 1 ) = t_x / qp.detJ;
        }
      }
      else if constexpr ( nDim == 3 && nDimLoc == 2 ) {
        Eigen::Vector3d t1    = J.col( 0 );
        Eigen::Vector3d t2    = J.col( 1 );
        Eigen::Vector3d cross = t1.cross( t2 );
        qp.detJ               = cross.norm();
        if ( qp.detJ > 1e-16 ) {
          qp.normal = cross / qp.detJ;
        }
      }

      Fastor::Tensor< double, nDim > n_tensor;
      for ( int d = 0; d < nDim; ++d )
        n_tensor( d ) = qp.normal( d );

      Fastor::Tensor< double, nDim, nDim > T_fastor;
      T_fastor.eye();
      // T_{ij} = I_{ij} - n_i n_j
      T_fastor = T_fastor -
                 Fastor::einsum< Fastor::Index< i_fastor >, Fastor::Index< j_fastor > >( n_tensor, n_tensor );

      // G = J^T * J => G_{kl} = J_{ik} J_{il}
      Fastor::Tensor< double, nDimLoc, nDimLoc >
                                                 G_fastor     = Fastor::einsum< Fastor::Index< i_fastor, k_fastor >,
                                   Fastor::Index< i_fastor, l_fastor > >( J_fastor, J_fastor );
      Fastor::Tensor< double, nDimLoc, nDimLoc > G_inv_fastor = Fastor::inverse( G_fastor );

      // Calculate spatial gradient completely in Fastor tensor math
      // pseudoinv_T_{il} = J_{ik} G_inv_{kl}
      Fastor::Tensor< double, nDim, nDimLoc >
        pseudoinv_T = Fastor::einsum< Fastor::Index< i_fastor, k_fastor >,
                                      Fastor::Index< k_fastor, l_fastor > >( J_fastor, G_inv_fastor );

      // grad_phi_s_{ia} = pseudoinv_T_{ik} dNdXi_{ka}
      Fastor::Tensor< double, nDim, nNodesLoc >
        grad_phi_s_fastor = Fastor::einsum< Fastor::Index< i_fastor, k_fastor >,
                                            Fastor::Index< k_fastor, a_fastor > >( pseudoinv_T, dNdXi_fastor );

      // grad_A_{ia} = T_{ij} grad_phi_s_{ja}
      Fastor::Tensor< double, nDim, nNodesLoc >
        grad_A_fastor = Fastor::einsum< Fastor::Index< i_fastor, j_fastor >,
                                        Fastor::Index< j_fastor, a_fastor > >( T_fastor, grad_phi_s_fastor );

      qp.N_jump.setZero();
      qp.B_surface.setZero();

      qp.N_local.zeros();
      qp.B_local.zeros();

      for ( int node_a = 0; node_a < nNodesLoc; ++node_a ) {
        for ( int mdof = 0; mdof < nDim; ++mdof ) {
          qp.N_jump( mdof, node_a * nDim + mdof )                 = -N_fastor( node_a );
          qp.N_jump( mdof, ( node_a + nNodesLoc ) * nDim + mdof ) = N_fastor( node_a );
          qp.N_local( mdof, node_a )                              = N_fastor( node_a );
        }

        for ( int mdof = 0; mdof < nDim; ++mdof ) {
          for ( int idim = 0; idim < nDim; ++idim ) {
            for ( int kdim = 0; kdim < nDim; ++kdim ) {
              int    row          = idim * nDim + kdim;
              double val          = 0.0;
              double val_unscaled = 0.0;

              if constexpr ( nDim == 2 ) {
                // Project only on spatial gradient side
                if ( idim == mdof ) {
                  val          = 0.5 * grad_A_fastor( kdim, node_a );
                  val_unscaled = grad_A_fastor( kdim, node_a );
                }
              }
              else {
                // Project on both component side (T_fastor) and spatial side (grad_A_fastor)
                val          = 0.5 * T_fastor( idim, mdof ) * grad_A_fastor( kdim, node_a );
                val_unscaled = T_fastor( idim, mdof ) * grad_A_fastor( kdim, node_a );
              }

              // Convert the results generated by FastorTensors to C++ Matrix formats used by the Marmot solver
              if ( val != 0.0 ) {
                qp.B_surface( row, node_a * nDim + mdof ) += val;
                qp.B_surface( row, ( node_a + nNodesLoc ) * nDim + mdof ) += val;
              }
              qp.B_local( idim, kdim, mdof, node_a ) = val_unscaled;
            }
          }
        }
      }

      if constexpr ( nDim == 3 ) {
        qp.J0xW = qp.weight * qp.detJ;
      }
      if constexpr ( nDim == 2 ) {
        const double& thickness = elementProperties[0];
        qp.J0xW                 = qp.weight * qp.detJ * thickness;
      }
      if constexpr ( nDim == 1 ) {
        const double& crossSection = elementProperties[0];
        qp.J0xW                    = qp.weight * qp.detJ * crossSection;
      }
    }
  }

#include <Fastor/Fastor.h>
  namespace InterfaceElementMatrices {
    using Fastor::Index;
    enum { i_f, j_f, k_f, a_f, b_f, m_f, n_f, p_f };

    template < int nDim, int nNodes >
    Eigen::Matrix< double, nNodes * nDim, nNodes * nDim > assign_K_jumpu_jumpv(
      const Eigen::Matrix< double, nDim, nNodes * nDim >& N_jump,
      const Eigen::Matrix< double, nDim, nDim >&          H_inv_ij )
    {
      constexpr int                     nd       = nNodes * nDim / 2;
      Eigen::Matrix< double, nDim, nd > N_matrix = N_jump.block( 0, nd, nDim, nd );

      Fastor::Tensor< double, nd, nDim >   N_T;
      Fastor::Tensor< double, nDim, nDim > H;
      Fastor::Tensor< double, nDim, nd >   N;

      for ( int r = 0; r < nDim; r++ )
        for ( int c = 0; c < nd; c++ ) {
          N( r, c )   = N_matrix( r, c );
          N_T( c, r ) = N_matrix( r, c );
        }
      for ( int r = 0; r < nDim; r++ )
        for ( int c = 0; c < nDim; c++ )
          H( r, c ) = H_inv_ij( r, c );

      Fastor::Tensor< double, nd, nd >
        f_K = Fastor::einsum< Index< i_f, m_f >, Index< m_f, n_f >, Index< n_f, j_f > >( N_T, H, N );

      Eigen::Matrix< double, nNodes * nDim, nNodes * nDim > K_out;
      for ( int r = 0; r < nd; r++ )
        for ( int c = 0; c < nd; c++ ) {
          double val              = f_K( r, c );
          K_out( r, c )           = val;
          K_out( r, c + nd )      = -val;
          K_out( r + nd, c )      = -val;
          K_out( r + nd, c + nd ) = val;
        }
      return K_out;
    }

    template < int nDim, int nNodes >
    Eigen::Matrix< double, nNodes * nDim, 1 > assign_P_jumpv(
      const Eigen::Matrix< double, nDim, nNodes * nDim >& N_jump,
      const Eigen::Matrix< double, nDim, 1 >&             force )
    {
      constexpr int                     nd       = nNodes * nDim / 2;
      Eigen::Matrix< double, nDim, nd > N_matrix = N_jump.block( 0, nd, nDim, nd );

      Fastor::Tensor< double, nd, nDim > N_T;
      Fastor::Tensor< double, nDim, 1 >  F;

      for ( int r = 0; r < nDim; r++ )
        for ( int c = 0; c < nd; c++ )
          N_T( c, r ) = N_matrix( r, c );
      for ( int r = 0; r < nDim; r++ )
        F( r, 0 ) = force( r, 0 );

      Fastor::Tensor< double, nd, 1 > f_P = Fastor::einsum< Index< i_f, m_f >, Index< m_f, j_f > >( N_T, F );

      Eigen::Matrix< double, nNodes * nDim, 1 > P_out;
      for ( int r = 0; r < nd; r++ ) {
        double val      = f_P( r, 0 );
        P_out( r )      = -val;
        P_out( r + nd ) = val;
      }
      return P_out;
    }

    template < int nDim, int nNodes >
    Eigen::Matrix< double, nNodes * nDim, nNodes * nDim > assign_K_grad_s_u_grad_s_v(
      const Eigen::Matrix< double, nDim * nDim, nNodes * nDim >& B_surface,
      const Eigen::Matrix< double, nDim * nDim, nDim * nDim >&   Z_ijkl )
    {
      constexpr int                   nd       = nNodes * nDim / 2;
      constexpr int                   ds       = nDim * nDim;
      Eigen::Matrix< double, ds, nd > B_matrix = 2.0 * B_surface.block( 0, nd, ds, nd );

      Fastor::Tensor< double, nd, ds > B_T;
      Fastor::Tensor< double, ds, ds > Z;
      Fastor::Tensor< double, ds, nd > B;

      for ( int r = 0; r < ds; r++ )
        for ( int c = 0; c < nd; c++ ) {
          B( r, c )   = B_matrix( r, c );
          B_T( c, r ) = B_matrix( r, c );
        }
      for ( int r = 0; r < ds; r++ )
        for ( int c = 0; c < ds; c++ )
          Z( r, c ) = Z_ijkl( r, c );

      Fastor::Tensor< double, nd, nd >
        f_K = Fastor::einsum< Index< i_f, m_f >, Index< m_f, n_f >, Index< n_f, j_f > >( B_T, Z, B );

      Eigen::Matrix< double, nNodes * nDim, nNodes * nDim > K_out;
      for ( int r = 0; r < nd; r++ )
        for ( int c = 0; c < nd; c++ ) {
          double val              = 0.25 * f_K( r, c );
          K_out( r, c )           = val;
          K_out( r, c + nd )      = val;
          K_out( r + nd, c )      = val;
          K_out( r + nd, c + nd ) = val;
        }
      return K_out;
    }

    template < int nDim, int nNodes >
    Eigen::Matrix< double, nNodes * nDim, 1 > assign_P_grad_s_v(
      const Eigen::Matrix< double, nDim * nDim, nNodes * nDim >& B_surface,
      const Eigen::Matrix< double, nDim * nDim, 1 >&             surface_stress )
    {
      constexpr int                   nd       = nNodes * nDim / 2;
      constexpr int                   ds       = nDim * nDim;
      Eigen::Matrix< double, ds, nd > B_matrix = 2.0 * B_surface.block( 0, nd, ds, nd );

      Fastor::Tensor< double, nd, ds > B_T;
      Fastor::Tensor< double, ds, 1 >  S;

      for ( int r = 0; r < ds; r++ )
        for ( int c = 0; c < nd; c++ )
          B_T( c, r ) = B_matrix( r, c );
      for ( int r = 0; r < ds; r++ )
        S( r, 0 ) = surface_stress( r, 0 );

      Fastor::Tensor< double, nd, 1 > f_P = Fastor::einsum< Index< i_f, m_f >, Index< m_f, j_f > >( B_T, S );

      Eigen::Matrix< double, nNodes * nDim, 1 > P_out;
      for ( int r = 0; r < nd; r++ ) {
        double val      = 0.5 * f_P( r, 0 );
        P_out( r )      = val;
        P_out( r + nd ) = val;
      }
      return P_out;
    }

    template < int nDim, int nNodes >
    Eigen::Matrix< double, nNodes * nDim, nNodes * nDim > assign_K_grad_s_u_jump_v(
      const Eigen::Matrix< double, nDim, nNodes * nDim >&        N_jump,
      const Eigen::Matrix< double, nDim * nDim, nNodes * nDim >& B_surface,
      const Eigen::Matrix< double, nDim, nDim * nDim >&          H_inv_nF )
    {
      constexpr int                            nd       = nNodes * nDim / 2;
      Eigen::Matrix< double, nDim, nd >        N_matrix = N_jump.block( 0, nd, nDim, nd );
      Eigen::Matrix< double, nDim * nDim, nd > B_matrix = 2.0 * B_surface.block( 0, nd, nDim * nDim, nd );

      Fastor::Tensor< double, nDim, nd >         N;
      Fastor::Tensor< double, nDim, nDim, nDim > t;
      Fastor::Tensor< double, nDim, nDim, nd >   B;

      for ( int r = 0; r < nDim; r++ )
        for ( int c = 0; c < nd; c++ )
          N( r, c ) = N_matrix( r, c );
      for ( int m = 0; m < nDim; m++ )
        for ( int n = 0; n < nDim; n++ )
          for ( int p = 0; p < nDim; p++ )
            t( m, n, p ) = H_inv_nF( m, n * nDim + p );
      for ( int m = 0; m < nDim; m++ )
        for ( int n = 0; n < nDim; n++ )
          for ( int d = 0; d < nd; d++ )
            B( m, n, d ) = B_matrix( m * nDim + n, d );

      Fastor::Tensor< double, nd, nd >
        f_A = Fastor::einsum< Index< i_f, a_f >, Index< i_f, j_f, k_f >, Index< j_f, k_f, b_f > >( N, t, B );

      Eigen::Matrix< double, nNodes * nDim, nNodes * nDim > K_out;
      for ( int r = 0; r < nd; r++ )
        for ( int c = 0; c < nd; c++ ) {
          double val              = f_A( r, c );
          K_out( r, c )           = -0.5 * val;
          K_out( r, c + nd )      = -0.5 * val;
          K_out( r + nd, c )      = 0.5 * val;
          K_out( r + nd, c + nd ) = 0.5 * val;
        }
      return K_out;
    }

    template < int nDim, int nNodes >
    Eigen::Matrix< double, nNodes * nDim, nNodes * nDim > assign_K_jump_u_grad_s_v(
      const Eigen::Matrix< double, nDim, nNodes * nDim >&        N_jump,
      const Eigen::Matrix< double, nDim * nDim, nNodes * nDim >& B_surface,
      const Eigen::Matrix< double, nDim, nDim * nDim >&          H_inv_nF )
    {
      constexpr int                            nd       = nNodes * nDim / 2;
      Eigen::Matrix< double, nDim, nd >        N_matrix = N_jump.block( 0, nd, nDim, nd );
      Eigen::Matrix< double, nDim * nDim, nd > B_matrix = 2.0 * B_surface.block( 0, nd, nDim * nDim, nd );

      Fastor::Tensor< double, nDim, nd >         N;
      Fastor::Tensor< double, nDim, nDim, nDim > t;
      Fastor::Tensor< double, nDim, nDim, nd >   B;

      for ( int r = 0; r < nDim; r++ )
        for ( int c = 0; c < nd; c++ )
          N( r, c ) = N_matrix( r, c );
      for ( int m = 0; m < nDim; m++ )
        for ( int n = 0; n < nDim; n++ )
          for ( int p = 0; p < nDim; p++ )
            t( m, n, p ) = H_inv_nF( m, n * nDim + p );
      for ( int m = 0; m < nDim; m++ )
        for ( int n = 0; n < nDim; n++ )
          for ( int d = 0; d < nd; d++ )
            B( m, n, d ) = B_matrix( m * nDim + n, d );

      Fastor::Tensor< double, nd, nd >
        f_A = Fastor::einsum< Index< i_f, a_f >, Index< i_f, j_f, k_f >, Index< j_f, k_f, b_f > >( N, t, B );

      Eigen::Matrix< double, nNodes * nDim, nNodes * nDim > K_out;
      for ( int r = 0; r < nd; r++ )
        for ( int c = 0; c < nd; c++ ) {
          double val_AT           = f_A( c, r );
          K_out( r, c )           = -0.5 * val_AT;
          K_out( r, c + nd )      = 0.5 * val_AT;
          K_out( r + nd, c )      = -0.5 * val_AT;
          K_out( r + nd, c + nd ) = 0.5 * val_AT;
        }
      return K_out;
    }
  } // namespace InterfaceElementMatrices

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::computeYourself( const double* QTotal_,
                                                                const double* dQ_,
                                                                double*       Pe_,
                                                                double*       Ke_,
                                                                const double* time,
                                                                double        dT,
                                                                double&       pNewDT )
  {
    using namespace Marmot;
    using namespace ContinuumMechanics::VoigtNotation;

    Map< const RhsSized > QTotal( QTotal_ );
    Map< const RhsSized > dQ( dQ_ );
    Map< KeSizedMatrix >  Ke( Ke_ );
    Map< RhsSized >       Pe( Pe_ );

    for ( QuadraturePoint& qp : qps ) {

      const NMatrixSized&  N_jump    = qp.N_jump;
      const BSurfaceSized& B_surface = qp.B_surface;

      Eigen::Matrix< double, nDim, 1 >        force          = Eigen::Matrix< double, nDim, 1 >::Zero();
      Eigen::Matrix< double, nDim * nDim, 1 > surface_stress = Eigen::Matrix< double, nDim * nDim, 1 >::Zero();

      Eigen::Matrix< double, nDim, nDim > Q_ij = Eigen::Matrix< double, nDim, nDim >::Zero();
      Eigen::Matrix< double, nDim * nDim, nDim* nDim >
                                                Z_ijkl = Eigen::Matrix< double, nDim * nDim, nDim * nDim >::Zero();
      Eigen::Matrix< double, nDim, nDim* nDim > H_ijk  = Eigen::Matrix< double, nDim, nDim * nDim >::Zero();
      Eigen::Matrix< double, nDim * nDim, nDim* nDim >
        Y_ijkl = Eigen::Matrix< double, nDim * nDim, nDim * nDim >::Zero();

      Eigen::Matrix< double, nDim * 2, 1 >        dU_GPs;
      Eigen::Matrix< double, nDim * nDim * 2, 1 > dSurface_strain_GPs;

      constexpr int                                   half_size = nNodes * nDim / 2;
      Eigen::Matrix< double, nDim, half_size >        N_matrix  = qp.N_jump.block( 0, half_size, nDim, half_size );
      Eigen::Matrix< double, nDim * nDim, half_size > B_matrix  = 2.0 * qp.B_surface.block( 0,
                                                                                           half_size,
                                                                                           nDim * nDim,
                                                                                           half_size );

      dU_GPs.segment( 0, nDim )    = N_matrix * dQ.segment( half_size, half_size );
      dU_GPs.segment( nDim, nDim ) = N_matrix * dQ.segment( 0, half_size );

      dSurface_strain_GPs.segment( 0, nDim * nDim )           = B_matrix * dQ.segment( half_size, half_size );
      dSurface_strain_GPs.segment( nDim * nDim, nDim * nDim ) = B_matrix * dQ.segment( 0, half_size );

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

      if ( pNewDT < 1.0 )
        return;

      // Integration weight tracking J0xW = weight * detJ * thickness
      Ke += InterfaceElementMatrices::assign_K_jumpu_jumpv< nDim, nNodes >( qp.N_jump, Q_ij ) * qp.J0xW;
      Pe -= InterfaceElementMatrices::assign_P_jumpv< nDim, nNodes >( qp.N_jump, force ) * qp.J0xW;

      Ke += InterfaceElementMatrices::assign_K_grad_s_u_grad_s_v< nDim, nNodes >( qp.B_surface, Z_ijkl ) * qp.J0xW;
      Pe -= InterfaceElementMatrices::assign_P_grad_s_v< nDim, nNodes >( qp.B_surface, surface_stress ) * qp.J0xW;

      Ke += InterfaceElementMatrices::assign_K_grad_s_u_grad_s_v< nDim, nNodes >( qp.B_surface, Y_ijkl ) * qp.J0xW;

      Ke += InterfaceElementMatrices::assign_K_jump_u_grad_s_v< nDim, nNodes >( qp.N_jump, qp.B_surface, H_ijk ) *
            qp.J0xW;
      Ke += InterfaceElementMatrices::assign_K_grad_s_u_jump_v< nDim, nNodes >( qp.N_jump, qp.B_surface, H_ijk ) *
            qp.J0xW;
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
    case MarmotElement::GeostaticStress: {
      if ( nDim >= 2 )
        for ( QuadraturePoint& qp : qps ) {
          throw std::invalid_argument( "Geostatic stress not supported for InterfaceFiniteElement" );
          Eigen::VectorXd coordAtGauss( 2 ); // dummy

          const double sigY1 = values[0];
          const double sigY2 = values[2];
          const double y1    = values[1];
          const double y2    = values[3];

          using namespace Math;
          qp.managedStateVars->stress( 1 ) = linearInterpolation( coordAtGauss[1], y1, y2, sigY1, sigY2 ); // sigma_y
          qp.managedStateVars->stress( 0 ) = values[4] * qp.managedStateVars->stress( 1 );                 // sigma_x
          qp.managedStateVars->stress( 2 ) = values[5] * qp.managedStateVars->stress( 1 );
        }
      break;
    }
    case MarmotElement::MarmotMaterialStateVars: {
      throw std::invalid_argument( "Please use initializeStateVars directly on material" );
    }
    default: throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid initial condition" );
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
    Map< RhsSized > fU( P );
    // Interface elements typically don't apply boundary conditions directly through faces in Marmot standard structure.
    throw std::invalid_argument(
      "computeDistributedLoad not strictly supported natively for Cohesive/Interface Elements" );
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::computeBodyForce( double*       P_,
                                                                 double*       K,
                                                                 const double* load,
                                                                 const double* QTotal,
                                                                 const double* time,
                                                                 double        dT )
  {
    Map< RhsSized >                              Pe( P_ );
    const Map< const Matrix< double, nDim, 1 > > f( load );

    // Average body force across upper and lower surface?
    for ( const auto& qp : qps ) {
      // Map body force equally to top and bottom nodes? Or only on base if defined properly.
      // Usually, interface elements don't get body forces natively distinct from the solid domains.
      // But we map using standard shape functions:
      const auto N_jump = qp.N_jump;
      // Note: N_jump is (nDim x sizeLoadVector), typically f is nDim x 1.
      Pe += N_jump.transpose() * f * qp.J0xW; // Simplification for cohesive element
    }
  }

  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::computeConsistentInertia( double* M )
  {
    Map< KeSizedMatrix > Me( M );
    Me.setZero();

    for ( const auto& qp : qps ) {
      const auto   N_  = qp.N_jump; // Use the appropriately sized matrix
      const double rho = qp.material->getDensity();
      Me += N_.transpose() * N_ * qp.detJ * qp.weight * rho;
    }
  }
  template < int nDim, int nNodes >
  void InterfaceFiniteElement< nDim, nNodes >::computeLumpedInertia( double* M )
  {
    Map< RhsSized > Me( M );
    Me.setZero();

    KeSizedMatrix CMM;
    CMM.setZero();
    computeConsistentInertia( CMM.data() );

    Me = CMM.rowwise().sum();
  }

  template < int nDim, int nNodes >
  std::vector< double > InterfaceFiniteElement< nDim, nNodes >::getCoordinatesAtCenter()
  {
    std::vector< double > coords( nDim );
    Eigen::Map< XiSized > coordsMap( &coords[0] );
    const auto            centerXi = XiSized::Zero();
    throw std::invalid_argument( "getCoordinatesAtCenter not implemented" );
    return coords;
  }

  template < int nDim, int nNodes >
  std::vector< std::vector< double > > InterfaceFiniteElement< nDim, nNodes >::getCoordinatesAtQuadraturePoints()
  {
    std::vector< std::vector< double > > listedCoords;
    std::vector< double >                coords( nDim );
    Eigen::Map< XiSized >                coordsMap( &coords[0] );

    for ( const auto& qp : qps ) {
      throw std::invalid_argument( "getCoordinatesAtQuadraturePoints not implemented" );
      listedCoords.push_back( coords );
    }

    return listedCoords;
  }

  template < int nDim, int nNodes >
  int InterfaceFiniteElement< nDim, nNodes >::getNumberOfQuadraturePoints()
  {
    return qps.size();
  }

  template class InterfaceFiniteElement< 2, 4 >;
  template class InterfaceFiniteElement< 2, 6 >;
  template class InterfaceFiniteElement< 3, 8 >;
  template class InterfaceFiniteElement< 3, 16 >;

} // namespace Marmot::Elements
