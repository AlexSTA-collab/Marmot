#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotTypedefs.h"

using namespace Eigen;
namespace Marmot {
  namespace FiniteElement {
    ElementShapes getElementShapeByMetric( int nDim, int nNodes )
    {

      switch ( nDim ) {
      case ( 1 ): {
        switch ( nNodes ) {
        case 2: return ElementShapes::Bar2;
        default: throw std::invalid_argument( "Invalid number of nodes for nDim=1" );
        }
      }

      case ( 2 ): {
        switch ( nNodes ) {
        case 4: return ElementShapes::Quad4;
        case 8: return ElementShapes::Quad8;
        default: throw std::invalid_argument( "Invalid number of nodes for nDim=2" );
        }
      }
      case ( 3 ): {
        switch ( nNodes ) {
        case 4: return ElementShapes::Tetra4;
        case 10: return ElementShapes::Tetra10;
        case 8: return ElementShapes::Hexa8;
        case 20: return ElementShapes::Hexa20;
        default: throw std::invalid_argument( "Invalid number of nodes for nDim=3" );
        }
      }
      default: throw std::invalid_argument( "Invalid dimension specified" );
      }
    }
    // create a large N matrix dependent on the degrees of freedom per node
    MatrixXd NB( const VectorXd& N, const int nDoFPerNode )
    {

      MatrixXd N_( nDoFPerNode, N.size() * nDoFPerNode );
      N_ = MatrixXd::Zero( nDoFPerNode, N.size() * nDoFPerNode );

      for ( int i = 0; i < N.size(); i++ ) {
        for ( int j = 0; j < nDoFPerNode; j++ ) {
          N_( j, nDoFPerNode * i + j ) = N( i );
        }
      }

      return N_;
    }

    MatrixXd Jacobian( const MatrixXd& dNdXi, const VectorXd& coordinates )
    {

      /* Dynamic version of jacobian, and not necessarily square!
       *
       * Notation:
       *
       * /                          \
       * | x1,xi1,  x1,xi2,  x1,xi3 |
       * |                          |
       * | x2,xi1,  x2,xi2,  x2,xi3 |
       * |                          |
       * | x3,xi1,  x3,xi2,  x3,xi3 |
       * \                          /
       *
       * */

      const int nDimXi = dNdXi.rows();
      const int nNodes = dNdXi.cols();
      const int nDimX  = coordinates.size() / nNodes;

      MatrixXd J_ = MatrixXd::Zero( nDimX, nDimXi );

      for ( int i = 0; i < nDimX; i++ )      // loop over global dimensions
        for ( int j = 0; j < nDimXi; j++ )   // loop over local dimensions
          for ( int k = 0; k < nNodes; k++ ) // Loop over nodes
            J_( i, j ) += dNdXi( j, k ) * coordinates( i + k * nDimX );
      return J_;
    }

    VectorXi expandNodeIndicesToCoordinateIndices( const VectorXi& nodeIndices, int nDim )
    {
      VectorXi indices( nDim * nodeIndices.size() );

      for ( int i = 0; i < nodeIndices.size(); i++ )
        for ( int k = 0; k < nDim; k++ )
          indices( ( i * nDim ) + k ) = ( nodeIndices( i ) * nDim ) + k;

      return indices;
    }
  } // namespace FiniteElement

  namespace FiniteElement::Quadrature {
    const std::vector< QuadraturePointInfo >& getGaussPointInfo( Marmot::FiniteElement::ElementShapes shape,
                                                                 IntegrationTypes                     integrationType )
    {
      using Marmot::FiniteElement::ElementShapes;
      switch ( shape ) {
      case ( ElementShapes::Bar2 ): {
        if ( integrationType == IntegrationTypes::FullIntegration )
          return Spatial1D::gaussPointList2;
        else
          return Spatial1D::gaussPointList1;
      }

      case ( ElementShapes::Bar3 ): {
        if ( integrationType == IntegrationTypes::FullIntegration )
          return Spatial1D::gaussPointList3;
        else
          return Spatial1D::gaussPointList2;
      }

      case ( ElementShapes::Quad4 ): {
        if ( integrationType == IntegrationTypes::FullIntegration )
          return Spatial2D::gaussPointList2x2;
        else
          return Spatial2D::gaussPointList1x1;
      }

      case ( ElementShapes::Quad8 ): {
        if ( integrationType == IntegrationTypes::FullIntegration )
          return Spatial2D::gaussPointList3x3;
        else
          return Spatial2D::gaussPointList2x2;
      }

      case ( ElementShapes::Tetra4 ): {
        return Spatial3D::gaussPointListTetra4;
      }
      case ( ElementShapes::Tetra10 ): {
        return Spatial3D::gaussPointListTetra10;
      }

      case ( ElementShapes::Hexa8 ): {
        if ( integrationType == IntegrationTypes::FullIntegration )
          return Spatial3D::gaussPointList2x2x2;
        else
          return Spatial3D::gaussPointList1x1x1;
      }

      case ( ElementShapes::Hexa20 ): {
        if ( integrationType == IntegrationTypes::FullIntegration )
          return Spatial3D::gaussPointList3x3x3;
        else
          return Spatial3D::gaussPointList2x2x2;
      }

      default: {
        throw std::invalid_argument( "Invalid shape/integrationType combination" );
      }
      }
    }

    int getNumGaussPoints( Marmot::FiniteElement::ElementShapes shape, IntegrationTypes integrationType )
    {
      return getGaussPointInfo( shape, integrationType ).size();
    }

  } // end of namespace FiniteElement::Quadrature
} // end of namespace Marmot
