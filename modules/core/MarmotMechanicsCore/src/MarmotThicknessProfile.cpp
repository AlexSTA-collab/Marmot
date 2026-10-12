#include "Marmot/MarmotThicknessProfile.h"

#include "Marmot/MarmotUtils.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace Marmot::Materials {

  std::pair< std::vector< double >, std::vector< double > > gaussLegendre( int n )
  {
    std::vector< double > nodes( n ), weights( n );
    for ( int i = 0; i < n; i++ ) {
      // Newton's method on the Legendre polynomial, started from the usual estimate of the root
      double x  = std::cos( M_PI * ( i + 0.75 ) / ( n + 0.5 ) );
      double dP = 1.0;
      for ( int iteration = 0; iteration < 100; iteration++ ) {
        double p0 = 1.0, p1 = x;
        for ( int k = 2; k <= n; k++ ) {
          const double p2 = ( ( 2. * k - 1. ) * x * p1 - ( k - 1. ) * p0 ) / k;
          p0              = p1;
          p1              = p2;
        }
        const double pm = ( n == 1 ) ? 1.0 : p0;
        dP              = n * ( x * p1 - pm ) / ( x * x - 1. );
        const double dx = p1 / dP;
        x -= dx;
        if ( std::abs( dx ) < 1e-15 )
          break;
      }
      nodes[i]   = x;
      weights[i] = 2. / ( ( 1. - x * x ) * dP * dP );
    }
    return { nodes, weights };
  }

  ThicknessProfile ThicknessProfile::fromMaterialProperties( const double* matProperties, int nMaterialProperties )
  {
    ThicknessProfile profile;
    profile.thickness = checkedMaterialProperty( matProperties, nMaterialProperties, 0 );

    if ( profile.thickness <= 0. )
      throw std::invalid_argument( "ThicknessProfile: the thickness must be positive" );

    const int mode   = static_cast< int >( checkedMaterialProperty( matProperties, nMaterialProperties, 1 ) );
    const int nGauss = static_cast< int >( checkedMaterialProperty( matProperties, nMaterialProperties, 2 ) );
    const int K      = static_cast< int >( checkedMaterialProperty( matProperties, nMaterialProperties, 3 ) );
    const int m      = static_cast< int >( checkedMaterialProperty( matProperties, nMaterialProperties, 4 ) );

    if ( mode != 0 && mode != 1 )
      throw std::invalid_argument( "ThicknessProfile: the mode must be 0 (layers) or 1 (interpolated stations)" );
    if ( nGauss < 1 || m < 1 || K < ( mode == 0 ? 1 : 2 ) )
      throw std::invalid_argument( "ThicknessProfile: nGauss and m must be positive, and K at least one layer or two "
                                   "stations" );
    if ( nMaterialProperties != 5 + K * ( 1 + m ) )
      throw std::invalid_argument( "ThicknessProfile: " + std::to_string( 5 + K * ( 1 + m ) ) +
                                   " material properties are required for K = " + std::to_string( K ) +
                                   " and m = " + std::to_string( m ) + ", but " +
                                   std::to_string( nMaterialProperties ) + " are given" );

    const auto position = [&]( int k ) { return matProperties[5 + k * ( 1 + m )]; };
    const auto property = [&]( int k, int i ) { return matProperties[5 + k * ( 1 + m ) + 1 + i]; };

    // segments [a, b] of the thickness, with the stations at both ends (equal for a layer)
    struct Segment {
      double a, b;
      int    kA, kB;
    };
    std::vector< Segment > segments;

    if ( mode == 0 ) {
      double sum = 0.;
      for ( int k = 0; k < K; k++ ) {
        if ( position( k ) <= 0. )
          throw std::invalid_argument( "ThicknessProfile: the layer thicknesses must be positive" );
        sum += position( k );
      }
      double a = 0.;
      for ( int k = 0; k < K; k++ ) {
        const double b = ( k == K - 1 ) ? 1.0 : a + position( k ) / sum;
        segments.push_back( { a, b, k, k } );
        a = b;
      }
    }
    else {
      if ( position( 0 ) != 0. || position( K - 1 ) != 1. )
        throw std::invalid_argument( "ThicknessProfile: the first station must be at 0 and the last at 1" );
      for ( int k = 0; k < K - 1; k++ ) {
        if ( position( k + 1 ) <= position( k ) )
          throw std::invalid_argument( "ThicknessProfile: the stations must be increasing" );
        segments.push_back( { position( k ), position( k + 1 ), k, k + 1 } );
      }
    }

    const auto [nodes, weights] = gaussLegendre( nGauss );

    for ( const auto& segment : segments ) {
      for ( int g = 0; g < nGauss; g++ ) {
        Point point;
        point.weight   = 0.5 * weights[g] * ( segment.b - segment.a );
        point.position = segment.a + 0.5 * ( nodes[g] + 1. ) * ( segment.b - segment.a );

        const double s = ( segment.kA == segment.kB ) ? 0. : ( point.position - segment.a ) / ( segment.b - segment.a );
        point.properties.resize( m );
        for ( int i = 0; i < m; i++ )
          point.properties[i] = ( 1. - s ) * property( segment.kA, i ) + s * property( segment.kB, i );

        profile.points.push_back( std::move( point ) );
      }
    }
    return profile;
  }

} // namespace Marmot::Materials
