#include "Marmot/MarmotTwoSidedInterfaceGeometry.h"

#include <stdexcept>

namespace Marmot::Materials::TwoSidedInterface {

  namespace {
    /// tolerance of the length of the connector
    constexpr double geometryTolerance = 1.0e-12;
  } // namespace

  InterfaceGeometry evaluateInterfaceGeometry( const Eigen::Vector3d& normal,
                                               const Eigen::Vector3d& separationVector,
                                               double                 thickness )
  {
    if ( separationVector.norm() <= geometryTolerance )
      return { thickness, Eigen::Vector3d::Zero() };

    const double normalSeparation = separationVector.dot( normal );
    if ( normalSeparation <= geometryTolerance )
      throw std::invalid_argument( "The top-bottom connector of the faces must have a positive normal component." );

    return { normalSeparation, separationVector - normalSeparation * normal };
  }

  Matrix9x3 jumpToGradient( const Eigen::Vector3d& normal, double normalSeparation )
  {
    Matrix9x3 map = Matrix9x3::Zero();
    for ( int k = 0; k < 3; k++ )
      for ( int l = 0; l < 3; l++ )
        map( flatIndex( k, l ), k ) = normal( l ) / normalSeparation;
    return map;
  }

  Matrix9x9 surfaceToGradient( const Eigen::Vector3d& normal,
                               const Eigen::Vector3d& tangentialSeparation,
                               double                 normalSeparation )
  {
    Matrix9x9 map = Matrix9x9::Zero();
    for ( int k = 0; k < 3; k++ )
      for ( int l = 0; l < 3; l++ )
        for ( int b = 0; b < 3; b++ )
          map( flatIndex( k, l ), flatIndex( k, b ) ) = ( l == b ? 1.0 : 0.0 ) -
                                                        tangentialSeparation( b ) * normal( l ) / normalSeparation;
    return map;
  }

} // namespace Marmot::Materials::TwoSidedInterface
