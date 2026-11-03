#include "Marmot/MarmotEquilibratedXFacePropsInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotEquilibratedXInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotTesting.h"

#include <Eigen/Dense>

#include <cmath>
#include <functional>
#include <stdexcept>
#include <vector>

using namespace Marmot::Testing;

namespace {

  using Base  = MarmotEquilibratedXInterfaceMaterialHypoElastic;
  using Faces = MarmotEquilibratedXFacePropsInterfaceMaterialHypoElastic;

  const double normalZ[3]      = { 0., 0., 1. };
  const double noSeparation[3] = { 0., 0., 0. };

  /// Generalized force of one update with the jump w and vanishing surface gradients.
  template < class Material >
  Eigen::Vector3d forceOfJump( Material& material, const Eigen::Vector3d& w )
  {
    std::vector< double > stateVars( material.getNumberOfRequiredStateVars(), 0.0 );
    material.initializeYourself( stateVars.data(), static_cast< int >( stateVars.size() ) );

    Base::State state{ Base::Tensor3d( 0.0 ), Base::Tensor33d( 0.0 ), Base::Tensor33d( 0.0 ), stateVars.data() };
    Base::Tangents tangents{ Base::Tensor33d( 0.0 ),   Base::Tensor333d( 0.0 ),  Base::Tensor333d( 0.0 ),
                             Base::Tensor333d( 0.0 ),  Base::Tensor3333d( 0.0 ), Base::Tensor3333d( 0.0 ),
                             Base::Tensor333d( 0.0 ),  Base::Tensor3333d( 0.0 ), Base::Tensor3333d( 0.0 ) };

    double dU[6]   = { w( 0 ), w( 1 ), w( 2 ), 0., 0., 0. };
    double dA[18]  = { 0. };
    Base::Deformation deformation{ Base::Tensor6d( dU ),
                                   Base::Tensor18d( dA ),
                                   Base::Tensor3d( normalZ ),
                                   Base::Tensor3d( noSeparation ) };

    material.computeStress( state, tangents, deformation, Base::TimeIncrement{ 0.0, 1.0 } );
    return Eigen::Map< const Eigen::Vector3d >( state.generalizedForce.data() );
  }

  double constrainedModulus( double E, double nu )
  {
    return E * ( 1. - nu ) / ( ( 1. + nu ) * ( 1. - 2. * nu ) );
  }

  double shearModulus( double E, double nu )
  {
    return E / ( 2. * ( 1. + nu ) );
  }

  void testEqualFacesGiveTheBaseMaterial()
  {
    const std::vector< double > base  = { 2e5, 0.3, 0.1 };
    const std::vector< double > faces = { 2e5, 0.3, 0.1, 2e5, 0.3 };
    Base                        reference( "LINEARELASTIC", base.data(), 3, 1 );
    Faces                       biphasic( "LINEARELASTIC", faces.data(), 5, 1 );

    const Eigen::Vector3d w( 1e-3, -2e-3, 3e-3 );
    const Eigen::Vector3d f0 = forceOfJump( reference, w );
    const Eigen::Vector3d f1 = forceOfJump( biphasic, w );
    throwExceptionOnFailure( ( f1 - f0 ).norm() / f0.norm() < 1e-12,
                             "Equal faces have to reproduce the base material." );
  }

  /// The two halves are in series, so the layer carries the harmonic mean of the two moduli.
  void testSeriesLaw()
  {
    const double Et = 4e5, Eb = 4e3, nu = 0.3;
    const std::vector< double > base  = { 2e5, nu, 0.1 };
    const std::vector< double > faces = { Et, nu, 0.1, Eb, nu };
    Base                        reference( "LINEARELASTIC", base.data(), 3, 1 );
    Faces                       biphasic( "LINEARELASTIC", faces.data(), 5, 1 );

    const Eigen::Vector3d normalJump( 0., 0., 1e-3 );
    const Eigen::Vector3d shearJump( 1e-3, 0., 0. );

    const double Mmean = 2. / ( 1. / constrainedModulus( Et, nu ) + 1. / constrainedModulus( Eb, nu ) );
    const double Gmean = 2. / ( 1. / shearModulus( Et, nu ) + 1. / shearModulus( Eb, nu ) );

    const double normalRatio = forceOfJump( biphasic, normalJump )( 2 ) / forceOfJump( reference, normalJump )( 2 );
    const double shearRatio  = forceOfJump( biphasic, shearJump )( 0 ) / forceOfJump( reference, shearJump )( 0 );

    throwExceptionOnFailure( std::abs( normalRatio - Mmean / constrainedModulus( 2e5, nu ) ) < 1e-9 * normalRatio,
                             "The normal stiffness of a biphasic layer is the harmonic mean of the two moduli." );
    throwExceptionOnFailure( std::abs( shearRatio - Gmean / shearModulus( 2e5, nu ) ) < 1e-9 * shearRatio,
                             "The shear stiffness of a biphasic layer is the harmonic mean of the two moduli." );
  }

  void testRejectsTooFewProperties()
  {
    const std::vector< double > faces = { 2e5, 0.3, 0.1, 2e5 };
    bool                        thrown = false;
    try {
      Faces material( "LINEARELASTIC", faces.data(), 4, 1 );
    }
    catch ( const std::invalid_argument& ) {
      thrown = true;
    }
    throwExceptionOnFailure( thrown, "Fewer than five properties must be rejected." );
  }

} // namespace

int main()
{
  std::vector< std::function< void() > > tests = { testEqualFacesGiveTheBaseMaterial,
                                                   testSeriesLaw,
                                                   testRejectsTooFewProperties };
  executeTestsAndCollectExceptions( tests );
  return 0;
}
