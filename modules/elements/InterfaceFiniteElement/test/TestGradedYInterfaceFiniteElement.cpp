#include "Marmot/InterfaceFiniteElement.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/YInterfaceFiniteElement.h"

#include <array>
#include <cmath>
#include <functional>
#include <memory>
#include <string>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::Elements;

namespace {

  constexpr double thickness = 0.1;

  /// Material properties of the graded layer: [h, mode, nGauss, K, m, (position, bulk properties) x K, flag].
  std::vector< double > gradedCard( const std::vector< std::pair< double, std::vector< double > > >& layers,
                                    int                                                              mode   = 0,
                                    int                                                              nGauss = 2,
                                    int tangentialGradient                                                  = 1 )
  {
    std::vector< double > properties = { thickness,
                                         double( mode ),
                                         double( nGauss ),
                                         double( layers.size() ),
                                         double( layers[0].second.size() ) };
    for ( const auto& [position, bulk] : layers ) {
      properties.push_back( position );
      properties.insert( properties.end(), bulk.begin(), bulk.end() );
    }
    properties.push_back( double( tangentialGradient ) );
    return properties;
  }

  /// Coordinates of an inclined, non-planar-looking interface: first the nodes of one side, then the other side.
  const std::array< double, 24 > coordinates3d = { -0.5, -0.5, 0.088163, 0.5,  -0.5, -0.088163, 0.5, 0.5,  -0.088163,
                                                   -0.5, 0.5,  0.088163, -0.5, -0.5, 0.188163,  0.5, -0.5, 0.011837,
                                                   0.5,  0.5,  0.011837, -0.5, 0.5,  0.188163 };

  /// the top nodes lie at the gap of 0.1 above the bottom nodes
  const std::array< double, 24 > flatCoordinates3d = { -0.5, -0.5, 0.0, 0.5, -0.5, 0.0, 0.5, 0.5, 0.0, -0.5, 0.5, 0.0,
                                                       -0.5, -0.5, 0.1, 0.5, -0.5, 0.1, 0.5, 0.5, 0.1, -0.5, 0.5, 0.1 };

  /// Response of one element evaluation: the residual, and the tangent in column-major order.
  struct Evaluation {
    std::vector< double > Pe;
    std::vector< double > Ke;
  };

  template < class Element >
  std::unique_ptr< Element > makeElement( const std::string&           materialName,
                                          const std::vector< double >& materialProperties,
                                          const double*                coordinates,
                                          std::vector< double >&       stateVars )
  {
    auto element = std::make_unique< Element >( 1,
                                                FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                                Element::SectionType::Interface );
    element->assignNodeCoordinates( coordinates );

    static std::array< double, 1 > elementProperties = { 1.0 };
    ElementProperties              properties( elementProperties.data(), 1 );
    element->assignProperty( properties );

    static std::vector< std::vector< double > > keepAlive;
    keepAlive.push_back( materialProperties );
    element->assignMaterial( materialName, keepAlive.back().data(), keepAlive.back().size() );

    stateVars.assign( element->getNumberOfRequiredStateVars(), 0.0 );
    element->assignStateVars( stateVars.data(), static_cast< int >( stateVars.size() ) );
    element->initializeYourself();
    element->setInitialConditions( MarmotElement::MarmotMaterialInitialization, nullptr );
    return element;
  }

  /// One increment from the initial state; dU is the increment of the nodal displacements.
  template < class Element, int nDof >
  Evaluation evaluateIncrement( const std::string&           materialName,
                                const std::vector< double >& materialProperties,
                                const double*                coordinates,
                                const std::vector< double >& dU )
  {
    std::vector< double > stateVars;
    auto                  element = makeElement< Element >( materialName, materialProperties, coordinates, stateVars );

    Evaluation evaluation;
    evaluation.Pe.assign( nDof, 0.0 );
    evaluation.Ke.assign( nDof * nDof, 0.0 );
    const std::vector< double > U( nDof, 0.0 );

    element->computeKernels( U.data(), dU.data(), evaluation.Pe.data(), evaluation.Ke.data(), 0.0, 1.0 );
    return evaluation;
  }

  std::vector< double > arbitraryIncrement( int nDof, double scale )
  {
    std::vector< double > dU( nDof );
    for ( int i = 0; i < nDof; i++ )
      dU[i] = scale * std::sin( 0.9 * i + 0.3 );
    return dU;
  }

  /// With two points at the faces and without the tangential gradient, the element is the equilibrated element.
  template < int nDim, int nNodes >
  void checkTwoFacesAreEquilibratedElement( const std::string&           materialName,
                                            const std::vector< double >& bulk,
                                            const double*                coordinates,
                                            double                       scale )
  {
    constexpr int nDof = nDim * nNodes;

    std::vector< double > equilibrated = { bulk[0], bulk[1], thickness };
    equilibrated.insert( equilibrated.end(), bulk.begin() + 2, bulk.end() );
    const auto graded = gradedCard( { { 0.0, bulk }, { 1.0, bulk } }, 2, 1, 0 );

    const auto dU = arbitraryIncrement( nDof, scale );

    const auto expected = evaluateIncrement< YInterfaceFiniteElement< nDim, nNodes >, nDof >( materialName,
                                                                                              equilibrated,
                                                                                              coordinates,
                                                                                              dU );
    const auto actual   = evaluateIncrement< GradedYInterfaceFiniteElement< nDim, nNodes >, nDof >( materialName,
                                                                                                  graded,
                                                                                                  coordinates,
                                                                                                  dU );

    const double
      scaleOfKe = Eigen::Map< const Eigen::VectorXd >( expected.Ke.data(), nDof * nDof ).cwiseAbs().maxCoeff();
    throwExceptionOnFailure( checkIfEqual< double >( Eigen::Map< const Eigen::VectorXd >( actual.Pe.data(), nDof ),
                                                     Eigen::Map< const Eigen::VectorXd >( expected.Pe.data(), nDof ),
                                                     1e-8 *
                                                       ( 1. +
                                                         Eigen::Map< const Eigen::VectorXd >( expected.Pe.data(), nDof )
                                                           .norm() ) ),
                             "The residual of the two-faced graded element is not the one of the equilibrated "
                             "element." );
    throwExceptionOnFailure( checkIfEqual< double >( Eigen::Map< const Eigen::VectorXd >( actual.Ke.data(),
                                                                                          nDof * nDof ),
                                                     Eigen::Map< const Eigen::VectorXd >( expected.Ke.data(),
                                                                                          nDof * nDof ),
                                                     1e-8 * scaleOfKe ),
                             "The tangent of the two-faced graded element is not the one of the equilibrated "
                             "element." );
  }

  void TestTwoFacesAreEquilibratedElement3D()
  {
    const std::vector< double > elastic  = { 4000., 0.3 };
    const std::vector< double > vonMises = { 1e4, 0.3, 40., 10., 0., 1., 2400. };
    for ( const double* coordinates : { flatCoordinates3d.data(), coordinates3d.data() } ) {
      checkTwoFacesAreEquilibratedElement< 3, 8 >( "LINEARELASTIC", elastic, coordinates, 1e-3 );
      checkTwoFacesAreEquilibratedElement< 3, 8 >( "VONMISES", vonMises, coordinates, 2e-3 );
    }
  }

  void TestTwoFacesAreEquilibratedElement2D()
  {
    const std::array< double, 8 > coordinates2d = { 0., 0., 1., 0.01, 0., 0.11, 1., 0.1 };
    checkTwoFacesAreEquilibratedElement< 2, 4 >( "LINEARELASTIC", { 4000., 0.3 }, coordinates2d.data(), 1e-3 );
    checkTwoFacesAreEquilibratedElement< 2, 4 >( "VONMISES",
                                                 { 1e4, 0.3, 40., 10., 0., 1., 2400. },
                                                 coordinates2d.data(),
                                                 2e-3 );
  }

  /// The tangent of the element is the derivative of its residual, also with the tangential gradient of the layers.
  template < class Element, int nDof >
  void checkTangentIsDerivativeOfResidual( const std::string&           materialName,
                                           const std::vector< double >& materialProperties,
                                           const double*                coordinates,
                                           double                       scale,
                                           double                       tolerance )
  {
    const auto dU   = arbitraryIncrement( nDof, scale );
    const auto base = evaluateIncrement< Element, nDof >( materialName, materialProperties, coordinates, dU );

    const double perturbation = 1e-8;
    double       scaleOfKe    = 0.;
    for ( const double entry : base.Ke )
      scaleOfKe = std::max( scaleOfKe, std::abs( entry ) );

    for ( int j = 0; j < nDof; j++ ) {
      auto dUPlus = dU;
      dUPlus[j] += perturbation;
      const auto perturbed = evaluateIncrement< Element, nDof >( materialName,
                                                                 materialProperties,
                                                                 coordinates,
                                                                 dUPlus );

      for ( int i = 0; i < nDof; i++ ) {
        const double numerical = ( perturbed.Pe[i] - base.Pe[i] ) / perturbation;
        throwExceptionOnFailure( std::abs( numerical - base.Ke[i + j * nDof] ) < tolerance * scaleOfKe,
                                 "The tangent of the graded two-sided element is not the derivative of its residual "
                                 "at (" +
                                   std::to_string( i ) + "," + std::to_string( j ) + ")." );
      }
    }
  }

  void TestTangentOfStackIsConsistent()
  {
    const auto layers = gradedCard( { { 0.3, { 4000., 0.3 } }, { 0.7, { 20000., 0.2 } } } );
    checkTangentIsDerivativeOfResidual< GradedYInterfaceFiniteElement< 3, 8 >, 24 >( "LINEARELASTIC",
                                                                                     layers,
                                                                                     coordinates3d.data(),
                                                                                     1e-3,
                                                                                     1e-5 );

    const std::array< double, 8 > coordinates2d = { 0., 0., 1., 0.01, 0., 0.11, 1., 0.1 };
    checkTangentIsDerivativeOfResidual< GradedYInterfaceFiniteElement< 2, 4 >, 8 >( "LINEARELASTIC",
                                                                                    layers,
                                                                                    coordinates2d.data(),
                                                                                    1e-3,
                                                                                    1e-5 );
  }

  /**
   * The stiffness of a patch of elastic graded elements is symmetric: the test function of the variational form is
   * the derivative of the strain of the layer, which contains the response of the microproblem and the tangential
   * gradient of the correction. This holds with the closure of the tangential gradient, and between all nodes, also
   * between the nodes of the two faces.
   */
  void TestStiffnessOfPatchIsSymmetric()
  {
    constexpr int nElements = 3;
    constexpr int nSide     = nElements + 1;
    constexpr int nDof      = 2 * 2 * nSide;

    const auto layers = gradedCard( { { 0.4, { 4000., 0.3 } }, { 0.6, { 20000., 0.2 } } } );

    Eigen::MatrixXd K = Eigen::MatrixXd::Zero( nDof, nDof );
    for ( int e = 0; e < nElements; e++ ) {
      const double x0 = 0.5 * e, x1 = 0.5 * ( e + 1 );
      // nodes of the bottom side, then of the top side, which lies at the thickness above
      const std::array< double, 8 > coordinates = { x0, 0., x1, 0., x0, thickness, x1, thickness };

      std::vector< double > dU( 8 );
      for ( int i = 0; i < 8; i++ )
        dU[i] = 1e-3 * std::sin( 0.9 * i + 0.3 + e );

      const auto evaluation = evaluateIncrement< GradedYInterfaceFiniteElement< 2, 4 >, 8 >( "LINEARELASTIC",
                                                                                             layers,
                                                                                             coordinates.data(),
                                                                                             dU );

      const int global[4] = { e, e + 1, nSide + e, nSide + e + 1 };
      for ( int a = 0; a < 4; a++ )
        for ( int b = 0; b < 4; b++ )
          for ( int i = 0; i < 2; i++ )
            for ( int j = 0; j < 2; j++ )
              K( 2 * global[a] + i, 2 * global[b] + j ) += evaluation.Ke[( 2 * b + j ) * 8 + 2 * a + i];
    }

    throwExceptionOnFailure( ( K - K.transpose() ).norm() < 1e-10 * K.norm(),
                             "The stiffness of a patch of graded elements is not symmetric." );
  }

  /**
   * A strip of two elastic layers, whose faces are displaced such that the top face has u_z = z0 + g x and u_x = 0, and
   * the bottom face is fixed, has an exact solution of plane strain, which is polynomial in the layers. It was derived
   * with the computer algebra system sympy for the thickness 0.02, z0 = 3e-5, g = 2e-3, the lower layer E = 4000,
   * nu = 0.3 of 40 % of the thickness and the upper layer E = 20000, nu = 0.2. The tractions of the faces do not depend
   * on the element size, so that the nodal forces of an interior node of a patch, divided by the width of the
   * elements, are the exact tractions. The tangential gradient of the correction of the microproblem is needed for the
   * shear traction: without it the traction of the bottom face is off by 3 %, and the one of the top face by a factor.
   */
  void TestPatchReproducesExactPlaneStrainTractions()
  {
    constexpr int    nElements = 6;
    constexpr int    nSide     = nElements + 1;
    constexpr int    nDof      = 4 * nSide;
    constexpr double h         = 0.02;
    constexpr double z0        = 3e-5;
    constexpr double g         = 2e-3;

    std::vector< double > card = { h, 0, 3, 2, 2, 0.4, 4000., 0.3, 0.6, 20000., 0.2, 1.0 };

    for ( const double width : { 0.2, 0.05, 0.01 } ) {
      Eigen::VectorXd dUglobal = Eigen::VectorXd::Zero( nDof );
      for ( int i = 0; i < nSide; i++ )
        dUglobal( 2 * ( nSide + i ) + 1 ) = z0 + g * width * i;

      Eigen::VectorXd P = Eigen::VectorXd::Zero( nDof );
      for ( int e = 0; e < nElements; e++ ) {
        const double                  x0 = width * e, x1 = width * ( e + 1 );
        const std::array< double, 8 > coordinates = { x0, 0., x1, 0., x0, h, x1, h };

        const int             global[4] = { e, e + 1, nSide + e, nSide + e + 1 };
        std::vector< double > dU( 8 );
        for ( int a = 0; a < 4; a++ )
          for ( int i = 0; i < 2; i++ )
            dU[2 * a + i] = dUglobal( 2 * global[a] + i );

        std::vector< double > stateVars;
        auto                  element = makeElement< GradedYInterfaceFiniteElement< 2, 4 > >( "LINEARELASTIC",
                                                                             card,
                                                                             coordinates.data(),
                                                                             stateVars );
        std::vector< double > U( 8, 0. ), Pe( 8, 0. ), Ke( 64, 0. );
        element->computeKernels( U.data(), dU.data(), Pe.data(), Ke.data(), 0., 1. );

        for ( int a = 0; a < 4; a++ )
          for ( int i = 0; i < 2; i++ )
            P( 2 * global[a] + i ) += Pe[2 * a + i];
      }

      // interior node 3 of the bottom and of the top side
      const double bottomShear = P( 2 * 3 ) / width;
      const double topShear    = P( 2 * ( nSide + 3 ) ) / width;
      const double normal      = P( 2 * ( nSide + 3 ) + 1 ) / width;

      const double expectedNormal = 14.809590973201692 * ( z0 + g * 3 * width ) / z0;

      throwExceptionOnFailure( std::abs( bottomShear + 6.397097558074329 ) < 1e-5 * 6.4,
                               "The shear traction of the bottom face differs from the exact plane strain solution." );
      throwExceptionOnFailure( std::abs( topShear - 0.05012999813074583 ) < 1e-5 * 6.4,
                               "The shear traction of the top face differs from the exact plane strain solution." );
      throwExceptionOnFailure( std::abs( normal - expectedNormal ) < 1e-8 * expectedNormal,
                               "The normal traction differs from the exact plane strain solution." );
    }
  }

  void TestElementsAreRegistered()
  {
    for ( const std::string name : { "HYLINE2", "HYQUAD4", "HYQUAD8" } ) {
      std::unique_ptr< MarmotElement > element( MarmotLibrary::MarmotElementFactory::createElement( name, 1 ) );
      throwExceptionOnFailure( element != nullptr, name + " is not registered with the element factory." );
    }

    std::unique_ptr< MarmotElement > hyquad4( MarmotLibrary::MarmotElementFactory::createElement( "HYQUAD4", 1 ) );
    throwExceptionOnFailure( hyquad4->getElementShape() == "hexa8", "HYQUAD4 reports the wrong result-file shape." );
  }

} // namespace

int main()
{
  auto tests = std::vector< std::function< void() > >{ TestTwoFacesAreEquilibratedElement3D,
                                                       TestTwoFacesAreEquilibratedElement2D,
                                                       TestTangentOfStackIsConsistent,
                                                       TestStiffnessOfPatchIsSymmetric,
                                                       TestPatchReproducesExactPlaneStrainTractions,
                                                       TestElementsAreRegistered };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
