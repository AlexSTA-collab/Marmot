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
                                                       TestElementsAreRegistered };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
