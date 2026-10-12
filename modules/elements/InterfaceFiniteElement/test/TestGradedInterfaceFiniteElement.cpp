#include "Marmot/InterfaceFiniteElement.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotTesting.h"

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

  constexpr double thickness = 0.18;

  /// Material properties of the graded layer: [h, mode, nGauss, K, m, (position, bulk properties) x K].
  std::vector< double > gradedCard( const std::vector< std::pair< double, std::vector< double > > >& layers,
                                    int                                                              mode   = 0,
                                    int                                                              nGauss = 2 )
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
    return properties;
  }

  /// Coordinates of an inclined, non-planar-looking interface: first the nodes of one side, then the other side.
  const std::array< double, 24 > coordinates3d = { -0.5, -0.5, 0.088163, 0.5,  -0.5, -0.088163, 0.5, 0.5,  -0.088163,
                                                   -0.5, 0.5,  0.088163, -0.5, -0.5, 0.188163,  0.5, -0.5, 0.011837,
                                                   0.5,  0.5,  0.011837, -0.5, 0.5,  0.188163 };

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

  /// The element with a homogeneous graded layer is the plain element, in 2D and 3D.
  template < int nDim, int nNodes >
  void checkHomogeneousLayerIsPlainElement( const double* coordinates )
  {
    constexpr int nDof = nDim * nNodes;

    const std::vector< double > bulk   = { 4000., 0.3 };
    const auto                  graded = gradedCard( { { 0.4, bulk }, { 0.6, bulk } } );
    const std::vector< double > plain  = { bulk[0], bulk[1], thickness };

    for ( const double scale : { 1e-4, 1e-3 } ) {
      const auto dU = arbitraryIncrement( nDof, scale );

      const auto expected = evaluateIncrement< InterfaceFiniteElement< nDim, nNodes >, nDof >( "LINEARELASTIC",
                                                                                               plain,
                                                                                               coordinates,
                                                                                               dU );
      const auto actual   = evaluateIncrement< GradedInterfaceFiniteElement< nDim, nNodes >, nDof >( "LINEARELASTIC",
                                                                                                   graded,
                                                                                                   coordinates,
                                                                                                   dU );

      throwExceptionOnFailure( checkIfEqual< double >( Eigen::Map< const Eigen::VectorXd >( actual.Pe.data(), nDof ),
                                                       Eigen::Map< const Eigen::VectorXd >( expected.Pe.data(), nDof ),
                                                       1e-9 * scale ),
                               "The residual of a homogeneous graded layer is not the one of the plain element." );
      throwExceptionOnFailure( checkIfEqual< double >( Eigen::Map< const Eigen::VectorXd >( actual.Ke.data(),
                                                                                            nDof * nDof ),
                                                       Eigen::Map< const Eigen::VectorXd >( expected.Ke.data(),
                                                                                            nDof * nDof ),
                                                       1e-9 ),
                               "The tangent of a homogeneous graded layer is not the one of the plain element." );
    }
  }

  void TestHomogeneousLayerIsPlainElement3D()
  {
    checkHomogeneousLayerIsPlainElement< 3, 8 >( coordinates3d.data() );
  }

  void TestHomogeneousLayerIsPlainElement2D()
  {
    const std::array< double, 8 > coordinates2d = { 0., 0., 1., 0., 0., 0.1, 1., 0.1 };
    checkHomogeneousLayerIsPlainElement< 2, 4 >( coordinates2d.data() );
  }

  /// The tangent of the element is the derivative of its residual, also for plastic layers with different history.
  template < class Element, int nDof >
  void checkTangentIsDerivativeOfResidual( const std::string&           materialName,
                                           const std::vector< double >& materialProperties,
                                           const double*                coordinates,
                                           double                       scale )
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
        throwExceptionOnFailure( std::abs( numerical - base.Ke[i + j * nDof] ) < 1e-4 * scaleOfKe,
                                 "The tangent of the graded element is not the derivative of its residual at (" +
                                   std::to_string( i ) + "," + std::to_string( j ) + ")." );
      }
    }
  }

  void TestTangentOfElasticStackIsConsistent()
  {
    const auto layers = gradedCard( { { 0.3, { 4000., 0.3 } }, { 0.7, { 20000., 0.2 } } } );
    checkTangentIsDerivativeOfResidual< GradedInterfaceFiniteElement< 3, 8 >, 24 >( "LINEARELASTIC",
                                                                                    layers,
                                                                                    coordinates3d.data(),
                                                                                    1e-3 );
  }

  void TestTangentOfPlasticStackIsConsistent()
  {
    // the soft layer yields in the increment, the hard one does not
    const auto layers = gradedCard(
      { { 0.4, { 1e4, 0.3, 40., 10., 0., 1., 2400. } }, { 0.6, { 2e4, 0.3, 4000., 10., 0., 1., 2400. } } } );
    checkTangentIsDerivativeOfResidual< GradedInterfaceFiniteElement< 3, 8 >, 24 >( "VONMISES",
                                                                                    layers,
                                                                                    coordinates3d.data(),
                                                                                    2e-3 );

    const std::array< double, 8 > coordinates2d = { 0., 0., 1., 0., 0., 0.1, 1., 0.1 };
    checkTangentIsDerivativeOfResidual< GradedInterfaceFiniteElement< 2, 4 >, 8 >( "VONMISES",
                                                                                   layers,
                                                                                   coordinates2d.data(),
                                                                                   2e-3 );
  }

  void TestTangentOfGradedLayerIsConsistent()
  {
    const auto layers = gradedCard( { { 0.0, { 4000., 0.3 } }, { 1.0, { 16000., 0.3 } } }, 1, 4 );
    checkTangentIsDerivativeOfResidual< GradedInterfaceFiniteElement< 3, 8 >, 24 >( "LINEARELASTIC",
                                                                                    layers,
                                                                                    coordinates3d.data(),
                                                                                    1e-3 );
  }

  void TestElementsAreRegistered()
  {
    for ( const std::string name : { "HILINE2", "HIQUAD4", "HIQUAD8" } ) {
      std::unique_ptr< MarmotElement > element( MarmotLibrary::MarmotElementFactory::createElement( name, 1 ) );
      throwExceptionOnFailure( element != nullptr, name + " is not registered with the element factory." );
    }

    std::unique_ptr< MarmotElement > hiquad4( MarmotLibrary::MarmotElementFactory::createElement( "HIQUAD4", 1 ) );
    throwExceptionOnFailure( hiquad4->getElementShape() == "hexa8", "HIQUAD4 reports the wrong result-file shape." );
  }

} // namespace

int main()
{
  auto tests = std::vector< std::function< void() > >{ TestHomogeneousLayerIsPlainElement3D,
                                                       TestHomogeneousLayerIsPlainElement2D,
                                                       TestTangentOfElasticStackIsConsistent,
                                                       TestTangentOfPlasticStackIsConsistent,
                                                       TestTangentOfGradedLayerIsConsistent,
                                                       TestElementsAreRegistered };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
