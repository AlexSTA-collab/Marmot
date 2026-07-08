#include "Marmot/InterfaceFiniteElementFiniteStrain.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotTesting.h"

#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::Elements;

namespace {

  template < typename DerivedA, typename DerivedB >
  void assertMatrixNear( const Eigen::MatrixBase< DerivedA >& actual,
                         const Eigen::MatrixBase< DerivedB >& expected,
                         double                               tol,
                         const std::string&                   message )
  {
    throwExceptionOnFailure( actual.rows() == expected.rows() && actual.cols() == expected.cols(),
                             message + ": matrix shape mismatch." );

    const double err = ( actual - expected ).template lpNorm< Eigen::Infinity >();
    throwExceptionOnFailure( err < tol, message + ": max error = " + std::to_string( err ) );
  }

  std::unique_ptr< InterfaceFiniteElementFiniteStrain< 2, 4 > > makeTwoDimensionalFiniteStrainInterfaceElement()
  {
    constexpr int nDim   = 2;
    constexpr int nNodes = 4;

    const int  elId    = 4;
    const auto intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
    const auto secType = InterfaceFiniteElementFiniteStrain< nDim, nNodes >::SectionType::Interface;

    auto element = std::make_unique< InterfaceFiniteElementFiniteStrain< nDim, nNodes > >( elId, intType, secType );

    static std::array< double, nDim* nNodes > coordinates = {
      0.0,
      0.0,
      1.0,
      0.0,
      0.0,
      0.0,
      1.0,
      0.0,
    };
    element->assignNodeCoordinates( coordinates.data() );

    static std::array< double, 1 > elPropsVec = { 1.0 };
    ElementProperties              elProps( elPropsVec.data(), static_cast< int >( elPropsVec.size() ) );
    element->assignProperty( elProps );

    static std::array< double, 3 > materialProperties = { 0.1, 3500.0, 1500.0 };
    element->assignMaterial( "COMPRESSIBLENEOHOOKE", materialProperties.data(), materialProperties.size() );

    return element;
  }

  template < int nDim, int nNodes >
  void initializeStateAndMaterial( InterfaceFiniteElementFiniteStrain< nDim, nNodes >& element,
                                   std::vector< double >&                              stateVars )
  {
    stateVars.assign( element.getNumberOfRequiredStateVars(), 0.0 );
    element.assignStateVars( stateVars.data(), static_cast< int >( stateVars.size() ) );
    element.initializeYourself();
    element.setInitialConditions( MarmotElement::MarmotMaterialInitialization, nullptr );
  }

} // namespace

void TestTwoDimensionalFiniteStrainInterfaceElementComputesConsistentTangent()
{
  std::cout << "\n--- TestTwoDimensionalFiniteStrainInterfaceElementComputesConsistentTangent ---\n";

  constexpr int nDim      = 2;
  constexpr int nNodes    = 4;
  constexpr int totalNDof = nDim * nNodes;

  auto element = makeTwoDimensionalFiniteStrainInterfaceElement();

  std::vector< double > stateVars;
  initializeStateAndMaterial( *element, stateVars );

  Eigen::Matrix< double, totalNDof, 1 > dUEigen;
  dUEigen << 0.000, 0.000, 0.006, 0.001, 0.001, 0.003, 0.009, 0.006;

  std::vector< double > U( dUEigen.data(), dUEigen.data() + dUEigen.size() );
  std::vector< double > dU( dUEigen.data(), dUEigen.data() + dUEigen.size() );
  std::vector< double > Pe( totalNDof, 0.0 );
  std::vector< double > Ke( totalNDof * totalNDof, 0.0 );

  double time   = 0.0;
  double dT     = 0.1;
  double pNewDT = 1.0;

  element->computeYourself( U.data(), dU.data(), Pe.data(), Ke.data(), &time, dT, pNewDT );

  Eigen::Map< Eigen::Matrix< double, totalNDof, 1 > >                          PeActual( Pe.data() );
  Eigen::Map< Eigen::Matrix< double, totalNDof, totalNDof, Eigen::RowMajor > > KeActual( Ke.data() );

  throwExceptionOnFailure( pNewDT == 1.0, "Finite-strain 2D interface element requested a smaller time step." );
  throwExceptionOnFailure( PeActual.allFinite(), "Finite-strain 2D interface residual contains nan or inf." );
  throwExceptionOnFailure( KeActual.allFinite(), "Finite-strain 2D interface tangent contains nan or inf." );
  throwExceptionOnFailure( PeActual.template lpNorm< Eigen::Infinity >() > 0.0,
                           "Finite-strain 2D interface residual should be nonzero." );
  throwExceptionOnFailure( KeActual.template lpNorm< Eigen::Infinity >() > 0.0,
                           "Finite-strain 2D interface tangent should be nonzero." );

  auto computePeForTotalDisplacement = [&]( const Eigen::Matrix< double, totalNDof, 1 >& UTotal ) {
    auto fdElement = makeTwoDimensionalFiniteStrainInterfaceElement();

    std::vector< double > fdStateVars;
    initializeStateAndMaterial( *fdElement, fdStateVars );

    std::vector< double > Ulocal( UTotal.data(), UTotal.data() + UTotal.size() );
    std::vector< double > dUlocal( UTotal.data(), UTotal.data() + UTotal.size() );
    std::vector< double > Pelocal( totalNDof, 0.0 );
    std::vector< double > Kelocal( totalNDof * totalNDof, 0.0 );

    double timeLocal   = 0.0;
    double dTLocal     = 0.1;
    double pNewDTLocal = 1.0;

    fdElement->computeYourself( Ulocal.data(),
                                dUlocal.data(),
                                Pelocal.data(),
                                Kelocal.data(),
                                &timeLocal,
                                dTLocal,
                                pNewDTLocal );

    return Eigen::Map< Eigen::Matrix< double, totalNDof, 1 > >( Pelocal.data() ).eval();
  };

  Eigen::Matrix< double, totalNDof, totalNDof > KeFiniteDifference;
  const double                                  eps = 1e-7;

  for ( int j = 0; j < totalNDof; ++j ) {
    Eigen::Matrix< double, totalNDof, 1 > dUPlus  = dUEigen;
    Eigen::Matrix< double, totalNDof, 1 > dUMinus = dUEigen;

    dUPlus( j ) += eps;
    dUMinus( j ) -= eps;

    KeFiniteDifference.col(
      j ) = ( computePeForTotalDisplacement( dUPlus ) - computePeForTotalDisplacement( dUMinus ) ) / ( 2.0 * eps );
  }

  const double fdNorm = KeFiniteDifference.template lpNorm< Eigen::Infinity >();
  const double relErr = ( KeActual + KeFiniteDifference ).template lpNorm< Eigen::Infinity >() /
                        std::max( 1.0, fdNorm );

  throwExceptionOnFailure( relErr < 1e-5, "Finite-strain 2D interface tangent is not consistent with Ke = -dPe/ddU." );
}

int main()
{
  auto tests = std::vector< std::function< void() > >{
    TestTwoDimensionalFiniteStrainInterfaceElementComputesConsistentTangent,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
