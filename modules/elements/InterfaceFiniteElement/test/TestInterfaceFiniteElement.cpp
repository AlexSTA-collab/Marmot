#include "Marmot/InterfaceFiniteElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotTesting.h"

#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::Elements;

void TestAngledInterfaceKinematics()
{
  std::cout << "\n--- TestAngledInterfaceKinematics ---\n";

  constexpr int nDim   = 2;
  constexpr int nNodes = 4;

  const int  elId    = 1;
  const auto intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
  const auto secType = InterfaceFiniteElement< nDim, nNodes >::SectionType::Interface;

  auto element = std::make_unique< InterfaceFiniteElement< nDim, nNodes > >( elId, intType, secType );

  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0 };

  element->assignNodeCoordinates( nodeCoordsVec.data() );

  const std::vector< double > elPropsVec = { 1.0 };
  ElementProperties           elProps( elPropsVec.data(), elPropsVec.size() );
  element->assignProperty( elProps );

  element->initializeYourself();

  for ( int q = 0; q < element->getNumberOfQuadraturePoints(); ++q ) {
    const auto& qp = element->qps[q];

    std::cout << "qp " << q << " xi = " << qp.xi.transpose() << " normal = " << qp.normal.transpose()
              << " norm = " << qp.normal.norm() << " sqrtDetG = " << qp.sqrtDetG << "\n";

    throwExceptionOnFailure( std::abs( qp.normal.norm() - 1.0 ) < 1e-12, "Interface normal is not normalized." );
  }
}

void TestLinearElasticInterfaceAssignment()
{
  std::cout << "\n--- TestLinearElasticInterfaceAssignment ---\n";

  constexpr int nDim   = 3;
  constexpr int nNodes = 8;

  const int  elId    = 1;
  const auto intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
  const auto secType = InterfaceFiniteElement< nDim, nNodes >::SectionType::Interface;

  auto element = std::make_unique< InterfaceFiniteElement< nDim, nNodes > >( elId, intType, secType );

  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0,

                                                0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0 };

  element->assignNodeCoordinates( nodeCoordsVec.data() );

  const std::vector< double > elPropsVec = { 1.0 };
  ElementProperties           elProps( elPropsVec.data(), elPropsVec.size() );
  element->assignProperty( elProps );

  std::vector< double > matProps = { 1000.0, 0.3, 10.0, 10.0 };

  std::cout << "assignMaterial...\n";

  element->assignMaterial( "LINEARELASTICINTERFACE", matProps.data(), static_cast< int >( matProps.size() ) );

  for ( size_t q = 0; q < element->qps.size(); ++q ) {
    throwExceptionOnFailure( static_cast< bool >( element->qps[q].material ), "Material assignment returned nullptr." );
  }

  std::cout << "getNumberOfRequiredStateVars...\n";

  const int nStateVarsTotal = element->getNumberOfRequiredStateVars();

  std::cout << "nStateVarsTotal = " << nStateVarsTotal << "\n";

  std::vector< double > stateVars( nStateVarsTotal, 0.0 );

  std::cout << "assignStateVars...\n";

  element->assignStateVars( stateVars.data(), nStateVarsTotal );

  std::cout << "initializeYourself...\n";

  element->initializeYourself();

  std::cout << "initialize material qps...\n";

  for ( auto& qp : element->qps ) {
    throwExceptionOnFailure( static_cast< bool >( qp.material ),
                             "Material pointer is null before initializeYourself." );

    qp.material->initializeYourself();
  }

  std::vector< double > Q( nNodes * nDim, 0.0 );
  std::vector< double > dQ( nNodes * nDim, 0.01 );
  std::vector< double > Pe( nNodes * nDim, 0.0 );
  std::vector< double > Ke( nNodes * nDim * nNodes * nDim, 0.0 );

  double time   = 0.0;
  double dT     = 0.1;
  double pNewDT = 1.0;

  std::cout << "computeYourself...\n";

  element->computeYourself( Q.data(), dQ.data(), Pe.data(), Ke.data(), &time, dT, pNewDT );

  Eigen::Map< Eigen::Matrix< double, nNodes * nDim, nNodes * nDim > > K_matrix( Ke.data() );

  throwExceptionOnFailure( K_matrix.allFinite(), "Ke contains nan or inf." );

  throwExceptionOnFailure( element->getNumberOfQuadraturePoints() == 4, "Incorrect number of quadrature points." );

  throwExceptionOnFailure( element->getNDofPerElement() == nNodes * nDim, "Incorrect number of DOFs." );
}

int main()
{
  auto tests = std::vector< std::function< void() > >{ TestAngledInterfaceKinematics,
                                                       TestLinearElasticInterfaceAssignment };

  executeTestsAndCollectExceptions( tests );

  return 0;
}