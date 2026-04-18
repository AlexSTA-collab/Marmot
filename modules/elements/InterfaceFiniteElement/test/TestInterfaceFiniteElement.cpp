#include "Marmot/InterfaceFiniteElement.h"
#include "Marmot/LinearElasticInterface.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotTesting.h"
#include <cmath>
#include <iostream>

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::Elements;

void TestAngledInterfaceKinematics()
{
  constexpr int nDim    = 2;
  constexpr int nNodes  = 4;
  const int     elId    = 1;
  const auto    intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
  const auto    secType = InterfaceFiniteElement< nDim, nNodes >::SectionType::Solid;

  auto element = std::make_unique< InterfaceFiniteElement< nDim, nNodes > >( elId, intType, secType );

  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0 };

  element->assignNodeCoordinates( nodeCoordsVec.data() );

  const std::vector< double > elPropsVec = { 1.0 };
  ElementProperties           elProps( elPropsVec.data(), elPropsVec.size() );
  element->assignProperty( elProps );

  element->initializeYourself();

  double invSqrt2 = 1.0 / std::sqrt( 2.0 );

  for ( int q = 0; q < element->getNumberOfQuadraturePoints(); ++q ) {
    auto normal = element->qps[q].normal;
    throwExceptionOnFailure( checkIfEqual( normal( 0 ), -invSqrt2 ), "Incorrect normal(0)" );
    throwExceptionOnFailure( checkIfEqual( normal( 1 ), invSqrt2 ), "Incorrect normal(1)" );
  }
}

void TestLinearElasticInterfaceAssignment()
{
  constexpr int nDim    = 3;
  constexpr int nNodes  = 8;
  const int     elId    = 1;
  const auto    intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
  const auto    secType = InterfaceFiniteElement< nDim, nNodes >::SectionType::Solid;

  auto element = std::make_unique< InterfaceFiniteElement< nDim, nNodes > >( elId, intType, secType );

  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0,
                                                0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0 };
  element->assignNodeCoordinates( nodeCoordsVec.data() );
  std::vector< double > matProps = { 1000.0, 0.3, 10.0, 10.0 };

  for ( auto& qp : element->qps ) {
    qp.material = std::unique_ptr< MarmotMaterialHypoElasticInterface >(
      dynamic_cast< MarmotMaterialHypoElasticInterface* >(
        MarmotLibrary::MarmotMaterialHypoElasticInterfaceFactory::createMaterial( "LINEARELASTICINTERFACE",
                                                                                  matProps.data(),
                                                                                  matProps.size(),
                                                                                  elId ) ) );
  }

  const int             nStateVarsTotal = element->getNumberOfRequiredStateVars();
  std::vector< double > stateVars( nStateVarsTotal, 0.0 );
  element->assignStateVars( stateVars.data(), nStateVarsTotal );

  element->initializeYourself();
  for ( auto& qp : element->qps ) {
    qp.material->setCharacteristicElementLength( std::sqrt( 4 * qp.detJ ) );
    qp.material->initializeYourself();
  }

  std::vector< double > Q( nNodes * nDim, 0.0 );
  std::vector< double > dQ( nNodes * nDim, 0.01 );
  std::vector< double > Pe( nNodes * nDim, 0.0 );
  std::vector< double > Ke( nNodes * nDim * nNodes * nDim, 0.0 );
  double                time   = 0.0;
  double                dT     = 0.1;
  double                pNewDT = 1.0;

  element->computeYourself( Q.data(), dQ.data(), Pe.data(), Ke.data(), &time, dT, pNewDT );

  Eigen::Map< Eigen::MatrixXd > K_matrix( Ke.data(), nNodes * nDim, nNodes * nDim );
  double                        toleranceSymmetry = 1e-10;
  bool                          isSymmetric       = K_matrix.isApprox( K_matrix.transpose(), toleranceSymmetry );
  throwExceptionOnFailure( isSymmetric, "Stiffness matrix Ke is not symmetric!" );

  int nQP_expected = ( nDim == 3 ) ? 4 : 2;
  throwExceptionOnFailure( element->getNumberOfQuadraturePoints() == nQP_expected,
                           "Incorrect number of quadrature points." );
  int nDof = element->getNDofPerElement();
  throwExceptionOnFailure( nDof == nNodes * nDim, "Incorrect number of DOFs." );
}

int main()
{
  auto tests = std::vector< std::function< void() > >{ TestAngledInterfaceKinematics,
                                                       TestLinearElasticInterfaceAssignment };
  executeTestsAndCollectExceptions( tests );
  return 0;
}
