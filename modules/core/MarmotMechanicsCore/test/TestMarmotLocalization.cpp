#include "Marmot/MarmotLocalization.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"
#include <Marmot/Marmot.h>
#include <Marmot/MarmotMaterial.h>
#include <Marmot/VonMises.h>

using namespace Marmot::Testing;

void testMarmotLocalization()
{
  using namespace Marmot::ContinuumMechanics::LocalizationAnalysis;

  // material properties
  const int                  vonMisesCode = 2;
  Eigen::Vector< double, 6 > materialProperties;
  materialProperties << 210000., 0.3, 200., 0., 0., 0.;

  // instantiate material
  auto material = std::unique_ptr< MarmotMaterialHypoElastic >( dynamic_cast< MarmotMaterialHypoElastic* >(
    MarmotLibrary::MarmotMaterialFactory::createMaterial( vonMisesCode,
                                                          materialProperties.data(),
                                                          materialProperties.size(),
                                                          0 ) ) );

  // set state vars
  if ( material->getNumberOfRequiredStateVars() > 1 ) {
    throw std::runtime_error( "Number of required state vars for Mises model changed!" );
  }
  double kappa = 0;
  material->assignStateVars( &kappa, 1 );

  Marmot::Vector6d stress;
  Marmot::Matrix6d dStress_dStrain;
  Marmot::Vector6d dStrain;
  const double     timeOld = 0.;
  const double     dT      = 1;
  double           pNewDt;

  // initialize stress on shear meridian close to yield surface
  stress << 115.4, -115.4, 0., 0., 0., 0.;

  // set strain increment to reach shear meridian
  dStrain << 1e-6, -1e-6, 0., 0., 0., 0.;

  // evaluate material
  material->computeStress( stress.data(), dStress_dStrain.data(), dStrain.data(), &timeOld, dT, pNewDt );

  // normal vector
  Marmot::Vector3d n;
  n << 1, -1, 0;
  n /= n.norm();

  throwExceptionOnFailure( checkIfEqual( computeAcousticTensor( dStress_dStrain, n ).determinant(), 0. ),
                           "for n = [1/√2, 1/√2, 0], determinant of acoustic tensor should be zero on shear meridian" );

  // normal vector
  n << 1, 0, 0;
  n /= n.norm();

  throwExceptionOnFailure( checkIfEqual( minimumDeterminantAcousticTensor( dStress_dStrain ) /
                                           computeAcousticTensor( dStress_dStrain, n ).determinant(),
                                         0.0301537,
                                         1e-6 ),
                           "test minimumDeterminantAcousticTensor" );
}

int main()
{
  testMarmotLocalization();

  return 0;
}
