#include "Marmot/Marmot.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>

// Use namespaces for brevity
using namespace Marmot::Testing;
using namespace Marmot::ContinuumMechanics::Elasticity::Isotropic;
using namespace Marmot::ContinuumMechanics::Elasticity::TransverseIsotropic;

// Function to create a MarmotMaterialHypoElastic object
// Inputs:
// - materialName: The name of the material (e.g., "LINEARELASTIC")
// - materialProperties: Array of material parameters (e.g., Young's modulus, Poisson's ratio)
// - nMaterialProperties: Number of parameters in the materialProperties array
std::unique_ptr< MarmotMaterialHypoElastic > createMarmotMaterialHypoElastic( const std::string& materialName,
                                                                              const double*      materialProperties,
                                                                              int                nMaterialProperties )
{
  // Element label (arbitrary value)
  const int elLabel = 1;

  // Create the material object using Marmot's factory method
  auto mat = std::unique_ptr< MarmotMaterialHypoElastic >( dynamic_cast< MarmotMaterialHypoElastic* >(
    MarmotLibrary::MarmotMaterialFactory::createMaterial( MarmotLibrary::MarmotMaterialFactory::getMaterialCodeFromName(
                                                            materialName ),
                                                          materialProperties,
                                                          nMaterialProperties,
                                                          elLabel ) ) );

  return mat; // Return the created material object
}

// Function to test the isotropic material response for a normal strain increment
void testMaterialResponse()
{
  // Define material parameters (Young's modulus and Poisson's ratio)
  const double materialProperties[2] = { 20000, 0.25 };
  const int    nMaterialProperties   = 2;

  // Create the material object
  auto mat = createMarmotMaterialHypoElastic( "LINEARELASTIC", materialProperties, nMaterialProperties );

  // Define initial stress state (set to zero) and strain increment
  Eigen::Matrix< double, 6, 1 > stress = Eigen::Matrix< double, 6, 1 >::Zero();
  Eigen::Matrix< double, 6, 1 > dStrain;

  // Apply a small strain increment
  dStrain << 1e-3, 0., 0., 0., 0., 0.;

  // Define a matrix to store the tangent stiffness (stress-strain relation)
  Eigen::Matrix< double, 6, 6 > dStressDDStrain;

  // Define time parameters for the material response calculation
  const double timeOld = 0.0; // Previous time step
  const double dT      = 1.0; // Time increment
  double       pNewDT;        // Placeholder for the new time increment

  // Compute the stress response of the material
  mat->computeStress( stress.data(), dStressDDStrain.data(), dStrain.data(), &timeOld, dT, pNewDT );

  // Define the expected stress values for the applied strain increment
  Eigen::Matrix< double, 6, 1 > stressTarget;
  stressTarget << 24., 8., 8., 0., 0., 0.;

  // Compare the computed stress to the expected stress and throw an exception if they differ
  throwExceptionOnFailure( checkIfEqual< double >( stress, stressTarget, 1e-10 ),
                           "Stress computation failed for normal strain in " + std::string( __PRETTY_FUNCTION__ ) );
}

// Function to test the isotropic material response for a shear strain increment
void testShearMaterialResponse()
{
  // Define material parameters (Young's modulus and Poisson's ratio)
  const double materialProperties[2] = { 20000, 0.25 };
  const int    nMaterialProperties   = 2;

  // Create the material object
  auto mat = createMarmotMaterialHypoElastic( "LINEARELASTIC", materialProperties, nMaterialProperties );

  // Define initial stress state (set to zero) and strain increment
  Eigen::Matrix< double, 6, 1 > stress = Eigen::Matrix< double, 6, 1 >::Zero();
  Eigen::Matrix< double, 6, 1 > dStrain;

  // Apply a small strain increment
  dStrain << 0., 0., 0., 1e-3, 1e-3, 1e-3;

  // Define a matrix to store the tangent stiffness (stress-strain relation)
  Eigen::Matrix< double, 6, 6 > dStressDDStrain;

  // Define time parameters for the material response calculation
  const double timeOld = 0.0; // Previous time step
  const double dT      = 1.0; // Time increment
  double       pNewDT;        // Placeholder for the new time increment

  // Compute the stress response of the material
  mat->computeStress( stress.data(), dStressDDStrain.data(), dStrain.data(), &timeOld, dT, pNewDT );

  // Define the expected stress values for the applied strain increment
  Eigen::Matrix< double, 6, 1 > stressTargetShear;
  stressTargetShear << 0., 0., 0., 8., 8., 8.;

  // Compare the computed stress to the expected stress and throw an exception if they differ
  throwExceptionOnFailure( checkIfEqual< double >( stress, stressTargetShear, 1e-10 ),
                           "Stress computation failed for shear strain in " + std::string( __PRETTY_FUNCTION__ ) );
}

// Function to test the transversely isotropic material response for a normal strain increment
void testTransverseIsotropicMaterialResponse()
{
  // Define material parameters for transverse isotropy:
  // - E1: Longitudinal Young's modulus
  // - E2: Transverse Young's modulus
  // - nu12: Poisson's ratio in the 1-2 plane
  // - nu23: Poisson's ratio in the 2-3 plane
  // - G12: Shear modulus in the 1-2 plane
  const double materialProperties[11] = { 20000, 10000, 0.25, 0.3, 4000, 1, 0, 0, 0, 0, 1 };
  const int    nMaterialProperties    = 11;

  // Create the material object
  auto mat = createMarmotMaterialHypoElastic( "LINEARELASTIC", materialProperties, nMaterialProperties );

  // Define initial stress state (set to zero) and strain increment
  Eigen::Matrix< double, 6, 1 > stress = Eigen::Matrix< double, 6, 1 >::Zero();
  Eigen::Matrix< double, 6, 1 > dStrain;

  // Apply a small strain increment
  dStrain << 1e-3, 0., 0., 0., 0., 0.;

  // Define a matrix to store the tangent stiffness (stress-strain relation)
  Eigen::Matrix< double, 6, 6 > dStressDDStrain;

  // Define time parameters for the material response calculation
  const double timeOld = 0.0; // Previous time step
  const double dT      = 1.0; // Time increment
  double       pNewDT;        // Placeholder for the new time increment

  // Compute the stress response of the material
  mat->computeStress( stress.data(), dStressDDStrain.data(), dStrain.data(), &timeOld, dT, pNewDT );

  // Define the expected stress values for the applied strain increment
  Eigen::Matrix< double, 6, 1 > stressTarget_transverseIsotropic;
  stressTarget_transverseIsotropic << 31.1111111111111111, 11.1111111111111111, 11.1111111111111111, 0., 0., 0.;

  // Compare the computed stress to the expected stress and throw an exception if they differ
  throwExceptionOnFailure( checkIfEqual< double >( stress, stressTarget_transverseIsotropic, 1e-10 ),
                           "Stress computation failed for transverse isotropic material in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// Function to test the transversely isotropic material response for a shear strain increment
void testTransverseIsotropicShearMaterialResponse()
{
  // Define material parameters for transverse isotropy:
  // - E1: Longitudinal Young's modulus
  // - E2: Transverse Young's modulus
  // - nu12: Poisson's ratio in the 1-2 plane
  // - nu23: Poisson's ratio in the 2-3 plane
  // - G12: Shear modulus in the 1-2 plane
  const double materialProperties[11] = { 20000, 10000, 0.25, 0.3, 4000, 1, 0, 0, 0, 1, 0 };
  const int    nMaterialProperties    = 11;

  // Create the material object
  auto mat = createMarmotMaterialHypoElastic( "LINEARELASTIC", materialProperties, nMaterialProperties );

  // Define initial stress state (set to zero) and strain increment
  Eigen::Matrix< double, 6, 1 > stress = Eigen::Matrix< double, 6, 1 >::Zero();
  Eigen::Matrix< double, 6, 1 > dStrain;

  // Apply small strain increment
  dStrain << 0., 0., 0., 1e-3, 1e-3, 1e-3;

  // Define a matrix to store the tangent stiffness (stress-strain relation)
  Eigen::Matrix< double, 6, 6 > dStressDDStrain;

  // Define time parameters for the material response calculation
  const double timeOld = 0.0; // Previous time step
  const double dT      = 1.0; // Time increment
  double       pNewDT;        // Placeholder for the new time increment

  // Compute the stress response of the material
  mat->computeStress( stress.data(), dStressDDStrain.data(), dStrain.data(), &timeOld, dT, pNewDT );

  // Define the expected stress values for the applied shear strain increment
  Eigen::Matrix< double, 6, 1 > stressTargetShear_transverseIsotropic;
  stressTargetShear_transverseIsotropic << 0., 0., 0., 4., 4., 3.846153846153847;

  // Compare the computed stress to the expected stress and throw an exception if they differ
  throwExceptionOnFailure( checkIfEqual< double >( stress, stressTargetShear_transverseIsotropic, 1e-10 ),
                           "Stress computation failed for transverse isotropic material in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// Function to test the orthotropic material response for a normal strain increment
void testOrthotropicMaterialResponse()
{
  // Define material parameters for orthotropy:
  // - E1: Longitudinal Young's modulus
  // - E2: Transverse Young's modulus
  // - E3: Young's modulus in the third direction
  // - nu12: Poisson's ratio in the 1-2 plane
  // - nu23: Poisson's ratio in the 2-3 plane
  // - nu13: Poisson's ratio in the 1-3 plane
  // - G12: Shear modulus in the 1-2 plane
  // - G23: Shear modulus in the 2-3 plane
  // - G13: Shear modulus in the 1-3 plane
  const double materialProperties[15] = { 20000, 10000, 15000, 0.25, 0.3, 0.35, 4000, 5000, 6000, 1, 0, 0, 0, 1, 0 };
  const int    nMaterialProperties    = 15;

  // Create the material object
  auto mat = createMarmotMaterialHypoElastic( "LINEARELASTIC", materialProperties, nMaterialProperties );

  // Define initial stress state (set to zero) and strain increment
  Eigen::Matrix< double, 6, 1 > stress = Eigen::Matrix< double, 6, 1 >::Zero();
  Eigen::Matrix< double, 6, 1 > dStrain;

  // Apply a small strain increment
  dStrain << 1e-3, 0., 0., 0., 0., 0.;

  // Define a matrix to store the tangent stiffness (stress-strain relation)
  Eigen::Matrix< double, 6, 6 > dStressDDStrain;

  // Define time parameters for the material response calculation
  const double timeOld = 0.0; // Previous time step
  const double dT      = 1.0; // Time increment
  double       pNewDT;        // Placeholder for the new time increment

  // Compute the stress response of the material
  mat->computeStress( stress.data(), dStressDDStrain.data(), dStrain.data(), &timeOld, dT, pNewDT );

  // Define the expected stress values for the applied strain increment
  Eigen::Matrix< double, 6, 1 > stressTarget_orthotropic;
  stressTarget_orthotropic << 32.32091690544413, 11.002865329512895, 14.613180515759312, 0., 0., 0.;

  // Compare the computed stress to the expected stress and throw an exception if they differ
  throwExceptionOnFailure( checkIfEqual< double >( stress, stressTarget_orthotropic, 1e-10 ),
                           "Stress computation failed for orthotropic material in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// Function to test the orthotropic material response for a shear strain increment
void testOrthotropicShearMaterialResponse()
{
  // Define material parameters for orthotropy:
  // - E1: Longitudinal Young's modulus
  // - E2: Transverse Young's modulus
  // - E3: Young's modulus in the third direction
  // - nu12: Poisson's ratio in the 1-2 plane
  // - nu23: Poisson's ratio in the 2-3 plane
  // - nu13: Poisson's ratio in the 1-3 plane
  // - G12: Shear modulus in the 1-2 plane
  // - G23: Shear modulus in the 2-3 plane
  // - G13: Shear modulus in the 1-3 plane
  const double materialProperties[15] = { 20000, 10000, 15000, 0.25, 0.3, 0.35, 4000, 5000, 6000, 1, 0, 0, 0, 1, 0 };
  const int    nMaterialProperties    = 15;

  // Create the material object
  auto mat = createMarmotMaterialHypoElastic( "LINEARELASTIC", materialProperties, nMaterialProperties );

  // Define initial stress state (set to zero) and strain increment
  Eigen::Matrix< double, 6, 1 > stress = Eigen::Matrix< double, 6, 1 >::Zero();
  Eigen::Matrix< double, 6, 1 > dStrain;

  // Apply a small strain increment
  dStrain << 0., 0., 0., 1e-3, 1e-3, 1e-3;

  // Define a matrix to store the tangent stiffness (stress-strain relation)
  Eigen::Matrix< double, 6, 6 > dStressDDStrain;

  // Define time parameters for the material response calculation
  const double timeOld = 0.0; // Previous time step
  const double dT      = 1.0; // Time increment
  double       pNewDT;        // Placeholder for the new time increment

  // Compute the stress response of the material
  mat->computeStress( stress.data(), dStressDDStrain.data(), dStrain.data(), &timeOld, dT, pNewDT );

  // Define the expected stress values for the applied strain increment
  Eigen::Matrix< double, 6, 1 > stressTargetShear_orthotropic;
  stressTargetShear_orthotropic << 0., 0., 0., 4., 6., 5.;

  // Compare the computed stress to the expected stress and throw an exception if they differ
  throwExceptionOnFailure( checkIfEqual< double >( stress, stressTargetShear_orthotropic, 1e-10 ),
                           "Stress computation failed for orthotropic material in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// Function to test the orthotropic material response for a normal strain increment
void testOrthotropicMaterialResponseRotation()
{
  // Define material parameters for orthotropy:
  // - E1: Longitudinal Young's modulus
  // - E2: Transverse Young's modulus
  // - E3: Young's modulus in the third direction
  // - nu12: Poisson's ratio in the 1-2 plane
  // - nu23: Poisson's ratio in the 2-3 plane
  // - nu13: Poisson's ratio in the 1-3 plane
  // - G12: Shear modulus in the 1-2 plane
  // - G23: Shear modulus in the 2-3 plane
  // - G13: Shear modulus in the 1-3 plane
  const double materialProperties[15] =
    { 1000, 30, 30, 0.009, 0, 0, 50, 50, 50, 0.906307787, 0.422618262, 0, -0.422618262, 0.906307787, 0 };
  const int nMaterialProperties = 15;

  // Create the material object
  auto mat = createMarmotMaterialHypoElastic( "LINEARELASTIC", materialProperties, nMaterialProperties );

  // Define initial stress state (set to zero) and strain increment
  Eigen::Matrix< double, 6, 1 > stress = Eigen::Matrix< double, 6, 1 >::Zero();
  Eigen::Matrix< double, 6, 1 > dStrain;

  // Apply a small strain increment
  dStrain << -2.4e-4, -3.9e-4, 0., 1.76e-3, 0., 0.;

  // Define a matrix to store the tangent stiffness (stress-strain relation)
  Eigen::Matrix< double, 6, 6 > dStressDDStrain;

  // Define time parameters for the material response calculation
  const double timeOld = 0.0; // Previous time step
  const double dT      = 1.0; // Time increment
  double       pNewDT;        // Placeholder for the new time increment

  // Compute the stress response of the material
  mat->computeStress( stress.data(), dStressDDStrain.data(), dStrain.data(), &timeOld, dT, pNewDT );

  // Define the expected st ress values for the applied strain increment
  Eigen::Matrix< double, 6, 1 > stressTarget_orthotropic;
  stressTarget_orthotropic << 0.28394632, 0.08759518, 0., 0.19606294, 0., 0.;

  // Compare the computed stress to the expected stress and throw an exception if they differ
  throwExceptionOnFailure( checkIfEqual< double >( stress, stressTarget_orthotropic, 1e-8 ),
                           "Stress computation failed for orthotropic material in " +
                             std::string( __PRETTY_FUNCTION__ ) );
}

// Function to test getDensity for isotropic material
void testGetDensityIsotropic()
{
  // Define material parameters (Young's modulus, Poisson's ratio and Density)
  const double materialProperties[3] = { 20000, 0.25, 10 };
  const int    nMaterialProperties   = 3;

  auto mat = createMarmotMaterialHypoElastic( "LINEARELASTIC", materialProperties, nMaterialProperties );

  // Get the density
  mat->getDensity();

  // Defined density at start (materialProperties[3])
  double expectedDensityI = 10;

  // Compare the retrieved density with the given number
  throwExceptionOnFailure( checkIfEqual( mat->getDensity(), expectedDensityI, 1e-10 ),
                           "Density retrieval failed for isotropic material in " + std::string( __PRETTY_FUNCTION__ ) );
}

// Function to test getDensity for transversely isotropic material
void testGetDensityTransverselyIsotropic()
{
  // Define material parameters for transverse isotropy:
  // - E1: Longitudinal Young's modulus
  // - E2: Transverse Young's modulus
  // - nu12: Poisson's ratio in the 1-2 plane
  // - nu23: Poisson's ratio in the 2-3 plane
  // - G12: Shear modulus in the 1-2 plane
  const double materialProperties[12] = { 20000, 10000, 0.25, 0.3, 4000, 1, 0, 0, 0, 0, 1, 10 };
  const int    nMaterialProperties    = 12;

  auto mat = createMarmotMaterialHypoElastic( "LINEARELASTIC", materialProperties, nMaterialProperties );

  // Get the density
  mat->getDensity();

  // Defined density at start (materialProperties[3])
  double expectedDensityTI = 10;

  // Compare the retrieved density with the given number
  throwExceptionOnFailure( checkIfEqual( mat->getDensity(), expectedDensityTI, 1e-10 ),
                           "Density retrieval failed for isotropic material in " + std::string( __PRETTY_FUNCTION__ ) );
}

// Function to test getDensity for transversely isotropic material
void testGetDensityOrthotropic()
{
  // Define material parameters for orthotropy:
  // - E1: Longitudinal Young's modulus
  // - E2: Transverse Young's modulus
  // - E3: Young's modulus in the third direction
  // - nu12: Poisson's ratio in the 1-2 plane
  // - nu23: Poisson's ratio in the 2-3 plane
  // - nu13: Poisson's ratio in the 1-3 plane
  // - G12: Shear modulus in the 1-2 plane
  // - G23: Shear modulus in the 2-3 plane
  // - G13: Shear modulus in the 1-3 plane
  const double materialProperties[16] =
    { 1000, 30, 30, 0.009, 0, 0, 50, 0, 0, 0.906307787, 0.422618262, 0, -0.422618262, 0.906307787, 0, 10 };
  const int nMaterialProperties = 16;

  auto mat = createMarmotMaterialHypoElastic( "LINEARELASTIC", materialProperties, nMaterialProperties );

  // Get the density
  mat->getDensity();

  // Defined density at start (materialProperties[3])
  double expectedDensityO = 10;

  // Compare the retrieved density with the given number
  throwExceptionOnFailure( checkIfEqual( mat->getDensity(), expectedDensityO, 1e-10 ),
                           "Density retrieval failed for isotropic material in " + std::string( __PRETTY_FUNCTION__ ) );
}

int main()
{
  // Run tests
  testMaterialResponse();                         // Test for normal strain
  testShearMaterialResponse();                    // Test for shear strain
  testTransverseIsotropicMaterialResponse();      // Test for transverse isotropic normal strain
  testTransverseIsotropicShearMaterialResponse(); // Test for transverse isotropic shear strain
  testOrthotropicMaterialResponse();              // Test for orthotropic normal strain
  testOrthotropicShearMaterialResponse();         // Test for orthotropic shear strain
  testOrthotropicMaterialResponseRotation();      // Test for orthotropic normal strain with rotation
  testGetDensityIsotropic();                      // Test for density retrieval for isotropic case
  testGetDensityTransverselyIsotropic();          // Test for density retrieval for transversely isotropic case
  testGetDensityOrthotropic();                    // Test for density retrieval for orthotropic case
  return 0;
}
