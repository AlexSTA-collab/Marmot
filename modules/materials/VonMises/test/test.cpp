#include "Marmot/MarmotMaterialPointSolverHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/VonMises.h"

using namespace Marmot::Testing;
using namespace Marmot::Solvers;

void testVonMisesIncrementalPotentialIsConsistent()
{
  const double                              materialProperties[6] = { 210000., 0.3, 200., 2100., 20., 20. };
  Marmot::Materials::VonMisesModel          material( materialProperties, 6, 1 );
  const MarmotMaterialHypoElastic::timeInfo timeInfo{ 1., 1. };

  struct Evaluation {
    Marmot::Vector6d stress;
    Marmot::Matrix6d tangent;
    double           potential;
  };

  const auto evaluate = [&]( const Marmot::Vector6d& strainIncrement ) {
    Eigen::VectorXd stateVars( material.getNumberOfRequiredStateVars() );
    material.initializeYourself( stateVars.data(), stateVars.size() );
    Evaluation                         result{ Marmot::Vector6d::Zero(), Marmot::Matrix6d::Zero(), 0.0 };
    MarmotMaterialHypoElastic::state3D state{ result.stress, 0.0, 0.0, stateVars.data() };
    const bool                         exact = material.computeStressAndIncrementalPotential( state,
                                                                      result.tangent,
                                                                      strainIncrement,
                                                                      timeInfo,
                                                                      result.potential );
    throwExceptionOnFailure( exact, "Von Mises did not provide its exact incremental potential." );
    result.stress = state.stress;
    return result;
  };

  Marmot::Vector6d strainIncrement;
  strainIncrement << 0.006, -0.0015, -0.002, 0.001, -0.0007, 0.0009;
  const Evaluation base = evaluate( strainIncrement );

  Marmot::Vector6d finiteDifferenceGradient = Marmot::Vector6d::Zero();
  Marmot::Matrix6d finiteDifferenceHessian  = Marmot::Matrix6d::Zero();
  constexpr double perturbation             = 1e-7;
  for ( int column = 0; column < 6; ++column ) {
    Marmot::Vector6d plus  = strainIncrement;
    Marmot::Vector6d minus = strainIncrement;
    plus[column] += perturbation;
    minus[column] -= perturbation;
    const Evaluation plusEvaluation  = evaluate( plus );
    const Evaluation minusEvaluation = evaluate( minus );
    finiteDifferenceGradient[column] = ( plusEvaluation.potential - minusEvaluation.potential ) /
                                       ( 2.0 * perturbation );
    finiteDifferenceHessian.col( column ) = ( plusEvaluation.stress - minusEvaluation.stress ) / ( 2.0 * perturbation );
  }

  const double gradientError = ( finiteDifferenceGradient - base.stress ).norm() / std::max( 1.0, base.stress.norm() );
  const double hessianError  = ( finiteDifferenceHessian - base.tangent ).norm() / std::max( 1.0, base.tangent.norm() );
  throwExceptionOnFailure( gradientError < 2e-7, "Von Mises incremental-potential gradient does not match stress." );
  throwExceptionOnFailure( hessianError < 2e-6,
                           "Von Mises incremental-potential Hessian does not match the algorithmic tangent." );

  const double     shearModulus             = materialProperties[0] / ( 2.0 * ( 1.0 + materialProperties[1] ) );
  const double     yieldShearStrain         = materialProperties[2] / ( std::sqrt( 3.0 ) * shearModulus );
  constexpr double transitionPerturbation   = 1e-9;
  Marmot::Vector6d belowYield               = Marmot::Vector6d::Zero();
  Marmot::Vector6d aboveYield               = Marmot::Vector6d::Zero();
  belowYield[3]                             = yieldShearStrain - transitionPerturbation;
  aboveYield[3]                             = yieldShearStrain + transitionPerturbation;
  const Evaluation belowYieldEvaluation     = evaluate( belowYield );
  const Evaluation aboveYieldEvaluation     = evaluate( aboveYield );
  const double     transitionPotentialSlope = ( aboveYieldEvaluation.potential - belowYieldEvaluation.potential ) /
                                          ( 2.0 * transitionPerturbation );
  const double yieldShearStress = materialProperties[2] / std::sqrt( 3.0 );
  throwExceptionOnFailure( std::abs( transitionPotentialSlope - yieldShearStress ) < 2e-5,
                           "Von Mises incremental potential is not continuously differentiable at first yield." );
}

void testVonMisesCoordinateInvariance()
{
  // material properties
  std::vector< double > materialProperties = { 210000., 0.3, 200., 2100., 20., 20 };
  auto                  solveropts         = MarmotMaterialPointSolverHypoElastic::SolverOptions();
  std::string           matName            = "VONMISES";

  // create material point solver instance
  auto solver = MarmotMaterialPointSolverHypoElastic( matName,
                                                      &materialProperties[0],
                                                      materialProperties.size(),
                                                      solveropts );

  // define a step
  MarmotMaterialPointSolverHypoElastic::Step step;

  // define step parameters
  step.isStrainComponentControlled = { true, true, true, true, true, true };
  step.isStressComponentControlled = { false, false, false, false, false, false };
  step.stressIncrementTarget       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  step.strainIncrementTarget       = { 0.01, 0., 0., 0., 0.03, 0.0 };
  step.dTStart                     = 1.0;

  // add step to solver
  solver.addStep( step );
  // check coordinate invariance with Turbokreisel
  throwExceptionOnFailure( spinTurbokreisel( solver, 1e-10, 1e-8 ), "Turbokreisel failed!" );
}

void testVonMises()
{
  // material properties
  std::vector< double > materialProperties = { 210000., 0.3, 200., 2100., 20., 20 };
  auto                  solveropts         = MarmotMaterialPointSolverHypoElastic::SolverOptions();
  std::string           matName            = "VONMISES";

  // create material point solver instance
  auto solver = MarmotMaterialPointSolverHypoElastic( matName,
                                                      &materialProperties[0],
                                                      materialProperties.size(),
                                                      solveropts );

  // define a step
  MarmotMaterialPointSolverHypoElastic::Step step;

  // define step parameters
  step.isStrainComponentControlled = { true, true, true, true, true, true };
  step.isStressComponentControlled = { false, false, false, false, false, false };
  step.stressIncrementTarget       = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  step.strainIncrementTarget       = { 0.00839244, 0.00089344, -0.00703916, 0.00013635, 0.00160548, 0.00572825 };
  step.dTStart                     = 1.0;

  // add step to solver
  solver.addStep( step );

  solver.solve();

  // read history
  auto history = solver.getHistory();
  // reference solution
  const Marmot::Vector6d stressTarget =
    { 511.262747695, 395.408929527, 272.856322779, 1.0532516407, 12.4017194288, 44.2485420671 };

  // compare solutions
  throwExceptionOnFailure( checkIfEqual< double >( history.back().stress, stressTarget, 1e-9 ),
                           "comparison with reference solution failed" );
}

int main()
{
  std::vector< std::function< void( void ) > > tests = {
    testVonMisesIncrementalPotentialIsConsistent,
    testVonMises,
    testVonMisesCoordinateInvariance,
  };

  executeTestsAndCollectExceptions( tests );
  return 0;
}
