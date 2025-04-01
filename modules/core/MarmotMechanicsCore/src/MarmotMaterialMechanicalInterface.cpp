#include "Marmot/MarmotMaterialMechanicalInterface.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotLowerDimensionalStress.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotVoigt.h"
#include <iostream>

using namespace Eigen;

void MarmotMaterialMechanicalInterface::computeStress( double*  force,
                                                        double*  surface_stress,
                                                        double* dStress_dStrain,
                                                        const double* dU,
                                                        const double* dSurface_strain,
                                                        const double* normal,
                                                        const double* timeOld,
                                                        const double  dT,
                                                        double&       pNewDT)
{};
// void MarmotMaterialHypoElasticInterface::computeStress( Tensor1D&  force,
//                        Tensor2D&  surface_stress,
//                        Fastor::Tensor<double, 21,21>& dStress_dStrain,
//                        const Fastor::Tensor<double, 6,1>& dU,
//                        const Fastor::Tensor<double, 18,1>& dSurface_strain,
//                        const Tensor1D& normal,
//                        const double* timeOld,
//                        const double  dT,
//                        double&       pNewDT  )
//{
// const Map<const Matrix3d> FNew( FNew_ );
// const Map<const Matrix3d> FOld( FOld_ );

////Marmot::Vector6d dEps = 1./2 * ( Marmot::ContinuumMechanics::VoigtNotation::StrainToVoigt( H + H.tranpose() ) );

////computeStress (stress_, dStressDDStrain_, dEps.data(), timeOld, dT, pNewDT);
//}
