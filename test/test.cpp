#include "Marmot/CompressibleNeoHooke.h"
#include "Marmot/Marmot.h"
#include "Marmot/MarmotMaterialFiniteStrain.h"
#include "Marmot/MarmotMicromorphicTensorBasics.h"

int main()
{

  using namespace Marmot::Materials;
  std::array< double, 2 > materialProperties_ = { 3500, 1500 };
  const double                  nMaterialProperties = 2;
  const int                     elLabel             = 1;

  CompressibleNeoHooke mat = CompressibleNeoHooke( &materialProperties_[0], nMaterialProperties, elLabel );

  CompressibleNeoHooke::DeformationIncrement< 3 > def;
  CompressibleNeoHooke::TimeIncrement             timeInc = { 0, 0.1 };
  double                                          pNewDT;

  CompressibleNeoHooke::ConstitutiveResponse< 3 > response;
  CompressibleNeoHooke::AlgorithmicModuli< 3 >    tangent;

  def.F_np = Marmot::FastorStandardTensors::Spatial3D::I;
  def.F_np( 0, 0 ) += 1e-4;
  def.F_np( 1, 1 ) += 1e-4;
  def.F_np( 2, 2 ) += 1e-4;


  mat.computeStress( response, tangent, def, timeInc, pNewDT );

  std::cout << "deformation=\n" << def.F_np << std::endl;
  std::cout << "stress = \n" << response.S << std::endl;
  std::cout << "tangent = \n" << tangent.dS_dF << std::endl;

  return 0;
}
