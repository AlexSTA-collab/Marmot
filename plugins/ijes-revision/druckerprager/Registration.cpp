#include "Marmot/DruckerPrager.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

namespace Marmot::Materials::Registration {
  const static bool DruckerPragerIsRegistered = MarmotLibrary::MarmotMaterialHypoElasticFactory::registerMaterial<
    DruckerPragerModel >( "DRUCKERPRAGER" );
}
