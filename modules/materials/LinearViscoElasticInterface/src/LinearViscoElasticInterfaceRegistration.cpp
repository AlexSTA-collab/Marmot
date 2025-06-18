#include "../include/Marmot/LinearViscoElasticInterface.h"
#include "Marmot/MarmotMaterialRegistrationHelper.h"

namespace Marmot::Materials {

  namespace Registration {

    constexpr int LinearViscoElasticInterfaceCode = 1;

    using namespace MarmotLibrary;

    const static bool LinearViscoElasticIsRegistered = MarmotMaterialFactory::
      registerMaterial( LinearViscoElasticInterfaceCode,
                        "LINEARVISCOELASTICINTERFACE",
                        makeDefaultMarmotMaterialFactoryFunction< class LinearViscoElasticInterface >() );

  } // namespace Registration
} // namespace Marmot::Materials
