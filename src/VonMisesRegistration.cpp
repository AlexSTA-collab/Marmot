#include "Marmot/MarmotMaterialRegistrationHelper.h"
#include "Marmot/VonMises.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool VonMisesIsRegistered = MarmotMaterialFactory::
      registerMaterial( MaterialCode::VonMises,
                        "VONMISES",
                        makeDefaultMarmotMaterialFactoryFunction< class VonMisesModel >() );

  } // namespace Registration

} // namespace Marmot::Materials
