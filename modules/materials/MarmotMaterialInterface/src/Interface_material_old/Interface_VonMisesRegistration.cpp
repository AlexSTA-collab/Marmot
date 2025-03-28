#include "Marmot/MarmotMaterialRegistrationHelper.h"
#include "Marmot/InterfaceVonMises.h"

namespace Marmot::Materials {

  namespace Registration {

    constexpr int vonMisesCode = 2;

    using namespace MarmotLibrary;

    const static bool InterfaceVonMisesIsRegistered = MarmotMaterialFactory::
      registerMaterial( InterfacevonMisesCode, "INTVONMISES", makeDefaultMarmotMaterialFactoryFunction< class InterfaceVonMisesModel >() );

  } // namespace Registration

} // namespace Marmot::Materials
