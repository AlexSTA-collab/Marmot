#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/WarpingStabPressureMiniInterfaceFiniteElement.h"

namespace Marmot::Elements::Registration {
  using namespace MarmotLibrary;
  using namespace Marmot::FiniteElement::Quadrature;

  namespace {
    MarmotElement* makeWarpingMini( int id )
    {
      return new WarpingStabPressureMiniInterfaceFiniteElement( id,
                                                                FullIntegration,
                                                                WarpingStabPressureMiniInterfaceFiniteElement::
                                                                  Interface );
    }
  } // namespace

  const static bool WIQUAD4_STABP_MINI_isRegistered = MarmotElementFactory::registerElement( "WIQUAD4_STABP_MINI",
                                                                                             makeWarpingMini );
} // namespace Marmot::Elements::Registration
