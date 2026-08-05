#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/YStabPressureMiniInterfaceFiniteElement.h"

namespace Marmot::Elements::Registration {
  using namespace MarmotLibrary;
  using namespace Marmot::FiniteElement::Quadrature;
  namespace {
    MarmotElement* makeMini( int id )
    {
      return new YStabPressureMiniInterfaceFiniteElement( id,
                                                          FullIntegration,
                                                          YStabPressureMiniInterfaceFiniteElement::Interface );
    }
  } // namespace
  const static bool YIQUAD4_STABP_MINI_isRegistered = MarmotElementFactory::registerElement( "YIQUAD4_STABP_MINI",
                                                                                             makeMini );
} // namespace Marmot::Elements::Registration
