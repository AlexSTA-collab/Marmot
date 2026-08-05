#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/YStabPressureInterfaceFiniteElement.h"

namespace Marmot::Elements::Registration {
  using namespace MarmotLibrary;
  using namespace Marmot::FiniteElement::Quadrature;
  namespace {
    MarmotElement* make( int id )
    {
      return new YStabPressureInterfaceFiniteElement( id,
                                                      FullIntegration,
                                                      YStabPressureInterfaceFiniteElement::Interface );
    }
  } // namespace
  const static bool YIQUAD4_STABP_isRegistered = MarmotElementFactory::registerElement( "YIQUAD4_STABP", make );
} // namespace Marmot::Elements::Registration
