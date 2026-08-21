#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/ZStabPressureInterfaceFiniteElement.h"

namespace Marmot::Elements::Registration {

  using namespace MarmotLibrary;
  using namespace Marmot::FiniteElement::Quadrature;

  MarmotElement* makeZStabPressureInterfaceFiniteElement( int elementID )
  {
    return new ZStabPressureInterfaceFiniteElement( elementID,
                                                    FullIntegration,
                                                    ZStabPressureInterfaceFiniteElement::SectionType::Interface );
  }

  const static bool
    ZIQUAD4_STABP_MINI_isRegistered = MarmotElementFactory::registerElement( "ZIQUAD4_STABP_MINI",
                                                                             makeZStabPressureInterfaceFiniteElement );

} // namespace Marmot::Elements::Registration
