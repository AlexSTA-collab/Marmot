#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/ZInterfaceFiniteElement.h"

namespace Marmot::Elements::Registration {

  using namespace MarmotLibrary;
  using namespace Marmot::FiniteElement::Quadrature;

  MarmotElement* makeZInterfaceFiniteElement( int elementID )
  {
    return new ZInterfaceFiniteElement( elementID, FullIntegration, ZInterfaceFiniteElement::SectionType::Interface );
  }

  const static bool ZIQUAD4_isRegistered = MarmotElementFactory::registerElement( "ZIQUAD4",
                                                                                  makeZInterfaceFiniteElement );

} // namespace Marmot::Elements::Registration
