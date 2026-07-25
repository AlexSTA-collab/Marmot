#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/YInterfaceFiniteElement.h"

namespace Marmot::Elements::Registration {

  template < class T,
             Marmot::FiniteElement::Quadrature::IntegrationTypes integrationType,
             typename T::SectionType                             sectionType >
  MarmotLibrary::MarmotElementFactory::elementFactoryFunction makeFactoryFunction()
  {
    return []( int elementID ) -> MarmotElement* { return new T( elementID, integrationType, sectionType ); };
  }

  using namespace MarmotLibrary;
  using namespace Marmot::FiniteElement::Quadrature;

  const static bool YILINE2_isRegistered = MarmotElementFactory::
    registerElement( "YILINE2",
                     makeFactoryFunction< YInterfaceFiniteElement< 2, 4 >,
                                          FullIntegration,
                                          YInterfaceFiniteElement< 2, 4 >::SectionType::Interface >() );

  const static bool YIQUAD4_isRegistered = MarmotElementFactory::
    registerElement( "YIQUAD4",
                     makeFactoryFunction< YInterfaceFiniteElement< 3, 8 >,
                                          FullIntegration,
                                          YInterfaceFiniteElement< 3, 8 >::SectionType::Interface >() );

} // namespace Marmot::Elements::Registration
