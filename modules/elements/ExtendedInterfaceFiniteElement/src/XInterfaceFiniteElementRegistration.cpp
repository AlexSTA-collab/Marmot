#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/XInterfaceFiniteElement.h"

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

  const static bool XILINE2_isRegistered = MarmotElementFactory::
    registerElement( "XILINE2",
                     makeFactoryFunction< XInterfaceFiniteElement< 2, 4 >,
                                          FullIntegration,
                                          XInterfaceFiniteElement< 2, 4 >::SectionType::Interface >() );

  const static bool XIQUAD4_isRegistered = MarmotElementFactory::
    registerElement( "XIQUAD4",
                     makeFactoryFunction< XInterfaceFiniteElement< 3, 8 >,
                                          FullIntegration,
                                          XInterfaceFiniteElement< 3, 8 >::SectionType::Interface >() );

} // namespace Marmot::Elements::Registration
