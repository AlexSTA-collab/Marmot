#include "Marmot/GaussLobattoBBarInterfaceFiniteElement.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotFiniteElement.h"

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

  const static bool GLILINE2_BBAR_isRegistered = MarmotElementFactory::
    registerElement( "GLILINE2_BBAR",
                     makeFactoryFunction< GaussLobattoBBarInterfaceFiniteElement< 2, 4 >,
                                          FullIntegration,
                                          GaussLobattoBBarInterfaceFiniteElement< 2, 4 >::SectionType::Interface >() );

  const static bool GLIQUAD4_BBAR_isRegistered = MarmotElementFactory::
    registerElement( "GLIQUAD4_BBAR",
                     makeFactoryFunction< GaussLobattoBBarInterfaceFiniteElement< 3, 8 >,
                                          FullIntegration,
                                          GaussLobattoBBarInterfaceFiniteElement< 3, 8 >::SectionType::Interface >() );

} // namespace Marmot::Elements::Registration
