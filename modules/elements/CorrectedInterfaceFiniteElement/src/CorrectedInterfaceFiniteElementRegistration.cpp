#include "Marmot/CorrectedInterfaceFiniteElement.h"
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

  const static bool EILINE2_isRegistered = MarmotElementFactory::
    registerElement( "EILINE2",
                     makeFactoryFunction< CorrectedInterfaceFiniteElement< 2, 4 >,
                                          FullIntegration,
                                          CorrectedInterfaceFiniteElement< 2, 4 >::SectionType::Interface >() );

  const static bool EIQUAD4_isRegistered = MarmotElementFactory::
    registerElement( "EIQUAD4",
                     makeFactoryFunction< CorrectedInterfaceFiniteElement< 3, 8 >,
                                          FullIntegration,
                                          CorrectedInterfaceFiniteElement< 3, 8 >::SectionType::Interface >() );

} // namespace Marmot::Elements::Registration
