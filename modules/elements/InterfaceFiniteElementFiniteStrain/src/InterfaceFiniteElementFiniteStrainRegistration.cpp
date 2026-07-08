#include "Marmot/InterfaceFiniteElementFiniteStrain.h"
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

  const static bool ILINE2UL_isRegistered = MarmotElementFactory::
    registerElement( "ILINE2UL",
                     makeFactoryFunction< InterfaceFiniteElementFiniteStrain< 2, 4 >,
                                          FullIntegration,
                                          InterfaceFiniteElementFiniteStrain< 2, 4 >::SectionType::Interface >() );

  const static bool IQUAD4UL_isRegistered = MarmotElementFactory::
    registerElement( "IQUAD4UL",
                     makeFactoryFunction< InterfaceFiniteElementFiniteStrain< 3, 8 >,
                                          FullIntegration,
                                          InterfaceFiniteElementFiniteStrain< 3, 8 >::SectionType::Interface >() );

} // namespace Marmot::Elements::Registration
