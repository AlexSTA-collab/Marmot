#include "Marmot/GaussLobattoInterfaceFiniteElement.h"
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

  const static bool GLILINE2_isRegistered = MarmotElementFactory::
    registerElement( "GLILINE2",
                     makeFactoryFunction< GaussLobattoInterfaceFiniteElement< 2, 4 >,
                                          FullIntegration,
                                          GaussLobattoInterfaceFiniteElement< 2, 4 >::SectionType::Interface >() );

  const static bool GLIQUAD4_isRegistered = MarmotElementFactory::
    registerElement( "GLIQUAD4",
                     makeFactoryFunction< GaussLobattoInterfaceFiniteElement< 3, 8 >,
                                          FullIntegration,
                                          GaussLobattoInterfaceFiniteElement< 3, 8 >::SectionType::Interface >() );

  // Internal station-count variants for the Section 13 convergence study
  // ONLY -- not a user-facing production input. GLIQUAD4 (five stations,
  // above) remains the validated, registered production element.
  const static bool GLIQUAD4_N3_isRegistered = MarmotElementFactory::
    registerElement( "GLIQUAD4_N3",
                     makeFactoryFunction< GaussLobattoInterfaceFiniteElement< 3, 8, 3 >,
                                          FullIntegration,
                                          GaussLobattoInterfaceFiniteElement< 3, 8, 3 >::SectionType::Interface >() );

  const static bool GLIQUAD4_N4_isRegistered = MarmotElementFactory::
    registerElement( "GLIQUAD4_N4",
                     makeFactoryFunction< GaussLobattoInterfaceFiniteElement< 3, 8, 4 >,
                                          FullIntegration,
                                          GaussLobattoInterfaceFiniteElement< 3, 8, 4 >::SectionType::Interface >() );

  const static bool GLIQUAD4_N7_isRegistered = MarmotElementFactory::
    registerElement( "GLIQUAD4_N7",
                     makeFactoryFunction< GaussLobattoInterfaceFiniteElement< 3, 8, 7 >,
                                          FullIntegration,
                                          GaussLobattoInterfaceFiniteElement< 3, 8, 7 >::SectionType::Interface >() );

} // namespace Marmot::Elements::Registration
