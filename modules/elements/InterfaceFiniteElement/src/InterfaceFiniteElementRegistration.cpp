#include "Marmot/InterfaceFiniteElement.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/YInterfaceFiniteElement.h"

namespace Marmot::Elements::Registration {

  /**
   * @brief Create the factory function of an element type.
   *
   * @tparam T element type
   * @tparam integrationType integration type of the created elements
   * @tparam sectionType section type of the created elements
   * @return a function that creates an element of type T from an element number
   */
  template < class T,
             Marmot::FiniteElement::Quadrature::IntegrationTypes integrationType,
             typename T::SectionType                             sectionType >
  MarmotLibrary::MarmotElementFactory::elementFactoryFunction makeFactoryFunction()
  {
    return []( int elementID ) -> MarmotElement* { return new T( elementID, integrationType, sectionType ); };
  }

  using namespace MarmotLibrary;
  using namespace Marmot::FiniteElement::Quadrature;

  const static bool ILINE2_isRegistered = MarmotElementFactory::
    registerElement( "ILINE2",
                     makeFactoryFunction< InterfaceFiniteElement< 2, 4 >,
                                          FullIntegration,
                                          InterfaceFiniteElement< 2, 4 >::SectionType::Interface >() );

  const static bool IQUAD4_isRegistered = MarmotElementFactory::
    registerElement( "IQUAD4",
                     makeFactoryFunction< InterfaceFiniteElement< 3, 8 >,
                                          FullIntegration,
                                          InterfaceFiniteElement< 3, 8 >::SectionType::Interface >() );

  const static bool IQUAD8_isRegistered = MarmotElementFactory::
    registerElement( "IQUAD8",
                     makeFactoryFunction< InterfaceFiniteElement< 3, 16 >,
                                          FullIntegration,
                                          InterfaceFiniteElement< 3, 16 >::SectionType::Interface >() );

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

  const static bool YIQUAD8_isRegistered = MarmotElementFactory::
    registerElement( "YIQUAD8",
                     makeFactoryFunction< YInterfaceFiniteElement< 3, 16 >,
                                          FullIntegration,
                                          YInterfaceFiniteElement< 3, 16 >::SectionType::Interface >() );

  const static bool HILINE2_isRegistered = MarmotElementFactory::
    registerElement( "HILINE2",
                     makeFactoryFunction< GradedInterfaceFiniteElement< 2, 4 >,
                                          FullIntegration,
                                          GradedInterfaceFiniteElement< 2, 4 >::SectionType::Interface >() );

  const static bool HIQUAD4_isRegistered = MarmotElementFactory::
    registerElement( "HIQUAD4",
                     makeFactoryFunction< GradedInterfaceFiniteElement< 3, 8 >,
                                          FullIntegration,
                                          GradedInterfaceFiniteElement< 3, 8 >::SectionType::Interface >() );

  const static bool HIQUAD8_isRegistered = MarmotElementFactory::
    registerElement( "HIQUAD8",
                     makeFactoryFunction< GradedInterfaceFiniteElement< 3, 16 >,
                                          FullIntegration,
                                          GradedInterfaceFiniteElement< 3, 16 >::SectionType::Interface >() );

} // namespace Marmot::Elements::Registration