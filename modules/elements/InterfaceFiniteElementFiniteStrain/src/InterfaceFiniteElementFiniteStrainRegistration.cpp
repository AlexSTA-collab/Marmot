#include "Marmot/GYInterfaceFiniteElementFiniteStrain.h"
#include "Marmot/InterfaceFiniteElementFiniteStrain.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotFiniteElement.h"

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

  const static bool IQUAD8UL_isRegistered = MarmotElementFactory::
    registerElement( "IQUAD8UL",
                     makeFactoryFunction< InterfaceFiniteElementFiniteStrain< 3, 16 >,
                                          FullIntegration,
                                          InterfaceFiniteElementFiniteStrain< 3, 16 >::SectionType::Interface >() );

  const static bool YILINE2UL_isRegistered = MarmotElementFactory::
    registerElement( "YILINE2UL",
                     makeFactoryFunction< YInterfaceFiniteElementFiniteStrain< 2, 4 >,
                                          FullIntegration,
                                          YInterfaceFiniteElementFiniteStrain< 2, 4 >::SectionType::Interface >() );

  const static bool YIQUAD4UL_isRegistered = MarmotElementFactory::
    registerElement( "YIQUAD4UL",
                     makeFactoryFunction< YInterfaceFiniteElementFiniteStrain< 3, 8 >,
                                          FullIntegration,
                                          YInterfaceFiniteElementFiniteStrain< 3, 8 >::SectionType::Interface >() );

  const static bool YIQUAD8UL_isRegistered = MarmotElementFactory::
    registerElement( "YIQUAD8UL",
                     makeFactoryFunction< YInterfaceFiniteElementFiniteStrain< 3, 16 >,
                                          FullIntegration,
                                          YInterfaceFiniteElementFiniteStrain< 3, 16 >::SectionType::Interface >() );

  const static bool GYILINE2UL_isRegistered = MarmotElementFactory::
    registerElement( "GYILINE2UL",
                     makeFactoryFunction< GYInterfaceFiniteElementFiniteStrain< 2, 4 >,
                                          FullIntegration,
                                          GYInterfaceFiniteElementFiniteStrain< 2, 4 >::SectionType::Interface >() );

  const static bool GYIQUAD4UL_isRegistered = MarmotElementFactory::
    registerElement( "GYIQUAD4UL",
                     makeFactoryFunction< GYInterfaceFiniteElementFiniteStrain< 3, 8 >,
                                          FullIntegration,
                                          GYInterfaceFiniteElementFiniteStrain< 3, 8 >::SectionType::Interface >() );

  const static bool GYIQUAD8UL_isRegistered = MarmotElementFactory::
    registerElement( "GYIQUAD8UL",
                     makeFactoryFunction< GYInterfaceFiniteElementFiniteStrain< 3, 16 >,
                                          FullIntegration,
                                          GYInterfaceFiniteElementFiniteStrain< 3, 16 >::SectionType::Interface >() );

  const static bool GY2ILINE2UL_isRegistered = MarmotElementFactory::
    registerElement( "GY2ILINE2UL",
                     makeFactoryFunction< GY2InterfaceFiniteElementFiniteStrain< 2, 4 >,
                                          FullIntegration,
                                          GY2InterfaceFiniteElementFiniteStrain< 2, 4 >::SectionType::Interface >() );

  const static bool GY2IQUAD4UL_isRegistered = MarmotElementFactory::
    registerElement( "GY2IQUAD4UL",
                     makeFactoryFunction< GY2InterfaceFiniteElementFiniteStrain< 3, 8 >,
                                          FullIntegration,
                                          GY2InterfaceFiniteElementFiniteStrain< 3, 8 >::SectionType::Interface >() );

  const static bool GY2IQUAD8UL_isRegistered = MarmotElementFactory::
    registerElement( "GY2IQUAD8UL",
                     makeFactoryFunction< GY2InterfaceFiniteElementFiniteStrain< 3, 16 >,
                                          FullIntegration,
                                          GY2InterfaceFiniteElementFiniteStrain< 3, 16 >::SectionType::Interface >() );

  const static bool HILINE2UL_isRegistered = MarmotElementFactory::
    registerElement( "HILINE2UL",
                     makeFactoryFunction<
                       GradedInterfaceFiniteElementFiniteStrain< 2, 4 >,
                       FullIntegration,
                       GradedInterfaceFiniteElementFiniteStrain< 2, 4 >::SectionType::Interface >() );

  const static bool HIQUAD4UL_isRegistered = MarmotElementFactory::
    registerElement( "HIQUAD4UL",
                     makeFactoryFunction<
                       GradedInterfaceFiniteElementFiniteStrain< 3, 8 >,
                       FullIntegration,
                       GradedInterfaceFiniteElementFiniteStrain< 3, 8 >::SectionType::Interface >() );

  const static bool HIQUAD8UL_isRegistered = MarmotElementFactory::
    registerElement( "HIQUAD8UL",
                     makeFactoryFunction<
                       GradedInterfaceFiniteElementFiniteStrain< 3, 16 >,
                       FullIntegration,
                       GradedInterfaceFiniteElementFiniteStrain< 3, 16 >::SectionType::Interface >() );

} // namespace Marmot::Elements::Registration
