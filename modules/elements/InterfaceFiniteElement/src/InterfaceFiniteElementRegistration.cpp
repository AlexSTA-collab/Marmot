#include "Marmot/InterfaceFiniteElement.h"
#include "Marmot/Marmot.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotFiniteElementSpatialWrapper.h"

namespace Marmot::Elements::InterfaceRegistration {

  enum InterfaceElementCode {

    // Plane stress 2D - Line 2
    ILINE2_S  = 10402,
    ILINE2R_S = 10805,

    // Plane Strain 2D
    ILINE2  = 10407,
    ILINE3R = 10808,
    ILINE3  = 10807,

    // Solid
    IQUAD4   = 10803,
    IQUAD4_R = 10806,
    IQUAD8   = 12003,
    IQUAD8_R = 12006
  };

  template < class T,
             Marmot::FiniteElement::Quadrature::IntegrationTypes integrationType,
             typename T::SectionType                             sectionType >
  MarmotLibrary::MarmotElementFactory::elementFactoryFunction makeFactoryFunction()
  {
    return []( int elementID ) -> MarmotElement* { return new T( elementID, integrationType, sectionType ); };
  }

  using namespace MarmotLibrary;
  using namespace Marmot::FiniteElement::Quadrature;

  const static bool ILine2S_isRegistered = MarmotElementFactory::
    registerElement( "ILINE2S",
                     InterfaceElementCode::ILINE2_S,
                     makeFactoryFunction< InterfaceFiniteElement< 2, 4 >,
                                          FullIntegration,
                                          InterfaceFiniteElement< 2, 4 >::PlaneStress >() );

  const static bool ILine2_isRegistered = MarmotElementFactory::
    registerElement( "ILINE2",
                     InterfaceElementCode::ILINE2,
                     makeFactoryFunction< InterfaceFiniteElement< 2, 4 >,
                                          FullIntegration,
                                          InterfaceFiniteElement< 2, 4 >::PlaneStrain >() );

  const static bool ILine3S_isRegistered = MarmotElementFactory::
    registerElement( "ILINE3S",
                     InterfaceElementCode::ILINE3R,
                     makeFactoryFunction< InterfaceFiniteElement< 2, 6 >,
                                          ReducedIntegration,
                                          InterfaceFiniteElement< 2, 6 >::PlaneStress >() );

  const static bool ILine3_isRegistered = MarmotElementFactory::
    registerElement( "ILINE3",
                     InterfaceElementCode::ILINE3,
                     makeFactoryFunction< InterfaceFiniteElement< 2, 6 >,
                                          FullIntegration,
                                          InterfaceFiniteElement< 2, 6 >::PlaneStrain >() );

  const static bool IQuad4_isRegistered = MarmotLibrary::MarmotElementFactory::
    registerElement( "IQUAD4", InterfaceElementCode::IQUAD4, []( int elementID ) -> MarmotElement* {
      return new InterfaceFiniteElement< 3, 8 >( elementID,
                                                 Marmot::FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                                 InterfaceFiniteElement< 3, 8 >::SectionType::Solid );
    } );

  const static bool IQuad8_isRegistered = MarmotLibrary::MarmotElementFactory::
    registerElement( "IQUAD8", InterfaceElementCode::IQUAD8, []( int elementID ) -> MarmotElement* {
      return new InterfaceFiniteElement< 3, 16 >( elementID,
                                                  Marmot::FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                                  InterfaceFiniteElement< 3, 16 >::SectionType::Solid );
    } );

  const static bool IQuad8R_isRegistered = MarmotLibrary::MarmotElementFactory::
    registerElement( "IQUAD8R", InterfaceElementCode::IQUAD8_R, []( int elementID ) -> MarmotElement* {
      return new InterfaceFiniteElement< 3,
                                         16 >( elementID,
                                               Marmot::FiniteElement::Quadrature::IntegrationTypes::ReducedIntegration,
                                               InterfaceFiniteElement< 3, 16 >::SectionType::Solid );
    } );

} // namespace Marmot::Elements::InterfaceRegistration
