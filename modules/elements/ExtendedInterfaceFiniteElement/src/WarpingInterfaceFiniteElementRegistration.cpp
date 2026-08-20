#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/WarpingInterfaceFiniteElement.h"

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

  const static bool WILINE2_isRegistered = MarmotElementFactory::
    registerElement( "WILINE2",
                     makeFactoryFunction< WarpingInterfaceFiniteElement< 2, 4 >,
                                          FullIntegration,
                                          WarpingInterfaceFiniteElement< 2, 4 >::SectionType::Interface >() );

  const static bool WIQUAD4_isRegistered = MarmotElementFactory::
    registerElement( "WIQUAD4",
                     makeFactoryFunction< WarpingInterfaceFiniteElement< 3, 8 >,
                                          FullIntegration,
                                          WarpingInterfaceFiniteElement< 3, 8 >::SectionType::Interface >() );

  // Station-count variants, for the same through-thickness convergence checks
  // the Gauss-Lobatto element exposes. Four stations is the MINIMUM the warping
  // formulation admits: the three-point Lobatto rule samples xi = {-1,0,1},
  // where the antisymmetric profile phi_a = xi - xi^3 vanishes identically.
  const static bool WIQUAD4_N4_isRegistered = MarmotElementFactory::
    registerElement( "WIQUAD4_N4",
                     makeFactoryFunction< WarpingInterfaceFiniteElement< 3, 8, 4 >,
                                          FullIntegration,
                                          WarpingInterfaceFiniteElement< 3, 8, 4 >::SectionType::Interface >() );

  const static bool WIQUAD4_N7_isRegistered = MarmotElementFactory::
    registerElement( "WIQUAD4_N7",
                     makeFactoryFunction< WarpingInterfaceFiniteElement< 3, 8, 7 >,
                                          FullIntegration,
                                          WarpingInterfaceFiniteElement< 3, 8, 7 >::SectionType::Interface >() );

} // namespace Marmot::Elements::Registration
