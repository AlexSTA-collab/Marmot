// YIQUAD4_TWOPHASE: the extended interface element YIQUAD4 whose two faces carry their own complete property sets,
// see MarmotEquilibratedXTwoPhaseInterfaceMaterialHypoElastic.
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/YInterfaceFiniteElement.h"
#include "Marmot/MarmotEquilibratedXTwoPhaseInterfaceMaterialHypoElastic.h"

namespace {
  using Element = Marmot::Elements::InterfaceFiniteElementBase<
    3, 8,
    Marmot::Elements::EquilibratedTwoSidedFormulation< 3, 8, MarmotEquilibratedXTwoPhaseInterfaceMaterialHypoElastic > >;
  const static bool isRegistered = MarmotLibrary::MarmotElementFactory::registerElement(
    "YIQUAD4_TWOPHASE", []( int elementID ) -> MarmotElement* {
      return new Element( elementID, Marmot::FiniteElement::Quadrature::IntegrationTypes::FullIntegration, Element::SectionType::Interface );
    } );
}
