// YIQUAD4_STABP_MINI_SQ: the PSIE element with the squeeze-flow term in the pressure equation (see the header).
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/YStabPressureMiniSqueezeInterfaceFiniteElement.h"

namespace {
  using Element = Marmot::Elements::YStabPressureMiniSqueezeInterfaceFiniteElement< 3, 8 >;
  const static bool isRegistered = MarmotLibrary::MarmotElementFactory::registerElement(
    "YIQUAD4_STABP_MINI_SQ", []( int elementID ) -> MarmotElement* {
      return new Element( elementID, Marmot::FiniteElement::Quadrature::IntegrationTypes::FullIntegration, Element::SectionType::Interface );
    } );
}
