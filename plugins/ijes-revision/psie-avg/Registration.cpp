// YIQUAD4_STABP_MINI_AVG: the two-face PSIE element with the top-bottom differences neglected (mean surface gradient on both faces).
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/YStabPressureMiniAvgInterfaceFiniteElement.h"

namespace {
  using Element = Marmot::Elements::YStabPressureMiniAvgInterfaceFiniteElement< 3, 8 >;
  const static bool isRegistered = MarmotLibrary::MarmotElementFactory::registerElement(
    "YIQUAD4_STABP_MINI_AVG", []( int elementID ) -> MarmotElement* {
      return new Element( elementID, Marmot::FiniteElement::Quadrature::IntegrationTypes::FullIntegration, Element::SectionType::Interface );
    } );
}
