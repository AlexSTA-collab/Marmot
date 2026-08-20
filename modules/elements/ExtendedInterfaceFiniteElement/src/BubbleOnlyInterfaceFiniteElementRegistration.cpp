#include "Marmot/BubbleOnlyInterfaceFiniteElement.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotFiniteElement.h"

namespace Marmot::Elements::Registration {
  using namespace MarmotLibrary;
  using namespace Marmot::FiniteElement::Quadrature;

  namespace {
    MarmotElement* makeBubbleOnly( int id )
    {
      return new BubbleOnlyInterfaceFiniteElement( id, FullIntegration, BubbleOnlyInterfaceFiniteElement::Interface );
    }
  } // namespace

  const static bool GLIQUAD4_BUBBLE_isRegistered = MarmotElementFactory::registerElement( "GLIQUAD4_BUBBLE",
                                                                                          makeBubbleOnly );
} // namespace Marmot::Elements::Registration
