#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/YNodalGradientInterfaceFiniteElement.h"

namespace Marmot::Elements::Registration {

  using namespace MarmotLibrary;
  using namespace Marmot::FiniteElement::Quadrature;

  namespace {
    MarmotElement* makeYNodalGradientInterfaceFiniteElement( int elementID )
    {
      return new YNodalGradientInterfaceFiniteElement( elementID,
                                                       FullIntegration,
                                                       YNodalGradientInterfaceFiniteElement::SectionType::Interface );
    }
  } // namespace

  const static bool YIQUAD4_NODALGRADIENT_isRegistered = MarmotElementFactory::
    registerElement( "YIQUAD4_NODALGRADIENT", makeYNodalGradientInterfaceFiniteElement );

} // namespace Marmot::Elements::Registration
