// C3D8B: the eight-node displacement hexahedron with the B-bar (constant-dilatation) strain-displacement matrix
// of Hughes, B = B - 1/3 m m^T B + 1/3 m <m^T B>, where <.> is the volume average over the element and m = (1 1 1 0 0 0)^T.
// It differs from C3D8 only in initializeYourself(), which rewrites the stored B of every quadrature point; the residual,
// the tangent and the strain increment all read that matrix, so nothing else changes.
#include "Marmot/DisplacementFiniteElement.h"
#include "Marmot/MarmotElementFactory.h"

namespace {
  using Base = Marmot::Elements::DisplacementFiniteElement< 3, 8 >;

  class DisplacementBBarFiniteElement : public Base {
  public:
    using Base::Base;

    void initializeYourself() override
    {
      Base::initializeYourself();
      using Row = Eigen::Matrix< double, 1, 24 >;
      Row    vbar = Row::Zero();
      double V    = 0.0;
      for ( const auto& qp : qps ) {
        vbar += ( qp.B.row( 0 ) + qp.B.row( 1 ) + qp.B.row( 2 ) ) * qp.J0xW;
        V += qp.J0xW;
      }
      vbar /= V;
      Eigen::Matrix< double, 6, 1 > m;
      m << 1, 1, 1, 0, 0, 0;
      for ( auto& qp : qps ) {
        const Row v = qp.B.row( 0 ) + qp.B.row( 1 ) + qp.B.row( 2 );
        qp.B += ( 1.0 / 3.0 ) * m * ( vbar - v );
      }
    }
  };

  const static bool C3D8B_isRegistered = MarmotLibrary::MarmotElementFactory::registerElement(
    "C3D8B", []( int elementID ) -> MarmotElement* {
      return new DisplacementBBarFiniteElement( elementID,
                                               Marmot::FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                               Base::SectionType::Solid );
    } );
} // namespace
