// Small-strain hypoelastic Drucker-Prager plasticity with non-associated flow and linear hardening of the cohesion,
// for the interface elements (through provider=marmotmaterialhypoelasticinterface) and the resolved layer (C3D8, C3D8B).
//
// de Souza Neto, Peric & Owen (2008), Computational Methods for Plasticity, Sec. 8.3.3 and Box 8.9, with the cone
// matched to Mohr-Coulomb in plane strain (their Eq. 6.122). Stress positive in tension, p = tr(sigma) / 3:
//
//   yield      f = sqrt(J2) + eta    p - xi c(kappa)        eta    = 3 tan(phi) / sqrt(9 + 12 tan^2(phi))
//   potential  g = sqrt(J2) + etaBar p                       etaBar = 3 tan(psi) / sqrt(9 + 12 tan^2(psi))
//   hardening  c(kappa) = c0 + H kappa,  kappa' = xi lambda'  xi     = 3 / sqrt(9 + 12 tan^2(phi))
//
// psi != phi gives non-associated flow and a non-symmetric tangent; phi = psi = 0 is J2 plasticity with
// sigma_y = sqrt(3) c0 and H_J2 = 3 H, and phi = 0 < psi the non-associated von Mises model with dilatant flow.
// psi = 0 < phi has no apex return (the flow has no volumetric part): beyond the apex a tension cut-off puts the
// stress on the apex of the current cone, kappa unchanged, with a zero tangent.
//
// Optional tension cap (7th property pt, below the apex xi c0 / eta): f2 = p - pt <= 0 with associated, purely
// volumetric flow. Three plastic regimes, each with a closed-form return: cone, cap (p = pt, deviator elastic) and
// the corner of both (p = pt, radial return of the deviator on the cone, kappa hardens). The cap removes the apex, and
// its tangent keeps the shear stiffness. Without pt the apex return (psi > 0) or the cut-off (psi = 0) applies.
//
// Properties: E, nu, c0, phi [deg], psi [deg], H [, pt [, density]].   State: kappa.
#pragma once
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotStateHelpers.h"

namespace Marmot::Materials {

  class DruckerPragerModel : public MarmotMaterialHypoElastic {

  public:
    DruckerPragerModel( const double* materialProperties, const int nMaterialProperties, const int materialLabel );

    void computeStress( state3D&                state,
                        Marmot::Matrix6d&       dStressDDStrain,
                        const Marmot::Vector6d& dStrain,
                        const timeInfo&         timeInfo ) const override;

    double getDensity( const double* stateVars ) const override;

    // cone parameters of a friction (or dilatancy) angle in degrees, plane-strain match
    static double eta( double angleInDegrees );
    static double xi( double angleInDegrees );
  };

} // namespace Marmot::Materials
