# Plug-ins of the IJES revision

Elements and materials used in the analyses of the revised manuscript *A generalised interface model for
inelastic interphases* (IJES), built as shared libraries against an installed Marmot of this branch and
loaded into EdelweissFE with `LD_PRELOAD`.

| directory | registers | used for |
|---|---|---|
| `bbar` | `C3D8B` | resolved $J_2$ layer: hexahedron with the $\bar B$ treatment of the volumetric strain |
| `mixed` | `C3D8P` | mixed displacement-pressure hexahedron (locking checks) |
| `druckerprager` | `DRUCKERPRAGER` | Drucker--Prager material with tension cap, linear hardening, closed-form return on every branch |
| `psie-avg` | `YIQUAD4_STABP_MINI_AVG` | PSIE variant with the averaged pressure |
| `psie-squeeze` | `YIQUAD4_STABP_MINI_SQ` | PSIE with the squeeze-flow term (Section 4.4.4 of the manuscript) |
| `twophase` | `YIQUAD4_TWOPHASE` | EIE of a biphasic viscoelastic layer, each face with its own phase (Section 4.3.3) |

`druckerprager/test_dp.cpp` checks the return mapping and its consistent tangent against finite differences.

Build: `./build.sh <Marmot install prefix> <conda env prefix> [out dir]` (C++20, OpenMP; Eigen and autodiff from the
environment). The Drucker--Prager card is `E, nu, c0, phi, psi, H, p_t`; as interface material
(`provider=marmotmaterialhypoelasticinterface`) the normal thickness follows nu: `E, nu, h, c0, phi, psi, H, p_t`.
