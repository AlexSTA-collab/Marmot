#pragma once
#include "bftMaterial.h"

class BftMaterialMechanical : public BftMaterial {

    /*
       Abstract basic class for Mechanical materials.
       'Mechanical' is meant in the 'most general sense', i.e., any material which describes a mechanical (cauchy)
       stress - deformation relationship (e.g, hyperelastic, hypoelastic, elasto-plastic, visco-elastic materials)

       σ = f (σ, dxdX, t, .. ),

       formulated incrementally as σ_np = f (σ_n, dxdX_n, dxdX_np, Δt, t_n, .. )

       Algorithmic tangent: dσdF = d σ_np d (dxdX_np)

        Format:

        / d σ_11 d F_00,   d σ_11 d F_10,   d σ_11 d F_20,   d σ_11 d F_01, \
        | d σ_22 d F_00,                                                    |
        | d σ_33 d F_00,                                                    |
        | ...                                                               |
        | ...                                                               |
        \ ...                                                               /

        such that it can be interpreted as a column-major 6x3x3 tensor (4th order, voigt left tensor)
    */

    public:

        using BftMaterial::BftMaterial;

        virtual void computeStress( double*       stress,
                double*       dStressDDFNew,
                const double* FOld,
                const double* FNew,
                const double* timeOld,
                const double  dT,
                double&       pNewDT ) = 0;

        virtual void computePlaneStress( double*       stress,
                double*       dStressDDFNew,
                const double* FOld,
                double*       FNew,
                const double* timeOld,
                const double  dT,
                double&       pNewDT );

        virtual void computeUniaxialStress( double*       stress,
                double*       dStressDDFNew,
                const double* FOld,
                double*       FNew,
                const double* timeOld,
                const double  dT,
                double&       pNewDT );
};
