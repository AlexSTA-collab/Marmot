# MechanicsCore

\page continuummechanics Continuum Mechanics

Theoretical Background for Continuum Mechanics.


## Additional Implementations for Material Models Based on Plasticity Theory
### Duvaut Lions Viscosity

**Implementation:** \ref DuvautLionsViscosity.h

### Menetrey Willam Yield Surfaces

**Implementation:** \ref MenetreyWillam.h


\page mechanicalmaterials Mechanical Materials

Abstract base class for mechanical materials with scalar nonlocal interaction. 

**Implementation:** \ref MarmotMaterialMechanical

\page hypoelastic Hypoelastic Material Models

## Basic Theory 

**Implementation:** \ref MarmotMaterialHypoElastic

Derived abstract base class for elastic materials expressed purely in rate form. In general, the nominal stress rate tensor \f$ \dot{\boldsymbol{\sigma}} \f$ can be written as a function of the nominal stress tensor \f$ \boldsymbol{\sigma} \f$, the stretching rate tensor \f$ \dot{\boldsymbol{\varepsilon}} \f$ and the time \f$ t \f$.

\f[
  \displaystyle \dot{\boldsymbol{\sigma}} = f( \boldsymbol{\sigma}, \dot{\boldsymbol{\varepsilon}}, t, ...)
\f]

In course of numerical time integration, this relation will be formulated incrementally as 

\f[
  \displaystyle \Delta\boldsymbol{\sigma} = f ( \boldsymbol{\sigma}_n, \Delta\boldsymbol{\varepsilon}, \Delta t, t_n, ...)
\f]

with 

\f[
  \displaystyle \Delta\boldsymbol{\varepsilon} =  \dot{\boldsymbol{\varepsilon}}\cdot \Delta t
\f]

and the algorithmic tangent 

\f[
  \displaystyle \frac{d \boldsymbol{\sigma}}{d \boldsymbol{\varepsilon}} =  \frac{d \Delta \boldsymbol{\sigma}}{d \Delta \boldsymbol{\varepsilon}}
\f]

This formulation is compatible with an Abaqus interface.



\page hyperelastic Hyperelastic Material Models

## Basic Theory 

**Implementation:** \ref MarmotMaterialHyperElastic

Derived abstract base class for _simple_, purely hyperelastic materials to be used for finite elements based on the total lagrangian kinematic description (TL elements). The second Piola - Kirchhoff stress tensor \f$ S \f$ will be derived by

\f[
  \displaystyle S = \frac{\partial f(\boldsymbol{E},t )}{\partial \boldsymbol{E}}
\f]

with the Green - Lagrange strain tensor \f$ \boldsymbol{E} \f$

\f[
  \displaystyle E  = \frac{1}{2}\,\left(\boldsymbol{F}^T\cdot \boldsymbol{F} - \boldsymbol{I} \right)
\f]

as work conjugated measure and the variable \f$ \boldsymbol{F} \f$ denoting the deformation gradient. 
The algorithmic tangent will be calculated by 

\f[
  \displaystyle \frac{d \boldsymbol{S}}{d \boldsymbol{E}}
\f]

\page gradmechanicalmaterials Gradient Enhanced Mechanical Materials

## Basic Theory 

**Implementation:** \ref MarmotMaterialGradientEnhancedMechanical

Base class for mechanical materials with gradient enhanced regularization to assure mesh indepency in finite element simulations. 

\page gradhypoelastic Gradient Enhanced Hypoelastic Material Models

**Implementation:** \ref MarmotMaterialGradientEnhancedHypoElastic
  
\page substepper Substepping Algorithms

## Adaptive Substepper

### Substepper for Implicit Return Mapping Algorithms

**Implementation:** \ref AdaptiveSubstepper.h

Adaptive Substepper, employing an error estimation and Richardson extrapolation for an implicit return mapping algorithm.

### Substepper for Semi - Explicit Return Mapping Algorithms
 
**Implementation:** \ref AdaptiveSubstepperExplicit.h

Adaptive Substepper, employing an error estimation and Richardson extrapolation for an semi - explicit return mapping algorithm.

## Perez Fouget Substepper


### Substepper for Implicit Return Mapping Algorithms

**Implementation:** \ref PerezFougetSubstepperMarkII.h

Perez Fouget Substepper for elastoplastic materials with an implicit return mapping algorithm.

### Substepper for Semi - Explicit Return Mapping Algorithms

**Implementation:** \ref PerezFougetSubstepperExplicitMarkII.h

Perez Fouget Substepper for semi - explicit elastoplastic materials.

### Modified Version to Consider Time-Variant Elastic Stiffness Tensor

**Implementation:** \ref PerezFougetSubstepperTime.h


\page hugheswinget Hughes Winget

**Implementation:** \ref HughesWinget

\page others Others

## Voigt Notation

**Implementation:** \ref MarmotVoigt.h


## Yield Surface Combination Manager

Necessary for multisurface plasticity material models.

**Implementation:** \ref YieldSurfaceCombinationManager.h
