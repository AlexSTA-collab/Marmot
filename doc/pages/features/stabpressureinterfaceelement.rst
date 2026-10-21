.. _stabpressureinterfaceelement:

Interface element with a mixed pressure (``YIQUAD4_STABP_MINI``)
================================================================

Why a mixed pressure
--------------------

With a vanishing hardening modulus a von Mises layer reaches the incompressible limit,
:math:`K / \mu_\mathrm{p} \sim 10^{7}`. A displacement formulation evaluates the pressure
:math:`p = K \operatorname{tr}\boldsymbol{\varepsilon}` from a trace that is driven to zero, and so multiplies the
discretisation error by :math:`K`. The result is an element-to-element checkerboard of the pressure while the
deviatoric stress is uniform. A mean-dilatation fix inside the element cannot remove it, because the checkerboard is
an element-to-element mode.

Formulation
-----------

The pressure is an independent, :math:`C^0`-continuous field with one value per node of the midsurface, and the
equilibrium of the displacements uses

.. math::

   \tilde{\boldsymbol{\sigma}} = \operatorname{dev} \boldsymbol{\sigma}(\boldsymbol{\varepsilon}) - p_h \mathbf{I},
   \qquad p_h = \mathbf{N} \hat{\mathbf{p}} .

The volumetric law is enforced weakly, together with a Brezzi-Pitkaranta stabilisation of the pressure gradient,

.. math::

   R_u = \mathbf{W}^{T} \mathbf{f} + \mathbf{B}^{+T} \mathbf{S}^{+} + \mathbf{B}^{-T} \mathbf{S}^{-}, \qquad
   R_p = h \left[ \mathbf{N}^{T} r_p + \frac{\gamma h_e^{2}}{2\mu} \nabla \mathbf{N}^{T} \nabla \mathbf{N}\, \hat{\mathbf{p}} \right],

with :math:`r_p = \langle \operatorname{tr}\Delta\boldsymbol{\varepsilon} \rangle + \Delta p / K`. The stabilisation scale
is the mesh size :math:`h_e` and the term vanishes with it, so it does not change the converged solution. The unknowns
of the element are the 24 displacements and the 4 pressures.

The displacement of the midsurface is enriched with the bubble :math:`b = (1 - \xi^{2})(1 - \eta^{2})` and a vector
amplitude :math:`\boldsymbol{\beta}`. The amplitude enters the average surface gradient of both faces equally, and it
is condensed at every call by an inner Newton iteration with a backtracking line search, because it enters the
return map.

The material state is committed by every evaluation of the material. The element therefore restarts every trial of the
bubble iteration from the state at the beginning of the increment, and undoes a rejected trial before the next one.

Material
--------

``MarmotStabPressureXInterfaceMaterialHypoElastic`` provides the generalized stress from :math:`\tilde{\boldsymbol{\sigma}}`,
the volumetric residual and the tangent, which has purely geometric blocks for the pressure. Replacing the volumetric
part is admissible only for a yield surface that is independent of the pressure, as the one of von Mises.

Parameters
----------

The first element property is the thickness of the element, the second the coefficient :math:`\gamma` of the
stabilisation (default 0.2). ``MARMOT_STABP_GAMMA`` overrides the default in sweeps, and ``MARMOT_MINI_BUBBLE=OFF``
holds the bubble at zero, which gives the plain equal-order displacement-pressure pairing for ablation studies.
