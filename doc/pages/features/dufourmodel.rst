Dufour model
============

A viscoelastic-viscoplastic damage model for adhesive layers at finite strain, written against
``MarmotMaterialGradientEnhancedFiniteStrain``. It is registered as ``DUFOURMODEL``.

Theory
------

The Kirchhoff stress is :math:`\boldsymbol\tau = (1 - \omega)\,\boldsymbol\tau_\mathrm{eff}`. The effective stress follows
from a Pence-Gou potential (variant B) for the elastic part of a multiplicative split of the deformation gradient,
:math:`\boldsymbol F = \boldsymbol F^e \boldsymbol F^p`, optionally relaxed by a generalized Maxwell chain: the deviatoric
and the volumetric part of the second Piola-Kirchhoff stress have their own relative moduli and share the relaxation
times. The plastic flow is viscoplastic and non-associative. The yield function is pressure sensitive with the ratio of
the compressive to the tensile strength, and the flow potential depends on the plastic Poisson ratios in tension and in
compression. The return map is solved for the elastic deformation gradient, the plastic strain
:math:`\alpha_P` and the plastic multiplier.

Damage
^^^^^^

The damage variable has a softening and a failure part,

.. math::

   \omega = 1 - (1 - \omega_s)(1 - \omega_f), \qquad
   \omega_s = 1 - \exp\left(-\frac{\bar\alpha_P}{\varepsilon_F}\right), \qquad
   \bar\alpha_P = m\,\bar N + (1 - m)\,\alpha_P ,

where :math:`\bar N` is the nonlocal field and :math:`m` the weight of its part. The failure part is zero until a driver
:math:`D` reaches one, and then grows as :math:`\omega_f = 1 - \exp(-(D - 1)/d_F)`. The driver accumulates the plastic
strain with a stress-dependent weight, after Rice and Tracey and Smith et al.,

.. math::

   D = c_{SW} \int \left[ \mathrm{e}^{E(T)} - \frac{1}{b_{SW}}\,\mathrm{e}^{-E(T)} \right] w(\dot\kappa)\,
       \mathrm{d}\bar\alpha_P ,

with the stress triaxiality :math:`T` (see Abrari Vajari et al., CMAME 400 (2022) 115467). The exponent is
:math:`E'(T) = (b_0 + b_1 T + b_2 T^2)^2` for positive and :math:`E = s_0 T` for negative triaxialities, with the
Rice-Tracey slope :math:`s_0`, so that :math:`E` cannot decrease with the triaxiality. The factor
:math:`w = (\dot\kappa/\dot\kappa_\mathrm{ref})^{-s}` makes the driver rate dependent, and a Gaussian bump on
:math:`g = \mathrm{e}^{E}` and a cap on the triaxiality are optional. The triaxiality, the Lode parameter and the rate
factor are held fixed in the tangent.

For :math:`c_{SW} = 0` the failure part is absent. The growth of the driver in one increment is limited: a larger
increase makes the material throw ``Marmot::StressUpdateFailed``, so that the solver retakes the increment with a
smaller time step.

Energy and dissipation
^^^^^^^^^^^^^^^^^^^^^^

The elastic energy density of the response is :math:`(1 - \omega)\,\psi_\mathrm{eff}`. The dissipation is cumulative, as in
the other damage models: the incoming value is incremented by
:math:`(1 - \omega)\,\boldsymbol{M} : \Delta\boldsymbol{\varepsilon}^p + \psi_\mathrm{eff}\,\Delta\omega`, with the effective
Mandel stress :math:`\boldsymbol{M}` and the plastic increment of the return map. A decrease of the damage, which the nonlocal
field can cause, does not dissipate. Without Maxwell branches the work of the stress along a path is the stored energy plus
the dissipation. The viscous dissipation of the Maxwell branches is not included, because the energy stored in them is not
tracked.

Material parameters
-------------------

The card is positional.

.. list-table::
   :header-rows: 1

   * - Index
     - Description
   * - 0, 1
     - bulk and shear modulus :math:`K`, :math:`G`
   * - 2, 3
     - tensile and compressive strength :math:`f_t`, :math:`f_c`
   * - 4 to 10
     - hardening of the yield strength :math:`Q_1, Q_2, b_1 \dots b_5`
   * - 11, 12
     - viscosity :math:`\eta_{VP}` and rate exponent :math:`n`
   * - 13, 14
     - plastic Poisson ratios in tension and compression
   * - 15, 16
     - softening strain :math:`\varepsilon_F` and largest damage
   * - 17, 18
     - nonlocal radius :math:`l_d` and weight :math:`m` of the nonlocal field
   * - 19
     - density, optional
   * - 20, 21, 22
     - :math:`c_{SW}`, :math:`b_{SW}` and :math:`d_F`, optional
   * - 23
     - number of Maxwell branches, from 0 to 7
   * - 24 + 3i, 25 + 3i, 26 + 3i
     - deviatoric and volumetric relative modulus and relaxation time of branch i
   * - then
     - :math:`b_0, b_1, b_2`, the reference rate, the rate exponent, the cap of the triaxiality, and the amplitude, the
       centre and the width of the bump

A card of an older layout, or with entries beyond the layout, is refused.

State variables
---------------

``Fp``, ``alphaP``, ``omega``, ``damageDriver``, ``alphaPBar``, ``PK2Ref``, ``veDev``, ``veVol``, and the output
variables ``triax``, ``omegaS``, ``omegaF`` and ``lode``.

Use in an interface element
---------------------------

In the interface elements, the card is ``[K, G, h, ...]`` with the thickness ``h`` of the layer inserted after the first
two entries, see :ref:`equilibratedinterfacematerialfinitestrain`. The elements ``GYIQUAD4UL``, ``GYIQUAD8UL`` and
``GYILINE2UL`` carry the nonlocal field of the model. With the weight :math:`m = 0` the model is local, and it can also
be used in the elements without a nonlocal field.

Implementation
--------------

.. doxygenclass:: Marmot::Materials::DufourModel
   :members:
