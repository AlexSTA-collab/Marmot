.. _gradedinterfacematerial:

Graded Interface Materials
==========================

Theory
------

The interface materials of this page describe a layer of finite thickness :math:`h`, whose stiffness
:math:`\mathbb{C}(y)` varies through the thickness, with the coordinate :math:`y \in [0, 1]` from the bottom face to the
top face along the normal :math:`\mathbf{n}`. The stiffness may be given by layers of constant properties, or by
properties that vary continuously, and any bulk material of the library may be used.

The displacement of the layer is :math:`\mathbf{u}^0 + h\,\mathbf{u}^1`, with the macroscopic displacement
:math:`\mathbf{u}^0`, which interpolates the faces, and the displacement :math:`\mathbf{u}^1` of the microproblem,
which vanishes on both faces. With the normal projection :math:`\mathbf{N} = \mathbf{n}\otimes\mathbf{n}` and the
tangent projection :math:`\mathbf{T} = \mathbf{I} - \mathbf{N}`, the gradient of the layer at the height :math:`y` is

.. math::

   \mathbf{G}(y) = \underbrace{\frac{1}{\ell}\, [\![\mathbf{u}^0]\!] \otimes \mathbf{n} + (1-y)\,\mathbf{A}^- + y\,\mathbf{A}^+
   + \mathbf{B}_\tau \bar{\mathbf{A}}}_{\mathbf{u}^0} + \underbrace{\mathbf{w}(y) \otimes \mathbf{n}}_{u^1_{k|r} N_{rl}}
   + \underbrace{h\,\nabla_T \mathbf{u}^1(y)}_{h\,u^1_{k,r} T_{rl}}, \qquad \mathbf{w} = \partial_y \mathbf{u}^1,

with the jump :math:`[\![\mathbf{u}^0]\!]` of the faces, their surface gradients :math:`\mathbf{A}^\pm`, the normal
separation :math:`\ell` of the faces, and the correction :math:`\mathbf{B}_\tau \bar{\mathbf{A}}` of their tangential
separation. For finite strains the gradient is the one of the deformation gradient
:math:`\mathbf{F} = \mathbf{I} + \mathbf{G}` in the reference configuration.

The variational form is the stationarity of the energy :math:`\int_\Gamma h \int_0^1 \psi(\mathbf{G})\,\mathrm{d}y\,
\mathrm{d}\Gamma` for the test functions of both displacements,

.. math::

   \int_\Gamma h \int_0^1 \boldsymbol{\sigma}(\mathbf{G}) : \hat{\mathbf{G}} \,\mathrm{d}y\,\mathrm{d}\Gamma = 0, \qquad
   \hat{\mathbf{G}} = \mathbf{G}^0(\hat{\mathbf{u}}^0) + \partial_y \hat{\mathbf{u}}^1 \otimes \mathbf{n}
   + h\,\nabla_T \hat{\mathbf{u}}^1 .

* **Macro problem:** the variation of :math:`\mathbf{u}^0` with :math:`\hat{\mathbf{u}}^1 = \mathbf{0}`.
* **Microproblem:** the variation of :math:`\mathbf{u}^1` with :math:`\hat{\mathbf{u}}^0 = \mathbf{0}`. It is
  loaded by the macroscopic gradient, and its strong form is the equilibrium of the **total** traction
  :math:`\mathbf{t} = \mathbf{t}^0 + \mathbf{t}^1 = \boldsymbol{\sigma}(\mathbf{G}^0 + \mathbf{G}^1)\,\mathbf{n}` (small
  strains) or :math:`\mathbf{P}\,\mathbf{N}` (finite strains) through the thickness,

  .. math::

     \partial_y ( \boldsymbol{\sigma}\,\mathbf{n} ) + h\, \nabla_\Gamma \cdot ( \boldsymbol{\sigma}\,\mathbf{T} ) = \mathbf{0},
     \qquad \int_0^1 \mathbf{w}\,\mathrm{d}y = \mathbf{0}.

  The traction :math:`\mathbf{t}^0` of the macroscopic gradient alone is not constrained.

For two integration points at the faces, and without the terms of :math:`h\,\nabla_T \mathbf{u}^1`, the layer is the
equilibrated two-sided interface material.

The microproblem is solved at every stress update. With the weights :math:`\omega_q` and the acoustic tensors
:math:`\mathbf{Q}_q = \mathbf{n}\cdot\mathbb{C}_q\cdot\mathbf{n}` (finite strains: :math:`\mathbb{A}_q[\mathbf{N},
\mathbf{N}]`) of the integration points through the thickness, a Newton step for the corrections
:math:`\mathbf{w}_q` and the common traction :math:`\mathbf{t}` is

.. math::

   \delta \mathbf{t} = \Big( \sum_q \omega_q \mathbf{Q}_q^{-1} \Big)^{-1}
     \Big( \sum_q \omega_q \mathbf{Q}_q^{-1} ( \mathbf{a}_q - \mathbf{t} ) - \sum_q \omega_q \mathbf{w}_q \Big),
   \qquad
   \delta \mathbf{w}_q = \mathbf{Q}_q^{-1} \big( \delta \mathbf{t} - ( \mathbf{a}_q - \mathbf{t} ) \big),

where :math:`\mathbf{a}_q` is the traction of the point, including the change of the traction through the thickness
that the tangential gradient causes. A linear material converges in one step. The step is damped by a line search on
the norm of the residual, because a yielding point has a small tangent along its flow direction.

The tangential gradient of u¹
-----------------------------

The microproblem is solved locally at every surface point of the element, with the macroscopic gradient as its
loading, so that :math:`\mathbf{u}^1` is a function of the macroscopic gradient. Its tangential gradient follows from
the variation of this gradient along the surface, which the faces give: the jump varies along the surface with
:math:`\nabla_T [\![\mathbf{u}]\!] = \mathbf{A}^+ - \mathbf{A}^-`, so that

.. math::

   \partial_r \boldsymbol{\varepsilon} = \frac{1}{\ell}\, \mathrm{sym}\big( ( \mathbf{A}^+ - \mathbf{A}^- )\, \mathbf{e}_r
   \otimes \mathbf{n} \big), \qquad
   \partial_r \mathbf{w}_q = \frac{\partial \mathbf{w}_q}{\partial \boldsymbol{\varepsilon}}\, \partial_r \boldsymbol{\varepsilon},
   \qquad
   \mathbf{u}^1(y) = \int_0^y \mathbf{w}\,\mathrm{d}y,

with the derivative of the correction with respect to the strain of the layer, which the condensation of the
microproblem provides. This enters the two equations of the layer, both in proportion to :math:`\mathbf{A}^+ -
\mathbf{A}^-`:

* the strain of the point gains the tangential gradient :math:`h\,\nabla_T \mathbf{u}^1(y_q)`, which is the term
  :math:`h\,u^1_{k,r} T_{rl}` of the form, and
* the traction changes through the thickness by the tangential divergence of the stress, which is the same term of
  the test function integrated by parts along the surface.

The derivative of the correction is the one of the converged state of the increment, which is iterated until it does
not change. The second derivatives of the displacement of the faces (their curvature) are neglected.

Generalized stress and tangent
------------------------------

The generalized stress is the derivative of the energy of the layer with respect to the generalized strain
:math:`\mathbf{x} = ( [\![\mathbf{u}]\!], \mathbf{A}^+, \mathbf{A}^- )`,

.. math::

   \mathbf{p} = h \sum_q \omega_q\, \hat{\mathbf{B}}_q^\mathsf{T} \boldsymbol{\sigma}_q, \qquad
   \hat{\mathbf{B}}_q = \frac{\partial \mathbf{G}(y_q)}{\partial \mathbf{x}},

where the derivative contains, besides the direct dependence, the response :math:`\mathbf{w}(\mathbf{x})` of the
microproblem and the tangential gradient of the correction. The test function is therefore built from the same
operator as the strain, and the tangent

.. math::

   \mathbf{K} = h \sum_q \omega_q\, \hat{\mathbf{B}}_q^\mathsf{T} \mathbb{C}_q \hat{\mathbf{B}}_q

is symmetric for an elastic layer. The resultants of the stress through the thickness (the mean and the moments
:math:`h\int y\,\boldsymbol{\sigma}\,\mathrm{d}y`) are available as ``stressResultant()`` of the materials. They differ
from the generalized stress by a divergence along the surface, which does not change the solution of the assembled
problem.

Without the term of :math:`h\,\nabla_T \mathbf{u}^1` the work of the traction on the response of the microproblem
vanishes, because the traction is the same at all points, and the generalized stress reduces to the resultants.

Validation
----------

A strip of two elastic layers (:math:`E = 4000, \nu = 0.3` for 40 % of the thickness, :math:`E = 20000, \nu = 0.2`
for the rest, :math:`h = 0.02`) with the displacement of the top face :math:`u_z = z_0 + g x`, :math:`u_x = 0`, and the
bottom face fixed, has an exact polynomial solution of plane strain, derived with the computer algebra system sympy.
Its shear tractions on the faces are independent of :math:`x`, and a patch of elements reproduces them for every
element width:

.. list-table::
   :header-rows: 1

   * - Variant
     - Shear traction, bottom face
     - Shear traction, top face
   * - Exact solution
     - 6.3971
     - 0.05013
   * - ``HY*`` with the tangential gradient of u¹
     - 6.3971
     - 0.05013
   * - ``HY*`` without it (flag 0)
     - 6.609
     - 0.2617

The finite-strain element ``HY*UL`` reproduces the same solution for a strain of :math:`10^{-5}` times it, and the
elements ``GYH*`` the one of a homogeneous strip. The tangent of the elastic stack is symmetric to the precision of
the arithmetic, also for an assembled patch.

Material properties
--------------------

The properties of all graded materials are :math:`[h, \mathrm{mode}, n_G, K, m, (p_k, \text{bulk properties})\times K]`.

.. list-table::
   :header-rows: 1
   :widths: 10 20 70

   * - Index
     - Name
     - Meaning
   * - 0
     - ``h``
     - Thickness of the layer
   * - 1
     - ``mode``
     - 0: layers of constant properties, 1: properties interpolated linearly between stations, 2: integration points at
       the stations (trapezoidal rule)
   * - 2
     - ``nGauss``
     - Gauss points per layer (mode 0) or per interval between stations (mode 1), not used in mode 2
   * - 3
     - ``K``
     - Number of layers (mode 0) or of stations (modes 1 and 2, at least two)
   * - 4
     - ``m``
     - Number of properties of the bulk material

These are followed by ``K`` groups, each of the entry :math:`p_k` and the ``m`` properties of the bulk material.
In mode 0, :math:`p_k` is the relative thickness of the layer :math:`k`, counted from the bottom, normalized by the
sum of all of them. In modes 1 and 2, :math:`p_k` is the position of the station :math:`k`, with :math:`p_0 = 0`,
:math:`p_{K-1} = 1` and increasing positions in between.

One more property may follow the groups: the flag of the term :math:`h\,\nabla_T \mathbf{u}^1` of the variational
form. The term is kept when the flag is absent or not zero, and switched off when it is zero.

The materials
-------------

``MarmotGradedXInterfaceMaterialHypoElastic``
   Small-strain material for the elements ``HYLINE2``, ``HYQUAD4`` and ``HYQUAD8``, with the contract of
   ``MarmotEquilibratedXInterfaceMaterialHypoElastic``: the surface gradients of both faces are given, and a surface
   stress is returned for each of them. The bulk material is a hypoelastic material.

``MarmotGradedEquilibratedInterfaceMaterialFiniteStrain``
   Finite-strain material for the elements ``HYLINE2UL``, ``HYQUAD4UL`` and ``HYQUAD8UL``, with the contract of
   ``MarmotEquilibratedInterfaceMaterialFiniteStrain``, which receives the total kinematics. The bulk material is a
   finite-strain material.

``MarmotGradedEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain``
   Gradient-enhanced finite-strain material for the elements ``GYHnILINE2UL``, ``GYHnIQUAD4UL`` and
   ``GYHnIQUAD8UL`` with :math:`n = 2, 3, 4` stations. Its bulk material is gradient-enhanced (for example
   ``DUFOURMODEL``), the profile is in mode 2, and every station has a nonlocal field of its own, so that the damage
   varies through the thickness. The nonlocal equation of the station :math:`s` is, reduced to the midsurface,

   .. math::

      h\,\omega_s \Big[ \mathbf{N}^\mathsf{T} ( \bar N_s - L_s ) + l_s^2\, \nabla\mathbf{N}^\mathsf{T} \nabla \bar N_s \Big]
      + \frac{l_s^2}{h} \Big[ \frac{\mathbf{N}^\mathsf{T}( \bar N_s - \bar N_{s+1} )}{y_{s+1} - y_s}
      + \frac{\mathbf{N}^\mathsf{T}( \bar N_s - \bar N_{s-1} )}{y_s - y_{s-1}} \Big] = 0,

   where :math:`L_s` is the local driving field of the station, :math:`l_s` its nonlocal radius, and the last term the
   derivative of the field through the thickness. The momentum equation has the generalized stress of the layer, and
   the response of the microproblem depends on the nonlocal fields through the tangent of the points, which the
   tangent of the element includes. For two stations and without the term of :math:`h\,\nabla_T \mathbf{u}^1`, the
   element is ``GY2`` to the residual of the displacements and the full tangent, for the damage model and for the
   Dufour model, which is the way to test the Dufour model through a graded layer exactly. The fields of the bottom
   station and of the inner stations are degrees of freedom of the nodes of the bottom side, the field of the top
   station of the nodes of the top side.

All materials solve the microproblem and the condensation with the same solver
(``Marmot::Materials::TwoSidedGradedLayer``) for small and finite strains, with and without damage.

Limitations
-----------

* The tangential gradient of :math:`\mathbf{u}^1` neglects the second derivatives of the displacement of the faces,
  and the tangent neglects the dependence of the derivative of the correction on the state. The tangent is therefore
  exact for an elastic layer and without the term of :math:`h\,\nabla_T \mathbf{u}^1`, and quasi-Newton otherwise.
  With a damage model the dependence on the nonlocal field of the station is included, the remaining error of the
  tangent with respect to the fields is a few per cent when the faces have very different surface gradients.
* The equilibrated two-point material integrates the surface gradient, which is linear through the thickness, with
  the trapezoidal rule. The contribution of :math:`\mathbf{A}^+ - \mathbf{A}^-` to the moments of a homogeneous
  layer is then :math:`h/4` times the stiffness instead of the :math:`h/12` of the exact integral. The materials with
  Gauss points (modes 0 and 1) have the exact one, while the points at the stations (mode 2, required by the
  gradient-enhanced material) have the trapezoidal rule, as the faces have.
* Where the bulk is nearly perfectly plastic at all points that carry the traction, the acoustic tensor of each of
  them loses its stiffness along the flow direction, and the microproblem has no unique solution. The material then
  throws ``StressUpdateFailed``, which makes the element ask for a smaller time step.
* With softening, the strain localizes in one integration point. The dissipated energy then depends on the weight of
  that point, so on the number of integration points, unless the nonlocal field regularizes it, as in the
  gradient-enhanced material.
* The derivative of the nonlocal radius is not part of the tangents.

Implementation
--------------

.. doxygenstruct:: Marmot::Materials::ThicknessProfile
   :members:

.. doxygenclass:: MarmotGradedXInterfaceMaterialHypoElastic
   :members:

.. doxygenclass:: MarmotGradedEquilibratedInterfaceMaterialFiniteStrain
   :members:

.. doxygenclass:: MarmotGradedEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain
   :members:

.. doxygenclass:: Marmot::Materials::TwoSidedGradedLayer::LayerSolver
   :members:

.. doxygenstruct:: Marmot::Elements::GradientEnhancedTwoSidedStationFieldsFormulation
   :members:

.. doxygenstruct:: Marmot::Elements::GradientEnhancedStationFieldsFormulation
   :members:
