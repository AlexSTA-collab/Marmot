.. _gradedinterfacematerial:

Graded Interface Materials
==========================

Theory
------

The interface materials of this page describe a layer of finite thickness :math:`h`, whose stiffness
:math:`\mathbb{C}(y)` varies through the thickness, with the coordinate :math:`y \in [0, 1]` from the second side
of the interface to the first one along the normal :math:`\mathbf{N}`. The stiffness may be given by layers of
constant properties, or by properties that vary continuously, and any bulk material of the library may be used.

The kinematics are those of the other interface materials: the displacement jump :math:`[\![\mathbf{u}]\!]` and the
average surface gradient :math:`\langle \nabla_s \mathbf{u} \rangle` give the strain of the layer
:math:`\boldsymbol{\varepsilon}^0` (small strains), or the deformation gradient
:math:`\mathbf{F}^0 = \mathbf{I} + \langle \nabla_s \mathbf{u} \rangle + \bar{\mathbf{g}} \otimes \mathbf{N}` with the
reduced normal derivative :math:`\bar{\mathbf{g}}` read off the connector between the faces (finite strains). Because
the stiffness varies, the strain does not: the first-order displacement :math:`\mathbf{u}^1` of the microproblem
corrects it,

.. math::

   \boldsymbol{\varepsilon}(y) = \boldsymbol{\varepsilon}^0 + \mathrm{sym}\big( \mathbf{w}(y) \otimes \mathbf{n} \big),
   \qquad
   \mathbf{F}(y) = \mathbf{F}^0 + \mathbf{w}(y) \otimes \mathbf{N}, \qquad \mathbf{w} = \partial_y \mathbf{u}^1 .

The variational form with a vanishing macroscopic test function gives the microproblem: the traction
:math:`\mathbf{t} = \boldsymbol{\sigma}(y)\,\mathbf{n}` (small strains) or :math:`\mathbf{P}(y)\,\mathbf{N}` (finite
strains) is the same at every height, and :math:`\mathbf{u}^1` vanishes at both sides,

.. math::

   \int_0^1 \mathbf{w} \,\mathrm{d}y = \mathbf{0} .

The materials with one averaged surface gradient (the elements ``HI*`` and ``GYS*``) neglect the terms of
:math:`\mathbf{u}^1` with the tangential gradient, so that :math:`\mathbf{w}` is a function of the macroscopic strain
at the quadrature point. The two-sided materials (the elements ``HY*`` and ``GYH*``) keep them, see
:ref:`below <gradedinterfacematerial-twosided>`. In both cases :math:`\mathbf{u}^1` is condensed at every stress
update, and the element has the displacements as its only unknown.

With the weights :math:`\omega_q` and the acoustic tensors :math:`\mathbf{Q}_q = \mathbf{n}\cdot\mathbb{C}_q\cdot\mathbf{n}`
(finite strains: :math:`\mathbb{A}_q[\mathbf{N}, \mathbf{N}]`) of the integration points through the thickness, a
Newton step for the corrections :math:`\mathbf{w}_q` and the common traction :math:`\mathbf{t}` is

.. math::

   \delta \mathbf{t} = \Big( \sum_q \omega_q \mathbf{Q}_q^{-1} \Big)^{-1}
     \Big( \sum_q \omega_q \mathbf{Q}_q^{-1} ( \mathbf{a}_q - \mathbf{t} ) - \sum_q \omega_q \mathbf{w}_q \Big),
   \qquad
   \delta \mathbf{w}_q = \mathbf{Q}_q^{-1} \big( \delta \mathbf{t} - ( \mathbf{a}_q - \mathbf{t} ) \big),

where :math:`\mathbf{a}_q` is the traction of the point. A linear material converges in one step. The step is damped by
a line search on the norm of the residual, because a yielding point has a small tangent along its flow direction. The
tangent of the layer is condensed over the corrections,

.. math::

   \bar{\mathbb{C}} = \sum_q \omega_q \mathbb{C}_q \Big( \mathbb{I} + \frac{\partial\, \mathrm{sym}(\mathbf{w}_q \otimes \mathbf{n})}
   {\partial \boldsymbol{\varepsilon}^0} \Big),

and the force and the surface stress are those of the plain interface material for the through-thickness mean of the
stress. For a single homogeneous layer the response is that of the plain interface material.

.. _gradedinterfacematerial-twosided:

Two-sided materials and the tangential gradient of the correction
-----------------------------------------------------------------

The one-sided materials give the layer one surface gradient, the average of the two faces, so that they carry the
membrane terms of the layer, of the order of :math:`h`, but not the terms of the order of :math:`h^3` that come from the
variation of the surface gradient through the thickness, which the equilibrated two-sided materials (the elements
``YI*``, ``GYI*`` and ``GY2I*``) have. The two-sided graded materials generalize the equilibrated material, which
has the two faces as its only integration points through the thickness, to any number of integration points with
different bulk materials. They have the contract of the equilibrated material, and are used by the same element
formulation: the surface gradients :math:`\mathbf{A}^\pm` of both faces are given, and a surface stress is returned for
each of them. The gradient at the height :math:`y` is

.. math::

   \mathbf{G}(y) = \frac{1}{\ell}\, [\![\mathbf{u}]\!] \otimes \mathbf{n} + (1-y)\,\mathbf{A}^- + y\,\mathbf{A}^+
   + \mathbf{B}_\tau \bar{\mathbf{A}} + \mathbf{w}(y) \otimes \mathbf{n} + h\,\nabla_T \mathbf{u}^1(y),

with the normal separation :math:`\ell` of the faces and the correction :math:`\mathbf{B}_\tau \bar{\mathbf{A}}` of
their tangential separation, as for the equilibrated material. The generalized stress is
:math:`\mathbf{p} = h \sum_q \omega_q \mathbf{B}_q^\mathsf{T} \boldsymbol{\sigma}_q`, which for a homogeneous layer of
two points at the faces is that of the equilibrated material.

The last term of the gradient, which comes from the term :math:`h\,u^1_{k,r} T_{rl}` of the variational form, is of
the same order as the terms of the surface gradients, because :math:`\mathbf{u}^1` follows the strain of the layer and
the strain varies along the interface with the jump: :math:`\nabla_T [\![\mathbf{u}]\!] = \mathbf{A}^+ - \mathbf{A}^-`.
The tangential derivative of the correction of the microproblem is therefore closed with the derivative of the
correction with respect to the strain of the layer, which the condensation of the microproblem provides,

.. math::

   \partial_r \boldsymbol{\varepsilon} = \frac{1}{\ell}\, \mathrm{sym}\big( ( \mathbf{A}^+ - \mathbf{A}^- )\, \mathbf{e}_r
   \otimes \mathbf{n} \big), \qquad
   \partial_r \mathbf{w}_q = \frac{\partial \mathbf{w}_q}{\partial \boldsymbol{\varepsilon}}\, \partial_r \boldsymbol{\varepsilon},

and the second derivatives of the displacement of the faces are neglected. Two things follow, both proportional to
:math:`\mathbf{A}^+ - \mathbf{A}^-`:

* the strain of the point gains the tangential gradient :math:`h\,\nabla_T \mathbf{u}^1(y_q)` of
  :math:`\mathbf{u}^1(y) = \int_0^y \mathbf{w}\,\mathrm{d}y`, and
* the traction is no longer constant through the thickness, because the test function with the tangential gradient
  brings the tangential divergence of the stress into the microproblem,
  :math:`\partial_y ( \boldsymbol{\sigma}\,\mathbf{n} ) = - h\, \partial_r ( \boldsymbol{\sigma}\, \mathbf{e}_r )`.

The derivative of the correction is the one of the converged state of the increment, which is iterated until it does
not change. The tangent neglects the dependence of that derivative on the state, which is of the order of the second
derivative of the stress, so that the tangent is exact without the closure and to the order of the strain times the
closure terms with it. A layered strip with a jump varying linearly along the interface has a polynomial plane-strain
solution, in which the shear stress varies linearly through every layer; the material reproduces it, also the
moments :math:`h\int y\,\boldsymbol{\sigma}\,\mathrm{d}y` of the stress, which are the surface stresses of the
faces, to the precision of the arithmetic. The closure is switched off by the optional last material property, which
makes the material the equilibrated one for two points at the faces.

.. note::

   The equilibrated two-point material integrates the surface gradient, which is linear through the thickness, with
   the trapezoidal rule. The contribution of :math:`\mathbf{A}^+ - \mathbf{A}^-` to the moments
   :math:`h\int y\,\boldsymbol{\sigma}\,\mathrm{d}y` of a homogeneous layer is then :math:`h/4` times the stiffness,
   instead of the :math:`h/12` of the exact integral. The two-sided graded materials with Gauss points (modes 0 and
   1) have the exact one, while the points at the stations (mode 2) have the trapezoidal rule, as the faces have.

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

One more property may follow the groups, the flag of the tangential gradient of the correction: the two-sided
materials keep it when the flag is absent or not zero, and switch it off when it is zero. The one-sided materials
ignore it.

The materials
-------------

``MarmotGradedInterfaceMaterialHypoElastic``
   Small-strain material for the elements ``HILINE2``, ``HIQUAD4`` and ``HIQUAD8``. The bulk material is a
   hypoelastic material.

``MarmotGradedInterfaceMaterialFiniteStrain``
   Finite-strain material for the elements ``HILINE2UL``, ``HIQUAD4UL`` and ``HIQUAD8UL``. The bulk material is a
   finite-strain material. For a homogeneous layer the response is that of ``MarmotInterfaceMaterialFiniteStrain``.

``MarmotGradedGradientEnhancedInterfaceMaterialFiniteStrain``
   Gradient-enhanced finite-strain material for the elements ``GYSnILINE2UL``, ``GYSnIQUAD4UL`` and
   ``GYSnIQUAD8UL`` with :math:`n = 2, 3, 4` stations. Its bulk material is gradient-enhanced (for example
   ``DUFOURMODEL``), the profile is in mode 2, and every station has a nonlocal field of its own, so that the damage
   varies through the thickness. The nonlocal equation of the station :math:`s` is, reduced to the midsurface,

   .. math::

      h\,\omega_s \Big[ \mathbf{N}^\mathsf{T} ( \bar N_s - L_s ) + l_s^2\, \nabla\mathbf{N}^\mathsf{T} \nabla \bar N_s \Big]
      + \frac{l_s^2}{h} \Big[ \frac{\mathbf{N}^\mathsf{T}( \bar N_s - \bar N_{s+1} )}{y_{s+1} - y_s}
      + \frac{\mathbf{N}^\mathsf{T}( \bar N_s - \bar N_{s-1} )}{y_s - y_{s-1}} \Big] = 0,

   where :math:`L_s` is the local driving field of the station, :math:`l_s` its nonlocal radius, and the last term the
   derivative of the field through the thickness. For two stations it is the element with a nonlocal field for each
   face (``GY2``). Its nonlocal residuals are the same, and the virtual work of every virtual displacement with a
   common surface gradient of the faces is the same; the displacement residuals differ where the faces carry different
   surface stresses, because the layer has one tangential gradient. The fields of the bottom station and of the
   inner stations are degrees of freedom of the nodes of the bottom side, the field of the top station of the nodes of
   the top side.

``MarmotGradedXInterfaceMaterialHypoElastic``
   Two-sided small-strain material for the elements ``HYLINE2``, ``HYQUAD4`` and ``HYQUAD8``, with the contract of
   ``MarmotEquilibratedXInterfaceMaterialHypoElastic``.

``MarmotGradedEquilibratedInterfaceMaterialFiniteStrain``
   Two-sided finite-strain material for the elements ``HYLINE2UL``, ``HYQUAD4UL`` and ``HYQUAD8UL``, with the contract
   of ``MarmotEquilibratedInterfaceMaterialFiniteStrain``, which receives the total kinematics.

``MarmotGradedEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain``
   Two-sided gradient-enhanced finite-strain material for the elements ``GYHnILINE2UL``, ``GYHnIQUAD4UL`` and
   ``GYHnIQUAD8UL`` with :math:`n = 2, 3, 4` stations. The nonlocal fields and their equations are those of
   ``GYSn``, and the momentum equation has the surface stress of each face. For two stations and without the
   closure, the element is ``GY2`` to the residual of the displacements and the full tangent, for the damage model
   and for the Dufour model, which is the way to test the Dufour model through a graded layer exactly.

The two-sided materials solve the microproblem and the condensation with the same solver
(``Marmot::Materials::TwoSidedGradedLayer``) for small and finite strains.

Limitations
-----------

* The one-sided materials neglect the terms of :math:`\mathbf{u}^1` with the tangential gradient. The two-sided
  materials keep them in the closure above, which neglects the second derivatives of the displacement of the faces,
  and the dependence of the derivative of the correction on the state, which is not in the tangent. With a damage
  model that dependence includes the second derivative of the damage with respect to the nonlocal field, and the
  tangent with respect to the nonlocal fields has the corresponding error when the faces have very different
  surface gradients.
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

.. doxygenclass:: MarmotGradedInterfaceMaterialHypoElastic
   :members:

.. doxygenclass:: MarmotGradedInterfaceMaterialFiniteStrain
   :members:

.. doxygenclass:: MarmotGradedGradientEnhancedInterfaceMaterialFiniteStrain
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
