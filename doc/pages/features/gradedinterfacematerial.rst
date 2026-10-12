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

The terms of :math:`\mathbf{u}^1` with the tangential gradient, which are of higher order in :math:`h`, are
neglected, so that :math:`\mathbf{w}` is a function of the macroscopic strain at the quadrature point. It is
condensed at every stress update, and the element has the displacements as its only unknown.

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

Limitations
-----------

* The terms of :math:`\mathbf{u}^1` with the tangential gradient are neglected, which is of higher order in the
  thickness.
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

.. doxygenstruct:: Marmot::Elements::GradientEnhancedStationFieldsFormulation
   :members:
