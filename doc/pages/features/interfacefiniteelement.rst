.. _interfacefiniteelement:

Interface Finite Element
========================

Theory
------

The interface element models a **zero-thickness interface** between two bodies, or between two parts of one
body, by means of the displacement jump across it. The element has :math:`2 n` nodes: the first :math:`n` nodes
form the bottom side and the last :math:`n` nodes the top side of the interface, so that the element is
described by one surface of :math:`n` nodes in the reference configuration. The kinematics are
small-deformation (linearized).

At each quadrature point of the interface surface, the element evaluates the shape functions
:math:`\mathbf{N}` of one side, the surface metric, the unit normal :math:`\mathbf{n}`, and the surface
gradient operator :math:`\mathbf{B}_s` of one side. The two sides are combined into the operators

.. math::

   [\![\mathbf{u}]\!] = \mathbf{u}^{+} - \mathbf{u}^{-} = \mathbf{N}_{\mathrm{jump}}\, \mathbf{q}, \qquad
   \mathbf{N}_{\mathrm{jump}} = \big[\, -\mathbf{N} \;\; +\mathbf{N} \,\big],

.. math::

   \bar{\nabla}_s \mathbf{u} = \tfrac{1}{2}\big( \nabla_s \mathbf{u}^{-} + \nabla_s \mathbf{u}^{+} \big)
   = \mathbf{B}_{\mathrm{avg}}\, \mathbf{q}, \qquad
   \mathbf{B}_{\mathrm{avg}} = \tfrac{1}{2}\big[\, \mathbf{B}_s \;\; \mathbf{B}_s \,\big],

where :math:`\mathbf{q}` is the vector of nodal displacements.

The constitutive response is provided by the :ref:`interface material <interfacematerialhypoelastic>`.
Per quadrature point, the material receives the increments of the jump and of the surface gradient, and it returns
the force :math:`\mathbf{t}` (conjugate to the jump), the surface stress :math:`\mathbf{s}` (conjugate to
the average surface gradient), and the tangent operators :math:`\mathbf{Q}`, :math:`\mathbb{Z}`,
:math:`\mathbf{H}` and :math:`\mathbb{Y}`.

The element internal force vector and the tangent stiffness are

.. math::

   \mathbf{P}_e = \sum_{qp} \Big( \mathbf{N}_{\mathrm{jump}}^\mathsf{T}\, \mathbf{t}
                   + \mathbf{B}_{\mathrm{avg}}^\mathsf{T}\, \mathbf{s} \Big)\, J_0\, w,

.. math::

   \mathbf{K}_e = \sum_{qp} \Big( \mathbf{N}_{\mathrm{jump}}^\mathsf{T} \mathbf{Q}\, \mathbf{N}_{\mathrm{jump}}
     + \mathbf{B}_{\mathrm{avg}}^\mathsf{T} (\mathbb{Z} + \mathbb{Y})\, \mathbf{B}_{\mathrm{avg}}
     + \mathbf{N}_{\mathrm{jump}}^\mathsf{T} \mathbf{H}\, \mathbf{B}_{\mathrm{avg}}
     + \mathbf{B}_{\mathrm{avg}}^\mathsf{T} \mathbf{H}^\mathsf{T} \mathbf{N}_{\mathrm{jump}} \Big)\, J_0\, w,

where :math:`J_0 = \sqrt{\det \mathbf{G}}` is the surface Jacobian determinant (:math:`\mathbf{G}` is the
surface metric), :math:`w` the quadrature weight, and :math:`\sum_{qp}` the sum over quadrature points.
The sign convention is :math:`\mathbf{K}_e = +\partial \mathbf{P}_e / \partial \mathbf{q}`. The element has
no inertia terms.

Elements
--------

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Name
     - Nodes
     - Description
   * - ``IQUAD4``
     - 8
     - Three-dimensional interface element between two bilinear quadrilaterals (4 bottom and 4 top nodes),
       full integration
   * - ``IQUAD8``
     - 16
     - Three-dimensional interface element between two serendipity quadrilaterals (8 bottom and 8 top nodes),
       full integration
   * - ``YIQUAD4``
     - 8
     - As ``IQUAD4``, with the :ref:`equilibrated two-sided interface material <equilibratedinterfacematerialhypoelastic>`
   * - ``YIQUAD8``
     - 16
     - As ``IQUAD8``, with the equilibrated two-sided interface material
   * - ``YILINE2``
     - 4
     - Two-dimensional counterpart of ``YIQUAD4`` between two linear lines

The element is a template in the spatial dimension and the number of nodes. Besides the registered
three-dimensional ``IQUAD4``, a two-dimensional instantiation embeds its quadrature-point state in the
three-dimensional material layout.

The element is built from two parts. ``InterfaceFiniteElementBase`` holds what all interface elements have in
common: the geometry of the interface and of its quadrature points, the nodes, degrees of freedom and state
variables, the assignment of the material and the kinematic increments at the quadrature points. A *formulation*
supplies what differs between the elements: the interface material, the state variables stored at a quadrature
point, and the evaluation of the material response together with the assembly of the residual and the tangent.
``DisplacementJumpFormulation`` is the formulation described here, and ``InterfaceFiniteElement`` is the base with
this formulation.

The ``Y`` elements use ``EquilibratedTwoSidedFormulation``. Each side keeps its own surface gradient
:math:`\mathbf{A}^{\pm}`, and the material condenses the normal-gradient jump by enforcing traction equilibrium
:math:`\boldsymbol{\sigma}^{+}\mathbf{n} = \boldsymbol{\sigma}^{-}\mathbf{n}`. The faces may be separated: the
normal is oriented from the bottom to the top face, and nodes with a connector that has no component along the
normal are rejected when the element is initialized. The residual and the tangent are

.. math::

   \mathbf{P}_e = \sum_{qp} \Big( \mathbf{N}_{\mathrm{jump}}^\mathsf{T}\mathbf{f}
     + \mathbf{B}_{+}^\mathsf{T}\mathbf{S}^{+} + \mathbf{B}_{-}^\mathsf{T}\mathbf{S}^{-} \Big) J_0\, w,

with the tangent assembled from the nine blocks of the condensed material tangent over
:math:`(\mathbf{w}, \mathbf{A}^{+}, \mathbf{A}^{-})`.

Properties
----------

The interface material is defined by the section: its name is the name of the bulk material and its
properties are ``[E, nu, h, ...]``, see :ref:`interfacematerialhypoelastic`.
The first (optional) element property scales the integration measure :math:`J_0 w`; it defaults to 1.

State variables
---------------

Per quadrature point the element stores, in this order: the force :math:`\mathbf{t}`, the surface stress
:math:`\mathbf{s}`, the accumulated displacements :math:`[\mathbf{u}^{+}, \mathbf{u}^{-}]`, the accumulated
surface strains, followed by the state variables of the interface material (padded for alignment).

Implementation
--------------

.. doxygenclass:: Marmot::Elements::InterfaceFiniteElementBase
   :allow-dot-graphs:

.. doxygenstruct:: Marmot::Elements::DisplacementJumpFormulation
   :members:

.. doxygenstruct:: Marmot::Elements::EquilibratedTwoSidedFormulation
   :members:
