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
   * - ``IQUAD4UL``, ``IQUAD8UL``
     - 8, 16
     - As ``IQUAD4`` and ``IQUAD8`` at finite strain, with the finite-strain interface material
   * - ``ILINE2UL``
     - 4
     - Two-dimensional counterpart of ``IQUAD4UL`` between two linear lines
   * - ``YIQUAD4UL``, ``YIQUAD8UL``
     - 8, 16
     - As ``YIQUAD4`` and ``YIQUAD8`` at finite strain, with the :ref:`equilibrated finite-strain interface material
       <equilibratedinterfacematerialfinitestrain>`
   * - ``YILINE2UL``
     - 4
     - Two-dimensional counterpart of ``YIQUAD4UL``
   * - ``GYIQUAD4UL``, ``GYIQUAD8UL``
     - 8, 16
     - As ``YIQUAD4UL`` and ``YIQUAD8UL``, with a nonlocal field on the nodes of the bottom side
   * - ``GYILINE2UL``
     - 4
     - Two-dimensional counterpart of ``GYIQUAD4UL``
   * - ``GY2IQUAD4UL``, ``GY2IQUAD8UL``
     - 8, 16
     - As ``GYIQUAD4UL`` and ``GYIQUAD8UL``, with an independent nonlocal field on each face of the layer
   * - ``GY2ILINE2UL``
     - 4
     - Two-dimensional counterpart of ``GY2IQUAD4UL``
   * - ``YILINE2``
     - 4
     - Two-dimensional counterpart of ``YIQUAD4`` between two linear lines
   * - ``HIQUAD4``, ``HIQUAD8``
     - 8, 16
     - As ``IQUAD4`` and ``IQUAD8``, for a layer whose stiffness varies through the thickness, with the :ref:`graded
       interface material <gradedinterfacematerial>`
   * - ``HILINE2``
     - 4
     - Two-dimensional counterpart of ``HIQUAD4``
   * - ``HIQUAD4UL``, ``HIQUAD8UL``, ``HILINE2UL``
     - 8, 16, 4
     - As ``HIQUAD4``, ``HIQUAD8`` and ``HILINE2`` at finite strain
   * - ``GYS2IQUAD4UL``, ``GYS3IQUAD4UL``, ``GYS4IQUAD4UL``
     - 8
     - ``HIQUAD4UL`` with a gradient-enhanced bulk material and a nonlocal field at each of 2, 3 or 4 stations through
       the thickness. The same for ``IQUAD8`` (16 nodes) and ``ILINE2`` (4 nodes) as ``GYSnIQUAD8UL`` and
       ``GYSnILINE2UL``
   * - ``HYQUAD4``, ``HYQUAD8``, ``HYLINE2``
     - 8, 16, 4
     - As ``YIQUAD4``, ``YIQUAD8`` and ``YILINE2`` (equilibrated two-sided kinematics), for a layer whose stiffness
       varies through the thickness, with the two-sided :ref:`graded interface material <gradedinterfacematerial>`,
       which keeps the surface gradients of both faces and the tangential gradient of the correction of the
       microproblem
   * - ``HYQUAD4UL``, ``HYQUAD8UL``, ``HYLINE2UL``
     - 8, 16, 4
     - As ``HYQUAD4``, ``HYQUAD8`` and ``HYLINE2`` at finite strain
   * - ``GYH2IQUAD4UL``, ``GYH3IQUAD4UL``, ``GYH4IQUAD4UL``
     - 8
     - ``HYQUAD4UL`` with a gradient-enhanced bulk material and a nonlocal field at each of 2, 3 or 4 stations through
       the thickness. The same for ``IQUAD8`` (16 nodes) and ``ILINE2`` (4 nodes) as ``GYHnIQUAD8UL`` and
       ``GYHnILINE2UL``. For two stations and without the tangential gradient of the correction, the element is
       ``GY2IQUAD4UL``

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

The ``UL`` elements use ``FiniteStrainDisplacementJumpFormulation``. The interface material receives the *total*
displacements and surface gradients of the two sides and the connector between the paired faces, and evaluates a bulk
material of finite strain for the deformation gradient of the layer,

.. math::

   \mathbf{F} = \mathbf{I} + \langle\nabla_s \mathbf{u}\rangle
     + \frac{1}{\ell}\big( [\![\mathbf{u}]\!] - \langle\nabla_s \mathbf{u}\rangle\, \mathbf{d}_\tau \big)
     \otimes \mathbf{n}, \qquad \mathbf{X}^{+} - \mathbf{X}^{-} = \ell\, \mathbf{n} + \mathbf{d}_\tau .

The tangential part :math:`\mathbf{d}_\tau` of the connector keeps the layer free of stresses under a rigid rotation
of an inclined interface. The material properties are ``[h, properties of the bulk material...]``, where ``h`` is the
thickness of the layer. The force :math:`\mathbf{t} = (h/\ell)\, \mathbf{P}\mathbf{n}` and the surface stress
:math:`\mathbf{s} = h\, \mathbf{P}\mathbf{T} - (h/\ell)\, (\mathbf{P}\mathbf{n}) \otimes \mathbf{d}_\tau` follow
from the first Piola-Kirchhoff stress :math:`\mathbf{P}` of the bulk material and
:math:`\mathbf{T} = \mathbf{I} - \mathbf{n}\otimes\mathbf{n}`.

The ``GY`` elements use ``GradientEnhancedEquilibratedTwoSidedFormulation``. In addition to the displacements, they
carry a nonlocal field :math:`\bar{N}` with one degree of freedom on each node of the bottom side, which is the
midsurface of the layer, after the displacements in the vectors of the element. The bulk material is a
gradient-enhanced finite-strain material, such as the
:doc:`gradient-enhanced neo-Hookean damage material <gradientenhancedcompressibleneohookedamage>`. For a nonlocal
radius :math:`l` above the thickness :math:`h` of the layer, the field is flat across it, and the nonlocal equation
reduces to the midsurface,

.. math::

   h \left[ \mathbf{N}^\mathsf{T} \big( \bar{N} - \bar{L} \big) + l^2\, \nabla\mathbf{N}^\mathsf{T} \nabla\bar{N} \right] = 0,

where :math:`\bar{L}` is the average of the local driving fields of the two faces. The residual and the tangent of the
displacements, and of the nonlocal field, are coupled through the condensed derivatives of the
:ref:`equilibrated finite-strain interface material <equilibratedinterfacematerialfinitestrain>`.

The ``GY2`` elements use ``GradientEnhancedFaceFieldsFormulation``, in which the two faces of the layer have independent
nonlocal fields :math:`\bar N^{+}` and :math:`\bar N^{-}`, and so independent damage. Every node carries a nonlocal degree
of freedom: the top nodes the field of the top face and the bottom nodes the one of the bottom face. Each face has its own
nonlocal equation on the midsurface, with half of the thickness of the layer,

.. math::

   \frac{h}{2} \left[ \mathbf{N}^\mathsf{T} \big( \bar N^{\pm} - L^{\pm} \big)
   + l^2\, \nabla\mathbf{N}^\mathsf{T} \nabla\bar N^{\pm} \right] = 0 ,

where :math:`L^{\pm}` is the local driving field of the face, its plastic strain. The coupling through the thickness :math:`h`
of the layer comes from the derivative of the field through the thickness, for a field that varies linearly between the faces,

.. math::

   \pm \frac{(l^{\pm})^2}{h}\, \mathbf{N}^\mathsf{T} \big( \bar N^{+} - \bar N^{-} \big),

with the sign :math:`+` on the top and :math:`-` on the bottom face. Each face has the coefficient of its own nonlocal radius,
so that the radii of the faces may differ. The coupling vanishes for a small radius, where the faces are independent, and forces
:math:`\bar N^{+} = \bar N^{-}` for a radius that is large against the thickness, which is the element ``GYIQUAD4UL``. The
faces also interact through the traction equilibrium of the normal-gradient jump, which gives each field a derivative path to
the stress and to the local field of the other face.

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

.. doxygenstruct:: Marmot::Elements::FiniteStrainDisplacementJumpFormulation
   :members:

.. doxygenclass:: MarmotInterfaceMaterialFiniteStrain
   :members:

.. doxygenstruct:: Marmot::Elements::FiniteStrainEquilibratedTwoSidedFormulation
   :members:

.. doxygenstruct:: Marmot::Elements::GradientEnhancedEquilibratedTwoSidedFormulation
   :members:

.. doxygenstruct:: Marmot::Elements::GradientEnhancedFaceFieldsFormulation
   :members:
