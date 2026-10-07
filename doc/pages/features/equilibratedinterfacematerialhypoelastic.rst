.. _equilibratedinterfacematerialhypoelastic:

Equilibrated interface material (hypoelastic bulk material)
===========================================================

Theory
------

``MarmotEquilibratedXInterfaceMaterialHypoElastic`` turns **any registered hypoelastic bulk material** into the
constitutive law of a thin interface layer of thickness :math:`h`, like the
:ref:`interface material <interfacematerialhypoelastic>`. Instead of one material that sees the average of the two
interface sides, the layer is made of **two faces of thickness** :math:`h/2`, each with its own bulk material,
its own state variables and its own surface gradient :math:`\mathbf{A}^{+}` and :math:`\mathbf{A}^{-}`.

Normal-gradient jump
^^^^^^^^^^^^^^^^^^^^

With the displacement jump :math:`\mathbf{w} = [\![\mathbf{u}]\!]`, the average surface gradient
:math:`\bar{\mathbf{A}} = \tfrac{1}{2}(\mathbf{A}^{+} + \mathbf{A}^{-})` and the connector
:math:`\mathbf{d}` between the two faces, split into its normal part :math:`\ell` and its tangential part
:math:`\mathbf{d}_\tau`, the shared average of the normal gradient is

.. math::

   \bar{\mathbf{g}} = \frac{1}{\ell} \left( \mathbf{w} - \bar{\mathbf{A}}\, \mathbf{d}_\tau \right).

The two faces do **not** get the same normal gradient. The jump

.. math::

   \mathbf{z} = \mathbf{g}^{+} - \mathbf{g}^{-}, \qquad
   \mathbf{g}^{\pm} = \bar{\mathbf{g}} \pm \tfrac{1}{2}\mathbf{z},

is an internal variable. The deformation gradient increment of a face is
:math:`\mathbf{A}^{\pm} + \mathbf{g}^{\pm} \otimes \mathbf{n}`, and the strain increment passed to its bulk
material is the symmetric part of it. For a layer without a connector, i.e. coincident faces, the material uses
:math:`\ell = h` and :math:`\mathbf{d}_\tau = \mathbf{0}`.

Traction equilibrium
^^^^^^^^^^^^^^^^^^^^

At every stress update :math:`\mathbf{z}` is determined by a local Newton iteration that enforces the continuity of
the traction,

.. math::

   \mathbf{r}(\mathbf{z}) = \sig^{+}(\mathbf{z})\, \mathbf{n} - \sig^{-}(\mathbf{z})\, \mathbf{n} = \mathbf{0}.

Every trial evaluation restarts from the committed state of the two bulk materials, so the state is never
accumulated over the local iterations. The Jacobian is the average one-sided acoustic tensor
:math:`\langle \mathbf{Q} \rangle = \tfrac{1}{2}(\mathbf{n}\,\mathbb{C}^{+}\mathbf{n} + \mathbf{n}\,\mathbb{C}^{-}\mathbf{n})`.
It stays invertible where the acoustic tensor of a single face degenerates, for instance at yield, so only the
average is ever inverted. Near perfect plasticity the nonlinear condition may have no solution at all, because both
faces sit on the yield surface and their tractions stop responding to :math:`\mathbf{z}`. The iteration then
floors at a small residual and the best iterate is accepted.

Gauge
^^^^^

Where both faces flow in the same direction, for instance under pure sliding of a near-perfectly plastic layer,
:math:`\langle \mathbf{Q} \rangle` loses its transverse stiffness and the component of :math:`\mathbf{g}` along
:math:`\mathbf{m} = \tfrac{1}{2}(\mathbf{N}^{+} + \mathbf{N}^{-})\,\mathbf{n}` is undetermined, with
:math:`\mathbf{N} = \operatorname{dev}\sig / \lVert \operatorname{dev}\sig \rVert`. The local solve is then
bistable. The material detects this case, when both faces update their state and the ratio of the smallest to the
largest singular value of :math:`\langle \mathbf{Q} \rangle` falls below :math:`10^{-6}`, and fixes the tangential
part of :math:`\mathbf{m}`: the local solve and the condensation of the tangent are restricted to the
two-dimensional subspace that the layer still resists.

The decision whether the gauge is active has to be the same in all iterations of an increment, because otherwise a
few points whose decision keeps changing can hold the global Newton iteration in a limit cycle. The state variables
are restored at the beginning of every iteration, so they cannot carry the decision. The caller therefore passes an
``IncrementMemory`` per quadrature point to ``computeStress``. The first decision taken in an increment is kept for
the whole increment, including its cutbacks. Without a memory, the decision is taken anew at every call.

Generalized stress and condensed tangent
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The external generalized strain is :math:`\mathbf{x} = (\mathbf{w}, \mathbf{A}^{+}, \mathbf{A}^{-})` and the
generalized stress :math:`\mathbf{p} = (\mathbf{f}, \mathbf{S}^{+}, \mathbf{S}^{-})`, where :math:`\mathbf{f}` is
conjugate to the jump and :math:`\mathbf{S}^{\pm}` to the surface gradients. The stress of a face is recovered as
:math:`\sig^{\pm} = (\mathbf{S}^{\pm} + \mathbf{f} \otimes \mathbf{d}_\tau)/(h/2)`.

The internal variable :math:`\mathbf{z}` is condensed statically. The resulting tangent is a full :math:`3 \times 3`
block tangent over :math:`\mathbf{x}`,

.. math::

   \frac{\partial \mathbf{p}}{\partial \mathbf{x}} =
   \begin{pmatrix}
     \mathbf{Q}_{ww}   & \mathbf{Q}_{wA^+}   & \mathbf{Q}_{wA^-}   \\
     \mathbf{Q}_{A^+w} & \mathbf{Q}_{A^+A^+} & \mathbf{Q}_{A^+A^-} \\
     \mathbf{Q}_{A^-w} & \mathbf{Q}_{A^-A^+} & \mathbf{Q}_{A^-A^-}
   \end{pmatrix},

in general without major symmetry. The cross blocks between :math:`\mathbf{A}^{+}` and :math:`\mathbf{A}^{-}` do not
vanish, because the equilibrium constraint couples the two faces.

Properties
----------

The material is created with the **name of the bulk material** and the property array

.. code-block:: text

   [ E, nu, h, <remaining properties of the bulk material> ]

as the :ref:`interface material <interfacematerialhypoelastic>`. Both faces receive the same properties. The
thickness :math:`h` has to be positive.

State variables
---------------

Per quadrature point the material stores the state variables of the bulk materials of the top and of the bottom
face (``topMaterialStateVars``, ``bottomMaterialStateVars``), the internal normal-gradient jump
(``normalGradientJump``) as a warm start for the next increment, and the committed stress of each face
(``stressPlus``, ``stressMinus``) in Voigt notation. The stress of a face is stored separately, because recovering
it from the single generalized force is exact only where the tractions are identical.

Implementation
--------------

.. doxygenclass:: MarmotEquilibratedXInterfaceMaterialHypoElastic
   :allow-dot-graphs:
