.. _equilibratedinterfacematerialfinitestrain:

Equilibrated interface material (finite-strain bulk material)
=============================================================

Theory
------

``MarmotEquilibratedInterfaceMaterialFiniteStrain`` is the finite-strain counterpart of the
:ref:`equilibrated interface material <equilibratedinterfacematerialhypoelastic>`. The layer of thickness :math:`h`
consists of two faces of thickness :math:`h/2`, each with its own bulk material of finite strain, its own state
variables and its own surface gradient :math:`\mathbf{A}^{\pm}`. The material receives the **total** displacements and
surface gradients.

With the displacement jump :math:`\mathbf{w} = [\![\mathbf{u}]\!]`, the average surface gradient
:math:`\bar{\mathbf{A}}` and the connector :math:`\mathbf{d} = \ell\,\mathbf{N} + \mathbf{d}_\tau` between the faces,
the faces carry the normal gradients :math:`\mathbf{g}^{\pm} = \bar{\mathbf{g}} \pm \mathbf{z}/2`, with

.. math::

   \bar{\mathbf{g}} = \frac{1}{\ell}\left( \mathbf{w} - \bar{\mathbf{A}}\,\mathbf{d}_\tau \right), \qquad
   \mathbf{F}^{\pm} = \mathbf{I} + \mathbf{A}^{\pm} + \mathbf{g}^{\pm} \otimes \mathbf{N}.

The deformation gradients are affine in :math:`(\mathbf{w}, \mathbf{A}^{+}, \mathbf{A}^{-}, \mathbf{z})`. The internal
variable :math:`\mathbf{z}` is condensed at every stress update by a local Newton iteration, which enforces the
continuity of the nominal traction on the reference normal,

.. math::

   \mathbf{P}^{+}\mathbf{N} = \mathbf{P}^{-}\mathbf{N}.

Its Jacobian is the average acoustic tensor :math:`\tfrac{1}{2}\big(\mathbf{N}\,\mathbb{A}^{+}\mathbf{N} +
\mathbf{N}\,\mathbb{A}^{-}\mathbf{N}\big)`, with :math:`\mathbb{A} = \partial\mathbf{P}/\partial\mathbf{F}`. The
condensed tangent is the full :math:`3\times3` block tangent over :math:`(\mathbf{w}, \mathbf{A}^{+},
\mathbf{A}^{-})`, including the cross blocks.

The material does not have the gauge of the small-strain material, which fixes the undetermined component of
:math:`\mathbf{z}` for a perfectly plastic layer. For near-perfect plasticity the local residual stalls, and the best
iterate is accepted.

Properties
----------

The material properties are ``[p0, p1, h, ...]``: the first two properties of the bulk material, the thickness of
the layer normal to the interface, and the remaining properties of the bulk material.

The bulk material may also be a gradient-enhanced finite-strain material, such as the
:doc:`gradient-enhanced neo-Hookean damage material <gradientenhancedcompressibleneohookedamage>`. It is then driven
with its nonlocal field held at zero by ``MarmotMaterialFiniteStrainFromGradientEnhanced``. This is exact as long as
the stress does not depend on the field at that value, which the adapter checks on every call, and it refuses the
material otherwise. The nonlocal regularization itself is lost, because the interface carries no nonlocal field.

Implementation
--------------

.. doxygenclass:: MarmotEquilibratedInterfaceMaterialFiniteStrain
   :members:

.. doxygenclass:: MarmotMaterialFiniteStrainFromGradientEnhanced
   :members:
