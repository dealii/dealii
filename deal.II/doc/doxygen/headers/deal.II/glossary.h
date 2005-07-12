/**
 * @page DEALGlossary Glossary
 *
 * <dl>
 *
 * <dt>@anchor GlossActive Active cells</dt>
 * <dd>Mesh cells not refined any further in the hierarchy.</dd>
 *
 * <dt>@anchor GlossGeneralizedSupport Generalized support points</dt>
 * <dd>While @ref GlossSupport "support points" allow very simple interpolation
 * into the finite element space, their concept is restricted to
 * @ref GlossLagrange "Lagrange elements". For other elements, more general
 * interpolation operators can be defined, often relying on integral values
 * or moments. Since these integral values are again computed using a
 * Quadrature rule, we consider them a generalization of support
 * points.
 *
 * Note that there is no simple relation between
 * @ref GlossShape "shape functions" and generalized support points as for
 * regular @ref GlossSupport "support points". Instead, FiniteElement defines
 * a couple of interpolation functions doing the actual interpolation.
 *
 * @note If a finite element is Lagrangian, generalized support points
 * and support points coincide.
 *
 * </dd>
 *
 * <dt>@anchor GlossLagrange Lagrange elements</dt>
 * <dd>Finite elements based on Lagrangian interpolation at
 * @ref GlossSupport "support points".</dd>
 *
 * <dt>@anchor GlossReferenceCell Reference cell</dt>
 * <dd>The hypercube [0,1]<sup>dim</sup>, on which all parametric finite
 * element shape functions are defined.</dd>
 *
 * <dt>@anchor GlossSupport Support points</dt>
 * <dd>Support points are by definition those points <i>p<sub>i</sub></i>,
 * such that for the shape functions <i>v<sub>j</sub></i> holds
 * <i>v<sub>j</sub>(p<sub>i</sub>) = &delta;<sub>ij</sub></i>. Therefore, a
 * finite element interpolation can be defined uniquely by the values in the
 * support points.
 *
 * Lagrangian elements fill the vector accessed by
 * FiniteElementBase::get_unit_support_points(), such that the
 * function FiniteElementBase::has_support_points() returns
 * <tt>true</tt>. Naturally, these support points are on the
 * @ref GlossReferenceCell "reference cell".  Then, FEValues can be used
 * (in conjuncton with a Mapping) to access support points on the
 * actual grid cells.
 *
 * @note The concept of @ref GlossSupport "support points" is
 * restricted to the finite element families based on Lagrange
 * interpolation. For a more general concept, see
 * @ref GlossGeneralizedSupport "generalized support points".
 * </dd>
 *
 * <dt>@anchor GlossUnitCell Unit cell</dt>
 * <dd>See @ref GlossReferenceCell "Reference cell".</dd>
 *
 * <dt>@anchor GlossUnitSupport Unit support points</dt>
 * <dd>@ref GlossSupport "Support points" on the reference cell, defined in
 * FiniteElementBase.
 * </dd>
 *
 * </dl>
 */
