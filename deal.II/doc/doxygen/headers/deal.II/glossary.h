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
 * If a finite element is Lagrangian, generalized support points
 * and support points coincide.
 * </dd>
 *
 * <dt>@anchor GlossInterpolation Interpolation with finite elements</dt>
 * <dd>The purpose of interpolation with finite elements is computing
 * a vector of coefficients representing a finite element function,
 * such that the @ref GlossNodes "node values" of the original
 * function and the finite element function coincide. Therefore, the
 * interpolation process consists of evaluating all @ref GlossNodes
 * "node functionals" <i>N<sub>i</sub></i> for the given function
 * <i>f</i> and store the result as entry <i>i</i> in the coefficient
 * vector.
 *
 * <dt>@anchor GlossLagrange Lagrange elements</dt>
 * <dd>Finite elements based on Lagrangian interpolation at
 * @ref GlossSupport "support points".</dd>
 *
 * <dt>@anchor GlossNodes Node values or node functionals</dt>
 *
 * <dd>It is customary to define a FiniteElement as a pair consisting
 * of a local function space and a set of node values $N_i$ on the
 * mesh cells (usually defined on the @ref GlossReferenceCell
 * "reference cell"). Then, the basis of the local function space is
 * chosen such that $N_i(v_j) = \delta_{ij}$, the Kronecker delta.
 *
 * This splitting has several advantages, concerning analysis as well
 * as implementation. For the analysis, it means that conformity with
 * certain spaces (FiniteElementData::Conformity), e.g. continuity, is
 * up to the node values. In deal.II, it helps simplifying the
 * implementation of more complex elements like FE_RaviartThomas
 * considerably.
 *
 * Examples for node functionals are values in @ref GlossSupport
 * "support points" and moments with respect to Legendre
 * polynomials. Let us give some examples:
 *
 * <table><tr>
 * <th>Element</th>
 * <th>Function space</th>
 * <th>Node values</th></tr>
 * <tr><th>FE_Q, FE_DGQ</th>
 * <td><i>Q<sub>k</sub></i></td>
 * <td>values in support points</td></tr>
 * <tr><th>FE_DGP</th>
 * <td><i>P<sub>k</sub></i></td>
 * <td>moments with respect to Legendre polynomials</td></tr>
 * <tr><th>FE_RaviartThomas (2d)</th>
 * <td><i>Q<sub>k+1,k</sub> x Q<sub>k,k+1</sub></i></td>
 * <td>moments on edges and in the interior</td></tr>
 * <tr><th>FE_RaviartThomasNodal</th>
 * <td><i>Q<sub>k+1,k</sub> x Q<sub>k,k+1</sub></i></td>
 * <td>Gauss points on edges(faces) and anisotropic Gauss points in the interior</td></tr>
 * </table>
 *
 * <dt>@anchor GlossReferenceCell Reference cell</dt>
 * <dd>The hypercube [0,1]<sup>dim</sup>, on which all parametric finite
 * element shape functions are defined.</dd>
 *
 * <dt>@anchor GlossShape Shape functions</dt> <dd>The restriction of
 * the finite element basis functions to a single grid cell.</dd>
 *
 * <dt>@anchor GlossSupport Support points</dt> <dd>Support points are
 * by definition those points $p_i$, such that for the shape functions
 * $v_j$ holds $v_j(p_i) = \delta_{ij}$. Therefore, a finite element
 * interpolation can be defined uniquely by the values in the support
 * points.
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
