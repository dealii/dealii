/**
 * @page DEALGlossary Glossary
 *
 * This glossary explains a few terms that are frequently used in the
 * documentation of classes of deal.II. The glossary often only gives
 * a microscopic view of a particular concept; if you struggle with
 * the bigger picture, it may therefore also be worth to consult the
 * global overview of classes on the main/@ref index page.
 *
 * <dl>
 *
 * <dt>@anchor GlossActive <b>Active cells</b></dt>
 * <dd>Mesh cells not refined any further in the hierarchy.</dd>
 *
 *
 * <dt>@anchor GlossFaceOrientation <b>Face orientation</b></dt>
 * <dd>In a triangulation, the normal vector to a face
 * can be deduced from the face orientation by
 * applying the right hand side rule (x,y -> normal).  We note, that
 * in the standard orientation of faces in 2d, faces 0 and 2 have
 * normals that point into the cell, and faces 1 and 3 have normals
 * pointing outward. In 3d, faces 0, 2, and 4
 * have normals that point into the cell, while the normals of faces
 * 1, 3, and 5 point outward. This information, again, can be queried from
 * GeometryInfo<dim>::unit_normal_orientation.
 *
 * However, it turns out that a significant number of 3d meshes cannot
 * satisfy this convention. This is due to the fact that the face
 * convention for one cell already implies something for the
 * neighbor, since they share a common face and fixing it for the
 * first cell also fixes the normal vectors of the opposite faces of
 * both cells. It is easy to construct cases of loops of cells for
 * which this leads to cases where we cannot find orientations for
 * all faces that are consistent with this convention.
 *
 * For this reason, above convention is only what we call the
 * <em>standard orientation</em>. deal.II actually allows faces in 3d
 * to have either the standard direction, or its opposite, in which
 * case the lines that make up a cell would have reverted orders, and
 * the normal vector would have the opposite direction. You can ask a
 * cell whether a given face has standard orientation by calling
 * <tt>cell->face_orientation(face_no)</tt>: if the result is @p true,
 * then the face has standard orientation, otherwise its normal vector
 * is pointing the other direction. There are not very many places in
 * application programs where you need this information actually, but
 * a few places in the library make use of this. Note that in 2d, the
 * result is always @p true.
 *
 *
 * <dt>@anchor GlossGeneralizedSupport <b>Generalized support points</b></dt>
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
 *
 * <dt>@anchor GlossInterpolation <b>Interpolation with finite elements</b></dt>
 * <dd>The purpose of interpolation with finite elements is computing
 * a vector of coefficients representing a finite element function,
 * such that the @ref GlossNodes "node values" of the original
 * function and the finite element function coincide. Therefore, the
 * interpolation process consists of evaluating all @ref GlossNodes
 * "node functionals" <i>N<sub>i</sub></i> for the given function
 * <i>f</i> and store the result as entry <i>i</i> in the coefficient
 * vector.
 *
 *
 * <dt>@anchor GlossLagrange <b>Lagrange elements</b></dt>
 * <dd>Finite elements based on Lagrangian interpolation at
 * @ref GlossSupport "support points".</dd>
 *
 *
 * <dt>@anchor GlossNodes <b>Node values or node functionals</b></dt>
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
 *
 * <dt>@anchor GlossReferenceCell <b>Reference cell</b></dt>
 * <dd>The hypercube [0,1]<sup>dim</sup>, on which all parametric finite
 * element shape functions are defined.</dd>
 *
 *
 * <dt>@anchor GlossShape <b>Shape functions</b></dt> <dd>The restriction of
 * the finite element basis functions to a single grid cell.</dd>
 *
 *
 * <dt>@anchor GlossSupport <b>Support points</b></dt> <dd>Support points are
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
 *
 * <dt>@anchor GlossTargetComponent <b>Target component</b></dt> <dd>When
 * vectors and matrices are grouped into blocks by component, it is
 * often desirable to collect several of the original components into
 * a single one. This could be for instance, grouping the velocities
 * of a Stokes system into a single block.</dd>
 *
 *
 * <dt>@anchor GlossUnitCell <b>Unit cell</b></dt>
 * <dd>See @ref GlossReferenceCell "Reference cell".</dd>
 *
 *
 * <dt>@anchor GlossUnitSupport <b>Unit support points</b></dt>
 * <dd>@ref GlossSupport "Support points" on the reference cell, defined in
 * FiniteElementBase. For example, the usual Q1 element in 1d has support
 * points  at <tt>x=0</tt> and <tt>x=1</tt> (and similarly, in higher
 * dimensions at the vertices of the unit square or cube). On the other
 * hand, higher order Lagrangian elements have unit support points also
 * in the interior of the unit line, square, or cube.
 * </dd>
 *
 * </dl>
 */
