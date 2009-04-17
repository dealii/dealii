//-------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005, 2006, 2007, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-------------------------------------------------------------------------

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
 * <dt class="glossary">@anchor GlossActive Active cells</dt>
 * <dd>Mesh cells not refined any further in the hierarchy.</dd>
 *
 * <dt class="glossary">@anchor GlossBlock block</dt>
 * <dd>Blocks were introduced in BlockVector,
 * BlockSparseMatrix and related classes. These are used to reflect the
 * structure of a PDE system in linear algebra, in particular allowing
 * for modular solvers for problems with multiple solution components.
 * How to implement this is described in more detail in the
 * @ref vector_valued report and the tutorial programs referenced
 * therein.
 *
 * Originally, this concept was intermixed with the idea of the vector
 * @ref GlossComponent "component". Since the introduction of
 * non-@ref GlossPrimitive "primitive" elements, they became different. Take
 * for instance the solution of the mixed Laplacian system with
 * FE_RaviartThomas. There, the first <tt>dim</tt> components are the
 * directional derivatives. Since the shape functions are linear
 * combinations of those, they constitute only a single block. The
 * primal function <i>u</i> would be in the second block, but in the
 * <tt>dim+1</tt>st component.
 *
 * In most cases, when you subdivide a matrix or vector into blocks, you do so
 * by creating one block for each vector component. However, this is not
 * always so, and the DoFRenumbering::component_wise function allows to group
 * several vector components into the same block (see, for example, the @ref
 * step_22 "step-22" or @ref step_31 "step-31" tutorial programs, as opposed
 * to @ref step_20 "step-20").  </dd>
 *
 * <dt class="glossary">@anchor GlossComponent component</dt>
 *
 * <dd>For vector functions, component denotes the index in the
 * vector. For instance, in the mixed Laplacian system, the first
 * <tt>dim</tt> components are the derivatives in each coordinate
 * direction and the last component is the primal function <i>u</i>.
 *
 * Originally, components were not distinguished from @ref GlossBlock
 * "blocks", but since the introduction of non-@ref GlossPrimitive
 * "primitive" elements, they have to be distinguished. See
 * FiniteElementData::n_components() and the documentation of
 * FiniteElement</dd>
 *
 * <dt class="glossary">@anchor GlossFaceOrientation Face orientation</dt>
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
 * The only places in the library where face orientations play a
 * significant role are in the Triangulation and its accessors, and in
 * the QProjector class and its users.
 *
 *
 * <dt class="glossary">@anchor GlossGeneralizedSupport Generalized support points</dt>
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
 * <dt class="glossary">@anchor hp_paper hp paper</dt>
 * <dd>The "hp paper" is a paper by W. Bangerth and O. Kayser-Herold, titled
 * "Data Structures and Requirements for hp Finite Element Software", that
 * describes many of the algorithms and data structures used in the implementation
 * of the hp framework of deal.II. In particular, it summarizes many of the
 * tricky points that have to be considered for hp finite elements using continuous
 * elements.
 *
 * The full reference for this paper is as follows:
 * @code
@Article{BK07,
  author = 	 {Wolfgang Bangerth and Oliver Kayser-Herold},
  title = 	 {Data Structures and Requirements for $hp$ Finite Element
                  Software},
  journal = 	 {ACM Trans. Math. Softw.},
  year = 	 2009,
  volume = 	 36,
  number = 	 1,
  pages = 	 {4/1--4/31}
} 
 * @endcode
 * It is available as Technical Report ISC-07-04-MATH from the
 * <a href="http://www.isc.tamu.edu/publications-reports/technical_reports">Institute
 * for Scientific Computation, Texas A&amp;M University</a>, and also
 * from http://www.math.tamu.edu/~bangerth/publications.html .
 *
 * The numerical examples shown in that paper are generated with a slightly
 * modified version of @ref step_27 "step-27". The main difference to that
 * tutorial program is that various operations in the program were timed for
 * the paper to compare different options and show that $hp$ methods are
 * really not all that expensive.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossInterpolation Interpolation with finite elements</dt>
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
 * <dt class="glossary">@anchor GlossLagrange Lagrange elements</dt>
 * <dd>Finite elements based on Lagrangian interpolation at
 * @ref GlossSupport "support points".</dd>
 *
 *
 * <dt class="glossary">@anchor GlossNodes Node values or node functionals</dt>
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
 * <dt class="glossary">@anchor GlossPrimitive Primitive finite
 * elements</dt> <dd>Finite element shape function sets with a unique
 * relation from shape function number to vector @ref GlossComponent
 * "component". What this means is that each shape function of a
 * vector-valued element has exactly one-nonzero component if an
 * element is primitive. This includes, in particular, all scalar
 * elements as well as vector-valued elements assembled via the
 * FESystem class from other primitive (for example scalar) elements
 * as shown in @ref step_8 "step-8", @ref step_29 "step_29" or @ref
 * step_22 "step-22". On the other hand, the FE_RaviartThomas class used
 * in @ref step_20 "step-20" and @ref step_21 "step-21", or the
 * FE_Nedelec class provide non-primitive finite elements because
 * there each vector-value shape function may have several non-zero
 * components.</dd>
 *
 * <dt class="glossary">@anchor GlossReferenceCell Reference cell</dt>
 * <dd>The hypercube [0,1]<sup>dim</sup>, on which all parametric finite
 * element shape functions are defined.</dd>
 *
 *
 * <dt class="glossary">@anchor GlossShape Shape functions</dt> <dd>The restriction of
 * the finite element basis functions to a single grid cell.</dd>
 *
 *
 * <dt class="glossary">@anchor GlossSupport Support points</dt> <dd>Support points are
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
 * (in conjunction with a Mapping) to access support points on the
 * actual grid cells.
 *
 * @note The concept of @ref GlossSupport "support points" is
 * restricted to the finite element families based on Lagrange
 * interpolation. For a more general concept, see
 * @ref GlossGeneralizedSupport "generalized support points".
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossTargetComponent Target component</dt> <dd>When
 * vectors and matrices are grouped into blocks by component, it is
 * often desirable to collect several of the original components into
 * a single one. This could be for instance, grouping the velocities
 * of a Stokes system into a single block.</dd>
 *
 *
 * <dt class="glossary">@anchor GlossUnitCell Unit cell</dt>
 * <dd>See @ref GlossReferenceCell "Reference cell".</dd>
 *
 *
 * <dt class="glossary">@anchor GlossUnitSupport Unit support points</dt>
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
