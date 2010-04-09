//-------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005, 2006, 2007, 2008, 2009, 2010 by the deal.II authors
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
 * <dt class="glossary">@anchor GlossActive <b>Active cells</b></dt>
 * <dd>Mesh cells not refined any further in the hierarchy.</dd>
 *
 *
 * <dt class="glossary">@anchor GlossBlockLA <b>Block (linear algebra)</b></dt>

 * <dd>It is often convenient to treat a matrix or vector as a collection of
 * individual blocks. For example, in step-20 (and other tutorial
 * programs), we want to consider the global linear system $Ax=b$ in
 * the form
 * @f{eqnarray*}
  \left(\begin{array}{cc}
    M & B^T \\ B & 0
  \end{array}\right)
  \left(\begin{array}{cc}
    U \\ P
  \end{array}\right)
  =
  \left(\begin{array}{cc}
    F \\ G
  \end{array}\right),
 * @f}
 * where $U,P$ are the values of velocity and pressure degrees of freedom,
 * respectively, $M$ is the mass matrix on the velocity space, $B$ corresponds to
 * the negative divergence operator, and $B^T$ is its transpose and corresponds
 * to the negative gradient.
 *
 * Using such a decomposition into blocks, one can then define
 * preconditioners that are based on the individual operators that are
 * present in a system of equations (for example the Schur complement,
 * in the case of step-20), rather than the entire matrix. In essence,
 * blocks are used to reflect the structure of a PDE system in linear
 * algebra, in particular allowing for modular solvers for problems
 * with multiple solution components. On the other hand, the matrix
 * and right hand side vector can also treated as a unit, which is
 * convenient for example during assembly of the linear system when
 * one may not want to make a distinction between the individual
 * components, or for an outer Krylov space solver that doesn't care
 * about the block structure (e.g. if only the preconditioner needs
 * the block structure).
 *
 * Splitting matrices and vectors into blocks is supported by the
 * BlockSparseMatrix, BlockVector, and related classes. See the
 * overview of the various linear algebra classes in the @ref LAC
 * module. The objects present two interfaces: one that makes the
 * object look like a matrix or vector with global indexing
 * operations, and one that makes the object look like a collection of
 * sub-blocks that can be individually addressed. Depending on
 * context, one may wish to use one or the other interface.
 *
 * Typically, one defines the sub-structure of a matrix or vector by
 * grouping the degrees of freedom that make up groups of physical
 * quantities (for example all velocities) into individual blocks of
 * the linear system. This is defined in more detail below in the
 * glossary entry on @ref GlossBlock "Block (finite element)".
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossBlock <b>Block (finite element)</b></dt>
 * <dd>
 * <i>Intent:</i>
 * Blocks are a generalization of @ref GlossComponent "components" in that
 * they group together one or more components of a vector-valued finite
 * element that one would like to consider jointly. One often wants to do this
 * to define operators that correspond to the structure of a (part of a)
 * differential operator acting on the vector-valued solution, such as the
 * Schur complement solver in step-20, or the block solvers and
 * preconditioners of step-22.
 *
 * For the purpose of a discretization, blocks are the better concept to use
 * since it is not always possible to address individual components of a
 * solution. This is, in particular, the case for non-@ref GlossPrimitive
 * "primitive" elements. Take for instance the solution of the mixed Laplacian
 * system with the FE_RaviartThomas element (see step-20). There, the first
 * <tt>dim</tt> components are the directional velocities. Since the shape
 * functions are linear combinations of those, these <tt>dim</tt> components
 * constitute only a single block. On the other hand, the pressure variable is
 * scalar and would form a the second block, but in the <tt>dim+1</tt>st
 * component.
 *
 * The minimal size of each block is dictated by the underlying finite element
 * (a blocks consists of a single component for scalar elements, but in the
 * case of the FE_RaviartThomas, for example, a block consists of <tt>dim</tt>
 * components). However, several such minimal blocks can be grouped together
 * into user defined blocks at will, and in accordance with the
 * application. For instance, for the $Q_2^d\times Q^1$ (Taylor-Hood) Stokes
 * element, there are $d+1$ components each of which could in principle form
 * its own block. But we are typically more interested in having only two
 * blocks, one of which consists of all the velocity vector components
 * (i.e. this block would have $d$ components) and the other having only the
 * single pressure component.
 *
 * <i>Implementation:</i>
 * deal.II has a number of different finite element classes, all of which are
 * derived from the FiniteElement base class (see the @ref feall "module on
 * finite element classes"). With one exception, whether they are scalar or
 * vector valued, they all define a single block: all vector components the
 * finite element defines through its FiniteElement::n_components() function
 * form a single block, i.e. FiniteElement::n_blocks() returns one.
 *
 * The exception is the FESystem class that takes multiple simpler elements
 * and connects them into more complicated ones. Consequently, it can have
 * more than one block. A FESystem has as many blocks as it has base elements
 * times their multiplicity (see the constructors of FESystem to understand
 * this statement). In other words, it does not care how many blocks each base
 * element has, and consequently you can produce a Stokes element that has
 * only two blocks by creating the object
 * @code
 *    FESystem<dim> (FESystem<dim> (FE_Q<dim>(2), dim), 1,
 *                   FE_Q<dim>(1), 1);
 * @endcode
 * On the other hand, we could have produced a similar object with dim+1
 * blocks using
 * @code
 *    FESystem<dim> (FE_Q<dim>(2), dim,
 *                   FE_Q<dim>(1), 1);
 * @endcode
 * With the exception of the number of blocks, the two objects are the
 * same for all practical purposes, however.
 *
 * <i>Global degrees of freedom:</i>
 * While we have defined blocks above in terms of the vector components of a
 * vector-valued solution function (or, equivalently, in terms of the
 * vector-valued finite element space), every shape function of a finite
 * element is part of one block or another. Consequently, we can partition all
 * degrees of freedom defined on a DoFHandler into individual blocks. Since by
 * default the DoFHandler class enumerates degrees of freedom in a more or
 * less random way, you will first want to call the
 * DoFRenumbering::component_wise function to make sure that all degrees of
 * freedom that correspond to a single block are enumerated consecutively.
 *
 * If you do this, you naturally partition matrices and vectors into blocks as
 * well (see @ref GlossBlockLA "block (linear algebra)).  In most cases, when
 * you subdivide a matrix or vector into blocks, you do so by creating one
 * block for each block defined by the finite element (i.e. in most practical
 * cases the FESystem object). However, this needs not be so: the
 * DoFRenumbering::component_wise function allows to group several vector
 * components or finite element blocks into the same logical block (see, for
 * example, the @ref step_22 "step-22" or step-31 tutorial programs, as
 * opposed to step-20). As a consequence, using this feature, we can achieve
 * the same result, i.e. subdividing matrices into $2\times 2$ blocks and
 * vectors into 2 blocks, for the second way of creating a Stokes element
 * outlined above using an extra argument as we would have using the first way
 * of creating the Stokes element with two blocks right away.
 *
 * More information on this topic can be found in the documentation of
 * FESystem, the @ref vector_valued module and the tutorial programs
 * referenced therein.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossComponent <b>Component</b></dt>
 *
 * <dd> When considering systems of equations in which the solution is not
 * just a single scalar function, we say that we have a <i>vector system</i>
 * with a <i>vector-valued solution</i>. For example, the vector solution in
 * the elasticity equation considered in step-8 is $u=(u_x,u_y,u_z)^T$
 * consisting of the displacements in each of the three coordinate
 * directions. The solution then has three elements. Similarly, the 3d Stokes
 * equation considered in step-22 has four elements: $u=(v_x,v_y,v_z,p)^T$. We
 * call the elements of the vector-valued solution <i>components</i> in
 * deal.II. To be well-posed, for the solution to have $n$ components, there
 * need to be $n$ partial differential equations to describe them.
 *
 * In finite element programs, one frequently wants to address individual
 * elements (components) of this vector-valued solution, or sets of
 * components. For example, we do this extensively in step-8, and a lot
 * of documentation is also provided in the module on
 * @ref vector_valued "Handling vector valued problems". If you are thinking
 * only in terms of the partial differential equation (not in terms of
 * its discretization), then the concept of <i>components</i> is the natural
 * one.
 *
 * On the other hand, when talking about finite elements and degrees of
 * freedom, <i>components</i> are not always the correct concept because
 * components are not always individually addressable. In particular, this is
 * the case for @ref GlossPrimitive "non-primitive finite elements". Similarly,
 * one may not always <i>want</i> to address individual components but rather
 * sets of components &mdash; e.g. all velocity components together, and
 * separate from the pressure in the Stokes system, without further splitting
 * the velocities into their individual components. In either case, the
 * correct concept to think in is that of a @ref GlossBlock "block".  Since
 * each component, if individually addressable, is also a block, thinking in
 * terms of blocks is most frequently the better strategy.
 *
 * For a given finite element, the number of components can be queried using
 * the FiniteElementData::n_components() function. Individual components of a
 * shape function (if the element is primitive) can be queried using the
 * FiniteElement::shape_value_component() and
 * FiniteElement::shape_grad_component() functions on the reference cell. The
 * FEValues::shape_value_component() and FEValues::shape_grad_component()
 * functions do the same on a real cell. See also the documentation of the
 * FiniteElement and FEValues classes.</dd>
 *
 *
 * <dt class="glossary">@anchor GlossCompress <b>Compressing distributed
 *                                              vectors and matrices</b></dt>
 *
 * <dd>
 * For %parallel computations, deal.II uses the vector and matrix
 * classes defined in the PETScWrappers and TrilinosWrappers
 * namespaces. When running programs in %parallel using MPI, these
 * classes only store a certain number of rows or elements on the
 * current processor, whereas the rest of the vector or matrix is
 * stored on the other processors that belong to our MPI
 * universe. This presents a certain problem when you assemble linear
 * systems: we add elements to the matrix and right hand side vectors
 * that may or may not be stored locally. Sometimes, we may also want
 * to just <i>set</i> an element, not add to it.
 *
 * Both PETSc and Trilinos allow adding to or setting elements that
 * are not locally stored. In that case, they write the value that we
 * want to store or add into a cache, and we need to call one of the
 * functions TrilinosWrappers::VectorBase::compress(),
 * TrilinosWrappers::SparseMatrix::compress(),
 * PETScWrappers::VectorBase::compress() or
 * PETScWrappers::MatrixBase::compress() which will then ship the
 * values in the cache to the MPI process that owns the element to
 * which it is supposed to be added or written to. Due to the MPI
 * model that only allows to initiate communication from the sender
 * side (i.e. in particular, it is not a remote procedure call), these
 * functions are collective, i.e. they need to be called by all
 * processors.
 *
 * There is one snag, however: both PETSc and Trilinos need to know
 * whether the operation that these <code>compress()</code> functions
 * invoke applies to adding elements or setting them. Usually, you
 * will have written or added elements to the vector or matrix before
 * (and after <code>compress()</code> was last called), and in this
 * case the wrapper object knows that the global communication
 * operation is either an add or a set operation since it keeps track
 * of this sort of thing. However, there are cases where this isn't
 * so: for example, if you are working on a coarse grid and there are
 * more processors than coarse grid cells; in that case, some
 * processors will not assemble anything, and when they come to the
 * point where they call <code>compress()</code> on the system matrix
 * and right hand side, these objects are still in their pristine
 * state. In a case like this the wrapper object doesn't know whether
 * it is supposed to do a global exchange for add or set operations,
 * and in the worst case you end up with a deadlock (because those
 * processors that did assembly operations want to communicate, while
 * those that didn't assemble anything do not want to communicate).
 *
 * The way out of a situation like this is to use one of the two following
 * ways:
 * - You tell the object that you want to compress what operation is
 *   intended. The TrilinosWrappers::VectorBase::compress() can take such an
 *   additional argument. Or,
 * - You do a fake addition or set operation on the object in question.
 *
 * Some of the objects are also indifferent and can figure out what to
 * do without being told. The TrilinosWrappers::SparseMatrix can do that,
 * for example.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossDistorted <b>Distorted cells</b></dt>
 *
 * <dd>A <i>distorted cell</i> is a cell for which the mapping from
 * the reference cell to real cell has a Jacobian whose determinant is
 * non-positive somewhere in the cell. Typically, we only check the sign
 * of this determinant at the vertices of the cell. The function
 * GeometryInfo::alternating_form_at_vertices computes these
 * determinants at the vertices.
 *
 * By way of example, if all of the determinants are of roughly equal value
 * and on the order of $h^\text{dim}$ then the cell is well-shaped. For
 * example, a square cell or face has determinants equal to $h^\text{dim}$
 * whereas a strongly sheared parallelogram has a determinant much
 * smaller. Similarly, a cell with very unequal edge lengths will have widely
 * varying determinants. Conversely, a pinched cell in which the location of
 * two or more vertices is collapsed to a single point has a zero determinant
 * at this location. Finally, an inverted or twisted cell in which the
 * location of two vertices is out of order will have negative determinants.
 *
 * The following two images show a well-formed, a pinched, and a twisted cell
 * for both 2d and 3d:
 *
 * @image html distorted_2d.png "A well-formed, a pinched, and a twisted cell in 2d."
 *
 * @image html distorted_3d.png "A well-formed, a pinched, and a twisted cell in 3d."
 * </dd>
 *
 * Distorted cells can appear in two different ways: The original
 * coarse mesh can already contain such cells, or they can be created
 * as the result of mesh refinement if the boundary description in use
 * is sufficiently irregular.
 *
 * If the appropriate flag is given upon creation of a triangulation,
 * the function Triangulation::create_triangulation, which is called
 * by the various functions in GridGenerator and GridIn (but can also
 * be called from user code, see step-14, will signal
 * the creation of coarse meshes with distorted cells by throwing an
 * exception of type Triangulation::DistortedCellList. There are
 * legitimate cases for creating meshes with distorted cells (in
 * particular collapsed/pinched cells) if you don't intend to assemble
 * anything on these cells. For example, consider a case where one
 * would like to simulate the behavior of an elastic material with a
 * fluid-filled crack such as an oil reservoir. If the pressure
 * becomes too large, the crack is closed -- and the cells that
 * discretize the crack volume are collapsed to zero volume. As long
 * as you don't integrate over these cells to simulate the behavior of
 * the fluid (of which there isn't any if the crack has zero volume),
 * such meshes are perfectly legitimate. As a consequence,
 * Triangulation::create_triangulation does not simply abort the
 * program, but throws an exception that contains a list of cells that
 * are distorted; this exception can be caught and, if you believe
 * that you can ignore this condition, you can react by doing nothing
 * with the caught exception.
 *
 * The second case in which distorted cells can appear is through mesh
 * refinement when we have curved boundaries. Consider, for example, the
 * following case where the dashed line shows the exact boundary that the
 * lower edge of the cell is supposed to approximate (let's assume for
 * simplicity that the left, top and right edges are interior edges and
 * therefore will be considered as straight; in fact, for this particular case
 * in 2d where only one side of a cell is at the boundary we have special code
 * that avoids the situation depicted, but you will get the general idea of
 * the problem that holds in 3d or if more than one side of the cell is at the
 * boundary):
 *
 * @image html distorted_2d_refinement_01.png "One cell with an edge approximating a curved boundary"
 *
 * Now, if this cell is refined, we first split all edges and place
 * new mid-points on them. For the left, top and right edge, this is
 * trivial: because they are considered straight, we just take the
 * point in the middle between the two vertices. For the lower edge,
 * the Triangulation class asks the Boundary object associated with
 * this boundary (and in particular the Boundary::new_point_on_line
 * function) where the new point should lie. The four old vertices and
 * the four new points are shown here:
 *
 * @image html distorted_2d_refinement_02.png "Cell after edge refinement"
 *
 * The last step is to compute the location of the new point in the interior
 * of the cell. By default, it is chosen as the average location (arithmetic
 * mean of the coordinates) of these 8 points (in 3d, the 26 surrounding
 * points have different weights, but the idea is the same):
 *
 * @image html distorted_2d_refinement_03.png "Cell after edge refinement"
 *
 * The problem with that is, of course, that the bottom two child cells are
 * twisted, whereas the top two children are well-shaped. While such
 * meshes can happen with sufficiently irregular boundary descriptions
 * (and if the coarse mesh is entirely inadequate to resolve the
 * complexity of the boundary), the Triangulation class does not know
 * what to do in such situations. Consequently, the
 * Triangulation::execute_coarsening_and_refinement function does
 * create such meshes, but it keeps a list of cells whose children are
 * distorted. If this list is non-empty at the end of a refinement
 * step, it will throw an exception of type
 * Triangulation::DistortedCellList that contains those cells that
 * have distorted children. The caller of
 * Triangulation::execute_coarsening_and_refinement can then decide
 * what to do with this situation.
 *
 * One way to deal with this problem is to use the
 * GridTools::fix_up_distorted_child_cells function that attempts to
 * fix up exactly these cells if possible by moving around the node at
 * the center of the cell.
 *
 * Note that the Triangulation class does not test for the presence of
 * distorted cells by default, since the determination whether a cell
 * is distorted or not is not a cheap operation. If you want a
 * Triangulation object to test for distortion of cells, you need to
 * specify this upon creation of the object by passing the appropriate
 * flag.
 *
 *
 * <dt class="glossary">@anchor GlossFaceOrientation <b>Face orientation</b></dt>
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
 * <dt class="glossary">@anchor GlossGeneralizedSupport <b>Generalized support points</b></dt>
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
 * <dt class="glossary">@anchor hp_paper <b>%hp paper</b></dt>
 * <dd>The "hp paper" is a paper by W. Bangerth and O. Kayser-Herold, titled
 * "Data Structures and Requirements for hp Finite Element Software", that
 * describes many of the algorithms and data structures used in the implementation
 * of the hp framework of deal.II. In particular, it summarizes many of the
 * tricky points that have to be considered for %hp finite elements using continuous
 * elements.
 *
 * The full reference for this paper is as follows:
 * @code
Article{BK07,
  author = 	 {Wolfgang Bangerth and Oliver Kayser-Herold},
  title = 	 {Data Structures and Requirements for hp Finite Element
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
 * modified version of step-27. The main difference to that
 * tutorial program is that various operations in the program were timed for
 * the paper to compare different options and show that $hp$ methods are
 * really not all that expensive.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossInterpolation <b>Interpolation with finite elements</b></dt>
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
 * <dt class="glossary">@anchor GlossLagrange <b>Lagrange elements</b></dt>
 * <dd>Finite elements based on Lagrangian interpolation at
 * @ref GlossSupport "support points".</dd>
 *
 *
 * <dt class="glossary">@anchor mg_paper <b>%Multigrid paper</b></dt>
 * <dd>The "multigrid paper" is a paper by B. Janssen and G. Kanschat, titled
 * "Adaptive multilevel methods with local smoothing", that
 * describes many of the algorithms and data structures used in the implementation
 * of the multigrid framework of deal.II. It underlies the implementation of
 * the classes that are used in step-16 for multigrid
 * methods.
 *
 * The full reference for this paper is as follows:
 * @code
Article{JK10,
  author = 	 {B. Janssen and G. Kanschat},
  title = 	 {Adaptive multilevel methods with local smoothing},
  journal = 	 {submitted},
  year = 	 2010
}
 * @endcode
 * It is available as Technical Report IAMCS-2009-131 from the
 * <a href="http://iamcs.tamu.edu/research_sub.php?tab_sub=research&cms_id=8">Institute
 * for Applied Mathematics and Computational Science, Texas A&amp;M University</a>.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossNodes <b>Node values or node functionals</b></dt>
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
 * <dt class="glossary">@anchor GlossPrimitive <b>Primitive finite
 * elements</b></dt>
 * <dd>A finite element (described by its shape functions) is primitive if
 * there is a unique relation from shape function number to vector @ref
 * GlossComponent "component". What this means is that each shape function of
 * a vector-valued element has exactly one nonzero component if an element is
 * primitive. This includes, in particular, all scalar elements as well as
 * vector-valued elements assembled via the FESystem class from other
 * primitive (for example scalar) elements as shown in step-8, 
 * step-29, step-22 and several others. On the other hand,
 * the FE_RaviartThomas class used in step-20 and step-21, or the FE_Nedelec
 * class provide non-primitive finite elements because there, each
 * vector-value shape function may have several non-zero components.</dd>
 *
 * <dt class="glossary">@anchor GlossReferenceCell <b>Reference cell</b></dt>
 * <dd>The hypercube [0,1]<sup>dim</sup>, on which all parametric finite
 * element shape functions are defined.</dd>
 *
 *
 * <dt class="glossary">@anchor GlossShape <b>Shape functions</b></dt> <dd>The restriction of
 * the finite element basis functions to a single grid cell.</dd>
 *
 *
 * <dt class="glossary">@anchor GlossSubdomainId <b>Subdomain id</b></dt>
 * <dd>Each cell of a triangulation has associated with it a property called
 * the "subdomain id" that can be queried using a call like
 * <code>cell-@>subdomain_id()</code> and that can be set for example by using
 * <code>cell-@>set_subdomain_id(13)</code>. While in principle this property
 * can be used in any way application programs deem useful (it is simply an
 * integer associated with each cell that can indicate whatever you want), at
 * least for programs that run in %parallel it usually denotes the processor a
 * cell is associated with.
 *
 * For programs that are parallelized based on MPI but where each processor
 * stores the entire triangulation (as in, for example, step-18
 * or step-32, subdomain ids are assigned to cells by
 * partitioning a mesh, and each MPI process then only works on those cells it
 * "owns", i.e. that belong to a subdomain that it is associated with
 * (traditionally, this is the case for the subdomain id whose numerical value
 * coincides with the rank of the MPI process within the MPI
 * communicator). Partitioning is typically done using the
 * GridTools::partition() function, but any other method can also be used to
 * do this though most other ideas will likely lead to less well balanced
 * numbers of degrees of freedom on the various subdomains.
 * </dd>
 *
 *
 * <dt class="glossary">@anchor GlossSupport <b>Support points</b></dt> <dd>Support points are
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
 * <dt class="glossary">@anchor GlossTargetComponent <b>Target component</b></dt> <dd>When
 * vectors and matrices are grouped into blocks by component, it is
 * often desirable to collect several of the original components into
 * a single one. This could be for instance, grouping the velocities
 * of a Stokes system into a single block.</dd>
 *
 *
 * <dt class="glossary">@anchor GlossUnitCell <b>Unit cell</b></dt>
 * <dd>See @ref GlossReferenceCell "Reference cell".</dd>
 *
 *
 * <dt class="glossary">@anchor GlossUnitSupport <b>Unit support points</b></dt>
 * <dd>These are the @ref GlossSupport "support points" on the reference cell, defined in
 * FiniteElementBase. For example, the usual Q1 element in 1d has support
 * points  at <tt>x=0</tt> and <tt>x=1</tt> (and similarly, in higher
 * dimensions at the vertices of the unit square or cube). On the other
 * hand, higher order Lagrangian elements have unit support points also
 * in the interior of the unit line, square, or cube.
 * </dd>
 *
 * </dl>
 */
