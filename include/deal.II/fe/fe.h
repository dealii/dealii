// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef dealii_fe_h
#define dealii_fe_h

#include <deal.II/base/config.h>
#include <deal.II/fe/fe_base.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/block_mask.h>
#include <deal.II/fe/mapping.h>

#include <memory>


DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim> class FEValuesBase;
template <int dim, int spacedim> class FEValues;
template <int dim, int spacedim> class FEFaceValues;
template <int dim, int spacedim> class FESubfaceValues;
template <int dim, int spacedim> class FESystem;

/**
 * This is the base class for finite elements in arbitrary dimensions. It
 * declares the interface both in terms of member variables and public member
 * functions through which properties of a concrete implementation of a finite
 * element can be accessed. This interface generally consists of a number of
 * groups of variables and functions that can roughly be delineated as
 * follows:
 * - Basic information about the finite element, such as the number of
 * degrees of freedom per vertex, edge, or cell. This kind of data is stored
 * in the FiniteElementData base class. (Though the FiniteElement::get_name()
 * member function also falls into this category.)
 * - A description of the shape functions and their derivatives on the
 * reference cell $[0,1]^d$, if an element is indeed defined by mapping shape
 * functions from the reference cell to an actual cell.
 * - Matrices (and functions that access them) that describe how an
 * element's shape functions related to those on parent or child cells
 * (restriction or prolongation) or neighboring cells (for hanging node
 * constraints), as well as to other finite element spaces defined on the same
 * cell (e.g., when doing $p$ refinement).
 * - %Functions that describe the properties of individual shape functions,
 * for example which
 * @ref GlossComponent "vector components"
 * of a
 * @ref vector_valued "vector-valued finite element's"
 * shape function is nonzero, or whether an element is
 * @ref GlossPrimitive "primitive".
 * - For elements that are interpolatory, such as the common $Q_p$
 * Lagrange elements, data that describes where their
 * @ref GlossSupport "support points"
 * are located.
 * - %Functions that define the interface to the FEValues class that is
 * almost always used to access finite element shape functions from user code.
 *
 * The following sections discuss many of these concepts in more detail, and
 * outline strategies by which concrete implementations of a finite element
 * can provide the details necessary for a complete description of a finite
 * element space.
 *
 * As a general rule, there are three ways by which derived classes provide
 * this information:
 * - A number of fields that are generally easy to compute and that
 * are initialized by the constructor of this class (or the constructor of the
 * FiniteElementData base class) and derived classes therefore have to compute
 * in the process of calling this class's constructor. This is, specifically,
 * the case for the basic information and parts of the descriptive information
 * about shape functions mentioned above.
 * - Some common matrices that are widely used in the library and for
 * which this class provides protected member variables that the constructors
 * of derived classes need to fill. The purpose of providing these matrices in
 * this class is that (i) they are frequently used, and (ii) they are
 * expensive to compute. Consequently, it makes sense to only compute them
 * once, rather than every time they are used. In most cases, the constructor
 * of the current class already sets them to their correct size, and derived
 * classes therefore only have to fill them. Examples of this include the
 * matrices that relate the shape functions on one cell to the shape functions
 * on neighbors, children, and parents.
 * - Uncommon information, or information that depends on specific input
 * arguments, and that needs to be implemented by derived classes. For these,
 * this base class only declares abstract virtual member functions and derived
 * classes then have to implement them. Examples of this category would
 * include the functions that compute values and derivatives of shape
 * functions on the reference cell for which it is not possible to tabulate
 * values because there are infinitely many points at which one may want to
 * evaluate them. In some cases, derived classes may choose to simply not
 * implement <i>all</i> possible interfaces (or may not <i>yet</i> have a
 * complete implementation); for uncommon functions, there is then often a
 * member function derived classes can overload that describes whether a
 * particular feature is implemented. An example is whether an element
 * implements the information necessary to use it in the $hp$ finite element
 * context (see
 * @ref hp "hp finite element support").
 *
 *
 * <h3>Nomenclature</h3>
 *
 * Finite element classes have to define a large number of different
 * properties describing  a finite element space. The following subsections
 * describe some nomenclature that will be used in the documentation below.
 *
 * <h4>Components and blocks</h4>
 *
 * @ref vector_valued "Vector-valued finite element"
 * are elements used for systems of partial differential equations.
 * Oftentimes, they are composed via the FESystem class (which is itself
 * derived from the current class), but there are also non-composed elements
 * that have multiple components (for example the FE_Nedelec and
 * FE_RaviartThomas classes, among others). For any of these vector valued
 * elements, individual shape functions may be nonzero in one or several
 * @ref GlossComponent "components"
 * of the vector valued function. If the element is
 * @ref GlossPrimitive "primitive",
 * there is indeed a single component with a nonzero entry for each shape
 * function. This component can be determined using the
 * FiniteElement::system_to_component_index() function.
 *
 * On the other hand, if there is at least one shape function that is nonzero
 * in more than one vector component, then we call the entire element "non-
 * primitive". The FiniteElement::get_nonzero_components() can then be used to
 * determine which vector components of a shape function are nonzero. The
 * number of nonzero components of a shape function is returned by
 * FiniteElement::n_components(). Whether a shape function is non-primitive
 * can be queried by FiniteElement::is_primitive().
 *
 * Oftentimes, one may want to split linear system into blocks so that they
 * reflect the structure of the underlying operator. This is typically not
 * done based on vector components, but based on the use of
 * @ref GlossBlock "blocks",
 * and the result is then used to substructure objects of type BlockVector,
 * BlockSparseMatrix, BlockMatrixArray, and so on. If you use non-primitive
 * elements, you cannot determine the block number by
 * FiniteElement::system_to_component_index(). Instead, you can use
 * FiniteElement::system_to_block_index(). The number of blocks of a finite
 * element can be determined by FiniteElement::n_blocks().
 *
 * To better illustrate these concepts, let's consider the following
 * example of the multi-component system
 * @code
 * FESystem<dim> fe_basis(FE_Q<dim>(2), dim, FE_Q<dim>(1),1);
 * @endcode
 * with <code>dim=2</code>. The resulting finite element has 3 components:
 * two that come from the quadratic element and one from the linear element.
 * If, for example, this system were used to discretize a problem in fluid
 * dynamics then one could think of the first two components representing a
 * vector-valued velocity field whereas the last one corresponds to the scalar
 * pressure field. Without degree-of-freedom (DoF) renumbering this finite element will
 * produce the following distribution of local DoFs:
 *
 * @image html fe_system_example.png DoF indices
 *
 * Using the two functions FiniteElement::system_to_component_index() and
 * FiniteElement::system_to_base_index() one can get the
 * following information for each degree-of-freedom "i":
 * @code
 * const unsigned int component     = fe_basis.system_to_component_index(i).first;
 * const unsigned int within_base   = fe_basis.system_to_component_index(i).second;
 * const unsigned int base          = fe_basis.system_to_base_index(i).first.first;
 * const unsigned int multiplicity  = fe_basis.system_to_base_index(i).first.second;
 * const unsigned int within_base_  = fe_basis.system_to_base_index(i).second; // same as above
 * @endcode
 * which will result in:
 *
 * | DoF    | Component  | Base element | Shape function within base | Multiplicity |
 * | :----: | :--------: | :----------: | :------------------------: | :----------: |
 * |      0 |          0 |            0 |                          0 |            0 |
 * |      1 |          1 |            0 |                          0 |            1 |
 * |      2 |          2 |            1 |                          0 |            0 |
 * |      3 |          0 |            0 |                          1 |            0 |
 * |      4 |          1 |            0 |                          1 |            1 |
 * |      5 |          2 |            1 |                          1 |            0 |
 * |      6 |          0 |            0 |                          2 |            0 |
 * |      7 |          1 |            0 |                          2 |            1 |
 * |      8 |          2 |            1 |                          2 |            0 |
 * |      9 |          0 |            0 |                          3 |            0 |
 * |     10 |          1 |            0 |                          3 |            1 |
 * |     11 |          2 |            1 |                          3 |            0 |
 * |     12 |          0 |            0 |                          4 |            0 |
 * |     13 |          1 |            0 |                          4 |            1 |
 * |     14 |          0 |            0 |                          5 |            0 |
 * |     15 |          1 |            0 |                          5 |            1 |
 * |     16 |          0 |            0 |                          6 |            0 |
 * |     17 |          1 |            0 |                          6 |            1 |
 * |     18 |          0 |            0 |                          7 |            0 |
 * |     19 |          1 |            0 |                          7 |            1 |
 * |     20 |          0 |            0 |                          8 |            0 |
 * |     21 |          1 |            0 |                          8 |            1 |
 *
 * What we see is the following: there are a total of 22 degrees-of-freedom on this
 * element with components ranging from 0 to 2. Each DoF corresponds to
 * one of the two base elements used to build FESystem : $\mathbb Q_2$ or $\mathbb Q_1$.
 * Since FE_Q are primitive elements, we have a total of 9 distinct
 * scalar-valued shape functions for the quadratic element and 4 for the linear element.
 * Finally, for DoFs corresponding to the first base element multiplicity
 * is either zero or one, meaning that we use the same scalar valued $\mathbb Q_2$
 * for both $x$ and $y$ components of the velocity field $\mathbb Q_2 \otimes \mathbb Q_2$.
 * For DoFs corresponding to the second base element multiplicity is zero.
 *
 * <h4>Support points</h4>
 *
 * Finite elements are frequently defined by defining a polynomial space and a
 * set of dual functionals. If these functionals involve point evaluations,
 * then the element is "interpolatory" and it is possible to interpolate an
 * arbitrary (but sufficiently smooth) function onto the finite element space
 * by evaluating it at these points. We call these points "support points".
 *
 * Most finite elements are defined by mapping from the reference cell to a
 * concrete cell. Consequently, the support points are then defined on the
 * reference ("unit") cell, see
 * @ref GlossSupport "this glossary entry".
 * The support points on a concrete cell can then be computed by mapping the
 * unit support points, using the Mapping class interface and derived classes,
 * typically via the FEValues class.
 *
 * A typical code snippet to do so would look as follows:
 * @code
 * Quadrature<dim> dummy_quadrature (fe.get_unit_support_points());
 * FEValues<dim>   fe_values (mapping, fe, dummy_quadrature,
 *                            update_quadrature_points);
 * fe_values.reinit (cell);
 * Point<dim> mapped_point = fe_values.quadrature_point (i);
 * @endcode
 *
 * Alternatively, the points can be transformed one-by-one:
 * @code
 * const vector<Point<dim> > &unit_points =
 *    fe.get_unit_support_points();
 *
 * Point<dim> mapped_point =
 *    mapping.transform_unit_to_real_cell (cell, unit_points[i]);
 * @endcode
 *
 * @note Finite elements' implementation of the get_unit_support_points()
 * function returns these points in the same order as shape functions. As a
 * consequence, the quadrature points accessed above are also ordered in this
 * way. The order of shape functions is typically documented in the class
 * documentation of the various finite element classes.
 *
 *
 * <h3>Implementing finite element spaces in derived classes</h3>
 *
 * The following sections provide some more guidance for implementing concrete
 * finite element spaces in derived classes. This includes information that
 * depends on the dimension for which you want to provide something, followed
 * by a list of tools helping to generate information in concrete cases.
 *
 * It is important to note that there is a number of intermediate classes that
 * can do a lot of what is necessary for a complete description of finite
 * element spaces. For example, the FE_Poly, FE_PolyTensor, and FE_PolyFace
 * classes in essence build a complete finite element space if you only
 * provide them with an abstract description of the polynomial space upon
 * which you want to build an element. Using these intermediate classes
 * typically makes implementing finite element descriptions vastly simpler.
 *
 * As a general rule, if you want to implement an element, you will likely
 * want to look at the implementation of other, similar elements first. Since
 * many of the more complicated pieces of a finite element interface have to
 * do with how they interact with mappings, quadrature, and the FEValues
 * class, you will also want to read through the
 * @ref FE_vs_Mapping_vs_FEValues
 * documentation module.
 *
 *
 * <h4>Interpolation matrices in one dimension</h4>
 *
 * In one space dimension (i.e., for <code>dim==1</code> and any value of
 * <code>spacedim</code>), finite element classes implementing the interface
 * of the current base class need only set the #restriction and #prolongation
 * matrices that describe the interpolation of the finite element space on one
 * cell to that of its parent cell, and to that on its children, respectively.
 * The constructor of the current class in one dimension presets the
 * #interface_constraints matrix (used to describe hanging node constraints at
 * the interface between cells of different refinement levels) to have size
 * zero because there are no hanging nodes in 1d.
 *
 * <h4>Interpolation matrices in two dimensions</h4>
 *
 * In addition to the fields discussed above for 1D, a constraint matrix is
 * needed to describe hanging node constraints if the finite element has
 * degrees of freedom located on edges or vertices. These constraints are
 * represented by an $m\times n$-matrix #interface_constraints, where <i>m</i>
 * is the number of degrees of freedom on the refined side without the corner
 * vertices (those dofs on the middle vertex plus those on the two lines), and
 * <i>n</i> is that of the unrefined side (those dofs on the two vertices plus
 * those on the line). The matrix is thus a rectangular one. The $m\times n$
 * size of the #interface_constraints matrix can also be accessed through the
 * interface_constraints_size() function.
 *
 * The mapping of the dofs onto the indices of the matrix on the unrefined
 * side is as follows: let $d_v$ be the number of dofs on a vertex, $d_l$ that
 * on a line, then $n=0...d_v-1$ refers to the dofs on vertex zero of the
 * unrefined line, $n=d_v...2d_v-1$ to those on vertex one,
 * $n=2d_v...2d_v+d_l-1$ to those on the line.
 *
 * Similarly, $m=0...d_v-1$ refers to the dofs on the middle vertex of the
 * refined side (vertex one of child line zero, vertex zero of child line
 * one), $m=d_v...d_v+d_l-1$ refers to the dofs on child line zero,
 * $m=d_v+d_l...d_v+2d_l-1$ refers to the dofs on child line one.  Please note
 * that we do not need to reserve space for the dofs on the end vertices of
 * the refined lines, since these must be mapped one-to-one to the appropriate
 * dofs of the vertices of the unrefined line.
 *
 * Through this construction, the degrees of freedom on the child faces are
 * constrained to the degrees of freedom on the parent face. The information
 * so provided is typically consumed by the
 * DoFTools::make_hanging_node_constraints() function.
 *
 * @note The hanging node constraints described by these matrices are only
 * relevant to the case where the same finite element space is used on
 * neighboring (but differently refined) cells. The case that the finite
 * element spaces on different sides of a face are different, i.e., the $hp$
 * case (see
 * @ref hp "hp finite element support")
 * is handled by separate functions. See the
 * FiniteElement::get_face_interpolation_matrix() and
 * FiniteElement::get_subface_interpolation_matrix() functions.
 *
 *
 * <h4>Interpolation matrices in three dimensions</h4>
 *
 * For the interface constraints, the 3d case is similar to the 2d case. The
 * numbering for the indices $n$ on the mother face is obvious and keeps to
 * the usual numbering of degrees of freedom on quadrilaterals.
 *
 * The numbering of the degrees of freedom on the interior of the refined
 * faces for the index $m$ is as follows: let $d_v$ and $d_l$ be as above, and
 * $d_q$ be the number of degrees of freedom per quadrilateral (and therefore
 * per face), then $m=0...d_v-1$ denote the dofs on the vertex at the center,
 * $m=d_v...5d_v-1$ for the dofs on the vertices at the center of the bounding
 * lines of the quadrilateral, $m=5d_v..5d_v+4*d_l-1$ are for the degrees of
 * freedom on the four lines connecting the center vertex to the outer
 * boundary of the mother face, $m=5d_v+4*d_l...5d_v+4*d_l+8*d_l-1$ for the
 * degrees of freedom on the small lines surrounding the quad, and
 * $m=5d_v+12*d_l...5d_v+12*d_l+4*d_q-1$ for the dofs on the four child faces.
 * Note the direction of the lines at the boundary of the quads, as shown
 * below.
 *
 * The order of the twelve lines and the four child faces can be extracted
 * from the following sketch, where the overall order of the different dof
 * groups is depicted:
 * @verbatim
 *    *--15--4--16--*
 *    |      |      |
 *    10 19  6  20  12
 *    |      |      |
 *    1--7---0--8---2
 *    |      |      |
 *    9  17  5  18  11
 *    |      |      |
 *    *--13--3--14--*
 * @endverbatim
 * The numbering of vertices and lines, as well as the numbering of children
 * within a line is consistent with the one described in Triangulation.
 * Therefore, this numbering is seen from the outside and inside,
 * respectively, depending on the face.
 *
 * The three-dimensional case has a few pitfalls available for derived classes
 * that want to implement constraint matrices. Consider the following case:
 * @verbatim
 *          *-------*
 *         /       /|
 *        /       / |
 *       /       /  |
 *      *-------*   |
 *      |       |   *-------*
 *      |       |  /       /|
 *      |   1   | /       / |
 *      |       |/       /  |
 *      *-------*-------*   |
 *      |       |       |   *
 *      |       |       |  /
 *      |   2   |   3   | /
 *      |       |       |/
 *      *-------*-------*
 * @endverbatim
 * Now assume that we want to refine cell 2. We will end up with two faces
 * with hanging nodes, namely the faces between cells 1 and 2, as well as
 * between cells 2 and 3. Constraints have to be applied to the degrees of
 * freedom on both these faces. The problem is that there is now an edge (the
 * top right one of cell 2) which is part of both faces. The hanging node(s)
 * on this edge are therefore constrained twice, once from both faces. To be
 * meaningful, these constraints of course have to be consistent: both faces
 * have to constrain the hanging nodes on the edge to the same nodes on the
 * coarse edge (and only on the edge, as there can then be no constraints to
 * nodes on the rest of the face), and they have to do so with the same
 * weights. This is sometimes tricky since the nodes on the edge may have
 * different local numbers.
 *
 * For the constraint matrix this means the following: if a degree of freedom
 * on one edge of a face is constrained by some other nodes on the same edge
 * with some weights, then the weights have to be exactly the same as those
 * for constrained nodes on the three other edges with respect to the
 * corresponding nodes on these edges. If this isn't the case, you will get
 * into trouble with the ConstraintMatrix class that is the primary consumer
 * of the constraint information: while that class is able to handle
 * constraints that are entered more than once (as is necessary for the case
 * above), it insists that the weights are exactly the same.
 *
 * Using this scheme, child face degrees of freedom are constrained against
 * parent face degrees of freedom that contain those on the edges of the
 * parent face; it is possible that some of them are in turn constrained
 * themselves, leading to longer chains of constraints that the
 * ConstraintMatrix class will eventually have to sort out. (The constraints
 * described above are used by the DoFTools::make_hanging_node_constraints()
 * function that constructs a ConstraintMatrix object.) However, this is of no
 * concern for the FiniteElement and derived classes since they only act
 * locally on one cell and its immediate neighbor, and do not see the bigger
 * picture. The
 * @ref hp_paper
 * details how such chains are handled in practice.
 *
 *
 * <h4>Helper functions</h4>
 *
 * Construction of a finite element and computation of the matrices described
 * above is often a tedious task, in particular if it has to be performed for
 * several dimensions. Most of this work can be avoided by using the
 * intermediate classes already mentioned above (e.g., FE_Poly, FE_PolyTensor,
 * etc). Other tasks can be automated by some of the functions in namespace
 * FETools.
 *
 * <h5>Computing the correct basis from a set of linearly independent
 * functions</h5>
 *
 * First, it may already be difficult to compute the basis of shape functions
 * for arbitrary order and dimension. On the other hand, if the
 * @ref GlossNodes "node values"
 * are given, then the duality relation between node functionals and basis
 * functions defines the basis. As a result, the shape function space may be
 * defined from a set of linearly independent functions, such that the actual
 * finite element basis is computed from linear combinations of them. The
 * coefficients of these combinations are determined by the duality of node
 * values and form a matrix.
 *
 * Using this matrix allows the construction of the basis of shape functions
 * in two steps.
 * <ol>
 *
 * <li>Define the space of shape functions using an arbitrary basis
 * <i>w<sub>j</sub></i> and compute the matrix <i>M</i> of node functionals
 * <i>N<sub>i</sub></i> applied to these basis functions, such that its
 * entries are <i>m<sub>ij</sub> = N<sub>i</sub>(w<sub>j</sub>)</i>.
 *
 * <li>Compute the basis <i>v<sub>j</sub></i> of the finite element shape
 * function space by applying <i>M<sup>-1</sup></i> to the basis
 * <i>w<sub>j</sub></i>.
 * </ol>
 *
 * The matrix <i>M</i> may be computed using FETools::compute_node_matrix().
 * This function relies on the existence of #generalized_support_points and an
 * implementation of the FiniteElement::interpolate() function with
 * VectorSlice argument. (See the
 * @ref GlossGeneralizedSupport "glossary entry on generalized support points"
 * for more information.) With this, one can then use the following piece of
 * code in the constructor of a class derived from FiniteElement to compute the
 * $M$ matrix:
 * @code
 * FullMatrix<double> M(this->dofs_per_cell, this->dofs_per_cell);
 * FETools::compute_node_matrix(M, *this);
 * this->inverse_node_matrix.reinit(this->dofs_per_cell, this->dofs_per_cell);
 * this->inverse_node_matrix.invert(M);
 * @endcode
 * Don't forget to make sure that #unit_support_points or
 * #generalized_support_points are initialized before this!
 *
 * <h5>Computing prolongation matrices</h5>
 *
 * Once you have shape functions, you can define matrices that transfer data
 * from one cell to its children or the other way around. This is a common
 * operation in multigrid, of course, but is also used when interpolating the
 * solution from one mesh to another after mesh refinement, as well as in the
 * definition of some error estimators.
 *
 * To define the prolongation matrices, i.e., those matrices that describe the
 * transfer of a finite element field from one cell to its children,
 * implementations of finite elements can either fill the #prolongation array
 * by hand, or can call FETools::compute_embedding_matrices().
 *
 * In the latter case, all that is required is the following piece of code:
 * @code
 * for (unsigned int c=0; c<GeometryInfo<dim>::max_children_per_cell; ++c)
 *   this->prolongation[c].reinit (this->dofs_per_cell,
 *                                 this->dofs_per_cell);
 * FETools::compute_embedding_matrices (*this, this->prolongation);
 * @endcode
 * As in this example, prolongation is almost always implemented via
 * embedding, i.e., the nodal values of the function on the children may be
 * different from the nodal values of the function on the parent cell, but as
 * a function of $\mathbf x\in{\mathbb R}^\text{spacedim}$, the finite element
 * field on the child is the same as on the parent.
 *
 *
 * <h5>Computing restriction matrices</h5>
 *
 * The opposite operation, restricting a finite element function defined on
 * the children to the parent cell is typically implemented by interpolating
 * the finite element function on the children to the nodal values of the
 * parent cell. In deal.II, the restriction operation is implemented as a loop
 * over the children of a cell that each apply a matrix to the vector of
 * unknowns on that child cell (these matrices are stored in #restriction and
 * are accessed by get_restriction_matrix()). The operation that then needs to
 * be implemented turns out to be surprisingly difficult to describe, but is
 * instructive to describe because it also defines the meaning of the
 * #restriction_is_additive_flags array (accessed via the
 * restriction_is_additive() function).
 *
 * To give a concrete example, assume we use a $Q_1$ element in 1d, and that
 * on each of the parent and child cells degrees of freedom are (locally and
 * globally) numbered as follows:
 * @code
 * meshes:             *-------*                        *---*---*
 * local DoF numbers:  0       1                        0  1|0  1
 * global DoF numbers: 0       1                        0   1   2
 * @endcode
 * Then we want the restriction operation to take the value of the zeroth DoF
 * on child 0 as the value of the zeroth DoF on the parent, and take the value
 * of the first DoF on child 1 as the value of the first DoF on the parent.
 * Ideally, we would like to write this follows
 * @f[
 *   U^\text{coarse}|_\text{parent}
 *   = \sum_{\text{child}=0}^1 R_\text{child} U^\text{fine}|_\text{child}
 * @f]
 * where $U^\text{fine}|_\text{child=0}=(U^\text{fine}_0,U^\text{fine}_1)^T$
 * and $U^\text{fine}|_\text{child=1}=(U^\text{fine}_1,U^\text{fine}_2)^T$.
 * Writing the requested operation like this would here be possible by
 * choosing
 * @f[
 *   R_0 = \left(\begin{matrix}1 & 0 \\ 0 & 0\end{matrix}\right),
 *   \qquad\qquad
 *   R_1 = \left(\begin{matrix}0 & 0 \\ 0 & 1\end{matrix}\right).
 * @f]
 * However, this approach already fails if we go to a $Q_2$ element with the
 * following degrees of freedom:
 * @code
 * meshes:             *-------*                        *----*----*
 * local DoF numbers:  0   2   1                        0 2 1|0 2 1
 * global DoF numbers: 0   2   1                        0 2  1  4 3
 * @endcode
 * Writing things as the sum over matrix operations as above would not easily
 * work because we have to add nonzero values to $U^\text{coarse}_2$ twice,
 * once for each child.
 *
 * Consequently, restriction is typically implemented as a
 * <i>concatenation</i> operation. I.e., we first compute the individual
 * restrictions from each child,
 * @f[
 *   \tilde U^\text{coarse}_\text{child}
 *   = R_\text{child} U^\text{fine}|_\text{child},
 * @f]
 * and then compute the values of $U^\text{coarse}|_\text{parent}$ with the
 * following code:
 * @code
 * for (unsigned int child=0; child<cell->n_children(); ++child)
 *   for (unsigned int i=0; i<dofs_per_cell; ++i)
 *     if (U_tilde_coarse[child][i] != 0)
 *       U_coarse_on_parent[i] = U_tilde_coarse[child][i];
 * @endcode
 * In other words, each nonzero element of $\tilde
 * U^\text{coarse}_\text{child}$ <i>overwrites</i>, rather than adds to the
 * corresponding element of $U^\text{coarse}|_\text{parent}$. This typically
 * also implies that the restriction matrices from two different cells should
 * agree on a value for coarse degrees of freedom that they both want to touch
 * (otherwise the result would depend on the order in which we loop over
 * children, which would be unreasonable because the order of children is an
 * otherwise arbitrary convention). For example, in the example above, the
 * restriction matrices will be
 * @f[
 *   R_0 = \left(\begin{matrix}1 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 1 & 0 \end{matrix}\right),
 *   \qquad\qquad
 *   R_1 = \left(\begin{matrix}0 & 0 & 0 \\ 0 & 1 & 0 \\ 1 & 0 & 0 \end{matrix}\right),
 * @f]
 * and the compatibility condition is the $R_{0,21}=R_{1,20}$ because they
 * both indicate that $U^\text{coarse}|_\text{parent,2}$ should be set to one
 * times $U^\text{fine}|_\text{child=0,1}$ and
 * $U^\text{fine}|_\text{child=1,0}$.
 *
 * Unfortunately, not all finite elements allow to write the restriction
 * operation in this way. For example, for the piecewise constant FE_DGQ(0)
 * element, the value of the finite element field on the parent cell can not
 * be determined by interpolation from the children. Rather, the only
 * reasonable choice is to take it as the <i>average</i> value between the
 * children -- so we are back to the sum operation, rather than the
 * concatenation. Further thought shows that whether restriction should be
 * additive or not is a property of the individual shape function, not of the
 * finite element as a whole. Consequently, the
 * FiniteElement::restriction_is_additive() function returns whether a
 * particular shape function should act via concatenation (a return value of
 * @p false) or via addition (return value of @p true), and the correct code
 * for the overall operation is then as follows (and as, in fact, implemented
 * in DoFAccessor::get_interpolated_dof_values()):
 * @code
 * for (unsigned int child=0; child<cell->n_children(); ++child)
 *   for (unsigned int i=0; i<dofs_per_cell; ++i)
 *     if (fe.restriction_is_additive(i) == true)
 *       U_coarse_on_parent[i] += U_tilde_coarse[child][i];
 *     else
 *       if (U_tilde_coarse[child][i] != 0)
 *         U_coarse_on_parent[i] = U_tilde_coarse[child][i];
 * @endcode
 *
 *
 * <h5>Computing #interface_constraints</h5>
 *
 * Constraint matrices can be computed semi-automatically using
 * FETools::compute_face_embedding_matrices(). This function computes the
 * representation of the coarse mesh functions by fine mesh functions for each
 * child of a face separately. These matrices must be convoluted into a single
 * rectangular constraint matrix, eliminating degrees of freedom on common
 * vertices and edges as well as on the coarse grid vertices. See the
 * discussion above for details of this numbering.
 *
 * @ingroup febase fe
 *
 * @author Wolfgang Bangerth, Guido Kanschat, Ralf Hartmann, 1998, 2000, 2001,
 * 2005, 2015
 */
template <int dim, int spacedim=dim>
class FiniteElement : public Subscriptor,
  public FiniteElementData<dim>
{
public:
  /**
   * The dimension of the image space, corresponding to Triangulation.
   */
  static const unsigned int space_dimension = spacedim;

  /**
   * A base class for internal data that derived finite element classes may
   * wish to store.
   *
   * The class is used as follows: Whenever an FEValues (or FEFaceValues or
   * FESubfaceValues) object is initialized, it requests that the finite
   * element it is associated with creates an object of a class derived from
   * the current one here. This is done via each derived class's
   * FiniteElement::get_data() function. This object is then passed to the
   * FiniteElement::fill_fe_values(), FiniteElement::fill_fe_face_values(),
   * and FiniteElement::fill_fe_subface_values() functions as a constant
   * object. The intent of these objects is so that finite element classes can
   * pre-compute information once at the beginning (in the call to
   * FiniteElement::get_data() call) that can then be used on each cell that
   * is subsequently visited. An example for this is the values of shape
   * functions at the quadrature point of the reference cell, which remain the
   * same no matter the cell visited, and that can therefore be computed once
   * at the beginning and reused later on.
   *
   * Because only derived classes can know what they can pre-compute, each
   * derived class that wants to store information computed once at the
   * beginning, needs to derive its own InternalData class from this class,
   * and return an object of the derived type through its get_data() function.
   *
   * @author Guido Kanschat, 2001; Wolfgang Bangerth, 2015.
   */
  class InternalDataBase
  {
  private:
    /**
     * Copy construction is forbidden.
     */
    InternalDataBase (const InternalDataBase &) = delete;

  public:
    /**
     * Constructor. Sets update_flags to @p update_default and @p first_cell
     * to @p true.
     */
    InternalDataBase ();

    /**
     * Destructor. Made virtual to allow polymorphism.
     */
    virtual ~InternalDataBase () = default;

    /**
     * A set of update flags specifying the kind of information that an
     * implementation of the FiniteElement interface needs to compute on each
     * cell or face, i.e., in FiniteElement::fill_fe_values() and friends.
     *
     * This set of flags is stored here by implementations of
     * FiniteElement::get_data(), FiniteElement::get_face_data(), or
     * FiniteElement::get_subface_data(), and is that subset of the update
     * flags passed to those functions that require re-computation on every
     * cell. (The subset of the flags corresponding to information that can be
     * computed once and for all already at the time of the call to
     * FiniteElement::get_data() -- or an implementation of that interface --
     * need not be stored here because it has already been taken care of.)
     */
    UpdateFlags          update_each;

    /**
     * Return an estimate (in bytes) or the memory consumption of this object.
     */
    virtual std::size_t memory_consumption () const;
  };

public:
  /**
   * Constructor: initialize the fields of this base class of all finite
   * elements.
   *
   * @param[in] fe_data An object that stores identifying (typically integral)
   * information about the element to be constructed. In particular, this
   * object will contain data such as the number of degrees of freedom per
   * cell (and per vertex, line, etc), the number of vector components, etc.
   * This argument is used to initialize the base class of the current object
   * under construction.
   * @param[in] restriction_is_additive_flags A vector of size
   * <code>dofs_per_cell</code> (or of size one, see below) that for each
   * shape function states whether the shape function is additive or not. The
   * meaning of these flags is described in the section on restriction
   * matrices in the general documentation of this class.
   * @param[in] nonzero_components A vector of size <code>dofs_per_cell</code>
   * (or of size one, see below) that for each shape function provides a
   * ComponentMask (of size <code>fe_data.n_components()</code>) that
   * indicates in which vector components this shape function is nonzero
   * (after mapping the shape function to the real cell). For "primitive"
   * shape functions, this component mask will have a single entry (see
   * @ref GlossPrimitive
   * for more information about primitive elements). On the other hand, for
   * elements such as the Raviart-Thomas or Nedelec elements, shape functions
   * are nonzero in more than one vector component (after mapping to the real
   * cell) and the given component mask will contain more than one entry. (For
   * these two elements, all entries will in fact be set, but this would not
   * be the case if you couple a FE_RaviartThomas and a FE_Nedelec together
   * into a FESystem.)
   *
   * @pre <code>restriction_is_additive_flags.size() == dofs_per_cell</code>,
   * or <code>restriction_is_additive_flags.size() == 1</code>. In the latter
   * case, the array is simply interpreted as having size
   * <code>dofs_per_cell</code> where each element has the same value as the
   * single element given.
   *
   * @pre <code>nonzero_components.size() == dofs_per_cell</code>, or
   * <code>nonzero_components.size() == 1</code>. In the latter case, the
   * array is simply interpreted as having size <code>dofs_per_cell</code>
   * where each element equals the component mask provided in the single
   * element given.
   */
  FiniteElement (const FiniteElementData<dim>     &fe_data,
                 const std::vector<bool>          &restriction_is_additive_flags,
                 const std::vector<ComponentMask> &nonzero_components);

  /**
   * Move constructor.
   */
  FiniteElement (FiniteElement<dim, spacedim> &&) = default; // NOLINT

  /**
   * Copy constructor.
   */
  FiniteElement (const FiniteElement<dim, spacedim> &) = default;

  /**
   * Virtual destructor. Makes sure that pointers to this class are deleted
   * properly.
   */
  virtual ~FiniteElement () = default;

  /**
   * Creates information for creating a FESystem with this class as
   * base element and with multiplicity @p multiplicity. In particular,
   * the return type of this function can be used in the constructor
   * for a FESystem object.
   * This function calls clone() and hence creates a copy of the
   * current object.
   */
  std::pair<std::unique_ptr<FiniteElement<dim, spacedim> >, unsigned int>
  operator^ (const unsigned int multiplicity) const;

  /**
   * A sort of virtual copy constructor, this function returns a copy of
   * the finite element object. Derived classes need to override the function
   * here in this base class and return an object of the same type as the
   * derived class.
   *
   * Some places in the library, for
   * example the constructors of FESystem as well as the hp::FECollection
   * class, need to make copies of finite elements without knowing their exact
   * type. They do so through this function.
   */
  virtual
  std::unique_ptr<FiniteElement<dim,spacedim> >
  clone() const = 0;

  /**
   * Return a string that uniquely identifies a finite element. The general
   * convention is that this is the class name, followed by the dimension in
   * angle brackets, and the polynomial degree and whatever else is necessary
   * in parentheses. For example, <tt>FE_Q<2>(3)</tt> is the value returned
   * for a cubic element in 2d.
   *
   * Systems of elements have their own naming convention, see the FESystem
   * class.
   */
  virtual std::string get_name () const = 0;

  /**
   * This operator returns a reference to the present object if the argument
   * given equals to zero. While this does not seem particularly useful, it is
   * helpful in writing code that works with both ::DoFHandler and the hp
   * version hp::DoFHandler, since one can then write code like this:
   * @code
   *   dofs_per_cell
   *     = dof_handler->get_fe()[cell->active_fe_index()].dofs_per_cell;
   * @endcode
   *
   * This code doesn't work in both situations without the present operator
   * because DoFHandler::get_fe() returns a finite element, whereas
   * hp::DoFHandler::get_fe() returns a collection of finite elements that
   * doesn't offer a <code>dofs_per_cell</code> member variable: one first has
   * to select which finite element to work on, which is done using the
   * operator[]. Fortunately, <code>cell-@>active_fe_index()</code> also works
   * for non-hp classes and simply returns zero in that case. The present
   * operator[] accepts this zero argument, by returning the finite element
   * with index zero within its collection (that, of course, consists only of
   * the present finite element anyway).
   */
  const FiniteElement<dim,spacedim> &operator[] (const unsigned int fe_index) const;

  /**
   * @name Shape function access
   * @{
   */

  /**
   * Return the value of the @p ith shape function at the point @p p. @p p is
   * a point on the reference element. If the finite element is vector-valued,
   * then return the value of the only non-zero component of the vector value
   * of this shape function. If the shape function has more than one non-zero
   * component (which we refer to with the term non-primitive), then derived
   * classes implementing this function should throw an exception of type
   * ExcShapeFunctionNotPrimitive. In that case, use the
   * shape_value_component() function.
   *
   * Implementations of this function should throw an exception of type
   * ExcUnitShapeValuesDoNotExist if the shape functions of the FiniteElement
   * under consideration depend on the shape of the cell in real space, i.e.,
   * if the shape functions are not defined by mapping from the reference
   * cell. Some non-conforming elements are defined this way, as is the
   * FE_DGPNonparametric class, to name just one example.
   *
   * The default implementation of this virtual function does exactly this,
   * i.e., it simply throws an exception of type ExcUnitShapeValuesDoNotExist.
   */
  virtual double shape_value (const unsigned int  i,
                              const Point<dim>   &p) const;

  /**
   * Just like for shape_value(), but this function will be called when the
   * shape function has more than one non-zero vector component. In that case,
   * this function should return the value of the @p component-th vector
   * component of the @p ith shape function at point @p p.
   */
  virtual double shape_value_component (const unsigned int i,
                                        const Point<dim>   &p,
                                        const unsigned int component) const;

  /**
   * Return the gradient of the @p ith shape function at the point @p p. @p p
   * is a point on the reference element, and likewise the gradient is the
   * gradient on the unit cell with respect to unit cell coordinates. If the
   * finite element is vector-valued, then return the value of the only non-
   * zero component of the vector value of this shape function. If the shape
   * function has more than one non-zero component (which we refer to with the
   * term non-primitive), then derived classes implementing this function
   * should throw an exception of type ExcShapeFunctionNotPrimitive. In that
   * case, use the shape_grad_component() function.
   *
   * Implementations of this function should throw an exception of type
   * ExcUnitShapeValuesDoNotExist if the shape functions of the FiniteElement
   * under consideration depend on the shape of the cell in real space, i.e.,
   * if the shape functions are not defined by mapping from the reference
   * cell. Some non-conforming elements are defined this way, as is the
   * FE_DGPNonparametric class, to name just one example.
   *
   * The default implementation of this virtual function does exactly this,
   * i.e., it simply throws an exception of type ExcUnitShapeValuesDoNotExist.
   */
  virtual Tensor<1,dim> shape_grad (const unsigned int  i,
                                    const Point<dim>   &p) const;

  /**
   * Just like for shape_grad(), but this function will be called when the
   * shape function has more than one non-zero vector component. In that case,
   * this function should return the gradient of the @p component-th vector
   * component of the @p ith shape function at point @p p.
   */
  virtual Tensor<1,dim> shape_grad_component (const unsigned int i,
                                              const Point<dim>   &p,
                                              const unsigned int component) const;

  /**
   * Return the tensor of second derivatives of the @p ith shape function at
   * point @p p on the unit cell. The derivatives are derivatives on the unit
   * cell with respect to unit cell coordinates. If the finite element is
   * vector-valued, then return the value of the only non-zero component of
   * the vector value of this shape function. If the shape function has more
   * than one non-zero component (which we refer to with the term non-
   * primitive), then derived classes implementing this function should throw
   * an exception of type ExcShapeFunctionNotPrimitive. In that case, use the
   * shape_grad_grad_component() function.
   *
   * Implementations of this function should throw an exception of type
   * ExcUnitShapeValuesDoNotExist if the shape functions of the FiniteElement
   * under consideration depend on the shape of the cell in real space, i.e.,
   * if the shape functions are not defined by mapping from the reference
   * cell. Some non-conforming elements are defined this way, as is the
   * FE_DGPNonparametric class, to name just one example.
   *
   * The default implementation of this virtual function does exactly this,
   * i.e., it simply throws an exception of type ExcUnitShapeValuesDoNotExist.
   */
  virtual Tensor<2,dim> shape_grad_grad (const unsigned int  i,
                                         const Point<dim>   &p) const;

  /**
   * Just like for shape_grad_grad(), but this function will be called when
   * the shape function has more than one non-zero vector component. In that
   * case, this function should return the gradient of the @p component-th
   * vector component of the @p ith shape function at point @p p.
   */
  virtual Tensor<2,dim> shape_grad_grad_component (const unsigned int i,
                                                   const Point<dim>   &p,
                                                   const unsigned int component) const;

  /**
   * Return the tensor of third derivatives of the @p ith shape function at
   * point @p p on the unit cell. The derivatives are derivatives on the unit
   * cell with respect to unit cell coordinates. If the finite element is
   * vector-valued, then return the value of the only non-zero component of
   * the vector value of this shape function. If the shape function has more
   * than one non-zero component (which we refer to with the term non-
   * primitive), then derived classes implementing this function should throw
   * an exception of type ExcShapeFunctionNotPrimitive. In that case, use the
   * shape_3rd_derivative_component() function.
   *
   * Implementations of this function should throw an exception of type
   * ExcUnitShapeValuesDoNotExist if the shape functions of the FiniteElement
   * under consideration depend on the shape of the cell in real space, i.e.,
   * if the shape functions are not defined by mapping from the reference
   * cell. Some non-conforming elements are defined this way, as is the
   * FE_DGPNonparametric class, to name just one example.
   *
   * The default implementation of this virtual function does exactly this,
   * i.e., it simply throws an exception of type ExcUnitShapeValuesDoNotExist.
   */
  virtual Tensor<3,dim> shape_3rd_derivative (const unsigned int  i,
                                              const Point<dim>   &p) const;

  /**
   * Just like for shape_3rd_derivative(), but this function will be called
   * when the shape function has more than one non-zero vector component. In
   * that case, this function should return the gradient of the @p component-
   * th vector component of the @p ith shape function at point @p p.
   */
  virtual Tensor<3,dim> shape_3rd_derivative_component (const unsigned int i,
                                                        const Point<dim>   &p,
                                                        const unsigned int component) const;

  /**
   * Return the tensor of fourth derivatives of the @p ith shape function at
   * point @p p on the unit cell. The derivatives are derivatives on the unit
   * cell with respect to unit cell coordinates. If the finite element is
   * vector-valued, then return the value of the only non-zero component of
   * the vector value of this shape function. If the shape function has more
   * than one non-zero component (which we refer to with the term non-
   * primitive), then derived classes implementing this function should throw
   * an exception of type ExcShapeFunctionNotPrimitive. In that case, use the
   * shape_4th_derivative_component() function.
   *
   * Implementations of this function should throw an exception of type
   * ExcUnitShapeValuesDoNotExist if the shape functions of the FiniteElement
   * under consideration depend on the shape of the cell in real space, i.e.,
   * if the shape functions are not defined by mapping from the reference
   * cell. Some non-conforming elements are defined this way, as is the
   * FE_DGPNonparametric class, to name just one example.
   *
   * The default implementation of this virtual function does exactly this,
   * i.e., it simply throws an exception of type ExcUnitShapeValuesDoNotExist.
   */
  virtual Tensor<4,dim> shape_4th_derivative (const unsigned int  i,
                                              const Point<dim>   &p) const;

  /**
   * Just like for shape_4th_derivative(), but this function will be called
   * when the shape function has more than one non-zero vector component. In
   * that case, this function should return the gradient of the @p component-
   * th vector component of the @p ith shape function at point @p p.
   */
  virtual Tensor<4,dim> shape_4th_derivative_component (const unsigned int i,
                                                        const Point<dim>   &p,
                                                        const unsigned int component) const;
  /**
   * This function returns @p true, if the shape function @p shape_index has
   * non-zero function values somewhere on the face @p face_index. The
   * function is typically used to determine whether some matrix elements
   * resulting from face integrals can be assumed to be zero and may therefore
   * be omitted from integration.
   *
   * A default implementation is provided in this base class which always
   * returns @p true. This is the safe way to go.
   */
  virtual bool has_support_on_face (const unsigned int shape_index,
                                    const unsigned int face_index) const;

  //@}
  /**
   * @name Transfer and constraint matrices
   * @{
   */

  /**
   * Return the matrix that describes restricting a finite element field from
   * the given @p child (as obtained by the given @p refinement_case) to the
   * parent cell. The interpretation of the returned matrix depends on what
   * restriction_is_additive() returns for each shape function.
   *
   * Row and column indices are related to coarse grid and fine grid spaces,
   * respectively, consistent with the definition of the associated operator.
   *
   * If projection matrices are not implemented in the derived finite element
   * class, this function aborts with an exception of type
   * FiniteElement::ExcProjectionVoid. You can check whether this would happen
   * by first calling the restriction_is_implemented() or the
   * isotropic_restriction_is_implemented() function.
   */
  virtual const FullMatrix<double> &
  get_restriction_matrix (const unsigned int child,
                          const RefinementCase<dim> &refinement_case=RefinementCase<dim>::isotropic_refinement) const;

  /**
   * Prolongation/embedding matrix between grids.
   *
   * The identity operator from a coarse grid space into a fine grid space
   * (where both spaces are identified as functions defined on the parent and
   * child cells) is associated with a matrix @p P that maps the corresponding
   * representations of these functions in terms of their nodal values. The
   * restriction of this matrix @p P_i to a single child cell is returned
   * here.
   *
   * The matrix @p P is the concatenation, not the sum of the cell matrices @p
   * P_i. That is, if the same non-zero entry <tt>j,k</tt> exists in two
   * different child matrices @p P_i, the value should be the same in both
   * matrices and it is copied into the matrix @p P only once.
   *
   * Row and column indices are related to fine grid and coarse grid spaces,
   * respectively, consistent with the definition of the associated operator.
   *
   * These matrices are used by routines assembling the prolongation matrix
   * for multi-level methods.  Upon assembling the transfer matrix between
   * cells using this matrix array, zero elements in the prolongation matrix
   * are discarded and will not fill up the transfer matrix.
   *
   * If prolongation matrices are not implemented in the derived finite
   * element class, this function aborts with an exception of type
   * FiniteElement::ExcEmbeddingVoid. You can check whether this would happen
   * by first calling the prolongation_is_implemented() or the
   * isotropic_prolongation_is_implemented() function.
   */
  virtual const FullMatrix<double> &
  get_prolongation_matrix (const unsigned int child,
                           const RefinementCase<dim> &refinement_case=RefinementCase<dim>::isotropic_refinement) const;

  /**
   * Return whether this element implements its prolongation matrices. The
   * return value also indicates whether a call to the
   * get_prolongation_matrix() function will generate an error or not.
   *
   * Note, that this function returns <code>true</code> only if the
   * prolongation matrices of the isotropic and all anisotropic refinement
   * cases are implemented. If you are interested in the prolongation matrices
   * for isotropic refinement only, use the
   * isotropic_prolongation_is_implemented function instead.
   *
   * This function is mostly here in order to allow us to write more efficient
   * test programs which we run on all kinds of weird elements, and for which
   * we simply need to exclude certain tests in case something is not
   * implemented. It will in general probably not be a great help in
   * applications, since there is not much one can do if one needs these
   * features and they are not implemented. This function could be used to
   * check whether a call to <tt>get_prolongation_matrix()</tt> will succeed;
   * however, one then still needs to cope with the lack of information this
   * just expresses.
   */
  bool prolongation_is_implemented () const;

  /**
   * Return whether this element implements its prolongation matrices for
   * isotropic children. The return value also indicates whether a call to the
   * @p get_prolongation_matrix function will generate an error or not.
   *
   * This function is mostly here in order to allow us to write more efficient
   * test programs which we run on all kinds of weird elements, and for which
   * we simply need to exclude certain tests in case something is not
   * implemented. It will in general probably not be a great help in
   * applications, since there is not much one can do if one needs these
   * features and they are not implemented. This function could be used to
   * check whether a call to <tt>get_prolongation_matrix()</tt> will succeed;
   * however, one then still needs to cope with the lack of information this
   * just expresses.
   */
  bool isotropic_prolongation_is_implemented () const;

  /**
   * Return whether this element implements its restriction matrices. The
   * return value also indicates whether a call to the
   * get_restriction_matrix() function will generate an error or not.
   *
   * Note, that this function returns <code>true</code> only if the
   * restriction matrices of the isotropic and all anisotropic refinement
   * cases are implemented. If you are interested in the restriction matrices
   * for isotropic refinement only, use the
   * isotropic_restriction_is_implemented() function instead.
   *
   * This function is mostly here in order to allow us to write more efficient
   * test programs which we run on all kinds of weird elements, and for which
   * we simply need to exclude certain tests in case something is not
   * implemented. It will in general probably not be a great help in
   * applications, since there is not much one can do if one needs these
   * features and they are not implemented. This function could be used to
   * check whether a call to <tt>get_restriction_matrix()</tt> will succeed;
   * however, one then still needs to cope with the lack of information this
   * just expresses.
   */
  bool restriction_is_implemented () const;

  /**
   * Return whether this element implements its restriction matrices for
   * isotropic children. The return value also indicates whether a call to the
   * get_restriction_matrix() function will generate an error or not.
   *
   * This function is mostly here in order to allow us to write more efficient
   * test programs which we run on all kinds of weird elements, and for which
   * we simply need to exclude certain tests in case something is not
   * implemented. It will in general probably not be a great help in
   * applications, since there is not much one can do if one needs these
   * features and they are not implemented. This function could be used to
   * check whether a call to <tt>get_restriction_matrix()</tt> will succeed;
   * however, one then still needs to cope with the lack of information this
   * just expresses.
   */
  bool isotropic_restriction_is_implemented () const;


  /**
   * Access the #restriction_is_additive_flags field. See the discussion about
   * restriction matrices in the general class documentation for more
   * information.
   *
   * The index must be between zero and the number of shape functions of this
   * element.
   */
  bool restriction_is_additive (const unsigned int index) const;

  /**
   * Return a read only reference to the matrix that describes the constraints
   * at the interface between a refined and an unrefined cell.
   *
   * Some finite elements do not (yet) implement hanging node constraints. If
   * this is the case, then this function will generate an exception, since no
   * useful return value can be generated. If you should have a way to live
   * with this, then you might want to use the constraints_are_implemented()
   * function to check up front whether this function will succeed or generate
   * the exception.
   */
  const FullMatrix<double> &constraints (const dealii::internal::SubfaceCase<dim> &subface_case=dealii::internal::SubfaceCase<dim>::case_isotropic) const;

  /**
   * Return whether this element implements its hanging node constraints. The
   * return value also indicates whether a call to the constraints() function
   * will generate an error or not.
   *
   * This function is mostly here in order to allow us to write more efficient
   * test programs which we run on all kinds of weird elements, and for which
   * we simply need to exclude certain tests in case hanging node constraints
   * are not implemented. It will in general probably not be a great help in
   * applications, since there is not much one can do if one needs hanging
   * node constraints and they are not implemented. This function could be
   * used to check whether a call to <tt>constraints()</tt> will succeed;
   * however, one then still needs to cope with the lack of information this
   * just expresses.
   */
  bool constraints_are_implemented (const dealii::internal::SubfaceCase<dim> &subface_case=dealii::internal::SubfaceCase<dim>::case_isotropic) const;


  /**
   * Return whether this element implements its hanging node constraints in
   * the new way, which has to be used to make elements "hp compatible".  That
   * means, the element properly implements the get_face_interpolation_matrix
   * and get_subface_interpolation_matrix methods. Therefore the return value
   * also indicates whether a call to the get_face_interpolation_matrix()
   * method and the get_subface_interpolation_matrix() method will generate an
   * error or not.
   *
   * Currently the main purpose of this function is to allow the
   * make_hanging_node_constraints method to decide whether the new
   * procedures, which are supposed to work in the hp framework can be used,
   * or if the old well verified but not hp capable functions should be used.
   * Once the transition to the new scheme for computing the interface
   * constraints is complete, this function will be superfluous and will
   * probably go away.
   *
   * Derived classes should implement this function accordingly. The default
   * assumption is that a finite element does not provide hp capable face
   * interpolation, and the default implementation therefore returns @p false.
   */
  virtual bool hp_constraints_are_implemented () const;


  /**
   * Return the matrix interpolating from the given finite element to the
   * present one. The size of the matrix is then #dofs_per_cell times
   * <tt>source.#dofs_per_cell</tt>.
   *
   * Derived elements will have to implement this function. They may only
   * provide interpolation matrices for certain source finite elements, for
   * example those from the same family. If they don't implement interpolation
   * from a given element, then they must throw an exception of type
   * ExcInterpolationNotImplemented.
   */
  virtual void
  get_interpolation_matrix (const FiniteElement<dim,spacedim> &source,
                            FullMatrix<double>       &matrix) const;
  //@}

  /**
   * @name Functions to support hp
   * @{
   */


  /**
   * Return the matrix interpolating from a face of one element to the face
   * of the neighboring element.  The size of the matrix is then
   * <tt>source.#dofs_per_face</tt> times <tt>this->#dofs_per_face</tt>.
   *
   * Derived elements will have to implement this function. They may only
   * provide interpolation matrices for certain source finite elements, for
   * example those from the same family. If they don't implement interpolation
   * from a given element, then they must throw an exception of type
   * ExcInterpolationNotImplemented.
   */
  virtual void
  get_face_interpolation_matrix (const FiniteElement<dim,spacedim> &source,
                                 FullMatrix<double>       &matrix) const;


  /**
   * Return the matrix interpolating from a face of one element to the
   * subface of the neighboring element.  The size of the matrix is then
   * <tt>source.#dofs_per_face</tt> times <tt>this->#dofs_per_face</tt>.
   *
   * Derived elements will have to implement this function. They may only
   * provide interpolation matrices for certain source finite elements, for
   * example those from the same family. If they don't implement interpolation
   * from a given element, then they must throw an exception of type
   * ExcInterpolationNotImplemented.
   */
  virtual void
  get_subface_interpolation_matrix (const FiniteElement<dim,spacedim> &source,
                                    const unsigned int        subface,
                                    FullMatrix<double>       &matrix) const;
  //@}


  /**
   * @name Functions to support hp
   * @{
   */

  /**
   * If, on a vertex, several finite elements are active, the hp code first
   * assigns the degrees of freedom of each of these FEs different global
   * indices. It then calls this function to find out which of them should get
   * identical values, and consequently can receive the same global DoF index.
   * This function therefore returns a list of identities between DoFs of the
   * present finite element object with the DoFs of @p fe_other, which is a
   * reference to a finite element object representing one of the other finite
   * elements active on this particular vertex. The function computes which of
   * the degrees of freedom of the two finite element objects are equivalent,
   * both numbered between zero and the corresponding value of dofs_per_vertex
   * of the two finite elements. The first index of each pair denotes one of
   * the vertex dofs of the present element, whereas the second is the
   * corresponding index of the other finite element.
   */
  virtual
  std::vector<std::pair<unsigned int, unsigned int> >
  hp_vertex_dof_identities (const FiniteElement<dim,spacedim> &fe_other) const;

  /**
   * Same as hp_vertex_dof_indices(), except that the function treats degrees
   * of freedom on lines.
   */
  virtual
  std::vector<std::pair<unsigned int, unsigned int> >
  hp_line_dof_identities (const FiniteElement<dim,spacedim> &fe_other) const;

  /**
   * Same as hp_vertex_dof_indices(), except that the function treats degrees
   * of freedom on quads.
   */
  virtual
  std::vector<std::pair<unsigned int, unsigned int> >
  hp_quad_dof_identities (const FiniteElement<dim,spacedim> &fe_other) const;

  /**
   * Return whether this element dominates the one given as argument when they
   * meet at a common face, whether it is the other way around, whether
   * neither dominates, or if either could dominate.
   *
   * For a definition of domination, see FiniteElementDomination::Domination
   * and in particular the
   * @ref hp_paper "hp paper".
   */
  virtual
  FiniteElementDomination::Domination
  compare_for_face_domination (const FiniteElement<dim,spacedim> &fe_other) const;

  //@}

  /**
   * Comparison operator. We also check for equality of the constraint matrix,
   * which is quite an expensive operation.  Do therefore use this function
   * with care, if possible only for debugging purposes.
   *
   * Since this function is not that important, we avoid an implementational
   * question about comparing arrays and do not compare the matrix arrays
   * #restriction and #prolongation.
   */
  bool operator == (const FiniteElement<dim,spacedim> &) const;

  /**
   * @name Index computations
   * @{
   */
  /**
   * Compute vector component and index of this shape function within the
   * shape functions corresponding to this component from the index of a shape
   * function within this finite element.
   *
   * If the element is scalar, then the component is always zero, and the
   * index within this component is equal to the overall index.
   *
   * If the shape function referenced has more than one non-zero component,
   * then it cannot be associated with one vector component, and an exception
   * of type ExcShapeFunctionNotPrimitive will be raised.
   *
   * Note that if the element is composed of other (base) elements, and a base
   * element has more than one component but all its shape functions are
   * primitive (i.e. are non-zero in only one component), then this mapping
   * contains valid information. However, the index of a shape function of
   * this element within one component (i.e. the second number of the
   * respective entry of this array) does not indicate the index of the
   * respective shape function within the base element (since that has more
   * than one vector-component). For this information, refer to the
   * #system_to_base_table field and the system_to_base_index() function.
   *
   * See the class description above for an example of how this function is typically used.
   *
   * The use of this function is explained extensively in the step-8 and
   * @ref step_20 "step-20"
   * tutorial programs as well as in the
   * @ref vector_valued
   * module.
   */
  std::pair<unsigned int, unsigned int>
  system_to_component_index (const unsigned int index) const;

  /**
   * Compute the shape function for the given vector component and index.
   *
   * If the element is scalar, then the component must be zero, and the index
   * within this component is equal to the overall index.
   *
   * This is the opposite operation from the system_to_component_index()
   * function.
   */
  unsigned int component_to_system_index(const unsigned int component,
                                         const unsigned int index) const;

  /**
   * Same as system_to_component_index(), but do it for shape functions and
   * their indices on a face. The range of allowed indices is therefore
   * 0..#dofs_per_face.
   *
   * You will rarely need this function in application programs, since almost
   * all application codes only need to deal with cell indices, not face
   * indices. The function is mainly there for use inside the library.
   */
  std::pair<unsigned int, unsigned int>
  face_system_to_component_index (const unsigned int index) const;

  /**
   * For faces with non-standard face_orientation in 3D, the dofs on faces
   * (quads) have to be permuted in order to be combined with the correct
   * shape functions. Given a local dof @p index on a quad, return the local
   * index, if the face has non-standard face_orientation, face_flip or
   * face_rotation. In 2D and 1D there is no need for permutation and
   * consequently an exception is thrown.
   */
  unsigned int adjust_quad_dof_index_for_face_orientation (const unsigned int index,
                                                           const bool face_orientation,
                                                           const bool face_flip,
                                                           const bool face_rotation) const;

  /**
   * Given an index in the natural ordering of indices on a face, return the
   * index of the same degree of freedom on the cell.
   *
   * To explain the concept, consider the case where we would like to know
   * whether a degree of freedom on a face, for example as part of an FESystem
   * element, is primitive. Unfortunately, the is_primitive() function in the
   * FiniteElement class takes a cell index, so we would need to find the cell
   * index of the shape function that corresponds to the present face index.
   * This function does that.
   *
   * Code implementing this would then look like this:
   * @code
   * for (i=0; i<dofs_per_face; ++i)
   *  if (fe.is_primitive(fe.face_to_cell_index(i, some_face_no)))
   *   ... do whatever
   * @endcode
   * The function takes additional arguments that account for the fact that
   * actual faces can be in their standard ordering with respect to the cell
   * under consideration, or can be flipped, oriented, etc.
   *
   * @param face_dof_index The index of the degree of freedom on a face. This
   * index must be between zero and dofs_per_face.
   * @param face The number of the face this degree of freedom lives on. This
   * number must be between zero and GeometryInfo::faces_per_cell.
   * @param face_orientation One part of the description of the orientation of
   * the face. See
   * @ref GlossFaceOrientation.
   * @param face_flip One part of the description of the orientation of the
   * face. See
   * @ref GlossFaceOrientation.
   * @param face_rotation One part of the description of the orientation of
   * the face. See
   * @ref GlossFaceOrientation.
   * @return The index of this degree of freedom within the set of degrees of
   * freedom on the entire cell. The returned value will be between zero and
   * dofs_per_cell.
   *
   * @note This function exists in this class because that is where it was
   * first implemented. However, it can't really work in the most general case
   * without knowing what element we have. The reason is that when a face is
   * flipped or rotated, we also need to know whether we need to swap the
   * degrees of freedom on this face, or whether they are immune from this.
   * For this, consider the situation of a $Q_3$ element in 2d. If face_flip
   * is true, then we need to consider the two degrees of freedom on the edge
   * in reverse order. On the other hand, if the element were a $Q_1^2$, then
   * because the two degrees of freedom on this edge belong to different
   * vector components, they should not be considered in reverse order. What
   * all of this shows is that the function can't work if there are more than
   * one degree of freedom per line or quad, and that in these cases the
   * function will throw an exception pointing out that this functionality
   * will need to be provided by a derived class that knows what degrees of
   * freedom actually represent.
   */
  virtual
  unsigned int face_to_cell_index (const unsigned int face_dof_index,
                                   const unsigned int face,
                                   const bool face_orientation = true,
                                   const bool face_flip        = false,
                                   const bool face_rotation    = false) const;

  /**
   * For lines with non-standard line_orientation in 3D, the dofs on lines
   * have to be permuted in order to be combined with the correct shape
   * functions. Given a local dof @p index on a line, return the local index,
   * if the line has non-standard line_orientation. In 2D and 1D there is no
   * need for permutation, so the given index is simply returned.
   */
  unsigned int adjust_line_dof_index_for_line_orientation (const unsigned int index,
                                                           const bool line_orientation) const;

  /**
   * Return in which of the vector components of this finite element the @p
   * ith shape function is non-zero. The length of the returned array is equal
   * to the number of vector components of this element.
   *
   * For most finite element spaces, the result of this function will be a
   * vector with exactly one element being @p true, since for most spaces the
   * individual vector components are independent. In that case, the component
   * with the single zero is also the first element of what
   * system_to_component_index() returns.
   *
   * Only for those spaces that couple the components, for example to make a
   * shape function divergence free, will there be more than one @p true
   * entry.  Elements for which this is true are called non-primitive (see
   * @ref GlossPrimitive).
   */
  const ComponentMask &
  get_nonzero_components (const unsigned int i) const;

  /**
   * Return in how many vector components the @p ith shape function is non-
   * zero. This value equals the number of entries equal to @p true in the
   * result of the get_nonzero_components() function.
   *
   * For most finite element spaces, the result will be equal to one. It is
   * not equal to one only for those ansatz spaces for which vector-valued
   * shape functions couple the individual components, for example in order to
   * make them divergence-free.
   */
  unsigned int
  n_nonzero_components (const unsigned int i) const;

  /**
   * Return whether the entire finite element is primitive, in the sense that
   * all its shape functions are primitive. If the finite element is scalar,
   * then this is always the case.
   *
   * Since this is an extremely common operation, the result is cached and
   * returned by this function.
   */
  bool is_primitive () const;

  /**
   * Return whether the @p ith shape function is primitive in the sense that
   * the shape function is non-zero in only one vector component. Non-
   * primitive shape functions would then, for example, be those of divergence
   * free ansatz spaces, in which the individual vector components are
   * coupled.
   *
   * The result of the function is @p true if and only if the result of
   * <tt>n_nonzero_components(i)</tt> is equal to one.
   */
  bool
  is_primitive (const unsigned int i) const;

  /**
   * Number of base elements in a mixed discretization.
   *
   * Note that even for vector valued finite elements, the number of
   * components needs not coincide with the number of base elements, since
   * they may be reused. For example, if you create a FESystem with three
   * identical finite element classes by using the constructor that takes one
   * finite element and a multiplicity, then the number of base elements is
   * still one, although the number of components of the finite element is
   * equal to the multiplicity.
   */
  unsigned int n_base_elements () const;

  /**
   * Access to base element objects. If the element is atomic, then
   * <code>base_element(0)</code> is @p this.
   */
  virtual
  const FiniteElement<dim,spacedim> &
  base_element (const unsigned int index) const;

  /**
   * This index denotes how often the base element @p index is used in a
   * composed element. If the element is atomic, then the result is always
   * equal to one. See the documentation for the n_base_elements() function
   * for more details.
   */
  unsigned int
  element_multiplicity (const unsigned int index) const;

  /**
   * Return a reference to a contained finite element that matches the components
   * selected by the given ComponentMask @p mask.
   *
   * For an arbitrarily nested FESystem, this function returns the inner-most
   * FiniteElement that matches the given mask. The method fails if the @p mask
   * does not exactly match one of the contained finite elements. It is most
   * useful if the current object is an FESystem, as the return value can
   * only be @p this in all other cases.
   *
   * Note that the returned object can be an FESystem if the
   * mask matches it but not any of the contained objects.
   *
   * Let us illustrate the function with the an FESystem @p fe with 7 components:
   * @code
   * FESystem<2> fe_velocity(FE_Q<2>(2), 2);
   * FE_Q<2> fe_pressure(1);
   * FE_DGP<2> fe_dg(0);
   * FE_BDM<2> fe_nonprim(1);
   * FESystem<2> fe(fe_velocity, 1, fe_pressure, 1, fe_dg, 2, fe_nonprim, 1);
   * @endcode
   *
   * The following table lists all possible component masks you can use:
   * <table>
   * <tr>
   * <th>ComponentMask</th>
   * <th>Result</th>
   * <th>Description</th>
   * </tr>
   * <tr>
   * <td><code>[true,true,true,true,true,true,true]</code></td>
   * <td><code>FESystem<2>[FESystem<2>[FE_Q<2>(2)^2]-FE_Q<2>(1)-FE_DGP<2>(0)^2-FE_BDM<2>(1)]</code></td>
   * <td>@p fe itself, the whole @p FESystem</td>
   * </tr>
   * <tr>
   * <td><code>[true,true,false,false,false,false,false]</code></td>
   * <td><code>FESystem<2>[FE_Q<2>(2)^2]</code></td>
   * <td>just the @p fe_velocity</td>
   * </tr>
   * <tr>
   * <td><code>[true,false,false,false,false,false,false]</code></td>
   * <td><code>FE_Q<2>(2)</code></td>
   * <td>The first component in @p fe_velocity</td>
   * </tr>
   * <tr>
   * <td><code>[false,true,false,false,false,false,false]</code></td>
   * <td><code>FE_Q<2>(2)</code></td>
   * <td>The second component in @p fe_velocity</td>
   * </tr>
   * <tr>
   * <td><code>[false,false,true,false,false,false,false]</code></td>
   * <td><code>FE_Q<2>(1)</code></td>
   * <td>@p fe_pressure</td>
   * </tr>
   * <tr>
   * <td><code>[false,false,false,true,false,false,false]</code></td>
   * <td><code>FE_DGP<2>(0)</code></td>
   * <td>first copy of @p fe_dg</td>
   * </tr>
   * <tr>
   * <td><code>[false,false,false,false,true,false,false]</code></td>
   * <td><code>FE_DGP<2>(0)</code></td>
   * <td>second copy of @p fe_dg</td>
   * </tr>
   * <tr>
   * <td><code>[false,false,false,false,false,true,true]</code></td>
   * <td><code>FE_BDM<2>(1)</code></td>
   * <td>both components of @p fe_nonprim</td>
   * </tr>
   * </table>
   */
  const FiniteElement<dim,spacedim> &
  get_sub_fe (const ComponentMask &mask) const;

  /**
   * Return a reference to a contained finite element that matches the components
   * @p n_selected_components components starting at component with index
   * @p first_component.
   *
   * See the other get_sub_fe() function above for more details.
   */
  virtual
  const FiniteElement<dim,spacedim> &
  get_sub_fe (const unsigned int first_component,
              const unsigned int n_selected_components) const;

  /**
   * Return for shape function @p index the base element it belongs to, the
   * number of the copy of this base element (which is between zero and the
   * multiplicity of this element), and the index of this shape function
   * within this base element.
   *
   * If the element is not composed of others, then base and instance are
   * always zero, and the index is equal to the number of the shape function.
   * If the element is composed of single instances of other elements (i.e.
   * all with multiplicity one) all of which are scalar, then base values and
   * dof indices within this element are equal to the
   * #system_to_component_table. It differs only in case the element is
   * composed of other elements and at least one of them is vector-valued
   * itself.
   *
   * See the class documentation above for an example of how this function is typically used.
   *
   * This function returns valid values also in the case of vector-valued
   * (i.e. non-primitive) shape functions, in contrast to the
   * system_to_component_index() function.
   */
  std::pair<std::pair<unsigned int, unsigned int>, unsigned int>
  system_to_base_index (const unsigned int index) const;

  /**
   * Same as system_to_base_index(), but for degrees of freedom located on a
   * face. The range of allowed indices is therefore 0..#dofs_per_face.
   *
   * You will rarely need this function in application programs, since almost
   * all application codes only need to deal with cell indices, not face
   * indices. The function is mainly there for use inside the library.
   */
  std::pair<std::pair<unsigned int, unsigned int>, unsigned int>
  face_system_to_base_index (const unsigned int index) const;

  /**
   * Given a base element number, return the first block of a BlockVector it
   * would generate.
   */
  types::global_dof_index first_block_of_base (const unsigned int b) const;

  /**
   * For each vector component, return which base element implements this
   * component and which vector component in this base element this is. This
   * information is only of interest for vector-valued finite elements which
   * are composed of several sub-elements. In that case, one may want to
   * obtain information about the element implementing a certain vector
   * component, which can be done using this function and the
   * FESystem::base_element() function.
   *
   * If this is a scalar finite element, then the return value is always equal
   * to a pair of zeros.
   */
  std::pair<unsigned int, unsigned int>
  component_to_base_index (const unsigned int component) const;


  /**
   * Return the base element for this block and the number of the copy of the
   * base element.
   */
  std::pair<unsigned int,unsigned int>
  block_to_base_index (const unsigned int block) const;

  /**
   * The vector block and the index inside the block for this shape function.
   */
  std::pair<unsigned int,types::global_dof_index>
  system_to_block_index (const unsigned int component) const;

  /**
   * The vector block for this component.
   */
  unsigned int
  component_to_block_index (const unsigned int component) const;

  //@}

  /**
   * @name Component and block matrices
   * @{
   */

  /**
   * Return a component mask with as many elements as this object has vector
   * components and of which exactly the one component is true that
   * corresponds to the given argument. See
   * @ref GlossComponentMask "the glossary"
   * for more information.
   *
   * @param scalar An object that represents a single scalar vector component
   * of this finite element.
   * @return A component mask that is false in all components except for the
   * one that corresponds to the argument.
   */
  ComponentMask
  component_mask (const FEValuesExtractors::Scalar &scalar) const;

  /**
   * Return a component mask with as many elements as this object has vector
   * components and of which exactly the <code>dim</code> components are true
   * that correspond to the given argument. See
   * @ref GlossComponentMask "the glossary"
   * for more information.
   *
   * @param vector An object that represents dim vector components of this
   * finite element.
   * @return A component mask that is false in all components except for the
   * ones that corresponds to the argument.
   */
  ComponentMask
  component_mask (const FEValuesExtractors::Vector &vector) const;

  /**
   * Return a component mask with as many elements as this object has vector
   * components and of which exactly the <code>dim*(dim+1)/2</code> components
   * are true that correspond to the given argument. See
   * @ref GlossComponentMask "the glossary"
   * for more information.
   *
   * @param sym_tensor An object that represents dim*(dim+1)/2 components of
   * this finite element that are jointly to be interpreted as forming a
   * symmetric tensor.
   * @return A component mask that is false in all components except for the
   * ones that corresponds to the argument.
   */
  ComponentMask
  component_mask (const FEValuesExtractors::SymmetricTensor<2> &sym_tensor) const;

  /**
   * Given a block mask (see
   * @ref GlossBlockMask "this glossary entry"),
   * produce a component mask (see
   * @ref GlossComponentMask "this glossary entry")
   * that represents the components that correspond to the blocks selected in
   * the input argument. This is essentially a conversion operator from
   * BlockMask to ComponentMask.
   *
   * @param block_mask The mask that selects individual blocks of the finite
   * element
   * @return A mask that selects those components corresponding to the
   * selected blocks of the input argument.
   */
  ComponentMask
  component_mask (const BlockMask &block_mask) const;

  /**
   * Return a block mask with as many elements as this object has blocks and
   * of which exactly the one component is true that corresponds to the given
   * argument. See
   * @ref GlossBlockMask "the glossary"
   * for more information.
   *
   * @note This function will only succeed if the scalar referenced by the
   * argument encompasses a complete block. In other words, if, for example,
   * you pass an extractor for the single $x$ velocity and this object
   * represents an FE_RaviartThomas object, then the single scalar object you
   * selected is part of a larger block and consequently there is no block
   * mask that would represent it. The function will then produce an
   * exception.
   *
   * @param scalar An object that represents a single scalar vector component
   * of this finite element.
   * @return A component mask that is false in all components except for the
   * one that corresponds to the argument.
   */
  BlockMask
  block_mask (const FEValuesExtractors::Scalar &scalar) const;

  /**
   * Return a component mask with as many elements as this object has vector
   * components and of which exactly the <code>dim</code> components are true
   * that correspond to the given argument. See
   * @ref GlossBlockMask "the glossary"
   * for more information.
   *
   * @note The same caveat applies as to the version of the function above:
   * The extractor object passed as argument must be so that it corresponds to
   * full blocks and does not split blocks of this element.
   *
   * @param vector An object that represents dim vector components of this
   * finite element.
   * @return A component mask that is false in all components except for the
   * ones that corresponds to the argument.
   */
  BlockMask
  block_mask (const FEValuesExtractors::Vector &vector) const;

  /**
   * Return a component mask with as many elements as this object has vector
   * components and of which exactly the <code>dim*(dim+1)/2</code> components
   * are true that correspond to the given argument. See
   * @ref GlossBlockMask "the glossary"
   * for more information.
   *
   * @note The same caveat applies as to the version of the function above:
   * The extractor object passed as argument must be so that it corresponds to
   * full blocks and does not split blocks of this element.
   *
   * @param sym_tensor An object that represents dim*(dim+1)/2 components of
   * this finite element that are jointly to be interpreted as forming a
   * symmetric tensor.
   * @return A component mask that is false in all components except for the
   * ones that corresponds to the argument.
   */
  BlockMask
  block_mask (const FEValuesExtractors::SymmetricTensor<2> &sym_tensor) const;

  /**
   * Given a component mask (see
   * @ref GlossComponentMask "this glossary entry"),
   * produce a block mask (see
   * @ref GlossBlockMask "this glossary entry")
   * that represents the blocks that correspond to the components selected in
   * the input argument. This is essentially a conversion operator from
   * ComponentMask to BlockMask.
   *
   * @note This function will only succeed if the components referenced by the
   * argument encompasses complete blocks. In other words, if, for example,
   * you pass an component mask for the single $x$ velocity and this object
   * represents an FE_RaviartThomas object, then the single component you
   * selected is part of a larger block and consequently there is no block
   * mask that would represent it. The function will then produce an
   * exception.
   *
   * @param component_mask The mask that selects individual components of the
   * finite element
   * @return A mask that selects those blocks corresponding to the selected
   * blocks of the input argument.
   */
  BlockMask
  block_mask (const ComponentMask &component_mask) const;

  /**
   * Return a list of constant modes of the element. The number of rows in
   * the resulting table depends on the elements in use. For standard
   * elements, the table has as many rows as there are components in the
   * element and dofs_per_cell columns. To each component of the finite
   * element, the row in the returned table contains a basis representation of
   * the constant function 1 on the element. However, there are some scalar
   * elements where there is more than one constant mode, e.g. the element
   * FE_Q_DG0.
   *
   * In order to match the constant modes to the actual components in the
   * element, the returned data structure also returns a vector with as many
   * components as there are constant modes on the element that contains the
   * component number.
   */
  virtual std::pair<Table<2,bool>,std::vector<unsigned int> >
  get_constant_modes () const;

  //@}

  /**
   * @name Support points and interpolation
   * @{
   */

  /**
   * Return the support points of the trial functions on the unit cell, if the
   * derived finite element defines them.  Finite elements that allow some
   * kind of interpolation operation usually have support points. On the other
   * hand, elements that define their degrees of freedom by, for example,
   * moments on faces, or as derivatives, don't have support points. In that
   * case, the returned field is empty.
   *
   * If the finite element defines support points, then their number equals
   * the number of degrees of freedom of the element.  The order of points in
   * the array matches that returned by the <tt>cell->get_dof_indices</tt>
   * function.
   *
   * See the class documentation for details on support points.
   *
   * @note Finite elements' implementation of this function returns these
   * points in the same order as shape functions. The order of shape functions
   * is typically documented in the class documentation of the various finite
   * element classes. In particular, shape functions (and consequently the
   * mapped quadrature points discussed in the class documentation of this
   * class) will then traverse first those shape functions located on
   * vertices, then on lines, then on quads, etc.
   *
   * @note If this element implements support points, then it will return one
   * such point per shape function. Since multiple shape functions may be
   * defined at the same location, the support points returned here may be
   * duplicated. An example would be an element of the kind
   * <code>FESystem(FE_Q(1),3)</code> for which each support point would
   * appear three times in the returned array.
   */
  const std::vector<Point<dim> > &
  get_unit_support_points () const;

  /**
   * Return whether a finite element has defined support points. If the result
   * is true, then a call to the get_unit_support_points() yields a non-empty
   * array.
   *
   * The result may be false if an element is not defined by interpolating
   * shape functions, for example by P-elements on quadrilaterals. It will
   * usually only be true if the element constructs its shape functions by the
   * requirement that they be one at a certain point and zero at all the
   * points associated with the other shape functions.
   *
   * In composed elements (i.e. for the FESystem class), the result will be
   * true if all the base elements have defined support points. FE_Nothing
   * is a special case in FESystems, because it has 0 support points and
   * has_support_points() is false, but an FESystem containing an FE_Nothing
   * among other elements will return true.
   */
  bool has_support_points () const;

  /**
   * Return the position of the support point of the @p indexth shape
   * function. If it does not exist, raise an exception.
   *
   * The default implementation simply returns the respective element from the
   * array you get from get_unit_support_points(), but derived elements may
   * overload this function. In particular, note that the FESystem class
   * overloads it so that it can return the support points of individual base
   * elements, if not all the base elements define support points. In this
   * way, you can still ask for certain support points, even if
   * get_unit_support_points() only returns an empty array.
   */
  virtual
  Point<dim>
  unit_support_point (const unsigned int index) const;

  /**
   * Return the support points of the trial functions on the unit face, if the
   * derived finite element defines some.  Finite elements that allow some
   * kind of interpolation operation usually have support points. On the other
   * hand, elements that define their degrees of freedom by, for example,
   * moments on faces, or as derivatives, don't have support points. In that
   * case, the returned field is empty
   *
   * Note that elements that have support points need not necessarily have
   * some on the faces, even if the interpolation points are located
   * physically on a face. For example, the discontinuous elements have
   * interpolation points on the vertices, and for higher degree elements also
   * on the faces, but they are not defined to be on faces since in that case
   * degrees of freedom from both sides of a face (or from all adjacent
   * elements to a vertex) would be identified with each other, which is not
   * what we would like to have). Logically, these degrees of freedom are
   * therefore defined to belong to the cell, rather than the face or vertex.
   * In that case, the returned element would therefore have length zero.
   *
   * If the finite element defines support points, then their number equals
   * the number of degrees of freedom on the face (#dofs_per_face). The order
   * of points in the array matches that returned by the
   * <tt>cell->face(face)->get_dof_indices</tt> function.
   *
   * See the class documentation for details on support points.
   */
  const std::vector<Point<dim-1> > &
  get_unit_face_support_points () const;

  /**
   * Return whether a finite element has defined support points on faces. If
   * the result is true, then a call to the get_unit_face_support_points()
   * yields a non-empty vector.
   *
   * For more information, see the documentation for the has_support_points()
   * function.
   */
  bool has_face_support_points () const;

  /**
   * The function corresponding to the unit_support_point() function, but for
   * faces. See there for more information.
   */
  virtual
  Point<dim-1>
  unit_face_support_point (const unsigned int index) const;

  /**
   * Return a vector of generalized support points.
   *
   * @note The vector returned by this function is always a minimal set of
   * *unique* support points. This is in contrast to the behavior of
   * get_unit_support_points() that returns a repeated list of unit support
   * points for an FESystem of numerous (Lagrangian) base elements.
   *
   * See the
   * @ref GlossGeneralizedSupport "glossary entry on generalized support points"
   * for more information.
   */
  const std::vector<Point<dim> > &
  get_generalized_support_points () const;

  /**
   * Return whether a finite element has defined generalized support
   * points. If the result is true, then a call to the
   * get_generalized_support_points() yields a non-empty vector.
   *
   * See the
   * @ref GlossGeneralizedSupport "glossary entry on generalized support points"
   * for more information.
   */
  bool has_generalized_support_points () const;

  /**
   * Return the equivalent to get_generalized_support_points(), except
   * for faces.
   *
   * @deprecated In general, it is not possible to associate a unique
   * subset of generalized support points describing degrees of freedom for
   * a given face. Don't use this function
   */
  DEAL_II_DEPRECATED
  const std::vector<Point<dim-1> > &
  get_generalized_face_support_points () const;

  /**
   * Return whether a finite element has defined generalized support points on
   * faces. If the result is true, then a call to the
   * get_generalized_face_support_points() function yields a non-empty array.
   *
   * For more information, see the documentation for the has_support_points()
   * function.
   *
   * @deprecated In general, it is not possible to associate a unique
   * subset of generalized support points describing degrees of freedom for
   * a given face. Don't use this function
   */
  DEAL_II_DEPRECATED
  bool
  has_generalized_face_support_points () const;

  /**
   * For a given degree of freedom, return whether it is logically associated
   * with a vertex, line, quad or hex.
   *
   * For instance, for continuous finite elements this coincides with the
   * lowest dimensional object the support point of the degree of freedom lies
   * on. To give an example, for $Q_1$ elements in 3d, every degree of freedom
   * is defined by a shape function that we get by interpolating using support
   * points that lie on the vertices of the cell. The support of these points
   * of course extends to all edges connected to this vertex, as well as the
   * adjacent faces and the cell interior, but we say that logically the
   * degree of freedom is associated with the vertex as this is the lowest-
   * dimensional object it is associated with. Likewise, for $Q_2$ elements in
   * 3d, the degrees of freedom with support points at edge midpoints would
   * yield a value of GeometryPrimitive::line from this function, whereas
   * those on the centers of faces in 3d would return GeometryPrimitive::quad.
   *
   * To make this more formal, the kind of object returned by this function
   * represents the object so that the support of the shape function
   * corresponding to the degree of freedom, (i.e., that part of the domain
   * where the function "lives") is the union of all of the cells sharing this
   * object. To return to the example above, for $Q_2$ in 3d, the shape
   * function with support point at an edge midpoint has support on all cells
   * that share the edge and not only the cells that share the adjacent faces,
   * and consequently the function will return GeometryPrimitive::line.
   *
   * On the other hand, for discontinuous elements of type $DGQ_2$, a degree
   * of freedom associated with an interpolation polynomial that has its
   * support point physically located at a line bounding a cell, but is
   * nonzero only on one cell. Consequently, it is logically associated with
   * the interior of that cell (i.e., with a GeometryPrimitive::quad in 2d and
   * a GeometryPrimitive::hex in 3d).
   *
   * @param[in] cell_dof_index The index of a shape function or degree of
   * freedom. This index must be in the range <code>[0,dofs_per_cell)</code>.
   *
   * @note The integer value of the object returned by this function equals
   * the dimensionality of the object it describes, and can consequently be
   * used in generic programming paradigms. For example, if a degree of
   * freedom is associated with a vertex, then this function returns
   * GeometryPrimitive::vertex, which has a numeric value of zero (the
   * dimensionality of a vertex).
   */
  GeometryPrimitive
  get_associated_geometry_primitive (const unsigned int cell_dof_index) const;


  /**
   * Given the values of a function $f(\mathbf x)$ at the (generalized)
   * support points of the reference cell, this function then computes what
   * the nodal values of the element are, i.e., $\Psi_i[f]$, where $\Psi_i$
   * are the node functionals of the element
   * (see also @ref GlossNodes "Node values or node functionals").
   * The values $\Psi_i[f]$ are then the expansion coefficients
   * for the shape functions of the finite element function that
   * <i>interpolates</i> the given function $f(x)$, i.e.,
   * $ f_h(\mathbf x) = \sum_i \Psi_i[f] \varphi_i(\mathbf x)
   * $ is the finite element interpolant of $f$ with the current element.
   * The operation described here is used, for example, in the
   * FETools::compute_node_matrix() function.
   *
   * In more detail, let us assume that the generalized support points
   * (see
   * @ref GlossGeneralizedSupport "this glossary entry"
   * ) of the current
   * element are $\hat{\mathbf x}_i$ and that the node functionals associated
   * with the current element are $\Psi_i[\cdot]$. Then, the fact that the
   * element is based on generalized support points, implies that if we
   * apply $\Psi_i$ to a (possibly vector-valued) finite element function
   * $\varphi$, the result must have the form
   * $\Psi_i[\varphi] = f_i(\varphi(\hat{\mathbf x}_i))$ -- in other words,
   * the value of the node functional $\Psi_i$ applied to $\varphi$ <i>only</i>
   * depends on the <i>values of $\varphi$ at $\hat{\mathbf x}_i$</i> and not
   * on values anywhere else, or integrals of $\varphi$, or any other kind
   * of information.
   *
   * The exact form of $f_i$ depends on the element. For example, for scalar
   * @ref GlossLagrange "Lagrange elements", we have that in fact
   * $\Psi_i[\varphi] = \varphi(\hat{\mathbf x}_i)$. If you combine multiple
   * scalar Lagrange elements via an FESystem object, then
   * $\Psi_i[\varphi] = \varphi(\hat{\mathbf x}_i)_{c(i)}$ where $c(i)$
   * is the result of the FiniteElement::system_to_component_index()
   * function's return value's first component. In these two cases,
   * $f_i$ is therefore simply the identity (in the scalar case) or a
   * function that selects a particular vector component of its argument.
   * On the other hand, for Raviart-Thomas elements, one would have that
   * $f_i(\mathbf y) = \mathbf y \cdot \mathbf n_i$ where $\mathbf n_i$
   * is the normal vector of the face at which the shape function is
   * defined.
   *
   * Given all of this, what this function does is the following: If you
   * input a list of values of a function $\varphi$ at all generalized
   * support points (where each value is in fact a vector of values with
   * as many components as the element has), then this function returns
   * a vector of values obtained by applying the node functionals to
   * these values. In other words, if you pass in
   * $\{\varphi(\hat{\mathbf x}_i)\}_{i=0}^{N-1}$ then you
   * will get out a vector
   * $\{\Psi[\varphi]\}_{i=0}^{N-1}$ where $N$ equals @p dofs_per_cell.
   *
   * @param[in] support_point_values An array of size @p dofs_per_cell
   *   (which equals the number of points the get_generalized_support_points()
   *   function will return) where each element is a vector with as many entries
   *   as the element has vector components. This array should contain
   *   the values of a function at the generalized support points of the
   *   current element.
   * @param[out] nodal_values An array of size @p dofs_per_cell that contains
   *   the node functionals of the element applied to the given function.
   *
   * @note It is safe to call this function for (transformed) values on the
   * real cell only for elements with trivial MappingType. For all other
   * elements (for example for H(curl), or H(div) conforming elements)
   * vector values have to be transformed to the reference cell first.
   *
   * @note Given what the function is supposed to do, the function clearly
   * can only work for elements that actually implement (generalized) support
   * points. Elements that do not have generalized support points -- e.g.,
   * elements whose nodal functionals evaluate integrals or moments of
   * functions (such as FE_Q_Hierarchical) -- can in general not make
   * sense of the operation that is required for this function. They
   * consequently may not implement it.
   */
  virtual
  void
  convert_generalized_support_point_values_to_dof_values (const std::vector<Vector<double> > &support_point_values,
                                                          std::vector<double>                &nodal_values) const;

  //@}

  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object.
   *
   * This function is made virtual, since finite element objects are usually
   * accessed through pointers to their base class, rather than the class
   * itself.
   */
  virtual std::size_t memory_consumption () const;

  /**
   * Exception
   *
   * @ingroup Exceptions
   */
  DeclException1 (ExcShapeFunctionNotPrimitive,
                  int,
                  << "The shape function with index " << arg1
                  << " is not primitive, i.e. it is vector-valued and "
                  << "has more than one non-zero vector component. This "
                  << "function cannot be called for these shape functions. "
                  << "Maybe you want to use the same function with the "
                  << "_component suffix?");
  /**
   * Exception
   *
   * @ingroup Exceptions
   */
  DeclException0 (ExcFENotPrimitive);
  /**
   * Exception
   *
   * @ingroup Exceptions
   */
  DeclExceptionMsg (ExcUnitShapeValuesDoNotExist,
                    "You are trying to access the values or derivatives of shape functions "
                    "on the reference cell of an element that does not define its shape "
                    "functions through mapping from the reference cell. Consequently, "
                    "you cannot ask for shape function values or derivatives on the "
                    "reference cell.");

  /**
   * Attempt to access support points of a finite element that is not
   * Lagrangian.
   *
   * @ingroup Exceptions
   */
  DeclExceptionMsg (ExcFEHasNoSupportPoints,
                    "You are trying to access the support points of a finite "
                    "element that either has no support points at all, or for "
                    "which the corresponding tables have not been implemented.");

  /**
   * Attempt to access embedding matrices of a finite element that did not
   * implement these matrices.
   *
   * @ingroup Exceptions
   */
  DeclExceptionMsg (ExcEmbeddingVoid,
                    "You are trying to access the matrices that describe how "
                    "to embed a finite element function on one cell into the "
                    "finite element space on one of its children (i.e., the "
                    "'embedding' or 'prolongation' matrices). However, the "
                    "current finite element can either not define this sort of "
                    "operation, or it has not yet been implemented.");

  /**
   * Attempt to access restriction matrices of a finite element that did not
   * implement these matrices.
   *
   * Exception
   * @ingroup Exceptions
   */
  DeclExceptionMsg (ExcProjectionVoid,
                    "You are trying to access the matrices that describe how "
                    "to restrict a finite element function from the children "
                    "of one cell to the finite element space defined on their "
                    "parent (i.e., the 'restriction' or 'projection' matrices). "
                    "However, the current finite element can either not define "
                    "this sort of operation, or it has not yet been "
                    "implemented.");

  /**
   * Exception
   * @ingroup Exceptions
   */
  DeclException2 (ExcWrongInterfaceMatrixSize,
                  int, int,
                  << "The interface matrix has a size of " << arg1
                  << "x" << arg2
                  << ", which is not reasonable for the current element "
                  "in the present dimension.");
  /**
   * Exception
   * @ingroup Exceptions
   */
  DeclException0 (ExcInterpolationNotImplemented);

protected:

  /**
   * Reinit the vectors of restriction and prolongation matrices to the right
   * sizes: For every refinement case, except for
   * RefinementCase::no_refinement, and for every child of that refinement
   * case the space of one restriction and prolongation matrix is allocated,
   * see the documentation of the restriction and prolongation vectors for
   * more detail on the actual vector sizes.
   *
   * @param isotropic_restriction_only only the restriction matrices required
   * for isotropic refinement are reinited to the right size.
   * @param isotropic_prolongation_only only the prolongation matrices
   * required for isotropic refinement are reinited to the right size.
   */
  void reinit_restriction_and_prolongation_matrices(const bool isotropic_restriction_only=false,
                                                    const bool isotropic_prolongation_only=false);

  /**
   * Vector of projection matrices. See get_restriction_matrix() above. The
   * constructor initializes these matrices to zero dimensions, which can be
   * changed by derived classes implementing them.
   *
   * Note, that <code>restriction[refinement_case-1][child]</code> includes
   * the restriction matrix of child <code>child</code> for the RefinementCase
   * <code>refinement_case</code>. Here, we use <code>refinement_case-1</code>
   * instead of <code>refinement_case</code> as for
   * RefinementCase::no_refinement(=0) there are no restriction matrices
   * available.
   */
  std::vector<std::vector<FullMatrix<double> > > restriction;

  /**
   * Vector of embedding matrices. See <tt>get_prolongation_matrix()</tt>
   * above. The constructor initializes these matrices to zero dimensions,
   * which can be changed by derived classes implementing them.
   *
   * Note, that <code>prolongation[refinement_case-1][child]</code> includes
   * the prolongation matrix of child <code>child</code> for the
   * RefinementCase <code>refinement_case</code>. Here, we use
   * <code>refinement_case-1</code> instead of <code>refinement_case</code> as
   * for RefinementCase::no_refinement(=0) there are no prolongation matrices
   * available.
   */
  std::vector<std::vector<FullMatrix<double> > > prolongation;

  /**
   * Specify the constraints which the dofs on the two sides of a cell
   * interface underlie if the line connects two cells of which one is refined
   * once.
   *
   * For further details see the general description of the derived class.
   *
   * This field is obviously useless in one dimension and has there a zero
   * size.
   */
  FullMatrix<double> interface_constraints;

  /**
   * List of support points on the unit cell, in case the finite element has
   * any. The constructor leaves this field empty, derived classes may write
   * in some contents.
   *
   * Finite elements that allow some kind of interpolation operation usually
   * have support points. On the other hand, elements that define their
   * degrees of freedom by, for example, moments on faces, or as derivatives,
   * don't have support points. In that case, this field remains empty.
   */
  std::vector<Point<dim> > unit_support_points;

  /**
   * Same for the faces. See the description of the
   * get_unit_face_support_points() function for a discussion of what
   * contributes a face support point.
   */
  std::vector<Point<dim-1> > unit_face_support_points;

  /**
   * Support points used for interpolation functions of non-Lagrangian
   * elements.
   */
  std::vector<Point<dim> > generalized_support_points;

  /**
   * Face support points used for interpolation functions of non-Lagrangian
   * elements.
   */
  std::vector<Point<dim-1> > generalized_face_support_points;

  /**
   * For faces with non-standard face_orientation in 3D, the dofs on faces
   * (quads) have to be permuted in order to be combined with the correct
   * shape functions. Given a local dof @p index on a quad, return the shift
   * in the local index, if the face has non-standard face_orientation, i.e.
   * <code>old_index + shift = new_index</code>. In 2D and 1D there is no need
   * for permutation so the vector is empty. In 3D it has the size of <code>
   * #dofs_per_quad * 8 </code>, where 8 is the number of orientations, a face
   * can be in (all combinations of the three bool flags face_orientation,
   * face_flip and face_rotation).
   *
   * The standard implementation fills this with zeros, i.e. no permutation at
   * all. Derived finite element classes have to fill this Table with the
   * correct values.
   */
  Table<2,int> adjust_quad_dof_index_for_face_orientation_table;

  /**
   * For lines with non-standard line_orientation in 3D, the dofs on lines
   * have to be permuted in order to be combined with the correct shape
   * functions. Given a local dof @p index on a line, return the shift in the
   * local index, if the line has non-standard line_orientation, i.e.
   * <code>old_index + shift = new_index</code>. In 2D and 1D there is no need
   * for permutation so the vector is empty. In 3D it has the size of
   * #dofs_per_line.
   *
   * The standard implementation fills this with zeros, i.e. no permutation at
   * all. Derived finite element classes have to fill this vector with the
   * correct values.
   */
  std::vector<int> adjust_line_dof_index_for_line_orientation_table;

  /**
   * Store what system_to_component_index() will return.
   */
  std::vector<std::pair<unsigned int, unsigned int> > system_to_component_table;

  /**
   * Map between linear dofs and component dofs on face. This is filled with
   * default values in the constructor, but derived classes will have to
   * overwrite the information if necessary.
   *
   * By component, we mean the vector component, not the base element. The
   * information thus makes only sense if a shape function is non-zero in only
   * one component.
   */
  std::vector<std::pair<unsigned int, unsigned int> > face_system_to_component_table;

  /**
   * For each shape function, store to which base element and which instance
   * of this base element (in case its multiplicity is greater than one) it
   * belongs, and its index within this base element. If the element is not
   * composed of others, then base and instance are always zero, and the index
   * is equal to the number of the shape function. If the element is composed
   * of single instances of other elements (i.e. all with multiplicity one)
   * all of which are scalar, then base values and dof indices within this
   * element are equal to the #system_to_component_table. It differs only in
   * case the element is composed of other elements and at least one of them
   * is vector-valued itself.
   *
   * This array has valid values also in the case of vector-valued (i.e. non-
   * primitive) shape functions, in contrast to the
   * #system_to_component_table.
   */
  std::vector<std::pair<std::pair<unsigned int,unsigned int>,unsigned int> >
  system_to_base_table;

  /**
   * Likewise for the indices on faces.
   */
  std::vector<std::pair<std::pair<unsigned int,unsigned int>,unsigned int> >
  face_system_to_base_table;

  /**
   * For each base element, store the number of blocks generated by the base
   * and the first block in a block vector it will generate.
   */
  BlockIndices base_to_block_indices;

  /**
   * The base element establishing a component.
   *
   * For each component number <tt>c</tt>, the entries have the following
   * meaning: <dl> <dt><tt>table[c].first.first</tt></dt> <dd>Number of the
   * base element for <tt>c</tt>. This is the index you can pass to
   * base_element().</dd> <dt><tt>table[c].first.second</tt></dt>
   * <dd>Component within the base element for <tt>c</tt>. This value is
   * between 0 and the n_components() of this base element.</dd>
   * <dt><tt>table[c].second</tt></dt> <dd>Index of the multiple of the base
   * element that contains <tt>c</tt>. This value is between 0 and the
   * element_multiplicity() of this base element.</dd> </dl>
   *
   * This variable is set to the correct size by the constructor of this
   * class, but needs to be initialized by derived classes, unless its size is
   * one and the only entry is a zero, which is the case for scalar elements.
   * In that case, the initialization by the base class is sufficient.
   *
   * @note This table is filled by FETools::Compositing::build_cell_tables().
   */
  std::vector<std::pair<std::pair<unsigned int, unsigned int>, unsigned int> >
  component_to_base_table;

  /**
   * A flag determining whether restriction matrices are to be concatenated or
   * summed up. See the discussion about restriction matrices in the general
   * class documentation for more information.
   */
  const std::vector<bool> restriction_is_additive_flags;

  /**
   * For each shape function, give a vector of bools (with size equal to the
   * number of vector components which this finite element has) indicating in
   * which component each of these shape functions is non-zero.
   *
   * For primitive elements, there is only one non-zero component.
   */
  const std::vector<ComponentMask> nonzero_components;

  /**
   * This array holds how many values in the respective entry of the
   * #nonzero_components element are non-zero. The array is thus a short-cut
   * to allow faster access to this information than if we had to count the
   * non-zero entries upon each request for this information. The field is
   * initialized in the constructor of this class.
   */
  const std::vector<unsigned int> n_nonzero_components_table;

  /**
   * Store whether all shape functions are primitive. Since finding this out
   * is a very common operation, we cache the result, i.e. compute the value
   * in the constructor for simpler access.
   */
  const bool cached_primitivity;

  /**
   * Return the size of interface constraint matrices. Since this is needed in
   * every derived finite element class when initializing their size, it is
   * placed into this function, to avoid having to recompute the dimension-
   * dependent size of these matrices each time.
   *
   * Note that some elements do not implement the interface constraints for
   * certain polynomial degrees. In this case, this function still returns the
   * size these matrices should have when implemented, but the actual matrices
   * are empty.
   */
  TableIndices<2>
  interface_constraints_size () const;

  /**
   * Given the pattern of nonzero components for each shape function, compute
   * for each entry how many components are non-zero for each shape function.
   * This function is used in the constructor of this class.
   */
  static
  std::vector<unsigned int>
  compute_n_nonzero_components (const std::vector<ComponentMask> &nonzero_components);

  /**
   * Given a set of update flags, compute which other quantities <i>also</i>
   * need to be computed in order to satisfy the request by the given flags.
   * Then return the combination of the original set of flags and those just
   * computed.
   *
   * As an example, if @p update_flags contains update_gradients a finite
   * element class will typically require the computation of the inverse of
   * the Jacobian matrix in order to rotate the gradient of shape functions on
   * the reference cell to the real cell. It would then return not just
   * update_gradients, but also update_covariant_transformation, the flag that
   * makes the mapping class produce the inverse of the Jacobian matrix.
   *
   * An extensive discussion of the interaction between this function and
   * FEValues can be found in the
   * @ref FE_vs_Mapping_vs_FEValues
   * documentation module.
   *
   * @see UpdateFlags
   */
  virtual
  UpdateFlags
  requires_update_flags (const UpdateFlags update_flags) const = 0;

  /**
   * Create an internal data object and return a pointer to it of which the
   * caller of this function then assumes ownership. This object will then be
   * passed to the FiniteElement::fill_fe_values() every time the finite
   * element shape functions and their derivatives are evaluated on a concrete
   * cell. The object created here is therefore used by derived classes as a
   * place for scratch objects that are used in evaluating shape functions, as
   * well as to store information that can be pre-computed once and re-used on
   * every cell (e.g., for evaluating the values and gradients of shape
   * functions on the reference cell, for later re-use when transforming these
   * values to a concrete cell).
   *
   * This function is the first one called in the process of initializing a
   * FEValues object for a given mapping and finite element object. The
   * returned object will later be passed to FiniteElement::fill_fe_values()
   * for a concrete cell, which will itself place its output into an object of
   * type internal::FEValuesImplementation::FiniteElementRelatedData. Since there may be
   * data that can already be computed in its <i>final</i> form on the
   * reference cell, this function also receives a reference to the
   * internal::FEValuesImplementation::FiniteElementRelatedData object as its last argument.
   * This output argument is guaranteed to always be the same one when used
   * with the InternalDataBase object returned by this function. In other
   * words, the subdivision of scratch data and final data in the returned
   * object and the @p output_data object is as follows: If data can be pre-
   * computed on the reference cell in the exact form in which it will later
   * be needed on a concrete cell, then this function should already emplace
   * it in the @p output_data object. An example are the values of shape
   * functions at quadrature points for the usual Lagrange elements which on a
   * concrete cell are identical to the ones on the reference cell. On the
   * other hand, if some data can be pre-computed to make computations on a
   * concrete cell <i>cheaper</i>, then it should be put into the returned
   * object for later re-use in a derive class's implementation of
   * FiniteElement::fill_fe_values(). An example are the gradients of shape
   * functions on the reference cell for Lagrange elements: to compute the
   * gradients of the shape functions on a concrete cell, one has to multiply
   * the gradients on the reference cell by the inverse of the Jacobian of the
   * mapping; consequently, we cannot already compute the gradients on a
   * concrete cell at the time the current function is called, but we can at
   * least pre-compute the gradients on the reference cell, and store it in
   * the object returned.
   *
   * An extensive discussion of the interaction between this function and
   * FEValues can be found in the
   * @ref FE_vs_Mapping_vs_FEValues
   * documentation module. See also the documentation of the InternalDataBase
   * class.
   *
   * @param[in] update_flags A set of UpdateFlags values that describe what
   * kind of information the FEValues object requests the finite element to
   * compute. This set of flags may also include information that the finite
   * element can not compute, e.g., flags that pertain to data produced by the
   * mapping. An implementation of this function needs to set up all data
   * fields in the returned object that are necessary to produce the finite-
   * element related data specified by these flags, and may already pre-
   * compute part of this information as discussed above. Elements may want to
   * store these update flags (or a subset of these flags) in
   * InternalDataBase::update_each so they know at the time when
   * FiniteElement::fill_fe_values() is called what they are supposed to
   * compute
   * @param[in] mapping A reference to the mapping used for computing values
   * and derivatives of shape functions.
   * @param[in] quadrature A reference to the object that describes where the
   * shape functions should be evaluated.
   * @param[out] output_data A reference to the object that FEValues will use
   * in conjunction with the object returned here and where an implementation
   * of FiniteElement::fill_fe_values() will place the requested information.
   * This allows the current function to already pre-compute pieces of
   * information that can be computed on the reference cell, as discussed
   * above. FEValues guarantees that this output object and the object
   * returned by the current function will always be used together.
   * @return A pointer to an object of a type derived from InternalDataBase
   * and that derived classes can use to store scratch data that can be pre-
   * computed, or for scratch arrays that then only need to be allocated once.
   * The calling site assumes ownership of this object and will delete it when
   * it is no longer necessary.
   */
  virtual
  std::unique_ptr<InternalDataBase>
  get_data (const UpdateFlags                                                    update_flags,
            const Mapping<dim,spacedim>                                         &mapping,
            const Quadrature<dim>                                               &quadrature,
            dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim> &output_data) const = 0;

  /**
   * Like get_data(), but return an object that will later be used for
   * evaluating shape function information at quadrature points on faces of
   * cells. The object will then be used in calls to implementations of
   * FiniteElement::fill_fe_face_values(). See the documentation of get_data()
   * for more information.
   *
   * The default implementation of this function converts the face quadrature
   * into a cell quadrature with appropriate quadrature point locations, and
   * with that calls the get_data() function above that has to be implemented
   * in derived classes.
   *
   * @param[in] update_flags A set of UpdateFlags values that describe what
   * kind of information the FEValues object requests the finite element to
   * compute. This set of flags may also include information that the finite
   * element can not compute, e.g., flags that pertain to data produced by the
   * mapping. An implementation of this function needs to set up all data
   * fields in the returned object that are necessary to produce the finite-
   * element related data specified by these flags, and may already pre-
   * compute part of this information as discussed above. Elements may want to
   * store these update flags (or a subset of these flags) in
   * InternalDataBase::update_each so they know at the time when
   * FiniteElement::fill_fe_face_values() is called what they are supposed to
   * compute
   * @param[in] mapping A reference to the mapping used for computing values
   * and derivatives of shape functions.
   * @param[in] quadrature A reference to the object that describes where the
   * shape functions should be evaluated.
   * @param[out] output_data A reference to the object that FEValues will use
   * in conjunction with the object returned here and where an implementation
   * of FiniteElement::fill_fe_face_values() will place the requested
   * information. This allows the current function to already pre-compute
   * pieces of information that can be computed on the reference cell, as
   * discussed above. FEValues guarantees that this output object and the
   * object returned by the current function will always be used together.
   * @return A pointer to an object of a type derived from InternalDataBase
   * and that derived classes can use to store scratch data that can be pre-
   * computed, or for scratch arrays that then only need to be allocated once.
   * The calling site assumes ownership of this object and will delete it when
   * it is no longer necessary.
   */
  virtual
  std::unique_ptr<InternalDataBase>
  get_face_data (const UpdateFlags                                                    update_flags,
                 const Mapping<dim,spacedim>                                         &mapping,
                 const Quadrature<dim-1>                                             &quadrature,
                 dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim> &output_data) const;

  /**
   * Like get_data(), but return an object that will later be used for
   * evaluating shape function information at quadrature points on children of
   * faces of cells. The object will then be used in calls to implementations
   * of FiniteElement::fill_fe_subface_values(). See the documentation of
   * get_data() for more information.
   *
   * The default implementation of this function converts the face quadrature
   * into a cell quadrature with appropriate quadrature point locations, and
   * with that calls the get_data() function above that has to be implemented
   * in derived classes.
   *
   * @param[in] update_flags A set of UpdateFlags values that describe what
   * kind of information the FEValues object requests the finite element to
   * compute. This set of flags may also include information that the finite
   * element can not compute, e.g., flags that pertain to data produced by the
   * mapping. An implementation of this function needs to set up all data
   * fields in the returned object that are necessary to produce the finite-
   * element related data specified by these flags, and may already pre-
   * compute part of this information as discussed above. Elements may want to
   * store these update flags (or a subset of these flags) in
   * InternalDataBase::update_each so they know at the time when
   * FiniteElement::fill_fe_subface_values() is called what they are supposed
   * to compute
   * @param[in] mapping A reference to the mapping used for computing values
   * and derivatives of shape functions.
   * @param[in] quadrature A reference to the object that describes where the
   * shape functions should be evaluated.
   * @param[out] output_data A reference to the object that FEValues will use
   * in conjunction with the object returned here and where an implementation
   * of FiniteElement::fill_fe_subface_values() will place the requested
   * information. This allows the current function to already pre-compute
   * pieces of information that can be computed on the reference cell, as
   * discussed above. FEValues guarantees that this output object and the
   * object returned by the current function will always be used together.
   * @return A pointer to an object of a type derived from InternalDataBase
   * and that derived classes can use to store scratch data that can be pre-
   * computed, or for scratch arrays that then only need to be allocated once.
   * The calling site assumes ownership of this object and will delete it when
   * it is no longer necessary.
   */
  virtual
  std::unique_ptr<InternalDataBase>
  get_subface_data (const UpdateFlags                                                    update_flags,
                    const Mapping<dim,spacedim>                                         &mapping,
                    const Quadrature<dim-1>                                             &quadrature,
                    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim> &output_data) const;

  /**
   * Compute information about the shape functions on the cell denoted by the
   * first argument. Derived classes will have to implement this function
   * based on the kind of element they represent. It is called by
   * FEValues::reinit().
   *
   * Conceptually, this function evaluates shape functions and their
   * derivatives at the quadrature points represented by the mapped locations
   * of those described by the quadrature argument to this function. In many
   * cases, computing derivatives of shape functions (and in some cases also
   * computing values of shape functions) requires making use of the mapping
   * from the reference to the real cell; this information can either be taken
   * from the @p mapping_data object that has been filled for the current cell
   * before this function is called, or by calling the member functions of a
   * Mapping object with the @p mapping_internal object that also corresponds
   * to the current cell.
   *
   * The information computed by this function is used to fill the various
   * member variables of the output argument of this function. Which of the
   * member variables of that structure should be filled is determined by the
   * update flags stored in the FiniteElement::InternalDataBase::update_each
   * field of the object passed to this function. These flags are typically
   * set by FiniteElement::get_data(), FiniteElement::get_face_date() and
   * FiniteElement::get_subface_data() (or, more specifically, implementations
   * of these functions in derived classes).
   *
   * An extensive discussion of the interaction between this function and
   * FEValues can be found in the
   * @ref FE_vs_Mapping_vs_FEValues
   * documentation module.
   *
   * @param[in] cell The cell of the triangulation for which this function is
   * to compute a mapping from the reference cell to.
   * @param[in] cell_similarity Whether or not the cell given as first
   * argument is simply a translation, rotation, etc of the cell for which
   * this function was called the most recent time. This information is
   * computed simply by matching the vertices (as stored by the Triangulation)
   * between the previous and the current cell. The value passed here may be
   * modified by implementations of this function and should then be returned
   * (see the discussion of the return value of this function).
   * @param[in] quadrature A reference to the quadrature formula in use for
   * the current evaluation. This quadrature object is the same as the one
   * used when creating the @p internal_data object. The current object is
   * then responsible for evaluating shape functions at the mapped locations
   * of the quadrature points represented by this object.
   * @param[in] mapping A reference to the mapping object used to map from the
   * reference cell to the current cell. This object was used to compute the
   * information in the @p mapping_data object before the current function was
   * called. It is also the mapping object that created the @p
   * mapping_internal object via Mapping::get_data(). You will need the
   * reference to this mapping object most often to call Mapping::transform()
   * to transform gradients and higher derivatives from the reference to the
   * current cell.
   * @param[in] mapping_internal An object specific to the mapping object.
   * What the mapping chooses to store in there is of no relevance to the
   * current function, but you may have to pass a reference to this object to
   * certain functions of the Mapping class (e.g., Mapping::transform()) if
   * you need to call them from the current function.
   * @param[in] mapping_data The output object into which the
   * Mapping::fill_fe_values() function wrote the mapping information
   * corresponding to the current cell. This includes, for example, Jacobians
   * of the mapping that may be of relevance to the current function, as well
   * as other information that FEValues::reinit() requested from the mapping.
   * @param[in] fe_internal A reference to an object previously created by
   * get_data() and that may be used to store information the mapping can
   * compute once on the reference cell. See the documentation of the
   * FiniteElement::InternalDataBase class for an extensive description of the
   * purpose of these objects.
   * @param[out] output_data A reference to an object whose member variables
   * should be computed. Not all of the members of this argument need to be
   * filled; which ones need to be filled is determined by the update flags
   * stored inside the @p fe_internal object.
   *
   * @note FEValues ensures that this function is always called with the same
   * pair of @p fe_internal and @p output_data objects. In other words, if an
   * implementation of this function knows that it has written a piece of data
   * into the output argument in a previous call, then there is no need to
   * copy it there again in a later call if the implementation knows that this
   * is the same value.
   */
  virtual
  void
  fill_fe_values (const typename Triangulation<dim,spacedim>::cell_iterator           &cell,
                  const CellSimilarity::Similarity                                     cell_similarity,
                  const Quadrature<dim>                                               &quadrature,
                  const Mapping<dim,spacedim>                                         &mapping,
                  const typename Mapping<dim,spacedim>::InternalDataBase              &mapping_internal,
                  const dealii::internal::FEValuesImplementation::MappingRelatedData<dim, spacedim> &mapping_data,
                  const InternalDataBase                                              &fe_internal,
                  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim> &output_data) const = 0;

  /**
   * This function is the equivalent to FiniteElement::fill_fe_values(), but
   * for faces of cells. See there for an extensive discussion of its purpose.
   * It is called by FEFaceValues::reinit().
   *
   * @param[in] cell The cell of the triangulation for which this function is
   * to compute a mapping from the reference cell to.
   * @param[in] face_no The number of the face we are currently considering,
   * indexed among the faces of the cell specified by the previous argument.
   * @param[in] quadrature A reference to the quadrature formula in use for
   * the current evaluation. This quadrature object is the same as the one
   * used when creating the @p internal_data object. The current object is
   * then responsible for evaluating shape functions at the mapped locations
   * of the quadrature points represented by this object.
   * @param[in] mapping A reference to the mapping object used to map from the
   * reference cell to the current cell. This object was used to compute the
   * information in the @p mapping_data object before the current function was
   * called. It is also the mapping object that created the @p
   * mapping_internal object via Mapping::get_data(). You will need the
   * reference to this mapping object most often to call Mapping::transform()
   * to transform gradients and higher derivatives from the reference to the
   * current cell.
   * @param[in] mapping_internal An object specific to the mapping object.
   * What the mapping chooses to store in there is of no relevance to the
   * current function, but you may have to pass a reference to this object to
   * certain functions of the Mapping class (e.g., Mapping::transform()) if
   * you need to call them from the current function.
   * @param[in] mapping_data The output object into which the
   * Mapping::fill_fe_values() function wrote the mapping information
   * corresponding to the current cell. This includes, for example, Jacobians
   * of the mapping that may be of relevance to the current function, as well
   * as other information that FEValues::reinit() requested from the mapping.
   * @param[in] fe_internal A reference to an object previously created by
   * get_data() and that may be used to store information the mapping can
   * compute once on the reference cell. See the documentation of the
   * FiniteElement::InternalDataBase class for an extensive description of the
   * purpose of these objects.
   * @param[out] output_data A reference to an object whose member variables
   * should be computed. Not all of the members of this argument need to be
   * filled; which ones need to be filled is determined by the update flags
   * stored inside the @p fe_internal object.
   */
  virtual
  void
  fill_fe_face_values (const typename Triangulation<dim,spacedim>::cell_iterator           &cell,
                       const unsigned int                                                   face_no,
                       const Quadrature<dim-1>                                             &quadrature,
                       const Mapping<dim,spacedim>                                         &mapping,
                       const typename Mapping<dim,spacedim>::InternalDataBase              &mapping_internal,
                       const dealii::internal::FEValuesImplementation::MappingRelatedData<dim, spacedim> &mapping_data,
                       const InternalDataBase                                              &fe_internal,
                       dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim> &output_data) const = 0;

  /**
   * This function is the equivalent to FiniteElement::fill_fe_values(), but
   * for the children of faces of cells. See there for an extensive discussion
   * of its purpose. It is called by FESubfaceValues::reinit().
   *
   * @param[in] cell The cell of the triangulation for which this function is
   * to compute a mapping from the reference cell to.
   * @param[in] face_no The number of the face we are currently considering,
   * indexed among the faces of the cell specified by the previous argument.
   * @param[in] sub_no The number of the subface, i.e., the number of the
   * child of a face, that we are currently considering, indexed among the
   * children of the face specified by the previous argument.
   * @param[in] quadrature A reference to the quadrature formula in use for
   * the current evaluation. This quadrature object is the same as the one
   * used when creating the @p internal_data object. The current object is
   * then responsible for evaluating shape functions at the mapped locations
   * of the quadrature points represented by this object.
   * @param[in] mapping A reference to the mapping object used to map from the
   * reference cell to the current cell. This object was used to compute the
   * information in the @p mapping_data object before the current function was
   * called. It is also the mapping object that created the @p
   * mapping_internal object via Mapping::get_data(). You will need the
   * reference to this mapping object most often to call Mapping::transform()
   * to transform gradients and higher derivatives from the reference to the
   * current cell.
   * @param[in] mapping_internal An object specific to the mapping object.
   * What the mapping chooses to store in there is of no relevance to the
   * current function, but you may have to pass a reference to this object to
   * certain functions of the Mapping class (e.g., Mapping::transform()) if
   * you need to call them from the current function.
   * @param[in] mapping_data The output object into which the
   * Mapping::fill_fe_values() function wrote the mapping information
   * corresponding to the current cell. This includes, for example, Jacobians
   * of the mapping that may be of relevance to the current function, as well
   * as other information that FEValues::reinit() requested from the mapping.
   * @param[in] fe_internal A reference to an object previously created by
   * get_data() and that may be used to store information the mapping can
   * compute once on the reference cell. See the documentation of the
   * FiniteElement::InternalDataBase class for an extensive description of the
   * purpose of these objects.
   * @param[out] output_data A reference to an object whose member variables
   * should be computed. Not all of the members of this argument need to be
   * filled; which ones need to be filled is determined by the update flags
   * stored inside the @p fe_internal object.
   */
  virtual
  void
  fill_fe_subface_values (const typename Triangulation<dim,spacedim>::cell_iterator           &cell,
                          const unsigned int                                                   face_no,
                          const unsigned int                                                   sub_no,
                          const Quadrature<dim-1>                                             &quadrature,
                          const Mapping<dim,spacedim>                                         &mapping,
                          const typename Mapping<dim,spacedim>::InternalDataBase              &mapping_internal,
                          const dealii::internal::FEValuesImplementation::MappingRelatedData<dim, spacedim> &mapping_data,
                          const InternalDataBase                                              &fe_internal,
                          dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim> &output_data) const = 0;

  friend class InternalDataBase;
  friend class FEValuesBase<dim,spacedim>;
  friend class FEValues<dim,spacedim>;
  friend class FEFaceValues<dim,spacedim>;
  friend class FESubfaceValues<dim,spacedim>;
  friend class FESystem<dim,spacedim>;

  // explicitly check for sensible template arguments, but not on windows
  // because MSVC creates bogus warnings during normal compilation
#ifndef DEAL_II_MSVC
  static_assert (dim<=spacedim,
                 "The dimension <dim> of a FiniteElement must be less than or "
                 "equal to the space dimension <spacedim> in which it lives.");
#endif

};


//----------------------------------------------------------------------//


template <int dim, int spacedim>
inline
const FiniteElement<dim,spacedim> &
FiniteElement<dim,spacedim>::operator[] (const unsigned int fe_index) const
{
  (void)fe_index;
  Assert (fe_index == 0,
          ExcMessage ("A fe_index of zero is the only index allowed here"));
  return *this;
}



template <int dim, int spacedim>
inline
std::pair<unsigned int,unsigned int>
FiniteElement<dim,spacedim>::system_to_component_index (const unsigned int index) const
{
  Assert (index < system_to_component_table.size(),
          ExcIndexRange(index, 0, system_to_component_table.size()));
  Assert (is_primitive (index),
          ( typename FiniteElement<dim,spacedim>::ExcShapeFunctionNotPrimitive(index)) );
  return system_to_component_table[index];
}



template <int dim, int spacedim>
inline
unsigned int
FiniteElement<dim,spacedim>::n_base_elements () const
{
  return base_to_block_indices.size();
}



template <int dim, int spacedim>
inline
unsigned int
FiniteElement<dim,spacedim>::element_multiplicity (const unsigned int index) const
{
  return static_cast<unsigned int>(base_to_block_indices.block_size(index));
}



template <int dim, int spacedim>
inline
unsigned int
FiniteElement<dim,spacedim>::component_to_system_index (const unsigned int component,
                                                        const unsigned int index) const
{
  AssertIndexRange(component, this->n_components());
  const std::vector<std::pair<unsigned int, unsigned int> >::const_iterator
  it = std::find(system_to_component_table.begin(), system_to_component_table.end(),
                 std::pair<unsigned int, unsigned int>(component, index));

  Assert(it != system_to_component_table.end(),
         ExcMessage ("You are asking for the number of the shape function "
                     "within a system element that corresponds to vector "
                     "component " + Utilities::int_to_string(component) + " and within this to "
                     "index " + Utilities::int_to_string(index) + ". But no such "
                     "shape function exists."));
  return std::distance(system_to_component_table.begin(), it);
}



template <int dim, int spacedim>
inline
std::pair<unsigned int,unsigned int>
FiniteElement<dim,spacedim>::face_system_to_component_index (const unsigned int index) const
{
  Assert(index < face_system_to_component_table.size(),
         ExcIndexRange(index, 0, face_system_to_component_table.size()));

  // in debug mode, check whether the
  // function is primitive, since
  // otherwise the result may have no
  // meaning
  //
  // since the primitivity tables are
  // all geared towards cell dof
  // indices, rather than face dof
  // indices, we have to work a
  // little bit...
  //
  // in 1d, the face index is equal
  // to the cell index
  Assert (is_primitive(this->face_to_cell_index(index, 0)),
          (typename FiniteElement<dim,spacedim>::ExcShapeFunctionNotPrimitive(index)) );

  return face_system_to_component_table[index];
}




template <int dim, int spacedim>
inline
std::pair<std::pair<unsigned int,unsigned int>,unsigned int>
FiniteElement<dim,spacedim>::system_to_base_index (const unsigned int index) const
{
  Assert (index < system_to_base_table.size(),
          ExcIndexRange(index, 0, system_to_base_table.size()));
  return system_to_base_table[index];
}




template <int dim, int spacedim>
inline
std::pair<std::pair<unsigned int,unsigned int>,unsigned int>
FiniteElement<dim,spacedim>::face_system_to_base_index (const unsigned int index) const
{
  Assert(index < face_system_to_base_table.size(),
         ExcIndexRange(index, 0, face_system_to_base_table.size()));
  return face_system_to_base_table[index];
}



template <int dim, int spacedim>
inline
types::global_dof_index
FiniteElement<dim,spacedim>::first_block_of_base (const unsigned int index) const
{
  return base_to_block_indices.block_start(index);
}



template <int dim, int spacedim>
inline
std::pair<unsigned int,unsigned int>
FiniteElement<dim,spacedim>::component_to_base_index (const unsigned int index) const
{
  Assert(index < component_to_base_table.size(),
         ExcIndexRange(index, 0, component_to_base_table.size()));

  return component_to_base_table[index].first;
}



template <int dim, int spacedim>
inline
std::pair<unsigned int,unsigned int>
FiniteElement<dim,spacedim>::block_to_base_index (const unsigned int index) const
{
  return base_to_block_indices.global_to_local(index);
}



template <int dim, int spacedim>
inline
std::pair<unsigned int,types::global_dof_index>
FiniteElement<dim,spacedim>::system_to_block_index (const unsigned int index) const
{
  Assert (index < this->dofs_per_cell,
          ExcIndexRange(index, 0, this->dofs_per_cell));
  // The block is computed simply as
  // first block of this base plus
  // the index within the base blocks
  return std::pair<unsigned int, types::global_dof_index>(
           first_block_of_base(system_to_base_table[index].first.first)
           + system_to_base_table[index].first.second,
           system_to_base_table[index].second);
}



template <int dim, int spacedim>
inline
bool
FiniteElement<dim,spacedim>::restriction_is_additive (const unsigned int index) const
{
  Assert(index < this->dofs_per_cell,
         ExcIndexRange(index, 0, this->dofs_per_cell));
  return restriction_is_additive_flags[index];
}



template <int dim, int spacedim>
inline
const ComponentMask &
FiniteElement<dim,spacedim>::get_nonzero_components (const unsigned int i) const
{
  Assert (i < this->dofs_per_cell, ExcIndexRange (i, 0, this->dofs_per_cell));
  return nonzero_components[i];
}



template <int dim, int spacedim>
inline
unsigned int
FiniteElement<dim,spacedim>::n_nonzero_components (const unsigned int i) const
{
  Assert (i < this->dofs_per_cell, ExcIndexRange (i, 0, this->dofs_per_cell));
  return n_nonzero_components_table[i];
}



template <int dim, int spacedim>
inline
bool
FiniteElement<dim,spacedim>::is_primitive () const
{
  return cached_primitivity;
}



template <int dim, int spacedim>
inline
bool
FiniteElement<dim,spacedim>::is_primitive (const unsigned int i) const
{
  Assert (i < this->dofs_per_cell, ExcIndexRange (i, 0, this->dofs_per_cell));

  // return primitivity of a shape
  // function by checking whether it
  // has more than one non-zero
  // component or not. we could cache
  // this value in an array of bools,
  // but accessing a bit-vector (as
  // std::vector<bool> is) is
  // probably more expensive than
  // just comparing against 1
  //
  // for good measure, short circuit the test
  // if the entire FE is primitive
  return (is_primitive() ||
          (n_nonzero_components_table[i] == 1));
}



template <int dim, int spacedim>
inline
GeometryPrimitive
FiniteElement<dim,spacedim>::get_associated_geometry_primitive (const unsigned int cell_dof_index) const
{
  Assert (cell_dof_index < this->dofs_per_cell,
          ExcIndexRange (cell_dof_index, 0, this->dofs_per_cell));

  // just go through the usual cases, taking into account how DoFs
  // are enumerated on the reference cell
  if (cell_dof_index < this->first_line_index)
    return GeometryPrimitive::vertex;
  else if (cell_dof_index < this->first_quad_index)
    return GeometryPrimitive::line;
  else if (cell_dof_index < this->first_hex_index)
    return GeometryPrimitive::quad;
  else
    return GeometryPrimitive::hex;
}



DEAL_II_NAMESPACE_CLOSE

#endif
