//----------------------------  fe.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  fe.h  ---------------------------
#ifndef __deal2__fe_h
#define __deal2__fe_h


#include <base/exceptions.h>
#include <base/subscriptor.h>
#include <base/point.h>
#include <base/tensor.h>
#include <dofs/dof_handler.h>
#include <grid/geometry_info.h>
#include <lac/full_matrix.h>


/**
 * Dimension independent data for finite elements. See the #FiniteElementBase#
 * class for more information.
 */
template<int dim>
class FiniteElementData
{
  public:
				     /**
				      * Number of degrees of freedom on
				      * a vertex.
				      */
    const unsigned int dofs_per_vertex;

				     /** Number of degrees of freedom
				      *  on a line.
				      */
    const unsigned int dofs_per_line;

				     /** Number of degrees of freedom
				      *  on a quadrilateral.
				      */
    const unsigned int dofs_per_quad;

				     /** Number of degrees of freedom
				      *  on a hexahedron.
				      */
    const unsigned int dofs_per_hex;

				     /**
				      * First index of dof on a line.
				      */
    const unsigned int first_line_index;
    
				     /**
				      * First index of dof on a quad.
				      */
    const unsigned int first_quad_index;
    
				     /**
				      * First index of dof on a hexahedron.
				      */
    const unsigned int first_hex_index;
    
				     /**
				      * First index of dof on a line for face data.
				      */
    const unsigned int first_face_line_index;
    
				     /**
				      * First index of dof on a quad for face data.
				      */
    const unsigned int first_face_quad_index;

				     /**
				      * Number of degrees of freedom on a
				      * face. This information is
				      * redundant to some fields in the
				      * derived classes but makes
				      * writing dimension independant
				      * programs easier.
				      */
    const unsigned int dofs_per_face;
    
				     /**
				      * Total number of degrees of freedom
				      * on a cell. This information is
				      * redundant to some fields in the
				      * derived classes but makes
				      * writing dimension independant
				      * programs easier.
				      */
    const unsigned int dofs_per_cell;

				     /**
				      * Number of basis functions used for the
				      * transformation from unit cell to real
				      * cell. For a linear mapping, this number
				      * equals the number of vertices.
				      */
    const unsigned int transform_functions;


				     /**
				      * Number of components and dimension of
				      * the image space.
				      */
    const unsigned int components;

    				     /**
				      * Default constructor. Constructs
				      * an element
				      * which is not so useful. Checking
				      * #dofs_per_cell# is therefore a good way to
				      * check if something went wrong. 
				      */
    FiniteElementData ();

				     /**
				      * Constructor for a 1-dimensional
				      * object.
				      */
    FiniteElementData (const unsigned int dofs_per_vertex,
		       const unsigned int dofs_per_line,
		       const unsigned int n_transform_functions,
		       const unsigned int n_components);

				     /**
				      * Constructor for a 2-dimensional
				      * object.
				      */
    FiniteElementData (const unsigned int dofs_per_vertex,
		       const unsigned int dofs_per_line,
		       const unsigned int dofs_per_quad,
		       const unsigned int n_transform_functions,
		       const unsigned int n_components);

				     /**
				      * Constructor for a 3-dimensional
				      * object.
				      */
    FiniteElementData (const unsigned int dofs_per_vertex,
		       const unsigned int dofs_per_line,
		       const unsigned int dofs_per_quad,
		       const unsigned int dofs_per_hex,
		       const unsigned int n_transform_functions,
		       const unsigned int n_components);

				     /**
				      * Declare this destructor virtual in
				      * order to make the respective destructors
				      * in derived classes virtual as well.
				      */
    virtual ~FiniteElementData ();

				     /**
				      * Return the #dofs_per_vertex#.
				      */
    unsigned int n_dofs_per_vertex () const;

				     /**
				      * Return the #dofs_per_line#.
				      */
    unsigned int n_dofs_per_line () const;

    				     /**
				      * Return the #dofs_per_quad#.
				      */
    unsigned int n_dofs_per_quad () const;

    				     /**
				      * Return the #dofs_per_hex#.
				      */
    unsigned int n_dofs_per_hex () const;

    				     /**
				      * Return the #dofs_per_face#.
				      */
    unsigned int n_dofs_per_face () const;

    				     /**
				      * Return the #dofs_per_cell#.
				      */
    unsigned int n_dofs_per_cell () const;

    				     /**
				      * Return the #components#.
				      */
    unsigned int n_components () const;

    				     /**
				      * Return the #transform_functions#.
				      */
    unsigned int n_transform_functions () const;


				     /**
				      * Comparison operator. It is not clear to
				      * me (WB) why we have to declare and implement
				      * this one explicitely.
				      */
    bool operator == (const FiniteElementData<dim> &) const;

				     /**
				      * Exception
				      */
    DeclException2 (ExcDimensionMismatch, int, int,
		    << "used " << arg1 << "-d constructor for " << arg2 << "-d object");
};


/**
 * Base class for finite elements in arbitrary dimensions. This class provides
 * several fields which describe a specific finite element and which are filled
 * by derived classes. It more or less only offers the fields and access
 * functions which makes it possible to copy finite elements without knowledge
 * of the actual type (linear, quadratic, etc).
 *
 * The implementation of this base class is split into two parts: those fields
 * which are not common to all dimensions (#dofs_per_quad# for example are only
 * useful for #dim>=2#) are put into the #FiniteElementData<dim># class which
 * is explicitely specialized for all used dimensions, while those fields which
 * may be formulated in a dimension-independent way are put into the present
 * class.
 *
 * The different matrices are initialized with the correct size, such that in
 * the derived (concrete) finite element classes, their entries must only be
 * filled in; no resizing is needed.
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
class FiniteElementBase : public Subscriptor,
			  public FiniteElementData<dim>
{
  public:
				     /**
				      * Construct an object of this type.
				      * You have to set the
				      * matrices explicitely after calling
				      * this base class' constructor.
				      */
    FiniteElementBase (const FiniteElementData<dim> &fe_data,
		       const vector<bool> &restriction_is_additive_flags);
    
				     /**
				      * Return a readonly reference to the
				      * matrix which describes the transfer of a
				      * child with the given number to the
				      * mother cell. See the #restriction# array
				      * for more information.
				      */
    const FullMatrix<double> & restrict (const unsigned int child) const;

				     /**
				      * Return a readonly reference to the
				      * matrix which describes the transfer of a
				      * mother cell to the child with the given
				      * number.
				      */
    const FullMatrix<double> & prolongate (const unsigned int child) const;

				     /**
				      * Return a readonly reference to the
				      * matrix which describes the constraints
				      * at the interface between a refined and
				      * an unrefined cell.
				      *
				      * The matrix is obviously empty in only
				      * one space dimension, since there are no
				      * constraints then.
				      */
    const FullMatrix<double> & constraints () const;

				     /**
				      * Comparison operator. We also check for
				      * equality of the constraint matrix,
				      * which is quite an expensive operation.
				      * Do therefore
				      * use this function with care, if possible
				      * only for debugging purposes.
				      *
				      * Since this function is not that important,
				      * we avoid an implementational question
				      * about comparing arrays and do not compare
				      * the matrix arrays #restriction# and
				      * #prolongation#.
				      */
    bool operator == (const FiniteElementBase<dim> &) const;

				     /**
				      * Compute system index from components.
				      */
    unsigned int component_to_system_index (unsigned int component,
					    unsigned int component_index) const;
  
				     /**
				      * Compute component and index from
				      * system index.
				      *
				      * Return value contains first
				      * component and second index in
				      * component.
				      */
    pair<unsigned int,unsigned int> system_to_component_index (unsigned int index) const; 
    
				     /**
				      * Compute system index from components on a face.
				      */
    unsigned int face_component_to_system_index (unsigned int component,
						 unsigned int component_index) const;
  
				     /**
				      * Compute component and index from system
				      * index for a face.
				      *
				      * Return value contains first
				      * component and second index in
				      * component.
				      */
    pair<unsigned int,unsigned int> face_system_to_component_index (unsigned int index) const;
    
 				     /**
				      * The base element establishing a
				      * component.
				      *
				      * This table converts a
				      * component number to the
				      * #base_element# number. While
				      * component information contains
				      * multiplicity of base elements,
				      * the result allows access to
				      * shape functions of the base
				      * element.
				      */
    unsigned int component_to_base(unsigned int index) const;

				     /**
				      * Access the #restriction_is_additive_flag#
				      * field. See there for more information on 
				      * its contents.
				      */
    bool restriction_is_additive (const unsigned int component) const;

				     /**
				      * Exception
				      */
    DeclException2 (ExcWrongFieldDimension,
		    int, int,
		    << "The field has not the assumed dimension " << arg2
		    << ", but has " << arg1 << " elements.");
    DeclException2 (ExcWrongInterfaceMatrixSize,
		    int, int,
		    << "The interface matrix has a size of " << arg1
		    << "x" << arg2
		    << ", which is not reasonable in the present dimension.");
    
  protected:
				     /**
				      * Have #N=2^dim# matrices keeping the
				      * restriction constants for the transfer
				      * of the #i#th child to the mother cell.
				      * The numbering conventions for the
				      * degree of freedom indices are descibed
				      * in the derived classes.
				      * In this matrix, the row indices belong
				      * to the destination cell, i.e. the
				      * unrefined one, while the column indices
				      * are for the refined cell's degrees of
				      * freedom. The application of this matrix
				      * is therefore usually its being
				      * multiplied by the vector of nodal values
				      * on the child.
				      *
				      * In essence, using the matrices from the
				      * children to the mother cell amounts to
				      * computing the interpolation of the
				      * function on the refined to the coarse
				      * mesh. To get the vector of nodal values
				      * of the interpolant on the mother cell,
				      * you have to multiply the nodal value
				      * vectors of each of the child cell with
				      * the respective restriction matrix and
				      * clobber these contributions together.
				      * However, you must take care not to
				      * #add# these together, since nodes which
				      * belong to more than one child would then
				      * be counted more than once; rather, you
				      * have to overwrite the nonzero
				      * contributions of each child into the
				      * nodal value vector of the mother cell.
				      *
				      * While we could avoid this and rather add
				      * up the contributions of each child for
				      * nodes that are interior of the mother
				      * cell, we can't for nodes on the boundary
				      * of the mother cell. The reason for this
				      * is that we know how many children may
				      * contribute to the interpolated nodal
				      * value of an interior degree of freedom.
				      * However, we don't know for dofs on the
				      * boundary, for which we only know how many
				      * children from each side of the face
				      * contribute, but we would have to look out
				      * of the cell to know how many neighbors
				      * there are and then, still, we would have
				      * to have two different interpolation
				      * routines for local interpolation and for
				      * the contribution of a cell to a global
				      * interpolation if we wanted to compute
				      * that by adding up local contributions.
				      *
				      * Because of this problem, we chose to
				      * write rather than add the contributions
				      * of each cell to the interpolation onto
				      * the mother cell. However, there now is
				      * another problem which appears when using
				      * discontinuous finite elements. The
				      * process of 'writing' assumed that we
				      * get the same result for each degree of
				      * freedom from each of the children, such
				      * that 'over'writing would not destroy
				      * information; when using discontinuous
				      * elements, this assumption is violated.
				      *
				      * One way to get the whole thing working
				      * nonetheless would be to use a flag which
				      * tells us when the 'write' and when to
				      * 'add'. Adding is possible when we have
				      * only interior degrees of freedom, i.e.
				      * dofs that are not shared between cells.
				      * This is certainly the way to go for the
				      * DG(r) elements. However, this scheme
				      * does not work for discontinuous elements
				      * with degrees of freedom on the faces,
				      * such as the rotated bilinear
				      * (Rannacher-Turek) element or elements
				      * like the Crouzeix-Raviart one. The latter
				      * have degrees of freedom on the faces,
				      * e.g. mean values on a face, and these
				      * values are the same from both sides of
				      * the face, but the solutions are
				      * discontinuous there nevertheless, such
				      * that interpolation is not possible here.
				      * 
				      * At least for the first case, the flags
				      * #restriction_is_additive# was introduced,
				      * see there for more information. For the
				      * latter case, NO SOLUTION HAS BEEN MADE
				      * UP YET.
				      *
				      * To compute the interpolation of a
				      * finite element field to a cell, you
				      * may use the #get_interpolated_dof_values#
				      * function of the #DoFCellAccessor# class.
				      * See there for more information.
				      *
				      * Upon assembling the transfer matrix
				      * between cells using this matrix array,
				      * zero elements in the restriction
				      * matrix are discarded and will not fill
				      * up the transfer matrix.
				      */
#if !((__GNUC__==2) && (__GNUC_MINOR__==95))
    FullMatrix<double> restriction[GeometryInfo<dim>::children_per_cell];
#else
    FullMatrix<double> restriction[1 << dim];
#endif

    				     /**
				      * Have #N=2^dim# matrices keeping the
				      * prolongation constants for the transfer
				      * of the mother cell to the #i#th child.
				      * The numbering conventions for the
				      * degree of freedom indices are descibed
				      * in the derived classes.
				      * In this matrix, the row indices belong
				      * to the destination cell, i.e. the
				      * refined one, while the column indices
				      * are for the unrefined cell's degrees of
				      * freedom. Thus, if #u0# is the vector
				      * of values of degrees of freedom on the
				      * coarse cell, #prolongation[i]*u0#
				      * yields the vector of values of the
				      * degrees of freedom on the #i#th child
				      * cell.
				      *
				      * On the other hand, for finite elements
				      * with embedded spaces, the basis function
				      * phi0[i] on the coarse grid can be
				      * expressed by
				      * $\sum_c \sum_j p^c_{ji} phi1[j]$
				      * where the sum over c runs over the child
				      * cells and phi1[j] is the jth basis
				      * function on the cth child cell. Note
				      * that we need here the transpose of the
				      * matrix $p^c$ ($p^c$ is returned by this
				      * function with parameter c).
				      *
				      * Upon assembling the transfer matrix
				      * between cells using this matrix array,
				      * zero elements in the prolongation
				      * matrix are discarded and will not fill
				      * up the transfer matrix.
				      */
#if ! ((__GNUC__==2) && (__GNUC_MINOR__==95))
    FullMatrix<double> prolongation[GeometryInfo<dim>::children_per_cell];
#else
    FullMatrix<double> prolongation[1 << dim];
#endif

    				     /**
				      * Specify the constraints which
				      * the dofs on the two sides of a
				      * cell interface underly if the
				      * line connects two cells of
				      * which one is refined once.
				      *
				      * For further details see the
				      * general description of the
				      * derived class.
				      *
				      * This field is obviously
				      * useless in one space dimension
				      * and has there a zero size.
				      */
    FullMatrix<double> interface_constraints;

				     /**
				      * Map between linear dofs and
				      * component dofs.
				      */
    vector< pair<unsigned int, unsigned int> > system_to_component_table;

				     /**
				      * Map between linear dofs and
				      * component dofs on face.
				      */
    vector< pair<unsigned int, unsigned int> > face_system_to_component_table;

				     /**
				      * Map between component and
				      * linear dofs.
				      */
    vector< vector<unsigned int> > component_to_system_table;

				     /**
				      * Map between component and
				      * linear dofs on a face.
				      */
    vector< vector<unsigned int> > face_component_to_system_table;
    
				     /**
				      * The base element establishing
				      * a component.
				      *
				      * This table converts a
				      * component number to the
				      * #base_element# number. While
				      * component information contains
				      * multiplicity of base elements,
				      * the result allows access to
				      * shape functions of the base
				      * element.
				      */
    vector<unsigned int> component_to_base_table;

				     /**
				      * This flag determines how the
				      * restriction of data from child
				      * cells to its mother is to be
				      * done. In this, it also
				      * determines in which way the
				      * restriction matrices of the
				      * derived class are to be used.
				      *
				      * For most elements, the mode is
				      * the following. Consider a 1d
				      * linear element, with two
				      * children and nodal values 1
				      * and 2 on the first child, and
				      * 2 and 4 on the second
				      * child. The restriction to the
				      * mother child then yields the
				      * values 1 and four, i.e. the
				      * values on the mother cell can
				      * be obtained by pointwise
				      * interpolation, where for each
				      * nodal value on the mother
				      * child one point on exactly one
				      * child exists.  However,
				      * already on the quadratic
				      * element, the midpoint on the
				      * mother element can be obtained
				      * from any of the two children,
				      * which however would both yield
				      * the same value due to
				      * continuity. What we do in
				      * practice is to compute them
				      * from both sides and set them,
				      * rather than add them up.  This
				      * makes some things much
				      * easier. In practice, if a
				      * degree of freedom on one of
				      * the child cells yields a
				      * nonzero contribution to one of
				      * the degrees of freedom on the
				      * mother cell, we overwrite the
				      * value on the mother cell. This
				      * way, when setting up the
				      * restriction matrices, we do
				      * not have to track which child
				      * is responsible for setting a
				      * given value on the mother
				      * cell. We call this the
				      * non-additive mode.
				      *
				      * The other possibility would be
				      * to add up the contributions
				      * from the different
				      * children. This would mean that
				      * both of the inner endpoint of
				      * the quadratic child elements
				      * above would have a weight of
				      * 1/2 with respect to the
				      * midpoint value on the mother
				      * cell. However, this also means
				      * that we have to first compute
				      * the restriction to the mother
				      * cell by addition from the
				      * child cells, and afterwards
				      * set them to the global
				      * vector. The same process,
				      * adding up the local
				      * contributions to the global
				      * vector is not possible since
				      * we do not know how many coarse
				      * cells contribute to nodes on
				      * the boundary.
				      *
				      * In contrast to the
				      * non-additive mode described
				      * above, which is the simplest
				      * way for elements can be
				      * interpolated from its
				      * children, interpolation is not
				      * possible for piecewise
				      * constant elements, to name
				      * only one example.  Here, the
				      * value on the mother cell has
				      * to be taken the average of the
				      * values on the children,
				      * i.e. all children contribute
				      * alike to the one degree of
				      * freedom. Here, we have to sum
				      * up the contributions of all
				      * child cells with the same
				      * weight, and the non-additive
				      * mode of above would only set
				      * the value on the mother cell
				      * to the value of one of the
				      * child cell, irrespective of
				      * the values on the other cells.
				      *
				      * Similarly, for discontinuous
				      * linear elements, it might be
				      * better to not interpolate the
				      * values at the corners from the
				      * child cells, but to take a
				      * better average, for example
				      * interpolating at the centers
				      * of the child cells; in that
				      * case, the contributions of the
				      * child cells have to be
				      * additive as well.
				      *
				      * Given these notes, the flag
				      * under consideration has to be
				      * set to #false# for the usual
				      * continuous Lagrange elements,
				      * and #true# for the other cases
				      * mentioned above. The main
				      * function where it is used is
				      * #DoFAccessor::get_interpolated_dof_values#.
				      * There is one flag per
				      * component.
				      */
    const vector<bool> restriction_is_additive_flags;    
};


/**
 * Finite Element in any dimension. This class declares the
 * functionality to fill the fields of the #FiniteElementBase#
 * class. Since this is something that depends on the actual finite
 * element, the functions are declared virtual if it is not possible
 * to provide a reasonable standard implementation.
 *
 *
 * \subsection{Finite elements in one dimension}
 *
 * Finite elements in one dimension need only set the #restriction# and
 * #prolongation# matrices in #FiniteElementBase<1>#. The constructor of
 * this class in one dimension presets the #interface_constraints# matrix
 * by the unit matrix with dimension one. Changing this behaviour in
 * derived classes is generally not a reasonable idea and you risk getting
 * in terrible trouble.
 * 
 * 
 * \subsection{Finite elements in two dimensions}
 * 
 * In addition to the fields already present in 1D, a constraint matrix
 * is needed in case two quads meet at a common line of which one is refined
 * once more than the other one. Then there are constraints referring to the
 * hanging nodes on that side of the line which is refined. These constraints
 * are represented by a $m\times n$-matrix #interface_constraints#, where $n$ is the
 * number of degrees of freedom on the refined side (those dofs on the middle
 * vertex plus those on the two lines), and $m$ is that of the unrefined side
 * (those dofs on the two vertices plus those on the line). The matrix is thus
 * a rectangular one.
 *
 * The mapping of the dofs onto the indices of the matrix is as follows:
 * let $d_v$ be the number of dofs on a vertex, $d_l$ that on a line, then
 * $m=0...d_v-1$ refers to the dofs on vertex zero of the unrefined line,
 * $m=d_v...2d_v-1$ to those on vertex one,
 * $m=2d_v...2d_v+d_l-1$ to those on the line.
 *
 * Similarly, $n=0...d_v-1$ refers to the dofs on the middle vertex
 * (vertex one of child line zero, vertex zero of child line one),
 * $n=d_v...d_v+d_l-1$ refers to the dofs on child line zero,
 * $n=d_v+d_l...d_v+2d_l-1$ refers to the dofs on child line one.
 * Please note that we do not need to reserve space for the dofs on the
 * end vertices of the refined lines, since these must be mapped one-to-one
 * to the appropriate dofs of the vertices of the unrefined line.
 *
 * It should be noted that it is not possible to distribute a constrained
 * degree of freedom to other degrees of freedom which are themselves
 * constrained. Only one level of indirection is allowed. It is not known
 * at the time of this writing whether this is a constraint itself.
 *
 *
 * \subsection{Finite elements in three dimensions}
 *
 * For the interface constraints, almost the same holds as for the 2D case.
 * The numbering for the indices $m$ on the mother face is obvious and keeps
 * to the usual numbering of degrees of freedom on quadrilaterals.
 *
 * The numbering of the degrees of freedom on the interior of the refined
 * faces for the index $n$ is as follows: let $d_v$ and $d_l$ be as above,
 * and $d_q$ be the number of degrees of freedom per quadrilateral (and
 * therefore per face), then $n=0...d_v-1$ denote the dofs on the vertex at
 * the center, $n=d_v...5d_v-1$ for the dofs on the vertices at the center
 * of the bounding lines of the quadrilateral,
 * $n=5d_v..5d_v+4*d_l-1$ are for the degrees of freedom on
 * the four lines connecting the center vertex to the outer boundary of the
 * mother face, $n=5d_v+4*d_l...5d_v+4*d_l+8*d_l-1$ for the degrees of freedom
 * on the small lines surrounding the quad,
 * and $n=5d_v+12*d_l...5d_v+12*d_l+4*d_q-1$ for the dofs on the
 * four child faces. Note the direction of the lines at the boundary of the
 * quads, as shown below.
 *
 * The order of the twelve lines and the four child faces can be extracted
 * from the following sketch, where the overall order of the different
 * dof groups is depicted:
 * \begin{verbatim}
 *    *--13--3--14--*
 *    |      |      |
 *    16 20  7  19  12
 *    |      |      |
 *    4--8---0--6---2
 *    |      |      |
 *    15 17  5  18  11
 *    |      |      |
 *    *--9---1--10--*
 * \end{verbatim}
 * It should be noted that the face as shown here is in the standard form,
 * i.e. with vertex zero at the bottom left, and the other vertices numbered
 * counter clockwise. This explains the numbering of the lines labeled 13 and
 * 14, as well as those labeled 15 and 16. The dofs on the lines need to
 * be numbered in the direction of the lines, which is as follows:
 * \begin{verbatim}
 *    *-->---*-->---*
 *    |      |      |
 *    ^      ^      ^ 
 *    |      |      |
 *    *-->---*-->---*
 *    |      |      |
 *    ^      ^      ^ 
 *    |      |      |
 *    *-->---*-->---*
 * \end{verbatim}
 * The orientation of the quads should be obvious.
 *
 * The faces of a hexahedron are arranged in a way such that
 * some must be viewed from the inside and some from the outside of the cell to
 * show this order; refer to the documentation of the #Triangulation# class for
 * the definition of this.
 *
 * If of the cells adjacent to one line more than one is refined and there is
 * at least one unrefined cell, then the degrees of freedom on the refined line
 * are constrained from two cells. For example, consider the cell behind the
 * face shown above is refined, while the one in front of the face is not
 * refined; then the dofs on the lines numbered 9 and 10 are constrained. If
 * there are two more cells below the ones just introduced, with a common face
 * right below the one shown, and of these is one refined and one unrefined one,
 * then the degrees on the two mentioned small lines are constrained a second
 * time. Since these constraints must be unique, it follows that the constraints
 * for the degrees of freedom on refined lines may only be in terms of the
 * degrees of freedom on the unrefined line, not in terms of the other
 * degrees of freedom on a face.
 *
 * Since the handling of constraints on degrees of freedom is mostly done
 * by the #ConstraintMatrix# class, this class checks whether the constraints
 * introduced from the two sides are unique; it is able to handle the fact
 * that the constraints for some of the dofs are entered more than once.
 *
 *
 * \subsection{Notes on three space dimensions}
 *
 * In three space dimensions, using locally refined elements involves
 * a difficulty not found in one or two spatial dimensions: the common
 * face of two cells need not match exactly, if one of the cells is
 * refined and the two cells are at the boundary. To understand this,
 * look at the following sketch:
 * \begin{verbatim}
 *         *---------*---------*
 *        /         /         /|
 *       /         /         / |
 *      /         /         /  |
 *     *---------*---------*   |
 *     |         |         |   |
 *     |         |         |   *
 *     |         |         |  /
 *     |         |         | /
 *     |         |         |/
 *     *---------*---------*
 * \end{verbatim}
 *
 * Assume the two top faces represent the boundary of the
 * triangulation; assume further that the boundary of the original
 * domain is curved here. Then, refining one of the two cells will
 * lead to refinement of the top line of the common face of the two
 * cells with the new mid-point being raised or lowered, i.e. the two
 * children of this line will not take the same place as their mother
 * line did (this is not properly drawable using only ASCII characters,
 * use some imagination):
 * \begin{verbatim}
 *       ..*--.*---..*---------*
 *      *----*----* /         /|
 *      :    :    :/         / |
 *     :    :    :/         /  |
 *     *----*----*---------*   |
 *     |    |    |         |   |
 *     |    |    |         |   *
 *     *----*----*         |  /
 *     |    |    |         | /
 *     |    |    |         |/
 *     *----*----*---------*
 * \end{verbatim}
 * While this is the case with boundary faces in two spatial
 * dimensions also, it here leads to the fact that the four child
 * faces of the common face of the two cells will not coincide with
 * their mother face as well.
 *
 * Before proceeding to the consequences of this, we should note that
 * this problem occurs only for cells exactly at the boundary and if
 * exactly one of the two cells is refined once. If one of the two is
 * refined once and the other one twice, the problem again occurs only
 * for the outermost layer of cells, not for the others.
 *
 * Now for the consequences. Because most finite elements, at least
 * those implemented at present (February 1999) are implemented by
 * interpolation to certain Lagrange points, and because the Lagrange
 * points do not match, there is no easy way to obtain continuity or
 * any other constraint on the finite element functions at this
 * face. This is rather obvious since parts of the child faces of the
 * left, refined cell do not match any face of the right, unrefined
 * cell at all. This problem seems unsolvable using the usual finite
 * elements with trial functions computed on the unit cell without
 * taking into consideration the actual cell, so we do not even
 * attempt to solve it.
 *
 * A second, but related problem comes into play when trying to
 * compute integrals over faces which are refined from one side. For
 * this problem, the #FESubfaceValues# class exists, and it evaluates
 * certain functions of the finite element class involving the
 * Jacobian determinant of the mapping of unit face to real face,
 * restricted to a subface, and the normal vectors to the subfaces. We
 * should note that here, we talk only about evaluating the finite
 * element in the right cell, but on the common face; evaluating the
 * finite element in the small cells on the left poses no problem. The
 * question here is: what are the subfaces? It could either be the
 * four subfaces of the refined cell to the left, or the four subfaces
 * of the large face if it were refined with no curved boundary being
 * near it. In the first case, points where we evaluate jacobians as
 * well as normal vectors would match from both sides of the faces;
 * however, the points at which the finite element function is
 * evaluated, would not match, which needs to be that way because for
 * some points of the small faces of the left cell there are no
 * matching points on the right.
 *
 * The other possibility would be to totally ignore the existence of
 * the boundary and evaluate the finite element in the right cell at
 * subfaces which would be generated if the new vertex of the top line
 * of the common face was the midpoint of the line. This approach is
 * simpler from the implementational view point, but is also more
 * appropriate, since we evaluate on the right cell and do not want to
 * let this evaluation depend on the state of the left cell or its
 * children.
 *
 * Within this library, the present implementation uses the second way.
 *
 *
 * \subsection{Notes on extending the finite element library}
 *
 * The #deal.II# library was mainly made to use lagrange elements of arbitrary
 * order. For this reason, there may be places in the library where it uses
 * features of finite elements which may not be as general as desirable as may
 * be. Most of these restrictions don't come to mind and may cause problems
 * if someone wanted to implement a finite element which does not satisfy these
 * restrictions, leading to strange problems in places one does not expect.
 *
 * This section tries to collect some of these restrictions which are known.
 * There is no guarantee that this list is complete; in fact, doubts are in
 * place that that be so.
 *
 * \begin{itemize}
 * \item Lagrange elements: at several places in the library, use is made of the
 *   assumption that the basis functions of a finite element corresponds to a
 *   function value (as opposed to derivatives or the like, as used in the
 *   Hermitean finite element class or in the quintic Argyris element). It is
 *   further assumed that a basis function takes its nominal value at a
 *   certain point (e.g. linear trial functions take their value in the
 *   corners of the element; this restriction rules out spectral elements for
 *   the present library).
 *
 *   Both these assumptions are used when interpolation of a continuous
 *   function to the finite element space is applied. At present, only few
 *   places where this is used in the library come to mind to the author,
 *   namely the treating of boundary values in the #ProblemBase# class and
 *   the interpolation in the #VectorTools# collection. You should also
 *   look out for other places where explicit use of the support points is
 *   made if you want to use elements of other classes. A hint may be the
 *   use of the #get_support_points# and #get_face_support_points# functions
 *   of this class.
 *
 *   This also is used in some sense in the
 *   #DoFHandler::distribute_cell_to_dof_vector# where it is assumed that
 *   the degrees of freedom denote function values and not derivatives or
 *   the like.
 *
 * \item Vanishing of basis functions on faces: when projecting a function
 *   to the boundary, use if made of the assumption that all basis functions
 *   on a cell adjacent to the boundary vanish on the boundary except for those
 *   on the boundary face itself. For Lagrange elements this is true, but it
 *   may or may not be true in other cases.
 *
 *   This assumption is used in the #VectorTools::project_boundary_values#,
 *   #MatrixCreator::create_boundary_mass_matrix#,
 *   #DoFHandler::make_boundary_sparsity_pattern#,
 *   #DoFHandler::map_dof_to_boundary_indices# and may be a few other places.
 *   The places in the #DoFHandler# class are probably not that dangerous,
 *   since wrong results will most likely lead to internal errors through
 *   the #Assert# mechanism, but the first places will lead to undiscovered
 *   errors if not thought of properly.
 *
 *   This assumption also comes into play when computing the constraints of
 *   hanging nodes. If functions not located on a certain face vanish on
 *   that face (they do for Lagrange elements), then the distribution of
 *   constrained nodes happens with the face nodes on the large call. If
 *   the assumption does not hold, then the distribution has to happen
 *   with all nodes on the small and the large cells. This is not
 *   implemented in the #DoFHandler# class as of now.
 * \end{itemize}
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
class FiniteElement : public FiniteElementBase<dim>
{
  public:
				     /**
				      * Constructor
				      */
    FiniteElement (const FiniteElementData<dim> &fe_data,
		   const vector<bool> &restriction_is_additive_flags);

				     /**
				      * Destructor. Only declared to have a
				      * virtual destructor which the compiler
				      * wants to have.
				      */
    virtual ~FiniteElement () {};
    
				     /**
				      * Return the value of the #i#th shape
				      * function at the point #p#.
				      * #p# is a point on the reference element.
				      */
    virtual double shape_value (const unsigned int i,
			        const Point<dim> &p) const = 0;

				     /**
				      * Return the gradient of the #i#th shape
				      * function at the point #p#.
				      * #p# is a point on the reference element,
				      */
    virtual Tensor<1,dim> shape_grad (const unsigned int  i,
				      const Point<dim>   &p) const = 0;

				     /**
				      * Return the tensor of second derivatives
				      * of the #i#th shape function at
				      * point #p# on the unit cell.
				      */
    virtual Tensor<2,dim> shape_grad_grad (const unsigned int  i,
					   const Point<dim>   &p) const = 0;

				     /**
				      * Return the value of the #i#th shape
				      * function of the transformation mapping
				      * from unit cell to real cell. For
				      * isoparametric elements, this function
				      * is the same as the trial functions,
				      * but for sublinear or other mappings,
				      * they differ.
				      */
    virtual double shape_value_transform (const unsigned int i,
					  const Point<dim> &p) const = 0;

				     /**
				      * Same as above: return gradient of the
				      * #i#th shape function for the mapping
				      * from unit to real cell.
				      */
    virtual Tensor<1,dim> shape_grad_transform (const unsigned int i,
						const Point<dim> &p) const = 0;    
    
				     /**
				      * Compute the Jacobian matrix and the
				      * quadrature points as well as the trial
				      * function locations on the real cell in
				      * real space from the given cell
				      * and the given quadrature points on the
				      * unit cell. The Jacobian matrix is to
				      * be computed at every quadrature point.
				      * The derivative of the jacobian matrix
				      * is the derivative with respect to the
				      * unit cell coordinates.
				      * This function has to be in the finite
				      * element class, since different finite
				      * elements need different transformations
				      * of the unit cell to a real cell.
				      *
				      * The computation of these fields may
				      * share some common code, which is why we
				      * put it in one function. However, it may
				      * not always be necessary to really
				      * compute all fields, so there are
				      * bool flags which tell the function which
				      * of the fields to actually compute.
				      *
				      * Refer to the documentation of the
				      * \Ref{FEValues} class for a definition
				      * of the Jacobi matrix and of the various
				      * structures to be filled.
				      *
				      * This function is provided for
				      * the finite element class in
				      * one space dimension, but for
				      * higher dimensions, it depends
				      * on the present fe and needs
				      * reimplementation by the
				      * user. This is due to the fact
				      * that the user may want to use
				      * iso- or subparametric mappings
				      * of the unit cell to the real
				      * cell, which makes things much
				      * more complicated.
				      *
				      * The
				      * #shape_values/grads_transform#
				      * arrays store the values and
				      * gradients of the
				      * transformation basis
				      * functions.  While this
				      * information is not necessary
				      * for the computation of the
				      * other fields, it allows for
				      * significant speedups, since
				      * the values and gradients of
				      * the transform functions at the
				      * quadrature points need not be
				      * recomputed each time this
				      * function is called.
				      *
				      * The function assumes that the fields
				      * already have the right number of
				      * elements. It has to be
				      * guaranteed, that fields that are
				      * not requested for update are not changed.
				      * This also means, that these
				      * fields have to be filled with
				      * the correct values beforehand.
				      *
				      * This function is more or less an
				      * interface to the #FEValues# class and
				      * should not be used by users unless
				      * absolutely needed.
				      */
    virtual void fill_fe_values (const DoFHandler<dim>::cell_iterator &cell,
				 const vector<Point<dim> > &unit_points,
				 vector<Tensor<2,dim> >    &jacobians,
				 const bool                 compute_jacobians,
				 vector<Tensor<3,dim> >    &jacobians_grad,
				 const bool                 compute_jacobians_grad,
				 vector<Point<dim> >       &support_points,
				 const bool                 compute_support_points,
				 vector<Point<dim> >       &q_points,
				 const bool                 compute_q_points,
				 const FullMatrix<double>            &shape_values_transform,
				 const vector<vector<Tensor<1,dim> > > &shape_grads_transform) const;

				     /**
				      * Do the same thing that the other
				      * #fill_fe_values# function does,
				      * exception that a face rather than
				      * a cell is considered. The #face_no#
				      * parameter denotes the number of the
				      * face to the given cell to be
				      * considered.
				      *
				      * The unit points for the quadrature
				      * formula are given on the unit face
				      * which is a mannifold of dimension
				      * one less than the dimension of the
				      * cell. The #global_unit_points# 
				      * denote the position of the unit points
				      * on the selected face on the unit cell.
				      * This additional information is passed
				      * since the #FEFaceValues# class can
				      * compute them once and for all,
				      * eliminating the need to recompute it
				      * each time #FEFaceValues::reinit# is
				      * called.
				      *
				      * The jacobian matrix is evaluated at
				      * each of the quadrature points on the
				      * given face. The matrix is the
				      * transformation matrix of the unit cell
				      * to the real cell, not from the unit
				      * face to the real face. This is the
				      * necessary matrix to compute the real
				      * gradients.
				      *
				      * Conversely, the Jacobi determinants
				      * are the determinants of the
				      * transformation from the unit face to
				      * the real face. This information is
				      * needed to actually perform integrations
				      * along faces. Note that we here return
				      * the inverse of the determinant of the
				      * jacobi matrices as explained in the
				      * documentation of the #FEValues# class.
				      * 
				      * The support points are the
				      * off-points of those trial functions
				      * located on the given face; this
				      * information is taken over from the
				      * #get_face_support_points# function.
				      *
				      * The order of trial functions is the
				      * same as if it were a cell of dimension
				      * one less than the present. E.g. in
				      * two dimensions, the order is first
				      * the vertex functions (using the
				      * direction of the face induced by the
				      * given cell) then the interior functions.
				      * The same applies for the quadrature
				      * points which also use the standard
				      * direction of faces as laid down by
				      * the #Triangulation# class.
				      *
				      * There is a standard implementation for
				      * dimensions greater than one. It
				      * uses the #fill_fe_values()#
				      * function to retrieve the wanted
				      * information. Since this operation acts
				      * only on unit faces and cells it does
				      * not depend on a specific finite element
				      * transformation and is thus applicable
				      * for all finite elements and uses tha
				      * same mapping from the unit to the real
				      * cell as used for the other operations
				      * performed by the specific finite element
				      * class.
				      *
				      * Three fields remain to be finite element
				      * specific in this standard implementation:
				      * The jacobi determinants of the
				      * transformation from the unit face to the
				      * real face, the support points
				      * and the outward normal vectors. For
				      * these fields, there exist pure
				      * virtual functions, #get_face_jacobians#,
				      * #get_face_support_points# and
				      * #get_normal_vectors#.
				      *
				      * Though there is a standard
				      * implementation, there
				      * may be room for optimizations which is
				      * why this function is made virtual.
				      *
				      * Since any implementation for one
				      * dimension would be senseless, all
				      * derived classes should throw an error
				      * when called with #dim==1#.
				      *
				      * The function assumes that the fields
				      * already have the right number of
				      * elements.
				      *
				      * This function is more or less an
				      * interface to the #FEFaceValues# class
				      * and should not be used by users unless
				      * absolutely needed.
				      */
    virtual void fill_fe_face_values (const DoFHandler<dim>::cell_iterator &cell,
				      const unsigned int           face_no,
				      const vector<Point<dim-1> > &unit_points,
				      const vector<Point<dim> >   &global_unit_points,
				      vector<Tensor<2,dim> >      &jacobians,
				      const bool                   compute_jacobians,
				      vector<Tensor<3,dim> >      &jacobians_grad,
				      const bool                   compute_jacobians_grad,
				      vector<Point<dim> > &support_points,
				      const bool           compute_support_points,
				      vector<Point<dim> > &q_points,
				      const bool           compute_q_points,
				      vector<double>      &face_jacobi_determinants,
				      const bool           compute_face_jacobians,
				      vector<Point<dim> > &normal_vectors,
				      const bool           compute_normal_vectors,
				      const FullMatrix<double>      &shape_values_transform,
				      const vector<vector<Tensor<1,dim> > > &shape_grads_transform) const;

				     /**
				      * This function does almost the same as
				      * the above one, with the difference that
				      * it considers the restriction of a finite
				      * element to a subface (the child of a
				      * face) rather than to a face. The number
				      * of the subface in the face is given by
				      * the #subface_no# parameter. The meaning
				      * of the other parameters is the same as
				      * for the #fill_fe_face_values# function.
				      *
				      * Since the usage of support points on
				      * subfaces is not useful, it is excluded
				      * from the interface to this function.
				      *
				      * Like for the #fill_fe_face_values#
				      * function, there is a default
				      * implementation, using the
				      * #fill_fe_values# function. There may
				      * be better and more efficient solutions
				      * for a special finite element, which is
				      * why this function is made virtual.
				      *
				      * This function is more or less an
				      * interface to the #FESubfaceValues# class
				      * and should not be used by users unless
				      * absolutely needed.
				      */				       
    virtual void fill_fe_subface_values (const DoFHandler<dim>::cell_iterator &cell,
					 const unsigned int           face_no,
					 const unsigned int           subface_no,
					 const vector<Point<dim-1> > &unit_points,
					 const vector<Point<dim> >   &global_unit_points,
					 vector<Tensor<2,dim> >      &jacobians,
					 const bool                   compute_jacobians,
					 vector<Tensor<3,dim> >      &jacobians_grad,
					 const bool           compute_jacobians_grad,
					 vector<Point<dim> > &q_points,
					 const bool           compute_q_points,
					 vector<double>      &face_jacobi_determinants,
					 const bool           compute_face_jacobians,
					 vector<Point<dim> > &normal_vectors,
					 const bool           compute_normal_vectors,
					 const FullMatrix<double>      &shape_values_transform,
					 const vector<vector<Tensor<1,dim> > > &shape_grads_transform) const;

				     /**
				      * Return the support points of the
				      * trial functions on the unit cell.
				      *
				      * The function assumes that the
				      * #unit_points# array already has the
				      * right size. The order of points in
				      * the array matches that returned by
				      * the #cell->get_dof_indices# function.
				      *
				      * For one space dimension there is a
				      * standard implementation assuming
				      * equidistant off-points on the unit
				      * line. For all other dimensions, an
				      * overwritten function has to be provided.
				      */
    virtual void get_unit_support_points (vector<Point<dim> > &unit_points) const;
    
				     /**
				      * Compute the off-points of the finite
				      * element basis functions on the given
				      * cell in real space.
				      *
				      * This function implements a subset of
				      * the information delivered by the
				      * #fill_fe_values# function to the
				      * #FEValues# class. However, since it
				      * is useful to use information about
				      * off-points without using #FEValues#
				      * objects (e.g. in interpolating functions
				      * to the finite element space), this
				      * function is excluded from the
				      * abovementioned one.
				      *
				      * The function assumes that the
				      * #support_points# array already has the
				      * right size. The order of points in
				      * the array matches that returned by
				      * the #cell->get_dof_indices# function.
				      *
				      * For one space dimension there is a
				      * standard implementation assuming
				      * equidistant off-points on the unit
				      * line. For all other dimensions, an
				      * overwritten function has to be provided.
				      *
				      * For higher order transformations than
				      * the common (bi-, tri-)linear one,
				      * information about the boundary is
				      * needed, rather than only the readily
				      * available information on the location
				      * of the vertices. If necessary, we
				      * therefore rely on the boundary object
				      * of which a pointer is stored by the
				      * triangulation.
				      */
    virtual void get_support_points (const DoFHandler<dim>::cell_iterator &cell,
				     vector<Point<dim> > &support_points) const;
    
				     /**
				      * Compute the off-points of the finite
				      * element basis functions located on the
				      * face. It only returns the off-points
				      * of the trial functions which are
				      * located on the face, rather than of
				      * all basis functions, which is done by
				      * the #get_support_points# function.
				      *
				      * This function produces a subset of
				      * the information provided by the
				      * #fill_fe_face_values()# function.
				      * However, you should not try
				      * to implement this function using the
				      * abovementioned function, since usually
				      * that function uses this function to
				      * compute information.
				      *
				      * The function is excluded from the
				      * abovementioned one, since no information
				      * about the neighboring cell is needed,
				      * such that loops over faces alone are
				      * possible when using this function.
				      * This is useful for example if we want
				      * to interpolate boundary values to the
				      * finite element functions. If integration
				      * along faces is needed, we still need
				      * the #fill_fe_face_values# function.
				      *
				      * The function assumes that the
				      * #support_points# array already has the
				      * right size. The order of points in
				      * the array matches that returned by
				      * the #face->get_dof_indices# function.
				      *
				      * Since any implementation for one
				      * dimension would be senseless, all
				      * derived classes should throw an error
				      * when called with #dim==1#.
				      *
				      * Regarding information about the
				      * boundary, which is necessary for
				      * higher order transformations than
				      * the usual (bi-, tri-)linear ones,
				      * refer to the #get_support_points#
				      * function.
				      */
    virtual void get_face_support_points (const DoFHandler<dim>::face_iterator &face,
					  vector<Point<dim> > &support_points) const =0;

				     /**
				      * This is the second separated function
				      * described in the documentation of the
				      * #fill_fe_face_values# function. It
				      * returns the determinants of the
				      * transformation from the unit face to the
				      * real face at the
				      *
				      * Since any implementation for one
				      * dimension would be senseless, all
				      * derived classes should throw an error
				      * when called with #dim==1#.
				      *
				      * Regarding information about the
				      * boundary, which is necessary for
				      * higher order transformations than
				      * the usual (bi-, tri-)linear ones,
				      * refer to the #get_support_points#
				      * function.
				      */
    virtual void get_face_jacobians (const DoFHandler<dim>::face_iterator &face,
				     const vector<Point<dim-1> > &unit_points,
				     vector<double>      &face_jacobi_determinants) const =0;

				     /**
				      * Does the same as the above function,
				      * except that it computes the Jacobi
				      * determinant of the transformation from
				      * the unit face to the subface of #face#
				      * with number #subface_no#.
				      *
				      * The function needs not take special care
				      * about boundary approximation, since it
				      * must not be called for faces at the
				      * boundary.
				      */
    virtual void get_subface_jacobians (const DoFHandler<dim>::face_iterator &face,
					const unsigned int           subface_no,
					const vector<Point<dim-1> > &unit_points,
					vector<double>      &face_jacobi_determinants) const =0;

				     /**
				      * Compute the normal vectors to the cell
				      * at the quadrature points. See the
				      * documentation for the #fill_fe_face_values#
				      * function for more details. The function
				      * must guarantee that the length of the
				      * vectors be one.
				      *
				      * Since any implementation for one
				      * dimension would be senseless, all
				      * derived classes should throw an error
				      * when called with #dim==1#.
				      *
				      * Regarding information about the
				      * boundary, which is necessary for
				      * higher order transformations than
				      * the usual (bi-, tri-)linear ones,
				      * refer to the #get_support_points#
				      * function.
				      */
    virtual void get_normal_vectors (const DoFHandler<dim>::cell_iterator &cell,
				     const unsigned int          face_no,
				     const vector<Point<dim-1> > &unit_points,
				     vector<Point<dim> >         &normal_vectors) const =0;

				     /**
				      * Does the same as the above function,
				      * except that it refers to the
				      * subface #subface_no# of the given face.
				      *
				      * The function needs not take special care
				      * about boundary approximation, since it
				      * must not be called for faces at the
				      * boundary.
				      */
    virtual void get_normal_vectors (const DoFHandler<dim>::cell_iterator &cell,
				     const unsigned int           face_no,
				     const unsigned int           subface_no,
				     const vector<Point<dim-1> > &unit_points,
				     vector<Point<dim> >         &normal_vectors) const =0;

				     /**
				      * Fill in the given matrix with the local
				      * mass matrix. The mass matrix must be
				      * exactly computed, not using a
				      * quadrature, which may be done using
				      * an equation object and an assembler,
				      * as is done for the Laplace matrix
				      * in the #MatrixTools# class for example.
				      *
				      * The exact integration is possible since
				      * an exact representation for the Jacobi
				      * determinant exists in all known cases of
				      * iso- or subparametric mappings. For
				      * example, usually the point in real
				      * space $\vec x$ referring to the point
				      * $\vec \xi$ on the unit cell is given
				      * by $\vec x = \sum_i \vec p_i \phi_i(\vec \xi)$,
				      * where the sum is over all basis functions
				      * $\phi_i$ and $\vec p_i$ are the points
				      * in real space where the basis function
				      * $\phi_i$ is located. The Jacobi
				      * determinant is the given by
				      * $|det J| = |\frac{\partial\vec x}{\partial\vec\xi}$,
				      * which can be evaluated in closed form.
				      * The mass matrix then is given by
				      * $m_{ij} = \int_{\hat K} \phi_i(\vec\xi)
				      * \phi_j(\vec\xi) |det J| d\xi$, where
				      * $\hat K$ is the unit cell. The integrand
				      * obviously is a polynom and can thus
				      * easily be integrated analytically, so
				      * the computation of the local mass matrix
				      * is reduced to the computation of a
				      * weighted evaluation of a polynom in
				      * the coordinates of the support points
				      * in real space (for linear mappings,
				      * these are the corner points, for
				      * quadratic mappings also the center of
				      * mass and the edge and face centers).
				      * For example, in one space dimension,
				      * the Jacobi determinant simply is $h$,
				      * the size of the cell, and the integral
				      * over the two basis functions can easily
				      * be calculated with a pen and a sheet of
				      * paper. The actual computation on this
				      * matrix then is simply a scaling of a
				      * known and constant matrix by $h$.
				      *
				      * The functions which override this one
				      * may make assumptions on the sign of
				      * the determinant if stated in the
				      * documentation, but should check for
				      * them in debug mode. For that purpose,
				      * an exception with the longish name
				      * #ExcJacobiDeterminantHasWrongSign#
				      * is declared.
				      *
				      * The function takes a #DoFHandler#
				      * iterator, which provides a superset
				      * of information to the geometrical
				      * information needed for the computations.
				      * The additional data should not be
				      * used, however a #DoFHandler# iterator
				      * was preferred over a #Triangulation#
				      * iterator since this is what usually
				      * is available in places where this
				      * function is called.
				      *
				      * The cell matrix is assumed to be of
				      * the right size already. Functions
				      * of derived classes shall be implemented
				      * in a way as to overwrite the previous
				      * contents of the matrix, so it need not
				      * be necessary to clear the matrix before
				      * use with this function.
				      *
				      * Some finite elements, especially in
				      * higher dimensions, may chose not to
				      * implement this function because the
				      * computational effort is growing
				      * rapidly, for the in-time computation
				      * of the matrix as well as for the
				      * setting up using a script. For
				      * example, the size of the generated
				      * #C++# code for the local mass
				      * matrix in 3d is 4.383.656 bytes
				      * already for the trilinear element.
				      * Higher order elements would
				      * produce even larger code.
				      *
				      * In the case of a finite element chosing
				      * not to implement the functionality of
				      * this function, that function is supposed
				      * to throw an exception of class
				      * #ExcComputationNotUseful# declared
				      * in this class, for example through the
				      * #AssertThrow# mechanism; you can catch
				      * this exception and compute the mass matrix
				      * by quadrature instead. Finite element
				      * classes not implementing this function
				      * are assumed to state this in their
				      * documentation.
				      *
				      * Regarding information about the
				      * boundary, which is necessary for
				      * higher order transformations than
				      * the usual (bi-, tri-)linear ones,
				      * refer to the #get_support_points#
				      * function.
				      */
    virtual void get_local_mass_matrix (const DoFHandler<dim>::cell_iterator &cell,
					FullMatrix<double>            &local_mass_matrix) const =0;

				     /**
				      * Number of base elements in a mixed
				      * discretization. This function returns
				      * 1 for simple elements.
				      */
    virtual unsigned int n_base_elements() const;
    
				     /**
				      * Access to base element
				      * objects.  By default,
				      * #base_element(0)# is #this#.
				      * This function is overloaded by
				      * system elements to allow
				      * access to the different
				      * components of mixed
				      * discretizations.
				      */
    virtual const FiniteElement<dim>& base_element(unsigned int index) const;
    
				     /**
				      * Exception
				      */
    DeclException0 (ExcPureFunctionCalled);
				     /**
				      * Exception
				      */
    DeclException0 (ExcBoundaryFaceUsed);
				     /**
				      * Exception
				      */
    DeclException0 (ExcJacobiDeterminantHasWrongSign);
				     /**
				      * Exception
				      */
    DeclException1 (ExcComputationNotUseful,
		    int,
		    << "The computation you required from this function is not "
		    << "feasible or not probable in the present dimension ("
		    << arg1 << ") because it would be prohibitively expensive.");
};


/* ------------------------------- Inline functions ----------------------- */


template <int dim>
inline
unsigned int 
FiniteElementData<dim>::n_dofs_per_vertex () const
{
  return dofs_per_vertex;
};


template <int dim>
inline
unsigned int 
FiniteElementData<dim>::n_dofs_per_line () const
{
  return dofs_per_line;
};


template <int dim>
inline
unsigned int 
FiniteElementData<dim>::n_dofs_per_quad () const
{
  return dofs_per_quad;
};


template <int dim>
inline
unsigned int 
FiniteElementData<dim>::n_dofs_per_hex () const
{
  return dofs_per_hex;
};


template <int dim>
inline
unsigned int 
FiniteElementData<dim>::n_dofs_per_face () const
{
  return dofs_per_face;
};


template <int dim>
inline
unsigned int 
FiniteElementData<dim>::n_dofs_per_cell () const
{
  return dofs_per_cell;
};


template <int dim>
inline
unsigned int 
FiniteElementData<dim>::n_components () const
{
  return components;
};


template <int dim>
inline
unsigned int 
FiniteElementData<dim>::n_transform_functions () const
{
  return transform_functions;
};


template <int dim>
inline
unsigned int
FiniteElementBase<dim>::component_to_system_index (unsigned int component,
						   unsigned int component_index) const
{
  Assert(component<n_components(),
	 ExcIndexRange(component, 0, n_components()));
  Assert(component_index<component_to_system_table[component].size(),
	 ExcIndexRange(component_index, 0,
		       component_to_system_table[component].size()));
  return component_to_system_table[component][component_index];
}


template <int dim>  
inline
pair<unsigned int,unsigned int>
FiniteElementBase<dim>::system_to_component_index (unsigned int index) const
{
  Assert(index < system_to_component_table.size(),
	 ExcIndexRange(index, 0, system_to_component_table.size()));
  return system_to_component_table[index];
}


template <int dim>
inline
unsigned int
FiniteElementBase<dim>::face_component_to_system_index (unsigned int component,
							unsigned int component_index) const
{
  Assert(component<n_components(),
	 ExcIndexRange(component, 0, n_components()));
  Assert(component_index<face_component_to_system_table[component].size(),
	 ExcIndexRange(component_index, 0,
		       face_component_to_system_table[component].size()));
  return face_component_to_system_table[component][component_index];
}


template <int dim>  
inline
pair<unsigned int,unsigned int>
FiniteElementBase<dim>::face_system_to_component_index (unsigned int index) const
{
  Assert(index < face_system_to_component_table.size(),
	 ExcIndexRange(index, 0, face_system_to_component_table.size()));
  return face_system_to_component_table[index];
}


template <int dim>  
inline
unsigned int
FiniteElementBase<dim>::component_to_base (unsigned int index) const
{
  if (n_components() == 1)
    return 0;
  Assert(index < component_to_base_table.size(),
	 ExcIndexRange(index, 0, component_to_base_table.size()));
  return component_to_base_table[index];
}


template <int dim>
inline
bool
FiniteElementBase<dim>::restriction_is_additive (const unsigned int component) const
{
  Assert(component<n_components(),
	 ExcIndexRange(component, 0, n_components()));
  return restriction_is_additive_flags[component];
}


#endif
