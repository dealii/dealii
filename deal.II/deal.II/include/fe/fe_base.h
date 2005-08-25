//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__fe_base_h
#define __deal2__fe_base_h

#include <base/config.h>
#include <base/exceptions.h>
#include <base/subscriptor.h>
#include <base/point.h>
#include <base/tensor.h>
#include <base/table.h>
#include <base/vector_slice.h>
#include <base/geometry_info.h>
#include <lac/full_matrix.h>
#include <fe/fe_update_flags.h>
#include <fe/mapping.h>

#include <string>
#include <vector>

template<int dim> class FESystem;

/*!@addtogroup febase */
/*@{*/

/**
 * Dimension independent data for finite elements. See the derived
 * class FiniteElement class for information on its use. All
 * its data are available to the implementation in a concrete finite
 * element class.
 *
 * Remark on a change in implementation: it is now wrong to cast a
 * pointer to FiniteElement to a pointer to FiniteElementData and
 * delete it. The virtual destructor has been moved up. In a later
 * version, FiniteElementData and FiniteElement should be private
 * base classes of FiniteElement.
 *
 * @author Wolfgang Bangerth, Guido Kanschat, 1998, 1999, 2000, 2001, 2003, 2005
 */
template <int dim>
class FiniteElementData
{
  public:
				     /**
				      * Enumerator for the different
				      * types of continuity a finite
				      * element may have. Continuity
				      * is measured by the Sobolev
				      * space containing the
				      * constructed finite element
				      * space and is also called this
				      * way.
				      *
				      * Note that certain continuities
				      * may imply others. For
				      * instance, a function in
				      * <i>H<sup>1</sup></i> is in
				      * <i>H<sup>curl</sup></i> and
				      * <i>H<sup>div</sup></i> as
				      * well.
				      *
				      * If you are interested in
				      * continuity in the classical
				      * sense, then the following
				      * relations hold:
				      *
				      * <ol>
				      *
				      * <li> <i>H<sup>1</sup></i>
				      * implies that the function is
				      * continuous over cell
				      * boundaries.
				      *
				      * <li> <i>H<sup>2</sup></i>
				      * implies that the function is
				      * continuously differentiable
				      * over cell boundaries.
				      *
				      * <li> <i>L<sup>2</sup></i>
				      * indicates that the element is
				      * discontinuous. Since
				      * discontinuous elements have no
				      * topological couplings between
				      * grid cells and code may
				      * actually depend on this
				      * property, <i>L<sup>2</sup></i>
				      * conformity is handled in a
				      * special way in the sense that
				      * it is <b>not</b> implied by
				      * any higher conformity.
				      * </ol>
				      *
				      * In order to test if a finite
				      * element conforms to a certain
				      * space, use
				      * FiniteElementData<dim>::conforms().
				      */
    enum Conformity
    {
					   /**
					    * Indicates incompatible
					    * continuities of a
					    * system.
					    */
	  unknown = 0x00,
	  
					   /**
					    * Discontinuous
					    * elements. See above!
					    */
	  L2 = 0x01,
	  
					   /**
					    * Conformity with the
					    * space
					    * <i>H<sup>curl</sup></i>
					    * (continuous tangential
					    * component of a vector
					    * field)
					    */
	  Hcurl = 0x02,
	  
					   /**
					    * Conformity with the
					    * space
					    * <i>H<sup>div</sup></i>
					    * (continuous normal
					    * component of a vector
					    * field)
					    */
	  Hdiv = 0x04,
	  
					   /**
					    * Conformity with the
					    * space
					    * <i>H<sup>1</sup></i>
					    * (continuous)
					    */
	  H1 = Hcurl | Hdiv,
	  
					   /**
					    * Conformity with the
					    * space
					    * <i>H<sup>2</sup></i>
					    * (continuously
					    * differentiable)
					    */
	  H2 = 0x0e
    };
    
				     /**
				      * Number of degrees of freedom on
				      * a vertex.
				      */
    const unsigned int dofs_per_vertex;

				     /** Number of degrees of freedom
				      *  in a line; not including the
				      *  degrees of freedom on the
				      *  vertices of the line.
				      */
    const unsigned int dofs_per_line;

				     /** Number of degrees of freedom
				      *  in a quadrilateral; not
				      *  including the degrees of
				      *  freedom on the lines and
				      *  vertices of the
				      *  quadrilateral.
				      */
    const unsigned int dofs_per_quad;

				     /** Number of degrees of freedom
				      *  in a hexahedron; not
				      *  including the degrees of
				      *  freedom on the
				      *  quadrilaterals, lines and
				      *  vertices of the hecahedron.
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
				      * Number of degrees of freedom
				      * on a face. This is the
				      * accumulated number of degrees
				      * of freedom on all the objects
				      * of dimension up to
				      * <tt>dim-1</tt> constituting a
				      * face.
				      */
    const unsigned int dofs_per_face;
    
				     /**
				      * Total number of degrees of freedom
				      * on a cell. This is the
				      * accumulated number of degrees
				      * of freedom on all the objects
				      * of dimension up to
				      * <tt>dim</tt> constituting a
				      * cell.
				      */
    const unsigned int dofs_per_cell;

				     /**
				      * Number of vector components of
				      * this finite element, and
				      * dimension of the image
				      * space. For vector-valued
				      * finite elements (i.e. when
				      * this number is greater than
				      * one), the number of vector
				      * components is in many cases
				      * equal to the number of base
				      * elements glued together with
				      * the help of the FESystem
				      * class. However, for elements
				      * like the Nedelec element, the
				      * number is greater than one
				      * even though we only have one
				      * base element.
				      */
    const unsigned int components;

				     /**
				      * Maximal polynomial degree of a
				      * shape function in a single
				      * coordinate direction.
				      */
    const unsigned int degree;

				     /**
				      * Indicate the space this element conforms to.
				      */
    const Conformity conforming_space;
    
    				     /**
				      * Default
				      * constructor. Constructs an
				      * element with no dofs. Checking
				      * n_dofs_per_cell() is therefore
				      * a good way to check if
				      * something went wrong.
				      */
    FiniteElementData ();

				     /**
				      * Constructor, computing all
				      * necessary values from the
				      * distribution of dofs to
				      * geometrcal objects.
				      *
				      * @param dofs_per_object Number
				      * of dofs on geometrical objects
				      * for each dimension. In this
				      * vector, entry 0 refers to dofs
				      * on vertices, entry 1 on lines
				      * and so on. Its length must be
				      * <i>dim+1</i>.
				      * @param n_components Number of
				      * vector components of the
				      * element.
				      * @param degree
				      * Maximal polynomial degree in a
				      * single direction.
				      * @param conformity The finite
				      * element space has continuity
				      * of this Sobolev space.
				      */
    FiniteElementData (const std::vector<unsigned int> &dofs_per_object,
		       const unsigned int n_components,
		       const unsigned int degree = deal_II_numbers::invalid_unsigned_int,
		       const Conformity conformity = unknown);

				     /**
				      * Number of dofs per vertex.
				      */
    unsigned int n_dofs_per_vertex () const;

				     /**
				      * Number of dofs per line. Not
				      * including dofs on lower
				      * dimensional objects.
				      */
    unsigned int n_dofs_per_line () const;

    				     /**
				      * Number of dofs per quad. Not
				      * including dofs on lower
				      * dimensional objects.
				      */
    unsigned int n_dofs_per_quad () const;

    				     /**
				      * Number of dofs per hex. Not
				      * including dofs on lower
				      * dimensional objects.
				      */
    unsigned int n_dofs_per_hex () const;

    				     /**
				      * Number of dofs per face,
				      * accumulating degrees of
				      * freedom of all lower
				      * dimensional objects.
				      */
    unsigned int n_dofs_per_face () const;

    				     /**
				      * Number of dofs per cell,
				      * accumulating degrees of
				      * freedom of all lower
				      * dimensional objects.
				      */
    unsigned int n_dofs_per_cell () const;

    				     /**
				      * Number of components.
				      */
    unsigned int n_components () const;

				     /**
				      * Maximal polynomial degree of a
				      * shape function in a single
				      * coordinate direction.
				      *
				      * This function can be used to
				      * determine the optimal
				      * quadrature rule.
				      */
    unsigned int tensor_degree () const;

				     /**
				      * Test whether a finite element
				      * space conforms to a certain
				      * Sobolev space.
				      *
				      * @note This function will
				      * return a true value even if
				      * the finite element space has
				      * higher regularity than asked
				      * for.
				      */
    bool conforms (const Conformity) const;
    
				     /**
				      * Comparison operator.
				      */
    bool operator == (const FiniteElementData &) const;
};




/**
 * Base class for finite elements in arbitrary dimensions. This class
 * provides several fields which describe a specific finite element
 * and which are filled by derived classes. It more or less only
 * offers the fields and access functions which makes it possible to
 * copy finite elements without knowledge of the actual type (linear,
 * quadratic, etc).
 *
 * The different matrices are initialized with the correct size, such
 * that in the derived (concrete) finite element classes, their
 * entries only have to be filled in; no resizing is needed. If the
 * matrices are not defined by a concrete finite element, they should
 * be resized to zero. This way functions using them can find out,
 * that they are missing. On the other hand, it is possible to use
 * finite element classes without implementation of the full
 * functionality, if only part of it is needed. The functionality
 * under consideration here is hanging nodes constraints and grid
 * transfer, respectively.
 *
 *
 * <h3>Support points</h3>
 *
 * Since a FiniteElement does not have information on the actual
 * grid cell, it can only provide support points on the unit
 * cell. Support points on the actual grid cell must be computed by
 * mapping these points. The class used for this kind of operation is
 * FEValues. In most cases, code of the following type will
 * serve to provide the mapped support points.
 *
 * @code
 * Quadrature<dim> dummy_quadrature (fe.get_unit_support_points());
 * FEValues<dim>   fe_values (mapping, fe, dummy_quadrature,
 *                            update_q_points);
 * fe_values.reinit (cell);
 * Point<dim>& mapped_point = fe_values.quadrature_point (i);
 * @endcode
 *
 * Alternatively, the points can be transformed one-by-one:
 * @code
 * const vector<Point<dim> >& unit_points =
 *    fe.get_unit_support_points();
 *
 * Point<dim> mapped_point =
 *    mapping.transform_unit_to_real_cell (cell, unit_points[i]);
 * @endcode
 * This is a shortcut, and as all shortcuts should be used cautiously.
 * If the mapping of all support points is needed, the first variant should
 * be preferred for efficiency.
 *
 * <h3>Finite elements in one dimension</h3>
 *
 * Finite elements in one dimension need only set the #restriction
 * and #prolongation matrices. The constructor of this class in one
 * dimension presets the #interface_constraints matrix to have
 * dimension zero. Changing this behaviour in derived classes is
 * generally not a reasonable idea and you risk getting into trouble.
 * 
 * <h3>Finite elements in two dimensions</h3>
 * 
 * In addition to the fields already present in 1D, a constraint
 * matrix is needed, if the finite element has node values located on
 * edges or vertices. These constraints are represented by an $m\times
 * n$-matrix #interface_constraints, where <i>m</i> is the number of
 * degrees of freedom on the refined side without the corner vertices
 * (those dofs on the middle vertex plus those on the two lines), and
 * <i>n</i> is that of the unrefined side (those dofs on the two
 * vertices plus those on the line). The matrix is thus a rectangular
 * one. The $m\times n$ size of the #interface_constraints matrix can
 * also be accessed through the interface_constraints_size() function.
 *
 * The mapping of the dofs onto the indices of the matrix on the
 * unrefined side is as follows: let $d_v$ be the number of dofs on a
 * vertex, $d_l$ that on a line, then $n=0...d_v-1$ refers to the dofs
 * on vertex zero of the unrefined line, $n=d_v...2d_v-1$ to those on
 * vertex one, $n=2d_v...2d_v+d_l-1$ to those on the line.
 *
 * Similarly, $m=0...d_v-1$ refers to the dofs on the middle vertex of
 * the refined side (vertex one of child line zero, vertex zero of
 * child line one), $m=d_v...d_v+d_l-1$ refers to the dofs on child
 * line zero, $m=d_v+d_l...d_v+2d_l-1$ refers to the dofs on child
 * line one.  Please note that we do not need to reserve space for the
 * dofs on the end vertices of the refined lines, since these must be
 * mapped one-to-one to the appropriate dofs of the vertices of the
 * unrefined line.
 *
 * It should be noted that it is not possible to distribute a constrained
 * degree of freedom to other degrees of freedom which are themselves
 * constrained. Only one level of indirection is allowed. It is not known
 * at the time of this writing whether this is a constraint itself.
 *
 * 
 * <h3>Finite elements in three dimensions</h3>
 *
 * For the interface constraints, almost the same holds as for the 2D case.
 * The numbering for the indices $n$ on the mother face is obvious and keeps
 * to the usual numbering of degrees of freedom on quadrilaterals.
 *
 * The numbering of the degrees of freedom on the interior of the refined
 * faces for the index $m$ is as follows: let $d_v$ and $d_l$ be as above,
 * and $d_q$ be the number of degrees of freedom per quadrilateral (and
 * therefore per face), then $m=0...d_v-1$ denote the dofs on the vertex at
 * the center, $m=d_v...5d_v-1$ for the dofs on the vertices at the center
 * of the bounding lines of the quadrilateral,
 * $m=5d_v..5d_v+4*d_l-1$ are for the degrees of freedom on
 * the four lines connecting the center vertex to the outer boundary of the
 * mother face, $m=5d_v+4*d_l...5d_v+4*d_l+8*d_l-1$ for the degrees of freedom
 * on the small lines surrounding the quad,
 * and $m=5d_v+12*d_l...5d_v+12*d_l+4*d_q-1$ for the dofs on the
 * four child faces. Note the direction of the lines at the boundary of the
 * quads, as shown below.
 *
 * The order of the twelve lines and the four child faces can be extracted
 * from the following sketch, where the overall order of the different
 * dof groups is depicted:
 * @verbatim
 *    *--13--3--14--*
 *    |      |      |
 *    16 20  7  19  12
 *    |      |      |
 *    4--8---0--6---2
 *    |      |      |
 *    15 17  5  18  11
 *    |      |      |
 *    *--9---1--10--*
 * @endverbatim
 * The numbering of vertices and lines, as well as the numbering of
 * children within a line is consistent with the one described in
 * Triangulation. Therefore, this numbering is seen from the
 * outside and inside, respectively, depending on the face.
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
 * freedom on both these faces. The problem is that there is now an edge
 * (the top right one of cell 2) which is part of both faces. The hanging
 * node(s) on this edge are therefore constrained twice, once from both
 * faces. To be meaningful, these constraints of course have to be
 * consistent: both faces have to constrain the hanging nodes on the edge to
 * the same nodes on the coarse edge (and only on the edge, as there can
 * then be no constraints to nodes on the rest of the face), and they have
 * to do so with the same weights. This is sometimes tricky since the nodes
 * on the edge may have different local numbers.
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
 * @author Wolfgang Bangerth, 1998, 2002; Ralf Hartmann, Guido Kanschat, 2001
 */
/*@}*/


// --------- inline and template functions ---------------

template <int dim>
inline
unsigned int 
FiniteElementData<dim>::n_dofs_per_vertex () const
{
  return dofs_per_vertex;
}


template <int dim>
inline
unsigned int 
FiniteElementData<dim>::n_dofs_per_line () const
{
  return dofs_per_line;
}


template <int dim>
inline
unsigned int 
FiniteElementData<dim>::n_dofs_per_quad () const
{
  return dofs_per_quad;
}


template <int dim>
inline
unsigned int 
FiniteElementData<dim>::n_dofs_per_hex () const
{
  return dofs_per_hex;
}


template <int dim>
inline
unsigned int 
FiniteElementData<dim>::n_dofs_per_face () const
{
  return dofs_per_face;
}


template <int dim>
inline
unsigned int 
FiniteElementData<dim>::n_dofs_per_cell () const
{
  return dofs_per_cell;
}


template <int dim>
inline
unsigned int 
FiniteElementData<dim>::n_components () const
{
  return components;
}



template <int dim>
inline
unsigned int 
FiniteElementData<dim>::tensor_degree () const
{
  return degree;
}


template <int dim>
inline
bool
FiniteElementData<dim>::conforms (Conformity space) const
{
  return ((space & conforming_space) != 0);
}




#endif
