//----------------------------  fe_base.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  fe_base.h  ---------------------------
#ifndef __deal2__fe_base_h
#define __deal2__fe_base_h

#include <base/config.h>
#include <base/exceptions.h>
#include <base/subscriptor.h>
#include <base/point.h>
#include <base/tensor.h>
#include <grid/geometry_info.h>
#include <lac/full_matrix.h>
#include <fe/fe_update_flags.h>
#include <fe/mapping.h>

template<int dim> class FESystem;



/**
 * Dimension independent data for finite elements. See the derived
 * class @ref{FiniteElementBase} class for information on its use. All
 * its data are available to the implementation in a concrete finite
 * element class.
 *
 * Remark on a change in implementation: it is now wrong to cast a
 * pointer to @ref{FiniteElement} to a pointer to
 * @p{FiniteElementData} and delete it. The virtual destructor has
 * been moved up. In a later version, @p{FiniteElementData} and
 * @ref{FiniteElementBase} should be private base classes of
 * @ref{FiniteElement}.
 *
 * @author Wolfgang Bangerth, Guido Kanschat, 1998, 1999, 2000, 2001
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
				      * writing dimension independent
				      * programs easier.
				      */
    const unsigned int dofs_per_face;
    
				     /**
				      * Total number of degrees of freedom
				      * on a cell. This information is
				      * redundant to some fields in the
				      * derived classes but makes
				      * writing dimension independent
				      * programs easier.
				      */
    const unsigned int dofs_per_cell;

				     /**
				      * Number of components and dimension of
				      * the image space.
				      */
    const unsigned int components;

    				     /**
				      * Default
				      * constructor. Constructs an
				      * element which is not so
				      * useful. Checking
				      * @p{dofs_per_cell} is therefore
				      * a good way to check if
				      * something went wrong.
				      */
    FiniteElementData ();

				     /**
				      * Constructor for
				      * all-dimensional objects. The
				      * numbers in @p{dofs_per_object}
				      * represent the numbers of DoFs
				      * of grid objects in
				      * dim-ascending order. That is,
				      * @p{dofs_per_object[0]=dofs_per_vertex},
				      * @p{dofs_per_object[1]=dofs_per_line},
				      * @p{dofs_per_object[2]=dofs_per_quad},
				      * @p{dofs_per_object[3]=dofs_per_hex}.
				      *
				      * Hence this constructor requires
				      * @p{dofs_per_object.size()==dim+1}.
				      */
    FiniteElementData (const std::vector<unsigned int> &dofs_per_object,
		       const unsigned int n_components);

				     /**
				      * Return the @p{dofs_per_vertex}.
				      */
    unsigned int n_dofs_per_vertex () const;

				     /**
				      * Return the @p{dofs_per_line}.
				      */
    unsigned int n_dofs_per_line () const;

    				     /**
				      * Return the @p{dofs_per_quad}.
				      */
    unsigned int n_dofs_per_quad () const;

    				     /**
				      * Return the @p{dofs_per_hex}.
				      */
    unsigned int n_dofs_per_hex () const;

    				     /**
				      * Return the @p{dofs_per_face}.
				      */
    unsigned int n_dofs_per_face () const;

    				     /**
				      * Return the @p{dofs_per_cell}.
				      */
    unsigned int n_dofs_per_cell () const;

    				     /**
				      * Return the @p{components}.
				      */
    unsigned int n_components () const;

				     /**
				      * Comparison operator.
				      */
    bool operator == (const FiniteElementData<dim> &) const;

				     /**
				      * Exception
				      */
    DeclException2 (ExcSpaceDimensionMismatch, int, int,
		    << "used " << arg1 << "-d function for " << arg2 << "-d object");
};




/**
 * Base class for finite elements in arbitrary dimensions. This class
 * provides several fields which describe a specific finite element
 * and which are filled by derived classes. It more or less only
 * offers the fields and access functions which makes it possible to
 * copy finite elements without knowledge of the actual type (linear,
 * quadratic, etc).
 *
 * The implementation of this base class is split into two parts:
 * those fields which are not common to all dimensions
 * (@p{dofs_per_quad} for example are only useful for @p{dim>=2}) are put
 * into the @p{FiniteElementData<dim>} class which is explicitly
 * specialized for all used dimensions, while those fields which may
 * be formulated in a dimension-independent way are put into the
 * present class.
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
 * @sect3{Support points}
 *
 * Since a @ref{FiniteElement} does not have information on the actual
 * grid cell, it can only provide support points on the unit
 * cell. Support points on the actual grid cell must be computed by
 * mapping these points. The class used for this kind of operation is
 * @ref{FEValues}. In most cases, code of the following type will
 * serve to provide the mapped support points.
 *
 * @begin{verbatim}
 * Quadrature<dim> dummy_quadrature (fe.get_unit_support_points());
 * FEValues<dim>   fe_values (mapping, fe, dummy_quadrature,
 *                            update_q_points);
 * fe_values.reinit (cell);
 * Point<dim>& mapped_point = fe_values.quadrature_point (i);
 * @end{verbatim}
 *
 * Alternatively, the points can be transformed one-by-one:
 * @begin{verbatim}
 * const vector<Point<dim> >& unit_points =
 *    fe.get_unit_support_points();
 *
 * Point<dim> mapped_point =
 *    mapping.transform_unit_to_real_cell (cell, unit_points[i]);
 * @end{verbatim}
 * This is a shortcut, and as all shortcuts should be used cautiously.
 * If the mapping of all support points is needed, the first variant should
 * be preferred for efficiency.
 *
 * @sect3{Finite elements in one dimension}
 *
 * Finite elements in one dimension need only set the @p{restriction}
 * and @p{prolongation} matrices. The constructor of this class in one
 * dimension presets the @p{interface_constraints} matrix to have
 * dimension zero. Changing this behaviour in derived classes is
 * generally not a reasonable idea and you risk getting into trouble.
 * 
 * @sect3{Finite elements in two dimensions}
 * 
 * In addition to the fields already present in 1D, a constraint
 * matrix is needed, if the finite element has node values located on
 * edges or vertices. These constraints are represented by a $m\times
 * n$-matrix @p{interface_constraints}, where $n$ is the number of
 * degrees of freedom on the refined side without the corner vertices
 * (those dofs on the middle vertex plus those on the two lines), and
 * $m$ is that of the unrefined side (those dofs on the two vertices
 * plus those on the line). The matrix is thus a rectangular one.
 *
 * The mapping of the dofs onto the indices of the matrix on the
 * unrefined side is as follows: let $d_v$ be the number of dofs on a
 * vertex, $d_l$ that on a line, then $m=0...d_v-1$ refers to the dofs
 * on vertex zero of the unrefined line, $m=d_v...2d_v-1$ to those on
 * vertex one, $m=2d_v...2d_v+d_l-1$ to those on the line.
 *
 * Similarly, $n=0...d_v-1$ refers to the dofs on the middle vertex of
 * the refined side (vertex one of child line zero, vertex zero of
 * child line one), $n=d_v...d_v+d_l-1$ refers to the dofs on child
 * line zero, $n=d_v+d_l...d_v+2d_l-1$ refers to the dofs on child
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
 * @sect3{Finite elements in three dimensions}
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
 * The numbering of vertices and lines, as well as the numbering of
 * children within a line is consistent with the one described in
 * @ref{Triangulation}. Therefore, this numbering is seen from the
 * outside and inside, respectively, depending on the face.
 *
 * If of the cells adjacent to one line more than one is refined and
 * there is at least one unrefined cell, then the degrees of freedom
 * on the refined line are constrained from two cells. For example,
 * consider the cell behind the face shown above is refined, while the
 * one in front of the face is not refined; then the dofs on the lines
 * numbered 9 and 10 are constrained. If there are two more cells
 * below the ones just introduced, with a common face right below the
 * one shown, and of these is one refined and one unrefined one, then
 * the degrees on the two mentioned small lines are constrained a
 * second time. Since these constraints must be unique, it follows
 * that the constraints for the degrees of freedom on refined lines
 * may only be in terms of the degrees of freedom on the unrefined
 * line, not in terms of the other degrees of freedom on a face.
 *
 * Since the handling of constraints on degrees of freedom is mostly done
 * by the @p{ConstraintMatrix} class, this class checks whether the constraints
 * introduced from the two sides are unique; it is able to handle the fact
 * that the constraints for some of the dofs are entered more than once.
 *
 * @author Wolfgang Bangerth, 1998, 2002, Ralf Hartmann, Guido Kanschat, 2001
 */
template <int dim>
class FiniteElementBase : public Subscriptor,
			  public FiniteElementData<dim>
{
  public:
				   /**
				    * Base class for internal data.
				    * Adds data for second derivatives to
				    * @ref{Mapping::InternalDataBase}
				    *
				    * For information about the
				    * general purpose of this class,
				    * see the documentation of the
				    * base class.
				    *
				    * @author Guido Kanschat, 2001
				    */
  class InternalDataBase : public Mapping<dim>::InternalDataBase
    {
      public:      
					 /**
					  * Destructor. Needed to
					  * avoid memory leaks with
					  * difference quotients.
					  */
	virtual ~InternalDataBase ();

					 /**
					  * Initialize some pointers
					  * used in the computation of
					  * second derivatives by
					  * finite differencing of
					  * gradients.
					  */
	void initialize_2nd (const FiniteElement<dim> *element,
			     const Mapping<dim>       &mapping,
			     const Quadrature<dim>    &quadrature);
	
					 /**
					  * Storage for @p{FEValues}
					  * objects needed to
					  * approximate second
					  * derivatives.
					  *
					  * The ordering is @p{p+hx},
					  * @p{p+hy}, @p{p+hz},
					  * @p{p-hx}, @p{p-hy},
					  * @p{p-hz}, where unused
					  * entries in lower dimensions
					  * are missing.
					  */
	std::vector<FEValues<dim>*> differences;
    };
  
				     /**
				      * Construct an object of this
				      * type. You have to set some
				      * member variables, for example
				      * some matrices, explicitly
				      * after calling this base class'
				      * constructor. For this see the
				      * existing finite element
				      * classes. For the second and
				      * third parameter of this
				      * constructor, see the document
				      * of the respective member
				      * variables.
				      */
    FiniteElementBase (const FiniteElementData<dim> &fe_data,
		       const std::vector<bool> &restriction_is_additive_flags,
		       const std::vector<std::vector<bool> > &nonzero_components);

				     /**
				      * Return the value of the
				      * @p{i}th shape function at the
				      * point @p{p}. @p{p} is a point
				      * on the reference element. If
				      * the finite element is
				      * vector-valued, then return the
				      * value of the only non-zero
				      * component of the vector value
				      * of this shape function. If the
				      * shape function has more than
				      * one non-zero component (which
				      * we refer to with the term
				      * non-primitive), then derived
				      * classes implementing this
				      * function should throw an
				      * exception of type
				      * @p{ExcShapeFunctionNotPrimitive}. In
				      * that case, use the
				      * @ref{shape_value_component}
				      * function.
				      *
				      * An
				      * @p{ExcUnitShapeValuesDoNotExist}
				      * is thrown if the shape values
				      * of the @p{FiniteElement} under
				      * consideration depends on the
				      * shape of the cell in real
				      * space.
				      */
    virtual double shape_value (const unsigned int  i,
			        const Point<dim>   &p) const;

				     /**
				      * Just like for @p{shape_value},
				      * but this function will be
				      * called when the shape function
				      * has more than one non-zero
				      * vector component. In that
				      * case, this function should
				      * return the value of the
				      * @p{component}-th vector
				      * component of the @p{i}th shape
				      * function at point @p{p}.
				      */
    virtual double shape_value_component (const unsigned int i,
					  const Point<dim>   &p,
					  const unsigned int component) const;
    
				     /**
				      * Return the gradient of the
				      * @p{i}th shape function at the
				      * point @p{p}. @p{p} is a point
				      * on the reference element, and
				      * likewise the gradient is the
				      * gradient on the unit cell with
				      * respect to unit cell
				      * coordinates. If
				      * the finite element is
				      * vector-valued, then return the
				      * value of the only non-zero
				      * component of the vector value
				      * of this shape function. If the
				      * shape function has more than
				      * one non-zero component (which
				      * we refer to with the term
				      * non-primitive), then derived
				      * classes implementing this
				      * function should throw an
				      * exception of type
				      * @p{ExcShapeFunctionNotPrimitive}. In
				      * that case, use the
				      * @ref{shape_grad_component}
				      * function.
				      *
				      * An
				      * @p{ExcUnitShapeValuesDoNotExist}
				      * is thrown if the shape values
				      * of the @p{FiniteElement} under
				      * consideration depends on the
				      * shape of the cell in real
				      * space.
				      */
    virtual Tensor<1,dim> shape_grad (const unsigned int  i,
				      const Point<dim>   &p) const;

				     /**
				      * Just like for @p{shape_grad},
				      * but this function will be
				      * called when the shape function
				      * has more than one non-zero
				      * vector component. In that
				      * case, this function should
				      * return the gradient of the
				      * @p{component}-th vector
				      * component of the @p{i}th shape
				      * function at point @p{p}.
				      */
    virtual Tensor<1,dim> shape_grad_component (const unsigned int i,
						const Point<dim>   &p,
						const unsigned int component) const;

				     /**
				      * Return the tensor of second
				      * derivatives of the @p{i}th
				      * shape function at point @p{p}
				      * on the unit cell. The
				      * derivatives are derivatives on
				      * the unit cell with respect to
				      * unit cell coordinates. If
				      * the finite element is
				      * vector-valued, then return the
				      * value of the only non-zero
				      * component of the vector value
				      * of this shape function. If the
				      * shape function has more than
				      * one non-zero component (which
				      * we refer to with the term
				      * non-primitive), then derived
				      * classes implementing this
				      * function should throw an
				      * exception of type
				      * @p{ExcShapeFunctionNotPrimitive}. In
				      * that case, use the
				      * @ref{shape_grad_grad_component}
				      * function.
				      *
				      * An
				      * @p{ExcUnitShapeValuesDoNotExist}
				      * is thrown if the shape values
				      * of the @p{FiniteElement} under
				      * consideration depends on the
				      * shape of the cell in real
				      * space.
				      */
    virtual Tensor<2,dim> shape_grad_grad (const unsigned int  i,
					   const Point<dim>   &p) const;

				     /**
				      * Just like for @p{shape_grad_grad},
				      * but this function will be
				      * called when the shape function
				      * has more than one non-zero
				      * vector component. In that
				      * case, this function should
				      * return the gradient of the
				      * @p{component}-th vector
				      * component of the @p{i}th shape
				      * function at point @p{p}.
				      */
    virtual Tensor<2,dim> shape_grad_grad_component (const unsigned int i,
						     const Point<dim>   &p,
						     const unsigned int component) const;

				     /**
				      * Projection from a fine grid
				      * space onto a coarse grid
				      * space. If this projection
				      * operator is associated with a
				      * matrix @p{P}, then the
				      * restriction of this matrix
				      * @p{P_i} to a single child cell
				      * is returned here.
				      *
				      * The matrix @p{P} is the
				      * concatenation or the sum of
				      * the cell matrices @p{P_i},
				      * depending on the
				      * @p{restriction_is_additive_flags}
				      * given to the constructor. This
				      * distinguishes interpolation
				      * (concatenation) and projection
				      * with respect to scalar
				      * products (summation).
				      *
				      * Row and column indices are
				      * related to coarse grid and
				      * fine grid spaces,
				      * respectively, consistent with
				      * the definition of the
				      * associated operator.
				      *
				      * If projection matrices are not
				      * implemented in the derived
				      * finite element class, this
				      * function aborts with
				      * @p{ExcProjectionVoid}.
				      */
    const FullMatrix<double> & restrict (const unsigned int child) const;

				     /**
				      * Embedding matrix between grids.
				      * 
				      * The identity operator from a
				      * coarse grid space into a fine
				      * grid space is associated with
				      * a matrix @p{P}. The
				      * restriction of this matrix @p{P_i} to
				      * a single child cell is
				      * returned here.
				      *
				      * The matrix @p{P} is the
				      * concatenation, not the sum of
				      * the cell matrices
				      * @p{P_i}. That is, if the same
				      * non-zero entry @p{j,k} exists
				      * in in two different child
				      * matrices @p{P_i}, the value
				      * should be the same in both
				      * matrices and it is copied into
				      * the matrix @p{P} only once.
				      *
				      * Row and column indices are
				      * related to fine grid and
				      * coarse grid spaces,
				      * respectively, consistent with
				      * the definition of the
				      * associated operator.
				      *
				      * These matrices are used by
				      * routines assembling the
				      * prolongation matrix for
				      * multi-level methods.  Upon
				      * assembling the transfer matrix
				      * between cells using this
				      * matrix array, zero elements in
				      * the prolongation matrix are
				      * discarded and will not fill up
				      * the transfer matrix.
				      *
				      * If projection matrices are not
				      * implemented in the derived
				      * finite element class, this
				      * function aborts with
				      * @p{ExcEmbeddingVoid}.
				      */
    const FullMatrix<double> & prolongate (const unsigned int child) const;

				     /**
				      * Return a readonly reference to
				      * the matrix which describes the
				      * constraints at the interface
				      * between a refined and an
				      * unrefined cell.
				      *
				      * The matrix is obviously empty
				      * in only one space dimension,
				      * since there are no constraints
				      * then.
				      */
    const FullMatrix<double> & constraints () const;

				     /**
				      * Comparison operator. We also
				      * check for equality of the
				      * constraint matrix, which is
				      * quite an expensive operation.
				      * Do therefore use this function
				      * with care, if possible only
				      * for debugging purposes.
				      *
				      * Since this function is not
				      * that important, we avoid an
				      * implementational question
				      * about comparing arrays and do
				      * not compare the matrix arrays
				      * @p{restriction} and
				      * @p{prolongation}.
				      */
    bool operator == (const FiniteElementBase<dim> &) const;

				     /**
				      * Given a vector component and
				      * an index of a shape function
				      * within the shape functions
				      * corresponding to this vector
				      * component, return the index of
				      * this shape function within the
				      * shape functions of this
				      * element. If this is a scalar
				      * element, then the given
				      * component may only be zero,
				      * and the given component index
				      * is also the return value.
				      *
				      * If the finite element is
				      * vector-valued and has
				      * non-primitive shape functions,
				      * i.e. some of its shape
				      * functions are non-zero in more
				      * than just one vector
				      * component, then this function
				      * cannot be used since shape
				      * functions are no more
				      * associated with individual
				      * vector components, and an
				      * exception of type
				      * @p{ExcFENotPrimitive} is
				      * thrown.
				      */
    unsigned int component_to_system_index (const unsigned int component,
					    const unsigned int component_index) const;
  
				     /**
				      * Same as above, but compute the
				      * data from the index of a shape
				      * function on a face.
				      */
    unsigned int face_component_to_system_index (const unsigned int component,
						 const unsigned int component_index) const;

				     /**
				      * Compute vector component and
				      * index of this shape function
				      * within the shape functions
				      * corresponding to this
				      * component from the index of a
				      * shape function within this
				      * finite element.
				      *
				      * If the element is scalar, then
				      * the component is always zero,
				      * and the index within this
				      * component is equal to the
				      * overall index.
				      *
				      * If the shape function
				      * referenced has more than one
				      * non-zero component, then it
				      * cannot be associated with one
				      * vector component, and an
				      * exception of type
				      * @p{ExcShapeFunctionNotPrimitive}
				      * will be raised.
				      *
				      * Note that if the element is
				      * composed of other (base)
				      * elements, and a base element
				      * has more than one component
				      * but all its shape functions
				      * are primitive (i.e. are
				      * non-zero in only one
				      * component), then this mapping
				      * contains valid
				      * information. However, the
				      * index of a shape function of
				      * this element within one
				      * component (i.e. the second
				      * number of the respective entry
				      * of this array) does not
				      * indicate the index of the
				      * respective shape function
				      * within the base element (since
				      * that has more than one
				      * vector-component). For this
				      * information, refer to the
				      * @p{system_to_base_table}
				      * field.
				      */
    std::pair<unsigned int,unsigned int>
    system_to_component_index (const unsigned int index) const;    
  
				     /**
				      * Same as above, but do it for
				      * shape functions and their
				      * indices on a face.
				      */
    std::pair<unsigned int,unsigned int>
    face_system_to_component_index (const unsigned int index) const;
    
 				     /**
				      * Given a vector component,
				      * return an index which base
				      * element implements this
				      * component, and which vector
				      * component is this base element
				      * this is. This information is
				      * only of interest for
				      * vector-valued finite elements
				      * which are composed of several
				      * sub-elements. In that case,
				      * one may want to obtain
				      * information about the element
				      * implementing a certain vector
				      * component, which can be done
				      * using this function and the
				      * @ref{FESystem}::@p{base_element}
				      * function.
				      *
				      * If this is a scalar finite
				      * element, then the return value
				      * is always equalt to zero.
				      */
    std::pair<unsigned int,unsigned int>
    component_to_base (unsigned int component) const;

				     /**
				      * Access the
				      * @p{restriction_is_additive_flag}
				      * field. See there for more
				      * information on its contents.
				      *
				      * The index must be between zero
				      * and the number of shape
				      * functions of this element.
				      */
    bool restriction_is_additive (const unsigned int index) const;

				     /**
				      * Return the support points of
				      * the trial functions on the
				      * unit cell, if the derived
				      * finite element defines some.
				      * Finite elements that allow
				      * some kind of interpolation
				      * operation usually have support
				      * points. On the other hand,
				      * elements that define their
				      * degrees of freedom by, for
				      * example, moments on faces, or
				      * as derivatives, don't have
				      * support points. In that case,
				      * the returned field is empty.
				      *
				      * If the finite element defines
				      * support points, then their
				      * number equals the number of
				      * degrees of freedom of the
				      * element.  The order of points
				      * in the array matches that
				      * returned by the
				      * @p{cell->get_dof_indices}
				      * function.
				      *
				      * See the class documentation
				      * for details on support points.
				      */
    const std::vector<Point<dim> > &
    get_unit_support_points () const;    

				     /**
				      * Return whether a finite
				      * element has defined support
				      * points. If the result is true,
				      * then a call to the
				      * @p{get_unit_support_points}
				      * yields a non-empty array.
				      */
    bool has_support_points () const;
    
				     /**
				      * Return the support points of
				      * the trial functions on the
				      * unit face, if the derived
				      * finite element defines some.
				      * Finite elements that allow
				      * some kind of interpolation
				      * operation usually have support
				      * points. On the other hand,
				      * elements that define their
				      * degrees of freedom by, for
				      * example, moments on faces, or
				      * as derivatives, don't have
				      * support points. In that case,
				      * the returned field is empty
				      *
				      * Note that elements that have
				      * support points need not
				      * necessarily have some on the
				      * faces, even if the
				      * interpolation points are
				      * located physically on a
				      * face. For example, the
				      * discontinuous elements have
				      * interpolation points on the
				      * vertices, and for higher
				      * degree elements also on the
				      * faces, but they are not
				      * defined to be on faces since
				      * in that case degrees of
				      * freedom from both sides of a
				      * face (or from all adjacent
				      * elements to a vertex) would be
				      * identified with each other,
				      * which is not what we would
				      * like to have). Logically,
				      * these degrees of freedom are
				      * therefore defined to belong to
				      * the cell, rather than the face
				      * or vertex. In that case, the
				      * returned element would
				      * therefore have length zero.
				      *
				      * If the finite element defines
				      * support points, then their
				      * number equals the number of
				      * degrees of freedom on the face
				      * (@p{dofs_per_face}). The order
				      * of points in the array matches
				      * that returned by the
				      * @p{cell->get_dof_indices}
				      * function.
				      *
				      * See the class documentation
				      * for details on support points.
				      */
    const std::vector<Point<dim-1> > &
    get_unit_face_support_points () const;    

				     /**
				      * Return whether a finite
				      * element has defined support
				      * points on faces. If the result
				      * is true, then a call to the
				      * @p{get_unit_support_points}
				      * yields a non-empty array.
				      */
    bool has_face_support_points () const;

				     /**
				      * Return in which of the vector
				      * components of this finite
				      * element the @p{i}the shape
				      * function is non-zero. The
				      * length of the returned array
				      * is equal to the number of
				      * vector components of this
				      * element.
				      *
				      * For most finite element
				      * spaces, the result of this
				      * function will be a vector with
				      * exactly one element being
				      * @p{true}, since for most
				      * spaces the individual vector
				      * components are
				      * independent. Only for those
				      * spaces that couple the
				      * components, for example to
				      * make a shape function
				      * divergence free, will there be
				      * more than one @p{true} entry.
				      */
    const std::vector<bool> &
    get_nonzero_components (const unsigned int i) const;

				     /**
				      * Return in how many vector
				      * components the @p{i}th shape
				      * function is non-zero. This
				      * value equals the number of
				      * entries equal to @p{true} in
				      * the result of the
				      * @p{get_nonzero_components}
				      * function.
				      *
				      * For most finite element
				      * spaces, the result will be
				      * equal to one. It is not equal
				      * to one only for those ansatz
				      * spaces for which vector-valued
				      * shape functions couple the
				      * individual components, for
				      * example in order to make them
				      * divergence-free.
				      */
    unsigned int
    n_nonzero_components (const unsigned int i) const;

				     /**
				      * Return whether the @p{i}th
				      * shape function is primitive in
				      * the sense that the shape
				      * function is non-zero in only
				      * one vector
				      * component. Non-primitive shape
				      * functions would then, for
				      * example, be those of
				      * divergence free ansatz spaces,
				      * in which the individual vector
				      * components are coupled.
				      *
				      * The result of the function is
				      * @p{true} if and only if the
				      * result of
				      * @p{n_nonzero_components(i)} is
				      * equal to one.
				      */
    bool
    is_primitive (const unsigned int i) const;

				     /**
				      * Return whether the entire
				      * finite element is primitive,
				      * in the sense that all its
				      * shape functions are
				      * primitive. If the finite
				      * element is scalar, then this
				      * is always the case.
				      *
				      * Since this is an extremely
				      * common operations, the result
				      * is cached in the
				      * @p{cached_primitivity}
				      * variable which is computed in
				      * the constructor.
				      */
    bool
    is_primitive () const;
    
				     /**
				      * Determine an estimate for the
				      * memory consumption (in bytes)
				      * of this object.
				      *
				      * This function is not
				      * virtual. Use a
				      * @ref{FiniteElement} object to
				      * get the actual size of a
				      * concrete element.
				      */
    unsigned int memory_consumption () const;

				     /**
				      * Exception
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
				      */
    DeclException0 (ExcFENotPrimitive);
				     /**
				      * Exception
				      */
    DeclException0 (ExcUnitShapeValuesDoNotExist);

				     /**
				      * Exception
				      */
    DeclException0 (ExcFEHasNoSupportPoints);

				     /**
				      * Exception
				      */
    DeclException0 (ExcEmbeddingVoid);
    
				     /**
				      * Exception
				      */
    DeclException0 (ExcProjectionVoid);
    
				     /**
				      * Exception
				      */
    DeclException0 (ExcConstraintsVoid);
    
				     /**
				      * Exception
				      */
    DeclException2 (ExcWrongInterfaceMatrixSize,
		    int, int,
		    << "The interface matrix has a size of " << arg1
		    << "x" << arg2
		    << ", which is not reasonable in the present dimension.");
    
  protected:  
 				     /**
				      * Array of projection
				      * matrices. See @p{restrict()}
				      * above.
				      *
				      * Matrices in this array are
				      * automatically initialized to
				      * correct size. If the derived
				      * finite element class does not
				      * implement these matrices, they
				      * should be resized to zero
				      * size.
				      */
    FullMatrix<double> restriction[GeometryInfo<dim>::children_per_cell];

    				     /**
				      * Array of embedding
				      * matrices. See @p{prolongate()}
				      * above.
				      *
				      * Matrices in this array are
				      * automatically initialized to
				      * correct size. If the derived
				      * finite element class does not
				      * implement these matrices, they
				      * should be resized to zero
				      * size.
				      */
    FullMatrix<double> prolongation[GeometryInfo<dim>::children_per_cell];

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
				      * Store what
				      * @p{system_to_component_index}
				      * will return.
				      */
    std::vector< std::pair<unsigned int, unsigned int> > system_to_component_table;

				     /**
  				      * Map between linear dofs and
 				      * component dofs on face. This
 				      * is filled with default values
 				      * in the constructor, but
 				      * derived classes will have to
 				      * overwrite the information if
 				      * necessary.
 				      *
 				      * By component, we mean the
 				      * vector component, not the base
 				      * element. The information thus
 				      * makes only sense if a shape
 				      * function is non-zero in only
 				      * one component.
				      */
    std::vector< std::pair<unsigned int, unsigned int> > face_system_to_component_table;

				     /**
				      * For each shape function, store
				      * to which base element and
				      * which instance of this base
				      * element (in case its
				      * multiplicity is greater than
				      * one) it belongs, and its index
				      * within this base element. If
				      * the element is not composed of
				      * others, then base and instance
				      * are always zero, and the index
				      * is equal to the number of the
				      * shape function. If the element
				      * is composed of single
				      * instances of other elements
				      * (i.e. all with multiplicity
				      * one) all of which are scalar,
				      * then base values and dof
				      * indices within this element
				      * are equal to the
				      * @p{system_to_component_table}. It
				      * differs only in case the
				      * element is composed of other
				      * elements and at least one of
				      * them is vector-valued itself.
				      *
				      * This array has valid values
				      * also in the case of
				      * vector-valued
				      * (i.e. non-primitive) shape
				      * functions, in contrast to the
				      * @p{system_to_component_table}.
				      */
    std::vector<std::pair<std::pair<unsigned int,unsigned int>,unsigned int> >
    system_to_base_table;

				     /**
				      * Likewise for the indices on
				      * faces.
				      */
    std::vector<std::pair<std::pair<unsigned int,unsigned int>,unsigned int> >
    face_system_to_base_table;
    
				     /**
				      * Map between component and
				      * linear dofs: For each pair of
				      * vector component and index
				      * within this component, store
				      * the global dof number in the
				      * composed element. If the
				      * element is scalar, then the
				      * outer (component) index can
				      * only be zero, and the inner
				      * index is equal to the stored
				      * value.
				      *
				      * If the element is not
				      * primitive, i.e. there are
				      * shape functions that are
				      * non-zero in more than one
				      * vector-component, then this
				      * function is obviously useless,
				      * and all entries will be
				      * invalid.
				      */
    std::vector< std::vector<unsigned int> > component_to_system_table;

				     /**
				      * Map between component and
				      * linear dofs on a face. Same
				      * applies as above.
				      */
    std::vector< std::vector<unsigned int> > face_component_to_system_table;
    
				     /**
				      * The base element establishing
				      * a component.
				      *
				      * This table converts a
				      * component number to a pair
				      * consisting of the
				      * @p{base_element} number, and
				      * the component within this base
				      * element. While component
				      * information contains
				      * multiplicity of base elements,
				      * the result allows access to
				      * shape functions of the base
				      * element.
				      *
				      * This variable is set to the
				      * correct size by the
				      * constructor of this class, but
				      * needs to be initialized by
				      * derived classes, unless its
				      * size is one and the only entry
				      * is a zero, which is the case
				      * for scalar elements. In that
				      * case, the initialization by
				      * the base class is sufficient.
				      */
    std::vector<std::pair<unsigned int, unsigned int> > component_to_base_table;
    
				     /**
				      * Projection matrices are
				      * concatenated or summed up.
				      *
				      * This flags decides on how the
				      * projection matrices of the
				      * children of the same father
				      * are put together to one
				      * operator. The possible modes
				      * are concatenation and
				      * summation.
				      *
				      * If the projection is defined
				      * by an interpolation operator,
				      * the child matrices are
				      * concatenated, i.e. values
				      * belonging to the same node
				      * functional are identified and
				      * enter the interpolated value
				      * only once. In this case, the
				      * flag must be @p{false}.
				      *
				      * For projections with respect
				      * to scalar products, the child
				      * matrices must be summed up to
				      * build the complete matrix. The
				      * flag should be @p{true}.
				      *
				      * For examples of use of these
				      * flags, see the places in the
				      * library where it is queried.
				      * 
				      * There is one flag per shape
				      * function, indicating whether
				      * it belongs to the class of
				      * shape functions that are
				      * additive in the restriction or
				      * not.
				      *
				      * Note that in previous versions
				      * of the library, there was one
				      * flag per vector component of
				      * the element. This is based on
				      * the fact that all the shape
				      * functions that belong to the
				      * same vector component must
				      * necessarily behave in the same
				      * way, to make things
				      * reasonable. However, the
				      * problem is that it is
				      * sometimes impossible to query
				      * this flag in the vector-valued
				      * case: this used to be done
				      * with the
				      * @p{system_to_component_index}
				      * function that returns which
				      * vector component a shape
				      * function is associated
				      * with. The point is that since
				      * we now support shape functions
				      * that are associated with more
				      * than one vector component (for
				      * example the shape functions of
				      * Raviart-Thomas, or Nedelec
				      * elements), that function can
				      * no more be used, so it can be
				      * difficult to find out which
				      * for vector component we would
				      * like to query the
				      * restriction-is-additive flags.
				      */
    const std::vector<bool> restriction_is_additive_flags;

				     /**
				      * List of support points on the
				      * unit cell, in case the finite
				      * element has any. The
				      * constructor leaves this field
				      * empty, derived classes may
				      * write in some contents.
				      *
				      * Finite elements that allow
				      * some kind of interpolation
				      * operation usually have support
				      * points. On the other hand,
				      * elements that define their
				      * degrees of freedom by, for
				      * example, moments on faces, or
				      * as derivatives, don't have
				      * support points. In that case,
				      * this field remains empty.
				      */
    std::vector<Point<dim> > unit_support_points;

				     /**
				      * Same for the faces. See the
				      * description of the
				      * @p{get_unit_face_support_points}
				      * function for a discussion of
				      * what contributes a face
				      * support point.
				      */
    std::vector<Point<dim-1> > unit_face_support_points;

				     /**
				      * For each shape function, give
				      * a vector of bools (with size
				      * equal to the number of vector
				      * components which this finite
				      * element has) indicating in
				      * which component each of these
				      * shape functions is non-zero.
				      *
				      * For primitive elements, there
				      * is only one non-zero
				      * component.
				      */
    const std::vector<std::vector<bool> > nonzero_components;

				     /**
				      * This array holds how many
				      * values in the respective entry
				      * of the @p{nonzero_components}
				      * element are non-zero. The
				      * array is thus a short-cut to
				      * allow faster access to this
				      * information than if we had to
				      * count the non-zero entries
				      * upon each request for this
				      * information. The field is
				      * initialized in the constructor
				      * of this class.
				      */
    const std::vector<unsigned int> n_nonzero_components_table;

				     /**
				      * Store whether all shape
				      * functions are primitive. Since
				      * finding this out is a very
				      * common operation, we cache the
				      * result, i.e. compute the value
				      * in the constructor for simpler
				      * access.
				      */
    const bool cached_primitivity;

                                     /**
				      * Compute second derivatives by
				      * finite differences of
				      * gradients.
				      */
    void compute_2nd (const Mapping<dim>                      &mapping,
		      const typename DoFHandler<dim>::cell_iterator    &cell,
		      const unsigned int                       offset,
		      typename Mapping<dim>::InternalDataBase &mapping_internal,
		      InternalDataBase                        &fe_internal,
		      FEValuesData<dim>                       &data) const;

  private:
				     /**
				      * Second derivatives of shapes
				      * functions are not computed
				      * analytically, but by finite
				      * differences of the
				      * gradients. This static
				      * variable denotes the step
				      * length to be used for
				      * that. It's value is set to
				      * 1e-6.
				      */
    static const double fd_step_length;

				     /**
				      * Given the pattern of nonzero
				      * components for each shape
				      * function, compute for each
				      * entry how many components are
				      * non-zero for each shape
				      * function. This function is
				      * used in the constructor of
				      * this class.
				      */
    static
    std::vector<unsigned int>
    compute_n_nonzero_components (const std::vector<std::vector<bool> > &nonzero_components);
    
				     /**
				      * Allow the FESystem class to
				      * access the restriction and
				      * prolongation matrices
				      * directly. Hence, FESystem has
				      * the possibility to see if
				      * these matrices are initialized
				      * without accessing these
				      * matrices through the
				      * @p{restrict} and
				      * @p{prolongate} functions. This
				      * is important as these
				      * functions include assertions
				      * that throw if the matrices are
				      * not already initialized.
				      */
    template <int dim_> friend class FESystem;
};


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
FiniteElementBase<dim>::component_to_system_index (const unsigned int component,
						   const unsigned int component_index) const
{
  Assert(component<this->n_components(),
	 ExcIndexRange(component, 0, this->n_components()));
  Assert(component_index<component_to_system_table[component].size(),
	 ExcIndexRange(component_index, 0,
		       component_to_system_table[component].size()));
  Assert (is_primitive(),
	  typename FiniteElementBase<dim>::ExcFENotPrimitive());
  
  return component_to_system_table[component][component_index];
}


template <int dim>  
inline
std::pair<unsigned int,unsigned int>
FiniteElementBase<dim>::system_to_component_index (const unsigned int index) const
{
  Assert (index < system_to_component_table.size(),
	 ExcIndexRange(index, 0, system_to_component_table.size()));
  Assert (is_primitive (index),
	  typename FiniteElementBase<dim>::ExcShapeFunctionNotPrimitive(index));
  return system_to_component_table[index];
}


template <int dim>
inline
unsigned int
FiniteElementBase<dim>::face_component_to_system_index (const unsigned int component,
							const unsigned int component_index) const
{
  Assert(component<this->n_components(),
	 ExcIndexRange(component, 0, this->n_components()));
  Assert(component_index<face_component_to_system_table[component].size(),
	 ExcIndexRange(component_index, 0,
		       face_component_to_system_table[component].size()));
  Assert (is_primitive(),
	  typename FiniteElementBase<dim>::ExcFENotPrimitive());

  return face_component_to_system_table[component][component_index];
}


template <int dim>  
inline
std::pair<unsigned int,unsigned int>
FiniteElementBase<dim>::face_system_to_component_index (const unsigned int index) const
{
  Assert(index < face_system_to_component_table.size(),
	 ExcIndexRange(index, 0, face_system_to_component_table.size()));
//TODO: check for primitivity of this shape function. this needs the global dof index
//    Assert (is_primitive (face_to_cell_index(index)),
//  	  typename FiniteElementBase<dim>::ExcShapeFunctionNotPrimitive(index));
  return face_system_to_component_table[index];
}


template <int dim>  
inline
std::pair<unsigned int,unsigned int>
FiniteElementBase<dim>::component_to_base (unsigned int index) const
{
  Assert(index < component_to_base_table.size(),
	 ExcIndexRange(index, 0, component_to_base_table.size()));

  return component_to_base_table[index];
}


template <int dim>
inline
bool
FiniteElementBase<dim>::restriction_is_additive (const unsigned int index) const
{
  Assert(index < this->dofs_per_cell,
	 ExcIndexRange(index, 0, this->dofs_per_cell));
  return restriction_is_additive_flags[index];
}


template <int dim>
inline
const std::vector<bool> &
FiniteElementBase<dim>::get_nonzero_components (const unsigned int i) const
{
  Assert (i < this->dofs_per_cell, ExcIndexRange (i, 0, this->dofs_per_cell));
  return nonzero_components[i];
};



template <int dim>
inline
unsigned int
FiniteElementBase<dim>::n_nonzero_components (const unsigned int i) const
{
  Assert (i < this->dofs_per_cell, ExcIndexRange (i, 0, this->dofs_per_cell));
  return n_nonzero_components_table[i];
};



template <int dim>
inline
bool
FiniteElementBase<dim>::is_primitive (const unsigned int i) const
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
  return (n_nonzero_components_table[i] == 1);
};


template <int dim>
inline
bool
FiniteElementBase<dim>::is_primitive () const
{
  return cached_primitivity;
};




#endif
