//----------------------------  fe_base.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  fe_base.h  ---------------------------
#ifndef __deal2__fe_base_h
#define __deal2__fe_base_h

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
				      * Number of components and dimension of
				      * the image space.
				      */
    const unsigned int components;

    				     /**
				      * Default constructor. Constructs
				      * an element
				      * which is not so useful. Checking
				      * @p{dofs_per_cell} is therefore a good way to
				      * check if something went wrong. 
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
				      * Comparison operator. It is not clear to
				      * me (WB) why we have to declare and implement
				      * this one explicitely.
				      */
//TODO:[WB] (compiler) remove operator and let the compiler generate it as soon as it is willing to do so    
    bool operator == (const FiniteElementData<dim> &) const;

				     /**
				      * Exception
				      */
    DeclException2 (ExcSpaceDimensionMismatch, int, int,
		    << "used " << arg1 << "-d function for " << arg2 << "-d object");
};


//TODO:[WB] Unify FiniteElementBase with FiniteElement? There is no clear distinction between the two

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
 * into the @p{FiniteElementData<dim>} class which is explicitely
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
 * \subsection{Finite elements in one dimension}
 *
 * Finite elements in one dimension need only set the @p{restriction}
 * and @p{prolongation} matrices. The constructor of this class in one
 * dimension presets the @p{interface_constraints} matrix to have
 * dimension zero. Changing this behaviour in derived classes is
 * generally not a reasonable idea and you risk getting into trouble.
 * 
 * \subsection{Finite elements in two dimensions}
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
 * @author Wolfgang Bangerth, 1998, Guido Kanschat, 2001
 */
template <int dim>
class FiniteElementBase : public Subscriptor,
			  public FiniteElementData<dim>
{
  public:
				   /**
				    * Basis class for internal data.
				    * Adds data for second derivatives to
				    * @ref{Mapping::InternalDataBase}
				    *
				    * @author Guido Kanschat, 2001
				    */
  class InternalDataBase : public Mapping<dim>::InternalDataBase
    {
      public:
				       /**
					* Initialize some pointers
					* used in the computation of
					* second derivatives by finite
					* differencing of gradients.
					*/
      void initialize_2nd (const FiniteElement<dim> *element,
			   const Mapping<dim>       &mapping,
			   const Quadrature<dim>    &quadrature);
      
      				       /**
					* Destructor. Needed to avoid
					* memory leaks with difference
					* quotients.
					*/
      ~InternalDataBase ();

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
				      * Construct an object of this type.
				      * You have to set the
				      * matrices explicitely after calling
				      * this base class' constructor.
				      */
    FiniteElementBase (const FiniteElementData<dim> &fe_data,
		       const std::vector<bool> &restriction_is_additive_flags);

				     /**
				      * Return the value of the
				      * @p{i}th shape function at the
				      * point @p{p}.  @p{p} is a point
				      * on the reference element.
				      *
				      * An
				      * @p{ExcUnitShapeValuesDoNotExist}
				      * is thrown if the shape values
				      * of the @p{FiniteElement} under
				      * consideration depend on the
				      * shape of the cell in real
				      * space.
				      */
    virtual double shape_value (const unsigned int i,
			        const Point<dim> &p) const;
    
				     /**
				      * Return the gradient of the
				      * @p{i}th shape function at the
				      * point @p{p}. @p{p} is a point
				      * on the reference element, and
				      * likewise the gradient is the
				      * gradient on the unit cell with
				      * respect to unit cell
				      * coordinates.
				      *
				      * An
				      * @p{ExcUnitShapeValuesDoNotExist}
				      * is thrown if the shape values
				      * of the @p{FiniteElement} under
				      * consideration depend on the
				      * shape of the cell in real
				      * space.
				      */
    virtual Tensor<1,dim> shape_grad (const unsigned int  i,
				      const Point<dim>   &p) const;

				     /**
				      * Return the tensor of second
				      * derivatives of the @p{i}th
				      * shape function at point @p{p}
				      * on the unit cell. The
				      * derivatives are derivatives on
				      * the unit cell with respect to
				      * unit cell coordinates.
				      *
				      * An
				      * @p{ExcUnitShapeValuesDoNotExist}
				      * is thrown if the shape values
				      * of the @p{FiniteElement} under
				      * consideration depend on the
				      * shape of the cell in real
				      * space.
				      */
    virtual Tensor<2,dim> shape_grad_grad (const unsigned int  i,
					   const Point<dim> &p) const;

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
				      * the matrix arrays @p{restriction} and
				      * @p{prolongation}.
				      */
    bool operator == (const FiniteElementBase<dim> &) const;

				     /**
				      * Compute system index from components.
				      */
    unsigned int component_to_system_index (const unsigned int component,
					    const unsigned int component_index) const;
  
				     /**
				      * Compute component and index from
				      * system index.
				      *
				      * Return value contains first
				      * component and second index in
				      * component.
				      */
    std::pair<unsigned int,unsigned int>
    system_to_component_index (const unsigned int index) const; 
    
				     /**
				      * Compute system index from components on a face.
				      */
    unsigned int face_component_to_system_index (const unsigned int component,
						 const unsigned int component_index) const;
  
				     /**
				      * Compute component and index from system
				      * index for a face.
				      *
				      * Return value contains first
				      * component and second index in
				      * component.
				      */
    std::pair<unsigned int,unsigned int>
    face_system_to_component_index (const unsigned int index) const;
    
 				     /**
				      * The base element establishing a
				      * component.
				      *
				      * This table converts a
				      * component number to the
				      * @p{base_element} number. While
				      * component information contains
				      * multiplicity of base elements,
				      * the result allows access to
				      * shape functions of the base
				      * element.
				      */
    unsigned int component_to_base(unsigned int index) const;

				     /**
				      * Access the @p{restriction_is_additive_flag}
				      * field. See there for more information on 
				      * its contents.
				      */
    bool restriction_is_additive (const unsigned int component) const;

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
				      */
    const std::vector<Point<dim> > & get_unit_support_points () const;    

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
				      */
    const std::vector<Point<dim-1> > & get_unit_face_support_points () const;    

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
    DeclException0 (ExcUnitShapeValuesDoNotExist);

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
#if !((__GNUC__==2) && (__GNUC_MINOR__==95))
    FullMatrix<double> restriction[GeometryInfo<dim>::children_per_cell];
#else
    FullMatrix<double> restriction[1 << dim];
#endif

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
    std::vector< std::pair<unsigned int, unsigned int> > system_to_component_table;

				     /**
				      * Map between linear dofs and
				      * component dofs on face.
				      */
    std::vector< std::pair<unsigned int, unsigned int> > face_system_to_component_table;

				     /**
				      * Map between component and
				      * linear dofs.
				      */
    std::vector< std::vector<unsigned int> > component_to_system_table;

				     /**
				      * Map between component and
				      * linear dofs on a face.
				      */
    std::vector< std::vector<unsigned int> > face_component_to_system_table;
    
				     /**
				      * The base element establishing
				      * a component.
				      *
				      * This table converts a
				      * component number to the
				      * @p{base_element} number. While
				      * component information contains
				      * multiplicity of base elements,
				      * the result allows access to
				      * shape functions of the base
				      * element.
				      */
    std::vector<unsigned int> component_to_base_table;

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
				      * There is one flag per
				      * component in vector valued
				      * elements.
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
    friend class FESystem<dim>;
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
std::pair<unsigned int,unsigned int>
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
std::pair<unsigned int,unsigned int>
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
