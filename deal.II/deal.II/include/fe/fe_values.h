//----------------------------  fe_values.h  ---------------------------
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
//----------------------------  fe_values.h  ---------------------------
#ifndef __deal2__fe_values_h
#define __deal2__fe_values_h


#include <base/exceptions.h>
#include <base/subscriptor.h>
#include <base/point.h>
#include <base/tensor.h>
#include <base/quadrature.h>
#include <lac/full_matrix.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <fe/fe_update_flags.h>
#include <fe/mapping.h>

template <int dim> class Quadrature;


/**
 * Contains all data vectors for @p{FEValues}.
 *
 * This class has been extracted from @p{FEValuesBase} to be handed
 * over to the fill functions of @p{Mapping} and
 * @p{FiniteElement}. All data fields are public, but this is not
 * critical, because access to this object is private in @p{FEValues}.
 *
 * @author Guido Kanschat, 2000
 */
template <int dim>
class FEValuesData
{
  public:
				     /**
				      * Initialize all vectors to correct size.
				      */
    void initialize (const unsigned int n_quadrature_points,
		     const unsigned int n_shapes,
		     const UpdateFlags  flags);

				     /**
				      * Storage type for shape values.
				      */
    typedef FullMatrix<double> ShapeVector;

				     /**
				      * Storage type for gradients.
				      */
    typedef typename std::vector<typename std::vector<Tensor<1,dim> > > GradientVector;
    
				     /**
				      * Store the values of the shape functions
				      * at the quadrature points. Rows in the
				      * matrices denote the values of a single
				      * shape function at the different points,
				      * columns are for a single point with the
				      * different shape functions.
				      */
    ShapeVector shape_values;

    				     /**
				      * Store the gradients of the shape
				      * functions at the quadrature points.
				      * For each shape function, there is a vector
				      * of the gradients in each quadrature point.
				      */
    GradientVector shape_gradients;

				     /**
				      * Store the 2nd derivatives of the shape
				      * functions at the quadrature points.
				      * For each shape function, there is a vector
				      * of the 2nd gradients in each quadrature point.
				      */
    typename std::vector<typename std::vector<Tensor<2,dim> > >  shape_2nd_derivatives;

				     /**
				      * Store an array of weights times the
				      * Jacobi determinant at the quadrature
				      * points. This function is reset each time
				      * @p{reinit} is called. The Jacobi determinant
				      * is actually the reciprocal value of the
				      * Jacobi matrices stored in this class,
				      * see the general documentation of this
				      * class for more information.
				      */
    std::vector<double>       JxW_values;

				     /**
				      * Array of quadrature points. This array
				      * is set up upon calling @p{reinit} and
				      * contains the quadrature points on the
				      * real element, rather than on the
				      * reference element.
				      */
    std::vector<Point<dim> >  quadrature_points;

				     /**
				      * List of outward normal vectors at the
				      * quadrature points. This field is filled
				      * in by the finite element class.
				      */
    std::vector<Point<dim> >  normal_vectors;

                                     /**
				      * List of boundary forms at the
				      * quadrature points. This field is filled
				      * in by the finite element class.
				      */
    std::vector<Tensor<1,dim> >  boundary_forms;

                                     /**
				      * Original update flags handed
				      * to the constructor of
				      * @p{FEValues}.
				      */
    UpdateFlags          update_flags;
};


/**
 * Common features of @p{FEValues*} classes.
 *
 * @p{FEValues*} objects are programming interfaces to finite element
 * and mapping classes on the one hand side, to cells and quadrature
 * rules on the other side. The reason for their existence is possible
 * optimization. Depending on the type of finite element and mapping,
 * some values can be computed once on the unit cell. Others must be
 * computed on each cell, but maybe computation of several values at
 * the same time offers ways for optimization. Since this interlay may
 * be complex and depends on the actual finite element, it cannot be
 * left to the applications programmer.
 *
 * @p{FEValues*} provides only data handling: computations are left to
 * objects of type @ref{Mapping} and @ref{FiniteElement}. These
 * provide functions @p{get_*_data} and @p{fill_*_values} which are
 * called by the constructor and @p{reinit} functions of
 * @p{FEValues*}, respectively.
 *
 * \subsection{General usage}
 *
 * Usually, an object of @p{FEValues*} is used in integration loops
 * over all cells of a triangulation. To take full advantage of the
 * optimization features, it should be constructed before the
 * loop. Then, it must be re-initialized for each grid cell. This is
 * like a magnifying glass being used to look at one item after the
 * other. A typical piece of code looks like this:
 *
 * \begin{verbatim}
 * FEValues values (mapping, finite_element, quadrature, flags);
 * for (cell = dof_handler.begin_active();
 *      cell != dof_handler.end();
 *      ++cell)
 *   {
 *     values.reinit(cell);
 *     ...
 *   }
 * \end{verbatim}
 *
 *
 *  \subsection{Member functions}
 *
 *  The functions of this class fall into different cathegories:
 *  \begin{itemize}
 *  \item @p{shape_value}, @p{shape_grad}, etc: return one of the values
 *    of this object at a time. These functions are inlined, so this
 *    is the suggested access to all finite element values. There should be
 *    no loss in performance with an optimizing compiler.
 *   
 *  \item @p{get_shape_values}, @p{get_shape_grads}, etc: these return
 *    a reference to a whole field. Usually these fields contain
 *    the values of all trial functions at all quadrature points.
 *
 *  \item @p{get_function_values}, @p{get_function_grads}, @p{...}:
 *    Compute a finite element function or its derivative
 *    in quadrature points.
 *
 *  \item @p{reinit}: initialize the @p{FEValues} object for a certain cell.
 *    This function is not in the present class but only in the derived
 *    classes and has a variable call syntax. 
 *    See the docs for the derived classes for more information.
 * \end{itemize}
 *
 *
 * \subsection{UpdateFlags}
 *
 * The @ref{UpdateFlags} object handed to the constructor is used to
 * determine, which of the data fields to compute. This way, it is
 * possible to avoid expensive computations of useless derivatives.
 * In the beginning, these flags are processed through the functions
 * @p{update_once} and @p{update_each} of @ref{Mapping} and
 * @p{FiniteElement}. All the results are bit-wise or'd and determine
 * the fields actually computed. This enables @ref{Mapping} and
 * @p{FiniteElement} to schedule auxiliary data fields for
 * updating. Still, it is recommended to give ALL needed update flags
 * to @p{FEValues}.
 *
 * @author Wolfgang Bangerth, 1998, Guido Kanschat, 2001
 */
template <int dim>
class FEValuesBase : protected FEValuesData<dim>
{
  public:
				     /**
				      * Number of quadrature points.
				      */
    const unsigned int n_quadrature_points;

				     /**
				      * Number of shape functions
				      * per cell. If we use this base class
				      * to evaluate a finite element on
				      * faces of cells, this is still the
				      * number of degrees of freedom per
				      * cell, not per face.
				      */
    const unsigned int dofs_per_cell;

    
				     /**
				      * Constructor. Set up the array
				      * sizes with @p{n_q_points}
				      * quadrature points, @p{n_dof}
				      * trial functions per cell and
				      * with the given pattern to
				      * update the fields when the
				      * @p{reinit} function of the
				      * derived classes is called. The
				      * fields themselves are not set
				      * up, this must happen in the
				      * constructor of the derived
				      * class.
				      */
    FEValuesBase (const unsigned int n_q_points,
		  const unsigned int dofs_per_cell,
		  const UpdateFlags         update_flags,
		  const Mapping<dim>       &mapping,
		  const FiniteElement<dim> &fe);


				     /**
				      * Destructor.
				      */
    ~FEValuesBase ();
    
				     /**
				      * Value of the @p{function_no}th shape
				      * function at the @p{point_no}th quadrature
				      * point on the cell, face or
				      * subface selected the last time
				      * the @p{reinit} function of the
				      * derived class was called.
				      */
    double shape_value (const unsigned int function_no,
			const unsigned int point_no) const;

				     /**
				      * Pointer to the matrix holding
				      * all values of shape functions at all
				      * integration points, on the present cell,
				      * face or subface selected
				      * the last time the @p{reinit} function
				      * of the derived class was called.
				      * For the format of this matrix, see the
				      * documentation for the matrix itself.
				      *
				      * Use of this function is deprecated.
				      */
    const typename FEValuesData<dim>::ShapeVector & get_shape_values () const;

				     /**
				      * Values of the finite
				      * element function characterized
				      * by @p{fe_function} restricted to
				      * the cell, face or subface selected
				      * the last time the @p{reinit} function
				      * of the derived class was called,
				      * at the quadrature points.
				      *
				      * If the present cell is not an active
				      * one the interpolated function values
				      * are returned.
				      *
				      * To get values of
				      * multi-component elements,
				      * there is another
				      * @p{get_function_values}
				      * returning a vector of vectors of
				      * results.
				      *
				      * The function assumes that the
				      * @p{values} object already has the
				      * correct size. 
				      *
				      * The actual data type of the
				      * input vector may be either a
				      * @p{Vector<double>},
				      * @p{Vector<float>}, or
				      * @p{BlockVector<double,...>}.
				      */
    template <class InputVector, typename number>
    void get_function_values (const InputVector& fe_function,
			      typename std::vector<number>& values) const;

				     /**
				      * Access to vector valued finite
				      * element functions.
				      *
				      * This function does the same as
				      * the other
				      * @p{get_function_values}, but
				      * applied to multi-component
				      * elements.
				      *
				      * The actual data type of the
				      * input vector may be either a
				      * @p{Vector<double>},
				      * @p{Vector<float>}, or
				      * @p{BlockVector<double,...>}.
				      */
    template <class InputVector, typename number>
    void get_function_values (const InputVector       &fe_function,
			      typename std::vector<Vector<number> > &values) const;

    				     /**
				      * Gradient of the @p{i}th shape
				      * function at the @p{j} quadrature
				      * point with respect to real
				      * cell coordinates.  If you want
				      * to get the derivative in one
				      * of the coordinate directions,
				      * use the appropriate function
				      * of the @ref{Tensor} class to
				      * extract one component. Since
				      * only a reference to the
				      * gradient's value is returned,
				      * there should be no major
				      * performance drawback.
				      */
    const Tensor<1,dim> & shape_grad (const unsigned int function,
				      const unsigned int quadrature_point) const;

				     /** 
				      * Pointer to the matrix holding
				      * all gradients of shape
				      * functions at all integration
				      * points, on the present cell.
				      * For the format of this matrix,
				      * see the documentation for the
				      * matrix itself.
				      *
				      * Use of this function is
				      * deprecated.
				     */
    const typename FEValuesData<dim>::GradientVector & get_shape_grads () const;

				     /**
				      * Gradients of the finite
				      * element function characterized
				      * by @p{fe_function} restricted to
				      * @p{cell} at the quadrature points.
				      *
				      * If the present cell is not an active
				      * one the interpolated function values
				      * are returned.
				      *
				      * The function assumes that the
				      * @p{gradients} object already has the
				      * right size.
				      *
				      * The actual data type of the
				      * input vector may be either a
				      * @p{Vector<double>},
				      * @p{Vector<float>}, or
				      * @p{BlockVector<double,...>}.
				      */
    template <class InputVector>
    void get_function_grads (const InputVector      &fe_function,
			     typename std::vector<Tensor<1,dim> > &gradients) const;

				     /**
				      * Return the gradients of the finite
				      * element function characterized
				      * by @p{fe_function} restricted to
				      * @p{cell} at the quadrature points.
				      *
				      * If the present cell is not an active
				      * one the interpolated function values
				      * are returned.
				      *
				      * The function assumes that the
				      * @p{gradients} object already has the
				      * right size.
				      *
				      * This function does the same as
				      * the other @p{get_function_values},
				      * but applied to multi-component
				      * elements.
				      *
				      * The actual data type of the
				      * input vector may be either a
				      * @p{Vector<double>},
				      * @p{Vector<float>}, or
				      * @p{BlockVector<double,...>}.
				      */
    template <class InputVector>
    void get_function_grads (const InputVector               &fe_function,
			     typename std::vector<typename std::vector<Tensor<1,dim> > > &gradients) const;

    				     /**
				      * 2nd derivatives of
				      * the @p{function_no}th shape function at
				      * the @p{point_no}th quadrature point
				      * with respect to real cell
				      * coordinates. If you want to
				      * get the derivatives in one of
				      * the coordinate directions, use
				      * the appropriate function of
				      * the @p{Tensor} class to
				      * extract one component. Since
				      * only a reference to the
				      * derivative values is returned,
				      * there should be no major
				      * performance drawback.  */
    const Tensor<2,dim> & shape_2nd_derivative (const unsigned int function_no,
						const unsigned int point_no) const;

				     /**
				      * Pointer to the
				      * matrix holding all 2nd
				      * derivatives of shape functions
				      * at all integration points, on
				      * the present cell.  For the
				      * format of this matrix, see the
				      * documentation for the matrix
				      * itself.
				      */
    const typename std::vector<typename std::vector<Tensor<2,dim> > > & get_shape_2nd_derivatives () const;
    
				     /**
				      * Tensor of second derivatives
				      * of the finite element function
				      * characterized by
				      * @p{fe_function} restricted to
				      * @p{cell} at the quadrature
				      * points.
				      *
				      * The function assumes that the
				      * @p{second_derivatives} object
				      * already has the correct size.
				      *
				      * The actual data type of the
				      * input vector may be either a
				      * @p{Vector<double>},
				      * @p{Vector<float>}, or
				      * @p{BlockVector<double,...>}.
				      */
    template <class InputVector>
    void get_function_2nd_derivatives (const InputVector& fe_function,
				       typename std::vector<Tensor<2,dim> >& second_derivatives) const;

    
				     /**
				      * Return the tensor of second
				      * derivatives of the finite
				      * element function characterized
				      * by @p{fe_function} restricted to
				      * @p{cell} at the quadrature points.
				      *
				      * The function assumes that the
				      * @p{second_derivatives} object already has
				      * the right size.
				      *
				      * This function does the same as
				      * the other one with the same
				      * name, but applies to
				      * vector-valued finite elements.
				      *
				      * The actual data type of the
				      * input vector may be either a
				      * @p{Vector<double>},
				      * @p{Vector<float>}, or
				      * @p{BlockVector<double,...>}.
				      */
    template <class InputVector>
    void get_function_2nd_derivatives (const InputVector      &fe_function,
				       typename std::vector<typename std::vector<Tensor<2,dim> > > &second_derivatives) const;
    
				     /**
				      * Position of the @p{i}th
				      * quadrature point in real space.
				      */
    const Point<dim> & quadrature_point (const unsigned int i) const;

				     /**
				      * Return a pointer to the vector of
				      * quadrature points.
				      */
    const typename std::vector<Point<dim> > & get_quadrature_points () const;

				     /**
				      * Mapped quadrature weight. This
				      * is the Jacobi determinant
				      * times the weight of the
				      * @p{i}th unit quadrature point.
				      *
				      * On faces, this is the mapped
				      * surface element.
				      */
    double JxW (const unsigned int quadrature_point) const;

				     /**
				      * Pointer to the array holding
				      * the Jacobi determinant times the
				      * quadrature weight at the different
				      * quadrature points.
				      */
    const std::vector<double> & get_JxW_values () const;
    
				     /**
				      * Return the present cell.
				      */
    const typename DoFHandler<dim>::cell_iterator & get_cell() const;

				     /**
				      * Constant reference to the
				      * selected mapping object.
				      */
    const Mapping<dim> & get_mapping () const;

				     /**
				      * Constant reference to the
				      * selected finite element
				      * object.
				      */
    const FiniteElement<dim> & get_fe () const;

				     /**
				      * Determine an estimate for the
				      * memory consumption (in bytes)
				      * of this object.
				      */
    unsigned int memory_consumption () const;

				     /**
				      * Exception
				      */
    DeclException0 (ExcAccessToUninitializedField);
				     /**
				      * Exception
				      */
    DeclException0 (ExcCannotInitializeField);
				     /**
				      * Exception
				      */
    DeclException0 (ExcWrongNoOfComponents);
				     /**
				      * Exception.
				      */
    DeclException2 (ExcWrongVectorSize,
		    int, int,
		    << "Vector has wrong size " << arg1
		    << ", expected size " << arg2);
				     /**
				      * Exception
				      */
    DeclException0 (ExcInvalidUpdateFlag);
				     /**
				      * Exception
				      */
    DeclException0 (ExcFEDontMatch);
    
  protected:
				     /**
				      * Store the cell selected last time
				      * the @p{reinit} function was called
				      * to make access
				      * to the @p{get_function_*} functions
				      * safer.
				      */
    typename DoFHandler<dim>::cell_iterator present_cell;

				     /**
				      * Storage for the mapping object.
				      */
    const SmartPointer<const Mapping<dim> > mapping;
    
				     /**
				      * Store the finite element for later use.
				      */
    const SmartPointer<const FiniteElement<dim> > fe;

    
				     /**
				      * Internal data of mapping.
				      */
    SmartPointer<typename Mapping<dim>::InternalDataBase> mapping_data;

				     /**
				      * Internal data of finite element.
				      */
    SmartPointer<typename Mapping<dim>::InternalDataBase> fe_data;

				     /**
				      * Initialize some update
				      * flags. Called from the
				      * @p{initialize} functions of
				      * derived classes, which are in
				      * turn called from their
				      * constructors.
				      *
				      * Basically, this function finds
				      * out using the finite element
				      * and mapping object already
				      * stored which flags need to be
				      * set to compute everything the
				      * user wants, as expressed
				      * through the flags passed as
				      * argument.
				      */
    UpdateFlags compute_update_flags (const UpdateFlags update_flags) const;

				     /**
				      * Returns reference to default
				      * @p{MappingQ1} object. Needed
				      * by constructors of derived
				      * classes that uses
				      * @p{MappingQ1} implicitly.
				      */
    static const Mapping<dim> &get_default_mapping();
};



/**
 * Finite element evaluated in quadrature points of a cell.
 *
 * This function implements the initialization routines for
 * @ref{FEValuesBase}, if values in quadrature points of a cell are
 * needed. For further documentation see this class.
 * 
 * @author Wolfgang Bangerth, 1998, Guido Kanschat, 2001
 */
template <int dim>
class FEValues : public FEValuesBase<dim>
{
  public:
				     /**
				      * Constructor. Gets cell
				      * independent data from mapping
				      * and finite element objects,
				      * matching the quadrature rule
				      * and update flags.
				      */
    FEValues (const Mapping<dim>       &mapping,
	      const FiniteElement<dim> &fe,
	      const Quadrature<dim>    &quadrature,
	      const UpdateFlags         update_flags);

                                     /**
				      * Constructor. Uses @ref{MappingQ1}
				      * implicitly.
				      */
    FEValues (const FiniteElement<dim> &fe,
	      const Quadrature<dim>    &quadrature,
	      const UpdateFlags         update_flags);
    
				     /**
				      * Reinitialize the gradients, Jacobi
				      * determinants, etc for the given cell
				      * and the given finite element.
				      */
    void reinit (const typename DoFHandler<dim>::cell_iterator &);

				     /**
				      * Return a reference to the copy
				      * of the quadrature formula
				      * stored by this object.
				      */
    const Quadrature<dim> & get_quadrature () const;
    
				     /**
				      * Determine an estimate for the
				      * memory consumption (in bytes)
				      * of this object.
				      */
    unsigned int memory_consumption () const;

  private:
				     /**
				      * Store a copy of the quadrature
				      * formula here.
				      */
    const Quadrature<dim> quadrature;

				     /**
				      * Do work common to the two
				      * constructors.
				      */
    void initialize (const UpdateFlags update_flags);
};


/**
 * Extend the interface of @ref{FEValuesBase} by surface values.
 *
 * On surfaces of mesh cells, normal vectors and boundary forms are
 * additional values that can be computed. This class provides the
 * interface to access those. Implementations are in derived classes
 * @p{FEFaceValues} and @p{FESubfaceValues}.
 *
 * @see{FEValuesBase}
 *
 *  @author Wolfgang Bangerth, 1998, Guido Kanschat, 2000, 2001
 */
template <int dim>
class FEFaceValuesBase : public FEValuesBase<dim>
{
  public:
				     /**
				      * Constructor. Call the constructor of
				      * the base class and set up the arrays
				      * of this class with the right sizes.
				      * Actually filling these arrays is a
				      * duty of the derived class's
				      * constructors.
				      *
				      * @p{n_faces_or_subfaces} is the number
				      * of faces or subfaces that this object
				      * is to store. The actual number depends
				      * on the derived class, for
				      * @p{FEFaceValues} it is @p{2*dim}, while for
				      * the @p{FESubfaceValues} class it is
				      * @p{2*dim*(1<<(dim-1))}, i.e. the number
				      * of faces times the number of subfaces
				      * per face.
				      */
    FEFaceValuesBase (const unsigned int n_q_points,
		      const unsigned int dofs_per_cell,
		      const UpdateFlags         update_flags,
		      const Mapping<dim>       &mapping,
		      const FiniteElement<dim> &fe,
		      const Quadrature<dim-1>& quadrature);

    				     /**
				      * Return the outward normal vector to
				      * the cell at the @p{i}th quadrature
				      * point. The length of the vector
				      * is normalized to one.
				      */
    const Point<dim> & normal_vector (const unsigned int i) const;
    
    				     /**
				      * Boundary form of the
				      * transformation of the cell at
				      * the @p{i}th quadrature point.
				      *
				      * The boundary form is the cross
				      * product of the images of the
				      * unit tangential
				      * vectors. Therefore, it is the
				      * unit normal vector multiplied
				      * with the surface
				      * element. Since it may be
				      * cheaper to compute the
				      * boundary form immediately, use
				      * this value to integrate
				      * @p{n.ds}.
				      */
    const Tensor<1,dim> & boundary_form (const unsigned int i) const;
    
				     /**
				      * Return the list of outward normal
				      * vectors to the cell at the
				      * quadrature points.
				      */
    const typename std::vector<Point<dim> > & get_normal_vectors () const;

				     /**
				      * Return the list of outward normal
				      * vectors times quadrature weights.
				      */
    const std::vector<Tensor<1,dim> > & get_boundary_forms () const;

				     /**
				      * Return the present
				      * face.
				      */
    typename DoFHandler<dim>::face_iterator get_face() const;

				     /**
				      * Return a reference to the copy
				      * of the quadrature formula
				      * stored by this object.
				      */
    const Quadrature<dim-1> & get_quadrature () const;
    
				     /**
				      * Determine an estimate for the
				      * memory consumption (in bytes)
				      * of this object.
				      */
    unsigned int memory_consumption () const;
    
  protected:
				     /**
				      * Store a copy of the quadrature
				      * formula here.
				      */
    const Quadrature<dim-1> quadrature;

				     /**
				      * Stores the face or subface,
				      * resp., that was selected the
				      * last time the @p{reinit}
				      * function was called. Is used
				      * by the @p{get_face} function.
				      */
    typename DoFHandler<dim>::face_iterator present_face;
};



/**
 * Finite element evaluated in quadrature points on a face.
 *
 * This class adds the functionality of @ref{FEFaceValuesBase} to
 * @ref{FEValues}; see there for more documentation.
 *
 * Since finite element functions and their derivatives may be
 * discontinuous at cell boundaries, there is no restriction of this
 * function to a mesh face. But, there are limits of these values
 * approaching the face from either of the neighboring cells.
 *
 * @author Wolfgang Bangerth, 1998, Guido Kanschat, 2000, 2001
 */
template <int dim>
class FEFaceValues : public FEFaceValuesBase<dim>
{
  public:
				     /**
				      * Constructor. Gets cell
				      * independent data from mapping
				      * and finite element objects,
				      * matching the quadrature rule
				      * and update flags.
				      */
    FEFaceValues (const Mapping<dim>& mapping,
		  const FiniteElement<dim>& fe,
		  const Quadrature<dim-1>& quadrature,
		  const UpdateFlags update_flags);

                                     /**
				      * Constructor. Uses @ref{MappingQ1}
				      * implicitly.
				      */
    FEFaceValues (const FiniteElement<dim>& fe,
		  const Quadrature<dim-1>& quadrature,
		  const UpdateFlags update_flags);

				     /**
				      * Reinitialize the gradients, Jacobi
				      * determinants, etc for the face with
				      * number @p{face_no} of @p{cell}
				      * and the given finite element.
				      */
    void reinit (const typename DoFHandler<dim>::cell_iterator &cell,
		 const unsigned int                    face_no);

  private:

				     /**
				      * Do work common to the two
				      * constructors.
				      */
    void initialize (const UpdateFlags update_flags);    
};


/**
 * Finite element evaluated in quadrature points on a face.
 *
 * This class adds the functionality of @ref{FEFaceValuesBase} to
 * @ref{FEValues}; see there for more documentation.
 *
 * This class is used for faces lying on a refinement edge. In this
 * case, the neighboring cell is refined. To be able to compute
 * differences between interior and exterior function values, the
 * refinement of the neighboring cell must be simulated on this
 * cell. This is achieved by applying a quadrature rule that simulates
 * the refinement. The resulting data fields are split up to reflect
 * the refinement structure of the neighbor: a subface number
 * corresponds to the number of the child of the neighboring face.
 *
 * @author Wolfgang Bangerth, 1998, Guido Kanschat, 2000, 2001
 */
template <int dim>
class FESubfaceValues : public FEFaceValuesBase<dim>
{
  public:
				     /**
				      * Constructor. Gets cell
				      * independent data from mapping
				      * and finite element objects,
				      * matching the quadrature rule
				      * and update flags.
				      */
    FESubfaceValues (const Mapping<dim>& mapping,
		     const FiniteElement<dim>& fe,
		     const Quadrature<dim-1>& face_quadrature,
		     const UpdateFlags update_flags);

                                     /**
				      * Constructor. Uses @ref{MappingQ1}
				      * implicitly.
				      */
    FESubfaceValues (const FiniteElement<dim>& fe,
		     const Quadrature<dim-1>& face_quadrature,
		     const UpdateFlags update_flags);

				     /**
				      * Reinitialize the gradients, Jacobi
				      * determinants, etc for the face with
				      * number @p{face_no} of @p{cell}
				      * and the given finite element.
				      */
    void reinit (const typename DoFHandler<dim>::cell_iterator &cell,
		 const unsigned int                    face_no,
		 const unsigned int                    subface_no);

				     /**
				      * Exception
				      */
    DeclException0 (ExcReinitCalledWithBoundaryFace);

				     /**
				      * Exception
				      */
    DeclException0 (ExcFaceHasNoSubfaces);

  private:

				     /**
				      * Do work common to the two
				      * constructors.
				      */
    void initialize (const UpdateFlags update_flags);
};



/*------------------------ Inline functions: FEValuesBase ------------------------*/


template <int dim>
inline
double
FEValuesBase<dim>::shape_value (const unsigned int i,
				const unsigned int j) const
{
  Assert (update_flags & update_values, ExcAccessToUninitializedField());
  return shape_values(i,j);
}


template <int dim>
inline
const Tensor<1,dim> &
FEValuesBase<dim>::shape_grad (const unsigned int i,
			       const unsigned int j) const
{
  Assert (i<shape_gradients.size(),
	  ExcIndexRange (i, 0, shape_gradients.size()));
  Assert (j<shape_gradients[i].size(),
	  ExcIndexRange (j, 0, shape_gradients[i].size()));
  Assert (update_flags & update_gradients, ExcAccessToUninitializedField());

  return shape_gradients[i][j];
};



template <int dim>
inline
const Tensor<2,dim> &
FEValuesBase<dim>::shape_2nd_derivative (const unsigned int i,
					 const unsigned int j) const
{
  Assert (i<shape_2nd_derivatives.size(),
	  ExcIndexRange (i, 0, shape_2nd_derivatives.size()));
  Assert (j<shape_2nd_derivatives[i].size(),
	  ExcIndexRange (j, 0, shape_2nd_derivatives[i].size()));
  Assert (update_flags & update_second_derivatives, ExcAccessToUninitializedField());

  return shape_2nd_derivatives[i][j];
}



template <int dim>
inline
const typename DoFHandler<dim>::cell_iterator &
FEValuesBase<dim>::get_cell() const
{
  return present_cell;
}



template <int dim>
inline
const FiniteElement<dim> & 
FEValuesBase<dim>::get_fe () const
{
  return *fe;
}


template <int dim>
inline
const Mapping<dim> & 
FEValuesBase<dim>::get_mapping () const
{
  return *mapping;
}


/*------------------------ Inline functions: FEValues ----------------------------*/


template <int dim>
const Quadrature<dim> &
FEValues<dim>::get_quadrature () const 
{
  return quadrature;
};


/*------------------------ Inline functions: FEFaceValuesBase --------------------*/


template <int dim>
const Point<dim> &
FEFaceValuesBase<dim>::normal_vector (const unsigned int i) const
{
  Assert (i<normal_vectors.size(), ExcIndexRange(i, 0, normal_vectors.size()));
  Assert (update_flags & update_normal_vectors,
	  FEValuesBase<dim>::ExcAccessToUninitializedField());
  
  return normal_vectors[i];
};


template <int dim>
inline
typename DoFHandler<dim>::face_iterator
FEFaceValuesBase<dim>::get_face() const 
{
  return present_face;
};




template <int dim>
const Quadrature<dim-1> &
FEFaceValuesBase<dim>::get_quadrature () const 
{
  return quadrature;
};



#endif
