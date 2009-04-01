//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name:  $
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__fe_values_h
#define __deal2__fe_values_h


#include <base/config.h>
#include <base/exceptions.h>
#include <base/subscriptor.h>
#include <base/point.h>
#include <base/tensor.h>
#include <base/symmetric_tensor.h>
#include <base/vector_slice.h>
#include <base/quadrature.h>
#include <base/table.h>
#include <lac/vector.h>
#include <lac/block_vector.h>
#include <lac/petsc_vector.h>
#include <lac/petsc_block_vector.h>
#include <lac/trilinos_vector.h>
#include <lac/trilinos_block_vector.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <hp/dof_handler.h>
#include <fe/fe.h>
#include <fe/fe_update_flags.h>
#include <fe/mapping.h>
#include <fe/mapping_q.h>
#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_dof_accessor.h>

#include <algorithm>
#include <memory>

DEAL_II_NAMESPACE_OPEN

template <int dim>   class Quadrature;
template <int dim, int spacedim=dim> class FEValuesBase;



/**
 * A namespace in which we declare "extractors", i.e. classes that when used
 * as subscripts in operator[] expressions on FEValues, FEFaceValues, and
 * FESubfaceValues objects extract certain components of a vector-valued
 * element. There are extractors for single scalar components as well as
 * vector components consisting of <code>dim</code> elements.
 *
 * See the description of the @ref vector_valued module for examples how to
 * use the features of this namespace.
 *
 * @ingroup feaccess vector_valued
 */
namespace FEValuesExtractors
{
				   /**
				    * Extractor for a single scalar component
				    * of a vector-valued element. The concept
				    * of extractors is defined in the
				    * documentation of the namespace
				    * FEValuesExtractors and in the @ref
				    * vector_valued module.
				    *
				    * @ingroup feaccess vector_valued
				    */
  struct Scalar
  {
				       /**
					* The selected scalar component of the
					* vector.
					*/
      unsigned int component;

				       /**
					* Constructor. Take the selected
					* vector component as argument.
					*/
      Scalar (const unsigned int component);
  };

  
				   /**
				    * Extractor for a vector of
				    * <code>dim</code> components of a
				    * vector-valued element. The value of
				    * <code>dim</code> is defined by the
				    * FEValues object the extractor is applied
				    * to.
				    *
				    * The concept of
				    * extractors is defined in the
				    * documentation of the namespace
				    * FEValuesExtractors and in the @ref
				    * vector_valued module.
				    *
				    * @ingroup feaccess vector_valued
				    */
  struct Vector
  {
				       /**
					* The first component of the vector
					* view.
					*/
      unsigned int first_vector_component;

				       /**
					* Constructor. Take the first
					* component of the selected vector
					* inside the FEValues object as
					* argument.
					*/
      Vector (const unsigned int first_vector_component);
  };
}



/**
 * A namespace for "views" on a FEValues, FEFaceValues, or FESubfaceValues
 * object. A view represents only a certain part of the whole: whereas the
 * FEValues object represents <i>all</i> values, gradients, or second
 * derivatives of all components of a vector-valued element, views restrict
 * the attention to only a single component or a subset of components.
 *
 * There are classes that present views for single scalar components as well
 * as vector components consisting of <code>dim</code> elements.
 *
 * See the description of the @ref vector_valued module for examples how to
 * use the features of this namespace.
 *
 * @ingroup feaccess vector_valued
 */
namespace FEValuesViews
{
				   /**
				    * A class representing a view to a single
				    * scalar component of a possibly
				    * vector-valued finite element. Views are
				    * discussed in the @ref vector_valued module.
				    *
				    * @ingroup feaccess vector_valued
				    */
  template <int dim, int spacedim=dim>  
  class Scalar
  {
    public:
				       /**
					* A typedef for the data type of
					* values of the view this class
					* represents. Since we deal with a
					* single components, the value type is
					* a scalar double.
					*/
      typedef double        value_type;

				       /**
					* A typedef for the type of gradients
					* of the view this class
					* represents. Here, for a scalar
					* component of the finite element, the
					* gradient is a
					* <code>Tensor@<1,dim@></code>.
					*/
      typedef Tensor<1,spacedim> gradient_type;

				       /**
					* A typedef for the type of second
					* derivatives of the view this class
					* represents. Here, for a scalar
					* component of the finite element, the
					* Hessian is a
					* <code>Tensor@<2,dim@></code>.
					*/
      typedef Tensor<2,spacedim> hessian_type;

				       /**
					* Default constructor. Creates an
					* invalid object.
					*/
      Scalar ();
      
				       /**
					* Constructor for an object that
					* represents a single scalar component
					* of a FEValuesBase object (or of one
					* of the classes derived from
					* FEValuesBase).
					*/
      Scalar (const FEValuesBase<dim,spacedim> &fe_values_base,
	      const unsigned int       component);

				       /**
					* Copy operator. This is not a
					* lightweight object so we don't allow
					* copying and generate an exception if
					* this function is called.
					*/
      Scalar & operator= (const Scalar<dim,spacedim> &);
      
				       /**
					* Return the value of the vector
					* component selected by this view, for
					* the shape function and quadrature
					* point selected by the arguments.
					*/
      value_type
      value (const unsigned int shape_function,
	     const unsigned int q_point) const;

				       /**
					* Return the gradient (a tensor of
					* rank 1) of the vector component
					* selected by this view, for the shape
					* function and quadrature point
					* selected by the arguments.
					*/
      gradient_type
      gradient (const unsigned int shape_function,
		const unsigned int q_point) const;

				       /**
					* Return the Hessian (the tensor of
					* rank 2 of all second derivatives) of
					* the vector component selected by
					* this view, for the shape function
					* and quadrature point selected by the
					* arguments.
					*/
      hessian_type
      hessian (const unsigned int shape_function,
	       const unsigned int q_point) const;

				       /**
					* Return the values of the selected
					* scalar component of the finite
					* element function characterized by
					* <tt>fe_function</tt> at the
					* quadrature points of the cell, face
					* or subface selected the last time
					* the <tt>reinit</tt> function of the
					* FEValues object was called, at the
					* quadrature points.
					*
					* This function is the equivalent of
					* the
					* FEValuesBase::get_function_values
					* function but it only works on the
					* selected scalar component.
					*/
      template <class InputVector>
      void get_function_values (const InputVector& fe_function,
				std::vector<value_type>& values) const;

				       /**
					* Return the gradients of the selected
					* scalar component of the finite
					* element function characterized by
					* <tt>fe_function</tt> at the
					* quadrature points of the cell, face
					* or subface selected the last time
					* the <tt>reinit</tt> function of the
					* FEValues object was called, at the
					* quadrature points.
					*
					* This function is the equivalent of
					* the
					* FEValuesBase::get_function_gradients
					* function but it only works on the
					* selected scalar component.
					*/
      template <class InputVector>
      void get_function_gradients (const InputVector& fe_function,
				   std::vector<gradient_type>& gradients) const;

				       /**
					* Return the Hessians of the selected
					* scalar component of the finite
					* element function characterized by
					* <tt>fe_function</tt> at the
					* quadrature points of the cell, face
					* or subface selected the last time
					* the <tt>reinit</tt> function of the
					* FEValues object was called, at the
					* quadrature points.
					*
					* This function is the equivalent of
					* the
					* FEValuesBase::get_function_hessians
					* function but it only works on the
					* selected scalar component.
					*/
      template <class InputVector>
      void get_function_hessians (const InputVector& fe_function,
				  std::vector<hessian_type>& hessians) const;

				       /**
					* Return the Laplacians of the
					* selected scalar component of
					* the finite element function
					* characterized by
					* <tt>fe_function</tt> at the
					* quadrature points of the
					* cell, face or subface
					* selected the last time the
					* <tt>reinit</tt> function of
					* the FEValues object was
					* called, at the quadrature
					* points. The Laplacians are
					* the trace of the Hessians.
					*
					* This function is the
					* equivalent of the
					* FEValuesBase::get_function_laplacians
					* function but it only works
					* on the selected scalar
					* component.
					*/
      template <class InputVector>
      void get_function_laplacians (const InputVector& fe_function,
				    std::vector<value_type>& laplacians) const;
      
      
    private:
				       /**
					* A reference to the FEValuesBase object
					* we operator on.
					*/
      const FEValuesBase<dim,spacedim> &fe_values;

				       /**
					* The single scalar component this
					* view represents of the FEValuesBase
					* object.
					*/
      const unsigned int component;

				       /**
					* A structure where for each shape
					* function we pre-compute a bunch of
					* data that will make later accesses
					* much cheaper.
					*/
      struct ShapeFunctionData
      {
					   /**
					    * For each shape function, store
					    * whether the selected vector
					    * component may be nonzero. For
					    * primitive shape functions we
					    * know for sure whether a certain
					    * scalar component of a given
					    * shape function is nonzero,
					    * whereas for non-primitive shape
					    * functions this may not be
					    * entirely clear (e.g. for RT
					    * elements it depends on the shape
					    * of a cell).
					    */
	  bool is_nonzero_shape_function_component;

					   /**
					    * For each shape function, store
					    * the row index within the
					    * shape_values, shape_gradients,
					    * and shape_hessians tables (the
					    * column index is the quadrature
					    * point index). If the shape
					    * function is primitive, then we
					    * can get this information from
					    * the shape_function_to_row_table
					    * of the FEValues object;
					    * otherwise, we have to work a bit
					    * harder to compute this
					    * information.
					    */
	  unsigned int row_index;
      };

				       /**
					* Store the data about shape
					* functions.
					*/
      std::vector<ShapeFunctionData> shape_function_data;
  };
  


				   /**
				    * A class representing a view to a set of
				    * <code>dim</code> components forming a
				    * vector part of a vector-valued finite
				    * element. Views are discussed in the
				    * @ref vector_valued module.
				    *
				    * @ingroup feaccess vector_valued
				    */
  template <int dim, int spacedim=dim>  
  class Vector
  {
    public:
				       /**
					* A typedef for the data type of
					* values of the view this class
					* represents. Since we deal with a set
					* of <code>dim</code> components, the
					* value type is a Tensor<1,spacedim>.
					*/
      typedef Tensor<1,spacedim>          value_type;

				       /**
					* A typedef for the type of gradients
					* of the view this class
					* represents. Here, for a set of
					* <code>dim</code> components of the
					* finite element, the gradient is a
					* <code>Tensor@<2,spacedim@></code>.
					*/
      typedef Tensor<2,spacedim>          gradient_type;

				       /**
					* A typedef for the type of
					* symmetrized gradients of the view
					* this class represents. Here, for a
					* set of <code>dim</code> components
					* of the finite element, the
					* symmetrized gradient is a
					* <code>SymmetricTensor@<2,spacedim@></code>.
					*/
      typedef SymmetricTensor<2,spacedim> symmetric_gradient_type;

				       /**
					* A typedef for the type of the
					* divergence of the view this class
					* represents. Here, for a set of
					* <code>dim</code> components of the
					* finite element, the divergence of
					* course is a scalar.
					*/
      typedef double                 divergence_type;
      
				       /**
					* A typedef for the type of second
					* derivatives of the view this class
					* represents. Here, for a set of
					* <code>dim</code> components of the
					* finite element, the Hessian is a
					* <code>Tensor@<3,dim@></code>.
					*/
      typedef Tensor<3,spacedim>          hessian_type;

				       /**
					* Default constructor. Creates an
					* invalid object.
					*/
      Vector ();

				       /**
					* Constructor for an object that
					* represents dim components of a
					* FEValuesBase object (or of one of
					* the classes derived from
					* FEValuesBase), representing a
					* vector-valued variable.
					*
					* The second argument denotes the
					* index of the first component of the
					* selected vector.
					*/
      Vector (const FEValuesBase<dim,spacedim> &fe_values_base,
	      const unsigned int first_vector_component);

				       /**
					* Copy operator. This is not a
					* lightweight object so we don't allow
					* copying and generate an exception if
					* this function is called.
					*/
      Vector & operator= (const Vector<dim,spacedim> &);

				       /**
					* Return the value of the vector
					* components selected by this view,
					* for the shape function and
					* quadrature point selected by the
					* arguments. Here, since the view
					* represents a vector-valued part of
					* the FEValues object with
					* <code>dim</code> components, the
					* return type is a tensor of rank 1
					* with <code>dim</code> components.
					*/
      value_type
      value (const unsigned int shape_function,
	     const unsigned int q_point) const;

				       /**
					* Return the gradient (a tensor of
					* rank 2) of the vector component
					* selected by this view, for the shape
					* function and quadrature point
					* selected by the arguments.
					*/
      gradient_type
      gradient (const unsigned int shape_function,
		const unsigned int q_point) const;

				       /**
					* Return the symmetric gradient (a
					* symmetric tensor of rank 2) of the
					* vector component selected by this
					* view, for the shape function and
					* quadrature point selected by the
					* arguments.
					*
					* The symmetric gradient is defined as
					* $\frac 12 [(\nabla \phi_i(x_q)) +
					* (\nabla \phi_i(x_q))^T]$, where
					* $\phi_i$ represents the
					* <code>dim</code> components selected
					* from the FEValuesBase object, and
					* $x_q$ is the location of the $q$-th
					* quadrature point.
					*/
      symmetric_gradient_type
      symmetric_gradient (const unsigned int shape_function,
			  const unsigned int q_point) const;

				       /**
					* Return the scalar divergence of
					* the vector components selected by
					* this view, for the shape function
					* and quadrature point selected by the
					* arguments.
					*/
      divergence_type
      divergence (const unsigned int shape_function,
		  const unsigned int q_point) const;

				       /**
					* Return the Hessian (the tensor of
					* rank 2 of all second derivatives) of
					* the vector components selected by
					* this view, for the shape function
					* and quadrature point selected by the
					* arguments.
					*/
      hessian_type
      hessian (const unsigned int shape_function,
	       const unsigned int q_point) const;

				       /**
					* Return the values of the selected
					* vector components of the finite
					* element function characterized by
					* <tt>fe_function</tt> at the
					* quadrature points of the cell, face
					* or subface selected the last time
					* the <tt>reinit</tt> function of the
					* FEValues object was called, at the
					* quadrature points.
					*
					* This function is the equivalent of
					* the
					* FEValuesBase::get_function_values
					* function but it only works on the
					* selected vector components.
					*/
      template <class InputVector>
      void get_function_values (const InputVector& fe_function,
				std::vector<value_type>& values) const;

				       /**
					* Return the gradients of the selected
					* vector components of the finite
					* element function characterized by
					* <tt>fe_function</tt> at the
					* quadrature points of the cell, face
					* or subface selected the last time
					* the <tt>reinit</tt> function of the
					* FEValues object was called, at the
					* quadrature points.
					*
					* This function is the equivalent of
					* the
					* FEValuesBase::get_function_gradients
					* function but it only works on the
					* selected vector components.
					*/
      template <class InputVector>
      void get_function_gradients (const InputVector& fe_function,
				   std::vector<gradient_type>& gradients) const;

				       /**
					* Return the symmetrized gradients of
					* the selected vector components of
					* the finite element function
					* characterized by
					* <tt>fe_function</tt> at the
					* quadrature points of the cell, face
					* or subface selected the last time
					* the <tt>reinit</tt> function of the
					* FEValues object was called, at the
					* quadrature points.
					*
					* There is no equivalent function such
					* as
					* FEValuesBase::get_function_gradients
					* in the FEValues classes but the
					* information can be obtained from
					* FEValuesBase::get_function_gradients,
					* of course.
					*/
      template <class InputVector>
      void get_function_symmetric_gradients (const InputVector& fe_function,
		    std::vector<symmetric_gradient_type>& symmetric_gradients) const;
      
				       /**
					* Return the divergence of the selected
					* vector components of the finite
					* element function characterized by
					* <tt>fe_function</tt> at the
					* quadrature points of the cell, face
					* or subface selected the last time
					* the <tt>reinit</tt> function of the
					* FEValues object was called, at the
					* quadrature points.
					*
					* There is no equivalent function such
					* as
					* FEValuesBase::get_function_gradients
					* in the FEValues classes but the
					* information can be obtained from
					* FEValuesBase::get_function_gradients,
					* of course.
					*/
      template <class InputVector>
      void get_function_divergences (const InputVector& fe_function,
				     std::vector<divergence_type>& divergences) const;
      
				       /**
					* Return the Hessians of the selected
					* vector components of the finite
					* element function characterized by
					* <tt>fe_function</tt> at the
					* quadrature points of the cell, face
					* or subface selected the last time
					* the <tt>reinit</tt> function of the
					* FEValues object was called, at the
					* quadrature points.
					*
					* This function is the equivalent of
					* the
					* FEValuesBase::get_function_hessians
					* function but it only works on the
					* selected vector components.
					*/
      template <class InputVector>
      void get_function_hessians (const InputVector& fe_function,
				  std::vector<hessian_type>& hessians) const;

				       /**
					* Return the Laplacians of the
					* selected vector components
					* of the finite element
					* function characterized by
					* <tt>fe_function</tt> at the
					* quadrature points of the
					* cell, face or subface
					* selected the last time the
					* <tt>reinit</tt> function of
					* the FEValues object was
					* called, at the quadrature
					* points. The Laplacians are
					* the trace of the Hessians.
					*
					* This function is the
					* equivalent of the
					* FEValuesBase::get_function_laplacians
					* function but it only works
					* on the selected vector
					* components.
					*/
      template <class InputVector>
      void get_function_laplacians (const InputVector& fe_function,
				    std::vector<value_type>& laplacians) const;
      
    private:
				       /**
					* A reference to the FEValuesBase object
					* we operator on.
					*/
      const FEValuesBase<dim,spacedim> &fe_values;

				       /**
					* The first component of the vector
					* this view represents of the
					* FEValuesBase object.
					*/
      const unsigned int first_vector_component;

				       /**
					* A structure where for each shape
					* function we pre-compute a bunch of
					* data that will make later accesses
					* much cheaper.
					*/
      struct ShapeFunctionData
      {      
					   /**
					    * For each pair (shape
					    * function,component within
					    * vector), store whether the
					    * selected vector component may be
					    * nonzero. For primitive shape
					    * functions we know for sure
					    * whether a certain scalar
					    * component of a given shape
					    * function is nonzero, whereas for
					    * non-primitive shape functions
					    * this may not be entirely clear
					    * (e.g. for RT elements it depends
					    * on the shape of a cell).
					    */
	  bool is_nonzero_shape_function_component[dim];

					   /**
					    * For each pair (shape function,
					    * component within vector), store
					    * the row index within the
					    * shape_values, shape_gradients,
					    * and shape_hessians tables (the
					    * column index is the quadrature
					    * point index). If the shape
					    * function is primitive, then we
					    * can get this information from
					    * the shape_function_to_row_table
					    * of the FEValues object;
					    * otherwise, we have to work a bit
					    * harder to compute this
					    * information.
					    */
	  unsigned int row_index[dim];

					   /**
					    * For each shape function say the
					    * following: if only a single
					    * entry in
					    * is_nonzero_shape_function_component
					    * for this shape function is
					    * nonzero, then store the
					    * corresponding value of row_index
					    * and
					    * single_nonzero_component_index
					    * represents the index between 0
					    * and dim for which it is
					    * attained. If multiple components
					    * are nonzero, then store -1. If
					    * no components are nonzero then
					    * store -2.
					    */
	  int          single_nonzero_component;
	  unsigned int single_nonzero_component_index;
      };

				       /**
					* Store the data about shape
					* functions.
					*/
      std::vector<ShapeFunctionData> shape_function_data;
  };
}


namespace internal
{
  namespace FEValuesViews
  {
				     /**
				      * A class objects of which store a
				      * collection of FEValuesViews::Scalar,
				      * FEValuesViews::Vector, etc object. The
				      * FEValuesBase class uses it to generate
				      * all possible Views classes upon
				      * construction time; we do this at
				      * construction time since the Views
				      * classes cache some information and are
				      * therefore relatively expensive to
				      * create.
				      */
    template <int dim, int spacedim>
    struct Cache
    {
					 /**
					  * Caches for scalar and
					  * vector-valued views.
					  */
	std::vector<dealii::FEValuesViews::Scalar<dim,spacedim> > scalars;
	std::vector<dealii::FEValuesViews::Vector<dim,spacedim> > vectors;

					 /**
					  * Constructor.
					  */
	Cache (const FEValuesBase<dim,spacedim> &fe_values);
    };
  }
}



  
/*!@addtogroup feaccess */
/*@{*/



//TODO: Add access to mapping values to FEValuesBase
//TODO: Several FEValuesBase of a system should share Mapping

/**
 * Contains all data vectors for FEValues.
 * This class has been extracted from FEValuesBase to be handed
 * over to the fill functions of Mapping and
 * FiniteElement.
 *
 * @note All data fields are public, but this is not
 * critical, because access to this object is private in FEValues.
 *
 * The purpose of this class is discussed
 * on the page on @ref UpdateFlagsEssay.
 * 
 * @author Guido Kanschat, 2000
 */
template <int dim, int spacedim=dim>
class FEValuesData
{
  public:
				     /**
				      * Initialize all vectors to
				      * correct size.
				      */
    void initialize (const unsigned int        n_quadrature_points,
		     const FiniteElement<dim,spacedim> &fe,
		     const UpdateFlags         flags);

				     /**
				      * Storage type for shape
				      * values. Each row in the matrix
				      * denotes the values of a single
				      * shape function at the
				      * different points, columns are
				      * for a single point with the
				      * different shape functions.
				      *
				      * If a shape function has more
				      * than one non-zero component
				      * (in deal.II diction: it is
				      * non-primitive), then we
				      * allocate one row per non-zero
				      * component, and shift
				      * subsequent rows backward.
				      * Lookup of the correct row for
				      * a shape function is thus
				      * simple in case the entire
				      * finite element is primitive
				      * (i.e. all shape functions are
				      * primitive), since then the
				      * shape function number equals
				      * the row number. Otherwise, use
				      * the
				      * #shape_function_to_row_table
				      * array to get at the first row
				      * that belongs to this
				      * particular shape function, and
				      * navigate among all the rows
				      * for this shape function using
				      * the
				      * FiniteElement::get_nonzero_components()
				      * function which tells us which
				      * components are non-zero and
				      * thus have a row in the array
				      * presently under discussion.
				      */
    typedef Table<2,double> ShapeVector;

				     /**
				      * Storage type for
				      * gradients. The layout of data
				      * is the same as for the
				      * #ShapeVector data type.
				      */
    typedef std::vector<std::vector<Tensor<1,spacedim> > > GradientVector;

				     /**
				      * Likewise for second order
				      * derivatives.
				      */
    typedef std::vector<std::vector<Tensor<2,spacedim> > > HessianVector;

				     /**
				      * Store the values of the shape
				      * functions at the quadrature
				      * points. See the description of
				      * the data type for the layout
				      * of the data in this field.
				      */
    ShapeVector shape_values;

    				     /**
				      * Store the gradients of the
				      * shape functions at the
				      * quadrature points. See the
				      * description of the data type
				      * for the layout of the data in
				      * this field.
				      */
    GradientVector shape_gradients;

				     /**
				      * Store the 2nd derivatives of
				      * the shape functions at the
				      * quadrature points.  See the
				      * description of the data type
				      * for the layout of the data in
				      * this field.
				      */
    HessianVector shape_hessians;

				     /**
				      * Store an array of weights
				      * times the Jacobi determinant
				      * at the quadrature points. This
				      * function is reset each time
				      * reinit() is called. The
				      * Jacobi determinant is actually
				      * the reciprocal value of the
				      * Jacobi matrices stored in this
				      * class, see the general
				      * documentation of this class
				      * for more information.
				      *
				      * However, if this object refers
				      * to an FEFaceValues or
				      * FESubfaceValues object, then
				      * the JxW_values correspond to
				      * the Jacobian of the
				      * transformation of the face,
				      * not the cell, i.e. the
				      * dimensionality is that of a
				      * surface measure, not of a
				      * volume measure. In this case,
				      * it is computed from the
				      * boundary forms, rather than
				      * the Jacobian matrix.
				      */
    std::vector<double>       JxW_values;

				     /**
				      * Array of the Jacobian matrices at the
				      * quadrature points.
				      */
    std::vector<Tensor<2,spacedim> > jacobians;

				     /**
				      * Array of the derivatives of the Jacobian
				      * matrices at the quadrature points.
				      */
    std::vector<Tensor<3,spacedim> > jacobian_grads;

				     /**
				      * Array of the inverse Jacobian matrices 
				      * at the quadrature points.
				      */
    std::vector<Tensor<2,spacedim> > inverse_jacobians;

				     /**
				      * Array of quadrature points. This array
				      * is set up upon calling reinit() and
				      * contains the quadrature points on the
				      * real element, rather than on the
				      * reference element.
				      */
    std::vector<Point<spacedim> >  quadrature_points;

				     /**
				      * List of outward normal vectors at the
				      * quadrature points. This field is filled
				      * in by the finite element class.
				      */
    std::vector<Point<spacedim> >  normal_vectors;

				     /**
				      * List of outward vectors normal to the cell  
				      * surface (line) at the quadrature points
				      * for the codimension 1 case,
				      * when spacedim=3 (=2).
				      */
    std::vector<Point<spacedim> >  cell_normal_vectors;

                                     /**
				      * List of boundary forms at the
				      * quadrature points. This field is filled
				      * in by the finite element class.
				      */
    std::vector<Tensor<1,spacedim> >  boundary_forms;

				     /**
				      * Array of the mapped support
				      * points, filled by Mapping.
				      */
    std::vector<Point<spacedim> > support_points;

				     /**
				      * Array of the Jacobian of the
				      * mapping in the support points,
				      * filled by Mapping.
				      */
    HessianVector support_jacobians;
    
				     /**
				      * Array of the inverse Jacobian of the
				      * mapping in the support points,
				      * filled by Mapping.
				      */
    HessianVector support_inverse_jacobians;
    
				     /**
				      * Indicate the first row which a
				      * given shape function occupies
				      * in the #shape_values,
				      * #shape_gradients and
				      * #shape_hessians
				      * arrays. If all shape functions
				      * are primitive, then this is
				      * the identity mapping. If, on
				      * the other hand some shape
				      * functions have more than one
				      * non-zero vector components,
				      * then they may occupy more than
				      * one row, and this array
				      * indicates which is the first
				      * one.
				      *
				      * The questions which particular
				      * vector component occupies
				      * which row for a given shape
				      * function is answered as
				      * follows: we allocate one row
				      * for each non-zero component as
				      * indicated by the
				      * FiniteElement::get_nonzero_components()
				      * function, and the rows are in
				      * ascending order exactly those
				      * non-zero components.
				      */
    std::vector<unsigned int> shape_function_to_row_table;
    
                                     /**
				      * Original update flags handed
				      * to the constructor of
				      * FEValues.
				      */
    UpdateFlags          update_flags;
};


/**
 * FEValues, FEFaceValues and FESubfaceValues objects are interfaces to finite
 * element and mapping classes on the one hand side, to cells and quadrature
 * rules on the other side. They allow to evaluate values or derivatives of
 * shape functions at the quadrature points of a quadrature formula when
 * projected by a mapping from the unit cell onto a cell in real space. The
 * reason for this abstraction is possible optimization: Depending on the type
 * of finite element and mapping, some values can be computed once on the unit
 * cell. Others must be computed on each cell, but maybe computation of
 * several values at the same time offers ways for optimization. Since this
 * interlay may be complex and depends on the actual finite element, it cannot
 * be left to the applications programmer.
 *
 * FEValues, FEFaceValues and FESubfaceValues provide only data handling:
 * computations are left to objects of type Mapping and FiniteElement. These
 * provide functions <tt>get_*_data</tt> and <tt>fill_*_values</tt> which are
 * called by the constructor and <tt>reinit</tt> functions of
 * <tt>FEValues*</tt>, respectively.
 *
 * <h3>General usage</h3>
 *
 * Usually, an object of <tt>FEValues*</tt> is used in integration loops over
 * all cells of a triangulation (or faces of cells). To take full advantage of
 * the optimization features, it should be constructed before the loop so that
 * information that does not depend on the location and shape of cells can be
 * computed once and for all (this includes, for example, the values of shape
 * functions at quadrature points for the most common elements: we can
 * evaluate them on the unit cell and they will be the same when mapped to the
 * real cell). Then, in the loop over all cells, it must be re-initialized for
 * each grid cell to compute that part of the information that changes
 * depending on the actual cell (for example, the gradient of shape functions
 * equals the gradient on the unit cell -- which can be computed once and for
 * all -- times the Jacobian matrix of the mapping between unit and real cell,
 * which needs to be recomputed for each cell).
 *
 * A typical piece of code, adding up local contributions to the Laplace
 * matrix looks like this:
 *
 * @code
 * FEValues values (mapping, finite_element, quadrature, flags);
 * for (cell = dof_handler.begin_active();
 *      cell != dof_handler.end();
 *      ++cell)
 *   {
 *     values.reinit(cell);
 *     for (unsigned int q=0; q<quadrature.size(); ++q)
 *       for (unsigned int i=0; i<finite_element.dofs_per_cell; ++i)
 *         for (unsigned int j=0; j<finite_element.dofs_per_cell; ++j)
 *         A(i,j) += fe_values.shape_value(i,q) *
 *                   fe_values.shape_value(j,q) *
 *                   fe_values.JxW(q);
 *     ...
 *   }
 * @endcode
 *
 * The individual functions used here are described below. Note that by
 * design, the order of quadrature points used inside the FEValues object is
 * the same as defined by the quadrature formula passed to the constructor of
 * the FEValues object above.
 *
 *  <h3>Member functions</h3>
 *
 *  The functions of this class fall into different cathegories:
 *  <ul>
 *  <li> shape_value(), shape_grad(), etc: return one of the values 
 *    of this object at a time. These functions are inlined, so this
 *    is the suggested access to all finite element values. There
 *    should be no loss in performance with an optimizing compiler. If
 *    the finite element is vector valued, then these functions return
 *    the only non-zero component of the requested shape
 *    function. However, some finite elements have shape functions
 *    that have more than one non-zero component (we call them
 *    non-"primitive"), and in this case this set of functions will
 *    throw an exception since they cannot generate a useful
 *    result. Rather, use the next set of functions.
 *
 *  <li> shape_value_component(), shape_grad_component(), etc:
 *    This is the same set of functions as above, except that for vector
 *    valued finite elements they return only one vector component. This
 *    is useful for elements of which shape functions have more than one
 *    non-zero component, since then the above functions cannot be used,
 *    and you have to walk over all (or only the non-zero) components of
 *    the shape function using this set of functions.
 *   
 *  <li> get_function_values(), get_function_gradients(), etc.: Compute a
 *    finite element function or its derivative in quadrature points.
 *
 *  <li> reinit: initialize the FEValues object for a certain cell.
 *    This function is not in the present class but only in the derived
 *    classes and has a variable call syntax. 
 *    See the docs for the derived classes for more information.
 * </ul>
 *
 *
 * <h3>UpdateFlags</h3>
 *
 * The UpdateFlags object handed to the constructor is used to
 * determine, which of the data fields to compute. This way, it is
 * possible to avoid expensive computations of useless derivatives.
 * In the beginning, these flags are processed through the functions
 * Mapping::update_once(), Mapping::update_each(),
 * FiniteElement::update_once() FiniteElement::update_each(). All the
 * results are bit-wise or'd and determine the fields actually
 * computed. This enables Mapping and FiniteElement to schedule
 * auxiliary data fields for updating. Still, it is recommended to
 * give <b>all</b> needed update flags to FEValues.
 *
 * The mechanisms by which this class works is also discussed
 * on the page on @ref UpdateFlagsEssay. 
 *
 * @author Wolfgang Bangerth, 1998, 2003, Guido Kanschat, 2001
 */
template <int dim, int spacedim>
class FEValuesBase : protected FEValuesData<dim,spacedim>,
                     public Subscriptor
{
  public:
				     /**
				      * Number of quadrature points.
				      */
    const unsigned int n_quadrature_points;

				     /**
				      * Number of shape functions per
				      * cell. If we use this base
				      * class to evaluate a finite
				      * element on faces of cells,
				      * this is still the number of
				      * degrees of freedom per cell,
				      * not per face.
				      */
    const unsigned int dofs_per_cell;

    
				     /**
				      * Constructor. Set up the array
				      * sizes with <tt>n_q_points</tt>
				      * quadrature points, <tt>dofs_per_cell</tt>
				      * trial functions per cell and
				      * with the given pattern to
				      * update the fields when the
				      * <tt>reinit</tt> function of the
				      * derived classes is called. The
				      * fields themselves are not set
				      * up, this must happen in the
				      * constructor of the derived
				      * class.
				      */
    FEValuesBase (const unsigned int n_q_points,
		  const unsigned int dofs_per_cell,
		  const UpdateFlags update_flags,
		  const Mapping<dim,spacedim> &mapping,
		  const FiniteElement<dim,spacedim> &fe);


				     /**
				      * Destructor.
				      */
    ~FEValuesBase ();

				     /// @name Extractors Methods to extract individual components
				     //@{

				     /**
				      * Create a view of the current FEValues
				      * object that represents a particular
				      * scalar component of the possibly
				      * vector-valued finite element. The
				      * concept of views is explained in the
				      * documentation of the namespace
				      * FEValuesViews and in particular
				      * in the @ref vector_valued module.
				      */
    const FEValuesViews::Scalar<dim,spacedim> &
    operator[] (const FEValuesExtractors::Scalar &scalar) const;

				     /**
				      * Create a view of the current FEValues
				      * object that represents a set of
				      * <code>dim</code> scalar components
				      * (i.e. a vector) of the vector-valued
				      * finite element. The concept of views
				      * is explained in the documentation of
				      * the namespace FEValuesViews and in particular
				      * in the @ref vector_valued module.
				      */
    const FEValuesViews::Vector<dim,spacedim> &
    operator[] (const FEValuesExtractors::Vector &vector) const;

    				     //@}

    
				     /// @name ShapeAccess Access to shape function values
				     //@{
    
				     /**
				      * Value of a shape function at a
				      * quadrature point on the cell,
				      * face or subface selected the
				      * last time the <tt>reinit</tt>
				      * function of the derived class
				      * was called.
				      *
				      * If the shape function is
				      * vector-valued, then this
				      * returns the only non-zero
				      * component. If the shape
				      * function has more than one
				      * non-zero component (i.e. it is
				      * not primitive), then throw an
				      * exception of type
				      * ExcShapeFunctionNotPrimitive. In
				      * that case, use the
				      * shape_value_component()
				      * function.
				      *
				      * @arg function_no Number
				      * of the shape function to be
				      * computed
				      * @arg point_no Number of
				      * the quadrature point at which
				      * function is to be computed
				      */
    double shape_value (const unsigned int function_no,
			const unsigned int point_no) const;

				     /**
				      * Compute one vector component of
				      * the value of a shape function
				      * at a quadrature point. If the
				      * finite element is scalar, then
				      * only component zero is allowed
				      * and the return value equals
				      * that of the shape_value()
				      * function. If the finite
				      * element is vector valued but
				      * all shape functions are
				      * primitive (i.e. they are
				      * non-zero in only one
				      * component), then the value
				      * returned by shape_value()
				      * equals that of this function
				      * for exactly one
				      * component. This function is
				      * therefore only of greater
				      * interest if the shape function
				      * is not primitive, but then it
				      * is necessary since the other
				      * function cannot be used.
				      *
				      * @arg function_no Number
				      * of the shape function to be
				      * computed
				      * @arg point_no Number of
				      * the quadrature point at which
				      * function is to be computed
				      * @arg component vector component to be computed
				      */
    double shape_value_component (const unsigned int function_no,
				  const unsigned int point_no,
				  const unsigned int component) const;

    				     /**
				      * Compute the gradient of the
				      * <tt>i</tt>th shape function at the
				      * <tt>j</tt>th quadrature point with
				      * respect to real cell
				      * coordinates.  If you want to
				      * get the derivative in one of
				      * the coordinate directions, use
				      * the appropriate function of
				      * the Tensor class to
				      * extract one component. Since
				      * only a reference to the
				      * gradient's value is returned,
				      * there should be no major
				      * performance drawback.
				      *
				      * If the shape function is
				      * vector-valued, then this
				      * returns the only non-zero
				      * component. If the shape
				      * function has more than one
				      * non-zero component (i.e. it is
				      * not primitive), then throw an
				      * exception of type
				      * ExcShapeFunctionNotPrimitive. In
				      * that case, use the
				      * shape_grad_component()
				      * function.
				      */
    const Tensor<1,spacedim> &
    shape_grad (const unsigned int function,
		const unsigned int quadrature_point) const;

				     /**
				      * Return one vector component of
				      * the gradient of a shape function
				      * at a quadrature point. If the
				      * finite element is scalar, then
				      * only component zero is allowed
				      * and the return value equals
				      * that of the shape_grad()
				      * function. If the finite
				      * element is vector valued but
				      * all shape functions are
				      * primitive (i.e. they are
				      * non-zero in only one
				      * component), then the value
				      * returned by shape_grad()
				      * equals that of this function
				      * for exactly one
				      * component. This function is
				      * therefore only of greater
				      * interest if the shape function
				      * is not primitive, but then it
				      * is necessary since the other
				      * function cannot be used.
				      */
    Tensor<1,spacedim>
    shape_grad_component (const unsigned int function_no,
			  const unsigned int point_no,
			  const unsigned int component) const;

    				     /**
				      * Second derivatives of
				      * the <tt>function_no</tt>th shape function at
				      * the <tt>point_no</tt>th quadrature point
				      * with respect to real cell
				      * coordinates. If you want to
				      * get the derivatives in one of
				      * the coordinate directions, use
				      * the appropriate function of
				      * the Tensor class to
				      * extract one component. Since
				      * only a reference to the
				      * derivative values is returned,
				      * there should be no major
				      * performance drawback.
				      *
				      * If the shape function is
				      * vector-valued, then this
				      * returns the only non-zero
				      * component. If the shape
				      * function has more than one
				      * non-zero component (i.e. it is
				      * not primitive), then throw an
				      * exception of type
				      * ExcShapeFunctionNotPrimitive. In
				      * that case, use the
				      * shape_grad_grad_component()
				      * function.
				      */
    const Tensor<2,spacedim> &
    shape_hessian (const unsigned int function_no,
		   const unsigned int point_no) const;

				     /**
				      * @deprecated Wrapper for shape_hessian()
				      */
    const Tensor<2,spacedim> &
    shape_2nd_derivative (const unsigned int function_no,
			  const unsigned int point_no) const;


				     /**
				      * Return one vector component of
				      * the gradient of a shape
				      * function at a quadrature
				      * point. If the finite element
				      * is scalar, then only component
				      * zero is allowed and the return
				      * value equals that of the
				      * shape_hessian()
				      * function. If the finite
				      * element is vector valued but
				      * all shape functions are
				      * primitive (i.e. they are
				      * non-zero in only one
				      * component), then the value
				      * returned by
				      * shape_hessian()
				      * equals that of this function
				      * for exactly one
				      * component. This function is
				      * therefore only of greater
				      * interest if the shape function
				      * is not primitive, but then it
				      * is necessary since the other
				      * function cannot be used.
				      */
    Tensor<2,spacedim>
    shape_hessian_component (const unsigned int function_no,
			     const unsigned int point_no,
			     const unsigned int component) const;
    
				     /**
				      * @deprecated Wrapper for shape_hessian_component()
				      */
    Tensor<2,spacedim>
    shape_2nd_derivative_component (const unsigned int function_no,
				    const unsigned int point_no,
				    const unsigned int component) const;
    

				     //@}
				     /// @name FunctionAccess Access to values of global finite element functions
				     //@{
    
				     /**
				      * Returns the values of the
				      * finite element function
				      * characterized by
				      * <tt>fe_function</tt> restricted to
				      * the cell, face or subface
				      * selected the last time the
				      * <tt>reinit</tt> function of the
				      * derived class was called, at
				      * the quadrature points.
				      *
				      * If the present cell is not an
				      * active one the interpolated
				      * function values are returned.
				      *
				      * To get values of
				      * multi-component elements,
				      * there is another
				      * get_function_values() below,
				      * returning a vector of vectors
				      * of results.
				      *
				      * This function may only be used if the
				      * finite element in use is a scalar one,
				      * i.e. has only one vector component. If
				      * it is a vector-valued one, then use
				      * the other get_function_values()
				      * function.
				      * 
				      * The function assumes that the
				      * <tt>values</tt> object already has the
				      * correct size. 
				      *
				      * The actual data type of the
				      * input vector may be either a
				      * Vector&lt;T&gt;,
				      * BlockVector&lt;T&gt;, or one
				      * of the sequential PETSc or
				      * Trilinos vector wrapper
				      * classes. It represents a
				      * global vector of DoF values
				      * associated with the DofHandler
				      * object with which this
				      * FEValues object was last
				      * initialized.
				      */
    template <class InputVector, typename number>
    void get_function_values (const InputVector& fe_function,
			      std::vector<number>& values) const;

				     /**
				      * Access to vector valued finite
				      * element functions.
				      *
				      * This function does the same as
				      * the other
				      * get_function_values(), but
				      * applied to multi-component
				      * elements.
				      *
				      * The actual data type of the
				      * input vector may be either a
				      * Vector&lt;T&gt;,
				      * BlockVector&lt;T&gt;, or one
				      * of the sequential PETSc or
				      * Trilinos vector wrapper
				      * classes. It represents a
				      * global vector of DoF values
				      * associated with the DofHandler
				      * object with which this
				      * FEValues object was last
				      * initialized.
				      */
    template <class InputVector, typename number>
    void get_function_values (const InputVector       &fe_function,
			      std::vector<Vector<number> > &values) const;

				     /**
				      * Generate function values from
				      * an arbitrary vector.
				      * 
				      * This function offers the
				      * possibility to extract
				      * function values in quadrature
				      * points from vectors not
				      * corresponding to a whole
				      * discretization.
				      *
				      * The vector <tt>indices</tt>
				      * corresponds to the degrees of
				      * freedom on a single cell. Its
				      * length may even be a multiple
				      * of the number of dofs per
				      * cell. Then, the vectors in
				      * <tt>value</tt> should allow
				      * for the same multiple of the
				      * components of the finite
				      * element.
				      *
				      * You may want to use this
				      * function, if you want to
				      * access just a single block
				      * from a BlockVector, if you
				      * have a multi-level vector or
				      * if you already have a local
				      * representation of your finite
				      * element data.
				      */
    template <class InputVector, typename number>
    void get_function_values (const InputVector& fe_function,
			      const VectorSlice<const std::vector<unsigned int> >& indices,
			      std::vector<number>& values) const;

				     /**
				      * Generate vector function
				      * values from an arbitrary
				      * vector.
				      *
				      * This function offers the
				      * possibility to extract
				      * function values in quadrature
				      * points from vectors not
				      * corresponding to a whole
				      * discretization.
				      *
				      * The vector <tt>indices</tt>
				      * corresponds to the degrees of
				      * freedom on a single cell. Its
				      * length may even be a multiple
				      * of the number of dofs per
				      * cell. Then, the vectors in
				      * <tt>value</tt> should allow
				      * for the same multiple of the
				      * components of the finite
				      * element.
				      *
				      * You may want to use this
				      * function, if you want to
				      * access just a single block
				      * from a BlockVector, if you
				      * have a multi-level vector or
				      * if you already have a local
				      * representation of your finite
				      * element data.
				      *
				      * Since this function allows for
				      * fairly general combinations of
				      * argument sizes, be aware that
				      * the checks on the arguments
				      * may not detect errors.
				      */
    template <class InputVector, typename number>
    void get_function_values (const InputVector& fe_function,
			      const VectorSlice<const std::vector<unsigned int> >& indices,
			      std::vector<Vector<number> >& values) const;


				     /**
				      * Generate vector function
				      * values from an arbitrary
				      * vector.
				      *
				      * This function offers the
				      * possibility to extract
				      * function values in quadrature
				      * points from vectors not
				      * corresponding to a whole
				      * discretization.
				      *
				      * The vector <tt>indices</tt>
				      * corresponds to the degrees of
				      * freedom on a single cell. Its
				      * length may even be a multiple
				      * of the number of dofs per
				      * cell. Then, the vectors in
				      * <tt>value</tt> should allow
				      * for the same multiple of the
				      * components of the finite
				      * element.
				      *
				      * Depending on the value of the last
				      * argument, the outer vector of
				      * <tt>values</tt> has either the
				      * length of the quadrature rule
				      * (<tt>quadrature_points_fastest
				      * == false</tt>) or the length
				      * of components to be filled
				      * <tt>quadrature_points_fastest
				      * == true</tt>. If <tt>p</tt> is
				      * the current quadrature point
				      * number and <tt>i</tt> is the
				      * vector component of the
				      * solution desired, the access
				      * to <tt>values</tt> is
				      * <tt>values[p][i]</tt> if
				      * <tt>quadrature_points_fastest
				      * == false</tt>, and
				      * <tt>values[i][p]</tt>
				      * otherwise.
				      *
				      * You may want to use this
				      * function, if you want to
				      * access just a single block
				      * from a BlockVector, if you
				      * have a multi-level vector or
				      * if you already have a local
				      * representation of your finite
				      * element data.
				      *
				      * Since this function allows for
				      * fairly general combinations of
				      * argument sizes, be aware that
				      * the checks on the arguments
				      * may not detect errors.
				      */
    template <class InputVector, typename number>
    void get_function_values (const InputVector& fe_function,
			      const VectorSlice<const std::vector<unsigned int> >& indices,
			      std::vector<std::vector<number> >& values,
			      const bool quadrature_points_fastest) const;

				     /**
				      * Compute the gradients of the finite
				      * element function characterized
				      * by @p fe_function restricted to
				      * @p cell at the quadrature points.
				      *
				      * If the present cell is not an active
				      * one the interpolated function values
				      * are returned.
				      *
				      * This function may only be used if the
				      * finite element in use is a scalar one,
				      * i.e. has only one vector component. If
				      * it is a vector-valued one, then use
				      * the other get_function_gradients()
				      * function.
				      * 
				      * The function assumes that the
				      * @p gradients object already has the
				      * right size.
				      *
				      * The actual data type of the
				      * input vector may be either a
				      * Vector&lt;T&gt;,
				      * BlockVector&lt;T&gt;, or one
				      * of the sequential PETSc or
				      * Trilinos vector wrapper
				      * classes. It represents a
				      * global vector of DoF values
				      * associated with the DofHandler
				      * object with which this
				      * FEValues object was last
				      * initialized.
				      *
				      * The output are the gradients
				      * of the function represented by
				      * these DoF values, as computed
				      * in real space (as opposed to
				      * on the unit cell).
				      */
    template <class InputVector>
    void get_function_gradients (const InputVector      &fe_function,
				 std::vector<Tensor<1,spacedim> > &gradients) const;

				     /**
				      * @deprecated Use
				      * get_function_gradients() instead.
				      */
    template <class InputVector>
    void get_function_grads (const InputVector      &fe_function,
			     std::vector<Tensor<1,spacedim> > &gradients) const;

				     /**
				      * Compute the gradients of the finite
				      * element function characterized
				      * by @p fe_function restricted to
				      * @p cell at the quadrature points.
				      *
				      * If the present cell is not an active
				      * one the interpolated function values
				      * are returned.
				      *
				      * The function assumes that the
				      * @p gradients object already has the
				      * right size.
				      *
				      * This function does the same as
				      * the other get_function_values(),
				      * but applied to multi-component
				      * elements.
				      *
				      * The actual data type of the
				      * input vector may be either a
				      * Vector&lt;T&gt;,
				      * BlockVector&lt;T&gt;, or one
				      * of the sequential PETSc or
				      * Trilinos vector wrapper
				      * classes. It represents a
				      * global vector of DoF values
				      * associated with the DofHandler
				      * object with which this
				      * FEValues object was last
				      * initialized.
				      *
				      * The output are the gradients
				      * of the function represented by
				      * these DoF values, as computed
				      * in real space (as opposed to
				      * on the unit cell).
				      */
    template <class InputVector>
    void get_function_gradients (const InputVector               &fe_function,
				 std::vector<std::vector<Tensor<1,spacedim> > > &gradients) const;

    
				     /**
				      * @deprecated Use
				      * get_function_gradients() instead.
				      */
    template <class InputVector>
    void get_function_grads (const InputVector               &fe_function,
			     std::vector<std::vector<Tensor<1,spacedim> > > &gradients) const;

				     /**
				      * Function gradient access with
				      * more flexibility. see
				      * get_function_values() with
				      * corresponding arguments.
				      */
    template <class InputVector>
    void get_function_gradients (const InputVector& fe_function,
				 const VectorSlice<const std::vector<unsigned int> >& indices,
				 std::vector<Tensor<1,spacedim> >& gradients) const;

				     /**
				      * @deprecated Use
				      * get_function_gradients() instead.
				      */
    template <class InputVector>
    void get_function_grads (const InputVector& fe_function,
			     const VectorSlice<const std::vector<unsigned int> >& indices,
			     std::vector<Tensor<1,spacedim> >& gradients) const;

				     /**
				      * Function gradient access with
				      * more flexibility. see
				      * get_function_values() with
				      * corresponding arguments.
				      */
    template <class InputVector>
    void get_function_gradients (const InputVector& fe_function,
				 const VectorSlice<const std::vector<unsigned int> >& indices,
				 std::vector<std::vector<Tensor<1,spacedim> > >& gradients,
				 bool quadrature_points_fastest = false) const;

				     /**
				      * @deprecated Use
				      * get_function_gradients() instead.
				      */
    template <class InputVector>
    void get_function_grads (const InputVector& fe_function,
			     const VectorSlice<const std::vector<unsigned int> >& indices,
			     std::vector<std::vector<Tensor<1,spacedim> > >& gradients,
			     bool quadrature_points_fastest = false) const;

				     /**
				      * Compute the tensor of second
				      * derivatives of the finite
				      * element function characterized
				      * by @p fe_function restricted
				      * to @p cell at the quadrature
				      * points.
				      *
				      * The function assumes that the
				      * @p hessians object
				      * already has the correct size.
				      *
				      * This function may only be used if the
				      * finite element in use is a scalar one,
				      * i.e. has only one vector component. If
				      * it is a vector-valued one, then use
				      * the other
				      * get_function_hessians()
				      * function.
				      * 
				      * The actual data type of the
				      * input vector may be either a
				      * Vector&lt;T&gt;,
				      * BlockVector&lt;T&gt;, or one
				      * of the sequential PETSc or
				      * Trilinos vector wrapper
				      * classes. It represents a
				      * global vector of DoF values
				      * associated with the DofHandler
				      * object with which this
				      * FEValues object was last
				      * initialized.
				      *
				      * The output are the second
				      * derivatives of the function
				      * represented by these DoF
				      * values, as computed in real
				      * space (as opposed to on the
				      * unit cell).
				      */
    template <class InputVector>
    void
    get_function_hessians (const InputVector& fe_function,
			   std::vector<Tensor<2,spacedim> >& hessians) const;

    
				     /**
				      * Compute the tensor of second
				      * derivatives of the finite
				      * element function characterized
				      * by @p fe_function restricted to
				      * @p cell at the quadrature points.
				      *
				      * The function assumes that the
				      * @p hessians object already has
				      * the right size.
				      *
				      * This function does the same as
				      * the other one with the same
				      * name, but applies to
				      * vector-valued finite elements.
				      *
				      * The actual data type of the
				      * input vector may be either a
				      * Vector&lt;T&gt;,
				      * BlockVector&lt;T&gt;, or one
				      * of the sequential PETSc or
				      * Trilinos vector wrapper
				      * classes. It represents a
				      * global vector of DoF values
				      * associated with the DofHandler
				      * object with which this
				      * FEValues object was last
				      * initialized.
				      *
				      * The output are the second derivatives
				      * of the function represented by
				      * these DoF values, as computed
				      * in real space (as opposed to
				      * on the unit cell).
				      */
    template <class InputVector>
    void
    get_function_hessians (const InputVector      &fe_function,
			   std::vector<std::vector<Tensor<2,spacedim> > > &hessians,
			   bool quadrature_points_fastest = false) const;

				     /**
				      * Access to the second
				      * derivatives of a function with
				      * more flexibility. see
				      * get_function_values() with
				      * corresponding arguments.
				      */
    template <class InputVector>
    void get_function_hessians (
      const InputVector& fe_function,
      const VectorSlice<const std::vector<unsigned int> >& indices,
      std::vector<Tensor<2,spacedim> >& hessians) const;

				     /**
				      * Access to the second
				      * derivatives of a function with
				      * more flexibility. see
				      * get_function_values() with
				      * corresponding arguments.
				      */
    template <class InputVector>
    void get_function_hessians (
      const InputVector& fe_function,
      const VectorSlice<const std::vector<unsigned int> >& indices,
      std::vector<std::vector<Tensor<2,spacedim> > >& hessians,
      bool quadrature_points_fastest = false) const;

				     /**
				      * @deprecated Wrapper for get_function_hessians()
				      */
    template <class InputVector>
    void
    get_function_2nd_derivatives (const InputVector&,
				  std::vector<Tensor<2,spacedim> >&) const;

				     /**
				      * @deprecated Wrapper for get_function_hessians()
				      */
    template <class InputVector>
    void
    get_function_2nd_derivatives (const InputVector&,
				  std::vector<std::vector<Tensor<2,spacedim> > >&,
				  bool = false) const;
    

				     /**
				      * Compute the (scalar) Laplacian
				      * of the finite element function
				      * characterized by @p
				      * fe_function restricted to @p
				      * cell at the quadrature
				      * points. The Laplacian output
				      * vector is equivalent to
				      * getting
				      * <tt>trace(hessians)</tt>,
				      * where <tt>hessian</tt> would
				      * be the output of the
				      * get_function_hessians()
				      * function.
				      *
				      * The function assumes that the
				      * @p laplacians object
				      * already has the correct size.
				      *
				      * This function may only be used if the
				      * finite element in use is a scalar one,
				      * i.e. has only one vector component. If
				      * it is a vector-valued one, then use
				      * the other
				      * get_function_laplacians()
				      * function.
				      * 
				      * The actual data type of the
				      * input vector may be either a
				      * Vector&lt;T&gt;,
				      * BlockVector&lt;T&gt;, or one
				      * of the sequential PETSc or
				      * Trilinos vector wrapper
				      * classes. It represents a
				      * global vector of DoF values
				      * associated with the DofHandler
				      * object with which this
				      * FEValues object was last
				      * initialized.
				      *
				      * The output are the traces of
				      * the second derivatives
				      * (i.e. Laplacians) of the
				      * function represented by these
				      * DoF values, as computed in
				      * real space (as opposed to on
				      * the unit cell).
				      */
    template <class InputVector, typename number>
    void
    get_function_laplacians (const InputVector& fe_function,
			     std::vector<number>& laplacians) const;

    
				     /**
				      * Compute the (scalar) Laplacian
				      * of the finite element function
				      * characterized by @p
				      * fe_function restricted to @p
				      * cell at the quadrature
				      * points. The Laplacian output
				      * vector is equivalent to
				      * getting
				      * <tt>trace(hessians)</tt>, with
				      * <tt>hessian</tt> corresponding
				      * to the output of the
				      * get_function_hessians()
				      * function.
				      *
				      * The function assumes that the
				      * @p laplacians object
				      * already has the correct size.
				      *
				      * This function does the same as
				      * the other one with the same
				      * name, but applies to
				      * vector-valued finite elements.
				      *
				      * The actual data type of the
				      * input vector may be either a
				      * Vector&lt;T&gt;,
				      * BlockVector&lt;T&gt;, or one
				      * of the sequential PETSc or
				      * Trilinos vector wrapper
				      * classes. It represents a
				      * global vector of DoF values
				      * associated with the DofHandler
				      * object with which this
				      * FEValues object was last
				      * initialized.
				      *
				      * The output are the traces of
				      * the second derivatives (i.e.
				      * Laplacians) of the function
				      * represented by these DoF
				      * values, as computed in real
				      * space (as opposed to on the
				      * unit cell).
				      */
    template <class InputVector, typename number>
    void
    get_function_laplacians (const InputVector      &fe_function,
			     std::vector<Vector<number> > &laplacians) const;

				     /**
				      * Access to the second
				      * derivatives of a function with
				      * more flexibility. see
				      * get_function_values() with
				      * corresponding arguments.
				      */
    template <class InputVector, typename number>
    void get_function_laplacians (
      const InputVector& fe_function,
      const VectorSlice<const std::vector<unsigned int> >& indices,
      std::vector<number>& laplacians) const;

				     /**
				      * Access to the second
				      * derivatives of a function with
				      * more flexibility. see
				      * get_function_values() with
				      * corresponding arguments.
				      */
    template <class InputVector, typename number>
    void get_function_laplacians (
      const InputVector& fe_function,
      const VectorSlice<const std::vector<unsigned int> >& indices,
      std::vector<Vector<number> >& laplacians) const;

				     /**
				      * Access to the second
				      * derivatives of a function with
				      * more flexibility. see
				      * get_function_values() with
				      * corresponding arguments.
				      */
    template <class InputVector, typename number>
    void get_function_laplacians (
      const InputVector& fe_function,
      const VectorSlice<const std::vector<unsigned int> >& indices,
      std::vector<std::vector<number> >& laplacians,
      bool quadrature_points_fastest = false) const;


				     //@}
    
				     /**
				      * Position of the <tt>i</tt>th
				      * quadrature point in real space.
				      */
    const Point<spacedim> & quadrature_point (const unsigned int i) const;

				     /**
				      * Return a pointer to the vector of
				      * quadrature points.
				      */
    const std::vector<Point<spacedim> > & get_quadrature_points () const;

				     /**
				      * Mapped quadrature weight. If
				      * this object refers to a volume
				      * evaluation (i.e. the derived
				      * class is of type FEValues),
				      * then this is the Jacobi
				      * determinant times the weight
				      * of the *<tt>i</tt>th unit
				      * quadrature point.
				      *
				      * For surface evaluations
				      * (i.e. classes FEFaceValues or
				      * FESubfaceValues), it is the
				      * mapped surface element times
				      * the weight of the quadrature
				      * point.
				      *
				      * You can think of the quantity returned
				      * by this function as the volume or
				      * surface element $dx, ds$ in the
				      * integral that we implement here by
				      * quadrature.
				      */
    double JxW (const unsigned int quadrature_point) const;

				     /**
				      * Pointer to the array holding
				      * the values returned by JxW().
				      */
    const std::vector<double> & get_JxW_values () const;
    
				     /**
				      * Return the Jacobian of the
				      * transformation at the specified
				      * quadrature point, i.e.
				      * $J_{ij}=dx_i/d\hat x_j$
				      */
    const Tensor<2,spacedim> & jacobian (const unsigned int quadrature_point) const;

				     /**
				      * Pointer to the array holding
				      * the values returned by jacobian().
				      */
    const std::vector<Tensor<2,spacedim> > & get_jacobians () const;

    				     /**
				      * Return the second derivative of the
				      * transformation from unit to real cell,
				      * i.e. the first derivative of the
				      * Jacobian, at the specified quadrature
				      * point, i.e. $G_{ijk}=dJ_{jk}/d\hat x_i$.
				      */
    const Tensor<3,spacedim> & jacobian_grad (const unsigned int quadrature_point) const;

				     /**
				      * Pointer to the array holding
				      * the values returned by
				      * jacobian_grads().
				      */
    const std::vector<Tensor<3,spacedim> > & get_jacobian_grads () const;

				     /**
				      * Return the inverse Jacobian of the
				      * transformation at the specified
				      * quadrature point, i.e.
				      * $J_{ij}=d\hat x_i/dx_j$
				      */
    const Tensor<2,spacedim> & inverse_jacobian (const unsigned int quadrature_point) const;

				     /**
				      * Pointer to the array holding
				      * the values returned by 
				      * inverse_jacobian().
				      */
    const std::vector<Tensor<2,spacedim> > & get_inverse_jacobians () const;

				     /**
				      * Constant reference to the
				      * selected mapping object.
				      */
    const Mapping<dim,spacedim> & get_mapping () const;

				     /**
				      * Constant reference to the
				      * selected finite element
				      * object.
				      */
    const FiniteElement<dim,spacedim> & get_fe () const;

				     /**
				      * Return the update flags set
				      * for this object.
				      */
    UpdateFlags get_update_flags () const;
    
				     /**
				      * Return a triangulation
				      * iterator to the current cell.
				      */
    const typename Triangulation<dim,spacedim>::cell_iterator get_cell () const;

				     /**
				      * Returns the vectors normal to
				      * the cell in each of the
				      * quadrature points.
				      */

    const std::vector<Point<spacedim> > & get_cell_normal_vectors () const;

    				     /**
				      * Return the outward normal vector to
				      * the cell at the <tt>i</tt>th quadrature
				      * point. The length of the vector
				      * is normalized to one.
				      */
    const Point<spacedim> & cell_normal_vector (const unsigned int i) const;

				     /**
				      * Return the relation of the current
				      * cell to the previous cell. This
				      * allows re-use of some cell data
				      * (like local matrices for equations
				      * with constant coefficients) if the
				      * result is
				      * <tt>CellSimilarity::translation</tt>.
				      */
    CellSimilarity::Similarity get_cell_similarity () const;

				     /**
				      * Determine an estimate for the
				      * memory consumption (in bytes)
				      * of this object.
				      */
    unsigned int memory_consumption () const;

				     /**
				      * This exception is thrown if
				      * FEValuesBase is asked to
				      * return the value of a field
				      * which was not required by the
				      * UpdateFlags for this
				      * FEValuesBase.
				      *
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcAccessToUninitializedField);
				     /**
				      * @todo Document this
				      *
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcCannotInitializeField);
				     /**
				      * @todo Document this
				      *
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcInvalidUpdateFlag);
				     /**
				      * @todo Document this
				      *
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcFEDontMatch);
				     /**
				      * @todo Document this
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
				      * @todo Document this
				      *
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcFENotPrimitive);
    
  protected:
                                     /**
                                      * Objects of the FEValues
                                      * class need to store a pointer
                                      * (i.e. an iterator) to the
                                      * present cell in order to be
                                      * able to extract the values of
                                      * the degrees of freedom on this
                                      * cell in the
                                      * get_function_values() and
                                      * assorted functions. On the
                                      * other hand, this class should
                                      * also work for different
                                      * iterators, as long as they
                                      * have the same interface to
                                      * extract the DoF values (i.e.,
                                      * for example, they need to have
                                      * a @p get_interpolated_dof_values
                                      * function).
                                      *
                                      * This calls for a common base
                                      * class of iterator classes, and
                                      * making the functions we need
                                      * here @p virtual. On the other
                                      * hand, this is the only place
                                      * in the library where we need
                                      * this, and introducing a base
                                      * class of iterators and making
                                      * a function virtual penalizes
                                      * <em>all</em> users of the
                                      * iterators, which are basically
                                      * intended as very fast accessor
                                      * functions. So we do not want
                                      * to do this. Rather, what we do
                                      * here is making the functions
                                      * we need virtual only for use
                                      * with <em>this class</em>. The idea
                                      * is the following: have a
                                      * common base class which
                                      * declares some pure virtual
                                      * functions, and for each
                                      * possible iterator type, we
                                      * have a derived class which
                                      * stores the iterator to the
                                      * cell and implements these
                                      * functions. Since the iterator
                                      * classes have the same
                                      * interface, we can make the
                                      * derived classes a template,
                                      * templatized on the iterator
                                      * type.
                                      *
                                      * This way, the use of virtual
                                      * functions is restricted to
                                      * only this class, and other
                                      * users of iterators do not have
                                      * to bear the negative effects.
                                      *
                                      * @author Wolfgang Bangerth, 2003
                                      */
    class CellIteratorBase 
    {
      public:
                                         /**
                                          * Destructor. Made virtual
                                          * since we store only
                                          * pointers to the base
                                          * class.
                                          */
        virtual ~CellIteratorBase ();
        
                                         /**
                                          * Conversion operator to an
                                          * iterator for
                                          * triangulations. This
                                          * conversion is implicit for
                                          * the original iterators,
                                          * since they are derived
                                          * classes. However, since
                                          * here we have kind of a
                                          * parallel class hierarchy,
                                          * we have to have a
                                          * conversion operator.
                                          */
        virtual
        operator const typename Triangulation<dim,spacedim>::cell_iterator () const = 0;

                                         /**
                                          * Return the number of
                                          * degrees of freedom the DoF
                                          * handler object has to
                                          * which the iterator belongs
                                          * to.
                                          */
        virtual
        unsigned int
        n_dofs_for_dof_handler () const = 0;

                                         /**
                                          * Call
                                          * @p get_interpolated_dof_values
                                          * of the iterator with the
                                          * given arguments.
                                          */
        virtual
        void
        get_interpolated_dof_values (const Vector<double> &in,
                                     Vector<double>       &out) const = 0;

                                         /**
                                          * Call
                                          * @p get_interpolated_dof_values
                                          * of the iterator with the
                                          * given arguments.
                                          */
        virtual
        void
        get_interpolated_dof_values (const Vector<float> &in,
                                     Vector<float>       &out) const = 0;

                                         /**
                                          * Call
                                          * @p get_interpolated_dof_values
                                          * of the iterator with the
                                          * given arguments.
                                          */
        virtual
        void
        get_interpolated_dof_values (const Vector<long double> &in,
                                     Vector<long double>       &out) const = 0;

                                         /**
                                          * Call
                                          * @p get_interpolated_dof_values
                                          * of the iterator with the
                                          * given arguments.
                                          */
        virtual
        void
        get_interpolated_dof_values (const BlockVector<double> &in,
                                     Vector<double>       &out) const = 0;

                                         /**
                                          * Call
                                          * @p get_interpolated_dof_values
                                          * of the iterator with the
                                          * given arguments.
                                          */
        virtual
        void
        get_interpolated_dof_values (const BlockVector<float> &in,
                                     Vector<float>       &out) const = 0;

                                         /**
                                          * Call
                                          * @p get_interpolated_dof_values
                                          * of the iterator with the
                                          * given arguments.
                                          */
        virtual
        void
        get_interpolated_dof_values (const BlockVector<long double> &in,
                                     Vector<long double>       &out) const = 0;

#ifdef DEAL_II_USE_PETSC
                                         /**
                                          * Call
                                          * @p get_interpolated_dof_values
                                          * of the iterator with the
                                          * given arguments.
                                          */
        virtual
        void
        get_interpolated_dof_values (const PETScWrappers::Vector &in,
                                     Vector<PetscScalar>         &out) const = 0;

                                         /**
                                          * Call
                                          * @p get_interpolated_dof_values
                                          * of the iterator with the
                                          * given arguments.
                                          */
        virtual
        void
        get_interpolated_dof_values (const PETScWrappers::BlockVector &in,
                                     Vector<PetscScalar>              &out) const = 0;
#endif

#ifdef DEAL_II_USE_TRILINOS
                                         /**
                                          * Call
                                          * @p get_interpolated_dof_values
                                          * of the iterator with the
                                          * given arguments.
                                          */
        virtual
        void
        get_interpolated_dof_values (const TrilinosWrappers::Vector &in,
                                     Vector<TrilinosScalar>         &out) const = 0;

                                         /**
                                          * Call
                                          * @p get_interpolated_dof_values
                                          * of the iterator with the
                                          * given arguments.
                                          */
        virtual
        void
        get_interpolated_dof_values (const TrilinosWrappers::BlockVector &in,
                                     Vector<TrilinosScalar>              &out) const = 0;
                                         /**
                                          * Call
                                          * @p get_interpolated_dof_values
                                          * of the iterator with the
                                          * given arguments.
                                          */

        virtual
        void
        get_interpolated_dof_values (const TrilinosWrappers::MPI::Vector &in,
                                     Vector<TrilinosScalar>              &out) const = 0;

                                         /**
                                          * Call
                                          * @p get_interpolated_dof_values
                                          * of the iterator with the
                                          * given arguments.
                                          */
        virtual
        void
        get_interpolated_dof_values (const TrilinosWrappers::MPI::BlockVector &in,
                                     Vector<TrilinosScalar>                   &out) const = 0;
#endif
    };

				     /**
				      * Forward declaration of classes derived
				      * from CellIteratorBase. Their
				      * definition and implementation is given
				      * in the .cc file.
				      */
    template <typename CI> class CellIterator;
    class TriaCellIterator;
    
				     /**
				      * Store the cell selected last time the
				      * reinit() function was called.  This is
				      * necessary for the
				      * <tt>get_function_*</tt> functions as
				      * well as the functions of same name in
				      * the extractor classes.
				      */
    std::auto_ptr<const CellIteratorBase> present_cell;
    
				     /**
				      * Storage for the mapping object.
				      */
    const SmartPointer<const Mapping<dim,spacedim> > mapping;
    
				     /**
				      * Store the finite element for later use.
				      */
    const SmartPointer<const FiniteElement<dim,spacedim> > fe;

    
				     /**
				      * Internal data of mapping.
				      */
    SmartPointer<typename Mapping<dim,spacedim>::InternalDataBase> mapping_data;

				     /**
				      * Internal data of finite element.
				      */
    SmartPointer<typename Mapping<dim,spacedim>::InternalDataBase> fe_data;

				     /**
				      * Initialize some update
				      * flags. Called from the
				      * @p initialize functions of
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
				      * An enum variable that can store
				      * different states of the current cell
				      * in comparison to the previously
				      * visited cell. If wanted, additional
				      * states can be checked here and used
				      * in one of the methods used during
				      * reinit.
				      */
    CellSimilarity::Similarity cell_similarity;

				     /**
				      * A function that checks whether the
				      * new cell is similar to the one
				      * previously used. Then, a significant
				      * amount of the data can be reused,
				      * e.g. the derivatives of the basis
				      * functions in real space, shape_grad.
				      */
    void
    check_cell_similarity (const typename Triangulation<dim,spacedim>::cell_iterator &cell);

  private:
                                     /**
                                      * Copy constructor. Since
                                      * objects of this class are not
                                      * copyable, we make it private,
                                      * and also do not implement it.
                                      */
    FEValuesBase (const FEValuesBase &);
    
                                     /**
                                      * Copy operator. Since
                                      * objects of this class are not
                                      * copyable, we make it private,
                                      * and also do not implement it.
                                      */
    FEValuesBase & operator= (const FEValuesBase &);

				     /**
				      * A cache for all possible FEValuesViews
				      * objects.
				      */
    internal::FEValuesViews::Cache<dim,spacedim> fe_values_views_cache;
    
				     /**
				      * Make the view classes friends of this
				      * class, since they access internal
				      * data.
				      */
    template <int, int> friend class FEValuesViews::Scalar;
    template <int, int> friend class FEValuesViews::Vector;
};



/**
 * Finite element evaluated in quadrature points of a cell.
 *
 * This function implements the initialization routines for
 * FEValuesBase, if values in quadrature points of a cell are
 * needed. For further documentation see this class.
 * 
 * @author Wolfgang Bangerth, 1998, Guido Kanschat, 2001
 */
template <int dim, int spacedim=dim>
class FEValues : public FEValuesBase<dim,spacedim>
{
  public:
                                     /**
                                      * Dimension in which this object
                                      * operates.
                                      */
    static const unsigned int dimension = dim;

                                     /**
                                      * Dimension of the space in
                                      * which this object operates.
                                      */
    static const unsigned int space_dimension = spacedim;

                                     /**
                                      * Dimension of the object over
                                      * which we integrate. For the
                                      * present class, this is equal
                                      * to <code>dim</code>.
                                      */
    static const unsigned int integral_dimension = dim;
    
				     /**
				      * Constructor. Gets cell
				      * independent data from mapping
				      * and finite element objects,
				      * matching the quadrature rule
				      * and update flags.
				      */
    FEValues (const Mapping<dim,spacedim>       &mapping,
	      const FiniteElement<dim,spacedim> &fe,
	      const Quadrature<dim>             &quadrature,
	      const UpdateFlags                  update_flags);

                                     /**
				      * Constructor. Uses MappingQ1
				      * implicitly.
				      */
    FEValues (const FiniteElement<dim,spacedim> &fe,
	      const Quadrature<dim>             &quadrature,
	      const UpdateFlags                  update_flags);
    
				     /**
				      * Reinitialize the gradients,
				      * Jacobi determinants, etc for
				      * the given cell of type
				      * "iterator into a DoFHandler
				      * object", and the finite
				      * element associated with this
				      * object. It is assumed that the
				      * finite element used by the
				      * given cell is also the one
				      * used by this FEValues
				      * object.
				      */
    void reinit (const typename DoFHandler<dim,spacedim>::cell_iterator &cell);

				     /**
				      * Reinitialize the gradients,
				      * Jacobi determinants, etc for
				      * the given cell of type
				      * "iterator into a hp::DoFHandler
				      * object", and the finite
				      * element associated with this
				      * object. It is assumed that the
				      * finite element used by the
				      * given cell is also the one
				      * used by this @p FEValues
				      * object.
				      */
    void reinit (const typename hp::DoFHandler<dim,spacedim>::cell_iterator &cell);

				     /**
				      * Reinitialize the gradients,
				      * Jacobi determinants, etc for
				      * the given cell of type
				      * "iterator into a MGDoFHandler
				      * object", and the finite
				      * element associated with this
				      * object. It is assumed that the
				      * finite element used by the
				      * given cell is also the one
				      * used by this FEValues
				      * object.
				      */
    void reinit (const typename MGDoFHandler<dim,spacedim>::cell_iterator &cell);

				     /**
				      * Reinitialize the gradients,
				      * Jacobi determinants, etc for
				      * the given cell of type
				      * "iterator into a Triangulation
				      * object", and the given finite
				      * element. Since iterators into
				      * triangulation alone only
				      * convey information about the
				      * geometry of a cell, but not
				      * about degrees of freedom
				      * possibly associated with this
				      * cell, you will not be able to
				      * call some functions of this
				      * class if they need information
				      * about degrees of
				      * freedom. These functions are,
				      * above all, the
				      * <tt>get_function_value/gradients/hessians/laplacians</tt>
				      * functions. If you want to call
				      * these functions, you have to
				      * call the @p reinit variants
				      * that take iterators into
				      * DoFHandler or other DoF
				      * handler type objects.
				      */
    void reinit (const typename Triangulation<dim,spacedim>::cell_iterator &cell);

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

                                     /**
                                      * Return a reference to this
                                      * very object.
                                      *
                                      * Though it seems that it is not very
                                      * useful, this function is there to
                                      * provide capability to the hpFEValues
                                      * class, in which case it provides the
                                      * FEValues object for the present cell
                                      * (remember that for hp finite elements,
                                      * the actual FE object used may change
                                      * from cell to cell, so we also need
                                      * different FEValues objects for
                                      * different cells; once you reinitialize
                                      * the hpFEValues object for a specific
                                      * cell, it retrieves the FEValues object
                                      * for the FE on that cell and returns it
                                      * through a function of the same name as
                                      * this one; this function here therefore
                                      * only provides the same interface so
                                      * that one can templatize on
                                      * FEValues/hpFEValues).
                                      */
    const FEValues<dim,spacedim> & get_present_fe_values () const;
    
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

                                     /**
                                      * The reinit() functions do
                                      * only that part of the work
                                      * that requires knowledge of the
                                      * type of iterator. After
                                      * setting present_cell(),
                                      * they pass on to this function,
                                      * which does the real work, and
                                      * which is independent of the
                                      * actual type of the cell
                                      * iterator.
                                      */
    void do_reinit ();
};


/**
 * Extend the interface of FEValuesBase to values that only make sense
 * when evaluating something on the surface of a cell. All the data
 * that is available in the interior of cells is also available here.
 *
 * On surfaces of mesh cells, normal vectors and boundary forms are
 * additional values that can be computed. This class provides the
 * interface to access those. Implementations are in derived classes
 * FEFaceValues and FESubfaceValues.
 *
 * The boundary form is the cross product of the images of the unit
 * tangential vectors. Therefore, it is the unit normal vector
 * multiplied with the surface element. Since it may be cheaper to
 * compute the boundary form immediately, use this value to integrate
 * <tt>n.ds</tt>.
 *
 * See FEValuesBase
 *
 *  @author Wolfgang Bangerth, 1998, Guido Kanschat, 2000, 2001
 */
template <int dim, int spacedim=dim>  
class FEFaceValuesBase : public FEValuesBase<dim,spacedim>
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
				      * @p n_faces_or_subfaces is the number
				      * of faces or subfaces that this object
				      * is to store. The actual number depends
				      * on the derived class, for
				      * FEFaceValues it is <tt>2*dim</tt>, while for
				      * the FESubfaceValues class it is
				      * <tt>2*dim*(1<<(dim-1))</tt>, i.e. the number
				      * of faces times the number of subfaces
				      * per face.
				      */
    FEFaceValuesBase (const unsigned int                 n_q_points,
		      const unsigned int                 dofs_per_cell,
		      const UpdateFlags                  update_flags,
		      const Mapping<dim,spacedim>       &mapping,
		      const FiniteElement<dim,spacedim> &fe,
		      const Quadrature<dim-1>&           quadrature);

    				     /**
				      * Return the outward normal vector to
				      * the cell at the <tt>i</tt>th quadrature
				      * point. The length of the vector
				      * is normalized to one.
				      */
    const Point<dim> & normal_vector (const unsigned int i) const;
    
    				     /**
				      * Boundary form of the
				      * transformation of the cell at
				      * the <tt>i</tt>th quadrature point.
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
				      * <tt>n.ds</tt>.
				      */
    const Tensor<1,spacedim> & boundary_form (const unsigned int i) const;
    
				     /**
				      * Return the list of outward normal
				      * vectors to the cell at the
				      * quadrature points.
				      */
    const std::vector<Point<dim> > & get_normal_vectors () const;

				     /**
				      * Return the list of outward
				      * normal vectors times the
				      * Jacobian of the surface
				      * mapping.
				      */
    const std::vector<Tensor<1,spacedim> > & get_boundary_forms () const;

				     /**
				      * Return the index of the face
				      * selected the last time the
				      * reinit() function was called.
				      */
    unsigned int get_face_index() const;

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
				      * Index of the face selected the
				      * last time the reinit()
				      * function was called.
				      */
    unsigned int present_face_index;

				     /**
				      * Store a copy of the quadrature
				      * formula here.
				      */
    const Quadrature<dim-1> quadrature;
};



/**
 * Finite element evaluated in quadrature points on a face.
 *
 * This class adds the functionality of FEFaceValuesBase to
 * FEValues; see there for more documentation.
 *
 * Since finite element functions and their derivatives may be
 * discontinuous at cell boundaries, there is no restriction of this
 * function to a mesh face. But, there are limits of these values
 * approaching the face from either of the neighboring cells.
 *
 * @author Wolfgang Bangerth, 1998, Guido Kanschat, 2000, 2001
 */
template <int dim, int spacedim=dim>  
class FEFaceValues : public FEFaceValuesBase<dim,spacedim>
{
  public:
                                     /**
                                      * Dimension in which this object
                                      * operates.
                                      */

    static const unsigned int dimension = dim;

    static const unsigned int space_dimension = spacedim;

                                     /**
                                      * Dimension of the object over
                                      * which we integrate. For the
                                      * present class, this is equal
                                      * to <code>dim-1</code>.
                                      */
    static const unsigned int integral_dimension = dim-1;

    				     /**
				      * Constructor. Gets cell
				      * independent data from mapping
				      * and finite element objects,
				      * matching the quadrature rule
				      * and update flags.
				      */
    FEFaceValues (const Mapping<dim,spacedim>       &mapping,
		  const FiniteElement<dim,spacedim> &fe,
		  const Quadrature<dim-1>           &quadrature,
		  const UpdateFlags                  update_flags);

                                     /**
				      * Constructor. Uses MappingQ1
				      * implicitly.
				      */
    FEFaceValues (const FiniteElement<dim,spacedim> &fe,
		  const Quadrature<dim-1>           &quadrature,
		  const UpdateFlags                  update_flags);

				     /**
				      * Reinitialize the gradients, Jacobi
				      * determinants, etc for the face with
				      * number @p face_no of @p cell
				      * and the given finite element.
				      */
    void reinit (const typename DoFHandler<dim,spacedim>::cell_iterator &cell,
		 const unsigned int                                      face_no);

				     /**
				      * Reinitialize the gradients,
				      * Jacobi determinants, etc for
				      * the given cell of type
				      * "iterator into a hp::DoFHandler
				      * object", and the finite
				      * element associated with this
				      * object. It is assumed that the
				      * finite element used by the
				      * given cell is also the one
				      * used by this FEValues
				      * object.
				      */
    void reinit (const typename hp::DoFHandler<dim,spacedim>::cell_iterator &cell,
		 const unsigned int                                          face_no);

				     /**
				      * Reinitialize the gradients,
				      * Jacobi determinants, etc for
				      * the given cell of type
				      * "iterator into a MGDoFHandler
				      * object", and the finite
				      * element associated with this
				      * object. It is assumed that the
				      * finite element used by the
				      * given cell is also the one
				      * used by this FEValues
				      * object.
				      */
    void reinit (const typename MGDoFHandler<dim,spacedim>::cell_iterator &cell,
		 const unsigned int                                        face_no);

				     /**
				      * Reinitialize the gradients,
				      * Jacobi determinants, etc for
				      * the given face on given cell
				      * of type "iterator into a
				      * Triangulation object", and the
				      * given finite element. Since
				      * iterators into triangulation
				      * alone only convey information
				      * about the geometry of a cell,
				      * but not about degrees of
				      * freedom possibly associated
				      * with this cell, you will not
				      * be able to call some functions
				      * of this class if they need
				      * information about degrees of
				      * freedom. These functions are,
				      * above all, the
				      * <tt>get_function_value/gradients/hessians</tt>
				      * functions. If you want to call
				      * these functions, you have to
				      * call the @p reinit variants
				      * that take iterators into
				      * DoFHandler or other DoF
				      * handler type objects.
				      */
    void reinit (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
		 const unsigned int                                         face_no);
    
                                     /**
                                      * Return a reference to this
                                      * very object.
                                      *
                                      * Though it seems that it is not very
                                      * useful, this function is there to
                                      * provide capability to the hpFEValues
                                      * class, in which case it provides the
                                      * FEValues object for the present cell
                                      * (remember that for hp finite elements,
                                      * the actual FE object used may change
                                      * from cell to cell, so we also need
                                      * different FEValues objects for
                                      * different cells; once you reinitialize
                                      * the hpFEValues object for a specific
                                      * cell, it retrieves the FEValues object
                                      * for the FE on that cell and returns it
                                      * through a function of the same name as
                                      * this one; this function here therefore
                                      * only provides the same interface so
                                      * that one can templatize on
                                      * FEValues/hpFEValues).
                                      */
    const FEFaceValues<dim,spacedim> & get_present_fe_values () const;
  private:

				     /**
				      * Do work common to the two
				      * constructors.
				      */
    void initialize (const UpdateFlags update_flags);    

                                     /**
                                      * The reinit() functions do
                                      * only that part of the work
                                      * that requires knowledge of the
                                      * type of iterator. After
                                      * setting present_cell(),
                                      * they pass on to this function,
                                      * which does the real work, and
                                      * which is independent of the
                                      * actual type of the cell
                                      * iterator.
                                      */
    void do_reinit (const unsigned int face_no);
};


/**
 * Finite element evaluated in quadrature points on a face.
 *
 * This class adds the functionality of FEFaceValuesBase to
 * FEValues; see there for more documentation.
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
template <int dim, int spacedim=dim>  
class FESubfaceValues : public FEFaceValuesBase<dim,spacedim>
{
  public:
                                     /**
                                      * Dimension in which this object
                                      * operates.
                                      */
    static const unsigned int dimension = dim;

                                     /**
                                      * Dimension of the space in
                                      * which this object operates.
                                      */
    static const unsigned int space_dimension = spacedim;

                                     /**
                                      * Dimension of the object over
                                      * which we integrate. For the
                                      * present class, this is equal
                                      * to <code>dim-1</code>.
                                      */
    static const unsigned int integral_dimension = dim-1;

				     /**
				      * Constructor. Gets cell
				      * independent data from mapping
				      * and finite element objects,
				      * matching the quadrature rule
				      * and update flags.
				      */
    FESubfaceValues (const Mapping<dim,spacedim>       &mapping,
		     const FiniteElement<dim,spacedim> &fe,
		     const Quadrature<dim-1>  &face_quadrature,
		     const UpdateFlags         update_flags);

                                     /**
				      * Constructor. Uses MappingQ1
				      * implicitly.
				      */
    FESubfaceValues (const FiniteElement<dim,spacedim> &fe,
		     const Quadrature<dim-1>  &face_quadrature,
		     const UpdateFlags         update_flags);

				     /**
				      * Reinitialize the gradients,
				      * Jacobi determinants, etc for
				      * the given cell of type
				      * "iterator into a DoFHandler
				      * object", and the finite
				      * element associated with this
				      * object. It is assumed that the
				      * finite element used by the
				      * given cell is also the one
				      * used by this
				      * FESubfaceValues object.
				      */
    void reinit (const typename DoFHandler<dim,spacedim>::cell_iterator &cell,
		 const unsigned int                    face_no,
		 const unsigned int                    subface_no);

				     /**
				      * Reinitialize the gradients,
				      * Jacobi determinants, etc for
				      * the given cell of type
				      * "iterator into a hp::DoFHandler
				      * object", and the finite
				      * element associated with this
				      * object. It is assumed that the
				      * finite element used by the
				      * given cell is also the one
				      * used by this FEValues
				      * object.
				      */
    void reinit (const typename hp::DoFHandler<dim,spacedim>::cell_iterator &cell,
		 const unsigned int                    face_no,
		 const unsigned int                    subface_no);
    
				     /**
				      * Reinitialize the gradients,
				      * Jacobi determinants, etc for
				      * the given cell of type
				      * "iterator into a MGDoFHandler
				      * object", and the finite
				      * element associated with this
				      * object. It is assumed that the
				      * finite element used by the
				      * given cell is also the one
				      * used by this FEValues
				      * object.
				      */
    void reinit (const typename MGDoFHandler<dim,spacedim>::cell_iterator &cell,
		 const unsigned int                    face_no,
		 const unsigned int                    subface_no);
    
				     /**
				      * Reinitialize the gradients,
				      * Jacobi determinants, etc for
				      * the given subface on given
				      * cell of type "iterator into a
				      * Triangulation object", and the
				      * given finite element. Since
				      * iterators into triangulation
				      * alone only convey information
				      * about the geometry of a cell,
				      * but not about degrees of
				      * freedom possibly associated
				      * with this cell, you will not
				      * be able to call some functions
				      * of this class if they need
				      * information about degrees of
				      * freedom. These functions are,
				      * above all, the
				      * <tt>get_function_value/gradients/hessians</tt>
				      * functions. If you want to call
				      * these functions, you have to
				      * call the @p reinit variants
				      * that take iterators into
				      * DoFHandler or other DoF
				      * handler type objects.
				      */
    void reinit (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
		 const unsigned int                    face_no,
		 const unsigned int                    subface_no);

                                     /**
                                      * Return a reference to this
                                      * very object.
                                      *
                                      * Though it seems that it is not very
                                      * useful, this function is there to
                                      * provide capability to the hpFEValues
                                      * class, in which case it provides the
                                      * FEValues object for the present cell
                                      * (remember that for hp finite elements,
                                      * the actual FE object used may change
                                      * from cell to cell, so we also need
                                      * different FEValues objects for
                                      * different cells; once you reinitialize
                                      * the hpFEValues object for a specific
                                      * cell, it retrieves the FEValues object
                                      * for the FE on that cell and returns it
                                      * through a function of the same name as
                                      * this one; this function here therefore
                                      * only provides the same interface so
                                      * that one can templatize on
                                      * FEValues/hpFEValues).
                                      */
    const FESubfaceValues<dim,spacedim> & get_present_fe_values () const;

                                     /**
				      * @todo Document this
				      *
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcReinitCalledWithBoundaryFace);

				     /**
				      * @todo Document this
				      *
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcFaceHasNoSubfaces);

  private:

				     /**
				      * Do work common to the two
				      * constructors.
				      */
    void initialize (const UpdateFlags update_flags);

                                     /**
                                      * The reinit() functions do
                                      * only that part of the work
                                      * that requires knowledge of the
                                      * type of iterator. After
                                      * setting present_cell(),
                                      * they pass on to this function,
                                      * which does the real work, and
                                      * which is independent of the
                                      * actual type of the cell
                                      * iterator.
                                      */
    void do_reinit (const unsigned int face_no,
                    const unsigned int subface_no);
};

/*@}*/

#ifndef DOXYGEN


/*------------------------ Inline functions: namespace FEValuesExtractors --------*/

namespace FEValuesExtractors
{
  inline
  Scalar::Scalar (const unsigned int component)
		  :
		  component (component)
  {}



  inline
  Vector::Vector (const unsigned int first_vector_component)
		  :
		  first_vector_component (first_vector_component)
  {}
}


/*------------------------ Inline functions: namespace FEValuesViews --------*/

namespace FEValuesViews
{
  template <int dim, int spacedim>  
  inline
  typename Scalar<dim,spacedim>::value_type
  Scalar<dim,spacedim>::value (const unsigned int shape_function,
			       const unsigned int q_point) const
  {
    typedef FEValuesBase<dim,spacedim> FVB;
    Assert (shape_function < fe_values.fe->dofs_per_cell,
	    ExcIndexRange (shape_function, 0, fe_values.fe->dofs_per_cell));
    Assert (fe_values.update_flags & update_values,
	    typename FVB::ExcAccessToUninitializedField());    

				     // an adaptation of the
				     // FEValuesBase::shape_value_component
				     // function except that here we know the
				     // component as fixed and we have
				     // pre-computed and cached a bunch of
				     // information. see the comments there
    if (shape_function_data[shape_function].is_nonzero_shape_function_component)
      return fe_values.shape_values(shape_function_data[shape_function]
				    .row_index,
				    q_point);
    else
      return 0;
  }
  



  template <int dim, int spacedim>  
  inline
  typename Scalar<dim,spacedim>::gradient_type
  Scalar<dim,spacedim>::gradient (const unsigned int shape_function,
				  const unsigned int q_point) const
  {
    typedef FEValuesBase<dim,spacedim> FVB;
    Assert (shape_function < fe_values.fe->dofs_per_cell,
	    ExcIndexRange (shape_function, 0, fe_values.fe->dofs_per_cell));
    Assert (fe_values.update_flags & update_gradients,
	    typename FVB::ExcAccessToUninitializedField());    

				     // an adaptation of the
				     // FEValuesBase::shape_grad_component
				     // function except that here we know the
				     // component as fixed and we have
				     // pre-computed and cached a bunch of
				     // information. see the comments there
    if (shape_function_data[shape_function].is_nonzero_shape_function_component)
      return fe_values.shape_gradients[shape_function_data[shape_function]
				       .row_index][q_point];
    else
      return gradient_type();
  }



  template <int dim, int spacedim>  
  inline
  typename Scalar<dim,spacedim>::hessian_type
  Scalar<dim,spacedim>::hessian (const unsigned int shape_function,
				 const unsigned int q_point) const
  {
    typedef FEValuesBase<dim,spacedim> FVB;
    Assert (shape_function < fe_values.fe->dofs_per_cell,
	    ExcIndexRange (shape_function, 0, fe_values.fe->dofs_per_cell));
    Assert (fe_values.update_flags & update_hessians,
	    typename FVB::ExcAccessToUninitializedField());    

				     // an adaptation of the
				     // FEValuesBase::shape_grad_component
				     // function except that here we know the
				     // component as fixed and we have
				     // pre-computed and cached a bunch of
				     // information. see the comments there
    if (shape_function_data[shape_function].is_nonzero_shape_function_component)
      return fe_values.shape_hessians[shape_function_data[shape_function].row_index][q_point];
    else
      return hessian_type();
  }



  template <int dim, int spacedim>  
  inline
  typename Vector<dim,spacedim>::value_type
  Vector<dim,spacedim>::value (const unsigned int shape_function,
			       const unsigned int q_point) const
  {
    typedef FEValuesBase<dim,spacedim> FVB;
    Assert (shape_function < fe_values.fe->dofs_per_cell,
	    ExcIndexRange (shape_function, 0, fe_values.fe->dofs_per_cell));
    Assert (fe_values.update_flags & update_values,
	    typename FVB::ExcAccessToUninitializedField());    

				     // same as for the scalar case except
				     // that we have one more index
    const int snc = shape_function_data[shape_function].single_nonzero_component;
    if (snc == -2)
      return value_type();
    else if (snc != -1)
      {
	value_type return_value;
	return_value[shape_function_data[shape_function].single_nonzero_component_index]
	  = fe_values.shape_values(snc,q_point);
	return return_value;
      }
    else
      {
	value_type return_value;
	for (unsigned int d=0; d<dim; ++d)    
	  if (shape_function_data[shape_function].is_nonzero_shape_function_component[d])
	    return_value[d]
	      = fe_values.shape_values(shape_function_data[shape_function].row_index[d],q_point);
	
	return return_value;
      }
  }

  

  template <int dim, int spacedim>  
  inline
  typename Vector<dim,spacedim>::gradient_type
  Vector<dim,spacedim>::gradient (const unsigned int shape_function,
				  const unsigned int q_point) const
  {
    typedef FEValuesBase<dim,spacedim> FVB;    
    Assert (shape_function < fe_values.fe->dofs_per_cell,
	    ExcIndexRange (shape_function, 0, fe_values.fe->dofs_per_cell));
    Assert (fe_values.update_flags & update_gradients,
	    typename FVB::ExcAccessToUninitializedField());    

				     // same as for the scalar case except
				     // that we have one more index
    const int snc = shape_function_data[shape_function].single_nonzero_component;
    if (snc == -2)
      return gradient_type();
    else if (snc != -1)
      {
	gradient_type return_value;
	return_value[shape_function_data[shape_function].single_nonzero_component_index]
	  = fe_values.shape_gradients[snc][q_point];
	return return_value;
      }
    else
      {
	gradient_type return_value;
	for (unsigned int d=0; d<dim; ++d)    
	  if (shape_function_data[shape_function].is_nonzero_shape_function_component[d])
	    return_value[d]
	      = fe_values.shape_gradients[shape_function_data[shape_function].row_index[d]][q_point];

	return return_value;
      }
  }



  template <int dim, int spacedim>  
  inline
  typename Vector<dim,spacedim>::divergence_type
  Vector<dim,spacedim>::divergence (const unsigned int shape_function,
				    const unsigned int q_point) const
  {
				     // this function works like in
				     // the case above
    typedef FEValuesBase<dim,spacedim> FVB; 
    Assert (shape_function < fe_values.fe->dofs_per_cell,
	    ExcIndexRange (shape_function, 0, fe_values.fe->dofs_per_cell));
    Assert (fe_values.update_flags & update_gradients,
	    typename FVB::ExcAccessToUninitializedField());    

				     // same as for the scalar case except
				     // that we have one more index
    const int snc = shape_function_data[shape_function].single_nonzero_component;
    if (snc == -2)
      return divergence_type();
    else if (snc != -1)
      return
	fe_values.shape_gradients[snc][q_point][shape_function_data[shape_function].single_nonzero_component_index];
    else
      {
	divergence_type return_value = 0;
	for (unsigned int d=0; d<dim; ++d)    
	  if (shape_function_data[shape_function].is_nonzero_shape_function_component[d])
	    return_value
	      += fe_values.shape_gradients[shape_function_data[shape_function].row_index[d]][q_point][d];

	return return_value;
      }
  }
  


  template <int dim, int spacedim>  
  inline
  typename Vector<dim,spacedim>::hessian_type
  Vector<dim,spacedim>::hessian (const unsigned int shape_function,
				 const unsigned int q_point) const
  {
				     // this function works like in
				     // the case above				
    typedef FEValuesBase<dim,spacedim> FVB;    
    Assert (shape_function < fe_values.fe->dofs_per_cell,
	    ExcIndexRange (shape_function, 0, fe_values.fe->dofs_per_cell));
    Assert (fe_values.update_flags & update_hessians,
	    typename FVB::ExcAccessToUninitializedField());    

				     // same as for the scalar case except
				     // that we have one more index
    const int snc = shape_function_data[shape_function].single_nonzero_component;
    if (snc == -2)
      return hessian_type();
    else if (snc != -1)
      {
	hessian_type return_value;
	return_value[shape_function_data[shape_function].single_nonzero_component_index]
	  = fe_values.shape_hessians[snc][q_point];
	return return_value;
      }
    else
      {
	hessian_type return_value;
	for (unsigned int d=0; d<dim; ++d)    
	  if (shape_function_data[shape_function].is_nonzero_shape_function_component[d])
	    return_value[d]
	      = fe_values.shape_hessians[shape_function_data[shape_function].row_index[d]][q_point];

	return return_value;
      }
  }


  namespace
  {
				     /**
				      * Return the symmetrized version of a
				      * tensor whose n'th row equals the
				      * second argument, with all other rows
				      * equal to zero.
				      */
    inline
    SymmetricTensor<2,1>
    symmetrize_single_row (const unsigned int n,
			   const Tensor<1,1> &t)
    {
      const unsigned int dim = 1;
      Assert (n < dim, ExcIndexRange (n, 0, dim));
      
      const double array[1] = { t[0] };
      return SymmetricTensor<2,1>(array);
    }


    inline
    SymmetricTensor<2,2>
    symmetrize_single_row (const unsigned int n,
			   const Tensor<1,2> &t)
    {
      switch (n)
	{
	  case 0:
	  {
	    const double array[3] = { t[0], 0, t[1]/2 };
	    return SymmetricTensor<2,2>(array);
	  }
	  case 1:
	  {
	    const double array[3] = { 0, t[1], t[0]/2 };
	    return SymmetricTensor<2,2>(array);
	  }
	  default:
	  {
	    Assert (false, ExcIndexRange (n, 0, 2));
	    return SymmetricTensor<2,2>();
	  }
	}
    }


    inline
    SymmetricTensor<2,3>
    symmetrize_single_row (const unsigned int n,
			   const Tensor<1,3> &t)
    {
      switch (n)
	{
	  case 0:
	  {
	    const double array[6] = { t[0], 0, 0, t[1]/2, t[2]/2, 0 };
	    return SymmetricTensor<2,3>(array);
	  }
	  case 1:
	  {
	    const double array[6] = { 0, t[1], 0, t[0]/2, 0, t[2]/2 };
	    return SymmetricTensor<2,3>(array);
	  }
	  case 2:
	  {
	    const double array[6] = { 0, 0, t[2], 0, t[0]/2, t[1]/2 };
	    return SymmetricTensor<2,3>(array);
	  }
	  default:
	  {
	    Assert (false, ExcIndexRange (n, 0, 3));
	    return SymmetricTensor<2,3>();
	  }
	}
    }
  }
  

  template <int dim, int spacedim>
  inline
  typename Vector<dim,spacedim>::symmetric_gradient_type
  Vector<dim,spacedim>::symmetric_gradient (const unsigned int shape_function,
					    const unsigned int q_point) const
  {
    typedef FEValuesBase<dim,spacedim> FVB;    
    Assert (shape_function < fe_values.fe->dofs_per_cell,
	    ExcIndexRange (shape_function, 0, fe_values.fe->dofs_per_cell));
    Assert (fe_values.update_flags & update_gradients,
	    typename FVB::ExcAccessToUninitializedField());    

				     // same as for the scalar case except
				     // that we have one more index
    const int snc = shape_function_data[shape_function].single_nonzero_component;
    if (snc == -2)
      return symmetric_gradient_type();
    else if (snc != -1)
      return symmetrize_single_row (shape_function_data[shape_function].single_nonzero_component_index,
				    fe_values.shape_gradients[snc][q_point]);
    else
      {
	gradient_type return_value;
	for (unsigned int d=0; d<dim; ++d)    
	  if (shape_function_data[shape_function].is_nonzero_shape_function_component[d])
	    return_value[d]
	      = fe_values.shape_gradients[shape_function_data[shape_function].row_index[d]][q_point];

	return symmetrize(return_value);
      }
  }
}



/*------------------------ Inline functions: FEValuesBase ------------------------*/



template <int dim, int spacedim>
inline
const FEValuesViews::Scalar<dim,spacedim> &
FEValuesBase<dim,spacedim>::
operator[] (const FEValuesExtractors::Scalar &scalar) const
{
  Assert (scalar.component < fe_values_views_cache.scalars.size(),
	  ExcIndexRange (scalar.component,
			 0, fe_values_views_cache.scalars.size()));
  
  return fe_values_views_cache.scalars[scalar.component];
}



template <int dim, int spacedim>  
inline
const FEValuesViews::Vector<dim,spacedim> &
FEValuesBase<dim,spacedim>::
operator[] (const FEValuesExtractors::Vector &vector) const
{
  Assert (vector.first_vector_component <
	  fe_values_views_cache.vectors.size(),
	  ExcIndexRange (vector.first_vector_component,
			 0, fe_values_views_cache.vectors.size()));
  
  return fe_values_views_cache.vectors[vector.first_vector_component];
}



template <int dim, int spacedim>  
inline
double
FEValuesBase<dim,spacedim>::shape_value (const unsigned int i,
					 const unsigned int j) const
{
  Assert (i < fe->dofs_per_cell,
	  ExcIndexRange (i, 0, fe->dofs_per_cell));
  Assert (this->update_flags & update_values,
	  ExcAccessToUninitializedField());
  Assert (fe->is_primitive (i),
	  ExcShapeFunctionNotPrimitive(i));

				   // if the entire FE is primitive,
				   // then we can take a short-cut:
  if (fe->is_primitive())
    return this->shape_values(i,j);
  else
				     // otherwise, use the mapping
				     // between shape function numbers
				     // and rows. note that by the
				     // assertions above, we know that
				     // this particular shape function
				     // is primitive, so there is no
				     // question to which vector
				     // component the call of this
				     // function refers
    return this->shape_values(this->shape_function_to_row_table[i], j);
}



template <int dim, int spacedim>  
inline
double
FEValuesBase<dim,spacedim>::shape_value_component (const unsigned int i,
						   const unsigned int j,
						   const unsigned int component) const
{
  Assert (i < fe->dofs_per_cell,
	  ExcIndexRange (i, 0, fe->dofs_per_cell));
  Assert (this->update_flags & update_values,
	  ExcAccessToUninitializedField());
  Assert (component < fe->n_components(),
	  ExcIndexRange(component, 0, fe->n_components()));
			
				   // if this particular shape
				   // function is primitive, then we
				   // can take a short-cut by checking
				   // whether the requested component
				   // is the only non-zero one (note
				   // that calling
				   // system_to_component_table only
				   // works if the shape function is
				   // primitive):
  if (fe->is_primitive(i))
    {
      if (component == fe->system_to_component_index(i).first)
	return this->shape_values(this->shape_function_to_row_table[i],j);
      else
	return 0;
    }
  else
    {
				       // no, this shape function is
				       // not primitive. then we have
				       // to loop over its components
				       // to find the corresponding
				       // row in the arrays of this
				       // object. before that check
				       // whether the shape function
				       // is non-zero at all within
				       // this component:
      if (fe->get_nonzero_components(i)[component] == false)
	return 0;

				       // count how many non-zero
				       // component the shape function
				       // has before the one we are
				       // looking for, and add this to
				       // the offset of the first
				       // non-zero component of this
				       // shape function in the arrays
				       // we index presently:
      const unsigned int
	row = (this->shape_function_to_row_table[i]
	       +
	       std::count (fe->get_nonzero_components(i).begin(),
			   fe->get_nonzero_components(i).begin()+component,
			   true));
      return this->shape_values(row, j);
    }
}



template <int dim, int spacedim>  
inline
const Tensor<1,spacedim> &
FEValuesBase<dim,spacedim>::shape_grad (const unsigned int i,
					const unsigned int j) const
{
  Assert (i < fe->dofs_per_cell,
	  ExcIndexRange (i, 0, fe->dofs_per_cell));
  Assert (this->update_flags & update_gradients,
	  ExcAccessToUninitializedField());
  Assert (fe->is_primitive (i),
	  ExcShapeFunctionNotPrimitive(i));
  Assert (i<this->shape_gradients.size(),
	  ExcIndexRange (i, 0, this->shape_gradients.size()));
  Assert (j<this->shape_gradients[0].size(),
	  ExcIndexRange (j, 0, this->shape_gradients[0].size()));

				   // if the entire FE is primitive,
				   // then we can take a short-cut:
  if (fe->is_primitive())
    return this->shape_gradients[i][j];
  else
				     // otherwise, use the mapping
				     // between shape function numbers
				     // and rows. note that by the
				     // assertions above, we know that
				     // this particular shape function
				     // is primitive, so there is no
				     // question to which vector
				     // component the call of this
				     // function refers
    return this->shape_gradients[this->shape_function_to_row_table[i]][j];
}



template <int dim, int spacedim> 
inline
Tensor<1,spacedim>
FEValuesBase<dim,spacedim>::shape_grad_component (const unsigned int i,
						  const unsigned int j,
						  const unsigned int component) const
{
  Assert (i < fe->dofs_per_cell,
	  ExcIndexRange (i, 0, fe->dofs_per_cell));
  Assert (this->update_flags & update_gradients,
	  ExcAccessToUninitializedField());
  Assert (component < fe->n_components(),
	  ExcIndexRange(component, 0, fe->n_components()));
			
				   // if this particular shape
				   // function is primitive, then we
				   // can take a short-cut by checking
				   // whether the requested component
				   // is the only non-zero one (note
				   // that calling
				   // system_to_component_table only
				   // works if the shape function is
				   // primitive):
  if (fe->is_primitive(i))
    {
      if (component == fe->system_to_component_index(i).first)
	return this->shape_gradients[this->shape_function_to_row_table[i]][j];
      else
	return Tensor<1,spacedim>();
    }
  else
    {
				       // no, this shape function is
				       // not primitive. then we have
				       // to loop over its components
				       // to find the corresponding
				       // row in the arrays of this
				       // object. before that check
				       // whether the shape function
				       // is non-zero at all within
				       // this component:
      if (fe->get_nonzero_components(i)[component] == false)
	return Tensor<1,spacedim>();

				       // count how many non-zero
				       // component the shape function
				       // has before the one we are
				       // looking for, and add this to
				       // the offset of the first
				       // non-zero component of this
				       // shape function in the arrays
				       // we index presently:
      const unsigned int
	row = (this->shape_function_to_row_table[i]
	       +
	       std::count (fe->get_nonzero_components(i).begin(),
			   fe->get_nonzero_components(i).begin()+component,
			   true));
      return this->shape_gradients[row][j];
    }
}



template <int dim, int spacedim>  
inline
const Tensor<2,spacedim> &
FEValuesBase<dim,spacedim>::shape_hessian (const unsigned int i,
					   const unsigned int j) const
{
  Assert (i < fe->dofs_per_cell,
	  ExcIndexRange (i, 0, fe->dofs_per_cell));
  Assert (this->update_flags & update_hessians,
	  ExcAccessToUninitializedField());
  Assert (fe->is_primitive (i),
	  ExcShapeFunctionNotPrimitive(i));
  Assert (i<this->shape_hessians.size(),
	  ExcIndexRange (i, 0, this->shape_hessians.size()));
  Assert (j<this->shape_hessians[0].size(),
	  ExcIndexRange (j, 0, this->shape_hessians[0].size()));

				   // if the entire FE is primitive,
				   // then we can take a short-cut:
  if (fe->is_primitive())
    return this->shape_hessians[i][j];
  else
				     // otherwise, use the mapping
				     // between shape function numbers
				     // and rows. note that by the
				     // assertions above, we know that
				     // this particular shape function
				     // is primitive, so there is no
				     // question to which vector
				     // component the call of this
				     // function refers
    return this->shape_hessians[this->shape_function_to_row_table[i]][j];
}



template <int dim, int spacedim>  
inline
const Tensor<2,spacedim> &
FEValuesBase<dim,spacedim>::shape_2nd_derivative (const unsigned int i,
						  const unsigned int j) const
{
  return shape_hessian(i,j);
}



template <int dim, int spacedim>  
inline
Tensor<2,spacedim>
FEValuesBase<dim,spacedim>::shape_hessian_component (const unsigned int i,
						     const unsigned int j,
						     const unsigned int component) const
{
  Assert (i < fe->dofs_per_cell,
	  ExcIndexRange (i, 0, fe->dofs_per_cell));
  Assert (this->update_flags & update_hessians,
	  ExcAccessToUninitializedField());
  Assert (component < fe->n_components(),
	  ExcIndexRange(component, 0, fe->n_components()));
			
				   // if this particular shape
				   // function is primitive, then we
				   // can take a short-cut by checking
				   // whether the requested component
				   // is the only non-zero one (note
				   // that calling
				   // system_to_component_table only
				   // works if the shape function is
				   // primitive):
  if (fe->is_primitive(i))
    {
      if (component == fe->system_to_component_index(i).first)
	return this->shape_hessians[this->shape_function_to_row_table[i]][j];
      else
	return Tensor<2,spacedim>();
    }
  else
    {
				       // no, this shape function is
				       // not primitive. then we have
				       // to loop over its components
				       // to find the corresponding
				       // row in the arrays of this
				       // object. before that check
				       // whether the shape function
				       // is non-zero at all within
				       // this component:
      if (fe->get_nonzero_components(i)[component] == false)
	return Tensor<2,spacedim>();

				       // count how many non-zero
				       // component the shape function
				       // has before the one we are
				       // looking for, and add this to
				       // the offset of the first
				       // non-zero component of this
				       // shape function in the arrays
				       // we index presently:
      const unsigned int
	row = (this->shape_function_to_row_table[i]
	       +
	       std::count (fe->get_nonzero_components(i).begin(),
			   fe->get_nonzero_components(i).begin()+component,
			   true));
      return this->shape_hessians[row][j];
    }
}



template <int dim, int spacedim>  
inline
Tensor<2,spacedim>
FEValuesBase<dim,spacedim>::shape_2nd_derivative_component (const unsigned int i,
							    const unsigned int j,
							    const unsigned int component) const
{
  return shape_hessian_component(i,j,component);
}



template <int dim, int spacedim>  
inline
const FiniteElement<dim,spacedim> & 
FEValuesBase<dim,spacedim>::get_fe () const
{
  return *fe;
}


template <int dim, int spacedim>  
inline
const Mapping<dim,spacedim> & 
FEValuesBase<dim,spacedim>::get_mapping () const
{
  return *mapping;
}



template <int dim, int spacedim>  
inline
UpdateFlags
FEValuesBase<dim,spacedim>::get_update_flags () const
{
  return this->update_flags;
}



template <int dim, int spacedim>  
const typename Triangulation<dim,spacedim>::cell_iterator
FEValuesBase<dim,spacedim>::get_cell () const
{
  return *present_cell;
}


template <int dim, int spacedim>  
inline
const std::vector<Point<spacedim> > &
FEValuesBase<dim,spacedim>::get_quadrature_points () const
{
  Assert (this->update_flags & update_quadrature_points, ExcAccessToUninitializedField());
  return this->quadrature_points;
}



template <int dim, int spacedim>  
inline
const std::vector<double> &
FEValuesBase<dim,spacedim>::get_JxW_values () const
{
  Assert (this->update_flags & update_JxW_values, ExcAccessToUninitializedField());
  return this->JxW_values;
}



template <int dim, int spacedim>  
inline
const std::vector<Tensor<2,spacedim> >&
FEValuesBase<dim,spacedim>::get_jacobians () const
{
  Assert (this->update_flags & update_jacobians, ExcAccessToUninitializedField());
  return this->jacobians;
}



template <int dim, int spacedim>
inline
const std::vector<Tensor<3,spacedim> >&
FEValuesBase<dim,spacedim>::get_jacobian_grads () const
{
  Assert (this->update_flags & update_jacobian_grads, ExcAccessToUninitializedField());
  return this->jacobian_grads;
}



template <int dim, int spacedim>
inline
const std::vector<Tensor<2,spacedim> >&
FEValuesBase<dim,spacedim>::get_inverse_jacobians () const
{
  Assert (this->update_flags & update_inverse_jacobians, ExcAccessToUninitializedField());
  return this->inverse_jacobians;
}



template <int dim, int spacedim>
const Point<spacedim> &
FEValuesBase<dim,spacedim>::quadrature_point (const unsigned int i) const
{
  Assert (this->update_flags & update_quadrature_points, ExcAccessToUninitializedField());
  Assert (i<this->quadrature_points.size(), ExcIndexRange(i, 0, this->quadrature_points.size()));
  
  return this->quadrature_points[i];
}




template <int dim, int spacedim>  
inline
double
FEValuesBase<dim,spacedim>::JxW (const unsigned int i) const
{
  Assert (this->update_flags & update_JxW_values, ExcAccessToUninitializedField());
  Assert (i<this->JxW_values.size(), ExcIndexRange(i, 0, this->JxW_values.size()));
  
  return this->JxW_values[i];
}



template <int dim, int spacedim>  
inline
const Tensor<2,spacedim> &
FEValuesBase<dim,spacedim>::jacobian (const unsigned int i) const
{
  Assert (this->update_flags & update_jacobians, ExcAccessToUninitializedField());
  Assert (i<this->jacobians.size(), ExcIndexRange(i, 0, this->jacobians.size()));
  
  return this->jacobians[i];
}



template <int dim, int spacedim>
inline
const Tensor<3,spacedim> &
FEValuesBase<dim,spacedim>::jacobian_grad (const unsigned int i) const
{
  Assert (this->update_flags & update_jacobian_grads, ExcAccessToUninitializedField());
  Assert (i<this->jacobian_grads.size(), ExcIndexRange(i, 0, this->jacobian_grads.size()));
  
  return this->jacobian_grads[i];
}



template <int dim, int spacedim>
inline
const Tensor<2,spacedim> &
FEValuesBase<dim,spacedim>::inverse_jacobian (const unsigned int i) const
{
  Assert (this->update_flags & update_inverse_jacobians, ExcAccessToUninitializedField());
  Assert (i<this->inverse_jacobians.size(), ExcIndexRange(i, 0, this->inverse_jacobians.size()));
  
  return this->inverse_jacobians[i];
}



template <int dim, int spacedim>
template <class InputVector>
inline
void
FEValuesBase<dim,spacedim>::get_function_grads (const InputVector           &fe_function,
						std::vector<Tensor<1,spacedim> > &gradients) const
{
  get_function_gradients(fe_function, gradients);
}



template <int dim, int spacedim>  
template <class InputVector>
inline
void
FEValuesBase<dim,spacedim>::get_function_grads (
  const InputVector& fe_function,
  const VectorSlice<const std::vector<unsigned int> >& indices,
  std::vector<Tensor<1,spacedim> > &values) const
{
  get_function_gradients(fe_function, indices, values);
}



template <int dim, int spacedim>  
template <class InputVector>
inline
void
FEValuesBase<dim,spacedim>::
get_function_grads (const InputVector                         &fe_function,
		    std::vector<std::vector<Tensor<1,spacedim> > > &gradients) const
{
  get_function_gradients(fe_function, gradients);
}



template <int dim, int spacedim>  
template <class InputVector>
inline
void
FEValuesBase<dim,spacedim>::get_function_grads (
  const InputVector& fe_function,
  const VectorSlice<const std::vector<unsigned int> >& indices,
  std::vector<std::vector<Tensor<1,spacedim> > >& values,
  bool q_points_fastest) const
{
  get_function_gradients(fe_function, indices, values, q_points_fastest);
}



template <int dim, int spacedim>  
template <class InputVector>
inline
void
FEValuesBase<dim,spacedim>::
get_function_2nd_derivatives (const InputVector           &fe_function,
			      std::vector<Tensor<2,spacedim> > &hessians) const
{
  get_function_hessians(fe_function, hessians);
}



template <int dim, int spacedim>  
template <class InputVector>
inline
void
FEValuesBase<dim,spacedim>::
get_function_2nd_derivatives (const InputVector                         &fe_function,
			      std::vector<std::vector<Tensor<2,spacedim> > > &hessians,
			      bool quadrature_points_fastest) const
{
  get_function_hessians(fe_function, hessians, quadrature_points_fastest);
}



template <int dim, int spacedim>  
inline
const Point<spacedim> &
FEValuesBase<dim,spacedim>::cell_normal_vector (const unsigned int i) const
{
  typedef FEValuesBase<dim,spacedim> FVB;
  Assert (i<this->cell_normal_vectors.size(),
	  ExcIndexRange(i, 0, this->cell_normal_vectors.size()));
  Assert (this->update_flags & update_cell_normal_vectors,
	  typename FVB::ExcAccessToUninitializedField());
  
  return this->cell_normal_vectors[i];
}




/*------------------------ Inline functions: FEValues ----------------------------*/


template <int dim, int spacedim>  
inline
const Quadrature<dim> &
FEValues<dim,spacedim>::get_quadrature () const 
{
  return quadrature;
}



template <int dim, int spacedim>  
inline
const FEValues<dim,spacedim> &
FEValues<dim,spacedim>::get_present_fe_values () const 
{
  return *this;
}


/*------------------------ Inline functions: FEFaceValuesBase --------------------*/


template <int dim, int spacedim>  
inline
const Point<dim> &
FEFaceValuesBase<dim,spacedim>::normal_vector (const unsigned int i) const
{
  typedef FEValuesBase<dim,spacedim> FVB;
  Assert (i<this->normal_vectors.size(),
	  ExcIndexRange(i, 0, this->normal_vectors.size()));
  Assert (this->update_flags & update_normal_vectors,
	  typename FVB::ExcAccessToUninitializedField());
  
  return this->normal_vectors[i];
}


template <int dim, int spacedim>  
inline
unsigned int
FEFaceValuesBase<dim,spacedim>::get_face_index () const
{
  return present_face_index;
}


/*------------------------ Inline functions: FE*FaceValues --------------------*/

template <int dim, int spacedim>  
inline
const Quadrature<dim-1> &
FEFaceValuesBase<dim,spacedim>::get_quadrature () const 
{
  return quadrature;
}



template <int dim, int spacedim>  
inline
const FEFaceValues<dim,spacedim> &
FEFaceValues<dim,spacedim>::get_present_fe_values () const 
{
  return *this;
}



template <int dim, int spacedim>  
inline
const FESubfaceValues<dim,spacedim> &
FESubfaceValues<dim,spacedim>::get_present_fe_values () const 
{
  return *this;
}



template <int dim, int spacedim>  
inline
const Tensor<1,spacedim> &
FEFaceValuesBase<dim,spacedim>::boundary_form (const unsigned int i) const
{
  typedef FEValuesBase<dim,spacedim> FVB;
  Assert (i<this->boundary_forms.size(),
	  ExcIndexRange(i, 0, this->boundary_forms.size()));
  Assert (this->update_flags & update_boundary_forms,
	  typename FVB::ExcAccessToUninitializedField());
  
  return this->boundary_forms[i];
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
