// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2014 by the deal.II authors
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


#ifndef __deal2__matrix_free_fe_evaluation_h
#define __deal2__matrix_free_fe_evaluation_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/shape_info.h>
#include <deal.II/matrix_free/mapping_fe_evaluation.h>


DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  namespace distributed
  {
    template <typename> class Vector;
  }
}



// forward declarations
namespace internal
{
  DeclException0 (ExcAccessToUninitializedField);
}



/**
 * This is the base class for the FEEvaluation classes. This class is a base
 * class and needs usually not be called in user code. It does not have any
 * public constructor. The usage is through the class FEEvaluation instead. It
 * implements a reinit method that is used to set pointers so that operations
 * on quadrature points can be performed quickly, access functions to vectors
 * for the @p read_dof_values, @p set_dof_values, and @p
 * distributed_local_to_global functions, as well as methods to access values
 * and gradients of finite element functions.
 *
 * This class has three template arguments:
 *
 * @param dim Dimension in which this class is to be used
 *
 * @param n_components Number of vector components when solving a system of
 *                  PDEs. If the same operation is applied to several
 *                  components of a PDE (e.g. a vector Laplace equation), they
 *                  can be applied simultaneously with one call (and often
 *                  more efficiently)
 *
 * @param Number Number format, usually @p double or @p float
 *
 * @author Katharina Kormann and Martin Kronbichler, 2010, 2011
 */
template <int dim, int n_components_, typename Number>
class FEEvaluationBase
{
public:
  typedef Number                            number_type;
  typedef Tensor<1,n_components_,VectorizedArray<Number> > value_type;
  typedef Tensor<1,n_components_,Tensor<1,dim,VectorizedArray<Number> > > gradient_type;
  static const unsigned int dimension     = dim;
  static const unsigned int n_components  = n_components_;

  /**
   * @name 1: General operations
   */
  //@{
  /**
   * Initializes the operation pointer to the current cell. Unlike the
   * FEValues::reinit function, where the information related to a particular
   * cell is generated in the reinit call, this function is very cheap since
   * all data is pre-computed in @p matrix_free, and only a few indices have
   * to be set appropriately.
   */
  void reinit (const unsigned int cell);

  /**
   * For the transformation information stored in MappingInfo, this function
   * returns the index which belongs to the current cell as specified in @p
   * reinit. Note that MappingInfo has different fields for Cartesian cells,
   * cells with affine mapping and with general mappings, so in order to
   * access the correct data, this interface must be used together with
   * get_cell_type.
   */
  unsigned int get_cell_data_number() const;

  /**
   * Returns the type of the cell the @p reinit function has been called
   * for. Valid values are @p cartesian for Cartesian cells (which allows for
   * considerable data compression), @p affine for cells with affine mappings,
   * and @p general for general cells without any compressed storage applied.
   */
  internal::MatrixFreeFunctions::CellType get_cell_type() const;

  /**
   * Returns a reference to the ShapeInfo object currently in use.
   */
  const internal::MatrixFreeFunctions::ShapeInfo<Number> &
  get_shape_info() const;

  /**
   * Fills the JxW values currently used.
   */
  void
  fill_JxW_values(AlignedVector<VectorizedArray<Number> > &JxW_values) const;

  //@}

  /**
   * @name 2: Reading from and writing to vectors
   */
  //@{
  /**
   * For the vector @p src, read out the values on the degrees of freedom of
   * the current cell, and store them internally. Similar functionality as the
   * function DoFAccessor::get_interpolated_dof_values when no constraints are
   * present, but it also includes constraints from hanging nodes, so one can
   * see it as a similar function to ConstraintMatrix::read_dof_values as
   * well. Note that if vectorization is enabled, the DoF values for several
   * cells are set.
   *
   * If some constraints on the vector are inhomogeneous, use the function
   * read_dof_values_plain instead and provide the vector with useful data
   * also in constrained positions by calling
   * ConstraintMatrix::distribute. When accessing vector entries during the
   * solution of linear systems, the temporary solution should always have
   * homogeneous constraints and this method is the correct one.
   *
   * If the class was constructed through a MappingFEEvaluation object, only
   * one single cell is used by this class and this function extracts the
   * values of the underlying components on this cell. This call is slower
   * than the ones done through a MatrixFree object and lead to a structure
   * that does not effectively use vectorization in the evaluate routines
   * based on these values (instead, VectorizedArray<Number>::n_array_elements
   * same copies are worked on).
   */
  template <typename VectorType>
  void read_dof_values (const VectorType &src);

  /**
   * For a collection of several vector @p src, read out the values on the
   * degrees of freedom of the current cell for @p n_components (template
   * argument), starting at @p first_index, and store them internally. Similar
   * functionality as the function ConstraintMatrix::read_dof_values.  Note
   * that if vectorization is enabled, the DoF values for several cells are
   * set.
   */
  template <typename VectorType>
  void read_dof_values (const std::vector<VectorType> &src,
                        const unsigned int             first_index=0);

  /**
   * Reads data from several vectors. Same as other function with std::vector,
   * but accepts a vector of pointers to vectors.
   */
  template <typename VectorType>
  void read_dof_values (const std::vector<VectorType *> &src,
                        const unsigned int              first_index=0);

  /**
   * For the vector @p src, read out the values on the degrees of freedom of
   * the current cell, and store them internally. Similar functionality as the
   * function DoFAccessor::get_interpolated_dof_values. As opposed to the
   * read_dof_values function, this function reads out the plain entries from
   * vectors, without taking stored constraints into account. This way of
   * access is appropriate when the constraints have been distributed on the
   * vector by a call to ConstraintMatrix::distribute previously. This
   * function is also necessary when inhomogeneous constraints are to be used,
   * as MatrixFree can only handle homogeneous constraints. Note that if
   * vectorization is enabled, the DoF values for several cells are set.
   *
   * If the class was constructed through a MappingFEEvaluation object, only
   * one single cell is used by this class and this function extracts the
   * values of the underlying components on this cell. This call is slower
   * than the ones done through a MatrixFree object and lead to a structure
   * that does not effectively use vectorization in the evaluate routines
   * based on these values (instead, VectorizedArray<Number>::n_array_elements
   * same copies are worked on). In that case, no constraints can be
   * processed as these are not available here.
   */
  template <typename VectorType>
  void read_dof_values_plain (const VectorType &src);

  /**
   * For a collection of several vector @p src, read out the values on the
   * degrees of freedom of the current cell for @p n_components (template
   * argument), starting at @p first_index, and store them internally. Similar
   * functionality as the function DoFAccessor::read_dof_values.  Note that if
   * vectorization is enabled, the DoF values for several cells are set.
   */
  template <typename VectorType>
  void read_dof_values_plain (const std::vector<VectorType> &src,
                              const unsigned int             first_index=0);

  /**
   * Reads data from several vectors without resolving constraints. Same as
   * other function with std::vector, but accepts a vector of pointers to
   * vectors.
   */
  template <typename VectorType>
  void read_dof_values_plain (const std::vector<VectorType *> &src,
                              const unsigned int              first_index=0);

  /**
   * Takes the values stored internally on dof values of the current cell and
   * sums them into the vector @p dst. The function also applies constraints
   * during the write operation. The functionality is hence similar to the
   * function ConstraintMatrix::distribute_local_to_global. If vectorization
   * is enabled, the DoF values for several cells are used.
   *
   * If the class was constructed through a MappingFEEvaluation object, only
   * one single cell is used by this class and this function extracts the
   * values of the underlying components on this cell. This call is slower
   * than the ones done through a MatrixFree object and lead to a structure
   * that does not effectively use vectorization in the evaluate routines
   * based on these values (instead, VectorizedArray<Number>::n_array_elements
   * same copies are worked on). In that case, no constraints can be
   * processed as these are not available here.
   */
  template<typename VectorType>
  void distribute_local_to_global (VectorType &dst) const;

  /**
   * Takes the values stored internally on dof values of the current cell for
   * a vector-valued problem consisting of @p n_components (template argument)
   * and sums them into the collection of vectors vector @p dst, starting at
   * index @p first_index. The function also applies constraints during the
   * write operation. The functionality is hence similar to the function
   * ConstraintMatrix::distribute_local_to_global. If vectorization is
   * enabled, the DoF values for several cells are used.
   */
  template<typename VectorType>
  void distribute_local_to_global (std::vector<VectorType> &dst,
                                   const unsigned int       first_index=0) const;

  /**
   * Writes data to several vectors. Same as other function with std::vector,
   * but accepts a vector of pointers to vectors.
   */
  template<typename VectorType>
  void distribute_local_to_global (std::vector<VectorType *> &dst,
                                   const unsigned int       first_index=0) const;

  /**
   * Takes the values stored internally on dof values of the current cell and
   * sums them into the vector @p dst. The function also applies constraints
   * during the write operation. The functionality is hence similar to the
   * function ConstraintMatrix::distribute_local_to_global.  Note that if
   * vectorization is enabled, the DoF values for several cells are used.
   *
   * If the class was constructed through a MappingFEEvaluation object, only
   * one single cell is used by this class and this function extracts the
   * values of the underlying components on this cell. This call is slower
   * than the ones done through a MatrixFree object and lead to a structure
   * that does not effectively use vectorization in the evaluate routines
   * based on these values (instead, VectorizedArray<Number>::n_array_elements
   * same copies are worked on). In that case, no constraints can be
   * processed as these are not available here.
   */
  template<typename VectorType>
  void set_dof_values (VectorType &dst) const;

  /**
   * Takes the values stored internally on dof values of the current cell for
   * a vector-valued problem consisting of @p n_components (template argument)
   * and sums them into the collection of vectors vector @p dst, starting at
   * index @p first_index. The function also applies constraints during the
   * write operation. The functionality is hence similar to the function
   * ConstraintMatrix::distribute_local_to_global.  Note that if vectorization
   * is enabled, the DoF values for several cells are used.
   */
  template<typename VectorType>
  void set_dof_values (std::vector<VectorType> &dst,
                       const unsigned int       first_index=0) const;

  /**
   * Writes data to several vectors. Same as other function with std::vector,
   * but accepts a vector of pointers to vectors.
   */
  template<typename VectorType>
  void set_dof_values (std::vector<VectorType *> &dst,
                       const unsigned int        first_index=0) const;

  //@}

  /**
   * @name 3: Data access
   */
  //@{
  /**
   * Returns the value stored for the local degree of freedom with index @p
   * dof. If the object is vector-valued, a vector-valued return argument is
   * given. Note that when vectorization is enabled, values from several cells
   * are grouped together. If @p set_dof_values was called last, the value
   * corresponds to the one set there. If @p integrate was called last, it
   * instead corresponds to the value of the integrated function with the test
   * function of the given index.
   *
   * Note that the derived class FEEvaluationAccess overloads this operation
   * with specializations for the scalar case (n_components == 1) and for the
   * vector-valued case (n_components == dim).
   */
  value_type get_dof_value (const unsigned int dof) const;

  /**
   * Write a value to the field containing the degrees of freedom with
   * component @p dof. Writes to the same field as is accessed through @p
   * get_dof_value. Therefore, the original data that was read from a vector
   * is overwritten as soon as a value is submitted.
   *
   * Note that the derived class FEEvaluationAccess overloads this operation
   * with specializations for the scalar case (n_components == 1) and for the
   * vector-valued case (n_components == dim).
   */
  void submit_dof_value (const value_type   val_in,
                         const unsigned int dof);

  /**
   * Returns the value of a finite element function at quadrature point number
   * @p q_point after a call to @p evaluate(true,...), or the value that has
   * been stored there with a call to @p submit_value. If the object is
   * vector-valued, a vector-valued return argument is given. Note that when
   * vectorization is enabled, values from several cells are grouped together.
   *
   * Note that the derived class FEEvaluationAccess overloads this operation
   * with specializations for the scalar case (n_components == 1) and for the
   * vector-valued case (n_components == dim).
   */
  value_type get_value (const unsigned int q_point) const;

  /**
   * Write a value to the field containing the values on quadrature points
   * with component @p q_point. Access to the same field as through @p
   * get_value. If applied before the function @p integrate(true,...) is
   * called, this specifies the value which is tested by all basis function on
   * the current cell and integrated over.
   *
   * Note that the derived class FEEvaluationAccess overloads this operation
   * with specializations for the scalar case (n_components == 1) and for the
   * vector-valued case (n_components == dim).
   */
  void submit_value (const value_type   val_in,
                     const unsigned int q_point);

  /**
   * Returns the gradient of a finite element function at quadrature point
   * number @p q_point after a call to @p evaluate(...,true,...), or the value
   * that has been stored there with a call to @p submit_gradient.
   *
   * Note that the derived class FEEvaluationAccess overloads this operation
   * with specializations for the scalar case (n_components == 1) and for the
   * vector-valued case (n_components == dim).
   */
  gradient_type get_gradient (const unsigned int q_point) const;

  /**
   * Write a contribution that is tested by the gradient to the field
   * containing the values on quadrature points with component @p
   * q_point. Access to the same field as through @p get_gradient. If applied
   * before the function @p integrate(...,true) is called, this specifies what
   * is tested by all basis function gradients on the current cell and
   * integrated over.
   *
   * Note that the derived class FEEvaluationAccess overloads this operation
   * with specializations for the scalar case (n_components == 1) and for the
   * vector-valued case (n_components == dim).
   */
  void submit_gradient(const gradient_type grad_in,
                       const unsigned int  q_point);

  /**
   * Returns the Hessian of a finite element function at quadrature point
   * number @p q_point after a call to @p evaluate(...,true). If only the
   * diagonal or even the trace of the Hessian, the Laplacian, is needed, use
   * the other functions below.
   *
   * Note that the derived class FEEvaluationAccess overloads this operation
   * with specializations for the scalar case (n_components == 1) and for the
   * vector-valued case (n_components == dim).
   */
  Tensor<1,n_components_,Tensor<2,dim,VectorizedArray<Number> > >
  get_hessian (const unsigned int q_point) const;

  /**
   * Returns the diagonal of the Hessian of a finite element function at
   * quadrature point number @p q_point after a call to @p evaluate(...,true).
   *
   * Note that the derived class FEEvaluationAccess overloads this operation
   * with specializations for the scalar case (n_components == 1) and for the
   * vector-valued case (n_components == dim).
   */
  gradient_type get_hessian_diagonal (const unsigned int q_point) const;

  /**
   * Returns the Laplacian (i.e., the trace of the Hessian) of a finite
   * element function at quadrature point number @p q_point after a call to @p
   * evaluate(...,true). Compared to the case when computing the full Hessian,
   * some operations can be saved when only the Laplacian is requested.
   *
   * Note that the derived class FEEvaluationAccess overloads this operation
   * with specializations for the scalar case (n_components == 1) and for the
   * vector-valued case (n_components == dim).
   */
  value_type get_laplacian (const unsigned int q_point) const;

  /**
   * Takes values on quadrature points, multiplies by the Jacobian determinant
   * and quadrature weights (JxW) and sums the values for all quadrature
   * points on the cell. The result is a scalar, representing the integral
   * over the function over the cell. If a vector-element is used, the
   * resulting components are still separated. Moreover, if vectorization is
   * enabled, the integral values of several cells are represented together.
   */
  value_type integrate_value () const;

  //@}

  /**
   * @name 4: Access to internal data
   */
  //@{
  /**
   * Returns a read-only pointer to the first field of the dof values. This is
   * the data field the read_dof_values() functions write into. First come the
   * the dof values for the first component, then all values for the second
   * component, and so on. This is related to the internal data structures
   * used in this class. In general, it is safer to use the get_dof_value()
   * function instead.
   */
  const VectorizedArray<Number> *begin_dof_values () const;

  /**
   * Returns a read and write pointer to the first field of the dof
   * values. This is the data field the read_dof_values() functions write
   * into. First come the the dof values for the first component, then all
   * values for the second component, and so on. This is related to the
   * internal data structures used in this class. In general, it is safer to
   * use the get_dof_value() function instead.
   */
  VectorizedArray<Number> *begin_dof_values ();

  /**
   * Returns a read-only pointer to the first field of function values on
   * quadrature points. First come the function values on all quadrature
   * points for the first component, then all values for the second component,
   * and so on. This is related to the internal data structures used in this
   * class. The raw data after a call to @p evaluate only contains unit cell
   * operations, so possible transformations, quadrature weights etc. must be
   * applied manually. In general, it is safer to use the get_value() function
   * instead, which does all the transformation internally.
   */
  const VectorizedArray<Number> *begin_values () const;

  /**
   * Returns a read and write pointer to the first field of function values on
   * quadrature points. First come the function values on all quadrature
   * points for the first component, then all values for the second component,
   * and so on. This is related to the internal data structures used in this
   * class. The raw data after a call to @p evaluate only contains unit cell
   * operations, so possible transformations, quadrature weights etc. must be
   * applied manually. In general, it is safer to use the get_value() function
   * instead, which does all the transformation internally.
   */
  VectorizedArray<Number> *begin_values ();

  /**
   * Returns a read-only pointer to the first field of function gradients on
   * quadrature points. First comes the x-component of the gradient for the
   * first component on all quadrature points, then the y-component, and so
   * on. Next comes the x-component of the second component, and so on. This
   * is related to the internal data structures used in this class. The raw
   * data after a call to @p evaluate only contains unit cell operations, so
   * possible transformations, quadrature weights etc. must be applied
   * manually. In general, it is safer to use the get_gradient() function
   * instead, which does all the transformation internally.
   */
  const VectorizedArray<Number> *begin_gradients () const;

  /**
   * Returns a read and write pointer to the first field of function gradients
   * on quadrature points. First comes the x-component of the gradient for the
   * first component on all quadrature points, then the y-component, and so
   * on. Next comes the x-component of the second component, and so on. This
   * is related to the internal data structures used in this class. The raw
   * data after a call to @p evaluate only contains unit cell operations, so
   * possible transformations, quadrature weights etc. must be applied
   * manually. In general, it is safer to use the get_gradient() function
   * instead, which does all the transformation internally.
   */
  VectorizedArray<Number> *begin_gradients ();

  /**
   * Returns a read-only pointer to the first field of function hessians on
   * quadrature points. First comes the xx-component of the hessian for the
   * first component on all quadrature points, then the yy-component,
   * zz-component in (3D), then the xy-component, and so on. Next comes the
   * xx-component of the second component, and so on. This is related to the
   * internal data structures used in this class. The raw data after a call to
   * @p evaluate only contains unit cell operations, so possible
   * transformations, quadrature weights etc. must be applied manually. In
   * general, it is safer to use the get_laplacian() or get_hessian()
   * functions instead, which does all the transformation internally.
   */
  const VectorizedArray<Number> *begin_hessians () const;

  /**
   * Returns a read and write pointer to the first field of function hessians
   * on quadrature points. First comes the xx-component of the hessian for the
   * first component on all quadrature points, then the yy-component,
   * zz-component in (3D), then the xy-component, and so on. Next comes the
   * xx-component of the second component, and so on. This is related to the
   * internal data structures used in this class. The raw data after a call to
   * @p evaluate only contains unit cell operations, so possible
   * transformations, quadrature weights etc. must be applied manually. In
   * general, it is safer to use the get_laplacian() or get_hessian()
   * functions instead, which does all the transformation internally.
   */
  VectorizedArray<Number> *begin_hessians ();

  //@}

protected:

  /**
   * Constructor. Made protected to prevent users from directly using this
   * class. Takes all data stored in MatrixFree. If applied to problems with
   * more than one finite element or more than one quadrature formula selected
   * during construction of @p matrix_free, @p fe_no and @p quad_no allow to
   * select the appropriate components.
   */
  FEEvaluationBase (const MatrixFree<dim,Number> &matrix_free,
                    const unsigned int            fe_no,
                    const unsigned int            quad_no,
                    const unsigned int            fe_degree,
                    const unsigned int            n_q_points);

  /**
   * Constructor that comes with reduced functionality and works similar as
   * FEValues. The user has to provide a structure of type MappingFEEvaluation
   * and a DoFHandler in order to allow for reading out the finite element
   * data. It uses the data provided by dof_handler.get_fe(). If the element
   * is vector-valued, the optional argument allows to specify the index of
   * the base element (as long as the element is primitive, non-primitive are
   * not supported currently).
   *
   * With initialization from a FEValues object, no call to a reinit method of
   * this class is necessary. Instead, it is enough if the geometry is
   * initialized to a given cell iterator. It can also read from or write to
   * vectors in the standard way for DoFHandler<dim>::active_cell_iterator
   * types (which is less efficient with MPI since index translation has to be
   * done), but of course only for one cell at a time. Hence, a kernel using
   * this method does not vectorize over several elements (which is most
   * efficient for vector operations), but only possibly within the element if
   * the evaluate/integrate routines are combined (e.g. for computing cell
   * matrices).
   */
  FEEvaluationBase (const MappingFEEvaluation<dim,Number> &geometry,
                    const DoFHandler<dim>                 &dof_handler,
                    const unsigned int                     first_selected_component = 0);

  /**
   * Copy constructor
   */
  FEEvaluationBase (const FEEvaluationBase &other);

  /**
   * A unified function to read from and write into vectors based on the given
   * template operation. It can perform the operation for @p read_dof_values,
   * @p distribute_local_to_global, and @p set_dof_values. It performs the
   * operation for several vectors at a time.
   */
  template<typename VectorType, typename VectorOperation>
  void read_write_operation (const VectorOperation &operation,
                             VectorType            *vectors[]) const;

  /**
   * For a collection of several vector @p src, read out the values on the
   * degrees of freedom of the current cell for @p n_components (template
   * argument), and store them internally. Similar functionality as the
   * function DoFAccessor::read_dof_values. Note that if vectorization is
   * enabled, the DoF values for several cells are set.
   */
  template<typename VectorType>
  void read_dof_values_plain (const VectorType *src_data[]);

  /**
   * Internal data fields that store the values. Derived classes will know the
   * length of all arrays at compile time and allocate the memory on the
   * stack. This makes it possible to cheaply set up a FEEvaluation object and
   * write thread-safe programs by letting each thread own a private object of
   * this type. In this base class, only pointers to the actual data are
   * stored.
   *
   * This field stores the values for local degrees of freedom (e.g. after
   * reading out from a vector but before applying unit cell transformations
   * or before distributing them into a result vector). The methods
   * get_dof_value() and submit_dof_value() read from or write to this field.
   */
  VectorizedArray<Number> *values_dofs[n_components];

  /**
   * This field stores the values of the finite element function on quadrature
   * points after applying unit cell transformations or before
   * integrating. The methods get_value() and submit_value() access this
   * field.
   */
  VectorizedArray<Number> *values_quad[n_components];

  /**
   * This field stores the gradients of the finite element function on
   * quadrature points after applying unit cell transformations or before
   * integrating. The methods get_gradient() and submit_gradient() (as well as
   * some specializations like get_symmetric_gradient() or get_divergence())
   * access this field.
   */
  VectorizedArray<Number> *gradients_quad[n_components][dim];

  /**
   * This field stores the Hessians of the finite element function on
   * quadrature points after applying unit cell transformations. The methods
   * get_hessian(), get_laplacian(), get_hessian_diagonal() access this field.
   */
  VectorizedArray<Number> *hessians_quad[n_components][(dim*(dim+1))/2];

  /**
   * Stores the number of the quadrature formula of the present cell.
   */
  const unsigned int quad_no;

  /**
   * Stores the number of components in the finite element as detected in the
   * MatrixFree storage class for comparison with the template argument.
   */
  const unsigned int n_fe_components;

  /**
   * Stores the active fe index for this class for efficient indexing in the
   * hp case.
   */
  const unsigned int active_fe_index;

  /**
   * Stores the active quadrature index for this class for efficient indexing
   * in the hp case.
   */
  const unsigned int active_quad_index;

  /**
   * Stores a pointer to the underlying data.
   */
  const MatrixFree<dim,Number> *matrix_info;

  /**
   * Stores a pointer to the underlying DoF indices and constraint
   * description for the component specified at construction. Also contained
   * in matrix_info, but it simplifies code if we store a reference to it.
   */
  const internal::MatrixFreeFunctions::DoFInfo *dof_info;

  /**
   * Stores a pointer to the underlying transformation data from unit to
   * real cells for the given quadrature formula specified at construction.
   * Also contained in matrix_info, but it simplifies code if we store a
   * reference to it.
   */
  const internal::MatrixFreeFunctions::MappingInfo<dim,Number> *mapping_info;

  /**
   * In case the class is initialized from MappingFEEvaluation instead of
   * MatrixFree, this data structure holds the evaluated shape data.
   */
  std_cxx11::shared_ptr<internal::MatrixFreeFunctions::ShapeInfo<Number> > stored_shape_info;

  /**
   * Stores a pointer to the unit cell shape data, i.e., values, gradients and
   * Hessians in 1D at the quadrature points that constitute the tensor
   * product. Also contained in matrix_info, but it simplifies code if we
   * store a reference to it.
   */
  const internal::MatrixFreeFunctions::ShapeInfo<Number> *data;

  /**
   * A pointer to the Cartesian Jacobian information of the present cell. Only
   * set to a useful value if on a Cartesian cell, otherwise zero.
   */
  const Tensor<1,dim,VectorizedArray<Number> > *cartesian_data;

  /**
   * A pointer to the Jacobian information of the present cell. Only set to a
   * useful value if on a non-Cartesian cell.
   */
  const Tensor<2,dim,VectorizedArray<Number> > *jacobian;

  /**
   * A pointer to the Jacobian determinant of the present cell. If on a
   * Cartesian cell or on a cell with constant Jacobian, this is just the
   * Jacobian determinant, otherwise the Jacobian determinant times the
   * quadrature weight.
   */
  const VectorizedArray<Number> *J_value;

  /**
   * A pointer to the quadrature weights of the underlying quadrature formula.
   */
  const VectorizedArray<Number> *quadrature_weights;

  /**
   * A pointer to the quadrature points on the present cell.
   */
  const Point<dim,VectorizedArray<Number> > *quadrature_points;

  /**
   * A pointer to the diagonal part of the Jacobian gradient on the present
   * cell. Only set to a useful value if on a general cell with non-constant
   * Jacobian.
   */
  const Tensor<2,dim,VectorizedArray<Number> > *jacobian_grad;

  /**
   * A pointer to the upper diagonal part of the Jacobian gradient on the
   * present cell. Only set to a useful value if on a general cell with
   * non-constant Jacobian.
   */
  const Tensor<1,(dim>1?dim*(dim-1)/2:1),Tensor<1,dim,VectorizedArray<Number> > > * jacobian_grad_upper;

  /**
   * After a call to reinit(), stores the number of the cell we are currently
   * working with.
   */
  unsigned int cell;

  /**
   * Stores the type of the cell we are currently working with after a call to
   * reinit(). Valid values are @p cartesian, @p affine and @p general, which
   * have different implications on how the Jacobian transformations are
   * stored internally in MappingInfo.
   */
  internal::MatrixFreeFunctions::CellType cell_type;

  /**
   * The stride to access the correct data in MappingInfo.
   */
  unsigned int cell_data_number;

  /**
   * Debug information to track whether dof values have been initialized
   * before accessed. Used to control exceptions when uninitialized data is
   * used.
   */
  bool dof_values_initialized;

  /**
   * Debug information to track whether values on quadrature points have been
   * initialized before accessed. Used to control exceptions when
   * uninitialized data is used.
   */
  bool values_quad_initialized;

  /**
   * Debug information to track whether gradients on quadrature points have
   * been initialized before accessed. Used to control exceptions when
   * uninitialized data is used.
   */
  bool gradients_quad_initialized;

  /**
   * Debug information to track whether Hessians on quadrature points have
   * been initialized before accessed. Used to control exceptions when
   * uninitialized data is used.
   */
  bool hessians_quad_initialized;

  /**
   * Debug information to track whether values on quadrature points have been
   * submitted for integration before the integration is actually stared. Used
   * to control exceptions when uninitialized data is used.
   */
  bool values_quad_submitted;

  /**
   * Debug information to track whether gradients on quadrature points have
   * been submitted for integration before the integration is actually
   * stared. Used to control exceptions when uninitialized data is used.
   */
  bool gradients_quad_submitted;

  /**
   * Geometry data generated by FEValues on the fly.
   */
  SmartPointer<const MappingFEEvaluation<dim,Number> > mapped_geometry;

  /**
   * A pointer to the underlying DoFHandler.
   */
  const DoFHandler<dim> *dof_handler;

  /**
   * For a DoFHandler with more than one finite element, select at which
   * component this data structure should start.
   */
  const unsigned int first_selected_component;

  /**
   * A temporary data structure necessary to read degrees of freedom when no
   * MatrixFree object was given at initialization.
   */
  mutable std::vector<types::global_dof_index> local_dof_indices;
};



/**
 * This class provides access to the data fields of the FEEvaluation
 * classes. Generic access is achieved through the base class, and
 * specializations for scalar and vector-valued elements are defined
 * separately.
 *
 * @author Katharina Kormann and Martin Kronbichler, 2010, 2011
 */
template <int dim, int n_components_, typename Number>
class FEEvaluationAccess : public FEEvaluationBase<dim,n_components_,Number>
{
public:
  typedef Number                            number_type;
  typedef Tensor<1,n_components_,VectorizedArray<Number> > value_type;
  typedef Tensor<1,n_components_,Tensor<1,dim,VectorizedArray<Number> > > gradient_type;
  static const unsigned int dimension     = dim;
  static const unsigned int n_components  = n_components_;
  typedef FEEvaluationBase<dim,n_components_, Number> BaseClass;

protected:
  /**
   * Constructor. Made protected to prevent initialization in user code. Takes
   * all data stored in MatrixFree. If applied to problems with more than one
   * finite element or more than one quadrature formula selected during
   * construction of @p matrix_free, @p fe_no and @p quad_no allow to select
   * the appropriate components.
   */
  FEEvaluationAccess (const MatrixFree<dim,Number> &matrix_free,
                      const unsigned int            fe_no,
                      const unsigned int            quad_no,
                      const unsigned int            fe_degree,
                      const unsigned int            n_q_points);

  /**
   * Constructor that comes with reduced functionality and works similar as
   * FEValues. The user has to provide a structure of type MappingFEEvaluation
   * and a DoFHandler in order to allow for reading out the finite element
   * data. It uses the data provided by dof_handler.get_fe(). If the element
   * is vector-valued, the optional argument allows to specify the index of
   * the base element (as long as the element is primitive, non-primitive are
   * not supported currently).
   *
   * With initialization from a FEValues object, no call to a reinit method of
   * this class is necessary. Instead, it is enough if the geometry is
   * initialized to a given cell iterator. It can also read from or write to
   * vectors in the standard way for DoFHandler<dim>::active_cell_iterator
   * types (which is less efficient with MPI since index translation has to be
   * done), but of course only for one cell at a time. Hence, a kernel using
   * this method does not vectorize over several elements (which is most
   * efficient for vector operations), but only possibly within the element if
   * the evaluate/integrate routines are combined (e.g. for computing cell
   * matrices).
   * With this initialization, no call to a reinit method of this
   * class. Instead, it is enough if the geometry is initialized to a given
   * cell iterator. Moreover, beware that a kernel using this method does not
   * vectorize over several elements (which is most efficient for vector
   * operations), but only possibly within the element if the
   * evaluate/integrate routines are combined (e.g. for matrix assembly).
   */
  FEEvaluationAccess (const MappingFEEvaluation<dim,Number> &geometry,
                      const DoFHandler<dim>                 &dof_handler,
                      const unsigned int                     first_selected_component = 0);

  /**
   * Copy constructor
   */
  FEEvaluationAccess (const FEEvaluationAccess &other);
};




/**
 * This class provides access to the data fields of the FEEvaluation
 * classes. Partial specialization for scalar fields that defines access with
 * simple data fields, i.e., scalars for the values and Tensor<1,dim> for the
 * gradients.
 *
 * @author Katharina Kormann and Martin Kronbichler, 2010, 2011
 */
template <int dim, typename Number>
class FEEvaluationAccess<dim,1,Number> : public FEEvaluationBase<dim,1,Number>
{
public:
  typedef Number                                 number_type;
  typedef VectorizedArray<Number>                value_type;
  typedef Tensor<1,dim,VectorizedArray<Number> > gradient_type;
  static const unsigned int dimension          = dim;
  typedef FEEvaluationBase<dim,1,Number>         BaseClass;

  /**
   * Returns the value stored for the local degree of freedom with index @p
   * dof. If the object is vector-valued, a vector-valued return argument is
   * given. Note that when vectorization is enabled, values from several cells
   * are grouped together. If @p set_dof_values was called last, the value
   * corresponds to the one set there. If @p integrate was called last, it
   * instead corresponds to the value of the integrated function with the test
   * function of the given index.
   */
  value_type get_dof_value (const unsigned int dof) const;

  /**
   * Write a value to the field containing the degrees of freedom with
   * component @p dof. Access to the same field as through @p get_dof_value.
   */
  void submit_dof_value (const value_type   val_in,
                         const unsigned int dof);

  /**
   * Returns the value of a finite element function at quadrature point number
   * @p q_point after a call to @p evaluate(true,...), or the value that has
   * been stored there with a call to @p submit_value. If the object is
   * vector-valued, a vector-valued return argument is given. Note that when
   * vectorization is enabled, values from several cells are grouped together.
   */
  value_type get_value (const unsigned int q_point) const;

  /**
   * Write a value to the field containing the values on quadrature points
   * with component @p q_point. Access to the same field as through @p
   * get_value. If applied before the function @p integrate(true,...) is
   * called, this specifies the value which is tested by all basis function on
   * the current cell and integrated over.
   */
  void submit_value (const value_type   val_in,
                     const unsigned int q_point);

  /**
   * Returns the gradient of a finite element function at quadrature point
   * number @p q_point after a call to @p evaluate(...,true,...), or the value
   * that has been stored there with a call to @p submit_gradient.
   */
  gradient_type get_gradient (const unsigned int q_point) const;

  /**
   * Write a contribution that is tested by the gradient to the field
   * containing the values on quadrature points with component @p
   * q_point. Access to the same field as through @p get_gradient. If applied
   * before the function @p integrate(...,true) is called, this specifies what
   * is tested by all basis function gradients on the current cell and
   * integrated over.
   */
  void submit_gradient(const gradient_type grad_in,
                       const unsigned int  q_point);

  /**
   * Returns the Hessian of a finite element function at quadrature point
   * number @p q_point after a call to @p evaluate(...,true). If only the
   * diagonal part of the Hessian or its trace, the Laplacian, are needed, use
   * the respective functions below.
   */
  Tensor<2,dim,VectorizedArray<Number> >
  get_hessian (unsigned int q_point) const;

  /**
   * Returns the diagonal of the Hessian of a finite element function at
   * quadrature point number @p q_point after a call to @p evaluate(...,true).
   */
  gradient_type get_hessian_diagonal (const unsigned int q_point) const;

  /**
   * Returns the Laplacian of a finite element function at quadrature point
   * number @p q_point after a call to @p evaluate(...,true).
   */
  value_type get_laplacian (const unsigned int q_point) const;

  /**
   * Takes values on quadrature points, multiplies by the Jacobian determinant
   * and quadrature weights (JxW) and sums the values for all quadrature
   * points on the cell. The result is a scalar, representing the integral
   * over the function over the cell. If a vector-element is used, the
   * resulting components are still separated. Moreover, if vectorization is
   * enabled, the integral values of several cells are represented together.
   */
  value_type integrate_value () const;

protected:
  /**
   * Constructor. Made protected to avoid initialization in user code. Takes
   * all data stored in MatrixFree. If applied to problems with more than one
   * finite element or more than one quadrature formula selected during
   * construction of @p matrix_free, @p fe_no and @p quad_no allow to select
   * the appropriate components.
   */
  FEEvaluationAccess (const MatrixFree<dim,Number> &matrix_free,
                      const unsigned int            fe_no,
                      const unsigned int            quad_no,
                      const unsigned int            fe_degree,
                      const unsigned int            n_q_points);

  /**
   * Constructor that comes with reduced functionality and works similar as
   * FEValues. The user has to provide a structure of type MappingFEEvaluation
   * and a DoFHandler in order to allow for reading out the finite element
   * data. It uses the data provided by dof_handler.get_fe(). If the element
   * is vector-valued, the optional argument allows to specify the index of
   * the base element (as long as the element is primitive, non-primitive are
   * not supported currently).
   *
   * With initialization from a FEValues object, no call to a reinit method of
   * this class is necessary. Instead, it is enough if the geometry is
   * initialized to a given cell iterator. It can also read from or write to
   * vectors in the standard way for DoFHandler<dim>::active_cell_iterator
   * types (which is less efficient with MPI since index translation has to be
   * done), but of course only for one cell at a time. Hence, a kernel using
   * this method does not vectorize over several elements (which is most
   * efficient for vector operations), but only possibly within the element if
   * the evaluate/integrate routines are combined (e.g. for computing cell
   * matrices).
   * With this initialization, no call to a reinit method of this
   * class. Instead, it is enough if the geometry is initialized to a given
   * cell iterator. Moreover, beware that a kernel using this method does not
   * vectorize over several elements (which is most efficient for vector
   * operations), but only possibly within the element if the
   * evaluate/integrate routines are combined (e.g. for matrix assembly).
   */
  FEEvaluationAccess (const MappingFEEvaluation<dim,Number> &geometry,
                      const DoFHandler<dim>                 &dof_handler,
                      const unsigned int                     first_selected_component = 0);

  /**
   * Copy constructor
   */
  FEEvaluationAccess (const FEEvaluationAccess &other);
};



/**
 * This class provides access to the data fields of the FEEvaluation
 * classes. Partial specialization for fields with as many components as the
 * underlying space dimension, i.e., values are of type Tensor<1,dim> and
 * gradients of type Tensor<2,dim>. Provides some additional functions for
 * access, like the symmetric gradient and divergence.
 *
 * @author Katharina Kormann and Martin Kronbichler, 2010, 2011
 */
template <int dim, typename Number>
class FEEvaluationAccess<dim,dim,Number> : public FEEvaluationBase<dim,dim,Number>
{
public:
  typedef Number                            number_type;
  typedef Tensor<1,dim,VectorizedArray<Number> > value_type;
  typedef Tensor<2,dim,VectorizedArray<Number> > gradient_type;
  static const unsigned int dimension     = dim;
  static const unsigned int n_components  = dim;
  typedef FEEvaluationBase<dim,dim,Number> BaseClass;

  /**
   * Returns the gradient of a finite element function at quadrature point
   * number @p q_point after a call to @p evaluate(...,true,...).
   */
  gradient_type get_gradient (const unsigned int q_point) const;

  /**
   * Returns the divergence of a vector-valued finite element at quadrature
   * point number @p q_point after a call to @p evaluate(...,true,...).
   */
  VectorizedArray<Number> get_divergence (const unsigned int q_point) const;

  /**
   * Returns the symmetric gradient of a vector-valued finite element at
   * quadrature point number @p q_point after a call to @p
   * evaluate(...,true,...). It corresponds to <tt>0.5
   * (grad+grad<sup>T</sup>)</tt>.
   */
  SymmetricTensor<2,dim,VectorizedArray<Number> >
  get_symmetric_gradient (const unsigned int q_point) const;

  /**
   * Returns the curl of the vector field, $nabla \times v$ after a call to @p
   * evaluate(...,true,...).
   */
  Tensor<1,dim==2?1:dim,VectorizedArray<Number> >
  get_curl (const unsigned int q_point) const;

  /**
   * Returns the Hessian of a finite element function at quadrature point
   * number @p q_point after a call to @p evaluate(...,true). If only the
   * diagonal of the Hessian or its trace, the Laplacian, is needed, use the
   * respective functions.
   */
  Tensor<3,dim,VectorizedArray<Number> >
  get_hessian (const unsigned int q_point) const;

  /**
   * Returns the diagonal of the Hessian of a finite element function at
   * quadrature point number @p q_point after a call to @p evaluate(...,true).
   */
  gradient_type get_hessian_diagonal (const unsigned int q_point) const;

  /**
   * Write a contribution that is tested by the gradient to the field
   * containing the values on quadrature points with component @p
   * q_point. Access to the same field as through @p get_gradient. If applied
   * before the function @p integrate(...,true) is called, this specifies what
   * is tested by all basis function gradients on the current cell and
   * integrated over.
   */
  void submit_gradient(const gradient_type grad_in,
                       const unsigned int  q_point);

  /**
   * Write a contribution that is tested by the gradient to the field
   * containing the values on quadrature points with component @p
   * q_point. This function is an alternative to the other submit_gradient
   * function when using a system of fixed number of equations which happens
   * to coincide with the dimension for some dimensions, but not all. To allow
   * for dimension-independent programming, this function can be used instead.
   */
  void submit_gradient(const Tensor<1,dim,Tensor<1,dim,VectorizedArray<Number> > > grad_in,
                       const unsigned int q_point);

  /**
   * Write a constribution that is tested by the divergence to the field
   * containing the values on quadrature points with component @p
   * q_point. Access to the same field as through @p get_gradient. If applied
   * before the function @p integrate(...,true) is called, this specifies what
   * is tested by all basis function gradients on the current cell and
   * integrated over.
   */
  void submit_divergence (const VectorizedArray<Number> div_in,
                          const unsigned int q_point);

  /**
   * Write a contribution that is tested by the gradient to the field
   * containing the values on quadrature points with component @p
   * q_point. Access to the same field as through @p get_gradient. If applied
   * before the function @p integrate(...,true) is called, this specifies the
   * gradient which is tested by all basis function gradients on the current
   * cell and integrated over.
   */
  void submit_symmetric_gradient(const SymmetricTensor<2,dim,VectorizedArray<Number> > grad_in,
                                 const unsigned int      q_point);

  /**
   * Write the components of a curl containing the values on quadrature point
   * @p q_point. Access to the same data field as through @p get_gradient.
   */
  void submit_curl (const Tensor<1,dim==2?1:dim,VectorizedArray<Number> > curl_in,
                    const unsigned int q_point);

protected:
  /**
   * Constructor. Made protected to avoid initialization in user code. Takes
   * all data stored in MatrixFree. If applied to problems with more than one
   * finite element or more than one quadrature formula selected during
   * construction of @p matrix_free, @p fe_no and @p quad_no allow to select
   * the appropriate components.
   */
  FEEvaluationAccess (const MatrixFree<dim,Number> &matrix_free,
                      const unsigned int            fe_no,
                      const unsigned int            quad_no,
                      const unsigned int            dofs_per_cell,
                      const unsigned int            n_q_points);

  /**
   * Constructor that comes with reduced functionality and works similar as
   * FEValues. The user has to provide a structure of type MappingFEEvaluation
   * and a DoFHandler in order to allow for reading out the finite element
   * data. It uses the data provided by dof_handler.get_fe(). If the element
   * is vector-valued, the optional argument allows to specify the index of
   * the base element (as long as the element is primitive, non-primitive are
   * not supported currently).
   *
   * With initialization from a FEValues object, no call to a reinit method of
   * this class is necessary. Instead, it is enough if the geometry is
   * initialized to a given cell iterator. It can also read from or write to
   * vectors in the standard way for DoFHandler<dim>::active_cell_iterator
   * types (which is less efficient with MPI since index translation has to be
   * done), but of course only for one cell at a time. Hence, a kernel using
   * this method does not vectorize over several elements (which is most
   * efficient for vector operations), but only possibly within the element if
   * the evaluate/integrate routines are combined (e.g. for computing cell
   * matrices).
   * With this initialization, no call to a reinit method of this
   * class. Instead, it is enough if the geometry is initialized to a given
   * cell iterator. Moreover, beware that a kernel using this method does not
   * vectorize over several elements (which is most efficient for vector
   * operations), but only possibly within the element if the
   * evaluate/integrate routines are combined (e.g. for matrix assembly).
   */
  FEEvaluationAccess (const MappingFEEvaluation<dim,Number> &geometry,
                      const DoFHandler<dim>                 &dof_handler,
                      const unsigned int                     first_selected_component = 0);

  /**
   * Copy constructor
   */
  FEEvaluationAccess (const FEEvaluationAccess &other);
};



/**
 * The class that provides all functions necessary to evaluate functions at
 * quadrature points and cell integrations. In functionality, this class is
 * similar to FEValues<dim>, however, it includes a lot of specialized
 * functions that make it much faster (between 5 and 500, depending on the
 * polynomial order).
 *
 * This class can be used in two different ways. The first way is to
 * initialize it from a MatrixFree object that caches everything related to
 * the degrees of freedom and the mapping information. This way, it is
 * possible to use vectorization for applying a vector operation for several
 * cells at once. The second form of usage is to initialize it from geometry
 * information generated by FEValues, which is stored in the class
 * MappingFEEvaluation. Here, the operations can only work on a single cell, but
 * possibly be vectorized by combining several operations (e.g. when
 * performing matrix assembly).
 *
 * This class contains specialized evaluation routines for several elements,
 * including standard FE_Q or FE_DGQ elements and quadrature points symmetric
 * around 0.5 (like Gauss quadrature), FE_DGP elements based on truncated
 * tensor products as well as the faster case of Gauss-Lobatto elements with
 * Gauss-Lobatto quadrature which give diagonal mass matrices and quicker
 * evaluation internally. Note that many of the operations available through
 * this class are inherited from the base class FEEvaluationBase, in
 * particular reading from and writing to vectors. Also, the class inherits
 * from FEEvaluationAccess that implements access to values, gradients and
 * Hessians of the finite element function on quadrature points.
 *
 * This class assumes that shape functions of the FiniteElement under
 * consideration do <em>not</em> depend on the actual shape of the cells in
 * real space. Currently, other finite elements cannot be treated with the
 * matrix-free concept.
 *
 * This class has five template arguments:
 *
 * @param dim Dimension in which this class is to be used
 *
 * @param fe_degree Degree of the tensor product finite element with
 *                  fe_degree+1 degrees of freedom per coordinate direction
 *
 * @param n_q_points_1d Number of points in the quadrature formula in 1D,
 *                   usually chosen as fe_degree+1
 *
 * @param n_components Number of vector components when solving a system of
 *                  PDEs. If the same operation is applied to several
 *                  components of a PDE (e.g. a vector Laplace equation), they
 *                  can be applied simultaneously with one call (and often
 *                  more efficiently)
 *
 * @param Number Number format, usually @p double or @p float
 *
 * @author Katharina Kormann and Martin Kronbichler, 2010, 2011
 */
template <int dim, int fe_degree, int n_q_points_1d = fe_degree+1,
          int n_components_ = 1, typename Number = double >
class FEEvaluation : public FEEvaluationAccess<dim,n_components_,Number>
{
public:
  typedef FEEvaluationAccess<dim,n_components_,Number> BaseClass;
  typedef Number                            number_type;
  typedef typename BaseClass::value_type    value_type;
  typedef typename BaseClass::gradient_type gradient_type;
  static const unsigned int dimension     = dim;
  static const unsigned int n_components  = n_components_;
  static const unsigned int n_q_points    = Utilities::fixed_int_power<n_q_points_1d,dim>::value;
  static const unsigned int tensor_dofs_per_cell = Utilities::fixed_int_power<fe_degree+1,dim>::value;

  /**
   * Constructor. Takes all data stored in MatrixFree. If applied to problems
   * with more than one finite element or more than one quadrature formula
   * selected during construction of @p matrix_free, @p fe_no and @p quad_no
   * allow to select the appropriate components.
   */
  FEEvaluation (const MatrixFree<dim,Number> &matrix_free,
                const unsigned int            fe_no   = 0,
                const unsigned int            quad_no = 0);

  /**
   * Copy constructor
   */
  FEEvaluation (const FEEvaluation &other);

  /**
   * Constructor that comes with reduced functionality and works similar as
   * FEValues. The user has to provide a structure of type MappingFEEvaluation
   * and a DoFHandler in order to allow for reading out the finite element
   * data. It uses the data provided by dof_handler.get_fe(). If the element
   * is vector-valued, the optional argument allows to specify the index of
   * the base element (as long as the element is primitive, non-primitive are
   * not supported currently).
   *
   * With initialization from a FEValues object, no call to a reinit method of
   * this class is necessary. Instead, it is enough if the geometry is
   * initialized to a given cell iterator. It can also read from or write to
   * vectors in the standard way for DoFHandler<dim>::active_cell_iterator
   * types (which is less efficient with MPI since index translation has to be
   * done), but of course only for one cell at a time. Hence, a kernel using
   * this method does not vectorize over several elements (which is most
   * efficient for vector operations), but only possibly within the element if
   * the evaluate/integrate routines are combined (e.g. for computing cell
   * matrices).
   */
  FEEvaluation (const MappingFEEvaluation<dim,Number> &geometry,
                const DoFHandler<dim>                 &dof_handler,
                const unsigned int                     first_selected_component = 0);

  /**
   * Evaluates the function values, the gradients, and the Laplacians of the
   * FE function given at the DoF values in the input vector at the quadrature
   * points on the unit cell.  The function arguments specify which parts
   * shall actually be computed. Needs to be called before the functions @p
   * get_value(), @p get_gradient() or @p get_laplacian give useful
   * information (unless these values have been set manually).
   */
  void evaluate (const bool evaluate_val,
                 const bool evaluate_grad,
                 const bool evaluate_hess = false);

  /**
   * This function takes the values and/or gradients that are stored on
   * quadrature points, tests them by all the basis functions/gradients on the
   * cell and performs the cell integration. The two function arguments @p
   * integrate_val and @p integrate_grad are used to enable/disable some of
   * values or gradients.
   */
  void integrate (const bool integrate_val,
                  const bool integrate_grad);

  /**
   * Returns the q-th quadrature point stored in MappingInfo.
   */
  Point<dim,VectorizedArray<Number> >
  quadrature_point (const unsigned int q_point) const;

  /**
   * The number of scalar degrees of freedom on the cell. Usually close to
   * tensor_dofs_per_cell, but depends on the actual element selected and is
   * thus not static.
   */
  const unsigned int dofs_per_cell;

private:
  /**
   * Internally stored variables for the different data fields.
   */
  VectorizedArray<Number> my_data_array[n_components*(tensor_dofs_per_cell+1+(dim*dim+2*dim+1)*n_q_points)];

  /**
   * Checks if the template arguments regarding degree of the element
   * corresponds to the actual element used at initialization.
   */
  void check_template_arguments(const unsigned int fe_no);

  /**
   * Sets the pointers of the base class to my_data_array.
   */
  void set_data_pointers();

  /**
   * Function pointer for the evaluate function
   */
  void (*evaluate_funct) (const internal::MatrixFreeFunctions::ShapeInfo<Number> &,
                          VectorizedArray<Number> *values_dofs_actual[],
                          VectorizedArray<Number> *values_quad[],
                          VectorizedArray<Number> *gradients_quad[][dim],
                          VectorizedArray<Number> *hessians_quad[][(dim*(dim+1))/2],
                          const bool               evaluate_val,
                          const bool               evaluate_grad,
                          const bool               evaluate_lapl);

  /**
   * Function pointer for the integrate function
   */
  void (*integrate_funct)(const internal::MatrixFreeFunctions::ShapeInfo<Number> &,
                          VectorizedArray<Number> *values_dofs_actual[],
                          VectorizedArray<Number> *values_quad[],
                          VectorizedArray<Number> *gradients_quad[][dim],
                          const bool               evaluate_val,
                          const bool               evaluate_grad);
};



/**
 * Deprecated. Functionality has been merged into FEEvaluation. Use class
 * FEEvaluation instead.
 */
template <int dim, int fe_degree, int n_q_points_1d = fe_degree+1,
          int n_components_ = 1, typename Number = double >
class FEEvaluationGeneral : public FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_,Number>
{
public:
  typedef FEEvaluation<dim, fe_degree, n_q_points_1d, n_components_, Number> BaseClass;
  typedef Number                            number_type;
  typedef typename BaseClass::value_type    value_type;
  typedef typename BaseClass::gradient_type gradient_type;
  static const unsigned int dimension     = dim;
  static const unsigned int n_components  = n_components_;
  static const unsigned int dofs_per_cell = Utilities::fixed_int_power<fe_degree+1,dim>::value;
  static const unsigned int n_q_points    = Utilities::fixed_int_power<n_q_points_1d,dim>::value;

  /**
   * Constructor.
   */
  FEEvaluationGeneral (const MatrixFree<dim,Number> &matrix_free,
                       const unsigned int            fe_no   = 0,
                       const unsigned int            quad_no = 0) DEAL_II_DEPRECATED
:
  BaseClass (matrix_free, fe_no, quad_no)
  {}

  /**
   * Constructor.
   */
  FEEvaluationGeneral (const MappingFEEvaluation<dim,Number> &geometry,
                       const DoFHandler<dim>                 &dof_handler,
                       const unsigned int                     first_selected_component = 0) DEAL_II_DEPRECATED
:
  BaseClass (geometry, dof_handler, first_selected_component)
  {}
};



/**
 * Deprecated. Functionality has been merged into FEEvaluation. Use class
 * FEEvaluation instead.
 */
template <int dim, int fe_degree, int n_components_ = 1, typename Number = double >
class FEEvaluationGL :
  public FEEvaluation<dim,fe_degree,fe_degree+1,n_components_,Number>
{
public:
  typedef FEEvaluation<dim,fe_degree,fe_degree+1,n_components_,Number> BaseClass;
  typedef Number                            number_type;
  typedef typename BaseClass::value_type    value_type;
  typedef typename BaseClass::gradient_type gradient_type;
  static const unsigned int dimension     = dim;
  static const unsigned int n_components  = n_components_;
  static const unsigned int dofs_per_cell = Utilities::fixed_int_power<fe_degree+1,dim>::value;
  static const unsigned int n_q_points    = BaseClass::n_q_points;

  /**
   * Constructor.
   */
  FEEvaluationGL (const MatrixFree<dim,Number> &matrix_free,
                  const unsigned int          fe_no   = 0,
                  const unsigned int          quad_no = 0) DEAL_II_DEPRECATED
:
  BaseClass (matrix_free, fe_no, quad_no)
  {}

  /**
   * Constructor.
   */
  FEEvaluationGL (const MappingFEEvaluation<dim,Number> &geometry,
                  const DoFHandler<dim>                 &dof_handler,
                  const unsigned int                     first_selected_component = 0) DEAL_II_DEPRECATED
:
  BaseClass (geometry, dof_handler, first_selected_component)
  {}
};



namespace internal
{
  namespace MatrixFreeFunctions
  {
    // a helper function to compute the number of DoFs of a DGP element at compile
    // time, depending on the degree
    template <int dim, int degree>
    struct DGP_dofs_per_cell
    {
      // this division is always without remainder
      static const unsigned int value =
        (DGP_dofs_per_cell<dim-1,degree>::value * (degree+dim)) / dim;
    };

    // base specialization: 1d elements have 'degree+1' degrees of freedom
    template <int degree>
    struct DGP_dofs_per_cell<1,degree>
    {
      static const unsigned int value = degree+1;
    };
  }
}



/**
 * Functionality has been merged into FEEvaluation. Use class FEEvaluation
 * instead.
 */
template <int dim, int fe_degree, int n_q_points_1d = fe_degree+1,
          int n_components_ = 1, typename Number = double >
class FEEvaluationDGP :
  public FEEvaluation<dim,fe_degree,n_q_points_1d,n_components_,Number>
{
public:
  typedef FEEvaluation<dim,fe_degree,n_q_points_1d,n_components_,Number> BaseClass;
  typedef Number                            number_type;
  typedef typename BaseClass::value_type    value_type;
  typedef typename BaseClass::gradient_type gradient_type;
  static const unsigned int dimension     = dim;
  static const unsigned int n_components  = n_components_;
  static const unsigned int dofs_per_cell = internal::MatrixFreeFunctions::DGP_dofs_per_cell<dim,fe_degree>::value;
  static const unsigned int n_q_points    = BaseClass::n_q_points;

  /**
   * Constructor.
   */
  FEEvaluationDGP (const MatrixFree<dim,Number> &matrix_free,
                   const unsigned int            fe_no   = 0,
                   const unsigned int            quad_no = 0) DEAL_II_DEPRECATED
:
  BaseClass (matrix_free, fe_no, quad_no)
  {}

  /**
   * Constructor.
   */
  FEEvaluationDGP (const MappingFEEvaluation<dim,Number> &geometry,
                   const DoFHandler<dim>                 &dof_handler,
                   const unsigned int                     first_selected_component = 0) DEAL_II_DEPRECATED
:
  BaseClass (geometry, dof_handler, first_selected_component)
  {}
};



/*----------------------- Inline functions ----------------------------------*/

#ifndef DOXYGEN



/*----------------------- FEEvaluationBase ----------------------------------*/

template <int dim, int n_components_, typename Number>
inline
FEEvaluationBase<dim,n_components_,Number>
::FEEvaluationBase (const MatrixFree<dim,Number> &data_in,
                    const unsigned int fe_no_in,
                    const unsigned int quad_no_in,
                    const unsigned int fe_degree,
                    const unsigned int n_q_points)
  :
  quad_no            (quad_no_in),
  n_fe_components    (data_in.get_dof_info(fe_no_in).n_components),
  active_fe_index    (data_in.get_dof_info(fe_no_in).fe_index_from_degree
                      (fe_degree)),
  active_quad_index  (data_in.get_mapping_info().
                      mapping_data_gen[quad_no_in].
                      quad_index_from_n_q_points(n_q_points)),
  matrix_info        (&data_in),
  dof_info           (&data_in.get_dof_info(fe_no_in)),
  mapping_info       (&data_in.get_mapping_info()),
  data               (&data_in.get_shape_info
                      (fe_no_in, quad_no_in, active_fe_index,
                       active_quad_index)),
  cartesian_data     (0),
  jacobian           (0),
  J_value            (0),
  quadrature_weights (mapping_info->mapping_data_gen[quad_no].
                      quadrature_weights[active_quad_index].begin()),
  quadrature_points  (0),
  jacobian_grad      (0),
  jacobian_grad_upper(0),
  cell               (numbers::invalid_unsigned_int),
  cell_type          (internal::MatrixFreeFunctions::undefined),
  cell_data_number   (0),
  dof_handler        (0),
  first_selected_component (0)
{
  for (unsigned int c=0; c<n_components_; ++c)
    {
      values_dofs[c] = 0;
      values_quad[c] = 0;
      for (unsigned int d=0; d<dim; ++d)
        gradients_quad[c][d] = 0;
      for (unsigned int d=0; d<(dim*dim+dim)/2; ++d)
        hessians_quad[c][d] = 0;
    }
  Assert (matrix_info->mapping_initialized() == true,
          ExcNotInitialized());
  AssertDimension (matrix_info->get_size_info().vectorization_length,
                   VectorizedArray<Number>::n_array_elements);
  AssertDimension (data->dofs_per_cell,
                   dof_info->dofs_per_cell[active_fe_index]/n_fe_components);
  AssertDimension (data->n_q_points,
                   mapping_info->mapping_data_gen[quad_no].n_q_points[active_quad_index]);
  Assert (n_fe_components == 1 ||
          n_components == 1 ||
          n_components == n_fe_components,
          ExcMessage ("The underlying FE is vector-valued. In this case, the "
                      "template argument n_components must be a the same "
                      "as the number of underlying vector components."));


  // do not check for correct dimensions of data fields here, should be done
  // in derived classes
}



template <int dim, int n_components_, typename Number>
inline
FEEvaluationBase<dim,n_components_,Number>
::FEEvaluationBase (const MappingFEEvaluation<dim,Number> &geometry,
                    const DoFHandler<dim>                 &dof_handler_in,
                    const unsigned int                     first_selected_component)
  :
  quad_no            (-1),
  n_fe_components    (n_components_),
  active_fe_index    (-1),
  active_quad_index  (-1),
  matrix_info        (0),
  dof_info           (0),
  mapping_info       (0),
  stored_shape_info  (new internal::MatrixFreeFunctions::ShapeInfo<Number>(geometry.get_quadrature(), dof_handler_in.get_fe(), dof_handler_in.get_fe().component_to_base_index(first_selected_component).first)),
  data               (stored_shape_info.get()),
  cartesian_data     (0),
  jacobian           (geometry.get_inverse_jacobians().begin()),
  J_value            (geometry.get_JxW_values().begin()),
  quadrature_weights (0),
  quadrature_points  (geometry.get_quadrature_points().begin()),
  jacobian_grad      (0),
  jacobian_grad_upper(0),
  cell               (0),
  cell_type          (internal::MatrixFreeFunctions::general),
  cell_data_number   (0),
  mapped_geometry    (&geometry),
  dof_handler        (&dof_handler_in),
  first_selected_component (first_selected_component)
{
  const unsigned int base_element_number =
    dof_handler_in.get_fe().component_to_base_index(first_selected_component).first;
  for (unsigned int c=0; c<n_components_; ++c)
    {
      values_dofs[c] = 0;
      values_quad[c] = 0;
      for (unsigned int d=0; d<dim; ++d)
        gradients_quad[c][d] = 0;
      for (unsigned int d=0; d<(dim*dim+dim)/2; ++d)
        hessians_quad[c][d] = 0;
    }
  Assert(dof_handler->get_fe().element_multiplicity(base_element_number) == 1 ||
         dof_handler->get_fe().element_multiplicity(base_element_number)-first_selected_component >= n_components_,
         ExcMessage("The underlying element must at least contain as many "
                    "components as requested by this class"));
}



template <int dim, int n_components_, typename Number>
inline
FEEvaluationBase<dim,n_components_,Number>
::FEEvaluationBase (const FEEvaluationBase<dim,n_components_,Number> &other)
  :
  quad_no            (other.quad_no),
  n_fe_components    (other.n_fe_components),
  active_fe_index    (other.active_fe_index),
  active_quad_index  (other.active_quad_index),
  matrix_info        (other.matrix_info),
  dof_info           (other.dof_info),
  mapping_info       (other.mapping_info),
  stored_shape_info  (other.stored_shape_info),
  data               (other.data),
  cartesian_data     (other.cartesian_data),
  jacobian           (other.jacobian),
  J_value            (other.J_value),
  quadrature_weights (other.quadrature_weights),
  quadrature_points  (other.quadrature_points),
  jacobian_grad      (other.jacobian_grad),
  jacobian_grad_upper(other.jacobian_grad_upper),
  cell               (other.cell),
  cell_type          (other.cell_type),
  cell_data_number   (other.cell_data_number),
  mapped_geometry    (other.mapped_geometry),
  dof_handler        (other.dof_handler),
  first_selected_component (other.first_selected_component)
{
  for (unsigned int c=0; c<n_components_; ++c)
    {
      values_dofs[c] = 0;
      values_quad[c] = 0;
      for (unsigned int d=0; d<dim; ++d)
        gradients_quad[c][d] = 0;
      for (unsigned int d=0; d<(dim*dim+dim)/2; ++d)
        hessians_quad[c][d] = 0;
    }
}



template <int dim, int n_components_, typename Number>
inline
void
FEEvaluationBase<dim,n_components_,Number>::reinit (const unsigned int cell_in)
{
  Assert (mapped_geometry == 0, ExcMessage("FEEvaluation was initialized without a matrix-free object. Integer indexing is not possible"));
  if (mapped_geometry != 0)
    return;
  Assert (dof_info != 0, ExcNotInitialized());
  Assert (mapping_info != 0, ExcNotInitialized());
  AssertIndexRange (cell_in, dof_info->row_starts.size()-1);
  AssertDimension (((dof_info->cell_active_fe_index.size() > 0) ?
                    dof_info->cell_active_fe_index[cell_in] : 0),
                   active_fe_index);
  cell = cell_in;
  cell_type = mapping_info->get_cell_type(cell);
  cell_data_number = mapping_info->get_cell_data_index(cell);

  if (mapping_info->quadrature_points_initialized == true)
    {
      AssertIndexRange (cell_data_number, mapping_info->
                        mapping_data_gen[quad_no].rowstart_q_points.size());
      const unsigned int index = mapping_info->mapping_data_gen[quad_no].
                                 rowstart_q_points[cell];
      AssertIndexRange (index, mapping_info->mapping_data_gen[quad_no].
                        quadrature_points.size());
      quadrature_points =
        &mapping_info->mapping_data_gen[quad_no].quadrature_points[index];
    }

  if (cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      cartesian_data = &mapping_info->cartesian_data[cell_data_number].first;
      J_value        = &mapping_info->cartesian_data[cell_data_number].second;
    }
  else if (cell_type == internal::MatrixFreeFunctions::affine)
    {
      jacobian  = &mapping_info->affine_data[cell_data_number].first;
      J_value   = &mapping_info->affine_data[cell_data_number].second;
    }
  else
    {
      const unsigned int rowstart = mapping_info->
                                    mapping_data_gen[quad_no].rowstart_jacobians[cell_data_number];
      AssertIndexRange (rowstart, mapping_info->
                        mapping_data_gen[quad_no].jacobians.size());
      jacobian =
        &mapping_info->mapping_data_gen[quad_no].jacobians[rowstart];
      if (mapping_info->JxW_values_initialized == true)
        {
          AssertIndexRange (rowstart, mapping_info->
                            mapping_data_gen[quad_no].JxW_values.size());
          J_value = &(mapping_info->mapping_data_gen[quad_no].
                      JxW_values[rowstart]);
        }
      if (mapping_info->second_derivatives_initialized == true)
        {
          AssertIndexRange(rowstart, mapping_info->
                           mapping_data_gen[quad_no].jacobians_grad_diag.size());
          jacobian_grad = &mapping_info->mapping_data_gen[quad_no].
                          jacobians_grad_diag[rowstart];
          AssertIndexRange(rowstart, mapping_info->
                           mapping_data_gen[quad_no].jacobians_grad_upper.size());
          jacobian_grad_upper = &mapping_info->mapping_data_gen[quad_no].
                                jacobians_grad_upper[rowstart];
        }
    }

#ifdef DEBUG
  dof_values_initialized      = false;
  values_quad_initialized     = false;
  gradients_quad_initialized  = false;
  hessians_quad_initialized   = false;
#endif
}



template <int dim, int n_components_, typename Number>
inline
unsigned int
FEEvaluationBase<dim,n_components_,Number>
::get_cell_data_number () const
{
  Assert (cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  return cell_data_number;
}



template <int dim, int n_components_, typename Number>
inline
internal::MatrixFreeFunctions::CellType
FEEvaluationBase<dim,n_components_,Number>::get_cell_type () const
{
  Assert (cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  return cell_type;
}



template <int dim, int n_components_, typename Number>
inline
const internal::MatrixFreeFunctions::ShapeInfo<Number> &
FEEvaluationBase<dim,n_components_,Number>::get_shape_info() const
{
  Assert(data != 0, ExcInternalError());
  return *data;
}



template <int dim, int n_components_, typename Number>
inline
void
FEEvaluationBase<dim,n_components_,Number>
::fill_JxW_values(AlignedVector<VectorizedArray<Number> > &JxW_values) const
{
  AssertDimension(JxW_values.size(), data->n_q_points);
  Assert (this->J_value != 0, ExcNotImplemented());
  if (this->cell_type == internal::MatrixFreeFunctions::cartesian ||
      this->cell_type == internal::MatrixFreeFunctions::affine)
    {
      Assert (this->mapping_info != 0, ExcNotImplemented());
      VectorizedArray<Number> J = this->J_value[0];
      for (unsigned int q=0; q<this->data->n_q_points; ++q)
        JxW_values[q] = J * this->quadrature_weights[q];
    }
  else
    for (unsigned int q=0; q<data->n_q_points; ++q)
      JxW_values[q] = this->J_value[q];
}



namespace internal
{
  // write access to generic vectors that have operator ().
  template <typename VectorType>
  inline
  typename VectorType::value_type &
  vector_access (VectorType         &vec,
                 const unsigned int  entry)
  {
    return vec(entry);
  }



  // read access to generic vectors that have operator ().
  template <typename VectorType>
  inline
  typename VectorType::value_type
  vector_access (const VectorType   &vec,
                 const unsigned int  entry)
  {
    return vec(entry);
  }



  // write access to distributed MPI vectors that have a local_element(uint)
  // method to access data in local index space, which is what we use in
  // DoFInfo and hence in read_dof_values etc.
  template <typename Number>
  inline
  Number &
  vector_access (parallel::distributed::Vector<Number> &vec,
                 const unsigned int                     entry)
  {
    return vec.local_element(entry);
  }



  // read access to distributed MPI vectors that have a local_element(uint)
  // method to access data in local index space, which is what we use in
  // DoFInfo and hence in read_dof_values etc.
  template <typename Number>
  inline
  Number
  vector_access (const parallel::distributed::Vector<Number> &vec,
                 const unsigned int                           entry)
  {
    return vec.local_element(entry);
  }



  // this is to make sure that the parallel partitioning in the
  // parallel::distributed::Vector is really the same as stored in MatrixFree
  template <typename VectorType>
  inline
  void check_vector_compatibility (const VectorType                             &vec,
                                   const internal::MatrixFreeFunctions::DoFInfo &dof_info)
  {
    AssertDimension (vec.size(),
                     dof_info.vector_partitioner->size());
  }

  template <typename Number>
  inline
  void check_vector_compatibility (const parallel::distributed::Vector<Number>  &vec,
                                   const internal::MatrixFreeFunctions::DoFInfo &dof_info)
  {
    Assert (vec.partitioners_are_compatible(*dof_info.vector_partitioner),
            ExcMessage("The parallel layout of the given vector is not "
                       "compatible with the parallel partitioning in MatrixFree. "
                       "Use MatrixFree::initialize_dof_vector to get a "
                       "compatible vector."));
  }

  // A class to use the same code to read from and write to vector
  template <typename Number>
  struct VectorReader
  {
    template <typename VectorType>
    void process_dof (const unsigned int  index,
                      VectorType         &vec,
                      Number             &res) const
    {
      res = vector_access (const_cast<const VectorType &>(vec), index);
    }

    template <typename VectorType>
    void process_dof_global (const types::global_dof_index index,
                             VectorType         &vec,
                             Number             &res) const
    {
      res = const_cast<const VectorType &>(vec)(index);
    }

    void pre_constraints (const Number &,
                          Number       &res) const
    {
      res = Number();
    }

    template <typename VectorType>
    void process_constraint (const unsigned int index,
                             const Number       weight,
                             VectorType        &vec,
                             Number            &res) const
    {
      res += weight * vector_access (const_cast<const VectorType &>(vec), index);
    }

    void post_constraints (const Number &sum,
                           Number       &write_pos) const
    {
      write_pos = sum;
    }

    void process_empty (Number &res) const
    {
      res = Number();
    }
  };

  // A class to use the same code to read from and write to vector
  template <typename Number>
  struct VectorDistributorLocalToGlobal
  {
    template <typename VectorType>
    void process_dof (const unsigned int  index,
                      VectorType         &vec,
                      Number             &res) const
    {
      vector_access (vec, index) += res;
    }

    template <typename VectorType>
    void process_dof_global (const types::global_dof_index index,
                             VectorType         &vec,
                             Number             &res) const
    {
      vec(index) += res;
    }

    void pre_constraints (const Number &input,
                          Number       &res) const
    {
      res = input;
    }

    template <typename VectorType>
    void process_constraint (const unsigned int index,
                             const Number       weight,
                             VectorType        &vec,
                             Number            &res) const
    {
      vector_access (vec, index) += weight * res;
    }

    void post_constraints (const Number &,
                           Number &) const
    {
    }

    void process_empty (Number &) const
    {
    }
  };


  // A class to use the same code to read from and write to vector
  template <typename Number>
  struct VectorSetter
  {
    template <typename VectorType>
    void process_dof (const unsigned int  index,
                      VectorType         &vec,
                      Number             &res) const
    {
      vector_access (vec, index) = res;
    }

    template <typename VectorType>
    void process_dof_global (const types::global_dof_index index,
                             VectorType         &vec,
                             Number             &res) const
    {
      vec(index) = res;
    }

    void pre_constraints (const Number &,
                          Number &) const
    {
    }

    template <typename VectorType>
    void process_constraint (const unsigned int,
                             const Number,
                             VectorType &,
                             Number &) const
    {
    }

    void post_constraints (const Number &,
                           Number &) const
    {
    }

    void process_empty (Number &) const
    {
    }
  };

  // allows to select between block vectors and non-block vectors, which
  // allows to use a unified interface for extracting blocks on block vectors
  // and doing nothing on usual vectors
  template <typename VectorType, bool>
  struct BlockVectorSelector {};

  template <typename VectorType>
  struct BlockVectorSelector<VectorType,true>
  {
    typedef typename VectorType::BlockType BaseVectorType;

    static BaseVectorType *get_vector_component (VectorType &vec,
                                                 const unsigned int component)
    {
      AssertIndexRange (component, vec.n_blocks());
      return &vec.block(component);
    }
  };

  template <typename VectorType>
  struct BlockVectorSelector<VectorType,false>
  {
    typedef VectorType BaseVectorType;

    static BaseVectorType *get_vector_component (VectorType &vec,
                                                 const unsigned int)
    {
      return &vec;
    }
  };
}



template <int dim, int n_components_, typename Number>
template<typename VectorType, typename VectorOperation>
inline
void
FEEvaluationBase<dim,n_components_,Number>
::read_write_operation (const VectorOperation &operation,
                        VectorType            *src[]) const
{
  // This functions processes all the functions read_dof_values,
  // distribute_local_to_global, and set_dof_values with the same code. The
  // distinction between these three cases is made by the input
  // VectorOperation that either reads values from a vector and puts the data
  // into the local data field or write local data into the vector. Certain
  // operations are no-ops for the given use case.

  // Case 1: No MatrixFree object given, simple case because we do not need to
  // process constraints and need not care about vectorization
  if (matrix_info == 0)
    {
      Assert (dof_handler != 0, ExcNotInitialized());
      typename DoFHandler<dim>::cell_iterator cell (&dof_handler->get_tria(),
                                                    mapped_geometry->get_cell()->level(),
                                                    mapped_geometry->get_cell()->index(),
                                                    dof_handler);
      local_dof_indices.resize(dof_handler->get_fe().dofs_per_cell);
      cell->get_dof_indices(local_dof_indices);

      unsigned int index = first_selected_component * this->data->dofs_per_cell;
      for (unsigned int comp = 0; comp<n_components; ++comp)
        {
          for (unsigned int i=0; i<this->data->dofs_per_cell; ++i, ++index)
            {
              operation.process_dof_global(local_dof_indices[this->data->lexicographic_numbering[index]],
                                           *src[0], values_dofs[comp][i][0]);
              for (unsigned int v=1; v<VectorizedArray<Number>::n_array_elements; ++v)
                operation.process_empty(values_dofs[comp][i][v]);
            }
        }
      return;
    }

  Assert (dof_info != 0, ExcNotInitialized());
  Assert (matrix_info->indices_initialized() == true,
          ExcNotInitialized());
  Assert (cell != numbers::invalid_unsigned_int, ExcNotInitialized());

  // loop over all local dofs. ind_local holds local number on cell, index
  // iterates over the elements of index_local_to_global and dof_indices
  // points to the global indices stored in index_local_to_global
  const unsigned int *dof_indices = dof_info->begin_indices(cell);
  const std::pair<unsigned short,unsigned short> *indicators =
    dof_info->begin_indicators(cell);
  const std::pair<unsigned short,unsigned short> *indicators_end =
    dof_info->end_indicators(cell);
  unsigned int ind_local = 0;
  const unsigned int dofs_per_cell = this->data->dofs_per_cell;

  const unsigned int n_irreg_components_filled = dof_info->row_starts[cell][2];
  const bool at_irregular_cell = n_irreg_components_filled > 0;

  // scalar case (or case when all components have the same degrees of freedom
  // and sit on a different vector each)
  if (n_fe_components == 1)
    {
      const unsigned int n_local_dofs =
        VectorizedArray<Number>::n_array_elements * dofs_per_cell;
      for (unsigned int comp=0; comp<n_components; ++comp)
        internal::check_vector_compatibility (*src[comp], *dof_info);
      Number *local_data [n_components];
      for (unsigned int comp=0; comp<n_components; ++comp)
        local_data[comp] =
          const_cast<Number *>(&values_dofs[comp][0][0]);

      // standard case where there are sufficiently many cells to fill all
      // vectors
      if (at_irregular_cell == false)
        {
          // check whether there is any constraint on the current cell
          if (indicators != indicators_end)
            {
              for ( ; indicators != indicators_end; ++indicators)
                {
                  // run through values up to next constraint
                  for (unsigned int j=0; j<indicators->first; ++j)
                    for (unsigned int comp=0; comp<n_components; ++comp)
                      operation.process_dof (dof_indices[j], *src[comp],
                                             local_data[comp][ind_local+j]);

                  ind_local += indicators->first;
                  dof_indices   += indicators->first;

                  // constrained case: build the local value as a linear
                  // combination of the global value according to constraints
                  Number value [n_components];
                  for (unsigned int comp=0; comp<n_components; ++comp)
                    operation.pre_constraints (local_data[comp][ind_local],
                                               value[comp]);

                  const Number *data_val =
                    matrix_info->constraint_pool_begin(indicators->second);
                  const Number *end_pool =
                    matrix_info->constraint_pool_end(indicators->second);
                  for ( ; data_val != end_pool; ++data_val, ++dof_indices)
                    for (unsigned int comp=0; comp<n_components; ++comp)
                      operation.process_constraint (*dof_indices, *data_val,
                                                    *src[comp], value[comp]);

                  for (unsigned int comp=0; comp<n_components; ++comp)
                    operation.post_constraints (value[comp],
                                                local_data[comp][ind_local]);

                  ind_local++;
                }

              // get the dof values past the last constraint
              for (; ind_local < n_local_dofs; ++dof_indices, ++ind_local)
                {
                  for (unsigned int comp=0; comp<n_components; ++comp)
                    operation.process_dof (*dof_indices, *src[comp],
                                           local_data[comp][ind_local]);
                }
            }
          else
            {
              // no constraint at all: compiler can unroll at least the
              // vectorization loop
              AssertDimension (dof_info->end_indices(cell)-dof_indices,
                               static_cast<int>(n_local_dofs));
              for (unsigned int j=0; j<n_local_dofs; j+=VectorizedArray<Number>::n_array_elements)
                for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v)
                  for (unsigned int comp=0; comp<n_components; ++comp)
                    operation.process_dof (dof_indices[j+v], *src[comp],
                                           local_data[comp][j+v]);
            }
        }

      // non-standard case: need to fill in zeros for those components that
      // are not present (a bit more expensive), but there is not more than
      // one such cell
      else
        {
          Assert (n_irreg_components_filled > 0, ExcInternalError());
          for ( ; indicators != indicators_end; ++indicators)
            {
              for (unsigned int j=0; j<indicators->first; ++j)
                {
                  // non-constrained case: copy the data from the global
                  // vector, src, to the local one, local_src.
                  for (unsigned int comp=0; comp<n_components; ++comp)
                    operation.process_dof (dof_indices[j], *src[comp],
                                           local_data[comp][ind_local]);

                  // here we jump over all the components that are artificial
                  ++ind_local;
                  while (ind_local % VectorizedArray<Number>::n_array_elements
                         >= n_irreg_components_filled)
                    {
                      for (unsigned int comp=0; comp<n_components; ++comp)
                        operation.process_empty (local_data[comp][ind_local]);
                      ++ind_local;
                    }
                }
              dof_indices += indicators->first;

              // constrained case: build the local value as a linear
              // combination of the global value according to constraint
              Number value [n_components];
              for (unsigned int comp=0; comp<n_components; ++comp)
                operation.pre_constraints (local_data[comp][ind_local],
                                           value[comp]);

              const Number *data_val =
                matrix_info->constraint_pool_begin(indicators->second);
              const Number *end_pool =
                matrix_info->constraint_pool_end(indicators->second);

              for ( ; data_val != end_pool; ++data_val, ++dof_indices)
                for (unsigned int comp=0; comp<n_components; ++comp)
                  operation.process_constraint (*dof_indices, *data_val,
                                                *src[comp], value[comp]);

              for (unsigned int comp=0; comp<n_components; ++comp)
                operation.post_constraints (value[comp],
                                            local_data[comp][ind_local]);
              ind_local++;
              while (ind_local % VectorizedArray<Number>::n_array_elements
                     >= n_irreg_components_filled)
                {
                  for (unsigned int comp=0; comp<n_components; ++comp)
                    operation.process_empty (local_data[comp][ind_local]);
                  ++ind_local;
                }
            }
          for (; ind_local<n_local_dofs; ++dof_indices)
            {
              Assert (dof_indices != dof_info->end_indices(cell),
                      ExcInternalError());

              // non-constrained case: copy the data from the global vector,
              // src, to the local one, local_dst.
              for (unsigned int comp=0; comp<n_components; ++comp)
                operation.process_dof (*dof_indices, *src[comp],
                                       local_data[comp][ind_local]);
              ++ind_local;
              while (ind_local % VectorizedArray<Number>::n_array_elements
                     >= n_irreg_components_filled)
                {
                  for (unsigned int comp=0; comp<n_components; ++comp)
                    operation.process_empty(local_data[comp][ind_local]);
                  ++ind_local;
                }
            }
        }
    }
  else
    // case with vector-valued finite elements where all components are
    // included in one single vector. Assumption: first come all entries to
    // the first component, then all entries to the second one, and so
    // on. This is ensured by the way MatrixFree reads out the indices.
    {
      internal::check_vector_compatibility (*src[0], *dof_info);
      Assert (n_fe_components == n_components_, ExcNotImplemented());
      const unsigned int n_local_dofs =
        dofs_per_cell*VectorizedArray<Number>::n_array_elements * n_components;
      Number   *local_data =
        const_cast<Number *>(&values_dofs[0][0][0]);
      if (at_irregular_cell == false)
        {
          // check whether there is any constraint on the current cell
          if (indicators != indicators_end)
            {
              for ( ; indicators != indicators_end; ++indicators)
                {
                  // run through values up to next constraint
                  for (unsigned int j=0; j<indicators->first; ++j)
                    operation.process_dof (dof_indices[j], *src[0],
                                           local_data[ind_local+j]);
                  ind_local += indicators->first;
                  dof_indices   += indicators->first;

                  // constrained case: build the local value as a linear
                  // combination of the global value according to constraints
                  Number value;
                  operation.pre_constraints (local_data[ind_local], value);

                  const Number *data_val =
                    matrix_info->constraint_pool_begin(indicators->second);
                  const Number *end_pool =
                    matrix_info->constraint_pool_end(indicators->second);

                  for ( ; data_val != end_pool; ++data_val, ++dof_indices)
                    operation.process_constraint (*dof_indices, *data_val,
                                                  *src[0], value);

                  operation.post_constraints (value, local_data[ind_local]);
                  ind_local++;
                }

              // get the dof values past the last constraint
              for (; ind_local<n_local_dofs; ++dof_indices, ++ind_local)
                operation.process_dof (*dof_indices, *src[0],
                                       local_data[ind_local]);
              Assert (dof_indices == dof_info->end_indices(cell),
                      ExcInternalError());
            }
          else
            {
              // no constraint at all: compiler can unroll at least the
              // vectorization loop
              AssertDimension (dof_info->end_indices(cell)-dof_indices,
                               static_cast<int>(n_local_dofs));
              for (unsigned int j=0; j<n_local_dofs; j+=VectorizedArray<Number>::n_array_elements)
                for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v)
                  operation.process_dof (dof_indices[j+v], *src[0],
                                         local_data[j+v]);
            }
        }

      // non-standard case: need to fill in zeros for those components that
      // are not present (a bit more expensive), but there is not more than
      // one such cell
      else
        {
          Assert (n_irreg_components_filled > 0, ExcInternalError());
          for ( ; indicators != indicators_end; ++indicators)
            {
              for (unsigned int j=0; j<indicators->first; ++j)
                {
                  // non-constrained case: copy the data from the global
                  // vector, src, to the local one, local_src.
                  operation.process_dof (dof_indices[j], *src[0],
                                         local_data[ind_local]);

                  // here we jump over all the components that are artificial
                  ++ind_local;
                  while (ind_local % VectorizedArray<Number>::n_array_elements
                         >= n_irreg_components_filled)
                    {
                      operation.process_empty (local_data[ind_local]);
                      ++ind_local;
                    }
                }
              dof_indices += indicators->first;

              // constrained case: build the local value as a linear
              // combination of the global value according to constraint
              Number value;
              operation.pre_constraints (local_data[ind_local], value);

              const Number *data_val =
                matrix_info->constraint_pool_begin(indicators->second);
              const Number *end_pool =
                matrix_info->constraint_pool_end(indicators->second);

              for ( ; data_val != end_pool; ++data_val, ++dof_indices)
                operation.process_constraint (*dof_indices, *data_val,
                                              *src[0], value);

              operation.post_constraints (value, local_data[ind_local]);
              ind_local++;
              while (ind_local % VectorizedArray<Number>::n_array_elements
                     >= n_irreg_components_filled)
                {
                  operation.process_empty (local_data[ind_local]);
                  ++ind_local;
                }
            }
          for (; ind_local<n_local_dofs; ++dof_indices)
            {
              Assert (dof_indices != dof_info->end_indices(cell),
                      ExcInternalError());

              // non-constrained case: copy the data from the global vector,
              // src, to the local one, local_dst.
              operation.process_dof (*dof_indices, *src[0],
                                     local_data[ind_local]);
              ++ind_local;
              while (ind_local % VectorizedArray<Number>::n_array_elements
                     >= n_irreg_components_filled)
                {
                  operation.process_empty (local_data[ind_local]);
                  ++ind_local;
                }
            }
        }
    }
}



template <int dim, int n_components_, typename Number>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,n_components_,Number>
::read_dof_values (const VectorType &src)
{
  // select between block vectors and non-block vectors. Note that the number
  // of components is checked in the internal data
  typename internal::BlockVectorSelector<VectorType,
           IsBlockVector<VectorType>::value>::BaseVectorType *src_data[n_components];
  for (unsigned int d=0; d<n_components; ++d)
    src_data[d] = internal::BlockVectorSelector<VectorType, IsBlockVector<VectorType>::value>::get_vector_component(const_cast<VectorType &>(src), d);

  internal::VectorReader<Number> reader;
  read_write_operation (reader, src_data);

#ifdef DEBUG
  dof_values_initialized = true;
#endif
}



template <int dim, int n_components_, typename Number>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,n_components_,Number>
::read_dof_values (const std::vector<VectorType> &src,
                   const unsigned int             first_index)
{
  AssertIndexRange (first_index, src.size());
  Assert (n_fe_components == 1, ExcNotImplemented());
  Assert ((n_fe_components == 1 ?
           (first_index+n_components <= src.size()) : true),
          ExcIndexRange (first_index + n_components_, 0, src.size()));

  VectorType *src_data [n_components];
  for (unsigned int comp=0; comp<n_components; ++comp)
    src_data[comp] = const_cast<VectorType *>(&src[comp+first_index]);

  internal::VectorReader<Number> reader;
  read_write_operation (reader, src_data);

#ifdef DEBUG
  dof_values_initialized = true;
#endif
}



template <int dim, int n_components_, typename Number>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,n_components_,Number>
::read_dof_values (const std::vector<VectorType *> &src,
                   const unsigned int              first_index)
{
  AssertIndexRange (first_index, src.size());
  Assert (n_fe_components == 1, ExcNotImplemented());
  Assert ((n_fe_components == 1 ?
           (first_index+n_components <= src.size()) : true),
          ExcIndexRange (first_index + n_components_, 0, src.size()));

  const VectorType *src_data [n_components];
  for (unsigned int comp=0; comp<n_components; ++comp)
    src_data[comp] = const_cast<VectorType *>(src[comp+first_index]);

  internal::VectorReader<Number> reader;
  read_write_operation (reader, src_data);

#ifdef DEBUG
  dof_values_initialized = true;
#endif
}



template <int dim, int n_components_, typename Number>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,n_components_,Number>
::read_dof_values_plain (const VectorType &src)
{
  // select between block vectors and non-block vectors. Note that the number
  // of components is checked in the internal data
  const typename internal::BlockVectorSelector<VectorType,
        IsBlockVector<VectorType>::value>::BaseVectorType *src_data[n_components];
  for (unsigned int d=0; d<n_components; ++d)
    src_data[d] = internal::BlockVectorSelector<VectorType, IsBlockVector<VectorType>::value>::get_vector_component(const_cast<VectorType &>(src), d);

  read_dof_values_plain (src_data);
}



template <int dim, int n_components_, typename Number>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,n_components_,Number>
::read_dof_values_plain (const std::vector<VectorType> &src,
                         const unsigned int             first_index)
{
  AssertIndexRange (first_index, src.size());
  Assert (n_fe_components == 1, ExcNotImplemented());
  Assert ((n_fe_components == 1 ?
           (first_index+n_components <= src.size()) : true),
          ExcIndexRange (first_index + n_components_, 0, src.size()));
  const VectorType *src_data [n_components];
  for (unsigned int comp=0; comp<n_components; ++comp)
    src_data[comp] = &src[comp+first_index];
  read_dof_values_plain (src_data);
}



template <int dim, int n_components_, typename Number>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,n_components_,Number>
::read_dof_values_plain (const std::vector<VectorType *> &src,
                         const unsigned int              first_index)
{
  AssertIndexRange (first_index, src.size());
  Assert (n_fe_components == 1, ExcNotImplemented());
  Assert ((n_fe_components == 1 ?
           (first_index+n_components <= src.size()) : true),
          ExcIndexRange (first_index + n_components_, 0, src.size()));
  const VectorType *src_data [n_components];
  for (unsigned int comp=0; comp<n_components; ++comp)
    src_data[comp] = src[comp+first_index];
  read_dof_values_plain (src_data);
}



template <int dim, int n_components_, typename Number>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,n_components_,Number>
::distribute_local_to_global (VectorType &dst) const
{
  Assert (dof_values_initialized==true,
          internal::ExcAccessToUninitializedField());

  // select between block vectors and non-block vectors. Note that the number
  // of components is checked in the internal data
  typename internal::BlockVectorSelector<VectorType,
           IsBlockVector<VectorType>::value>::BaseVectorType *dst_data[n_components];
  for (unsigned int d=0; d<n_components; ++d)
    dst_data[d] = internal::BlockVectorSelector<VectorType, IsBlockVector<VectorType>::value>::get_vector_component(dst, d);

  internal::VectorDistributorLocalToGlobal<Number> distributor;
  read_write_operation (distributor, dst_data);
}



template <int dim, int n_components_, typename Number>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,n_components_,Number>
::distribute_local_to_global (std::vector<VectorType>  &dst,
                              const unsigned int        first_index) const
{
  AssertIndexRange (first_index, dst.size());
  Assert (n_fe_components == 1, ExcNotImplemented());
  Assert ((n_fe_components == 1 ?
           (first_index+n_components <= dst.size()) : true),
          ExcIndexRange (first_index + n_components_, 0, dst.size()));
  Assert (dof_values_initialized==true,
          internal::ExcAccessToUninitializedField());

  VectorType *dst_data [n_components];
  for (unsigned int comp=0; comp<n_components; ++comp)
    dst_data[comp] = &dst[comp+first_index];

  internal::VectorDistributorLocalToGlobal<Number> distributor;
  read_write_operation (distributor, dst_data);
}



template <int dim, int n_components_, typename Number>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,n_components_,Number>
::distribute_local_to_global (std::vector<VectorType *>  &dst,
                              const unsigned int         first_index) const
{
  AssertIndexRange (first_index, dst.size());
  Assert (n_fe_components == 1, ExcNotImplemented());
  Assert ((n_fe_components == 1 ?
           (first_index+n_components <= dst.size()) : true),
          ExcIndexRange (first_index + n_components_, 0, dst.size()));
  Assert (dof_values_initialized==true,
          internal::ExcAccessToUninitializedField());

  VectorType *dst_data [n_components];
  for (unsigned int comp=0; comp<n_components; ++comp)
    dst_data[comp] = dst[comp+first_index];

  internal::VectorDistributorLocalToGlobal<Number> distributor;
  read_write_operation (distributor, dst_data);
}



template <int dim, int n_components_, typename Number>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,n_components_,Number>
::set_dof_values (VectorType &dst) const
{
  Assert (dof_values_initialized==true,
          internal::ExcAccessToUninitializedField());

  // select between block vectors and non-block vectors. Note that the number
  // of components is checked in the internal data
  typename internal::BlockVectorSelector<VectorType,
           IsBlockVector<VectorType>::value>::BaseVectorType *dst_data[n_components];
  for (unsigned int d=0; d<n_components; ++d)
    dst_data[d] = internal::BlockVectorSelector<VectorType, IsBlockVector<VectorType>::value>::get_vector_component(dst, d);

  internal::VectorSetter<Number> setter;
  read_write_operation (setter, dst_data);
}



template <int dim, int n_components_, typename Number>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,n_components_,Number>
::set_dof_values (std::vector<VectorType>  &dst,
                  const unsigned int        first_index) const
{
  AssertIndexRange (first_index, dst.size());
  Assert (n_fe_components == 1, ExcNotImplemented());
  Assert ((n_fe_components == 1 ?
           (first_index+n_components <= dst.size()) : true),
          ExcIndexRange (first_index + n_components_, 0, dst.size()));

  Assert (dof_values_initialized==true,
          internal::ExcAccessToUninitializedField());

  VectorType *dst_data [n_components];
  for (unsigned int comp=0; comp<n_components; ++comp)
    dst_data[comp] = &dst[comp+first_index];

  internal::VectorSetter<Number> setter;
  read_write_operation (setter, dst_data);
}



template <int dim, int n_components_, typename Number>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,n_components_,Number>
::set_dof_values (std::vector<VectorType *>  &dst,
                  const unsigned int         first_index) const
{
  AssertIndexRange (first_index, dst.size());
  Assert (n_fe_components == 1, ExcNotImplemented());
  Assert ((n_fe_components == 1 ?
           (first_index+n_components <= dst.size()) : true),
          ExcIndexRange (first_index + n_components_, 0, dst.size()));

  Assert (dof_values_initialized==true,
          internal::ExcAccessToUninitializedField());

  VectorType *dst_data [n_components];
  for (unsigned int comp=0; comp<n_components; ++comp)
    dst_data[comp] = dst[comp+first_index];

  internal::VectorSetter<Number> setter;
  read_write_operation (setter, dst_data);
}



template <int dim, int n_components_, typename Number>
template<typename VectorType>
inline
void
FEEvaluationBase<dim,n_components_,Number>
::read_dof_values_plain (const VectorType *src[])
{
  // Case without MatrixFree initialization object
  if (matrix_info == 0)
    {
      internal::VectorReader<Number> reader;
      read_write_operation (reader, src);
      return;
    }

  // this is different from the other three operations because we do not use
  // constraints here, so this is a separate function.
  Assert (dof_info != 0, ExcNotInitialized());
  Assert (matrix_info->indices_initialized() == true,
          ExcNotInitialized());
  Assert (cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  Assert (dof_info->store_plain_indices == true, ExcNotInitialized());

  // loop over all local dofs. ind_local holds local number on cell, index
  // iterates over the elements of index_local_to_global and dof_indices
  // points to the global indices stored in index_local_to_global
  const unsigned int *dof_indices = dof_info->begin_indices_plain(cell);
  const unsigned int dofs_per_cell = this->data->dofs_per_cell;

  const unsigned int n_irreg_components_filled = dof_info->row_starts[cell][2];
  const bool at_irregular_cell = n_irreg_components_filled > 0;

  // scalar case (or case when all components have the same degrees of freedom
  // and sit on a different vector each)
  if (n_fe_components == 1)
    {
      const unsigned int n_local_dofs =
        VectorizedArray<Number>::n_array_elements * dofs_per_cell;
      for (unsigned int comp=0; comp<n_components; ++comp)
        internal::check_vector_compatibility (*src[comp], *dof_info);
      Number *local_src_number [n_components];
      for (unsigned int comp=0; comp<n_components; ++comp)
        local_src_number[comp] = &values_dofs[comp][0][0];

      // standard case where there are sufficiently many cells to fill all
      // vectors
      if (at_irregular_cell == false)
        {
          for (unsigned int j=0; j<n_local_dofs; ++j)
            for (unsigned int comp=0; comp<n_components; ++comp)
              local_src_number[comp][j] =
                internal::vector_access (*src[comp], dof_indices[j]);
        }

      // non-standard case: need to fill in zeros for those components that
      // are not present (a bit more expensive), but there is not more than
      // one such cell
      else
        {
          Assert (n_irreg_components_filled > 0, ExcInternalError());
          for (unsigned int ind_local=0; ind_local<n_local_dofs;
               ++dof_indices)
            {
              // non-constrained case: copy the data from the global vector,
              // src, to the local one, local_dst.
              for (unsigned int comp=0; comp<n_components; ++comp)
                local_src_number[comp][ind_local] =
                  internal::vector_access (*src[comp], *dof_indices);
              ++ind_local;
              while (ind_local % VectorizedArray<Number>::n_array_elements >= n_irreg_components_filled)
                {
                  for (unsigned int comp=0; comp<n_components; ++comp)
                    local_src_number[comp][ind_local] = 0.;
                  ++ind_local;
                }
            }
        }
    }
  else
    // case with vector-valued finite elements where all components are
    // included in one single vector. Assumption: first come all entries to
    // the first component, then all entries to the second one, and so
    // on. This is ensured by the way MatrixFree reads out the indices.
    {
      internal::check_vector_compatibility (*src[0], *dof_info);
      Assert (n_fe_components == n_components_, ExcNotImplemented());
      const unsigned int n_local_dofs =
        dofs_per_cell * VectorizedArray<Number>::n_array_elements * n_components;
      Number *local_src_number = &values_dofs[0][0][0];
      if (at_irregular_cell == false)
        {
          for (unsigned int j=0; j<n_local_dofs; ++j)
            local_src_number[j] =
              internal::vector_access (*src[0], dof_indices[j]);
        }

      // non-standard case: need to fill in zeros for those components that
      // are not present (a bit more expensive), but there is not more than
      // one such cell
      else
        {
          Assert (n_irreg_components_filled > 0, ExcInternalError());
          for (unsigned int ind_local=0; ind_local<n_local_dofs; ++dof_indices)
            {
              // non-constrained case: copy the data from the global vector,
              // src, to the local one, local_dst.
              local_src_number[ind_local] =
                internal::vector_access (*src[0], *dof_indices);
              ++ind_local;
              while (ind_local % VectorizedArray<Number>::n_array_elements >= n_irreg_components_filled)
                {
                  local_src_number[ind_local] = 0.;
                  ++ind_local;
                }
            }
        }
    }

#ifdef DEBUG
  dof_values_initialized = true;
#endif
}




/*------------------------------ access to data fields ----------------------*/


template <int dim, int n_components, typename Number>
inline
const VectorizedArray<Number> *
FEEvaluationBase<dim,n_components,Number>::
begin_dof_values () const
{
  return &values_dofs[0][0];
}



template <int dim, int n_components, typename Number>
inline
VectorizedArray<Number> *
FEEvaluationBase<dim,n_components,Number>::
begin_dof_values ()
{
#ifdef DEBUG
  dof_values_initialized = true;
#endif
  return &values_dofs[0][0];
}



template <int dim, int n_components, typename Number>
inline
const VectorizedArray<Number> *
FEEvaluationBase<dim,n_components,Number>::
begin_values () const
{
  Assert (values_quad_initialized || values_quad_submitted,
          ExcNotInitialized());
  return &values_quad[0][0];
}



template <int dim, int n_components, typename Number>
inline
VectorizedArray<Number> *
FEEvaluationBase<dim,n_components,Number>::
begin_values ()
{
#ifdef DEBUG
  values_quad_submitted = true;
#endif
  return &values_quad[0][0];
}



template <int dim, int n_components, typename Number>
inline
const VectorizedArray<Number> *
FEEvaluationBase<dim,n_components,Number>::
begin_gradients () const
{
  Assert (gradients_quad_initialized || gradients_quad_submitted,
          ExcNotInitialized());
  return &gradients_quad[0][0][0];
}



template <int dim, int n_components, typename Number>
inline
VectorizedArray<Number> *
FEEvaluationBase<dim,n_components,Number>::
begin_gradients ()
{
#ifdef DEBUG
  gradients_quad_submitted = true;
#endif
  return &gradients_quad[0][0][0];
}



template <int dim, int n_components, typename Number>
inline
const VectorizedArray<Number> *
FEEvaluationBase<dim,n_components,Number>::
begin_hessians () const
{
  Assert (hessians_quad_initialized, ExcNotInitialized());
  return &hessians_quad[0][0][0];
}



template <int dim, int n_components, typename Number>
inline
VectorizedArray<Number> *
FEEvaluationBase<dim,n_components,Number>::
begin_hessians ()
{
  return &hessians_quad[0][0][0];
}



template <int dim, int n_components_, typename Number>
inline
Tensor<1,n_components_,VectorizedArray<Number> >
FEEvaluationBase<dim,n_components_,Number>
::get_dof_value (const unsigned int dof) const
{
  AssertIndexRange (dof, this->data->dofs_per_cell);
  Tensor<1,n_components_,VectorizedArray<Number> > return_value (false);
  for (unsigned int comp=0; comp<n_components; comp++)
    return_value[comp] = this->values_dofs[comp][dof];
  return return_value;
}



template <int dim, int n_components_, typename Number>
inline
Tensor<1,n_components_,VectorizedArray<Number> >
FEEvaluationBase<dim,n_components_,Number>
::get_value (const unsigned int q_point) const
{
  Assert (this->values_quad_initialized==true,
          internal::ExcAccessToUninitializedField());
  AssertIndexRange (q_point, this->data->n_q_points);
  Tensor<1,n_components_,VectorizedArray<Number> > return_value (false);
  for (unsigned int comp=0; comp<n_components; comp++)
    return_value[comp] = this->values_quad[comp][q_point];
  return return_value;
}



template <int dim, int n_components_, typename Number>
inline
Tensor<1,n_components_,Tensor<1,dim,VectorizedArray<Number> > >
FEEvaluationBase<dim,n_components_,Number>
::get_gradient (const unsigned int q_point) const
{
  Assert (this->gradients_quad_initialized==true,
          internal::ExcAccessToUninitializedField());
  AssertIndexRange (q_point, this->data->n_q_points);

  Tensor<1,n_components_,Tensor<1,dim,VectorizedArray<Number> > > grad_out (false);

  // Cartesian cell
  if (this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      for (unsigned int comp=0; comp<n_components; comp++)
        for (unsigned int d=0; d<dim; ++d)
          grad_out[comp][d] = (this->gradients_quad[comp][d][q_point] *
                               cartesian_data[0][d]);
    }
  // cell with general/affine Jacobian
  else
    {
      const Tensor<2,dim,VectorizedArray<Number> > &jac =
        this->cell_type == internal::MatrixFreeFunctions::general ?
        jacobian[q_point] : jacobian[0];
      for (unsigned int comp=0; comp<n_components; comp++)
        {
          for (unsigned int d=0; d<dim; ++d)
            {
              grad_out[comp][d] = (jac[d][0] *
                                   this->gradients_quad[comp][0][q_point]);
              for (unsigned int e=1; e<dim; ++e)
                grad_out[comp][d] += (jac[d][e] *
                                      this->gradients_quad[comp][e][q_point]);
            }
        }
    }
  return grad_out;
}



namespace internal
{
  // compute tmp = hess_unit(u) * J^T. do this manually because we do not
  // store the lower diagonal because of symmetry
  template <typename Number>
  inline
  void
  hessian_unit_times_jac (const Tensor<2,1,VectorizedArray<Number> > &jac,
                          const VectorizedArray<Number> *const hessians_quad[1],
                          const unsigned int             q_point,
                          VectorizedArray<Number>       (&tmp)[1][1])
  {
    tmp[0][0] = jac[0][0] * hessians_quad[0][q_point];
  }

  template <typename Number>
  inline
  void
  hessian_unit_times_jac (const Tensor<2,2,VectorizedArray<Number> > &jac,
                          const VectorizedArray<Number> *const hessians_quad[3],
                          const unsigned int             q_point,
                          VectorizedArray<Number>       (&tmp)[2][2])
  {
    for (unsigned int d=0; d<2; ++d)
      {
        tmp[0][d] = (jac[d][0] * hessians_quad[0][q_point] +
                     jac[d][1] * hessians_quad[2][q_point]);
        tmp[1][d] = (jac[d][0] * hessians_quad[2][q_point] +
                     jac[d][1] * hessians_quad[1][q_point]);
      }
  }

  template <typename Number>
  inline
  void
  hessian_unit_times_jac (const Tensor<2,3,VectorizedArray<Number> > &jac,
                          const VectorizedArray<Number> *const hessians_quad[6],
                          const unsigned int             q_point,
                          VectorizedArray<Number>       (&tmp)[3][3])
  {
    for (unsigned int d=0; d<3; ++d)
      {
        tmp[0][d] = (jac[d][0] * hessians_quad[0][q_point] +
                     jac[d][1] * hessians_quad[3][q_point] +
                     jac[d][2] * hessians_quad[4][q_point]);
        tmp[1][d] = (jac[d][0] * hessians_quad[3][q_point] +
                     jac[d][1] * hessians_quad[1][q_point] +
                     jac[d][2] * hessians_quad[5][q_point]);
        tmp[2][d] = (jac[d][0] * hessians_quad[4][q_point] +
                     jac[d][1] * hessians_quad[5][q_point] +
                     jac[d][2] * hessians_quad[2][q_point]);
      }
  }
}



template <int dim, int n_components_, typename Number>
inline
Tensor<1,n_components_,Tensor<2,dim,VectorizedArray<Number> > >
FEEvaluationBase<dim,n_components_,Number>
::get_hessian (const unsigned int q_point) const
{
  Assert (this->hessians_quad_initialized==true,
          internal::ExcAccessToUninitializedField());
  AssertIndexRange (q_point, this->data->n_q_points);

  Tensor<2,dim,VectorizedArray<Number> > hessian_out [n_components];

  // Cartesian cell
  if (this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      const Tensor<1,dim,VectorizedArray<Number> > &jac = cartesian_data[0];
      for (unsigned int comp=0; comp<n_components; comp++)
        for (unsigned int d=0; d<dim; ++d)
          {
            hessian_out[comp][d][d] = (this->hessians_quad[comp][d][q_point] *
                                       jac[d] * jac[d]);
            switch (dim)
              {
              case 1:
                break;
              case 2:
                hessian_out[comp][0][1] = (this->hessians_quad[comp][2][q_point] *
                                           jac[0] * jac[1]);
                break;
              case 3:
                hessian_out[comp][0][1] = (this->hessians_quad[comp][3][q_point] *
                                           jac[0] * jac[1]);
                hessian_out[comp][0][2] = (this->hessians_quad[comp][4][q_point] *
                                           jac[0] * jac[2]);
                hessian_out[comp][1][2] = (this->hessians_quad[comp][5][q_point] *
                                           jac[1] * jac[2]);
                break;
              default:
                Assert (false, ExcNotImplemented());
              }
            for (unsigned int e=d+1; e<dim; ++e)
              hessian_out[comp][e][d] = hessian_out[comp][d][e];
          }
    }
  // cell with general Jacobian
  else if (this->cell_type == internal::MatrixFreeFunctions::general)
    {
      Assert (this->mapping_info->second_derivatives_initialized == true,
              ExcNotInitialized());
      const Tensor<2,dim,VectorizedArray<Number> > &jac = jacobian[q_point];
      const Tensor<2,dim,VectorizedArray<Number> > &jac_grad = jacobian_grad[q_point];
      const Tensor<1,(dim>1?dim*(dim-1)/2:1),
            Tensor<1,dim,VectorizedArray<Number> > >
            & jac_grad_UT = jacobian_grad_upper[q_point];
      for (unsigned int comp=0; comp<n_components; comp++)
        {
          // compute laplacian before the gradient because it needs to access
          // unscaled gradient data
          VectorizedArray<Number> tmp[dim][dim];
          internal::hessian_unit_times_jac (jac, this->hessians_quad[comp],
                                            q_point, tmp);

          // compute first part of hessian, J * tmp = J * hess_unit(u) * J^T
          for (unsigned int d=0; d<dim; ++d)
            for (unsigned int e=d; e<dim; ++e)
              {
                hessian_out[comp][d][e] = jac[d][0] * tmp[0][e];
                for (unsigned int f=1; f<dim; ++f)
                  hessian_out[comp][d][e] += jac[d][f] * tmp[f][e];
              }

          // add diagonal part of J' * grad(u)
          for (unsigned int d=0; d<dim; ++d)
            for (unsigned int e=0; e<dim; ++e)
              hessian_out[comp][d][d] += (jac_grad[d][e] *
                                          this->gradients_quad[comp][e][q_point]);

          // add off-diagonal part of J' * grad(u)
          for (unsigned int d=0, count=0; d<dim; ++d)
            for (unsigned int e=d+1; e<dim; ++e, ++count)
              for (unsigned int f=0; f<dim; ++f)
                hessian_out[comp][d][e] += (jac_grad_UT[count][f] *
                                            this->gradients_quad[comp][f][q_point]);

          // take symmetric part
          for (unsigned int d=0; d<dim; ++d)
            for (unsigned int e=d+1; e<dim; ++e)
              hessian_out[comp][e][d] = hessian_out[comp][d][e];
        }
    }
  // cell with general Jacobian, but constant within the cell
  else // if (this->cell_type == internal::MatrixFreeFunctions::affine)
    {
      const Tensor<2,dim,VectorizedArray<Number> > &jac = jacobian[0];
      for (unsigned int comp=0; comp<n_components; comp++)
        {
          // compute laplacian before the gradient because it needs to access
          // unscaled gradient data
          VectorizedArray<Number> tmp[dim][dim];
          internal::hessian_unit_times_jac (jac, this->hessians_quad[comp],
                                            q_point, tmp);

          // compute first part of hessian, J * tmp = J * hess_unit(u) * J^T
          for (unsigned int d=0; d<dim; ++d)
            for (unsigned int e=d; e<dim; ++e)
              {
                hessian_out[comp][d][e] = jac[d][0] * tmp[0][e];
                for (unsigned int f=1; f<dim; ++f)
                  hessian_out[comp][d][e] += jac[d][f] * tmp[f][e];
              }

          // no J' * grad(u) part here because the Jacobian is constant
          // throughout the cell and hence, its derivative is zero

          // take symmetric part
          for (unsigned int d=0; d<dim; ++d)
            for (unsigned int e=d+1; e<dim; ++e)
              hessian_out[comp][e][d] = hessian_out[comp][d][e];
        }
    }
  return Tensor<1,n_components_,Tensor<2,dim,VectorizedArray<Number> > >(hessian_out);
}



template <int dim, int n_components_, typename Number>
inline
Tensor<1,n_components_,Tensor<1,dim,VectorizedArray<Number> > >
FEEvaluationBase<dim,n_components_,Number>
::get_hessian_diagonal (const unsigned int q_point) const
{
  Assert (this->hessians_quad_initialized==true,
          internal::ExcAccessToUninitializedField());
  AssertIndexRange (q_point, this->data->n_q_points);

  Tensor<1,n_components_,Tensor<1,dim,VectorizedArray<Number> > > hessian_out (false);

  // Cartesian cell
  if (this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      const Tensor<1,dim,VectorizedArray<Number> > &jac = cartesian_data[0];
      for (unsigned int comp=0; comp<n_components; comp++)
        for (unsigned int d=0; d<dim; ++d)
          hessian_out[comp][d] = (this->hessians_quad[comp][d][q_point] *
                                  jac[d] * jac[d]);
    }
  // cell with general Jacobian
  else if (this->cell_type == internal::MatrixFreeFunctions::general)
    {
      Assert (this->mapping_info->second_derivatives_initialized == true,
              ExcNotInitialized());
      const Tensor<2,dim,VectorizedArray<Number> > &jac = jacobian[q_point];
      const Tensor<2,dim,VectorizedArray<Number> > &jac_grad = jacobian_grad[q_point];
      for (unsigned int comp=0; comp<n_components; comp++)
        {
          // compute laplacian before the gradient because it needs to access
          // unscaled gradient data
          VectorizedArray<Number> tmp[dim][dim];
          internal::hessian_unit_times_jac (jac, this->hessians_quad[comp],
                                            q_point, tmp);

          // compute only the trace part of hessian, J * tmp = J *
          // hess_unit(u) * J^T
          for (unsigned int d=0; d<dim; ++d)
            {
              hessian_out[comp][d] = jac[d][0] * tmp[0][d];
              for (unsigned int f=1; f<dim; ++f)
                hessian_out[comp][d] += jac[d][f] * tmp[f][d];
            }

          for (unsigned int d=0; d<dim; ++d)
            for (unsigned int e=0; e<dim; ++e)
              hessian_out[comp][d] += (jac_grad[d][e] *
                                       this->gradients_quad[comp][e][q_point]);
        }
    }
  // cell with general Jacobian, but constant within the cell
  else // if (this->cell_type == internal::MatrixFreeFunctions::affine)
    {
      const Tensor<2,dim,VectorizedArray<Number> > &jac = jacobian[0];
      for (unsigned int comp=0; comp<n_components; comp++)
        {
          // compute laplacian before the gradient because it needs to access
          // unscaled gradient data
          VectorizedArray<Number> tmp[dim][dim];
          internal::hessian_unit_times_jac (jac, this->hessians_quad[comp],
                                            q_point, tmp);

          // compute only the trace part of hessian, J * tmp = J *
          // hess_unit(u) * J^T
          for (unsigned int d=0; d<dim; ++d)
            {
              hessian_out[comp][d] = jac[d][0] * tmp[0][d];
              for (unsigned int f=1; f<dim; ++f)
                hessian_out[comp][d] += jac[d][f] * tmp[f][d];
            }
        }
    }
  return hessian_out;
}



template <int dim, int n_components_, typename Number>
inline
Tensor<1,n_components_,VectorizedArray<Number> >
FEEvaluationBase<dim,n_components_,Number>
::get_laplacian (const unsigned int q_point) const
{
  Assert (this->hessians_quad_initialized==true,
          internal::ExcAccessToUninitializedField());
  AssertIndexRange (q_point, this->data->n_q_points);
  Tensor<1,n_components_,VectorizedArray<Number> > laplacian_out (false);
  const Tensor<1,n_components_,Tensor<1,dim,VectorizedArray<Number> > > hess_diag
    = get_hessian_diagonal(q_point);
  for (unsigned int comp=0; comp<n_components; ++comp)
    {
      laplacian_out[comp] = hess_diag[comp][0];
      for (unsigned int d=1; d<dim; ++d)
        laplacian_out[comp] += hess_diag[comp][d];
    }
  return laplacian_out;
}



template <int dim, int n_components_, typename Number>
inline
void
FEEvaluationBase<dim,n_components_,Number>
::submit_dof_value (const Tensor<1,n_components_,VectorizedArray<Number> > val_in,
                    const unsigned int dof)
{
#ifdef DEBUG
  this->dof_values_initialized = true;
#endif
  AssertIndexRange (dof, this->data->dofs_per_cell);
  for (unsigned int comp=0; comp<n_components; comp++)
    this->values_dofs[comp][dof] = val_in[comp];
}



template <int dim, int n_components_, typename Number>
inline
void
FEEvaluationBase<dim,n_components_,Number>
::submit_value (const Tensor<1,n_components_,VectorizedArray<Number> > val_in,
                const unsigned int q_point)
{
#ifdef DEBUG
  Assert (this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  AssertIndexRange (q_point, this->data->n_q_points);
  this->values_quad_submitted = true;
#endif
  if (this->cell_type == internal::MatrixFreeFunctions::general)
    {
      const VectorizedArray<Number> JxW = J_value[q_point];
      for (unsigned int comp=0; comp<n_components; ++comp)
        this->values_quad[comp][q_point] = val_in[comp] * JxW;
    }
  else //if (this->cell_type < internal::MatrixFreeFunctions::general)
    {
      const VectorizedArray<Number> JxW = J_value[0] * quadrature_weights[q_point];
      for (unsigned int comp=0; comp<n_components; ++comp)
        this->values_quad[comp][q_point] = val_in[comp] * JxW;
    }
}



template <int dim, int n_components_, typename Number>
inline
void
FEEvaluationBase<dim,n_components_,Number>
::submit_gradient (const Tensor<1,n_components_,
                   Tensor<1,dim,VectorizedArray<Number> > >grad_in,
                   const unsigned int q_point)
{
#ifdef DEBUG
  Assert (this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  AssertIndexRange (q_point, this->data->n_q_points);
  this->gradients_quad_submitted = true;
#endif
  if (this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      const VectorizedArray<Number> JxW = J_value[0] * quadrature_weights[q_point];
      for (unsigned int comp=0; comp<n_components; comp++)
        for (unsigned int d=0; d<dim; ++d)
          this->gradients_quad[comp][d][q_point] = (grad_in[comp][d] *
                                                    cartesian_data[0][d] * JxW);
    }
  else
    {
      const Tensor<2,dim,VectorizedArray<Number> > &jac =
        this->cell_type == internal::MatrixFreeFunctions::general ?
        jacobian[q_point] : jacobian[0];
      const VectorizedArray<Number> JxW =
        this->cell_type == internal::MatrixFreeFunctions::general ?
        J_value[q_point] : J_value[0] * quadrature_weights[q_point];
      for (unsigned int comp=0; comp<n_components; ++comp)
        for (unsigned int d=0; d<dim; ++d)
          {
            VectorizedArray<Number> new_val = jac[0][d] * grad_in[comp][0];
            for (unsigned int e=1; e<dim; ++e)
              new_val += (jac[e][d] * grad_in[comp][e]);
            this->gradients_quad[comp][d][q_point] = new_val * JxW;
          }
    }
}



template <int dim, int n_components_, typename Number>
inline
Tensor<1,n_components_,VectorizedArray<Number> >
FEEvaluationBase<dim,n_components_,Number>
::integrate_value () const
{
#ifdef DEBUG
  Assert (this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  Assert (this->values_quad_submitted == true,
          internal::ExcAccessToUninitializedField());
#endif
  Tensor<1,n_components_,VectorizedArray<Number> > return_value (false);
  for (unsigned int comp=0; comp<n_components; ++comp)
    return_value[comp] = this->values_quad[comp][0];
  const unsigned int n_q_points = this->data->n_q_points;
  for (unsigned int q=1; q<n_q_points; ++q)
    for (unsigned int comp=0; comp<n_components; ++comp)
      return_value[comp] += this->values_quad[comp][q];
  return (return_value);
}



/*----------------------- FEEvaluationAccess --------------------------------*/


template <int dim, int n_components_, typename Number>
inline
FEEvaluationAccess<dim,n_components_,Number>
::FEEvaluationAccess (const MatrixFree<dim,Number> &data_in,
                      const unsigned int fe_no,
                      const unsigned int quad_no_in,
                      const unsigned int fe_degree,
                      const unsigned int n_q_points)
  :
  FEEvaluationBase <dim,n_components_,Number>
  (data_in, fe_no, quad_no_in, fe_degree, n_q_points)
{}



template <int dim, int n_components_, typename Number>
inline
FEEvaluationAccess<dim,n_components_,Number>
::FEEvaluationAccess (const MappingFEEvaluation<dim,Number> &geometry,
                      const DoFHandler<dim>               &dof_handler,
                      const unsigned int                   first_selected_component)
  :
  FEEvaluationBase <dim,n_components_,Number> (geometry, dof_handler, first_selected_component)
{}



template <int dim, int n_components_, typename Number>
inline
FEEvaluationAccess<dim,n_components_,Number>
::FEEvaluationAccess (const FEEvaluationAccess<dim,n_components_,Number> &other)
  :
  FEEvaluationBase <dim,n_components_,Number>(other)
{}




/*-------------------- FEEvaluationAccess scalar ----------------------------*/


template <int dim, typename Number>
inline
FEEvaluationAccess<dim,1,Number>
::FEEvaluationAccess (const MatrixFree<dim,Number> &data_in,
                      const unsigned int fe_no,
                      const unsigned int quad_no_in,
                      const unsigned int fe_degree,
                      const unsigned int n_q_points)
  :
  FEEvaluationBase <dim,1,Number>
  (data_in, fe_no, quad_no_in, fe_degree, n_q_points)
{}



template <int dim, typename Number>
inline
FEEvaluationAccess<dim,1,Number>
::FEEvaluationAccess (const MappingFEEvaluation<dim,Number> &geometry,
                      const DoFHandler<dim>               &dof_handler,
                      const unsigned int                   first_selected_component)
  :
  FEEvaluationBase <dim,1,Number> (geometry, dof_handler, first_selected_component)
{}



template <int dim, typename Number>
inline
FEEvaluationAccess<dim,1,Number>
::FEEvaluationAccess (const FEEvaluationAccess<dim,1,Number> &other)
  :
  FEEvaluationBase <dim,1,Number>(other)
{}



template <int dim, typename Number>
inline
VectorizedArray<Number>
FEEvaluationAccess<dim,1,Number>
::get_dof_value (const unsigned int dof) const
{
  AssertIndexRange (dof, this->data->dofs_per_cell);
  return this->values_dofs[0][dof];
}



template <int dim, typename Number>
inline
VectorizedArray<Number>
FEEvaluationAccess<dim,1,Number>
::get_value (const unsigned int q_point) const
{
  Assert (this->values_quad_initialized==true,
          internal::ExcAccessToUninitializedField());
  AssertIndexRange (q_point, this->data->n_q_points);
  return this->values_quad[0][q_point];
}



template <int dim, typename Number>
inline
Tensor<1,dim,VectorizedArray<Number> >
FEEvaluationAccess<dim,1,Number>
::get_gradient (const unsigned int q_point) const
{
  // could use the base class gradient, but that involves too many inefficient
  // initialization operations on tensors

  Assert (this->gradients_quad_initialized==true,
          internal::ExcAccessToUninitializedField());
  AssertIndexRange (q_point, this->data->n_q_points);

  Tensor<1,dim,VectorizedArray<Number> > grad_out (false);

  // Cartesian cell
  if (this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      for (unsigned int d=0; d<dim; ++d)
        grad_out[d] = (this->gradients_quad[0][d][q_point] *
                       this->cartesian_data[0][d]);
    }
  // cell with general/constant Jacobian
  else
    {
      const Tensor<2,dim,VectorizedArray<Number> > &jac =
        this->cell_type == internal::MatrixFreeFunctions::general ?
        this->jacobian[q_point] : this->jacobian[0];
      for (unsigned int d=0; d<dim; ++d)
        {
          grad_out[d] = (jac[d][0] * this->gradients_quad[0][0][q_point]);
          for (unsigned int e=1; e<dim; ++e)
            grad_out[d] += (jac[d][e] * this->gradients_quad[0][e][q_point]);
        }
    }
  return grad_out;
}



template <int dim, typename Number>
inline
Tensor<2,dim,VectorizedArray<Number> >
FEEvaluationAccess<dim,1,Number>
::get_hessian (const unsigned int q_point) const
{
  return BaseClass::get_hessian(q_point)[0];
}



template <int dim, typename Number>
inline
Tensor<1,dim,VectorizedArray<Number> >
FEEvaluationAccess<dim,1,Number>
::get_hessian_diagonal (const unsigned int q_point) const
{
  return BaseClass::get_hessian_diagonal(q_point)[0];
}



template <int dim, typename Number>
inline
VectorizedArray<Number>
FEEvaluationAccess<dim,1,Number>
::get_laplacian (const unsigned int q_point) const
{
  return BaseClass::get_laplacian(q_point)[0];
}



template <int dim, typename Number>
inline
void
FEEvaluationAccess<dim,1,Number>
::submit_dof_value (const VectorizedArray<Number> val_in,
                    const unsigned int dof)
{
#ifdef DEBUG
  this->dof_values_initialized = true;
  AssertIndexRange (dof, this->data->dofs_per_cell);
#endif
  this->values_dofs[0][dof] = val_in;
}



template <int dim, typename Number>
inline
void
FEEvaluationAccess<dim,1,Number>
::submit_value (const VectorizedArray<Number> val_in,
                const unsigned int q_point)
{
#ifdef DEBUG
  Assert (this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  AssertIndexRange (q_point, this->data->n_q_points);
  this->values_quad_submitted = true;
#endif
  if (this->cell_type == internal::MatrixFreeFunctions::general)
    {
      const VectorizedArray<Number> JxW = this->J_value[q_point];
      this->values_quad[0][q_point] = val_in * JxW;
    }
  else //if (this->cell_type < internal::MatrixFreeFunctions::general)
    {
      const VectorizedArray<Number> JxW = this->J_value[0] * this->quadrature_weights[q_point];
      this->values_quad[0][q_point] = val_in * JxW;
    }
}



template <int dim, typename Number>
inline
void
FEEvaluationAccess<dim,1,Number>
::submit_gradient (const Tensor<1,dim,VectorizedArray<Number> > grad_in,
                   const unsigned int q_point)
{
#ifdef DEBUG
  Assert (this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  AssertIndexRange (q_point, this->data->n_q_points);
  this->gradients_quad_submitted = true;
#endif
  if (this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      const VectorizedArray<Number> JxW = this->J_value[0] * this->quadrature_weights[q_point];
      for (unsigned int d=0; d<dim; ++d)
        this->gradients_quad[0][d][q_point] = (grad_in[d] *
                                               this->cartesian_data[0][d] *
                                               JxW);
    }
  // general/affine cell type
  else
    {
      const Tensor<2,dim,VectorizedArray<Number> > &jac =
        this->cell_type == internal::MatrixFreeFunctions::general ?
        this->jacobian[q_point] : this->jacobian[0];
      const VectorizedArray<Number> JxW =
        this->cell_type == internal::MatrixFreeFunctions::general ?
        this->J_value[q_point] : this->J_value[0] * this->quadrature_weights[q_point];
      for (unsigned int d=0; d<dim; ++d)
        {
          VectorizedArray<Number> new_val = jac[0][d] * grad_in[0];
          for (unsigned int e=1; e<dim; ++e)
            new_val += jac[e][d] * grad_in[e];
          this->gradients_quad[0][d][q_point] = new_val * JxW;
        }
    }
}



template <int dim, typename Number>
inline
VectorizedArray<Number>
FEEvaluationAccess<dim,1,Number>
::integrate_value () const
{
  return BaseClass::integrate_value()[0];
}




/*----------------- FEEvaluationAccess vector-valued ------------------------*/


template <int dim, typename Number>
inline
FEEvaluationAccess<dim,dim,Number>
::FEEvaluationAccess (const MatrixFree<dim,Number> &data_in,
                      const unsigned int fe_no,
                      const unsigned int quad_no_in,
                      const unsigned int fe_degree,
                      const unsigned int n_q_points)
  :
  FEEvaluationBase <dim,dim,Number>
  (data_in, fe_no, quad_no_in, fe_degree, n_q_points)
{}



template <int dim, typename Number>
inline
FEEvaluationAccess<dim,dim,Number>
::FEEvaluationAccess (const MappingFEEvaluation<dim,Number> &geometry,
                      const DoFHandler<dim>               &dof_handler,
                      const unsigned int                   first_selected_component)
  :
  FEEvaluationBase <dim,dim,Number> (geometry, dof_handler, first_selected_component)
{}



template <int dim, typename Number>
inline
FEEvaluationAccess<dim,dim,Number>
::FEEvaluationAccess (const FEEvaluationAccess<dim,dim,Number> &other)
  :
  FEEvaluationBase <dim,dim,Number>(other)
{}



template <int dim, typename Number>
inline
Tensor<2,dim,VectorizedArray<Number> >
FEEvaluationAccess<dim,dim,Number>
::get_gradient (const unsigned int q_point) const
{
  return BaseClass::get_gradient (q_point);
}



template <int dim, typename Number>
inline
VectorizedArray<Number>
FEEvaluationAccess<dim,dim,Number>
::get_divergence (const unsigned int q_point) const
{
  Assert (this->gradients_quad_initialized==true,
          internal::ExcAccessToUninitializedField());
  AssertIndexRange (q_point, this->data->n_q_points);

  VectorizedArray<Number> divergence;

  // Cartesian cell
  if (this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      divergence = (this->gradients_quad[0][0][q_point] *
                    this->cartesian_data[0][0]);
      for (unsigned int d=1; d<dim; ++d)
        divergence += (this->gradients_quad[d][d][q_point] *
                       this->cartesian_data[0][d]);
    }
  // cell with general/constant Jacobian
  else
    {
      const Tensor<2,dim,VectorizedArray<Number> > &jac =
        this->cell_type == internal::MatrixFreeFunctions::general ?
        this->jacobian[q_point] : this->jacobian[0];
      divergence = (jac[0][0] * this->gradients_quad[0][0][q_point]);
      for (unsigned int e=1; e<dim; ++e)
        divergence += (jac[0][e] * this->gradients_quad[0][e][q_point]);
      for (unsigned int d=1; d<dim; ++d)
        for (unsigned int e=0; e<dim; ++e)
          divergence += (jac[d][e] * this->gradients_quad[d][e][q_point]);
    }
  return divergence;
}



template <int dim, typename Number>
inline
SymmetricTensor<2,dim,VectorizedArray<Number> >
FEEvaluationAccess<dim,dim,Number>
::get_symmetric_gradient (const unsigned int q_point) const
{
  // copy from generic function into dim-specialization function
  const Tensor<2,dim,VectorizedArray<Number> > grad = get_gradient(q_point);
  VectorizedArray<Number> symmetrized [(dim*dim+dim)/2];
  VectorizedArray<Number> half = make_vectorized_array (0.5);
  for (unsigned int d=0; d<dim; ++d)
    symmetrized[d] = grad[d][d];
  switch (dim)
    {
    case 1:
      break;
    case 2:
      symmetrized[2] = grad[0][1] + grad[1][0];
      symmetrized[2] *= half;
      break;
    case 3:
      symmetrized[3] = grad[0][1] + grad[1][0];
      symmetrized[3] *= half;
      symmetrized[4] = grad[0][2] + grad[2][0];
      symmetrized[4] *= half;
      symmetrized[5] = grad[1][2] + grad[2][1];
      symmetrized[5] *= half;
      break;
    default:
      Assert (false, ExcNotImplemented());
    }
  return SymmetricTensor<2,dim,VectorizedArray<Number> > (symmetrized);
}



template <int dim, typename Number>
inline
Tensor<1,dim==2?1:dim,VectorizedArray<Number> >
FEEvaluationAccess<dim,dim,Number>
::get_curl (const unsigned int q_point) const
{
  // copy from generic function into dim-specialization function
  const Tensor<2,dim,VectorizedArray<Number> > grad = get_gradient(q_point);
  Tensor<1,dim==2?1:dim,VectorizedArray<Number> > curl (false);
  switch (dim)
    {
    case 1:
      Assert (false,
              ExcMessage("Computing the curl in 1d is not a useful operation"));
      break;
    case 2:
      curl[0] = grad[1][0] - grad[0][1];
      break;
    case 3:
      curl[0] = grad[2][1] - grad[1][2];
      curl[1] = grad[0][2] - grad[2][0];
      curl[2] = grad[1][0] - grad[0][1];
      break;
    default:
      Assert (false, ExcNotImplemented());
    }
  return curl;
}



template <int dim, typename Number>
inline
Tensor<2,dim,VectorizedArray<Number> >
FEEvaluationAccess<dim,dim,Number>
::get_hessian_diagonal (const unsigned int q_point) const
{
  Assert (this->hessians_quad_initialized==true,
          internal::ExcAccessToUninitializedField());
  AssertIndexRange (q_point, this->data->n_q_points);

  return BaseClass::get_hessian_diagonal (q_point);
}



template <int dim, typename Number>
inline
Tensor<3,dim,VectorizedArray<Number> >
FEEvaluationAccess<dim,dim,Number>
::get_hessian (const unsigned int q_point) const
{
  Assert (this->hessians_quad_initialized==true,
          internal::ExcAccessToUninitializedField());
  AssertIndexRange (q_point, this->data->n_q_points);
  return BaseClass::get_hessian(q_point);
}



template <int dim, typename Number>
inline
void
FEEvaluationAccess<dim,dim,Number>
::submit_gradient (const Tensor<2,dim,VectorizedArray<Number> > grad_in,
                   const unsigned int q_point)
{
  BaseClass::submit_gradient (grad_in, q_point);
}



template <int dim, typename Number>
inline
void
FEEvaluationAccess<dim,dim,Number>
::submit_gradient (const Tensor<1,dim,Tensor<1,dim,VectorizedArray<Number> > >
                   grad_in,
                   const unsigned int q_point)
{
  BaseClass::submit_gradient(grad_in, q_point);
}



template <int dim, typename Number>
inline
void
FEEvaluationAccess<dim,dim,Number>
::submit_divergence (const VectorizedArray<Number> div_in,
                     const unsigned int q_point)
{
#ifdef DEBUG
  Assert (this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  AssertIndexRange (q_point, this->data->n_q_points);
  this->gradients_quad_submitted = true;
#endif
  if (this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      const VectorizedArray<Number> fac = this->J_value[0] *
                                          this->quadrature_weights[q_point] * div_in;
      for (unsigned int d=0; d<dim; ++d)
        {
          this->gradients_quad[d][d][q_point] = (fac *
                                                 this->cartesian_data[0][d]);
          for (unsigned int e=d+1; e<dim; ++e)
            {
              this->gradients_quad[d][e][q_point] = VectorizedArray<Number>();
              this->gradients_quad[e][d][q_point] = VectorizedArray<Number>();
            }
        }
    }
  else
    {
      const Tensor<2,dim,VectorizedArray<Number> > &jac =
        this->cell_type == internal::MatrixFreeFunctions::general ?
        this->jacobian[q_point] : this->jacobian[0];
      const VectorizedArray<Number> fac =
        (this->cell_type == internal::MatrixFreeFunctions::general ?
         this->J_value[q_point] : this->J_value[0] *
         this->quadrature_weights[q_point]) * div_in;
      for (unsigned int d=0; d<dim; ++d)
        {
          for (unsigned int e=0; e<dim; ++e)
            this->gradients_quad[d][e][q_point] = jac[d][e] * fac;
        }
    }
}



template <int dim, typename Number>
inline
void
FEEvaluationAccess<dim,dim,Number>
::submit_symmetric_gradient(const SymmetricTensor<2,dim,VectorizedArray<Number> >
                            sym_grad,
                            const unsigned int q_point)
{
  // could have used base class operator, but that involves some overhead
  // which is inefficient. it is nice to have the symmetric tensor because
  // that saves some operations
#ifdef DEBUG
  Assert (this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  AssertIndexRange (q_point, this->data->n_q_points);
  this->gradients_quad_submitted = true;
#endif
  if (this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      const VectorizedArray<Number> JxW = this->J_value[0] * this->quadrature_weights[q_point];
      for (unsigned int d=0; d<dim; ++d)
        this->gradients_quad[d][d][q_point] = (sym_grad.access_raw_entry(d) *
                                               JxW *
                                               this->cartesian_data[0][d]);
      for (unsigned int e=0, counter=dim; e<dim; ++e)
        for (unsigned int d=e+1; d<dim; ++d, ++counter)
          {
            const VectorizedArray<Number> value = sym_grad.access_raw_entry(counter) * JxW;
            this->gradients_quad[e][d][q_point] = (value *
                                                   this->cartesian_data[0][d]);
            this->gradients_quad[d][e][q_point] = (value *
                                                   this->cartesian_data[0][e]);
          }
    }
  // general/affine cell type
  else
    {
      const VectorizedArray<Number> JxW =
        this->cell_type == internal::MatrixFreeFunctions::general ?
        this->J_value[q_point] : this->J_value[0] * this->quadrature_weights[q_point];
      const Tensor<2,dim,VectorizedArray<Number> > &jac =
        this->cell_type == internal::MatrixFreeFunctions::general ?
        this->jacobian[q_point] : this->jacobian[0];
      VectorizedArray<Number> weighted [dim][dim];
      for (unsigned int i=0; i<dim; ++i)
        weighted[i][i] = sym_grad.access_raw_entry(i) * JxW;
      for (unsigned int i=0, counter=dim; i<dim; ++i)
        for (unsigned int j=i+1; j<dim; ++j, ++counter)
          {
            const VectorizedArray<Number> value = sym_grad.access_raw_entry(counter) * JxW;
            weighted[i][j] = value;
            weighted[j][i] = value;
          }
      for (unsigned int comp=0; comp<dim; ++comp)
        for (unsigned int d=0; d<dim; ++d)
          {
            VectorizedArray<Number> new_val = jac[0][d] * weighted[comp][0];
            for (unsigned int e=1; e<dim; ++e)
              new_val += jac[e][d] * weighted[comp][e];
            this->gradients_quad[comp][d][q_point] = new_val;
          }
    }
}



template <int dim, typename Number>
inline
void
FEEvaluationAccess<dim,dim,Number>
::submit_curl (const Tensor<1,dim==2?1:dim,VectorizedArray<Number> > curl,
               const unsigned int q_point)
{
  Tensor<2,dim,VectorizedArray<Number> > grad;
  switch (dim)
    {
    case 1:
      Assert (false,
              ExcMessage("Testing by the curl in 1d is not a useful operation"));
      break;
    case 2:
      grad[1][0] = curl[0];
      grad[0][1] = -curl[0];
      break;
    case 3:
      grad[2][1] = curl[0];
      grad[1][2] = -curl[0];
      grad[0][2] = curl[1];
      grad[2][0] = -curl[1];
      grad[1][0] = curl[2];
      grad[0][1] = -curl[2];
      break;
    default:
      Assert (false, ExcNotImplemented());
    }
  submit_gradient (grad, q_point);
}



namespace internal
{
  /**
   * In this namespace, the evaluator routines that evaluate the tensor
   * products are implemented.
   */
  enum EvaluatorVariant
  {
    evaluate_general,
    evaluate_symmetric,
    evaluate_evenodd
  };

  /**
   * Generic evaluator framework
   */
  template <EvaluatorVariant variant, int dim, int fe_degree, int n_q_points_1d,
            typename Number>
  struct EvaluatorTensorProduct
  {};

  /**
   * Internal evaluator for 1d-3d shape function using the tensor product form
   * of the basis functions
   */
  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  struct EvaluatorTensorProduct<evaluate_general,dim,fe_degree,n_q_points_1d,Number>
  {
    static const unsigned int dofs_per_cell = Utilities::fixed_int_power<fe_degree+1,dim>::value;
    static const unsigned int n_q_points = Utilities::fixed_int_power<n_q_points_1d,dim>::value;

    /**
     * Empty constructor. Does nothing. Be careful when using 'values' and
     * related methods because they need to be filled with the other pointer
     */
    EvaluatorTensorProduct ()
      :
      shape_values (0),
      shape_gradients (0),
      shape_hessians (0)
    {}

    /**
     * Constructor, taking the data from ShapeInfo
     */
    EvaluatorTensorProduct (const AlignedVector<Number> &shape_values,
                            const AlignedVector<Number> &shape_gradients,
                            const AlignedVector<Number> &shape_hessians)
      :
      shape_values (shape_values.begin()),
      shape_gradients (shape_gradients.begin()),
      shape_hessians (shape_hessians.begin())
    {}

    template <int direction, bool dof_to_quad, bool add>
    void
    values (const Number in [],
            Number       out[]) const
    {
      apply<direction,dof_to_quad,add>(shape_values, in, out);
    }

    template <int direction, bool dof_to_quad, bool add>
    void
    gradients (const Number in [],
               Number       out[]) const
    {
      apply<direction,dof_to_quad,add>(shape_gradients, in, out);
    }

    template <int direction, bool dof_to_quad, bool add>
    void
    hessians (const Number in [],
              Number       out[]) const
    {
      apply<direction,dof_to_quad,add>(shape_hessians, in, out);
    }

    template <int direction, bool dof_to_quad, bool add>
    static void apply (const Number *shape_data,
                       const Number in [],
                       Number       out []);

    const Number *shape_values;
    const Number *shape_gradients;
    const Number *shape_hessians;
  };

  // evaluates the given shape data in 1d-3d using the tensor product
  // form. does not use a particular layout of entries in the matrices
  // like the functions below and corresponds to a usual matrix-matrix
  // product
  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  template <int direction, bool dof_to_quad, bool add>
  inline
  void
  EvaluatorTensorProduct<evaluate_general,dim,fe_degree,n_q_points_1d,Number>
  ::apply (const Number *shape_data,
           const Number in [],
           Number       out [])
  {
    AssertIndexRange (direction, dim);
    const int mm     = dof_to_quad ? (fe_degree+1) : n_q_points_1d,
              nn     = dof_to_quad ? n_q_points_1d : (fe_degree+1);

    const int n_blocks1 = (dim > 1 ? (direction > 0 ? nn : mm) : 1);
    const int n_blocks2 = (dim > 2 ? (direction > 1 ? nn : mm) : 1);
    const int stride    = Utilities::fixed_int_power<nn,direction>::value;

    for (int i2=0; i2<n_blocks2; ++i2)
      {
        for (int i1=0; i1<n_blocks1; ++i1)
          {
            for (int col=0; col<nn; ++col)
              {
                Number val0;
                if (dof_to_quad == true)
                  val0 = shape_data[col];
                else
                  val0 = shape_data[col*n_q_points_1d];
                Number res0 = val0 * in[0];
                for (int ind=1; ind<mm; ++ind)
                  {
                    if (dof_to_quad == true)
                      val0 = shape_data[ind*n_q_points_1d+col];
                    else
                      val0 = shape_data[col*n_q_points_1d+ind];
                    res0 += val0 * in[stride*ind];
                  }
                if (add == false)
                  out[stride*col]  = res0;
                else
                  out[stride*col] += res0;
              }

            // increment: in regular case, just go to the next point in
            // x-direction. If we are at the end of one chunk in x-dir, need
            // to jump over to the next layer in z-direction
            switch (direction)
              {
              case 0:
                in += mm;
                out += nn;
                break;
              case 1:
              case 2:
                ++in;
                ++out;
                break;
              default:
                Assert (false, ExcNotImplemented());
              }
          }
        if (direction == 1)
          {
            in += nn*(mm-1);
            out += nn*(nn-1);
          }
      }
  }



  // This method applies the tensor product operation to produce face values
  // out from cell values. As opposed to the apply_tensor_product method, this
  // method assumes that the directions orthogonal to the face have
  // fe_degree+1 degrees of freedom per direction and not n_q_points_1d for
  // those directions lower than the one currently applied
  template <int dim, int fe_degree, typename Number, int face_direction,
            bool dof_to_quad, bool add>
  inline
  void
  apply_tensor_product_face (const Number *shape_data,
                             const Number in [],
                             Number       out [])
  {
    const int n_blocks1 = dim > 1 ? (fe_degree+1) : 1;
    const int n_blocks2 = dim > 2 ? (fe_degree+1) : 1;

    AssertIndexRange (face_direction, dim);
    const int mm     = dof_to_quad ? (fe_degree+1) : 1,
              nn     = dof_to_quad ? 1 : (fe_degree+1);

    const int stride = Utilities::fixed_int_power<fe_degree+1,face_direction>::value;

    for (int i2=0; i2<n_blocks2; ++i2)
      {
        for (int i1=0; i1<n_blocks1; ++i1)
          {
            if (dof_to_quad == true)
              {
                Number res0 = shape_data[0] * in[0];
                for (int ind=1; ind<mm; ++ind)
                  res0 += shape_data[ind] * in[stride*ind];
                if (add == false)
                  out[0]  = res0;
                else
                  out[0] += res0;
              }
            else
              {
                for (int col=0; col<nn; ++col)
                  if (add == false)
                    out[col*stride]  = shape_data[col] * in[0];
                  else
                    out[col*stride] += shape_data[col] * in[0];
              }

            // increment: in regular case, just go to the next point in
            // x-direction. If we are at the end of one chunk in x-dir, need
            // to jump over to the next layer in z-direction
            switch (face_direction)
              {
              case 0:
                in += mm;
                out += nn;
                break;
              case 1:
                ++in;
                ++out;
                // faces 2 and 3 in 3D use local coordinate system zx, which
                // is the other way around compared to the tensor
                // product. Need to take that into account.
                if (dim == 3)
                  {
                    if (dof_to_quad)
                      out += fe_degree;
                    else
                      in += fe_degree;
                  }
                break;
              case 2:
                ++in;
                ++out;
                break;
              default:
                Assert (false, ExcNotImplemented());
              }
          }
        if (face_direction == 1 && dim == 3)
          {
            in += mm*(mm-1);
            out += nn*(nn-1);
            // adjust for local coordinate system zx
            if (dof_to_quad)
              out -= (fe_degree+1)*(fe_degree+1)-1;
            else
              in -= (fe_degree+1)*(fe_degree+1)-1;
          }
      }
  }



  // This class specializes the general application of tensor-product based
  // elements for "symmetric" finite elements, i.e., when the shape functions
  // are symmetric about 0.5 and the quadrature points are, too.
  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  struct EvaluatorTensorProduct<evaluate_symmetric,dim,fe_degree,n_q_points_1d,Number>
  {
    static const unsigned int dofs_per_cell = Utilities::fixed_int_power<fe_degree+1,dim>::value;
    static const unsigned int n_q_points = Utilities::fixed_int_power<n_q_points_1d,dim>::value;

    /**
     * Constructor, taking the data from ShapeInfo
     */
    EvaluatorTensorProduct (const AlignedVector<Number> &shape_values,
                            const AlignedVector<Number> &shape_gradients,
                            const AlignedVector<Number> &shape_hessians)
      :
      shape_values (shape_values.begin()),
      shape_gradients (shape_gradients.begin()),
      shape_hessians (shape_hessians.begin())
    {}

    template <int direction, bool dof_to_quad, bool add>
    void
    values (const Number in [],
            Number       out[]) const;

    template <int direction, bool dof_to_quad, bool add>
    void
    gradients (const Number in [],
               Number       out[]) const;

    template <int direction, bool dof_to_quad, bool add>
    void
    hessians (const Number in [],
              Number       out[]) const;

    const Number *shape_values;
    const Number *shape_gradients;
    const Number *shape_hessians;
  };



  // In this case, the 1D shape values read (sorted lexicographically, rows
  // run over 1D dofs, columns over quadrature points):
  // Q2 --> [ 0.687  0 -0.087 ]
  //        [ 0.4    1  0.4   ]
  //        [-0.087  0  0.687 ]
  // Q3 --> [ 0.66   0.003  0.002  0.049 ]
  //        [ 0.521  1.005 -0.01  -0.230 ]
  //        [-0.230 -0.01   1.005  0.521 ]
  //        [ 0.049  0.002  0.003  0.66  ]
  // Q4 --> [ 0.658  0.022  0 -0.007 -0.032 ]
  //        [ 0.608  1.059  0  0.039  0.176 ]
  //        [-0.409 -0.113  1 -0.113 -0.409 ]
  //        [ 0.176  0.039  0  1.059  0.608 ]
  //        [-0.032 -0.007  0  0.022  0.658 ]
  //
  // In these matrices, we want to use avoid computations involving zeros and
  // ones and in addition use the symmetry in entries to reduce the number of
  // read operations.
  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  template <int direction, bool dof_to_quad, bool add>
  inline
  void
  EvaluatorTensorProduct<evaluate_symmetric,dim,fe_degree,n_q_points_1d,Number>
  ::values (const Number in [],
            Number       out []) const
  {
    AssertIndexRange (direction, dim);
    const int mm     = dof_to_quad ? (fe_degree+1) : n_q_points_1d,
              nn     = dof_to_quad ? n_q_points_1d : (fe_degree+1);
    const int n_cols = nn / 2;
    const int mid    = mm / 2;

    const int n_blocks1 = (dim > 1 ? (direction > 0 ? nn : mm) : 1);
    const int n_blocks2 = (dim > 2 ? (direction > 1 ? nn : mm) : 1);
    const int stride    = Utilities::fixed_int_power<nn,direction>::value;

    for (int i2=0; i2<n_blocks2; ++i2)
      {
        for (int i1=0; i1<n_blocks1; ++i1)
          {
            for (int col=0; col<n_cols; ++col)
              {
                Number val0, val1, in0, in1, res0, res1;
                if (dof_to_quad == true)
                  {
                    val0 = shape_values[col];
                    val1 = shape_values[nn-1-col];
                  }
                else
                  {
                    val0 = shape_values[col*n_q_points_1d];
                    val1 = shape_values[(col+1)*n_q_points_1d-1];
                  }
                if (mid > 0)
                  {
                    in0 = in[0];
                    in1 = in[stride*(mm-1)];
                    res0 = val0 * in0;
                    res1 = val1 * in0;
                    res0 += val1 * in1;
                    res1 += val0 * in1;
                    for (int ind=1; ind<mid; ++ind)
                      {
                        if (dof_to_quad == true)
                          {
                            val0 = shape_values[ind*n_q_points_1d+col];
                            val1 = shape_values[ind*n_q_points_1d+nn-1-col];
                          }
                        else
                          {
                            val0 = shape_values[col*n_q_points_1d+ind];
                            val1 = shape_values[(col+1)*n_q_points_1d-1-ind];
                          }
                        in0 = in[stride*ind];
                        in1 = in[stride*(mm-1-ind)];
                        res0 += val0 * in0;
                        res1 += val1 * in0;
                        res0 += val1 * in1;
                        res1 += val0 * in1;
                      }
                  }
                else
                  res0 = res1 = Number();
                if (dof_to_quad == true)
                  {
                    if (mm % 2 == 1)
                      {
                        val0 = shape_values[mid*n_q_points_1d+col];
                        val1 = val0 * in[stride*mid];
                        res0 += val1;
                        res1 += val1;
                      }
                  }
                else
                  {
                    if (mm % 2 == 1 && nn % 2 == 0)
                      {
                        val0 = shape_values[col*n_q_points_1d+mid];
                        val1 = val0 * in[stride*mid];
                        res0 += val1;
                        res1 += val1;
                      }
                  }
                if (add == false)
                  {
                    out[stride*col]         = res0;
                    out[stride*(nn-1-col)]  = res1;
                  }
                else
                  {
                    out[stride*col]        += res0;
                    out[stride*(nn-1-col)] += res1;
                  }
              }
            if ( dof_to_quad == true && nn%2==1 && mm%2==1 )
              {
                if (add==false)
                  out[stride*n_cols]  = in[stride*mid];
                else
                  out[stride*n_cols] += in[stride*mid];
              }
            else if (dof_to_quad == true && nn%2==1)
              {
                Number res0;
                Number val0  = shape_values[n_cols];
                if (mid > 0)
                  {
                    res0  = in[0] + in[stride*(mm-1)];
                    res0 *= val0;
                    for (int ind=1; ind<mid; ++ind)
                      {
                        val0  = shape_values[ind*n_q_points_1d+n_cols];
                        Number val1  = in[stride*ind] + in[stride*(mm-1-ind)];
                        val1 *= val0;
                        res0 += val1;
                      }
                  }
                else
                  res0 = Number();
                if (add == false)
                  out[stride*n_cols]  = res0;
                else
                  out[stride*n_cols] += res0;
              }
            else if (dof_to_quad == false && nn%2 == 1)
              {
                Number res0;
                if (mid > 0)
                  {
                    Number val0 = shape_values[n_cols*n_q_points_1d];
                    res0 = in[0] + in[stride*(mm-1)];
                    res0 *= val0;
                    for (int ind=1; ind<mid; ++ind)
                      {
                        val0  = shape_values[n_cols*n_q_points_1d+ind];
                        Number val1 = in[stride*ind] + in[stride*(mm-1-ind)];
                        val1 *= val0;
                        res0 += val1;
                      }
                    if (mm % 2)
                      res0 += in[stride*mid];
                  }
                else
                  res0 = in[0];
                if (add == false)
                  out[stride*n_cols]  = res0;
                else
                  out[stride*n_cols] += res0;
              }

            // increment: in regular case, just go to the next point in
            // x-direction. If we are at the end of one chunk in x-dir, need to
            // jump over to the next layer in z-direction
            switch (direction)
              {
              case 0:
                in += mm;
                out += nn;
                break;
              case 1:
              case 2:
                ++in;
                ++out;
                break;
              default:
                Assert (false, ExcNotImplemented());
              }
          }
        if (direction == 1)
          {
            in += nn*(mm-1);
            out += nn*(nn-1);
          }
      }
  }



  // For the specialized loop used for the gradient computation in
  // here, the 1D shape values read (sorted lexicographically, rows
  // run over 1D dofs, columns over quadrature points):
  // Q2 --> [-2.549 -1  0.549 ]
  //        [ 3.098  0 -3.098 ]
  //        [-0.549  1  2.549 ]
  // Q3 --> [-4.315 -1.03  0.5  -0.44  ]
  //        [ 6.07  -1.44 -2.97  2.196 ]
  //        [-2.196  2.97  1.44 -6.07  ]
  //        [ 0.44  -0.5   1.03  4.315 ]
  // Q4 --> [-6.316 -1.3    0.333 -0.353  0.413 ]
  //        [10.111 -2.76  -2.667  2.066 -2.306 ]
  //        [-5.688  5.773  0     -5.773  5.688 ]
  //        [ 2.306 -2.066  2.667  2.76 -10.111 ]
  //        [-0.413  0.353 -0.333 -0.353  0.413 ]
  //
  // In these matrices, we want to use avoid computations involving
  // zeros and ones and in addition use the symmetry in entries to
  // reduce the number of read operations.
  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  template <int direction, bool dof_to_quad, bool add>
  inline
  void
  EvaluatorTensorProduct<evaluate_symmetric,dim,fe_degree,n_q_points_1d,Number>
  ::gradients (const Number in [],
               Number       out []) const
  {
    AssertIndexRange (direction, dim);
    const int mm     = dof_to_quad ? (fe_degree+1) : n_q_points_1d,
              nn     = dof_to_quad ? n_q_points_1d : (fe_degree+1);
    const int n_cols = nn / 2;
    const int mid    = mm / 2;

    const int n_blocks1 = (dim > 1 ? (direction > 0 ? nn : mm) : 1);
    const int n_blocks2 = (dim > 2 ? (direction > 1 ? nn : mm) : 1);
    const int stride    = Utilities::fixed_int_power<nn,direction>::value;

    for (int i2=0; i2<n_blocks2; ++i2)
      {
        for (int i1=0; i1<n_blocks1; ++i1)
          {
            for (int col=0; col<n_cols; ++col)
              {
                Number val0, val1, in0, in1, res0, res1;
                if (dof_to_quad == true)
                  {
                    val0 = shape_gradients[col];
                    val1 = shape_gradients[nn-1-col];
                  }
                else
                  {
                    val0 = shape_gradients[col*n_q_points_1d];
                    val1 = shape_gradients[(nn-col-1)*n_q_points_1d];
                  }
                if (mid > 0)
                  {
                    in0 = in[0];
                    in1 = in[stride*(mm-1)];
                    res0 = val0 * in0;
                    res1 = val1 * in0;
                    res0 -= val1 * in1;
                    res1 -= val0 * in1;
                    for (int ind=1; ind<mid; ++ind)
                      {
                        if (dof_to_quad == true)
                          {
                            val0 = shape_gradients[ind*n_q_points_1d+col];
                            val1 = shape_gradients[ind*n_q_points_1d+nn-1-col];
                          }
                        else
                          {
                            val0 = shape_gradients[col*n_q_points_1d+ind];
                            val1 = shape_gradients[(nn-col-1)*n_q_points_1d+ind];
                          }
                        in0 = in[stride*ind];
                        in1 = in[stride*(mm-1-ind)];
                        res0 += val0 * in0;
                        res1 += val1 * in0;
                        res0 -= val1 * in1;
                        res1 -= val0 * in1;
                      }
                  }
                else
                  res0 = res1 = Number();
                if (mm % 2 == 1)
                  {
                    if (dof_to_quad == true)
                      val0 = shape_gradients[mid*n_q_points_1d+col];
                    else
                      val0 = shape_gradients[col*n_q_points_1d+mid];
                    val1 = val0 * in[stride*mid];
                    res0 += val1;
                    res1 -= val1;
                  }
                if (add == false)
                  {
                    out[stride*col]         = res0;
                    out[stride*(nn-1-col)]  = res1;
                  }
                else
                  {
                    out[stride*col]        += res0;
                    out[stride*(nn-1-col)] += res1;
                  }
              }
            if ( nn%2 == 1 )
              {
                Number val0, res0;
                if (dof_to_quad == true)
                  val0 = shape_gradients[n_cols];
                else
                  val0 = shape_gradients[n_cols*n_q_points_1d];
                res0  = in[0] - in[stride*(mm-1)];
                res0 *= val0;
                for (int ind=1; ind<mid; ++ind)
                  {
                    if (dof_to_quad == true)
                      val0 = shape_gradients[ind*n_q_points_1d+n_cols];
                    else
                      val0 = shape_gradients[n_cols*n_q_points_1d+ind];
                    Number val1  = in[stride*ind] - in[stride*(mm-1-ind)];
                    val1 *= val0;
                    res0 += val1;
                  }
                if (add == false)
                  out[stride*n_cols]  = res0;
                else
                  out[stride*n_cols] += res0;
              }

            // increment: in regular case, just go to the next point in
            // x-direction. for y-part in 3D and if we are at the end of one
            // chunk in x-dir, need to jump over to the next layer in
            // z-direction
            switch (direction)
              {
              case 0:
                in += mm;
                out += nn;
                break;
              case 1:
              case 2:
                ++in;
                ++out;
                break;
              default:
                Assert (false, ExcNotImplemented());
              }
          }

        if (direction == 1)
          {
            in  += nn * (mm-1);
            out += nn * (nn-1);
          }
      }
  }



  // evaluates the given shape data in 1d-3d using the tensor product
  // form assuming the symmetries of unit cell shape hessians for
  // finite elements in FEEvaluation
  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  template <int direction, bool dof_to_quad, bool add>
  inline
  void
  EvaluatorTensorProduct<evaluate_symmetric,dim,fe_degree,n_q_points_1d,Number>
  ::hessians (const Number in [],
              Number       out []) const
  {
    AssertIndexRange (direction, dim);
    const int mm     = dof_to_quad ? (fe_degree+1) : n_q_points_1d,
              nn     = dof_to_quad ? n_q_points_1d : (fe_degree+1);
    const int n_cols = nn / 2;
    const int mid    = mm / 2;

    const int n_blocks1 = (dim > 1 ? (direction > 0 ? nn : mm) : 1);
    const int n_blocks2 = (dim > 2 ? (direction > 1 ? nn : mm) : 1);
    const int stride    = Utilities::fixed_int_power<nn,direction>::value;

    for (int i2=0; i2<n_blocks2; ++i2)
      {
        for (int i1=0; i1<n_blocks1; ++i1)
          {
            for (int col=0; col<n_cols; ++col)
              {
                Number val0, val1, in0, in1, res0, res1;
                if (dof_to_quad == true)
                  {
                    val0 = shape_hessians[col];
                    val1 = shape_hessians[nn-1-col];
                  }
                else
                  {
                    val0 = shape_hessians[col*n_q_points_1d];
                    val1 = shape_hessians[(col+1)*n_q_points_1d-1];
                  }
                if (mid > 0)
                  {
                    in0 = in[0];
                    in1 = in[stride*(mm-1)];
                    res0 = val0 * in0;
                    res1 = val1 * in0;
                    res0 += val1 * in1;
                    res1 += val0 * in1;
                    for (int ind=1; ind<mid; ++ind)
                      {
                        if (dof_to_quad == true)
                          {
                            val0 = shape_hessians[ind*n_q_points_1d+col];
                            val1 = shape_hessians[ind*n_q_points_1d+nn-1-col];
                          }
                        else
                          {
                            val0 = shape_hessians[col*n_q_points_1d+ind];
                            val1 = shape_hessians[(col+1)*n_q_points_1d-1-ind];
                          }
                        in0 = in[stride*ind];
                        in1 = in[stride*(mm-1-ind)];
                        res0 += val0 * in0;
                        res1 += val1 * in0;
                        res0 += val1 * in1;
                        res1 += val0 * in1;
                      }
                  }
                else
                  res0 = res1 = Number();
                if (mm % 2 == 1)
                  {
                    if (dof_to_quad == true)
                      val0 = shape_hessians[mid*n_q_points_1d+col];
                    else
                      val0 = shape_hessians[col*n_q_points_1d+mid];
                    val1 = val0 * in[stride*mid];
                    res0 += val1;
                    res1 += val1;
                  }
                if (add == false)
                  {
                    out[stride*col]         = res0;
                    out[stride*(nn-1-col)]  = res1;
                  }
                else
                  {
                    out[stride*col]        += res0;
                    out[stride*(nn-1-col)] += res1;
                  }
              }
            if ( nn%2 == 1 )
              {
                Number val0, res0;
                if (dof_to_quad == true)
                  val0 = shape_hessians[n_cols];
                else
                  val0 = shape_hessians[n_cols*n_q_points_1d];
                if (mid > 0)
                  {
                    res0  = in[0] + in[stride*(mm-1)];
                    res0 *= val0;
                    for (int ind=1; ind<mid; ++ind)
                      {
                        if (dof_to_quad == true)
                          val0 = shape_hessians[ind*n_q_points_1d+n_cols];
                        else
                          val0 = shape_hessians[n_cols*n_q_points_1d+ind];
                        Number val1  = in[stride*ind] + in[stride*(mm-1-ind)];
                        val1 *= val0;
                        res0 += val1;
                      }
                  }
                else
                  res0 = Number();
                if (mm % 2 == 1)
                  {
                    if (dof_to_quad == true)
                      val0 = shape_hessians[mid*n_q_points_1d+n_cols];
                    else
                      val0 = shape_hessians[n_cols*n_q_points_1d+mid];
                    res0 += val0 * in[stride*mid];
                  }
                if (add == false)
                  out[stride*n_cols]  = res0;
                else
                  out[stride*n_cols] += res0;
              }

            // increment: in regular case, just go to the next point in
            // x-direction. If we are at the end of one chunk in x-dir, need to
            // jump over to the next layer in z-direction
            switch (direction)
              {
              case 0:
                in += mm;
                out += nn;
                break;
              case 1:
              case 2:
                ++in;
                ++out;
                break;
              default:
                Assert (false, ExcNotImplemented());
              }
          }
        if (direction == 1)
          {
            in += nn*(mm-1);
            out += nn*(nn-1);
          }
      }
  }



  // This class implements a different approach to the symmetric case for
  // values, gradients, and Hessians also treated with the above functions: It
  // is possible to reduce the cost per dimension from N^2 to N^2/2, where N
  // is the number of 1D dofs (there are only N^2/2 different entries in the
  // shape matrix, so this is plausible). The approach is based on the idea of
  // applying the operator on the even and odd part of the input vectors
  // separately, given that the shape functions evaluated on quadrature points
  // are symmetric. This method is presented e.g. in the book "Implementing
  // Spectral Methods for Partial Differential Equations" by David A. Kopriva,
  // Springer, 2009, section 3.5.3 (Even-Odd-Decomposition). Even though the
  // experiments in the book say that the method is not efficient for N<20, it
  // is more efficient in the context where the loop bounds are compile-time
  // constants (templates).
  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  struct EvaluatorTensorProduct<evaluate_evenodd,dim,fe_degree,n_q_points_1d,Number>
  {
    static const unsigned int dofs_per_cell = Utilities::fixed_int_power<fe_degree+1,dim>::value;
    static const unsigned int n_q_points = Utilities::fixed_int_power<n_q_points_1d,dim>::value;

    /**
     * Empty constructor. Does nothing. Be careful when using 'values' and
     * related methods because they need to be filled with the other pointer
     */
    EvaluatorTensorProduct ()
      :
      shape_values (0),
      shape_gradients (0),
      shape_hessians (0)
    {}

    /**
     * Constructor, taking the data from ShapeInfo (using the even-odd
     * variants stored there)
     */
    EvaluatorTensorProduct (const AlignedVector<Number> &shape_values,
                            const AlignedVector<Number> &shape_gradients,
                            const AlignedVector<Number> &shape_hessians)
      :
      shape_values (shape_values.begin()),
      shape_gradients (shape_gradients.begin()),
      shape_hessians (shape_hessians.begin())
    {}

    template <int direction, bool dof_to_quad, bool add>
    void
    values (const Number in [],
            Number       out[]) const
    {
      apply<direction,dof_to_quad,add,0>(shape_values, in, out);
    }

    template <int direction, bool dof_to_quad, bool add>
    void
    gradients (const Number in [],
               Number       out[]) const
    {
      apply<direction,dof_to_quad,add,1>(shape_gradients, in, out);
    }

    template <int direction, bool dof_to_quad, bool add>
    void
    hessians (const Number in [],
              Number       out[]) const
    {
      apply<direction,dof_to_quad,add,2>(shape_hessians, in, out);
    }

    template <int direction, bool dof_to_quad, bool add, int type>
    static void apply (const Number *shape_data,
                       const Number  in [],
                       Number        out []);

    const Number *shape_values;
    const Number *shape_gradients;
    const Number *shape_hessians;
  };



  template <int dim, int fe_degree, int n_q_points_1d, typename Number>
  template <int direction, bool dof_to_quad, bool add, int type>
  inline
  void
  EvaluatorTensorProduct<evaluate_evenodd,dim,fe_degree,n_q_points_1d,Number>
  ::apply (const Number *shapes,
           const Number  in [],
           Number        out [])
  {
    AssertIndexRange (type, 3);
    AssertIndexRange (direction, dim);
    const int mm     = dof_to_quad ? (fe_degree+1) : n_q_points_1d,
              nn     = dof_to_quad ? n_q_points_1d : (fe_degree+1);
    const int n_cols = nn / 2;
    const int mid    = mm / 2;

    const int n_blocks1 = (dim > 1 ? (direction > 0 ? nn : mm) : 1);
    const int n_blocks2 = (dim > 2 ? (direction > 1 ? nn : mm) : 1);
    const int stride    = Utilities::fixed_int_power<nn,direction>::value;

    const int offset = (n_q_points_1d+1)/2;

    // this code may look very inefficient at first sight due to the many
    // different cases with if's at the innermost loop part, but all of the
    // conditionals can be evaluated at compile time because they are
    // templates, so the compiler should optimize everything away
    for (int i2=0; i2<n_blocks2; ++i2)
      {
        for (int i1=0; i1<n_blocks1; ++i1)
          {
            Number xp[mid>0?mid:1], xm[mid>0?mid:1];
            for (int i=0; i<mid; ++i)
              {
                if (dof_to_quad == true && type == 1)
                  {
                    xp[i] = in[stride*i] - in[stride*(mm-1-i)];
                    xm[i] = in[stride*i] + in[stride*(mm-1-i)];
                  }
                else
                  {
                    xp[i] = in[stride*i] + in[stride*(mm-1-i)];
                    xm[i] = in[stride*i] - in[stride*(mm-1-i)];
                  }
              }
            for (int col=0; col<n_cols; ++col)
              {
                Number r0, r1;
                if (mid > 0)
                  {
                    if (dof_to_quad == true)
                      {
                        r0 = shapes[col]                    * xp[0];
                        r1 = shapes[fe_degree*offset + col] * xm[0];
                      }
                    else
                      {
                        r0 = shapes[col*offset]             * xp[0];
                        r1 = shapes[(fe_degree-col)*offset] * xm[0];
                      }
                    for (int ind=1; ind<mid; ++ind)
                      {
                        if (dof_to_quad == true)
                          {
                            r0 += shapes[ind*offset+col]             * xp[ind];
                            r1 += shapes[(fe_degree-ind)*offset+col] * xm[ind];
                          }
                        else
                          {
                            r0 += shapes[col*offset+ind]             * xp[ind];
                            r1 += shapes[(fe_degree-col)*offset+ind] * xm[ind];
                          }
                      }
                  }
                else
                  r0 = r1 = Number();
                if (mm % 2 == 1 && dof_to_quad == true)
                  {
                    if (type == 1)
                      r1 += shapes[mid*offset+col] * in[stride*mid];
                    else
                      r0 += shapes[mid*offset+col] * in[stride*mid];
                  }
                else if (mm % 2 == 1 && (nn % 2 == 0 || type > 0))
                  r0 += shapes[col*offset+mid] * in[stride*mid];

                if (add == false)
                  {
                    out[stride*col]         = r0 + r1;
                    if (type == 1 && dof_to_quad == false)
                      out[stride*(nn-1-col)]  = r1 - r0;
                    else
                      out[stride*(nn-1-col)]  = r0 - r1;
                  }
                else
                  {
                    out[stride*col]        += r0 + r1;
                    if (type == 1 && dof_to_quad == false)
                      out[stride*(nn-1-col)] += r1 - r0;
                    else
                      out[stride*(nn-1-col)] += r0 - r1;
                  }
              }
            if ( type == 0 && dof_to_quad == true && nn%2==1 && mm%2==1 )
              {
                if (add==false)
                  out[stride*n_cols]  = in[stride*mid];
                else
                  out[stride*n_cols] += in[stride*mid];
              }
            else if (dof_to_quad == true && nn%2==1)
              {
                Number r0;
                if (mid > 0)
                  {
                    r0  = shapes[n_cols] * xp[0];
                    for (int ind=1; ind<mid; ++ind)
                      r0 += shapes[ind*offset+n_cols] * xp[ind];
                  }
                else
                  r0 = Number();
                if (type != 1 && mm % 2 == 1)
                  r0 += shapes[mid*offset+n_cols] * in[stride*mid];

                if (add == false)
                  out[stride*n_cols]  = r0;
                else
                  out[stride*n_cols] += r0;
              }
            else if (dof_to_quad == false && nn%2 == 1)
              {
                Number r0;
                if (mid > 0)
                  {
                    if (type == 1)
                      {
                        r0 = shapes[n_cols*offset] * xm[0];
                        for (int ind=1; ind<mid; ++ind)
                          r0 += shapes[n_cols*offset+ind] * xm[ind];
                      }
                    else
                      {
                        r0 = shapes[n_cols*offset] * xp[0];
                        for (int ind=1; ind<mid; ++ind)
                          r0 += shapes[n_cols*offset+ind] * xp[ind];
                      }
                  }
                else
                  r0 = Number();

                if (type == 0 && mm % 2 == 1)
                  r0 += in[stride*mid];
                else if (type == 2 && mm % 2 == 1)
                  r0 += shapes[n_cols*offset+mid] * in[stride*mid];

                if (add == false)
                  out[stride*n_cols]  = r0;
                else
                  out[stride*n_cols] += r0;
              }

            // increment: in regular case, just go to the next point in
            // x-direction. If we are at the end of one chunk in x-dir, need to
            // jump over to the next layer in z-direction
            switch (direction)
              {
              case 0:
                in += mm;
                out += nn;
                break;
              case 1:
              case 2:
                ++in;
                ++out;
                break;
              default:
                Assert (false, ExcNotImplemented());
              }
          }
        if (direction == 1)
          {
            in += nn*(mm-1);
            out += nn*(nn-1);
          }
      }
  }



  // Select evaluator type from element shape function type
  template <MatrixFreeFunctions::ElementType element, bool is_long>
  struct EvaluatorSelector {};

  template <bool is_long>
  struct EvaluatorSelector<MatrixFreeFunctions::tensor_general,is_long>
  {
    static const EvaluatorVariant variant = evaluate_general;
  };

  template <>
  struct EvaluatorSelector<MatrixFreeFunctions::tensor_symmetric,false>
  {
    static const EvaluatorVariant variant = evaluate_symmetric;
  };

  template <> struct EvaluatorSelector<MatrixFreeFunctions::tensor_symmetric,true>
  {
    static const EvaluatorVariant variant = evaluate_evenodd;
  };

  template <bool is_long>
  struct EvaluatorSelector<MatrixFreeFunctions::truncated_tensor,is_long>
  {
    static const EvaluatorVariant variant = evaluate_general;
  };

  template <>
  struct EvaluatorSelector<MatrixFreeFunctions::tensor_symmetric_plus_dg0,false>
  {
    static const EvaluatorVariant variant = evaluate_symmetric;
  };

  template <>
  struct EvaluatorSelector<MatrixFreeFunctions::tensor_symmetric_plus_dg0,true>
  {
    static const EvaluatorVariant variant = evaluate_evenodd;
  };

  template <bool is_long>
  struct EvaluatorSelector<MatrixFreeFunctions::tensor_gausslobatto,is_long>
  {
    static const EvaluatorVariant variant = evaluate_evenodd;
  };



  // This struct performs the evaluation of function values, gradients and
  // Hessians for tensor-product finite elements. The operation is used for
  // both the symmetric and non-symmetric case, which use different apply
  // functions 'values', 'gradients' in the individual coordinate
  // directions. The apply functions for values are provided through one of
  // the template classes EvaluatorTensorProduct which in turn are selected
  // from the MatrixFreeFunctions::ElementType template argument.
  //
  // There is a specialization made for Gauss-Lobatto elements further down
  // where the 'values' operation is identity, which allows us to write
  // shorter code.
  template <MatrixFreeFunctions::ElementType type, int dim, int fe_degree,
            int n_q_points_1d, int n_components, typename Number>
  struct FEEvaluationImpl
  {
    static
    void evaluate (const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
                   VectorizedArray<Number> *values_dofs_actual[],
                   VectorizedArray<Number> *values_quad[],
                   VectorizedArray<Number> *gradients_quad[][dim],
                   VectorizedArray<Number> *hessians_quad[][(dim*(dim+1))/2],
                   const bool               evaluate_val,
                   const bool               evaluate_grad,
                   const bool               evaluate_lapl);

    static
    void integrate (const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
                    VectorizedArray<Number> *values_dofs_actual[],
                    VectorizedArray<Number> *values_quad[],
                    VectorizedArray<Number> *gradients_quad[][dim],
                    const bool               evaluate_val,
                    const bool               evaluate_grad);
  };


  template <MatrixFreeFunctions::ElementType type, int dim, int fe_degree,
            int n_q_points_1d, int n_components, typename Number>
  inline
  void
  FEEvaluationImpl<type,dim,fe_degree,n_q_points_1d,n_components,Number>
  ::evaluate (const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
              VectorizedArray<Number> *values_dofs_actual[],
              VectorizedArray<Number> *values_quad[],
              VectorizedArray<Number> *gradients_quad[][dim],
              VectorizedArray<Number> *hessians_quad[][(dim*(dim+1))/2],
              const bool               evaluate_val,
              const bool               evaluate_grad,
              const bool               evaluate_lapl)
  {
    if (evaluate_val == false && evaluate_grad == false && evaluate_lapl == false)
      return;

    const EvaluatorVariant variant =
      EvaluatorSelector<type,(fe_degree+n_q_points_1d>4)>::variant;
    typedef EvaluatorTensorProduct<variant, dim, fe_degree, n_q_points_1d,
            VectorizedArray<Number> > Eval;
    Eval eval (variant == evaluate_evenodd ? shape_info.shape_val_evenodd :
               shape_info.shape_values,
               variant == evaluate_evenodd ? shape_info.shape_gra_evenodd :
               shape_info.shape_gradients,
               variant == evaluate_evenodd ? shape_info.shape_hes_evenodd :
               shape_info.shape_hessians);

    const unsigned int temp_size = Eval::dofs_per_cell > Eval::n_q_points ?
                                   Eval::dofs_per_cell : Eval::n_q_points;

    VectorizedArray<Number> **values_dofs = values_dofs_actual;
    VectorizedArray<Number> data_array[type!=MatrixFreeFunctions::truncated_tensor ? 1 :
                                       n_components*Eval::dofs_per_cell];
    VectorizedArray<Number> *expanded_dof_values[n_components];
    if (type == MatrixFreeFunctions::truncated_tensor)
      {
        for (unsigned int c=0; c<n_components; ++c)
          expanded_dof_values[c] = &data_array[c*Eval::dofs_per_cell];
        values_dofs = expanded_dof_values;

        unsigned int count_p = 0, count_q = 0;
        for (unsigned int i=0; i<(dim>2?fe_degree+1:1); ++i)
          {
            for (unsigned int j=0; j<(dim>1?fe_degree+1-i:1); ++j)
              {
                for (unsigned int k=0; k<fe_degree+1-j-i; ++k, ++count_p, ++count_q)
                  for (unsigned int c=0; c<n_components; ++c)
                    expanded_dof_values[c][count_q] = values_dofs_actual[c][count_p];
                for (unsigned int k=fe_degree+1-j-i; k<fe_degree+1; ++k, ++count_q)
                  for (unsigned int c=0; c<n_components; ++c)
                    expanded_dof_values[c][count_q] = VectorizedArray<Number>();
              }
            for (unsigned int j=fe_degree+1-i; j<fe_degree+1; ++j)
              for (unsigned int k=0; k<fe_degree+1; ++k, ++count_q)
                for (unsigned int c=0; c<n_components; ++c)
                  expanded_dof_values[c][count_q] = VectorizedArray<Number>();
          }
        AssertDimension(count_q, Eval::dofs_per_cell);
      }

    // These avoid compiler errors; they are only used in sensible context but
    // compilers typically cannot detect when we access something like
    // gradients_quad[2] only for dim==3.
    const unsigned int d1 = dim>1?1:0;
    const unsigned int d2 = dim>2?2:0;
    const unsigned int d3 = dim>2?3:0;
    const unsigned int d4 = dim>2?4:0;
    const unsigned int d5 = dim>2?5:0;

    switch (dim)
      {
      case 1:
        for (unsigned int c=0; c<n_components; c++)
          {
            if (evaluate_val == true)
              eval.template values<0,true,false> (values_dofs[c], values_quad[c]);
            if (evaluate_grad == true)
              eval.template gradients<0,true,false>(values_dofs[c], gradients_quad[c][0]);
            if (evaluate_lapl == true)
              eval.template hessians<0,true,false> (values_dofs[c], hessians_quad[c][0]);
          }
        break;

      case 2:
        for (unsigned int c=0; c<n_components; c++)
          {
            VectorizedArray<Number> temp1[temp_size];
            VectorizedArray<Number> temp2[temp_size];

            // grad x
            if (evaluate_grad == true)
              {
                eval.template gradients<0,true,false> (values_dofs[c], temp1);
                eval.template values<1,true,false> (temp1, gradients_quad[c][0]);
              }
            if (evaluate_lapl == true)
              {
                // grad xy
                if (evaluate_grad == false)
                  eval.template gradients<0,true,false>(values_dofs[c], temp1);
                eval.template gradients<1,true,false>  (temp1, hessians_quad[c][d1+d1]);

                // grad xx
                eval.template hessians<0,true,false>(values_dofs[c], temp1);
                eval.template values<1,true,false>  (temp1, hessians_quad[c][0]);
              }

            // grad y
            eval.template values<0,true,false> (values_dofs[c], temp1);
            if (evaluate_grad == true)
              eval.template gradients<1,true,false> (temp1, gradients_quad[c][d1]);

            // grad yy
            if (evaluate_lapl == true)
              eval.template hessians<1,true,false> (temp1, hessians_quad[c][d1]);

            // val: can use values applied in x
            if (evaluate_val == true)
              eval.template values<1,true,false> (temp1, values_quad[c]);
          }
        break;

      case 3:
        for (unsigned int c=0; c<n_components; c++)
          {
            VectorizedArray<Number> temp1[temp_size];
            VectorizedArray<Number> temp2[temp_size];

            if (evaluate_grad == true)
              {
                // grad x
                eval.template gradients<0,true,false> (values_dofs[c], temp1);
                eval.template values<1,true,false> (temp1, temp2);
                eval.template values<2,true,false> (temp2, gradients_quad[c][0]);
              }

            if (evaluate_lapl == true)
              {
                // grad xz
                if (evaluate_grad == false)
                  {
                    eval.template gradients<0,true,false> (values_dofs[c], temp1);
                    eval.template values<1,true,false> (temp1, temp2);
                  }
                eval.template gradients<2,true,false> (temp2, hessians_quad[c][d4]);

                // grad xy
                eval.template gradients<1,true,false> (temp1, temp2);
                eval.template values<2,true,false> (temp2, hessians_quad[c][d3]);

                // grad xx
                eval.template hessians<0,true,false>(values_dofs[c], temp1);
                eval.template values<1,true,false>  (temp1, temp2);
                eval.template values<2,true,false>  (temp2, hessians_quad[c][0]);
              }

            // grad y
            eval.template values<0,true,false> (values_dofs[c], temp1);
            if (evaluate_grad == true)
              {
                eval.template gradients<1,true,false>(temp1, temp2);
                eval.template values<2,true,false>   (temp2, gradients_quad[c][d1]);
              }

            if (evaluate_lapl == true)
              {
                // grad yz
                if (evaluate_grad == false)
                  eval.template gradients<1,true,false>(temp1, temp2);
                eval.template gradients<2,true,false>  (temp2, hessians_quad[c][d5]);

                // grad yy
                eval.template hessians<1,true,false> (temp1, temp2);
                eval.template values<2,true,false> (temp2, hessians_quad[c][d1]);
              }

            // grad z: can use the values applied in x direction stored in temp1
            eval.template values<1,true,false> (temp1, temp2);
            if (evaluate_grad == true)
              eval.template gradients<2,true,false> (temp2, gradients_quad[c][d2]);

            // grad zz: can use the values applied in x and y direction stored
            // in temp2
            if (evaluate_lapl == true)
              eval.template hessians<2,true,false>(temp2, hessians_quad[c][d2]);

            // val: can use the values applied in x & y direction stored in temp2
            if (evaluate_val == true)
              eval.template values<2,true,false> (temp2, values_quad[c]);
          }
        break;

      default:
        AssertThrow(false, ExcNotImplemented());
      }

    // case additional dof for FE_Q_DG0: add values; gradients and second
    // derivatives evaluate to zero
    if (type == MatrixFreeFunctions::tensor_symmetric_plus_dg0 && evaluate_val)
      for (unsigned int c=0; c<n_components; ++c)
        for (unsigned int q=0; q<Eval::n_q_points; ++q)
          values_quad[c][q] += values_dofs[c][Eval::dofs_per_cell];
  }



  template <MatrixFreeFunctions::ElementType type, int dim, int fe_degree,
            int n_q_points_1d, int n_components, typename Number>
  inline
  void
  FEEvaluationImpl<type,dim,fe_degree,n_q_points_1d,n_components,Number>
  ::integrate (const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
               VectorizedArray<Number> *values_dofs_actual[],
               VectorizedArray<Number> *values_quad[],
               VectorizedArray<Number> *gradients_quad[][dim],
               const bool               integrate_val,
               const bool               integrate_grad)
  {
    const EvaluatorVariant variant =
      EvaluatorSelector<type,(fe_degree+n_q_points_1d>4)>::variant;
    typedef EvaluatorTensorProduct<variant, dim, fe_degree, n_q_points_1d,
            VectorizedArray<Number> > Eval;
    Eval eval (variant == evaluate_evenodd ? shape_info.shape_val_evenodd :
               shape_info.shape_values,
               variant == evaluate_evenodd ? shape_info.shape_gra_evenodd :
               shape_info.shape_gradients,
               variant == evaluate_evenodd ? shape_info.shape_hes_evenodd :
               shape_info.shape_hessians);

    const unsigned int temp_size = Eval::dofs_per_cell > Eval::n_q_points ?
                                   Eval::dofs_per_cell : Eval::n_q_points;
    VectorizedArray<Number> temp1[temp_size];
    VectorizedArray<Number> temp2[temp_size];

    // expand dof_values to tensor product for truncated tensor products
    VectorizedArray<Number> **values_dofs = values_dofs_actual;
    VectorizedArray<Number> data_array[type!=MatrixFreeFunctions::truncated_tensor ? 1 :
                                       n_components*Eval::dofs_per_cell];
    VectorizedArray<Number> *expanded_dof_values[n_components];
    if (type == MatrixFreeFunctions::truncated_tensor)
      {
        for (unsigned int c=0; c<n_components; ++c)
          expanded_dof_values[c] = &data_array[c*Eval::dofs_per_cell];
        values_dofs = expanded_dof_values;
      }

    // These avoid compiler errors; they are only used in sensible context but
    // compilers typically cannot detect when we access something like
    // gradients_quad[2] only for dim==3.
    const unsigned int d1 = dim>1?1:0;
    const unsigned int d2 = dim>2?2:0;

    switch (dim)
      {
      case 1:
        for (unsigned int c=0; c<n_components; c++)
          {
            if (integrate_val == true)
              eval.template values<0,false,false> (values_quad[c], values_dofs[c]);
            if (integrate_grad == true)
              {
                if (integrate_val == true)
                  eval.template gradients<0,false,true> (gradients_quad[c][0], values_dofs[c]);
                else
                  eval.template gradients<0,false,false> (gradients_quad[c][0], values_dofs[c]);
              }
          }
        break;

      case 2:
        for (unsigned int c=0; c<n_components; c++)
          {
            if (integrate_val == true)
              {
                // val
                eval.template values<0,false,false> (values_quad[c], temp1);
                //grad x
                if (integrate_grad == true)
                  eval.template gradients<0,false,true> (gradients_quad[c][0], temp1);
                eval.template values<1,false,false>(temp1, values_dofs[c]);
              }
            if (integrate_grad == true)
              {
                // grad y
                eval.template values<0,false,false>  (gradients_quad[c][d1], temp1);
                if (integrate_val == false)
                  {
                    eval.template gradients<1,false,false>(temp1, values_dofs[c]);
                    //grad x
                    eval.template gradients<0,false,false> (gradients_quad[c][0], temp1);
                    eval.template values<1,false,true> (temp1, values_dofs[c]);
                  }
                else
                  eval.template gradients<1,false,true>(temp1, values_dofs[c]);
              }
          }
        break;

      case 3:
        for (unsigned int c=0; c<n_components; c++)
          {
            if (integrate_val == true)
              {
                // val
                eval.template values<0,false,false> (values_quad[c], temp1);
                //grad x: can sum to temporary value in temp1
                if (integrate_grad == true)
                  eval.template gradients<0,false,true> (gradients_quad[c][0], temp1);
                eval.template values<1,false,false>(temp1, temp2);
                eval.template values<0,false,false> (gradients_quad[c][d1], temp1);
                if (integrate_grad == true)
                  eval.template gradients<1,false,true>(temp1, temp2);
                eval.template values<2,false,false> (temp2, values_dofs[c]);
              }
            else if (integrate_grad == true)
              {
                eval.template gradients<0,false,false>(gradients_quad[c][0], temp1);
                eval.template values<1,false,false> (temp1, temp2);
                eval.template values<0,false,false> (gradients_quad[c][d1], temp1);
                eval.template gradients<1,false,true>(temp1, temp2);
                eval.template values<2,false,false> (temp2, values_dofs[c]);
              }
            if (integrate_grad == true)
              {
                // grad z: can sum to temporary x and y value in output
                eval.template values<0,false,false> (gradients_quad[c][d2], temp1);
                eval.template values<1,false,false> (temp1, temp2);
                eval.template gradients<2,false,true> (temp2, values_dofs[c]);
              }
          }
        break;

      default:
        AssertThrow(false, ExcNotImplemented());
      }

    // case FE_Q_DG0: add values, gradients and second derivatives are zero
    if (type == MatrixFreeFunctions::tensor_symmetric_plus_dg0)
      {
        if (integrate_val)
          for (unsigned int c=0; c<n_components; ++c)
            {
              values_dofs[c][Eval::dofs_per_cell] = values_quad[c][0];
              for (unsigned int q=1; q<Eval::n_q_points; ++q)
                values_dofs[c][Eval::dofs_per_cell] += values_quad[c][q];
            }
        else
          for (unsigned int c=0; c<n_components; ++c)
            values_dofs[c][Eval::dofs_per_cell] = VectorizedArray<Number>();
      }

    if (type == MatrixFreeFunctions::truncated_tensor)
      {
        unsigned int count_p = 0, count_q = 0;
        for (unsigned int i=0; i<(dim>2?fe_degree+1:1); ++i)
          {
            for (unsigned int j=0; j<(dim>1?fe_degree+1-i:1); ++j)
              {
                for (unsigned int k=0; k<fe_degree+1-j-i; ++k, ++count_p, ++count_q)
                  {
                    for (unsigned int c=0; c<n_components; ++c)
                      values_dofs_actual[c][count_p] = expanded_dof_values[c][count_q];
                  }
                count_q += j+i;
              }
            count_q += i*(fe_degree+1);
          }
        AssertDimension(count_q, Eval::dofs_per_cell);
      }
  }

  // This a specialization for Gauss-Lobatto elements where the 'values'
  // operation is identity, which allows us to write shorter code.
  template <int dim, int fe_degree, int n_q_points_1d, int n_components, typename Number>
  struct FEEvaluationImpl<MatrixFreeFunctions::tensor_gausslobatto, dim,
    fe_degree, n_q_points_1d, n_components, Number>
  {
    static
    void evaluate (const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
                   VectorizedArray<Number> *values_dofs[],
                   VectorizedArray<Number> *values_quad[],
                   VectorizedArray<Number> *gradients_quad[][dim],
                   VectorizedArray<Number> *hessians_quad[][(dim*(dim+1))/2],
                   const bool               evaluate_val,
                   const bool               evaluate_grad,
                   const bool               evaluate_lapl);

    static
    void integrate (const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
                    VectorizedArray<Number> *values_dofs[],
                    VectorizedArray<Number> *values_quad[],
                    VectorizedArray<Number> *gradients_quad[][dim],
                    const bool               integrate_val,
                    const bool               integrate_grad);
  };

  template <int dim, int fe_degree, int n_q_points_1d, int n_components, typename Number>
  inline
  void
  FEEvaluationImpl<MatrixFreeFunctions::tensor_gausslobatto, dim,
                   fe_degree, n_q_points_1d, n_components, Number>
                   ::evaluate (const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
                               VectorizedArray<Number> *values_dofs[],
                               VectorizedArray<Number> *values_quad[],
                               VectorizedArray<Number> *gradients_quad[][dim],
                               VectorizedArray<Number> *hessians_quad[][(dim*(dim+1))/2],
                               const bool               evaluate_val,
                               const bool               evaluate_grad,
                               const bool               evaluate_lapl)
  {
    typedef EvaluatorTensorProduct<evaluate_evenodd, dim, fe_degree, fe_degree+1,
            VectorizedArray<Number> > Eval;
    Eval eval (shape_info.shape_val_evenodd, shape_info.shape_gra_evenodd,
               shape_info.shape_hes_evenodd);

    // These avoid compiler errors; they are only used in sensible context but
    // compilers typically cannot detect when we access something like
    // gradients_quad[2] only for dim==3.
    const unsigned int d1 = dim>1?1:0;
    const unsigned int d2 = dim>2?2:0;
    const unsigned int d3 = dim>2?3:0;
    const unsigned int d4 = dim>2?4:0;
    const unsigned int d5 = dim>2?5:0;

    switch (dim)
      {
      case 1:
        if (evaluate_val == true)
          std::memcpy (values_quad[0], values_dofs[0],
                       eval.dofs_per_cell * n_components *
                       sizeof (values_dofs[0][0]));
        for (unsigned int c=0; c<n_components; c++)
          {
            if (evaluate_grad == true)
              eval.template gradients<0,true,false>(values_dofs[c], gradients_quad[c][0]);
            if (evaluate_lapl == true)
              eval.template hessians<0,true,false> (values_dofs[c], hessians_quad[c][0]);
          }
        break;

      case 2:
        if (evaluate_val == true)
          {
            std::memcpy (values_quad[0], values_dofs[0],
                         Eval::dofs_per_cell * n_components *
                         sizeof (values_dofs[0][0]));
          }
        if (evaluate_grad == true)
          for (unsigned int comp=0; comp<n_components; comp++)
            {
              // grad x
              eval.template gradients<0,true,false> (values_dofs[comp],
                                                     gradients_quad[comp][0]);
              // grad y
              eval.template gradients<1,true,false> (values_dofs[comp],
                                                     gradients_quad[comp][d1]);
            }
        if (evaluate_lapl == true)
          for (unsigned int comp=0; comp<n_components; comp++)
            {
              // hess x
              eval.template hessians<0,true,false> (values_dofs[comp],
                                                    hessians_quad[comp][0]);
              // hess y
              eval.template hessians<1,true,false> (values_dofs[comp],
                                                    hessians_quad[comp][d1]);

              VectorizedArray<Number> temp1[Eval::dofs_per_cell];
              // grad x grad y
              eval.template gradients<0,true,false> (values_dofs[comp], temp1);
              eval.template gradients<1,true,false> (temp1, hessians_quad[comp][d1+d1]);
            }
        break;

      case 3:
        if (evaluate_val == true)
          {
            std::memcpy (values_quad[0], values_dofs[0],
                         Eval::dofs_per_cell * n_components *
                         sizeof (values_dofs[0][0]));
          }
        if (evaluate_grad == true)
          for (unsigned int comp=0; comp<n_components; comp++)
            {
              // grad x
              eval.template gradients<0,true,false> (values_dofs[comp],
                                                     gradients_quad[comp][0]);
              // grad y
              eval.template gradients<1,true,false> (values_dofs[comp],
                                                     gradients_quad[comp][d1]);
              // grad y
              eval.template gradients<2,true,false> (values_dofs[comp],
                                                     gradients_quad[comp][d2]);
            }
        if (evaluate_lapl == true)
          for (unsigned int comp=0; comp<n_components; comp++)
            {
              // grad x
              eval.template hessians<0,true,false> (values_dofs[comp],
                                                    hessians_quad[comp][0]);
              // grad y
              eval.template hessians<1,true,false> (values_dofs[comp],
                                                    hessians_quad[comp][d1]);
              // grad y
              eval.template hessians<2,true,false> (values_dofs[comp],
                                                    hessians_quad[comp][d2]);

              VectorizedArray<Number> temp1[Eval::dofs_per_cell];
              // grad xy
              eval.template gradients<0,true,false> (values_dofs[comp], temp1);
              eval.template gradients<1,true,false> (temp1, hessians_quad[comp][d3]);
              // grad xz
              eval.template gradients<2,true,false> (temp1, hessians_quad[comp][d4]);
              // grad yz
              eval.template gradients<1,true,false> (values_dofs[comp], temp1);
              eval.template gradients<2,true,false> (temp1, hessians_quad[comp][d5]);
            }
        break;
      default:
        AssertThrow(false, ExcNotImplemented());
      }
  }

  template <int dim, int fe_degree, int n_q_points_1d, int n_components, typename Number>
  inline
  void
  FEEvaluationImpl<MatrixFreeFunctions::tensor_gausslobatto, dim,
                   fe_degree, n_q_points_1d, n_components, Number>
                   ::integrate (const MatrixFreeFunctions::ShapeInfo<Number> &shape_info,
                                VectorizedArray<Number> *values_dofs[],
                                VectorizedArray<Number> *values_quad[],
                                VectorizedArray<Number> *gradients_quad[][dim],
                                const bool               integrate_val,
                                const bool               integrate_grad)
  {
    typedef EvaluatorTensorProduct<evaluate_evenodd, dim, fe_degree, fe_degree+1,
            VectorizedArray<Number> > Eval;
    Eval eval (shape_info.shape_val_evenodd, shape_info.shape_gra_evenodd,
               shape_info.shape_hes_evenodd);

    // These avoid compiler errors; they are only used in sensible context but
    // compilers typically cannot detect when we access something like
    // gradients_quad[2] only for dim==3.
    const unsigned int d1 = dim>1?1:0;
    const unsigned int d2 = dim>2?2:0;

    if (integrate_val == true)
      std::memcpy (values_dofs[0], values_quad[0],
                   Eval::dofs_per_cell * n_components *
                   sizeof (values_dofs[0][0]));
    switch (dim)
      {
      case 1:
        for (unsigned int c=0; c<n_components; c++)
          {
            if (integrate_grad == true)
              {
                if (integrate_val == true)
                  eval.template gradients<0,false,true> (gradients_quad[c][0],
                                                         values_dofs[c]);
                else
                  eval.template gradients<0,false,false> (gradients_quad[c][0],
                                                          values_dofs[c]);
              }
          }

        break;
      case 2:
        if (integrate_grad == true)
          for (unsigned int comp=0; comp<n_components; comp++)
            {
              // grad x: If integrate_val == true we have to add to the
              // previous output
              if (integrate_val == true)
                eval.template gradients<0, false, true> (gradients_quad[comp][0],
                                                         values_dofs[comp]);
              else
                eval.template gradients<0, false, false> (gradients_quad[comp][0],
                                                          values_dofs[comp]);

              // grad y
              eval.template gradients<1, false, true> (gradients_quad[comp][d1],
                                                       values_dofs[comp]);
            }
        break;

      case 3:
        if (integrate_grad == true)
          for (unsigned int comp=0; comp<n_components; comp++)
            {
              // grad x: If integrate_val == true we have to add to the
              // previous output
              if (integrate_val == true)
                eval.template gradients<0, false, true> (gradients_quad[comp][0],
                                                         values_dofs[comp]);
              else
                eval.template gradients<0, false, false> (gradients_quad[comp][0],
                                                          values_dofs[comp]);

              // grad y
              eval.template gradients<1, false, true> (gradients_quad[comp][d1],
                                                       values_dofs[comp]);

              // grad z
              eval.template gradients<2, false, true> (gradients_quad[comp][d2],
                                                       values_dofs[comp]);
            }
        break;

      default:
        AssertThrow(false, ExcNotImplemented());
      }
  }

} // end of namespace internal



/*-------------------------- FEEvaluation -----------------------------------*/


template <int dim, int fe_degree,  int n_q_points_1d, int n_components_,
          typename Number>
inline
FEEvaluation<dim,fe_degree,n_q_points_1d,n_components_,Number>
::FEEvaluation (const MatrixFree<dim,Number> &data_in,
                const unsigned int fe_no,
                const unsigned int quad_no)
  :
  BaseClass (data_in, fe_no, quad_no, fe_degree, n_q_points),
  dofs_per_cell (this->data->dofs_per_cell)
{
  check_template_arguments(fe_no);
  set_data_pointers();
}



template <int dim, int fe_degree,  int n_q_points_1d, int n_components_,
          typename Number>
inline
FEEvaluation<dim,fe_degree,n_q_points_1d,n_components_,Number>
::FEEvaluation (const MappingFEEvaluation<dim,Number> &geometry,
                const DoFHandler<dim>               &dof_handler,
                const unsigned int                   first_selected_component)
  :
  BaseClass (geometry, dof_handler, first_selected_component),
  dofs_per_cell (this->data->dofs_per_cell)
{
  check_template_arguments(numbers::invalid_unsigned_int);
  set_data_pointers();
}



template <int dim, int fe_degree,  int n_q_points_1d, int n_components_,
          typename Number>
inline
FEEvaluation<dim,fe_degree,n_q_points_1d,n_components_,Number>
::FEEvaluation (const FEEvaluation &other)
  :
  BaseClass (other),
  dofs_per_cell (this->data->dofs_per_cell)
{
  set_data_pointers();
}



template <int dim, int fe_degree,  int n_q_points_1d, int n_components_,
          typename Number>
inline
void
FEEvaluation<dim,fe_degree,n_q_points_1d,n_components_,Number>
::set_data_pointers()
{
  AssertIndexRange(this->data->dofs_per_cell, tensor_dofs_per_cell+2);
  const unsigned int desired_dofs_per_cell = this->data->dofs_per_cell;

  // set the pointers to the correct position in the data array
  for (unsigned int c=0; c<n_components_; ++c)
    {
      this->values_dofs[c] = &my_data_array[c*desired_dofs_per_cell];
      this->values_quad[c] = &my_data_array[n_components*desired_dofs_per_cell+c*n_q_points];
      for (unsigned int d=0; d<dim; ++d)
        this->gradients_quad[c][d] = &my_data_array[n_components*(desired_dofs_per_cell+
                                                                  n_q_points)
                                                    +
                                                    (c*dim+d)*n_q_points];
      for (unsigned int d=0; d<(dim*dim+dim)/2; ++d)
        this->hessians_quad[c][d] = &my_data_array[n_components*((dim+1)*n_q_points+
                                                                 desired_dofs_per_cell)
                                                   +
                                                   (c*(dim*dim+dim)+d)*n_q_points];
    }

  switch (this->data->element_type)
    {
    case internal::MatrixFreeFunctions::tensor_symmetric:
      evaluate_funct =
        internal::FEEvaluationImpl<internal::MatrixFreeFunctions::tensor_symmetric,
        dim, fe_degree, n_q_points_1d, n_components_,
        Number>::evaluate;
      integrate_funct =
        internal::FEEvaluationImpl<internal::MatrixFreeFunctions::tensor_symmetric,
        dim, fe_degree, n_q_points_1d, n_components_,
        Number>::integrate;
      break;

    case internal::MatrixFreeFunctions::tensor_symmetric_plus_dg0:
      evaluate_funct =
        internal::FEEvaluationImpl<internal::MatrixFreeFunctions::tensor_symmetric_plus_dg0,
        dim, fe_degree, n_q_points_1d, n_components_,
        Number>::evaluate;
      integrate_funct =
        internal::FEEvaluationImpl<internal::MatrixFreeFunctions::tensor_symmetric_plus_dg0,
        dim, fe_degree, n_q_points_1d, n_components_,
        Number>::integrate;
      break;

    case internal::MatrixFreeFunctions::tensor_general:
      evaluate_funct =
        internal::FEEvaluationImpl<internal::MatrixFreeFunctions::tensor_general,
        dim, fe_degree, n_q_points_1d, n_components_,
        Number>::evaluate;
      integrate_funct =
        internal::FEEvaluationImpl<internal::MatrixFreeFunctions::tensor_general,
        dim, fe_degree, n_q_points_1d, n_components_,
        Number>::integrate;
      break;

    case internal::MatrixFreeFunctions::tensor_gausslobatto:
      evaluate_funct =
        internal::FEEvaluationImpl<internal::MatrixFreeFunctions::tensor_gausslobatto,
        dim, fe_degree, n_q_points_1d, n_components_,
        Number>::evaluate;
      integrate_funct =
        internal::FEEvaluationImpl<internal::MatrixFreeFunctions::tensor_gausslobatto,
        dim, fe_degree, n_q_points_1d, n_components_,
        Number>::integrate;
      break;

    case internal::MatrixFreeFunctions::truncated_tensor:
      evaluate_funct =
        internal::FEEvaluationImpl<internal::MatrixFreeFunctions::truncated_tensor,
        dim, fe_degree, n_q_points_1d, n_components_,
        Number>::evaluate;
      integrate_funct =
        internal::FEEvaluationImpl<internal::MatrixFreeFunctions::truncated_tensor,
        dim, fe_degree, n_q_points_1d, n_components_,
        Number>::integrate;
      break;

    default:
      AssertThrow(false, ExcNotImplemented());
    }

}



template <int dim, int fe_degree,  int n_q_points_1d, int n_components_,
          typename Number>
inline
void
FEEvaluation<dim,fe_degree,n_q_points_1d,n_components_,Number>
::check_template_arguments(const unsigned int fe_no)
{
#ifdef DEBUG
  // print error message when the dimensions do not match. Propose a possible
  // fix
  if (fe_degree != this->data->fe_degree
      ||
      n_q_points != this->data->n_q_points)
    {
      std::string message =
        "-------------------------------------------------------\n";
      message += "Illegal arguments in constructor/wrong template arguments!\n";
      message += "    Called -->   FEEvaluation<dim,";
      message += Utilities::int_to_string(fe_degree) + ",";
      message += Utilities::int_to_string(n_q_points_1d);
      message += "," + Utilities::int_to_string(n_components);
      message += ",Number>(data, ";
      message += Utilities::int_to_string(fe_no) + ", ";
      message += Utilities::int_to_string(this->quad_no) + ")\n";

      // check whether some other vector component has the correct number of
      // points
      unsigned int proposed_dof_comp = numbers::invalid_unsigned_int,
                   proposed_quad_comp = numbers::invalid_unsigned_int;
      if (fe_no != numbers::invalid_unsigned_int)
        {
          if (fe_degree == this->data->fe_degree)
            proposed_dof_comp = fe_no;
          else
            for (unsigned int no=0; no<this->matrix_info->n_components(); ++no)
              if (this->matrix_info->get_shape_info(no,0,this->active_fe_index,0).fe_degree
                  == fe_degree)
                {
                  proposed_dof_comp = no;
                  break;
                }
          if (n_q_points ==
              this->mapping_info->mapping_data_gen[this->quad_no].n_q_points[this->active_quad_index])
            proposed_quad_comp = this->quad_no;
          else
            for (unsigned int no=0; no<this->mapping_info->mapping_data_gen.size(); ++no)
              if (this->mapping_info->mapping_data_gen[no].n_q_points[this->active_quad_index]
                  == n_q_points)
                {
                  proposed_quad_comp = no;
                  break;
                }
        }
      if (proposed_dof_comp  != numbers::invalid_unsigned_int &&
          proposed_quad_comp != numbers::invalid_unsigned_int)
        {
          if (proposed_dof_comp != fe_no)
            message += "Wrong vector component selection:\n";
          else
            message += "Wrong quadrature formula selection:\n";
          message += "    Did you mean FEEvaluation<dim,";
          message += Utilities::int_to_string(fe_degree) + ",";
          message += Utilities::int_to_string(n_q_points_1d);
          message += "," + Utilities::int_to_string(n_components);
          message += ",Number>(data, ";
          message += Utilities::int_to_string(proposed_dof_comp) + ", ";
          message += Utilities::int_to_string(proposed_quad_comp) + ")?\n";
          std::string correct_pos;
          if (proposed_dof_comp != fe_no)
            correct_pos = " ^ ";
          else
            correct_pos = "   ";
          if (proposed_quad_comp != this->quad_no)
            correct_pos += " ^\n";
          else
            correct_pos += "  \n";
          message += "                                                     " + correct_pos;
        }
      // ok, did not find the numbers specified by the template arguments in
      // the given list. Suggest correct template arguments
      const unsigned int proposed_n_q_points_1d = static_cast<unsigned int>(std::pow(1.001*this->data->n_q_points,1./dim));
      message += "Wrong template arguments:\n";
      message += "    Did you mean FEEvaluation<dim,";
      message += Utilities::int_to_string(this->data->fe_degree) + ",";
      message += Utilities::int_to_string(proposed_n_q_points_1d);
      message += "," + Utilities::int_to_string(n_components);
      message += ",Number>(data, ";
      message += Utilities::int_to_string(fe_no) + ", ";
      message += Utilities::int_to_string(this->quad_no) + ")?\n";
      std::string correct_pos;
      if (this->data->fe_degree != fe_degree)
        correct_pos = " ^";
      else
        correct_pos = "  ";
      if (proposed_n_q_points_1d != n_q_points_1d)
        correct_pos += " ^\n";
      else
        correct_pos += "  \n";
      message += "                                 " + correct_pos;

      Assert (fe_degree == this->data->fe_degree &&
              n_q_points == this->data->n_q_points,
              ExcMessage(message));
    }
  if (fe_no != numbers::invalid_unsigned_int)
    {
      AssertDimension (n_q_points,
                       this->mapping_info->mapping_data_gen[this->quad_no].
                       n_q_points[this->active_quad_index]);
      AssertDimension (this->data->dofs_per_cell * this->n_fe_components,
                       this->dof_info->dofs_per_cell[this->active_fe_index]);
    }
#endif
}



template <int dim, int fe_degree,  int n_q_points_1d, int n_components_,
          typename Number>
inline
Point<dim,VectorizedArray<Number> >
FEEvaluation<dim,fe_degree,n_q_points_1d,n_components_,Number>
::quadrature_point (const unsigned int q) const
{
  Assert (this->mapping_info->quadrature_points_initialized == true,
          ExcNotInitialized());
  AssertIndexRange (q, n_q_points);

  // Cartesian mesh: not all quadrature points are stored, only the
  // diagonal. Hence, need to find the tensor product index and retrieve the
  // value from that
  if (this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      Point<dim,VectorizedArray<Number> > point (false);
      switch (dim)
        {
        case 1:
          return this->quadrature_points[q];
        case 2:
          point[0] = this->quadrature_points[q%n_q_points_1d][0];
          point[1] = this->quadrature_points[q/n_q_points_1d][1];
          return point;
        case 3:
          point[0] = this->quadrature_points[q%n_q_points_1d][0];
          point[1] = this->quadrature_points[(q/n_q_points_1d)%n_q_points_1d][1];
          point[2] = this->quadrature_points[q/(n_q_points_1d*n_q_points_1d)][2];
          return point;
        default:
          Assert (false, ExcNotImplemented());
          return point;
        }
    }
  // all other cases: just return the respective data as it is fully stored
  else
    return this->quadrature_points[q];
}



template <int dim, int fe_degree,  int n_q_points_1d, int n_components_,
          typename Number>
inline
void
FEEvaluation<dim,fe_degree,n_q_points_1d,n_components_,Number>
::evaluate (const bool evaluate_val,
            const bool evaluate_grad,
            const bool evaluate_lapl)
{
  Assert (this->dof_values_initialized == true,
          internal::ExcAccessToUninitializedField());

  // Select algorithm matching the element type at run time (the function
  // pointer is easy to predict, so negligible in cost)
  evaluate_funct (*this->data, &this->values_dofs[0],
                  this->values_quad, this->gradients_quad, this->hessians_quad,
                  evaluate_val, evaluate_grad, evaluate_lapl);

#ifdef DEBUG
  if (evaluate_val == true)
    this->values_quad_initialized = true;
  if (evaluate_grad == true)
    this->gradients_quad_initialized = true;
  if (evaluate_lapl == true)
    this->hessians_quad_initialized  = true;
#endif
}



template <int dim, int fe_degree,  int n_q_points_1d, int n_components_,
          typename Number>
inline
void
FEEvaluation<dim,fe_degree,n_q_points_1d,n_components_,Number>
::integrate (bool integrate_val,bool integrate_grad)
{
  if (integrate_val == true)
    Assert (this->values_quad_submitted == true,
            internal::ExcAccessToUninitializedField());
  if (integrate_grad == true)
    Assert (this->gradients_quad_submitted == true,
            internal::ExcAccessToUninitializedField());

  // Select algorithm matching the element type at run time (the function
  // pointer is easy to predict, so negligible in cost)
  integrate_funct (*this->data, this->values_dofs, this->values_quad,
                   this->gradients_quad, integrate_val, integrate_grad);

#ifdef DEBUG
  this->dof_values_initialized = true;
#endif
}



#endif  // ifndef DOXYGEN


DEAL_II_NAMESPACE_CLOSE

#endif
