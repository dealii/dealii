// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2017 by the deal.II authors
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


#ifndef dealii_matrix_free_fe_evaluation_h
#define dealii_matrix_free_fe_evaluation_h


#include <deal.II/base/array_view.h>
#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/matrix_free/mapping_data_on_the_fly.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/shape_info.h>
#include <deal.II/matrix_free/evaluation_kernels.h>
#include <deal.II/matrix_free/tensor_product_kernels.h>
#include <deal.II/matrix_free/evaluation_selector.h>

#include <deal.II/lac/vector_operation.h>

DEAL_II_NAMESPACE_OPEN



// forward declarations
namespace LinearAlgebra
{
  namespace distributed
  {
    template <typename> class Vector;
  }
}
namespace internal
{
  DeclException0 (ExcAccessToUninitializedField);
}

template <int dim, int fe_degree, int n_q_points_1d = fe_degree+1,
          int n_components_ = 1, typename Number = double > class FEEvaluation;


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
 * PDEs. If the same operation is applied to several components of a PDE (e.g.
 * a vector Laplace equation), they can be applied simultaneously with one
 * call (and often more efficiently)
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
  static constexpr unsigned int dimension     = dim;
  static constexpr unsigned int n_components  = n_components_;

  /**
   * @name 1: General operations
   */
  //@{
  /**
   * Destructor.
   */
  ~FEEvaluationBase();

  /**
   * Initializes the operation pointer to the current cell. Unlike the reinit
   * functions taking a cell iterator as argument below and the
   * FEValues::reinit() methods, where the information related to a particular
   * cell is generated in the reinit call, this function is very cheap since
   * all data is pre-computed in @p matrix_free, and only a few indices have
   * to be set appropriately.
   */
  void reinit (const unsigned int cell);

  /**
   * Initialize the data to the current cell using a TriaIterator object as
   * usual in FEValues. The argument is either of type
   * DoFHandler::active_cell_iterator or DoFHandler::level_cell_iterator. This
   * option is only available if the FEEvaluation object was created with a
   * finite element, quadrature formula and correct update flags and
   * <b>without</b> a MatrixFree object. This initialization method loses the
   * ability to use vectorization, see also the description of the
   * FEEvaluation class. When this reinit method is used, FEEvaluation can
   * also read from vectors (but less efficient than with data coming from
   * MatrixFree).
   */
  template <typename DoFHandlerType, bool level_dof_access>
  void reinit (const TriaIterator<DoFCellAccessor<DoFHandlerType,level_dof_access> > &cell);

  /**
   * Initialize the data to the current cell using a TriaIterator object as
   * usual in FEValues. This option is only available if the FEEvaluation
   * object was created with a finite element, quadrature formula and correct
   * update flags and <b>without</b> a MatrixFree object. This initialization
   * method loses the ability to use vectorization, see also the description
   * of the FEEvaluation class. When this reinit method is used, FEEvaluation
   * can <b>not</b> read from vectors because no DoFHandler information is
   * available.
   */
  void reinit (const typename Triangulation<dim>::cell_iterator &cell);

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
   * Return the type of the cell the @p reinit function has been called for.
   * Valid values are @p cartesian for Cartesian cells (which allows for
   * considerable data compression), @p affine for cells with affine mappings,
   * and @p general for general cells without any compressed storage applied.
   */
  internal::MatrixFreeFunctions::CellType get_cell_type() const;

  /**
   * Return a reference to the ShapeInfo object currently in use.
   */
  const internal::MatrixFreeFunctions::ShapeInfo<VectorizedArray<Number>> &
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
   * also in constrained positions by calling ConstraintMatrix::distribute.
   * When accessing vector entries during the solution of linear systems, the
   * temporary solution should always have homogeneous constraints and this
   * method is the correct one.
   *
   * If this class was constructed without a MatrixFree object and the
   * information is acquired on the fly through a
   * DoFHandler<dim>::cell_iterator, only one single cell is used by this
   * class and this function extracts the values of the underlying components
   * on the given cell. This call is slower than the ones done through a
   * MatrixFree object and lead to a structure that does not effectively use
   * vectorization in the evaluate routines based on these values (instead,
   * VectorizedArray::n_array_elements same copies are worked on).
   *
   * If the given vector template class is a block vector (determined through
   * the template function 'IsBlockVector<VectorType>::value', which checks
   * for vectors derived from dealii::BlockVectorBase) or an
   * std::vector<VectorType> or std::vector<VectorType *>, this function reads
   * @p n_components blocks from the block vector starting at the index
   * @p first_index. For non-block vectors, @p first_index is ignored.
   */
  template <typename VectorType>
  void read_dof_values (const VectorType  &src,
                        const unsigned int first_index = 0);

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
   * If this class was constructed without a MatrixFree object and the
   * information is acquired on the fly through a
   * DoFHandler<dim>::cell_iterator, only one single cell is used by this
   * class and this function extracts the values of the underlying components
   * on the given cell. This call is slower than the ones done through a
   * MatrixFree object and lead to a structure that does not effectively use
   * vectorization in the evaluate routines based on these values (instead,
   * VectorizedArray::n_array_elements same copies are worked on).
   *
   * If the given vector template class is a block vector (determined through
   * the template function 'IsBlockVector<VectorType>::value', which checks
   * for vectors derived from dealii::BlockVectorBase) or an
   * std::vector<VectorType> or std::vector<VectorType *>, this function reads
   * @p n_components blocks from the block vector starting at the index
   * @p first_index. For non-block vectors, @p first_index is ignored.
   */
  template <typename VectorType>
  void read_dof_values_plain (const VectorType  &src,
                              const unsigned int first_index = 0);

  /**
   * Takes the values stored internally on dof values of the current cell and
   * sums them into the vector @p dst. The function also applies constraints
   * during the write operation. The functionality is hence similar to the
   * function ConstraintMatrix::distribute_local_to_global. If vectorization
   * is enabled, the DoF values for several cells are used.
   *
   * If this class was constructed without a MatrixFree object and the
   * information is acquired on the fly through a
   * DoFHandler<dim>::cell_iterator, only one single cell is used by this
   * class and this function extracts the values of the underlying components
   * on the given cell. This call is slower than the ones done through a
   * MatrixFree object and lead to a structure that does not effectively use
   * vectorization in the evaluate routines based on these values (instead,
   * VectorizedArray::n_array_elements same copies are worked on).
   *
   * If the given vector template class is a block vector (determined through
   * the template function 'IsBlockVector<VectorType>::value', which checks
   * for vectors derived from dealii::BlockVectorBase) or an
   * std::vector<VectorType> or std::vector<VectorType *>, this function
   * writes to @p n_components blocks of the block vector starting at the
   * index @p first_index. For non-block vectors, @p first_index is ignored.
   */
  template <typename VectorType>
  void distribute_local_to_global (VectorType        &dst,
                                   const unsigned int first_index = 0) const;

  /**
   * Takes the values stored internally on dof values of the current cell and
   * writes them into the vector @p dst. The function skips the degrees of
   * freedom which are constrained. As opposed to the
   * distribute_local_to_global method, the old values at the position given
   * by the current cell are overwritten. Thus, if a degree of freedom is
   * associated to more than one cell (as usual in continuous finite
   * elements), the values will be overwritten and only the value written last
   * is retained.
   *
   * If this class was constructed without a MatrixFree object and the
   * information is acquired on the fly through a
   * DoFHandler<dim>::cell_iterator, only one single cell is used by this
   * class and this function extracts the values of the underlying components
   * on the given cell. This call is slower than the ones done through a
   * MatrixFree object and lead to a structure that does not effectively use
   * vectorization in the evaluate routines based on these values (instead,
   * VectorizedArray::n_array_elements same copies are worked on).
   *
   * If the given vector template class is a block vector (determined through
   * the template function 'IsBlockVector<VectorType>::value', which checks
   * for vectors derived from dealii::BlockVectorBase) or an
   * std::vector<VectorType> or std::vector<VectorType *>, this function
   * writes to @p n_components blocks of the block vector starting at the
   * index @p first_index. For non-block vectors, @p first_index is ignored.
   */
  template <typename VectorType>
  void set_dof_values (VectorType        &dst,
                       const unsigned int first_index = 0) const;

  //@}

  /**
   * @name 3: Data access
   */
  //@{
  /**
   * Return the value stored for the local degree of freedom with index @p
   * dof. If the object is vector-valued, a vector-valued return argument is
   * given. Thus, the argument @p dof can at most run until @p
   * dofs_per_component rather than @p dofs_per_cell since the different
   * components of a vector-valued FE are return together. Note that when
   * vectorization is enabled, values from several cells are grouped
   * together. If @p set_dof_values was called last, the value corresponds to
   * the one set there. If @p integrate was called last, it instead
   * corresponds to the value of the integrated function with the test
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
   * Return the value of a finite element function at quadrature point number
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
   * Return the gradient of a finite element function at quadrature point
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
   * containing the values on quadrature points with component @p q_point.
   * Access to the same field as through @p get_gradient. If applied before
   * the function @p integrate(...,true) is called, this specifies what is
   * tested by all basis function gradients on the current cell and integrated
   * over.
   *
   * Note that the derived class FEEvaluationAccess overloads this operation
   * with specializations for the scalar case (n_components == 1) and for the
   * vector-valued case (n_components == dim).
   */
  void submit_gradient(const gradient_type grad_in,
                       const unsigned int  q_point);

  /**
   * Return the Hessian of a finite element function at quadrature point
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
   * Return the diagonal of the Hessian of a finite element function at
   * quadrature point number @p q_point after a call to @p evaluate(...,true).
   *
   * Note that the derived class FEEvaluationAccess overloads this operation
   * with specializations for the scalar case (n_components == 1) and for the
   * vector-valued case (n_components == dim).
   */
  gradient_type get_hessian_diagonal (const unsigned int q_point) const;

  /**
   * Return the Laplacian (i.e., the trace of the Hessian) of a finite
   * element function at quadrature point number @p q_point after a call to @p
   * evaluate(...,true). Compared to the case when computing the full Hessian,
   * some operations can be saved when only the Laplacian is requested.
   *
   * Note that the derived class FEEvaluationAccess overloads this operation
   * with specializations for the scalar case (n_components == 1) and for the
   * vector-valued case (n_components == dim).
   */
  value_type get_laplacian (const unsigned int q_point) const;

#ifdef DOXYGEN
  // doxygen does not anyhow mention functions coming from partial template
  // specialization of the base class, in this case FEEvaluationAccess<dim,dim>.
  // For now, hack-in those functions manually only to fix documentation:

  /** @copydoc FEEvaluationAccess<dim,dim,Number>::get_divergence()
   * @note Only available for n_components_==dim.
   */
  VectorizedArray<Number> get_divergence (const unsigned int q_point) const;

  /** @copydoc FEEvaluationAccess<dim,dim,Number>::get_symmetric_gradient()
   * @note Only available for n_components_==dim.
   */
  SymmetricTensor<2, dim, VectorizedArray<Number> > get_symmetric_gradient (const unsigned int q_point) const;

  /** @copydoc FEEvaluationAccess<dim,dim,Number>::get_curl()
   * @note Only available for n_components_==dim.
   */
  Tensor<1,(dim==2?1:dim), VectorizedArray<Number> > get_curl (const unsigned int q_point) const;

  /** @copydoc FEEvaluationAccess<dim,dim,Number>::submit_divergence()
   * @note Only available for n_components_==dim.
   */
  void submit_divergence (const VectorizedArray<Number> div_in, const unsigned int q_point);

  /** @copydoc FEEvaluationAccess<dim,dim,Number>::submit_symmetric_gradient()
   * @note Only available for n_components_==dim.
   */
  void submit_symmetric_gradient (const SymmetricTensor<2, dim, VectorizedArray<Number> > grad_in, const unsigned int q_point);

  /** @copydoc FEEvaluationAccess<dim,dim,Number>::submit_curl()
   * @note Only available for n_components_==dim.
   */
  void submit_curl (const Tensor<1, dim==2?1:dim, VectorizedArray<Number> > curl_in, const unsigned int q_point);

#endif

  /**
   * Takes values on quadrature points, multiplies by the Jacobian determinant
   * and quadrature weights (JxW) and sums the values for all quadrature
   * points on the cell. The result is a scalar, representing the integral
   * over the function over the cell. If a vector-element is used, the
   * resulting components are still separated. Moreover, if vectorization is
   * enabled, the integral values of several cells are contained in the slots
   * of the returned VectorizedArray field.
   *
   * @note In case the FEEvaluation object is initialized with a batch of
   * cells where not all lanes in the SIMD vector VectorizedArray are
   * representing actual data, this method performs computations on dummy data
   * (that is copied from the last valid lane) that will not make sense. Thus,
   * the user needs to make sure that it is not used in any computation
   * explicitly, like when summing the results of several cells.
   */
  value_type integrate_value () const;

  /**
   * Return the determinant of the Jacobian from the unit to the real cell
   * times the quadrature weight.
   */
  VectorizedArray<Number> JxW(const unsigned int q_point) const;

  //@}

  /**
   * @name 4: Access to internal data
   */
  //@{
  /**
   * Return a read-only pointer to the first field of the dof values. This is
   * the data field the read_dof_values() functions write into. First come the
   * the dof values for the first component, then all values for the second
   * component, and so on. This is related to the internal data structures
   * used in this class. In general, it is safer to use the get_dof_value()
   * function instead.
   */
  const VectorizedArray<Number> *begin_dof_values () const;

  /**
   * Return a read and write pointer to the first field of the dof values.
   * This is the data field the read_dof_values() functions write into. First
   * come the dof values for the first component, then all values for the
   * second component, and so on. This is related to the internal data
   * structures used in this class. In general, it is safer to use the
   * get_dof_value() function instead.
   */
  VectorizedArray<Number> *begin_dof_values ();

  /**
   * Return a read-only pointer to the first field of function values on
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
   * Return a read and write pointer to the first field of function values on
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
   * Return a read-only pointer to the first field of function gradients on
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
   * Return a read and write pointer to the first field of function gradients
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
   * Return a read-only pointer to the first field of function hessians on
   * quadrature points. First comes the xx-component of the hessian for the
   * first component on all quadrature points, then the yy-component, zz-
   * component in (3D), then the xy-component, and so on. Next comes the xx-
   * component of the second component, and so on. This is related to the
   * internal data structures used in this class. The raw data after a call to
   * @p evaluate only contains unit cell operations, so possible
   * transformations, quadrature weights etc. must be applied manually. In
   * general, it is safer to use the get_laplacian() or get_hessian()
   * functions instead, which does all the transformation internally.
   */
  const VectorizedArray<Number> *begin_hessians () const;

  /**
   * Return a read and write pointer to the first field of function hessians
   * on quadrature points. First comes the xx-component of the hessian for the
   * first component on all quadrature points, then the yy-component, zz-
   * component in (3D), then the xy-component, and so on. Next comes the xx-
   * component of the second component, and so on. This is related to the
   * internal data structures used in this class. The raw data after a call to
   * @p evaluate only contains unit cell operations, so possible
   * transformations, quadrature weights etc. must be applied manually. In
   * general, it is safer to use the get_laplacian() or get_hessian()
   * functions instead, which does all the transformation internally.
   */
  VectorizedArray<Number> *begin_hessians ();

  /**
   * Return the numbering of local degrees of freedom within the evaluation
   * routines of FEEvaluation in terms of the standard numbering on finite
   * elements.
   */
  const std::vector<unsigned int> &
  get_internal_dof_numbering() const;

  /**
   * Return an ArrayView to internal memory for temporary use. Note that some
   * of this memory is overwritten during evaluate() and integrate() calls so
   * do not assume it to be stable over those calls. The maximum size you can
   * write into is 3*dofs_per_cell+2*n_q_points.
   */
  ArrayView<VectorizedArray<Number> >
  get_scratch_data() const;

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
   * FEValues. The arguments are similar to the ones passed to the constructor
   * of FEValues, with the notable difference that FEEvaluation expects a one-
   * dimensional quadrature formula, Quadrature<1>, instead of a @p dim
   * dimensional one. The finite element can be both scalar or vector valued,
   * but this method always only selects a scalar base element at a time (with
   * @p n_components copies as specified by the class template argument). For
   * vector-valued elements, the optional argument @p first_selected_component
   * allows to specify the index of the base element to be used for
   * evaluation. Note that the internal data structures always assume that the
   * base element is primitive, non-primitive are not supported currently.
   *
   * As known from FEValues, a call to the reinit method with a
   * Triangulation::cell_iterator is necessary to make the geometry and
   * degrees of freedom of the current class known. If the iterator includes
   * DoFHandler information (i.e., it is a DoFHandler::cell_iterator or
   * similar), the initialization allows to also read from or write to vectors
   * in the standard way for DoFHandler::active_cell_iterator types for one
   * cell at a time. However, this approach is much slower than the path with
   * MatrixFree with MPI since index translation has to be done. As only one
   * cell at a time is used, this method does not vectorize over several
   * elements (which is most efficient for vector operations), but only
   * possibly within the element if the evaluate/integrate routines are
   * combined inside user code (e.g. for computing cell matrices).
   *
   * The optional FEEvaluationBase object allows several FEEvaluation objects
   * to share the geometry evaluation, i.e., the underlying mapping and
   * quadrature points do only need to be evaluated once. This only works if
   * the quadrature formulas are the same. Otherwise, a new evaluation object
   * is created. Make sure to not pass an optional object around when you
   * intend to use the FEEvaluation object in %parallel with another one
   * because otherwise the intended sharing may create race conditions.
   */
  template <int n_components_other>
  FEEvaluationBase (const Mapping<dim>       &mapping,
                    const FiniteElement<dim> &fe,
                    const Quadrature<1>      &quadrature,
                    const UpdateFlags         update_flags,
                    const unsigned int        first_selected_component,
                    const FEEvaluationBase<dim,n_components_other,Number> *other);

  /**
   * Copy constructor. If FEEvaluationBase was constructed from a mapping, fe,
   * quadrature, and update flags, the underlying geometry evaluation based on
   * FEValues will be deep-copied in order to allow for using in parallel with
   * threads.
   */
  FEEvaluationBase (const FEEvaluationBase &other);

  /**
   * Copy assignment operator. If FEEvaluationBase was constructed from a
   * mapping, fe, quadrature, and update flags, the underlying geometry
   * evaluation based on FEValues will be deep-copied in order to allow for
   * using in parallel with threads.
   */
  FEEvaluationBase &operator = (const FEEvaluationBase &other);

  /**
   * A unified function to read from and write into vectors based on the given
   * template operation. It can perform the operation for @p read_dof_values,
   * @p distribute_local_to_global, and @p set_dof_values. It performs the
   * operation for several vectors at a time.
   */
  template <typename VectorType, typename VectorOperation>
  void read_write_operation (const VectorOperation &operation,
                             VectorType            *vectors[]) const;

  /**
   * For a collection of several vector @p src, read out the values on the
   * degrees of freedom of the current cell for @p n_components (template
   * argument), and store them internally. Similar functionality as the
   * function DoFAccessor::read_dof_values. Note that if vectorization is
   * enabled, the DoF values for several cells are set.
   */
  template <typename VectorType>
  void read_dof_values_plain (const VectorType *src_data[]);

  /**
   * This is the general array for all data fields.
   */
  AlignedVector<VectorizedArray<Number> > *scratch_data_array;

  /**
   * This is the user-visible part of scratch_data_array, only showing the
   * last part of scratch_data_array. The first part is consumed by
   * values_dofs, values_quad, etc.
   */
  VectorizedArray<Number> *scratch_data;

  /**
   * This field stores the values for local degrees of freedom (e.g. after
   * reading out from a vector but before applying unit cell transformations
   * or before distributing them into a result vector). The methods
   * get_dof_value() and submit_dof_value() read from or write to this field.
   *
   * The values of this array are stored in the start section of
   * @p scratch_data_array. Due to its access as a thread local memory, the
   * memory can get reused between different calls. As opposed to requesting
   * memory on the stack, this approach allows for very large polynomial
   * degrees.
   */
  VectorizedArray<Number> *values_dofs[n_components];

  /**
   * This field stores the values of the finite element function on quadrature
   * points after applying unit cell transformations or before integrating.
   * The methods get_value() and submit_value() access this field.
   *
   * The values of this array are stored in the start section of
   * @p scratch_data_array. Due to its access as a thread local memory, the
   * memory can get reused between different calls. As opposed to requesting
   * memory on the stack, this approach allows for very large polynomial
   * degrees.
   */
  VectorizedArray<Number> *values_quad[n_components];

  /**
   * This field stores the gradients of the finite element function on
   * quadrature points after applying unit cell transformations or before
   * integrating. The methods get_gradient() and submit_gradient() (as well as
   * some specializations like get_symmetric_gradient() or get_divergence())
   * access this field.
   *
   * The values of this array are stored in the start section of
   * @p scratch_data_array. Due to its access as a thread local memory, the
   * memory can get reused between different calls. As opposed to requesting
   * memory on the stack, this approach allows for very large polynomial
   * degrees.
   */
  VectorizedArray<Number> *gradients_quad[n_components][dim];

  /**
   * This field stores the Hessians of the finite element function on
   * quadrature points after applying unit cell transformations. The methods
   * get_hessian(), get_laplacian(), get_hessian_diagonal() access this field.
   *
   * The values of this array are stored in the start section of
   * @p scratch_data_array. Due to its access as a thread local memory, the
   * memory can get reused between different calls. As opposed to requesting
   * memory on the stack, this approach allows for very large polynomial
   * degrees.
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
   * Stores a pointer to the underlying DoF indices and constraint description
   * for the component specified at construction. Also contained in
   * matrix_info, but it simplifies code if we store a reference to it.
   */
  const internal::MatrixFreeFunctions::DoFInfo *dof_info;

  /**
   * Stores a pointer to the underlying transformation data from unit to real
   * cells for the given quadrature formula specified at construction. Also
   * contained in matrix_info, but it simplifies code if we store a reference
   * to it.
   */
  const internal::MatrixFreeFunctions::MappingInfo<dim,Number> *mapping_info;

  /**
   * Stores a pointer to the unit cell shape data, i.e., values, gradients and
   * Hessians in 1D at the quadrature points that constitute the tensor
   * product. Also contained in matrix_info, but it simplifies code if we
   * store a reference to it.
   */
  const internal::MatrixFreeFunctions::ShapeInfo<VectorizedArray<Number>> *data;

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
   * present cell. Only set to a useful value if on a general cell with non-
   * constant Jacobian.
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
   * been submitted for integration before the integration is actually stared.
   * Used to control exceptions when uninitialized data is used.
   */
  bool gradients_quad_submitted;

  /**
   * Geometry data that can be generated FEValues on the fly with the
   * respective constructor.
   */
  std::shared_ptr<internal::MatrixFreeFunctions::MappingDataOnTheFly<dim,Number> > mapped_geometry;

  /**
   * For use with on-the-fly evaluation, provide a data structure to store the
   * global dof indices on the current cell from a reinit call.
   */
  std::vector<types::global_dof_index> old_style_dof_indices;

  /**
   * For a FiniteElement with more than one finite element, select at which
   * component this data structure should start.
   */
  const unsigned int first_selected_component;

  /**
   * A temporary data structure necessary to read degrees of freedom when no
   * MatrixFree object was given at initialization.
   */
  mutable std::vector<types::global_dof_index> local_dof_indices;

private:
  /**
   * Sets the pointers for values, gradients, hessians to the central
   * scratch_data_array.
   */
  void set_data_pointers();

  /**
   * Make other FEEvaluationBase as well as FEEvaluation objects friends.
   */
  template <int, int, typename> friend class FEEvaluationBase;
  template <int, int, int, int, typename> friend class FEEvaluation;
};



/**
 * This class provides access to the data fields of the FEEvaluation classes.
 * Generic access is achieved through the base class, and specializations for
 * scalar and vector-valued elements are defined separately.
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
  static constexpr unsigned int dimension     = dim;
  static constexpr unsigned int n_components  = n_components_;
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
   * Constructor with reduced functionality for similar usage of FEEvaluation
   * as FEValues, including matrix assembly.
   */
  template <int n_components_other>
  FEEvaluationAccess (const Mapping<dim>       &mapping,
                      const FiniteElement<dim> &fe,
                      const Quadrature<1>      &quadrature,
                      const UpdateFlags         update_flags,
                      const unsigned int        first_selected_component,
                      const FEEvaluationBase<dim,n_components_other,Number> *other);

  /**
   * Copy constructor
   */
  FEEvaluationAccess (const FEEvaluationAccess &other);

  /**
   * Copy assignment operator
   */
  FEEvaluationAccess &operator= (const FEEvaluationAccess &other);
};




/**
 * This class provides access to the data fields of the FEEvaluation classes.
 * Partial specialization for scalar fields that defines access with simple
 * data fields, i.e., scalars for the values and Tensor<1,dim> for the
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
  static constexpr unsigned int dimension          = dim;
  typedef FEEvaluationBase<dim,1,Number>         BaseClass;

  /** @copydoc FEEvaluationBase<dim,1,Number>::get_dof_value()
   */
  value_type get_dof_value (const unsigned int dof) const;

  /** @copydoc FEEvaluationBase<dim,1,Number>::submit_dof_value()
   */
  void submit_dof_value (const value_type   val_in,
                         const unsigned int dof);

  /** @copydoc FEEvaluationBase<dim,1,Number>::get_value()
   */
  value_type get_value (const unsigned int q_point) const;

  /** @copydoc FEEvaluationBase<dim,1,Number>::submit_value()
   */
  void submit_value (const value_type   val_in,
                     const unsigned int q_point);

  /** @copydoc FEEvaluationBase<dim,1,Number>::get_gradient()
   */
  gradient_type get_gradient (const unsigned int q_point) const;

  /** @copydoc FEEvaluationBase<dim,1,Number>::submit_gradient()
   */
  void submit_gradient(const gradient_type grad_in,
                       const unsigned int  q_point);

  /** @copydoc FEEvaluationBase<dim,1,Number>::get_hessian()
   */
  Tensor<2,dim,VectorizedArray<Number> >
  get_hessian (unsigned int q_point) const;

  /** @copydoc FEEvaluationBase<dim,1,Number>::get_hessian_diagonal()
   */
  gradient_type get_hessian_diagonal (const unsigned int q_point) const;

  /** @copydoc FEEvaluationBase<dim,1,Number>::get_laplacian()
   */
  value_type get_laplacian (const unsigned int q_point) const;

  /** @copydoc FEEvaluationBase<dim,1,Number>::integrate_value()
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
   * Constructor with reduced functionality for similar usage of FEEvaluation
   * as FEValues, including matrix assembly.
   */
  template <int n_components_other>
  FEEvaluationAccess (const Mapping<dim>       &mapping,
                      const FiniteElement<dim> &fe,
                      const Quadrature<1>      &quadrature,
                      const UpdateFlags         update_flags,
                      const unsigned int        first_selected_component,
                      const FEEvaluationBase<dim,n_components_other,Number> *other);

  /**
   * Copy constructor
   */
  FEEvaluationAccess (const FEEvaluationAccess &other);

  /**
   * Copy assignment operator
   */
  FEEvaluationAccess &operator= (const FEEvaluationAccess &other);
};



/**
 * This class provides access to the data fields of the FEEvaluation classes.
 * Partial specialization for fields with as many components as the underlying
 * space dimension, i.e., values are of type Tensor<1,dim> and gradients of
 * type Tensor<2,dim>. Provides some additional functions for access, like the
 * symmetric gradient and divergence.
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
  static constexpr unsigned int dimension     = dim;
  static constexpr unsigned int n_components  = dim;
  typedef FEEvaluationBase<dim,dim,Number> BaseClass;

  /** @copydoc FEEvaluationBase<dim,dim,Number>::get_gradient()
   */
  gradient_type get_gradient (const unsigned int q_point) const;

  /**
   * Return the divergence of a vector-valued finite element at quadrature
   * point number @p q_point after a call to @p evaluate(...,true,...).
   */
  VectorizedArray<Number> get_divergence (const unsigned int q_point) const;

  /**
   * Return the symmetric gradient of a vector-valued finite element at
   * quadrature point number @p q_point after a call to @p
   * evaluate(...,true,...). It corresponds to <tt>0.5
   * (grad+grad<sup>T</sup>)</tt>.
   */
  SymmetricTensor<2,dim,VectorizedArray<Number> >
  get_symmetric_gradient (const unsigned int q_point) const;

  /**
   * Return the curl of the vector field, $\nabla \times v$ after a call to @p
   * evaluate(...,true,...).
   */
  Tensor<1,(dim==2?1:dim),VectorizedArray<Number> >
  get_curl (const unsigned int q_point) const;

  /** @copydoc FEEvaluationBase<dim,dim,Number>::get_hessian()
   */
  Tensor<3,dim,VectorizedArray<Number> >
  get_hessian (const unsigned int q_point) const;

  /** @copydoc FEEvaluationBase<dim,dim,Number>::get_hessian_diagonal()
   */
  gradient_type get_hessian_diagonal (const unsigned int q_point) const;

  /** @copydoc FEEvaluationBase<dim,dim,Number>::submit_gradient()
   */
  void submit_gradient(const gradient_type grad_in,
                       const unsigned int  q_point);

  /**
   * Write a contribution that is tested by the gradient to the field
   * containing the values on quadrature points with component @p q_point.
   * This function is an alternative to the other submit_gradient function
   * when using a system of fixed number of equations which happens to
   * coincide with the dimension for some dimensions, but not all. To allow
   * for dimension-independent programming, this function can be used instead.
   */
  void submit_gradient(const Tensor<1,dim,Tensor<1,dim,VectorizedArray<Number> > > grad_in,
                       const unsigned int q_point);

  /**
   * Write a contribution that is tested by the divergence to the field
   * containing the values on quadrature points with component @p q_point.
   * Access to the same field as through @p get_gradient. If applied before
   * the function @p integrate(...,true) is called, this specifies what is
   * tested by all basis function gradients on the current cell and integrated
   * over.
   */
  void submit_divergence (const VectorizedArray<Number> div_in,
                          const unsigned int q_point);

  /**
   * Write a contribution that is tested by the symmetric gradient to the field
   * containing the values on quadrature points with component @p q_point.
   * Access to the same field as through @p get_symmetric_gradient. If applied before
   * the function @p integrate(...,true) is called, this specifies the
   * symmetric gradient which is tested by all basis function symmetric gradients on the current
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
   * Constructor with reduced functionality for similar usage of FEEvaluation
   * as FEValues, including matrix assembly.
   */
  template <int n_components_other>
  FEEvaluationAccess (const Mapping<dim>       &mapping,
                      const FiniteElement<dim> &fe,
                      const Quadrature<1>      &quadrature,
                      const UpdateFlags         update_flags,
                      const unsigned int        first_selected_component,
                      const FEEvaluationBase<dim,n_components_other,Number> *other);

  /**
   * Copy constructor
   */
  FEEvaluationAccess (const FEEvaluationAccess &other);

  /**
   * Copy assignment operator
   */
  FEEvaluationAccess &operator= (const FEEvaluationAccess &other);
};


/**
 * This class provides access to the data fields of the FEEvaluation classes.
 * Partial specialization for scalar fields in 1d that defines access with
 * simple data fields, i.e., scalars for the values and Tensor<1,1> for the
 * gradients.
 *
 * @author Katharina Kormann and Martin Kronbichler, 2010, 2011, Shiva
 * Rudraraju, 2014
 */
template <typename Number>
class FEEvaluationAccess<1,1,Number> : public FEEvaluationBase<1,1,Number>
{
public:
  typedef Number                                 number_type;
  typedef VectorizedArray<Number>                value_type;
  typedef Tensor<1,1,VectorizedArray<Number> >   gradient_type;
  static constexpr unsigned int dimension          = 1;
  typedef FEEvaluationBase<1,1,Number>           BaseClass;

  /** @copydoc FEEvaluationBase<1,1,Number>::get_dof_value()
   */
  value_type get_dof_value (const unsigned int dof) const;

  /** @copydoc FEEvaluationBase<1,1,Number>::submit_dof_value()
   */
  void submit_dof_value (const value_type   val_in,
                         const unsigned int dof);

  /** @copydoc FEEvaluationBase<1,1,Number>::get_value()
   */
  value_type get_value (const unsigned int q_point) const;

  /** @copydoc FEEvaluationBase<1,1,Number>::submit_value()
   */
  void submit_value (const value_type   val_in,
                     const unsigned int q_point);

  /** @copydoc FEEvaluationBase<1,1,Number>::get_gradient()
   */
  gradient_type get_gradient (const unsigned int q_point) const;

  /** @copydoc FEEvaluationBase<1,1,Number>::submit_gradient()
   */
  void submit_gradient(const gradient_type grad_in,
                       const unsigned int  q_point);

  /** @copydoc FEEvaluationBase<1,1,Number>::get_hessian()
   */
  Tensor<2,1,VectorizedArray<Number> >
  get_hessian (unsigned int q_point) const;

  /** @copydoc FEEvaluationBase<1,1,Number>::get_hessian_diagonal()
   */
  gradient_type get_hessian_diagonal (const unsigned int q_point) const;

  /** @copydoc FEEvaluationBase<1,1,Number>::get_laplacian()
   */
  value_type get_laplacian (const unsigned int q_point) const;

  /** @copydoc FEEvaluationBase<1,1,Number>::integrate_value()
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
  FEEvaluationAccess (const MatrixFree<1,Number> &matrix_free,
                      const unsigned int          fe_no,
                      const unsigned int          quad_no,
                      const unsigned int          fe_degree,
                      const unsigned int          n_q_points);

  /**
   * Constructor with reduced functionality for similar usage of FEEvaluation
   * as FEValues, including matrix assembly.
   */
  template <int n_components_other>
  FEEvaluationAccess (const Mapping<1>       &mapping,
                      const FiniteElement<1> &fe,
                      const Quadrature<1>    &quadrature,
                      const UpdateFlags       update_flags,
                      const unsigned int      first_selected_component,
                      const FEEvaluationBase<1,n_components_other,Number> *other);

  /**
   * Copy constructor
   */
  FEEvaluationAccess (const FEEvaluationAccess &other);

  /**
   * Copy assignment operator
   */
  FEEvaluationAccess &operator= (const FEEvaluationAccess &other);
};



/**
 * The class that provides all functions necessary to evaluate functions at
 * quadrature points and cell integrations. In functionality, this class is
 * similar to FEValues, however, it includes a lot of specialized functions
 * that make it much faster (between 5 and 500, depending on the polynomial
 * degree).
 *
 * <h3>Usage and initialization</h3>
 *
 * <h4>Fast usage in combination with MatrixFree</h4>
 *
 * The first and foremost way of usage is to initialize this class from a
 * MatrixFree object that caches everything related to the degrees of freedom
 * and the mapping information. This way, it is possible to use vectorization
 * for applying a differential operator for several cells at once.
 *
 * The capabilities of FEEvaluation span a large spectrum of integration tasks
 * for weak forms. In general, there are two classes of tasks that get
 * done. One is the @p evaluate path that interpolates from a solution vector
 * to quadrature points:
 *
 * @code
 * FEEvaluation<dim,fe_degree> phi(matrix_free);
 * for (unsigned int cell_index = cell_range.first;
 *      cell_index < cell_range.second; ++cell_index)
 *   {
 *     phi.reinit(cell_index);
 *     phi.read_dof_values(vector);
 *     phi.evaluate(true, false);   // interpolate values, but not gradients
 *     for (unsigned int q_index=0; q_index<phi.n_q_points; ++q_index)
 *       {
 *         VectorizedArray<double> val = phi.get_value(q_index);
 *         // do something with val
 *       }
 *   }
 * @endcode
 *
 * Likewise, a gradient of the finite element solution represented by @p
 * vector can be interpolated to the quadrature points by @p
 * phi.get_gradient(q). The combination of read_dof_values(), evaluate() and
 * get_value() is similar to what FEValues::get_function_values or
 * FEValues::get_function_gradients does, but it is in general much faster
 * because it makes use of the tensor product, see the description of the
 * evaluation routines below, and can do this operation for several cells at
 * once through vectorization.
 *
 * The second class of tasks done by FEEvaluation are integration tasks for
 * right hand sides. In finite element computations, these typically consist
 * of multiplying a quantity on quadrature points (a function value, or a
 * field interpolated by the finite element space itself) by a set of test
 * functions and integrating over the cell through summation of the values in
 * each quadrature point, multiplied by the quadrature weight and the Jacobian
 * determinant of the transformation. If a generic Function object is given
 * and we want to compute $v_i = \int_\Omega \varphi_i f dx$, this is done by
 * the following cell-wise integration:
 *
 * @code
 * FEEvaluation<dim,fe_degree> phi(matrix_free);
 * Function<dim> &function = ...;
 * for (unsigned int cell_index = cell_range.first;
 *      cell_index < cell_range.second; ++cell_index)
 *   {
 *     phi.reinit(cell_index);
 *     for (unsigned int q_index=0; q_index<phi.n_q_points; ++q_index)
 *       {
 *         Point<dim,VectorizedArray<double> > p_vect =
 *           phi.quadrature_point(q_index);
 *         // Need to evaluate function for each component in VectorizedArray
 *         VectorizedArray<double> f_value;
 *         for (unsigned int v=0; v<VectorizedArray<double>::n_array_elements; ++v)
 *           {
 *             Point<dim> p;
 *             for (unsigned int d=0; d<dim; ++d)
 *               p[d] = p_vect[d][v];
 *             f_value[v] = function.value(p);
 *           }
 *         phi.submit_value(f_value, q);
 *       }
 *     phi.integrate(true, false);
 *     phi.distribute_local_to_global(dst);
 *   }
 * @endcode
 *
 * In this code, the call to @p phi.submit_value() prepares for the
 * multiplication by the test function prior to the actual integration (inside
 * the submit call, the value to be tested is also multiplied by the
 * determinant of the Jacobian and the quadrature weight). In the
 * @p integrate() call, an integral contribution tested by each basis function
 * underlying the FEEvaluation object (e.g. the four linear shape functions of
 * FE_Q@<2@>(1) in 2D) is computed, which gives the vector entries to be
 * summed into the @p dst vector. Note that the above code needs to explicitly
 * loop over the components in the vectorized array for evaluating the
 * function, which is necessary for interfacing with a generic Function object
 * with double arguments. Simple functions can also be implemented in
 * VectorizedArray form directly as VectorizedArray provides the basic math
 * operations.
 *
 * For evaluating a bilinear form, the evaluation on a source vector is
 * combined with the integration involving test functions that get written
 * into a result vector. This setting is the context of matrix-free operator
 * evaluation and explained in the step-37 and step-48 tutorial programs.
 *
 * Note that the two vector accesses through FEEvaluation::read_dof_values and
 * FEEvaluation::distribute_local_to_global resolve constraints on the fly,
 * based on the ConstraintMatrix object specified at the MatrixFree::reinit()
 * call. In case the values in the degrees of freedom are of interest (usually
 * only the values in quadrature points are necessary), these can be accessed
 * through FEEvaluation::get_dof_value(i), where i is the index of the basis
 * function. Note that the numbering of the degrees of freedom for continuous
 * elements in FEEvaluation is different from the ordering in FE_Q (or
 * FEValues) because FEEvaluation needs to access them in lexicographic order,
 * which is the ordering used in FE_DGQ, for instance. Re-indexing would be
 * too expensive because the access inside evaluate() and integrate() is on
 * the critical path in the tensorial evaluation parts. An alternative to
 * filling the DoF values by read_dof_values() before an evaluate() call is to
 * manually assign a value by a set_dof_value() call. Likewise, if the local
 * result of integration should be further processed rather than scattered
 * into a vector by distribute_local_to_global(), one can access it by
 * get_dof_value() after an integrate() call. An example for using the values
 * of an integral in a different context is fast assembly of matrices as shown
 * in the next subsection.
 *
 * For most operator evaluation tasks that repeatedly go through the mesh, the
 * realization by MatrixFree that combines pre-computed data for the mapping
 * (Jacobian transformations for the geometry description) with on-the-fly
 * evaluation of basis functions is the most efficient way of doing things. In
 * other words, the framework selects a trade-off between memory usage and
 * initialization of objects that is suitable for replacement of matrix-vector
 * products or explicit time integration in a matrix-free way.
 *
 * <h4>Usage without pre-initialized MatrixFree object</h4>
 *
 * The second form of usage is to initialize FEEvaluation from geometry
 * information generated by FEValues. This allows to apply the integration
 * loops on the fly without prior initialization of MatrixFree objects. This
 * can be useful when the memory and initialization cost of MatrixFree is not
 * acceptable, e.g. when a different number of quadrature points should be
 * used for one single evaluation in error computation. Also, when using the
 * routines of this class to assemble matrices the trade-off implied by the
 * MatrixFree class may not be desired. In such a case, the cost to initialize
 * the necessary geometry data on the fly is comparably low and thus avoiding
 * a global object MatrixFree can be useful. When used in this way, reinit
 * methods reminiscent from FEValues with a cell iterator are used. However,
 * note that this model results in working on a single cell at a time, with
 * geometry data duplicated in all components of the vectorized array. Thus,
 * vectorization is only useful when it can apply the same operation on
 * different data, e.g. when performing matrix assembly.
 *
 * As an example, consider the following code to assemble the contributions to
 * the Laplace matrix:
 *
 * @code
 * FEEvaluation<dim,fe_degree> fe_eval (mapping, finite_element,
 *                                      QGauss<1>(fe_degree+1), flags);
 * for (cell = dof_handler.begin_active();
 *      cell != dof_handler.end();
 *      ++cell)
 *   {
 *     fe_eval.reinit(cell);
 *     for (unsigned int i=0; i<dofs_per_cell;
 *          i += VectorizedArray<double>::n_array_elements)
 *       {
 *         const unsigned int n_items =
 *           i+VectorizedArray<double>::n_array_elements > dofs_per_cell ?
 *           (dofs_per_cell - i) : VectorizedArray<double>::n_array_elements;
 *
 *         // Set n_items unit vectors
 *         for (unsigned int j=0; j<dofs_per_cell; ++j)
 *           fe_eval.set_dof_value(VectorizedArray<double>(), j);
 *         for (unsigned int v=0; v<n_items; ++v)
 *           {
 *             VectorizedArray<double> one_value = VectorizedArray<double>();
 *             one_value[v] = 1.;
 *             fe_eval.set_dof_value(one_value, i+v);
 *           }
 *
 *         // Apply operator on unit vector to generate the next few matrix
 *         // columns
 *         fe_eval.evaluate(true, true);
 *         for (unsigned int q=0; q<n_q_points; ++q)
 *           {
 *             fe_eval.submit_value(10.*fe_eval.get_value(q), q);
 *             fe_eval.submit_gradient(fe_eval.get_gradient(q), q);
 *           }
 *         fe_eval.integrate(true, true);
 *
 *         // Insert computed entries in matrix
 *         for (unsigned int v=0; v<n_items; ++v)
 *           for (unsigned int j=0; j<dofs_per_cell; ++j)
 *             cell_matrix(fe_eval.get_internal_dof_numbering()[j],
 *                         fe_eval.get_internal_dof_numbering()[i+v])
 *               = fe_eval.get_dof_value(j)[v];
 *       }
 *     ...
 *   }
 * @endcode
 *
 * This code generates the columns of the cell matrix with the loop over @p i
 * above. The way this is done is the following: FEEvaluation's routines focus
 * on the evaluation of finite element operators, so for computing a cell
 * matrix out of an operator evaluation it is applied to all the unit vectors
 * on the cell. Applying the operator on a unit vector might seem inefficient
 * but the evaluation routines used here are so quick that they still work
 * much faster than what is possible with FEValues. In particular, the
 * complexity is <code>(fe_degree+1)<sup>2*dim+1</sup> </code> rather than
 * <code>(fe_degree+1)<sup>3*dim</sup> </code>.
 *
 * Due to vectorization, we can generate matrix columns for several unit
 * vectors at a time (e.g. 4). The variable @p n_items make sure that we do
 * the last iteration where the number of cell dofs is not divisible by the
 * vectorization length correctly. Also note that we need to get the internal
 * dof numbering applied by fe_eval because FEEvaluation internally uses a
 * lexicographic numbering of degrees of freedom as explained above.
 *
 * <h4>Internal data organization</h4>
 *
 * The temporary data for holding the solution values on the local degrees of
 * freedom as well as the interpolated values, gradients, and Hessians on
 * quadrature points is a scratch array provided by
 * MatrixFree::acquire_scratch_data() that is re-used between different calls
 * to FEEvaluation. Therefore, constructing an FEEvaluation object is
 * typically cheap and does not involve any expensive operation. Only a few
 * dozen pointers to the actual data fields are set during
 * construction. Therefore, no negative performance impact arises when
 * creating an FEEvaluation several times per loop, such as at the top of a @p
 * local_cell_operation operation that is split in small chunks for a parallel
 * for loop, obviating a separate scratch data field for parallel loops as
 * necessary in the loop of @p WorkStream.
 *
 * When using the FEEvaluation class in multithreaded mode, the thread local
 * storage of the scratch data in MatrixFree automatically makes sure that
 * each thread gets it private data array. Note, however, that deal.II must be
 * compiled with thread support also when all the thread parallelization is
 * provided externally and not done via deal.II's routines, such as
 * OpenMP. This is because deal.II needs to know the notation of thread local
 * storage. The FEEvaluation kernels have been verified to work within OpenMP
 * loops.
 *
 * <h4>Vectorization scheme through VectorizedArray</h4>
 *
 * This class is designed to perform all arithmetics on single-instruction
 * multiple-data (SIMD) instructions present on modern CPUs by explicit
 * vectorization, which are made available in deal.II through the class
 * VectorizedArray, using the widest vector width available at
 * configure/compile time. In order to keep programs flexible, FEEvaluation
 * always applies vectorization over several elements. This is often the best
 * compromise because computations on different elements are usually
 * independent in the finite element method (except of course the process of
 * adding an integral contribution to a global residual vector), also in more
 * complicated scenarios: Stabilization parameter can e.g. be defined as the
 * maximum of some quantities on all quadrature points of a cell divided by
 * the cell's volume, but without locally mixing the results with
 * neighbors. Using the terminology from computer architecture, the design of
 * FEEvaluation relies on not doing any cross-lane data exchange when
 * operating on the cell in typical integration scenarios.
 *
 * When the number of cells in the problem is not a multiple of the number of
 * array elements in the SIMD vector, the implementation of FEEvaluation fills
 * in some dummy entries in the unused SIMD lanes and carries them around
 * nonetheless, a choice made necessary since the length of VectorizedArray is
 * fixed at compile time. Yet, this approach most often results in superior
 * code as compared to an auto-vectorization setup where an alternative
 * unvectorized code path would be necessary next to the vectorized version to
 * be used on fully populated lanes, together with a dispatch mechanism. In
 * @p read_dof_values, the empty lanes resulting from a reinit() call to an
 * incomplete batch of cells are set to zero, whereas
 * @p distribute_local_to_global or @p set_dof_values simply ignores the
 * content in the empty lanes. The number of actually filled SIMD lanes can by
 * queried by MatrixFree::n_components_filled().
 *
 * Obviously, the computations performed on the artificial lanes (without real
 * data) should never be mixed with valid results. The contract in using this
 * class is that the user makes sure that lanes are not crossed in user code,
 * in particular since it is not clear a priori which cells are going to be
 * put together in vectorization. For example, results on an element should
 * not be added to results on other elements except through the global vector
 * access methods or by access that is masked by
 * MatrixFree::n_components_filled(). No guarantee can be made that results on
 * artificial lanes will always be zero that can safely be added to other
 * results: The data on JxW or Jacobians is copied from the last valid lane in
 * order to avoid division by zero that could trigger floating point
 * exceptions or trouble in other situations.
 *
 * <h3>Description of evaluation routines</h3>
 *
 * This class contains specialized evaluation routines for elements based on
 * tensor-product quadrature formulas and tensor-product-like shape functions,
 * including standard FE_Q or FE_DGQ elements and quadrature points symmetric
 * around 0.5 (like Gauss quadrature), FE_DGP elements based on truncated
 * tensor products as well as the faster case of Gauss-Lobatto elements with
 * Gauss-Lobatto quadrature which give diagonal mass matrices and quicker
 * evaluation internally. The main benefit of this class is the evaluation of
 * all shape functions in all quadrature or integration over all shape
 * functions in <code>dim (fe_degree+1)<sup>dim+1</sup> </code> operations
 * instead of the slower <code> (fe_degree+1)<sup>2*dim</sup></code>
 * complexity in the evaluation routines of FEValues. This is done by an
 * algorithm called sum factorization which factors out constant factors
 * during the evaluation along a coordinate direction. This algorithm is the
 * basis of many spectral element algorithms.
 *
 * Note that many of the operations available through this class are inherited
 * from the base class FEEvaluationBase, in particular reading from and
 * writing to vectors. Also, the class inherits from FEEvaluationAccess that
 * implements access to values, gradients and Hessians of the finite element
 * function on quadrature points.
 *
 * This class assumes that the shape functions of the FiniteElement under
 * consideration do <em>not</em> depend on the geometry of the cells in real
 * space. Currently, other finite elements cannot be treated with the
 * matrix-free concept.
 *
 * <h4>Degree of finite element as a compile-time parameter</h4>
 *
 * The default usage model of FEEvaluation expects the polynomial degree to be
 * given as a template parameter. This requirement guarantees maximum
 * efficiency: The sum factorization evaluation performs a number of nested
 * short 1D loops of length equal to the polynomial degree plus one. If the
 * loop bounds are known at compile time, the compiler can unroll loops as
 * deemed most efficient by its heuristics. At least the innermost loop is
 * almost always completely unrolled, avoiding the loop overhead.
 *
 * However, carrying the polynomial degree (and the number of quadrature
 * points) as a template parameter makes things more complicated in codes
 * where different polynomial degrees should be considered, e.g. in
 * application codes where the polynomial degree is given through an input
 * file. The recommended approach for good performance is to create different
 * cell functions, possibly through different operator classes that are
 * derived from a common base class with virtual functions to access the
 * operator evaluation. The program then keeps only a pointer to the common
 * base class after initializing templated derived classes with fixed
 * polynomial degree that are selected from the detected polynomial
 * degree. This approach requires a-priori knowledge of the range of valid
 * degrees and can lead to rather long compile times in programs with many
 * apply functions.
 *
 * A flexible choice of the polynomial degree based on the information in the
 * element passed to the initialization is also supported by FEEvaluation,
 * even though it runs two to three times more slowly. For this, set the
 * template parameter for the polynomial degree to -1 (and choose an arbitrary
 * number for the number of quadrature points), which switches to the other
 * code path. Internally, an evaluator object with template specialization for
 * -1 is invoked that runs according to run-time bounds, whereas the default
 * case with compile-time bounds given through fe_degree selects another
 * template class without compromising efficiency.
 *
 * An overview of the performance of FEEvaluation is given in the following
 * figure. It considers the time spent per degree of freedom for evaluating
 * the Laplacian with continuous finite elements using a code similar to the
 * step-37 tutorial program for single-precision arithmetics. The time is
 * based on an experiment on a single core of an Intel Xeon E5-2687W v4,
 * running at 3.4 GHz and measured at problem sizes around 10 million. The
 * plot lists the computational time (around 0.1 seconds) divided by the
 * number of degrees freedom.
 *
 * @image html fe_evaluation_laplacian_time_per_dof.png
 *
 * The figure shows that the templated version is between 2.5 and 3 times
 * faster. The fastest turnaround on this setup is for polynomial degree 5 at
 * 7.4e-9 seconds per degree of freedom or 134 million degrees of freedom per
 * second - on a single core. The non-templated version is also fastest at
 * polynomial degree 5 with 2.1e-9 seconds per degree of freedom or 48 million
 * degrees of freedom per second.
 *
 * <h3>Handling multi-component systems</h3>
 *
 * FEEvaluation also allows for treating vector-valued problems through a
 * template parameter on the number of components:
 *
 * @code
 * FEEvaluation<dim,fe_degree,n_q_points_1d,n_components> phi(matrix_free);
 * @endcode
 *
 * If used this way, the components can be gathered from several components of
 * an @p std::vector<VectorType> through the call
 *
 * @code
 * phi.read_dof_values(src, 0);
 * @endcode
 *
 * where the 0 means that the vectors starting from the zeroth vector in the
 * @p std::vector should be used, <code>src[0], src[1], ...,
 * src[n_components-1]</code>.
 *
 * An alternative way for reading multi-component systems is possible if the
 * DoFHandler underlying the MatrixFree data is based on an FESystem of @p
 * n_components entries. In that case, a single vector is provided for the
 * read_dof_values() and distribute_local_to_global() calls.
 *
 * An important property of FEEvaluation in multi-component systems is the
 * layout of multiple components in the get_value(), get_gradient(), or
 * get_dof_value() calls. In this case, instead of a scalar return field
 * VectorizedArray@<double@> a tensor is returned,
 *
 * @code
 * get_value -> Tensor<1,n_components,VectorizedArray<double> >
 * get_gradient -> Tensor<1,n_components,Tensor<1,dim,VectorizedArray<double> > >
 * @endcode
 *
 * In a similar vein, the submit_value() and submit_gradient() calls take
 * tensors of values. Note that there exist specializations for @p
 * n_components=1 and @p n_components=dim, which are provided through the base
 * class FEEvaluationAccess. In the scalar case, these provide the scalar
 * return types described above. In the vector-valued case, the gradient is
 * converted from <code>Tensor@<1,dim,Tensor@<1,dim,VectorizedArray@<double@>
 * @> @></code> to <code>Tensor@<2,dim,VectorizedArray@<double@>
 * @></code>. Furthermore, additional operations such as the diveregence or
 * curl are available.
 *
 * In case different shape functions are combined, for example mixed finite
 * element formulations in Stokes flow, two FEEvaluation objects are created,
 * one for the velocity and one for the pressure. Those are then combined on
 * quadrature points:
 *
 * @code
 * FEEvaluation<dim,degree_p+1,degree_p+2,dim> velocity (data, 0);
 * FEEvaluation<dim,degree_p,  degree_p+2,1, > pressure (data, 1);
 *
 * for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
 *   {
 *     velocity.reinit (cell);
 *     velocity.read_dof_values (src.block(0));
 *     velocity.evaluate (false,true);
 *     pressure.reinit (cell);
 *     pressure.read_dof_values (src.block(1));
 *     pressure.evaluate (true,false);
 *
 *     for (unsigned int q=0; q<velocity.n_q_points; ++q)
 *       {
 *         SymmetricTensor<2,dim,VectorizedArray<double> > sym_grad_u =
 *           velocity.get_symmetric_gradient (q);
 *         VectorizedArray<double> pres = pressure.get_value(q);
 *         VectorizedArray<double> div = -trace(sym_grad_u);
 *         pressure.submit_value (div, q);
 *
 *         // subtract p * I
 *         for (unsigned int d=0; d<dim; ++d)
 *           sym_grad_u[d][d] -= pres;
 *
 *         velocity.submit_symmetric_gradient(sym_grad_u, q);
 *      }
 *
 *     velocity.integrate (false,true);
 *     velocity.distribute_local_to_global (dst.block(0));
 *     pressure.integrate (true,false);
 *     pressure.distribute_local_to_global (dst.block(1));
 *   }
 * @endcode
 *
 * This code assumes that a BlockVector of two components describes the
 * velocity and pressure components, respectively. For identifying the
 * different DoFHandler objects for velocity and pressure, the second argument
 * to the FEEvaluation objects specify the respective component 0 for velocity
 * and 1 for pressure. For further examples of vector-valued problems, the
 * deal.II test suite includes a few additional examples as well, e.g. the
 * Stokes operator described above is found at
 * https://github.com/dealii/dealii/blob/master/tests/matrix_free/matrix_vector_stokes_noflux.cc
 *
 * <h3>Handling several integration tasks and data storage in quadrature points</h3>
 *
 * The design of FEEvaluation and MatrixFree separates the geometry from the
 * basis functions. Therefore, several DoFHandler objects (or the same
 * DoFHandler equipped with different constraint objects) can share the same
 * geometry information like in the Stokes example above. All geometry is
 * cached once in MatrixFree, so FEEvaluation does not need to do expensive
 * initialization calls and rather sets a few pointers. This realization is
 * based on the idea that the geometry information is needed only once also
 * when several fields are evaluated, in a departure from FEValues which sets
 * up the internal mapping data for each field. If for example a
 * multi-component PDE involves the shape values on one component and the
 * shape gradient on the other, no efficiency is lost if both are based on the
 * same MatrixFree object where the update flags specify that both @p
 * update_values , @p update_gradients , and @p update_jxw_values are
 * given. The selection of desired quantities of shape values is through the
 * flags in the evaluate() or integrate calls and the access at quadrature
 * points:
 *
 * @code
 * phi1.evaluate(true, false);
 * phi2.evaluate(false, true);
 * for (unsigned int q_index=0; q_index<phi1.n_q_points; ++q_index)
 *   {
 *     VectorizedArray<double> val1 = phi1.get_value(q);
 *     Tensor<1,dim,VectorizedArray<double> > grad2 = phi2.get_gradient(q);
 *     Point<dim,VectorizedArray<double> > point = phi1.get_quadrature_point(q);
 *     // ... some complicated formula combining those three...
 *   }
 * @endcode
 *
 * In the loop over quadrature points, one can ask any of the two FEEvaluation
 * objects &mdash; it does not really matter which one because they only keep
 * pointers to the quadrature point data &mdash; to provide the quadrature
 * point location.
 *
 * This observation also translates to the case when different differential
 * operators are implemented in a program, for example the action of a mass
 * matrix for one phase of the algorithm and the action of a stiffness matrix
 * in another one. Only a single MatrixFree object is necessary, maintaining
 * full efficiency by using different local functions with the respective
 * implementation in separate FEEvaluation objects. In other words, a user
 * does not need to bother about being conservative when providing
 * update_flags to the initialization of MatrixFree for efficiency reasons -
 * no overhead incurs inside FEEvaluation, except for at most one or two more
 * @p if statements inside the FEEvaluation::reinit() call. Rather, the
 * largest set of flags necessary among all calls is perfectly fine from an
 * efficiency point of view.
 *
 * For the combination of different fields, including different solution
 * vectors that come from different time steps, it is mandatory that all
 * FEEvaluation objects share the same MatrixFree object. This is because the
 * way cells are looped by MatrixFree::cell_loop() can be different for
 * different DoFHandler or ConstraintMatrix arguments. More precisely, even
 * though the layout is going to be the same in serial, there is no guarantee
 * about the ordering for different DoFHandler/Constraints in the MPI
 * case. The reason is that the algorithm detects cells that need data
 * exchange with MPI and those can change for different elements &mdash; FE_Q
 * with hanging node constraints connects to more neighbors than a FE_DGQ
 * element, for instance, and cells which need data exchange are put in
 * different positions inside the cell loop. Of course, if the exact same
 * DoFHandler, ConstraintMatrix, and options (such as the setting for thread
 * parallelism) are set, then the order is going to be the same because the
 * algorithm is deterministic.
 *
 * @tparam dim Dimension in which this class is to be used
 *
 * @tparam fe_degree Degree of the tensor product finite element with
 * fe_degree+1 degrees of freedom per coordinate direction. Can be set to -1
 * if the degree is not known at compile time, but performance will usually be
 * worse by a factor of 2-3.
 *
 * @tparam n_q_points_1d Number of points in the quadrature formula in 1D,
 * defaults to fe_degree+1
 *
 * @tparam n_components Number of vector components when solving a system of
 * PDEs. If the same operation is applied to several components of a PDE (e.g.
 * a vector Laplace equation), they can be applied simultaneously with one
 * call (and often more efficiently). Defaults to 1.
 *
 * @tparam Number Number format, usually @p double or @p float. Defaults to @p
 * double
 *
 * @author Katharina Kormann and Martin Kronbichler, 2010, 2011
 */
template <int dim, int fe_degree, int n_q_points_1d, int n_components_,
          typename Number >
class FEEvaluation : public FEEvaluationAccess<dim,n_components_,Number>
{
public:
  typedef FEEvaluationAccess<dim,n_components_,Number> BaseClass;
  typedef Number                            number_type;
  typedef typename BaseClass::value_type    value_type;
  typedef typename BaseClass::gradient_type gradient_type;
  static constexpr unsigned int dimension     = dim;
  static constexpr unsigned int n_components  = n_components_;
  static constexpr unsigned int static_n_q_points    = Utilities::fixed_int_power<n_q_points_1d,dim>::value;
  static constexpr unsigned int static_dofs_per_component = Utilities::fixed_int_power<fe_degree+1,dim>::value;
  static constexpr unsigned int tensor_dofs_per_cell = static_dofs_per_component *n_components;
  static constexpr unsigned int static_dofs_per_cell = static_dofs_per_component *n_components;

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
   * Constructor that comes with reduced functionality and works similar as
   * FEValues. The arguments are similar to the ones passed to the constructor
   * of FEValues, with the notable difference that FEEvaluation expects a one-
   * dimensional quadrature formula, Quadrature<1>, instead of a @p dim
   * dimensional one. The finite element can be both scalar or vector valued,
   * but this method always only selects a scalar base element at a time (with
   * @p n_components copies as specified by the class template). For vector-
   * valued elements, the optional argument @p first_selected_component allows
   * to specify the index of the base element to be used for evaluation. Note
   * that the internal data structures always assume that the base element is
   * primitive, non-primitive are not supported currently.
   *
   * As known from FEValues, a call to the reinit method with a
   * Triangulation<dim>::cell_iterator is necessary to make the geometry and
   * degrees of freedom of the current class known. If the iterator includes
   * DoFHandler information (i.e., it is a DoFHandler<dim>::cell_iterator or
   * similar), the initialization allows to also read from or write to vectors
   * in the standard way for DoFHandler<dim>::active_cell_iterator types for
   * one cell at a time. However, this approach is much slower than the path
   * with MatrixFree with MPI since index translation has to be done. As only
   * one cell at a time is used, this method does not vectorize over several
   * elements (which is most efficient for vector operations), but only
   * possibly within the element if the evaluate/integrate routines are
   * combined inside user code (e.g. for computing cell matrices).
   */
  FEEvaluation (const Mapping<dim>       &mapping,
                const FiniteElement<dim> &fe,
                const Quadrature<1>      &quadrature,
                const UpdateFlags         update_flags,
                const unsigned int        first_selected_component = 0);

  /**
   * Constructor for the reduced functionality. This constructor is equivalent
   * to the other one except that it makes the object use a $Q_1$ mapping
   * (i.e., an object of type MappingQGeneric(1)) implicitly.
   */
  FEEvaluation (const FiniteElement<dim> &fe,
                const Quadrature<1>      &quadrature,
                const UpdateFlags         update_flags,
                const unsigned int        first_selected_component = 0);

  /**
   * Constructor for the reduced functionality. Similar to the other
   * constructor with FiniteElement argument but using another
   * FEEvaluationBase object to provide info about the geometry. This allows
   * several FEEvaluation objects to share the geometry evaluation, i.e., the
   * underlying mapping and quadrature points do only need to be evaluated
   * once. Make sure to not pass an optional object around when you intend to
   * use the FEEvaluation object in %parallel to the given one because
   * otherwise the intended sharing may create race conditions.
   */
  template <int n_components_other>
  FEEvaluation (const FiniteElement<dim> &fe,
                const FEEvaluationBase<dim,n_components_other,Number> &other,
                const unsigned int        first_selected_component = 0);

  /**
   * Copy constructor. If FEEvaluationBase was constructed from a mapping, fe,
   * quadrature, and update flags, the underlying geometry evaluation based on
   * FEValues will be deep-copied in order to allow for using in parallel with
   * threads.
   */
  FEEvaluation (const FEEvaluation &other);

  /**
   * Copy assignment operator. If FEEvaluationBase was constructed from a
   * mapping, fe, quadrature, and update flags, the underlying geometry
   * evaluation based on FEValues will be deep-copied in order to allow for
   * using in parallel with threads.
   */
  FEEvaluation &operator= (const FEEvaluation &other);

  /**
   * Evaluates the function values, the gradients, and the Hessians of the
   * FE function given at the DoF values in the input vector at the quadrature
   * points on the unit cell.  The function arguments specify which parts
   * shall actually be computed. Needs to be called before the functions @p
   * get_value(), @p get_gradient() or @p get_laplacian give useful
   * information (unless these values have been set manually).
   */
  void evaluate (const bool evaluate_values,
                 const bool evaluate_gradients,
                 const bool evaluate_hessians = false);

  /**
   * This function takes the values and/or gradients that are stored on
   * quadrature points, tests them by all the basis functions/gradients on the
   * cell and performs the cell integration. The two function arguments
   * @p integrate_values and @p integrate_gradients define which of the values
   * or gradients (or both) are summed together.
   */
  void integrate (const bool integrate_values,
                  const bool integrate_gradients);

  /**
   * Return the q-th quadrature point stored in MappingInfo.
   */
  Point<dim,VectorizedArray<Number> >
  quadrature_point (const unsigned int q_point) const;

  /**
   * The number of degrees of freedom of a single component on the cell for
   * the underlying evaluation object. Usually close to
   * static_dofs_per_component, but the number depends on the actual element
   * selected and is thus not static.
   */
  const unsigned int dofs_per_component;

  /**
   * The number of degrees of freedom on the cell accumulated over all
   * components in the current evaluation object. Usually close to
   * static_dofs_per_cell = static_dofs_per_component*n_components, but the
   * number depends on the actual element selected and is thus not static.
   */
  const unsigned int dofs_per_cell;

  /**
   * The number of quadrature points in use for the current evaluation
   * object. It equals static_n_q_points in case the element degree is not set
   * to -1, i.e., the worker kernels can be set according to the templates @p
   * fe_degree and @p n_q_points_1d rather than at run time.
   */
  const unsigned int n_q_points;

private:
  /**
   * Checks if the template arguments regarding degree of the element
   * corresponds to the actual element used at initialization.
   */
  void check_template_arguments(const unsigned int fe_no,
                                const unsigned int first_selected_component);
};



namespace internal
{
  namespace MatrixFreeFunctions
  {
    // a helper function to compute the number of DoFs of a DGP element at compile
    // time, depending on the degree
    template <int dim, int degree>
    struct DGP_dofs_per_component
    {
      // this division is always without remainder
      static constexpr unsigned int value =
        (DGP_dofs_per_component<dim-1,degree>::value * (degree+dim)) / dim;
    };

    // base specialization: 1d elements have 'degree+1' degrees of freedom
    template <int degree>
    struct DGP_dofs_per_component<1,degree>
    {
      static constexpr unsigned int value = degree+1;
    };
  }
}


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
  scratch_data_array (data_in.acquire_scratch_data()),
  quad_no            (quad_no_in),
  n_fe_components    (data_in.get_dof_info(fe_no_in).n_components),
  active_fe_index    (fe_degree != numbers::invalid_unsigned_int ?
                      data_in.get_dof_info(fe_no_in).fe_index_from_degree(fe_degree)
                      :
                      0),
  active_quad_index  (fe_degree != numbers::invalid_unsigned_int ?
                      data_in.get_mapping_info().
                      mapping_data_gen[quad_no_in].
                      quad_index_from_n_q_points(n_q_points)
                      :
                      0),
  matrix_info        (&data_in),
  dof_info           (&data_in.get_dof_info(fe_no_in)),
  mapping_info       (&data_in.get_mapping_info()),
  data               (&data_in.get_shape_info
                      (fe_no_in, quad_no_in, active_fe_index,
                       active_quad_index)),
  cartesian_data     (nullptr),
  jacobian           (nullptr),
  J_value            (nullptr),
  quadrature_weights (mapping_info->mapping_data_gen[quad_no].
                      quadrature_weights[active_quad_index].begin()),
  quadrature_points  (nullptr),
  jacobian_grad      (nullptr),
  jacobian_grad_upper(nullptr),
  cell               (numbers::invalid_unsigned_int),
  cell_type          (internal::MatrixFreeFunctions::undefined),
  cell_data_number   (numbers::invalid_unsigned_int),
  dof_values_initialized    (false),
  values_quad_initialized   (false),
  gradients_quad_initialized(false),
  hessians_quad_initialized (false),
  values_quad_submitted     (false),
  gradients_quad_submitted  (false),
  first_selected_component  (0)
{
  set_data_pointers();
  Assert (matrix_info->mapping_initialized() == true,
          ExcNotInitialized());
  AssertDimension (matrix_info->get_size_info().vectorization_length,
                   VectorizedArray<Number>::n_array_elements);
  AssertDimension (data->dofs_per_component_on_cell*n_fe_components,
                   dof_info->dofs_per_cell[active_fe_index]);
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
template <int n_components_other>
inline
FEEvaluationBase<dim,n_components_,Number>
::FEEvaluationBase (const Mapping<dim>       &mapping,
                    const FiniteElement<dim> &fe,
                    const Quadrature<1>      &quadrature,
                    const UpdateFlags         update_flags,
                    const unsigned int        first_selected_component,
                    const FEEvaluationBase<dim,n_components_other,Number> *other)
  :
  scratch_data_array (new AlignedVector<VectorizedArray<Number> >()),
  quad_no            (numbers::invalid_unsigned_int),
  n_fe_components    (n_components_),
  active_fe_index    (numbers::invalid_unsigned_int),
  active_quad_index  (numbers::invalid_unsigned_int),
  matrix_info        (nullptr),
  dof_info           (nullptr),
  mapping_info       (nullptr),
  // select the correct base element from the given FE component
  data               (new internal::MatrixFreeFunctions::ShapeInfo<VectorizedArray<Number>>(quadrature, fe, fe.component_to_base_index(first_selected_component).first)),
  cartesian_data     (nullptr),
  jacobian           (nullptr),
  J_value            (nullptr),
  quadrature_weights (nullptr),
  quadrature_points  (nullptr),
  jacobian_grad      (nullptr),
  jacobian_grad_upper(nullptr),
  cell               (0),
  cell_type          (internal::MatrixFreeFunctions::general),
  cell_data_number   (numbers::invalid_unsigned_int),
  dof_values_initialized    (false),
  values_quad_initialized   (false),
  gradients_quad_initialized(false),
  hessians_quad_initialized (false),
  values_quad_submitted     (false),
  gradients_quad_submitted  (false),
  // keep the number of the selected component within the current base element
  // for reading dof values
  first_selected_component  (fe.component_to_base_index(first_selected_component).second)
{
  const unsigned int base_element_number =
    fe.component_to_base_index(first_selected_component).first;
  set_data_pointers();

  Assert(other == nullptr || other->mapped_geometry.get() != nullptr, ExcInternalError());
  if (other != nullptr &&
      other->mapped_geometry->get_quadrature() == quadrature)
    mapped_geometry = other->mapped_geometry;
  else
    mapped_geometry
      = std::make_shared<internal::MatrixFreeFunctions::MappingDataOnTheFly<dim,Number> >
        (mapping, quadrature, update_flags);
  jacobian = mapped_geometry->get_inverse_jacobians().begin();
  J_value = mapped_geometry->get_JxW_values().begin();
  quadrature_points = mapped_geometry->get_quadrature_points().begin();

  Assert(fe.element_multiplicity(base_element_number) == 1 ||
         fe.element_multiplicity(base_element_number)-first_selected_component >= n_components_,
         ExcMessage("The underlying element must at least contain as many "
                    "components as requested by this class"));
  (void) base_element_number;
}



template <int dim, int n_components_, typename Number>
inline
FEEvaluationBase<dim,n_components_,Number>
::FEEvaluationBase (const FEEvaluationBase<dim,n_components_,Number> &other)
  :
  scratch_data_array (other.matrix_info == nullptr ?
                      new AlignedVector<VectorizedArray<Number> >() :
                      other.matrix_info->acquire_scratch_data()),
  quad_no            (other.quad_no),
  n_fe_components    (other.n_fe_components),
  active_fe_index    (other.active_fe_index),
  active_quad_index  (other.active_quad_index),
  matrix_info        (other.matrix_info),
  dof_info           (other.dof_info),
  mapping_info       (other.mapping_info),
  data               (other.matrix_info == nullptr ?
                      new internal::MatrixFreeFunctions::ShapeInfo<VectorizedArray<Number>>(*other.data) :
                      other.data),
  cartesian_data     (nullptr),
  jacobian           (nullptr),
  J_value            (nullptr),
  quadrature_weights (mapping_info != nullptr ?
                      mapping_info->mapping_data_gen[quad_no].
                      quadrature_weights[active_quad_index].begin()
                      :
                      nullptr),
  quadrature_points  (nullptr),
  jacobian_grad      (nullptr),
  jacobian_grad_upper(nullptr),
  cell               (numbers::invalid_unsigned_int),
  cell_type          (internal::MatrixFreeFunctions::general),
  cell_data_number   (numbers::invalid_unsigned_int),
  dof_values_initialized    (false),
  values_quad_initialized   (false),
  gradients_quad_initialized(false),
  hessians_quad_initialized (false),
  values_quad_submitted     (false),
  gradients_quad_submitted  (false),
  first_selected_component  (other.first_selected_component)
{
  set_data_pointers();

  // Create deep copy of mapped geometry for use in parallel...
  if (other.mapped_geometry.get() != nullptr)
    {
      mapped_geometry.reset
      (new internal::MatrixFreeFunctions::
       MappingDataOnTheFly<dim,Number>(other.mapped_geometry->get_fe_values().get_mapping(),
                                       other.mapped_geometry->get_quadrature(),
                                       other.mapped_geometry->get_fe_values().get_update_flags()));
      jacobian = mapped_geometry->get_inverse_jacobians().begin();
      J_value = mapped_geometry->get_JxW_values().begin();
      quadrature_points = mapped_geometry->get_quadrature_points().begin();
      cell = 0;
    }
}



template <int dim, int n_components_, typename Number>
inline
FEEvaluationBase<dim,n_components_,Number> &
FEEvaluationBase<dim,n_components_,Number>
::operator= (const FEEvaluationBase<dim,n_components_,Number> &other)
{
  AssertDimension(quad_no, other.quad_no);
  AssertDimension(n_fe_components, other.n_fe_components);
  AssertDimension(active_fe_index, other.active_fe_index);
  AssertDimension(active_quad_index, other.active_quad_index);
  AssertDimension(first_selected_component, other.first_selected_component);

  // release old memory
  if (matrix_info == nullptr)
    {
      delete data;
      delete scratch_data_array;
    }
  else
    {
      matrix_info->release_scratch_data(scratch_data_array);
    }

  matrix_info = other.matrix_info;
  dof_info = other.dof_info;
  mapping_info = other.mapping_info;
  if (other.matrix_info == nullptr)
    {
      data = new internal::MatrixFreeFunctions::ShapeInfo<VectorizedArray<Number>>(*other.data);
      scratch_data_array = new AlignedVector<VectorizedArray<Number> >();
    }
  else
    {
      data = other.data;
      scratch_data_array = matrix_info->acquire_scratch_data();
    }
  set_data_pointers();

  cartesian_data = nullptr;
  jacobian = nullptr;
  J_value = nullptr;
  quadrature_weights = mapping_info != nullptr ?
                       mapping_info->mapping_data_gen[quad_no].
                       quadrature_weights[active_quad_index].begin()
                       :
                       nullptr;
  quadrature_points = nullptr;
  jacobian_grad = nullptr;
  jacobian_grad_upper = nullptr;
  cell = numbers::invalid_unsigned_int;
  cell_type = internal::MatrixFreeFunctions::general;
  cell_data_number = numbers::invalid_unsigned_int;

  // Create deep copy of mapped geometry for use in parallel...
  if (other.mapped_geometry.get() != nullptr)
    {
      mapped_geometry.reset
      (new internal::MatrixFreeFunctions::
       MappingDataOnTheFly<dim,Number>(other.mapped_geometry->get_fe_values().get_mapping(),
                                       other.mapped_geometry->get_quadrature(),
                                       other.mapped_geometry->get_fe_values().get_update_flags()));
      jacobian = mapped_geometry->get_inverse_jacobians().begin();
      J_value = mapped_geometry->get_JxW_values().begin();
      quadrature_points = mapped_geometry->get_quadrature_points().begin();
      cell = 0;
    }

  return *this;
}



template <int dim, int n_components_, typename Number>
inline
FEEvaluationBase<dim,n_components_,Number>::~FEEvaluationBase ()
{
  if (matrix_info != nullptr)
    {
      try
        {
          matrix_info->release_scratch_data(scratch_data_array);
        }
      catch (...)
        {}
    }
  else
    {
      delete scratch_data_array;
      delete data;
    }
}



template <int dim, int n_components_, typename Number>
inline
void
FEEvaluationBase<dim,n_components_,Number>
::set_data_pointers()
{
  Assert(scratch_data_array != nullptr, ExcInternalError());

  const unsigned int tensor_dofs_per_component =
    Utilities::fixed_power<dim>(this->data->fe_degree+1);
  const unsigned int dofs_per_component = this->data->dofs_per_component_on_cell;
  const unsigned int n_quadrature_points = this->data->n_q_points;

  const unsigned int shift = std::max(tensor_dofs_per_component+1, dofs_per_component)*
                             n_components_*3 + 2 * n_quadrature_points;
  const unsigned int allocated_size = shift + n_components_ * dofs_per_component
                                      + (n_components_*(dim*dim+2*dim+1)*n_quadrature_points);
  scratch_data_array->resize_fast(allocated_size);

  // set the pointers to the correct position in the data array
  for (unsigned int c=0; c<n_components_; ++c)
    {
      this->values_dofs[c] = scratch_data_array->begin() + c*dofs_per_component;
      this->values_quad[c] = scratch_data_array->begin() +
                             n_components*dofs_per_component+c*n_quadrature_points;
      for (unsigned int d=0; d<dim; ++d)
        this->gradients_quad[c][d] = scratch_data_array->begin() +
                                     n_components*(dofs_per_component+n_quadrature_points) +
                                     (c*dim+d)*n_quadrature_points;
      for (unsigned int d=0; d<(dim*dim+dim)/2; ++d)
        this->hessians_quad[c][d] = scratch_data_array->begin() +
                                    n_components*((dim+1)*n_quadrature_points + dofs_per_component) +
                                    (c*(dim*dim+dim)+d)*n_quadrature_points;
    }
  scratch_data = scratch_data_array->begin() + n_components_ * dofs_per_component
                 + (n_components_*(dim*dim+2*dim+1)*n_quadrature_points);
}



template <int dim, int n_components_, typename Number>
inline
void
FEEvaluationBase<dim,n_components_,Number>::reinit (const unsigned int cell_in)
{
  Assert (mapped_geometry == nullptr,
          ExcMessage("FEEvaluation was initialized without a matrix-free object."
                     " Integer indexing is not possible"));
  if (mapped_geometry != nullptr)
    return;
  Assert (dof_info != nullptr, ExcNotInitialized());
  Assert (mapping_info != nullptr, ExcNotInitialized());
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
template <typename DoFHandlerType, bool level_dof_access>
inline
void
FEEvaluationBase<dim,n_components_,Number>
::reinit (const TriaIterator<DoFCellAccessor<DoFHandlerType,level_dof_access> > &cell)
{
  Assert(matrix_info == nullptr,
         ExcMessage("Cannot use initialization from cell iterator if "
                    "initialized from MatrixFree object. Use variant for "
                    "on the fly computation with arguments as for FEValues "
                    "instead"));
  Assert(mapped_geometry.get() != nullptr, ExcNotInitialized());
  mapped_geometry->reinit(static_cast<typename Triangulation<dim>::cell_iterator>(cell));
  local_dof_indices.resize(cell->get_fe().dofs_per_cell);
  if (level_dof_access)
    cell->get_mg_dof_indices(local_dof_indices);
  else
    cell->get_dof_indices(local_dof_indices);
}



template <int dim, int n_components_, typename Number>
inline
void
FEEvaluationBase<dim,n_components_,Number>
::reinit (const typename Triangulation<dim>::cell_iterator &cell)
{
  Assert(matrix_info == 0,
         ExcMessage("Cannot use initialization from cell iterator if "
                    "initialized from MatrixFree object. Use variant for "
                    "on the fly computation with arguments as for FEValues "
                    "instead"));
  Assert(mapped_geometry.get() != 0, ExcNotInitialized());
  mapped_geometry->reinit(cell);
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
const internal::MatrixFreeFunctions::ShapeInfo<VectorizedArray<Number>> &
    FEEvaluationBase<dim,n_components_,Number>::get_shape_info() const
{
  Assert(data != nullptr, ExcInternalError());
  return *data;
}



template <int dim, int n_components_, typename Number>
inline
void
FEEvaluationBase<dim,n_components_,Number>
::fill_JxW_values(AlignedVector<VectorizedArray<Number> > &JxW_values) const
{
  AssertDimension(JxW_values.size(), data->n_q_points);
  Assert (this->J_value != nullptr, ExcNotInitialized());
  if (this->cell_type == internal::MatrixFreeFunctions::cartesian ||
      this->cell_type == internal::MatrixFreeFunctions::affine)
    {
      Assert (this->mapping_info != nullptr, ExcNotImplemented());
      VectorizedArray<Number> J = this->J_value[0];
      for (unsigned int q=0; q<this->data->n_q_points; ++q)
        JxW_values[q] = J * this->quadrature_weights[q];
    }
  else
    for (unsigned int q=0; q<data->n_q_points; ++q)
      JxW_values[q] = this->J_value[q];
}



template <int dim, int n_components_, typename Number>
inline
VectorizedArray<Number>
FEEvaluationBase<dim,n_components_,Number>::JxW(const unsigned int q_point) const
{
  Assert (this->J_value != nullptr, ExcNotInitialized());
  if (this->cell_type == internal::MatrixFreeFunctions::cartesian ||
      this->cell_type == internal::MatrixFreeFunctions::affine)
    {
      Assert (this->mapping_info != nullptr, ExcInternalError());
      return this->J_value[0] * this->quadrature_weights[q_point];
    }
  else
    return this->J_value[q_point];
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
  vector_access (LinearAlgebra::distributed::Vector<Number> &vec,
                 const unsigned int                          entry)
  {
    return vec.local_element(entry);
  }



  // read access to distributed MPI vectors that have a local_element(uint)
  // method to access data in local index space, which is what we use in
  // DoFInfo and hence in read_dof_values etc.
  template <typename Number>
  inline
  Number
  vector_access (const LinearAlgebra::distributed::Vector<Number> &vec,
                 const unsigned int                                entry)
  {
    return vec.local_element(entry);
  }



  // this is to make sure that the parallel partitioning in the
  // LinearAlgebra::distributed::Vector is really the same as stored in
  // MatrixFree
  template <typename VectorType>
  inline
  void check_vector_compatibility (const VectorType                             &vec,
                                   const internal::MatrixFreeFunctions::DoFInfo &dof_info)
  {
    (void) vec;
    (void) dof_info;

    AssertDimension (vec.size(),
                     dof_info.vector_partitioner->size());
  }

  template <typename Number>
  inline
  void check_vector_compatibility (const LinearAlgebra::distributed::Vector<Number> &vec,
                                   const internal::MatrixFreeFunctions::DoFInfo     &dof_info)
  {
    (void) vec;
    (void) dof_info;
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

    // variant where VectorType::value_type is the same as Number -> can call
    // gather
    template <typename VectorType>
    void process_dof_gather (const unsigned int      *indices,
                             VectorType              &vec,
                             VectorizedArray<Number> &res,
                             std::integral_constant<bool, true>) const
    {
      res.gather(vec.begin(), indices);
    }

    // variant where VectorType::value_type is not the same as Number -> must
    // manually load the data
    template <typename VectorType>
    void process_dof_gather (const unsigned int      *indices,
                             VectorType              &vec,
                             VectorizedArray<Number> &res,
                             std::integral_constant<bool, false>) const
    {
      for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v)
        res[v] = vector_access(const_cast<const VectorType &>(vec), indices[v]);
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

    // variant where VectorType::value_type is the same as Number -> can call
    // scatter
    template <typename VectorType>
    void process_dof_gather (const unsigned int      *indices,
                             VectorType              &vec,
                             VectorizedArray<Number> &res,
                             std::integral_constant<bool, true>) const
    {
      // TODO: enable scatter path when indices are fixed
      //#if DEAL_II_COMPILER_VECTORIZATION_LEVEL < 3
#if 1
      for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v)
        vector_access(vec, indices[v]) += res[v];
#else
      // only use gather in case there is also scatter.
      VectorizedArray<Number> tmp;
      tmp.gather(vec.begin(), indices);
      tmp += res;
      tmp.scatter(indices, vec.begin());
#endif
    }

    // variant where VectorType::value_type is not the same as Number -> must
    // manually append all data
    template <typename VectorType>
    void process_dof_gather (const unsigned int      *indices,
                             VectorType              &vec,
                             VectorizedArray<Number> &res,
                             std::integral_constant<bool, false>) const
    {
      for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v)
        vector_access(vec, indices[v]) += res[v];
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
    void process_dof_gather (const unsigned int      *indices,
                             VectorType              &vec,
                             VectorizedArray<Number> &res,
                             std::integral_constant<bool, true>) const
    {
      res.scatter(indices, vec.begin());
    }

    template <typename VectorType>
    void process_dof_gather (const unsigned int      *indices,
                             VectorType              &vec,
                             VectorizedArray<Number> &res,
                             std::integral_constant<bool, false>) const
    {
      for (unsigned int v=0; v<VectorizedArray<Number>::n_array_elements; ++v)
        vector_access(vec, indices[v]) = res[v];
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
                                                 const unsigned int component)
    {
      // FEEvaluation allows to combine several vectors from a scalar
      // FiniteElement into a "vector-valued" FEEvaluation object with
      // multiple components. These components can be extracted with the other
      // get_vector_component functions. If we do not get a vector of vectors
      // (std::vector<VectorType>, std::vector<VectorType*>, BlockVector), we
      // must make sure that we do not duplicate the components in input
      // and/or duplicate the resulting integrals. In such a case, we should
      // only get the zeroth component in the vector contained set nullptr for
      // the others which allows us to catch unintended use in
      // read_write_operation.
      if (component == 0)
        return &vec;
      else
        return nullptr;
    }
  };

  template <typename VectorType>
  struct BlockVectorSelector<std::vector<VectorType>,false>
  {
    typedef VectorType BaseVectorType;

    static BaseVectorType *get_vector_component (std::vector<VectorType> &vec,
                                                 const unsigned int component)
    {
      AssertIndexRange (component, vec.size());
      return &vec[component];
    }
  };

  template <typename VectorType>
  struct BlockVectorSelector<std::vector<VectorType *>,false>
  {
    typedef VectorType BaseVectorType;

    static BaseVectorType *get_vector_component (std::vector<VectorType *> &vec,
                                                 const unsigned int component)
    {
      AssertIndexRange (component, vec.size());
      return vec[component];
    }
  };
}



template <int dim, int n_components_, typename Number>
template <typename VectorType, typename VectorOperation>
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
  if (matrix_info == nullptr)
    {
      Assert (!local_dof_indices.empty(), ExcNotInitialized());

      unsigned int index = first_selected_component * this->data->dofs_per_component_on_cell;
      for (unsigned int comp = 0; comp<n_components; ++comp)
        {
          for (unsigned int i=0; i<this->data->dofs_per_component_on_cell; ++i, ++index)
            {
              operation.process_dof_global(local_dof_indices[this->data->lexicographic_numbering[index]],
                                           *src[0], values_dofs[comp][i][0]);
              for (unsigned int v=1; v<VectorizedArray<Number>::n_array_elements; ++v)
                operation.process_empty(values_dofs[comp][i][v]);
            }
        }
      return;
    }

  Assert (dof_info != nullptr, ExcNotInitialized());
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
  const unsigned int dofs_per_component = this->data->dofs_per_component_on_cell;

  const unsigned int n_irreg_components_filled = dof_info->row_starts[cell][2];
  const bool at_irregular_cell = n_irreg_components_filled > 0;

  // scalar case (or case when all components have the same degrees of freedom
  // and sit on a different vector each)
  if (n_fe_components == 1)
    {
      for (unsigned int c=0; c<n_components; ++c)
        Assert(src[c] != nullptr,
               ExcMessage("The finite element underlying this FEEvaluation "
                          "object is scalar, but you requested " +
                          std::to_string(n_components) +
                          " components via the template argument in "
                          "FEEvaluation. In that case, you must pass an "
                          "std::vector<VectorType> or a BlockVector to " +
                          "read_dof_values and distribute_local_to_global."));

      const unsigned int n_local_dofs =
        VectorizedArray<Number>::n_array_elements * dofs_per_component;
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
              for (unsigned int j=0, ind=0; j<dofs_per_component; ++j, ind += VectorizedArray<Number>::n_array_elements)
                for (unsigned int comp=0; comp<n_components; ++comp)
                  operation.process_dof_gather(dof_indices+ind,
                                               *src[comp], values_dofs[comp][j],
                                               std::integral_constant<bool, std::is_same<typename VectorType::value_type,Number>::value>());
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
        dofs_per_component*VectorizedArray<Number>::n_array_elements * n_components;
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
              for (unsigned int comp=0, ind=0; comp<n_components; ++comp)
                for (unsigned int j=0; j<dofs_per_component; ++j, ind += VectorizedArray<Number>::n_array_elements)
                  operation.process_dof_gather(dof_indices+ind,
                                               *src[0], values_dofs[comp][j],
                                               std::integral_constant<bool, std::is_same<typename VectorType::value_type,Number>::value>());
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
template <typename VectorType>
inline
void
FEEvaluationBase<dim,n_components_,Number>
::read_dof_values (const VectorType  &src,
                   const unsigned int first_index)
{
  // select between block vectors and non-block vectors. Note that the number
  // of components is checked in the internal data
  typename internal::BlockVectorSelector<VectorType,
           IsBlockVector<VectorType>::value>::BaseVectorType *src_data[n_components];
  for (unsigned int d=0; d<n_components; ++d)
    src_data[d] = internal::BlockVectorSelector<VectorType, IsBlockVector<VectorType>::value>::get_vector_component(const_cast<VectorType &>(src), d+first_index);

  internal::VectorReader<Number> reader;
  read_write_operation (reader, src_data);

#ifdef DEBUG
  dof_values_initialized = true;
#endif
}



template <int dim, int n_components_, typename Number>
template <typename VectorType>
inline
void
FEEvaluationBase<dim,n_components_,Number>
::read_dof_values_plain (const VectorType  &src,
                         const unsigned int first_index)
{
  // select between block vectors and non-block vectors. Note that the number
  // of components is checked in the internal data
  const typename internal::BlockVectorSelector<VectorType,
        IsBlockVector<VectorType>::value>::BaseVectorType *src_data[n_components];
  for (unsigned int d=0; d<n_components; ++d)
    src_data[d] = internal::BlockVectorSelector<VectorType, IsBlockVector<VectorType>::value>::get_vector_component(const_cast<VectorType &>(src), d+first_index);

  read_dof_values_plain (src_data);
}



template <int dim, int n_components_, typename Number>
template <typename VectorType>
inline
void
FEEvaluationBase<dim,n_components_,Number>
::distribute_local_to_global (VectorType        &dst,
                              const unsigned int first_index) const
{
  Assert (dof_values_initialized==true,
          internal::ExcAccessToUninitializedField());

  // select between block vectors and non-block vectors. Note that the number
  // of components is checked in the internal data
  typename internal::BlockVectorSelector<VectorType,
           IsBlockVector<VectorType>::value>::BaseVectorType *dst_data[n_components];
  for (unsigned int d=0; d<n_components; ++d)
    dst_data[d] = internal::BlockVectorSelector<VectorType, IsBlockVector<VectorType>::value>::get_vector_component(dst, d+first_index);

  internal::VectorDistributorLocalToGlobal<Number> distributor;
  read_write_operation (distributor, dst_data);
}



template <int dim, int n_components_, typename Number>
template <typename VectorType>
inline
void
FEEvaluationBase<dim,n_components_,Number>
::set_dof_values (VectorType        &dst,
                  const unsigned int first_index) const
{
  Assert (dof_values_initialized==true,
          internal::ExcAccessToUninitializedField());

  // select between block vectors and non-block vectors. Note that the number
  // of components is checked in the internal data
  typename internal::BlockVectorSelector<VectorType,
           IsBlockVector<VectorType>::value>::BaseVectorType *dst_data[n_components];
  for (unsigned int d=0; d<n_components; ++d)
    dst_data[d] = internal::BlockVectorSelector<VectorType, IsBlockVector<VectorType>::value>::get_vector_component(dst, d+first_index);

  internal::VectorSetter<Number> setter;
  read_write_operation (setter, dst_data);
}



template <int dim, int n_components_, typename Number>
template <typename VectorType>
inline
void
FEEvaluationBase<dim,n_components_,Number>
::read_dof_values_plain (const VectorType *src[])
{
  // Case without MatrixFree initialization object
  if (matrix_info == nullptr)
    {
      internal::VectorReader<Number> reader;
      read_write_operation (reader, src);
      return;
    }

  // this is different from the other three operations because we do not use
  // constraints here, so this is a separate function.
  Assert (dof_info != nullptr, ExcNotInitialized());
  Assert (matrix_info->indices_initialized() == true,
          ExcNotInitialized());
  Assert (cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  Assert (dof_info->store_plain_indices == true, ExcNotInitialized());

  // loop over all local dofs. ind_local holds local number on cell, index
  // iterates over the elements of index_local_to_global and dof_indices
  // points to the global indices stored in index_local_to_global
  const unsigned int *dof_indices = dof_info->begin_indices_plain(cell);
  const unsigned int dofs_per_component = this->data->dofs_per_component_on_cell;

  const unsigned int n_irreg_components_filled = dof_info->row_starts[cell][2];
  const bool at_irregular_cell = n_irreg_components_filled > 0;

  // scalar case (or case when all components have the same degrees of freedom
  // and sit on a different vector each)
  if (n_fe_components == 1)
    {
      for (unsigned int c=0; c<n_components; ++c)
        Assert(src[c] != nullptr,
               ExcMessage("The finite element underlying this FEEvaluation "
                          "object is scalar, but you requested " +
                          std::to_string(n_components) +
                          " components via the template argument in "
                          "FEEvaluation. In that case, you must pass an "
                          "std::vector<VectorType> or a BlockVector to " +
                          "read_dof_values_plain."));

      const unsigned int n_local_dofs =
        VectorizedArray<Number>::n_array_elements * dofs_per_component;
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
        dofs_per_component * VectorizedArray<Number>::n_array_elements * n_components;
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
const std::vector<unsigned int> &
FEEvaluationBase<dim,n_components,Number>::
get_internal_dof_numbering() const
{
  return data->lexicographic_numbering;
}



template <int dim, int n_components, typename Number>
inline
ArrayView<VectorizedArray<Number> >
FEEvaluationBase<dim,n_components,Number>::
get_scratch_data() const
{
  return ArrayView<VectorizedArray<Number> >(const_cast<VectorizedArray<Number> *>(scratch_data),
                                             scratch_data_array->end()-
                                             scratch_data);
}



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
  AssertIndexRange (dof, this->data->dofs_per_component_on_cell);
  Tensor<1,n_components_,VectorizedArray<Number> > return_value;
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
  Tensor<1,n_components_,VectorizedArray<Number> > return_value;
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

  Tensor<1,n_components_,Tensor<1,dim,VectorizedArray<Number> > > grad_out;

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

  Tensor<1,n_components_,Tensor<1,dim,VectorizedArray<Number> > > hessian_out;

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
  Tensor<1,n_components_,VectorizedArray<Number> > laplacian_out;
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
  AssertIndexRange (dof, this->data->dofs_per_component_on_cell);
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
  Tensor<1,n_components_,VectorizedArray<Number> > return_value;
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
template <int n_components_other>
inline
FEEvaluationAccess<dim,n_components_,Number>
::FEEvaluationAccess (const Mapping<dim>       &mapping,
                      const FiniteElement<dim> &fe,
                      const Quadrature<1>      &quadrature,
                      const UpdateFlags         update_flags,
                      const unsigned int        first_selected_component,
                      const FEEvaluationBase<dim,n_components_other,Number> *other)
  :
  FEEvaluationBase <dim,n_components_,Number>(mapping, fe, quadrature, update_flags,
                                              first_selected_component, other)
{}



template <int dim, int n_components_, typename Number>
inline
FEEvaluationAccess<dim,n_components_,Number>
::FEEvaluationAccess (const FEEvaluationAccess<dim,n_components_,Number> &other)
  :
  FEEvaluationBase <dim,n_components_,Number>(other)
{}



template <int dim, int n_components_, typename Number>
inline
FEEvaluationAccess<dim,n_components_,Number> &
FEEvaluationAccess<dim,n_components_,Number>
::operator= (const FEEvaluationAccess<dim,n_components_,Number> &other)
{
  this->FEEvaluationBase<dim,n_components_,Number>::operator=(other);
  return *this;
}



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
template <int n_components_other>
inline
FEEvaluationAccess<dim,1,Number>
::FEEvaluationAccess (const Mapping<dim>       &mapping,
                      const FiniteElement<dim> &fe,
                      const Quadrature<1>      &quadrature,
                      const UpdateFlags         update_flags,
                      const unsigned int        first_selected_component,
                      const FEEvaluationBase<dim,n_components_other,Number> *other)
  :
  FEEvaluationBase <dim,1,Number> (mapping, fe, quadrature, update_flags,
                                   first_selected_component, other)
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
FEEvaluationAccess<dim,1,Number> &
FEEvaluationAccess<dim,1,Number>
::operator= (const FEEvaluationAccess<dim,1,Number> &other)
{
  this->FEEvaluationBase<dim,1,Number>::operator=(other);
  return *this;
}



template <int dim, typename Number>
inline
VectorizedArray<Number>
FEEvaluationAccess<dim,1,Number>
::get_dof_value (const unsigned int dof) const
{
  AssertIndexRange (dof, this->data->dofs_per_component_on_cell);
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

  Tensor<1,dim,VectorizedArray<Number> > grad_out;

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
  AssertIndexRange (dof, this->data->dofs_per_component_on_cell);
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
template <int n_components_other>
inline
FEEvaluationAccess<dim,dim,Number>
::FEEvaluationAccess (const Mapping<dim>       &mapping,
                      const FiniteElement<dim> &fe,
                      const Quadrature<1>      &quadrature,
                      const UpdateFlags         update_flags,
                      const unsigned int        first_selected_component,
                      const FEEvaluationBase<dim,n_components_other,Number> *other)
  :
  FEEvaluationBase <dim,dim,Number> (mapping, fe, quadrature, update_flags,
                                     first_selected_component, other)
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
FEEvaluationAccess<dim,dim,Number> &
FEEvaluationAccess<dim,dim,Number>
::operator= (const FEEvaluationAccess<dim,dim,Number> &other)
{
  this->FEEvaluationAccess<dim,dim,Number>::operator=(other);
  return *this;
}



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
  VectorizedArray<Number> half = make_vectorized_array<Number> (0.5);
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
Tensor<1,(dim==2?1:dim),VectorizedArray<Number> >
FEEvaluationAccess<dim,dim,Number>
::get_curl (const unsigned int q_point) const
{
  // copy from generic function into dim-specialization function
  const Tensor<2,dim,VectorizedArray<Number> > grad = get_gradient(q_point);
  Tensor<1,(dim==2?1:dim),VectorizedArray<Number> > curl;
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


/*-------------------- FEEvaluationAccess scalar for 1d ----------------------------*/


template <typename Number>
inline
FEEvaluationAccess<1,1,Number>
::FEEvaluationAccess (const MatrixFree<1,Number> &data_in,
                      const unsigned int fe_no,
                      const unsigned int quad_no_in,
                      const unsigned int fe_degree,
                      const unsigned int n_q_points)
  :
  FEEvaluationBase <1,1,Number>
  (data_in, fe_no, quad_no_in, fe_degree, n_q_points)
{}



template <typename Number>
template <int n_components_other>
inline
FEEvaluationAccess<1,1,Number>
::FEEvaluationAccess (const Mapping<1>       &mapping,
                      const FiniteElement<1> &fe,
                      const Quadrature<1>    &quadrature,
                      const UpdateFlags       update_flags,
                      const unsigned int      first_selected_component,
                      const FEEvaluationBase<1,n_components_other,Number> *other)
  :
  FEEvaluationBase <1,1,Number> (mapping, fe, quadrature, update_flags,
                                 first_selected_component, other)
{}



template <typename Number>
inline
FEEvaluationAccess<1,1,Number>
::FEEvaluationAccess (const FEEvaluationAccess<1,1,Number> &other)
  :
  FEEvaluationBase <1,1,Number>(other)
{}



template <typename Number>
inline
FEEvaluationAccess<1,1,Number> &
FEEvaluationAccess<1,1,Number>
::operator= (const FEEvaluationAccess<1,1,Number> &other)
{
  this->FEEvaluationBase<1,1,Number>::operator=(other);
  return *this;
}



template <typename Number>
inline
VectorizedArray<Number>
FEEvaluationAccess<1,1,Number>
::get_dof_value (const unsigned int dof) const
{
  AssertIndexRange (dof, this->data->dofs_per_component_on_cell);
  return this->values_dofs[0][dof];
}



template <typename Number>
inline
VectorizedArray<Number>
FEEvaluationAccess<1,1,Number>
::get_value (const unsigned int q_point) const
{
  Assert (this->values_quad_initialized==true,
          internal::ExcAccessToUninitializedField());
  AssertIndexRange (q_point, this->data->n_q_points);
  return this->values_quad[0][q_point];
}



template <typename Number>
inline
Tensor<1,1,VectorizedArray<Number> >
FEEvaluationAccess<1,1,Number>
::get_gradient (const unsigned int q_point) const
{
  // could use the base class gradient, but that involves too many inefficient
  // initialization operations on tensors

  Assert (this->gradients_quad_initialized==true,
          internal::ExcAccessToUninitializedField());
  AssertIndexRange (q_point, this->data->n_q_points);

  Tensor<1,1,VectorizedArray<Number> > grad_out;

  // Cartesian cell
  if (this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      grad_out[0] = (this->gradients_quad[0][0][q_point] *
                     this->cartesian_data[0][0]);
    }
  // cell with general/constant Jacobian
  else
    {
      const Tensor<2,1,VectorizedArray<Number> > &jac =
        this->cell_type == internal::MatrixFreeFunctions::general ?
        this->jacobian[q_point] : this->jacobian[0];

      grad_out[0] = (jac[0][0] * this->gradients_quad[0][0][q_point]);
    }
  return grad_out;
}



template <typename Number>
inline
Tensor<2,1,VectorizedArray<Number> >
FEEvaluationAccess<1,1,Number>
::get_hessian (const unsigned int q_point) const
{
  return BaseClass::get_hessian(q_point)[0];
}



template <typename Number>
inline
Tensor<1,1,VectorizedArray<Number> >
FEEvaluationAccess<1,1,Number>
::get_hessian_diagonal (const unsigned int q_point) const
{
  return BaseClass::get_hessian_diagonal(q_point)[0];
}



template <typename Number>
inline
VectorizedArray<Number>
FEEvaluationAccess<1,1,Number>
::get_laplacian (const unsigned int q_point) const
{
  return BaseClass::get_laplacian(q_point)[0];
}



template <typename Number>
inline
void
FEEvaluationAccess<1,1,Number>
::submit_dof_value (const VectorizedArray<Number> val_in,
                    const unsigned int dof)
{
#ifdef DEBUG
  this->dof_values_initialized = true;
  AssertIndexRange (dof, this->data->dofs_per_component_on_cell);
#endif
  this->values_dofs[0][dof] = val_in;
}



template <typename Number>
inline
void
FEEvaluationAccess<1,1,Number>
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



template <typename Number>
inline
void
FEEvaluationAccess<1,1,Number>
::submit_gradient (const Tensor<1,1,VectorizedArray<Number> > grad_in,
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
      this->gradients_quad[0][0][q_point] = (grad_in[0] *
                                             this->cartesian_data[0][0] *
                                             JxW);
    }
  // general/affine cell type
  else
    {
      const Tensor<2,1,VectorizedArray<Number> > &jac =
        this->cell_type == internal::MatrixFreeFunctions::general ?
        this->jacobian[q_point] : this->jacobian[0];
      const VectorizedArray<Number> JxW =
        this->cell_type == internal::MatrixFreeFunctions::general ?
        this->J_value[q_point] : this->J_value[0] * this->quadrature_weights[q_point];

      this->gradients_quad[0][0][q_point] = jac[0][0] * grad_in[0] * JxW;
    }
}



template <typename Number>
inline
VectorizedArray<Number>
FEEvaluationAccess<1,1,Number>
::integrate_value () const
{
  return BaseClass::integrate_value()[0];
}




/*-------------------------- FEEvaluation -----------------------------------*/


template <int dim, int fe_degree,  int n_q_points_1d, int n_components_,
          typename Number>
inline
FEEvaluation<dim,fe_degree,n_q_points_1d,n_components_,Number>
::FEEvaluation (const MatrixFree<dim,Number> &data_in,
                const unsigned int fe_no,
                const unsigned int quad_no)
  :
  BaseClass (data_in, fe_no, quad_no, fe_degree, static_n_q_points),
  dofs_per_component (this->data->dofs_per_component_on_cell),
  dofs_per_cell (this->data->dofs_per_component_on_cell *n_components_),
  n_q_points (this->data->n_q_points)
{
  check_template_arguments(fe_no, 0);
}



template <int dim, int fe_degree,  int n_q_points_1d, int n_components_,
          typename Number>
inline
FEEvaluation<dim,fe_degree,n_q_points_1d,n_components_,Number>
::FEEvaluation (const Mapping<dim>       &mapping,
                const FiniteElement<dim> &fe,
                const Quadrature<1>      &quadrature,
                const UpdateFlags         update_flags,
                const unsigned int        first_selected_component)
  :
  BaseClass (mapping, fe, quadrature, update_flags,
             first_selected_component,
             static_cast<FEEvaluationBase<dim,1,Number>*>(nullptr)),
  dofs_per_component (this->data->dofs_per_component_on_cell),
  dofs_per_cell (this->data->dofs_per_component_on_cell *n_components_),
  n_q_points (this->data->n_q_points)
{
  check_template_arguments(numbers::invalid_unsigned_int, 0);
}



template <int dim, int fe_degree,  int n_q_points_1d, int n_components_,
          typename Number>
inline
FEEvaluation<dim,fe_degree,n_q_points_1d,n_components_,Number>
::FEEvaluation (const FiniteElement<dim> &fe,
                const Quadrature<1>      &quadrature,
                const UpdateFlags         update_flags,
                const unsigned int        first_selected_component)
  :
  BaseClass (StaticMappingQ1<dim>::mapping, fe, quadrature, update_flags,
             first_selected_component,
             static_cast<FEEvaluationBase<dim,1,Number>*>(nullptr)),
  dofs_per_component (this->data->dofs_per_component_on_cell),
  dofs_per_cell (this->data->dofs_per_component_on_cell *n_components_),
  n_q_points (this->data->n_q_points)
{
  check_template_arguments(numbers::invalid_unsigned_int, 0);
}



template <int dim, int fe_degree,  int n_q_points_1d, int n_components_,
          typename Number>
template <int n_components_other>
inline
FEEvaluation<dim,fe_degree,n_q_points_1d,n_components_,Number>
::FEEvaluation (const FiniteElement<dim> &fe,
                const FEEvaluationBase<dim,n_components_other,Number> &other,
                const unsigned int        first_selected_component)
  :
  BaseClass (other.mapped_geometry->get_fe_values().get_mapping(),
             fe, other.mapped_geometry->get_quadrature(),
             other.mapped_geometry->get_fe_values().get_update_flags(),
             first_selected_component, &other),
  dofs_per_component (this->data->dofs_per_component_on_cell),
  dofs_per_cell (this->data->dofs_per_component_on_cell *n_components_),
  n_q_points (this->data->n_q_points)
{
  check_template_arguments(numbers::invalid_unsigned_int, 0);
}



template <int dim, int fe_degree,  int n_q_points_1d, int n_components_,
          typename Number>
inline
FEEvaluation<dim,fe_degree,n_q_points_1d,n_components_,Number>
::FEEvaluation (const FEEvaluation &other)
  :
  BaseClass (other),
  dofs_per_component (this->data->dofs_per_component_on_cell),
  dofs_per_cell (this->data->dofs_per_component_on_cell *n_components_),
  n_q_points (this->data->n_q_points)
{
  check_template_arguments(numbers::invalid_unsigned_int, 0);
}



template <int dim, int fe_degree,  int n_q_points_1d, int n_components_,
          typename Number>
inline
FEEvaluation<dim,fe_degree,n_q_points_1d,n_components_,Number> &
FEEvaluation<dim,fe_degree,n_q_points_1d,n_components_,Number>
::operator= (const FEEvaluation &other)
{
  BaseClass::operator=(other);
  check_template_arguments(numbers::invalid_unsigned_int, 0);
  return *this;
}



template <int dim, int fe_degree,  int n_q_points_1d, int n_components_,
          typename Number>
inline
void
FEEvaluation<dim,fe_degree,n_q_points_1d,n_components_,Number>
::check_template_arguments(const unsigned int fe_no,
                           const unsigned int first_selected_component)
{
  (void)fe_no;
  (void)first_selected_component;

#ifdef DEBUG
  // print error message when the dimensions do not match. Propose a possible
  // fix
  if ((fe_degree != -1 && static_cast<unsigned int>(fe_degree) != this->data->fe_degree)
      ||
      (fe_degree != -1 && static_n_q_points != this->data->n_q_points))
    {
      std::string message =
        "-------------------------------------------------------\n";
      message += "Illegal arguments in constructor/wrong template arguments!\n";
      message += "    Called -->   FEEvaluation<dim,";
      message += Utilities::int_to_string(fe_degree) + ",";
      message += Utilities::int_to_string(n_q_points_1d);
      message += "," + Utilities::int_to_string(n_components);
      message += ",Number>(data";
      if (fe_no != numbers::invalid_unsigned_int)
        {
          message += ", " + Utilities::int_to_string(fe_no) + ", ";
          message += Utilities::int_to_string(this->quad_no);
        }
      message += ")\n";

      // check whether some other vector component has the correct number of
      // points
      unsigned int proposed_dof_comp = numbers::invalid_unsigned_int,
                   proposed_quad_comp = numbers::invalid_unsigned_int;
      if (fe_no != numbers::invalid_unsigned_int)
        {
          if (static_cast<unsigned int>(fe_degree) == this->data->fe_degree)
            proposed_dof_comp = fe_no;
          else
            for (unsigned int no=0; no<this->matrix_info->n_components(); ++no)
              if (this->matrix_info->get_shape_info(no,0,this->active_fe_index,0).fe_degree
                  == static_cast<unsigned int>(fe_degree))
                {
                  proposed_dof_comp = no;
                  break;
                }
          if (static_n_q_points ==
              this->mapping_info->mapping_data_gen[this->quad_no].n_q_points[this->active_quad_index])
            proposed_quad_comp = this->quad_no;
          else
            for (unsigned int no=0; no<this->mapping_info->mapping_data_gen.size(); ++no)
              if (this->mapping_info->mapping_data_gen[no].n_q_points[this->active_quad_index]
                  == static_n_q_points)
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
          message += ",Number>(data";
          if (fe_no != numbers::invalid_unsigned_int)
            {
              message += ", " + Utilities::int_to_string(proposed_dof_comp) + ", ";
              message += Utilities::int_to_string(proposed_quad_comp);
            }
          message += ")?\n";
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
      message += ",Number>(data";
      if (fe_no != numbers::invalid_unsigned_int)
        {
          message += ", " + Utilities::int_to_string(fe_no) + ", ";
          message += Utilities::int_to_string(this->quad_no);
        }
      message += ")?\n";
      std::string correct_pos;
      if (this->data->fe_degree != static_cast<unsigned int>(fe_degree))
        correct_pos = " ^";
      else
        correct_pos = "  ";
      if (proposed_n_q_points_1d != n_q_points_1d)
        correct_pos += " ^\n";
      else
        correct_pos += "  \n";
      message += "                                 " + correct_pos;

      Assert (static_cast<unsigned int>(fe_degree) == this->data->fe_degree &&
              static_n_q_points == this->data->n_q_points,
              ExcMessage(message));
    }
  if (fe_no != numbers::invalid_unsigned_int)
    {
      AssertDimension (n_q_points,
                       this->mapping_info->mapping_data_gen[this->quad_no].
                       n_q_points[this->active_quad_index]);
      AssertDimension (this->data->dofs_per_component_on_cell * this->n_fe_components,
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
  Assert (this->quadrature_points != nullptr, ExcNotInitialized());
  AssertIndexRange (q, n_q_points);

  // Cartesian mesh: not all quadrature points are stored, only the
  // diagonal. Hence, need to find the tensor product index and retrieve the
  // value from that
  if (this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      Point<dim,VectorizedArray<Number> > point;
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
::evaluate (const bool evaluate_values,
            const bool evaluate_gradients,
            const bool evaluate_hessians)
{
  Assert (this->dof_values_initialized == true,
          internal::ExcAccessToUninitializedField());
  Assert(this->matrix_info != nullptr ||
         this->mapped_geometry->is_initialized(), ExcNotInitialized());

  SelectEvaluator<dim, fe_degree, n_q_points_1d, n_components, Number>
  ::evaluate (*this->data, &this->values_dofs[0], this->values_quad,
              this->gradients_quad, this->hessians_quad, this->scratch_data,
              evaluate_values, evaluate_gradients, evaluate_hessians);

#ifdef DEBUG
  if (evaluate_values == true)
    this->values_quad_initialized = true;
  if (evaluate_gradients == true)
    this->gradients_quad_initialized = true;
  if (evaluate_hessians == true)
    this->hessians_quad_initialized  = true;
#endif
}



template <int dim, int fe_degree,  int n_q_points_1d, int n_components_,
          typename Number>
inline
void
FEEvaluation<dim,fe_degree,n_q_points_1d,n_components_,Number>
::integrate (const bool integrate_values,
             const bool integrate_gradients)
{
  if (integrate_values == true)
    Assert (this->values_quad_submitted == true,
            internal::ExcAccessToUninitializedField());
  if (integrate_gradients == true)
    Assert (this->gradients_quad_submitted == true,
            internal::ExcAccessToUninitializedField());
  Assert(this->matrix_info != nullptr ||
         this->mapped_geometry->is_initialized(), ExcNotInitialized());

  SelectEvaluator<dim, fe_degree, n_q_points_1d, n_components, Number>
  ::integrate (*this->data, &this->values_dofs[0], this->values_quad,
               this->gradients_quad, this->scratch_data,
               integrate_values, integrate_gradients);

#ifdef DEBUG
  this->dof_values_initialized = true;
#endif
}



#endif  // ifndef DOXYGEN


DEAL_II_NAMESPACE_CLOSE

#endif
