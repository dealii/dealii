// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


#ifndef dealii_matrix_free_fe_evaluation_h
#define dealii_matrix_free_fe_evaluation_h


#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/lac/vector_operation.h>

#include <deal.II/matrix_free/evaluation_kernels.h>
#include <deal.II/matrix_free/evaluation_selector.h>
#include <deal.II/matrix_free/mapping_data_on_the_fly.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/shape_info.h>
#include <deal.II/matrix_free/tensor_product_kernels.h>
#include <deal.II/matrix_free/type_traits.h>
#include <deal.II/matrix_free/vector_access_internal.h>


DEAL_II_NAMESPACE_OPEN



namespace internal
{
  DeclException0(ExcAccessToUninitializedField);
}

template <int dim,
          int fe_degree,
          int n_q_points_1d            = fe_degree + 1,
          int n_components_            = 1,
          typename Number              = double,
          typename VectorizedArrayType = VectorizedArray<Number>>
class FEEvaluation;


/**
 * This is the base class for the FEEvaluation classes. This class is a base
 * class and needs usually not be called in user code. It does not have any
 * public constructor. The usage is through the class FEEvaluation instead. It
 * implements a reinit method that is used to set pointers so that operations
 * on quadrature points can be performed quickly, access functions to vectors
 * for the @p read_dof_values, @p set_dof_values, and @p
 * distribute_local_to_global functions, as well as methods to access values
 * and gradients of finite element functions.
 *
 * This class has three template arguments:
 *
 * @tparam dim Dimension in which this class is to be used
 *
 * @tparam n_components Number of vector components when solving a system of
 * PDEs. If the same operation is applied to several components of a PDE (e.g.
 * a vector Laplace equation), they can be applied simultaneously with one
 * call (and often more efficiently)
 *
 * @tparam Number Number format, usually @p double or @p float
 *
 * @tparam VectorizedArrayType Type of array to be woked on in a vectorized
 *                             fashion, defaults to VectorizedArray<Number>
 *
 * @note Currently only VectorizedArray<Number, width> is supported as
 *       VectorizedArrayType.
 *
 *
 * @ingroup matrixfree
 *
 * @author Katharina Kormann and Martin Kronbichler, 2010-2018
 */
template <int dim,
          int n_components_,
          typename Number,
          bool is_face                 = false,
          typename VectorizedArrayType = VectorizedArray<Number>>
class FEEvaluationBase
{
  static_assert(
    std::is_same<Number, typename VectorizedArrayType::value_type>::value,
    "Type of Number and of VectorizedArrayType do not match.");

public:
  using number_type = Number;
  using value_type  = Tensor<1, n_components_, VectorizedArrayType>;
  using gradient_type =
    Tensor<1, n_components_, Tensor<1, dim, VectorizedArrayType>>;
  static constexpr unsigned int dimension    = dim;
  static constexpr unsigned int n_components = n_components_;

  /**
   * @name 1: General operations
   */
  //@{
  /**
   * Destructor.
   */
  ~FEEvaluationBase();

  /**
   * @deprecated Use get_mapping_data_index_offset() instead.
   */
  DEAL_II_DEPRECATED
  unsigned int
  get_cell_data_number() const;

  /**
   * Return the index offset within the geometry fields for the cell the @p
   * reinit() function has been called for. This index can be used to access
   * an index into a field that has the same compression behavior as the
   * Jacobian of the geometry, e.g., to store an effective coefficient tensors
   * that combines a coefficient with the geometry for lower memory transfer
   * as the available data fields.
   */
  unsigned int
  get_mapping_data_index_offset() const;

  /**
   * Return the type of the cell the @p reinit() function has been called for.
   * Valid values are @p cartesian for Cartesian cells (which allows for
   * considerable data compression), @p affine for cells with affine mappings,
   * and @p general for general cells without any compressed storage applied.
   */
  internal::MatrixFreeFunctions::GeometryType
  get_cell_type() const;

  /**
   * Return a reference to the ShapeInfo object currently in use.
   */
  const internal::MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> &
  get_shape_info() const;

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
   * see it as a similar function to AffineConstraints::read_dof_values as
   * well. Note that if vectorization is enabled, the DoF values for several
   * cells are set.
   *
   * If some constraints on the vector are inhomogeneous, use the function
   * read_dof_values_plain instead and provide the vector with useful data
   * also in constrained positions by calling AffineConstraints::distribute.
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
  void
  read_dof_values(const VectorType &src, const unsigned int first_index = 0);

  /**
   * For the vector @p src, read out the values on the degrees of freedom of
   * the current cell, and store them internally. Similar functionality as the
   * function DoFAccessor::get_interpolated_dof_values. As opposed to the
   * read_dof_values function, this function reads out the plain entries from
   * vectors, without taking stored constraints into account. This way of
   * access is appropriate when the constraints have been distributed on the
   * vector by a call to AffineConstraints::distribute previously. This
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
  void
  read_dof_values_plain(const VectorType & src,
                        const unsigned int first_index = 0);

  /**
   * Takes the values stored internally on dof values of the current cell and
   * sums them into the vector @p dst. The function also applies constraints
   * during the write operation. The functionality is hence similar to the
   * function AffineConstraints::distribute_local_to_global. If vectorization
   * is enabled, the DoF values for several cells are used. The mask can be
   * used to suppress the write access for some of the cells contained in the
   * current cell vectorization batch, e.g. in case of local time stepping,
   * where some cells are excluded from a call. A value of `true` in the
   * bitset means that the respective lane index will be processed, whereas a
   * value of `false` skips this index. The default setting is a bitset that
   * contains all ones, which will write the accumulated integrals to all
   * cells in the batch.
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
  void
  distribute_local_to_global(
    VectorType &                                              dst,
    const unsigned int                                        first_index = 0,
    const std::bitset<VectorizedArrayType::n_array_elements> &mask =
      std::bitset<VectorizedArrayType::n_array_elements>().flip()) const;

  /**
   * Takes the values stored internally on dof values of the current cell and
   * writes them into the vector @p dst. The function skips the degrees of
   * freedom which are constrained. As opposed to the
   * distribute_local_to_global method, the old values at the position given
   * by the current cell are overwritten. Thus, if a degree of freedom is
   * associated to more than one cell (as usual in continuous finite
   * elements), the values will be overwritten and only the value written last
   * is retained. The mask can be used to suppress the write access for some
   * of the cells contained in the current cell vectorization batch, e.g. in
   * case of local time stepping, where some  cells are excluded from a call.
   * A value of `true` in the bitset means that the respective lane index will
   * be processed, whereas a value of `false` skips this index. The default
   * setting is a bitset that contains all ones, which will write the
   * accumulated integrals to all cells in the batch.
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
  void
  set_dof_values(
    VectorType &                                              dst,
    const unsigned int                                        first_index = 0,
    const std::bitset<VectorizedArrayType::n_array_elements> &mask =
      std::bitset<VectorizedArrayType::n_array_elements>().flip()) const;

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
  value_type
  get_dof_value(const unsigned int dof) const;

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
  void
  submit_dof_value(const value_type val_in, const unsigned int dof);

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
  value_type
  get_value(const unsigned int q_point) const;

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
  void
  submit_value(const value_type val_in, const unsigned int q_point);

  /**
   * Return the gradient of a finite element function at quadrature point
   * number @p q_point after a call to @p evaluate(...,true,...), or the value
   * that has been stored there with a call to @p submit_gradient.
   *
   * Note that the derived class FEEvaluationAccess overloads this operation
   * with specializations for the scalar case (n_components == 1) and for the
   * vector-valued case (n_components == dim).
   */
  gradient_type
  get_gradient(const unsigned int q_point) const;

  /**
   * Return the derivative of a finite element function at quadrature point
   * number @p q_point after a call to @p evaluate(...,true,...) in the
   * direction normal to the face:
   * $\boldsymbol \nabla u(\mathbf x_q) \cdot \mathbf n(\mathbf x_q)$
   *
   * This call is equivalent to calling `get_gradient() * get_normal_vector()`
   * but will use a more efficient internal representation of data.
   *
   * Note that the derived class FEEvaluationAccess overloads this operation
   * with specializations for the scalar case (n_components == 1) and for the
   * vector-valued case (n_components == dim).
   */
  value_type
  get_normal_derivative(const unsigned int q_point) const;

  /**
   * Write a contribution that is tested by the gradient to the field
   * containing the values on quadrature points with component @p q_point.
   * Access to the same field as through get_gradient(). If applied before
   * the function @p integrate(...,true) is called, this specifies what is
   * tested by all basis function gradients on the current cell and integrated
   * over.
   *
   * Note that the derived class FEEvaluationAccess overloads this operation
   * with specializations for the scalar case (n_components == 1) and for the
   * vector-valued case (n_components == dim).
   */
  void
  submit_gradient(const gradient_type grad_in, const unsigned int q_point);

  /**
   * Write a contribution that is tested by the gradient to the field
   * containing the values on quadrature points with component @p
   * q_point. Access to the same field as through get_gradient() or
   * get_normal_derivative(). If applied before the function @p
   * integrate(...,true) is called, this specifies what is tested by all basis
   * function gradients on the current cell and integrated over.
   *
   * @note This operation writes the data to the same field as
   * submit_gradient(). As a consequence, only one of these two can be
   * used. Usually, the contribution of a potential call to this function must
   * be added into the contribution for submit_gradient().
   *
   * @note The derived class FEEvaluationAccess overloads this operation
   * with specializations for the scalar case (n_components == 1) and for the
   * vector-valued case (n_components == dim).
   */
  void
  submit_normal_derivative(const value_type   grad_in,
                           const unsigned int q_point);

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
  Tensor<1, n_components_, Tensor<2, dim, VectorizedArrayType>>
  get_hessian(const unsigned int q_point) const;

  /**
   * Return the diagonal of the Hessian of a finite element function at
   * quadrature point number @p q_point after a call to @p evaluate(...,true).
   *
   * Note that the derived class FEEvaluationAccess overloads this operation
   * with specializations for the scalar case (n_components == 1) and for the
   * vector-valued case (n_components == dim).
   */
  gradient_type
  get_hessian_diagonal(const unsigned int q_point) const;

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
  value_type
  get_laplacian(const unsigned int q_point) const;

#ifdef DOXYGEN
  // doxygen does not anyhow mention functions coming from partial template
  // specialization of the base class, in this case FEEvaluationAccess<dim,dim>.
  // For now, hack-in those functions manually only to fix documentation:

  /** @copydoc FEEvaluationAccess<dim,dim,Number,is_face>::get_divergence()
   * @note Only available for n_components_==dim.
   */
  VectorizedArrayType
  get_divergence(const unsigned int q_point) const;

  /** @copydoc FEEvaluationAccess<dim,dim,Number,is_face>::get_symmetric_gradient()
   * @note Only available for n_components_==dim.
   */
  SymmetricTensor<2, dim, VectorizedArrayType>
  get_symmetric_gradient(const unsigned int q_point) const;

  /** @copydoc FEEvaluationAccess<dim,dim,Number,is_face>::get_curl()
   * @note Only available for n_components_==dim.
   */
  Tensor<1, (dim == 2 ? 1 : dim), VectorizedArrayType>
  get_curl(const unsigned int q_point) const;

  /** @copydoc FEEvaluationAccess<dim,dim,Number,is_face>::submit_divergence()
   * @note Only available for n_components_==dim.
   *
   * @note This operation writes the data to the same field as
   * submit_gradient(). As a consequence, only one of these two can be
   * used. Usually, the contribution of a potential call to this function must
   * be added into the diagonal of the contribution for submit_gradient().
   */
  void
  submit_divergence(const VectorizedArrayType div_in,
                    const unsigned int        q_point);

  /** @copydoc FEEvaluationAccess<dim,dim,Number,is_face>::submit_symmetric_gradient()
   * @note Only available for n_components_==dim.
   *
   * @note This operation writes the data to the same field as
   * submit_gradient(). As a consequence, only one of these two can be
   * used. Usually, the contribution of a potential call to this function must
   * be added to the respective entries of the rank-2 tensor for
   * submit_gradient().
   */
  void
  submit_symmetric_gradient(
    const SymmetricTensor<2, dim, VectorizedArrayType> grad_in,
    const unsigned int                                 q_point);

  /** @copydoc FEEvaluationAccess<dim,dim,Number,is_face>::submit_curl()
   * @note Only available for n_components_==dim.
   *
   * @note This operation writes the data to the same field as
   * submit_gradient(). As a consequence, only one of these two can be
   * used. Usually, the contribution of a potential call to this function must
   * be added to the respective entries of the rank-2 tensor for
   * submit_gradient().
   */
  void
  submit_curl(const Tensor<1, dim == 2 ? 1 : dim, VectorizedArrayType> curl_in,
              const unsigned int                                       q_point);

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
  value_type
  integrate_value() const;

  /**
   * Return the determinant of the Jacobian from the unit to the real cell
   * times the quadrature weight.
   */
  VectorizedArrayType
  JxW(const unsigned int q_index) const;

  /**
   * Fills the JxW values currently used into the given array.
   *
   * @deprecated Use JxW() instead.
   */
  DEAL_II_DEPRECATED void
  fill_JxW_values(AlignedVector<VectorizedArrayType> &JxW_values) const;

  /**
   * Return the inverse and transposed version of Jacobian of the mapping
   * between the unit to the real cell (representing the covariant
   * transformation). This is exactly the matrix used internally to transform
   * the unit cell gradients to gradients on the real cell.
   */
  Tensor<2, dim, VectorizedArrayType>
  inverse_jacobian(const unsigned int q_index) const;

  /**
   * Return the unit normal vector on a face. Note that both sides of a face
   * use the same orientation of the normal vector: For the faces enumerated
   * as `interior` in FaceToCellTopology and selected with the
   * `is_interior_face=true` flag of the constructor, this corresponds to the
   * outer normal vector, whereas for faces enumerated as `exterior` in
   * FaceToCellTopology and selected with the `is_interior_face=false` flag of
   * the constructor, the normal points into the element as a consequence of
   * the single normal vector.
   *
   * @note Only implemented in case `is_face == true`.
   */
  Tensor<1, dim, VectorizedArrayType>
  get_normal_vector(const unsigned int q_point) const;

  /**
   * Provides a unified interface to access data in a vector of
   * VectorizedArray fields of length MatrixFree::n_macro_cells() +
   * MatrixFree::n_macro_ghost_cells() for both cells (plain read) and faces
   * (indirect addressing).
   */
  VectorizedArrayType
  read_cell_data(const AlignedVector<VectorizedArrayType> &array) const;

  //@}

  /**
   * @name 4: Access to internal data
   */
  //@{
  /**
   * Return a read-only pointer to the first field of the dof values. This is
   * the data field the read_dof_values() functions write into. First come the
   * dof values for the first component, then all values for the second
   * component, and so on. This is related to the internal data structures
   * used in this class. In general, it is safer to use the get_dof_value()
   * function instead.
   */
  const VectorizedArrayType *
  begin_dof_values() const;

  /**
   * Return a read and write pointer to the first field of the dof values.
   * This is the data field the read_dof_values() functions write into. First
   * come the dof values for the first component, then all values for the
   * second component, and so on. This is related to the internal data
   * structures used in this class. In general, it is safer to use the
   * get_dof_value() function instead.
   */
  VectorizedArrayType *
  begin_dof_values();

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
  const VectorizedArrayType *
  begin_values() const;

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
  VectorizedArrayType *
  begin_values();

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
  const VectorizedArrayType *
  begin_gradients() const;

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
  VectorizedArrayType *
  begin_gradients();

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
  const VectorizedArrayType *
  begin_hessians() const;

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
  VectorizedArrayType *
  begin_hessians();

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
  ArrayView<VectorizedArrayType>
  get_scratch_data() const;

  //@}

protected:
  /**
   * Constructor. Made protected to prevent users from directly using this
   * class. Takes all data stored in MatrixFree. If applied to problems with
   * more than one finite element or more than one quadrature formula selected
   * during construction of @p matrix_free, @p dof_no, @p
   * first_selected_component and @p quad_no allow to select the appropriate
   * components.
   */
  FEEvaluationBase(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    const unsigned int                                  dof_no,
    const unsigned int first_selected_component,
    const unsigned int quad_no,
    const unsigned int fe_degree,
    const unsigned int n_q_points,
    const bool         is_interior_face);

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
  FEEvaluationBase(const Mapping<dim> &      mapping,
                   const FiniteElement<dim> &fe,
                   const Quadrature<1> &     quadrature,
                   const UpdateFlags         update_flags,
                   const unsigned int        first_selected_component,
                   const FEEvaluationBase<dim,
                                          n_components_other,
                                          Number,
                                          is_face,
                                          VectorizedArrayType> *other);

  /**
   * Copy constructor. If FEEvaluationBase was constructed from a mapping, fe,
   * quadrature, and update flags, the underlying geometry evaluation based on
   * FEValues will be deep-copied in order to allow for using in parallel with
   * threads.
   */
  FEEvaluationBase(const FEEvaluationBase &other);

  /**
   * Copy assignment operator. If FEEvaluationBase was constructed from a
   * mapping, fe, quadrature, and update flags, the underlying geometry
   * evaluation based on FEValues will be deep-copied in order to allow for
   * using in parallel with threads.
   */
  FEEvaluationBase &
  operator=(const FEEvaluationBase &other);

  /**
   * A unified function to read from and write into vectors based on the given
   * template operation. It can perform the operation for @p read_dof_values,
   * @p distribute_local_to_global, and @p set_dof_values. It performs the
   * operation for several vectors at a time.
   */
  template <typename VectorType, typename VectorOperation>
  void
  read_write_operation(
    const VectorOperation &                                   operation,
    VectorType *                                              vectors[],
    const std::bitset<VectorizedArrayType::n_array_elements> &mask,
    const bool apply_constraints = true) const;

  /**
   * A unified function to read from and write into vectors based on the given
   * template operation for DG-type schemes where all degrees of freedom on
   * cells are contiguous. It can perform the operation for read_dof_values(),
   * distribute_local_to_global(), and set_dof_values() for several vectors at
   * a time, depending on n_components.
   */
  template <typename VectorType, typename VectorOperation>
  void
  read_write_operation_contiguous(
    const VectorOperation &                                   operation,
    VectorType *                                              vectors[],
    const std::bitset<VectorizedArrayType::n_array_elements> &mask) const;

  /**
   * A unified function to read from and write into vectors based on the given
   * template operation for the case when we do not have an underlying
   * MatrixFree object. It can perform the operation for @p read_dof_values,
   * @p distribute_local_to_global, and @p set_dof_values. It performs the
   * operation for several vectors at a time, depending on n_components.
   */
  template <typename VectorType, typename VectorOperation>
  void
  read_write_operation_global(const VectorOperation &operation,
                              VectorType *           vectors[]) const;

  /**
   * This is the general array for all data fields.
   */
  AlignedVector<VectorizedArrayType> *scratch_data_array;

  /**
   * This is the user-visible part of scratch_data_array, only showing the
   * last part of scratch_data_array. The first part is consumed by
   * values_dofs, values_quad, etc.
   */
  VectorizedArrayType *scratch_data;

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
  VectorizedArrayType *values_dofs[n_components];

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
  VectorizedArrayType *values_quad[n_components];

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
  VectorizedArrayType *gradients_quad[n_components][dim];

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
  VectorizedArrayType *hessians_quad[n_components][(dim * (dim + 1)) / 2];

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
   * Stores the number of quadrature points in the current evaluation context.
   */
  const unsigned int n_quadrature_points;

  /**
   * Stores a pointer to the underlying data.
   */
  const MatrixFree<dim, Number, VectorizedArrayType> *matrix_info;

  /**
   * Stores a pointer to the underlying DoF indices and constraint description
   * for the component specified at construction. Also contained in
   * matrix_info, but it simplifies code if we store a reference to it.
   */
  const internal::MatrixFreeFunctions::DoFInfo *dof_info;

  /**
   * Stores a pointer to the underlying transformation data from unit to
   * real cells for the given quadrature formula specified at construction.
   * Also contained in matrix_info, but it simplifies code if we store a
   * reference to it.
   */
  const internal::MatrixFreeFunctions::MappingInfoStorage<
    (is_face ? dim - 1 : dim),
    dim,
    Number,
    VectorizedArrayType> *mapping_data;

  /**
   * Stores a pointer to the unit cell shape data, i.e., values, gradients and
   * Hessians in 1D at the quadrature points that constitute the tensor
   * product. Also contained in matrix_info, but it simplifies code if we
   * store a reference to it.
   */
  const internal::MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> *data;

  /**
   * A pointer to the Jacobian information of the present cell. Only set to a
   * useful value if on a non-Cartesian cell.
   */
  const Tensor<2, dim, VectorizedArrayType> *jacobian;

  /**
   * A pointer to the Jacobian determinant of the present cell. If on a
   * Cartesian cell or on a cell with constant Jacobian, this is just the
   * Jacobian determinant, otherwise the Jacobian determinant times the
   * quadrature weight.
   */
  const VectorizedArrayType *J_value;

  /**
   * A pointer to the normal vectors at faces.
   */
  const Tensor<1, dim, VectorizedArrayType> *normal_vectors;

  /**
   * A pointer to the normal vectors times the jacobian at faces.
   */
  const Tensor<1, dim, VectorizedArrayType> *normal_x_jacobian;

  /**
   * A pointer to the quadrature weights of the underlying quadrature formula.
   */
  const Number *quadrature_weights;

  /**
   * After a call to reinit(), stores the number of the cell we are currently
   * working with.
   */
  unsigned int cell;

  /**
   * Flag holding information whether a face is an interior or exterior face
   * according to the defined direction of the normal.  Not used for cells.
   */
  bool is_interior_face;

  /**
   * Stores the index an FEFaceEvaluation object is currently pointing into
   * (interior face, exterior face, data associated with cell).
   */
  internal::MatrixFreeFunctions::DoFInfo::DoFAccessIndex dof_access_index;

  /**
   * Stores the current number of a face within the given cell in case
   * `is_face==true`, using values between `0` and `2*dim`.
   */
  unsigned int face_no;

  /**
   * Stores the orientation of the given face with respect to the standard
   * orientation, 0 if in standard orientation.
   */
  unsigned int face_orientation;

  /**
   * Stores the subface index of the given face. Usually, this variable takes
   * the value numbers::invalid_unsigned_int to indicate integration over the
   * full face, but in case the current physical face has a neighbor that is
   * more refined, it is a subface and must scale the entries in ShapeInfo
   * appropriately.
   */
  unsigned int subface_index;

  /**
   * Stores the type of the cell we are currently working with after a call to
   * reinit(). Valid values are @p cartesian, @p affine and @p general, which
   * have different implications on how the Jacobian transformations are
   * stored internally in MappingInfo.
   */
  internal::MatrixFreeFunctions::GeometryType cell_type;

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
  std::shared_ptr<internal::MatrixFreeFunctions::
                    MappingDataOnTheFly<dim, Number, VectorizedArrayType>>
    mapped_geometry;

  /**
   * For a FiniteElement with more than one base element, select at which
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
  void
  set_data_pointers();

  // Make other FEEvaluationBase as well as FEEvaluation objects friends.
  template <int, int, typename, bool, typename>
  friend class FEEvaluationBase;
  template <int, int, int, int, typename, typename>
  friend class FEEvaluation;
};



/**
 * This class provides access to the data fields of the FEEvaluation classes.
 * Generic access is achieved through the base class, and specializations for
 * scalar and vector-valued elements are defined separately.
 *
 * @ingroup matrixfree
 *
 * @author Katharina Kormann and Martin Kronbichler, 2010, 2011
 */
template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType = VectorizedArray<Number>>
class FEEvaluationAccess : public FEEvaluationBase<dim,
                                                   n_components_,
                                                   Number,
                                                   is_face,
                                                   VectorizedArrayType>
{
  static_assert(
    std::is_same<Number, typename VectorizedArrayType::value_type>::value,
    "Type of Number and of VectorizedArrayType do not match.");

public:
  using number_type = Number;
  using value_type  = Tensor<1, n_components_, VectorizedArrayType>;
  using gradient_type =
    Tensor<1, n_components_, Tensor<1, dim, VectorizedArrayType>>;
  static constexpr unsigned int dimension    = dim;
  static constexpr unsigned int n_components = n_components_;
  using BaseClass =
    FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>;

protected:
  /**
   * Constructor. Made protected to prevent initialization in user code. Takes
   * all data stored in MatrixFree. If applied to problems with more than one
   * finite element or more than one quadrature formula selected during
   * construction of @p matrix_free, @p first_selected_component and @p
   * quad_no allow to select the appropriate components.
   */
  FEEvaluationAccess(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    const unsigned int                                  dof_no,
    const unsigned int first_selected_component,
    const unsigned int quad_no,
    const unsigned int fe_degree,
    const unsigned int n_q_points,
    const bool         is_interior_face = true);

  /**
   * Constructor with reduced functionality for similar usage of FEEvaluation
   * as FEValues, including matrix assembly.
   */
  template <int n_components_other>
  FEEvaluationAccess(const Mapping<dim> &      mapping,
                     const FiniteElement<dim> &fe,
                     const Quadrature<1> &     quadrature,
                     const UpdateFlags         update_flags,
                     const unsigned int        first_selected_component,
                     const FEEvaluationBase<dim,
                                            n_components_other,
                                            Number,
                                            is_face,
                                            VectorizedArrayType> *other);

  /**
   * Copy constructor
   */
  FEEvaluationAccess(const FEEvaluationAccess &other);

  /**
   * Copy assignment operator
   */
  FEEvaluationAccess &
  operator=(const FEEvaluationAccess &other);
};



/**
 * This class provides access to the data fields of the FEEvaluation classes.
 * Partial specialization for scalar fields that defines access with simple
 * data fields, i.e., scalars for the values and Tensor<1,dim> for the
 * gradients.
 *
 * @ingroup matrixfree
 *
 * @author Katharina Kormann and Martin Kronbichler, 2010, 2011
 */
template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
class FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>
  : public FEEvaluationBase<dim, 1, Number, is_face, VectorizedArrayType>
{
  static_assert(
    std::is_same<Number, typename VectorizedArrayType::value_type>::value,
    "Type of Number and of VectorizedArrayType do not match.");

public:
  using number_type                       = Number;
  using value_type                        = VectorizedArrayType;
  using gradient_type                     = Tensor<1, dim, VectorizedArrayType>;
  static constexpr unsigned int dimension = dim;
  using BaseClass =
    FEEvaluationBase<dim, 1, Number, is_face, VectorizedArrayType>;

  /** @copydoc FEEvaluationBase<dim,1,Number,is_face>::get_dof_value()
   */
  value_type
  get_dof_value(const unsigned int dof) const;

  /** @copydoc FEEvaluationBase<dim,1,Number,is_face>::submit_dof_value()
   */
  void
  submit_dof_value(const value_type val_in, const unsigned int dof);

  /** @copydoc FEEvaluationBase<dim,1,Number,is_face>::get_value()
   */
  value_type
  get_value(const unsigned int q_point) const;

  /** @copydoc FEEvaluationBase<dim,1,Number,is_face>::submit_value()
   */
  void
  submit_value(const value_type val_in, const unsigned int q_point);

  /** @copydoc FEEvaluationBase<dim,1,Number,is_face>::submit_value()
   */
  void
  submit_value(const Tensor<1, 1, VectorizedArrayType> val_in,
               const unsigned int                      q_point);

  /** @copydoc FEEvaluationBase<dim,1,Number,is_face>::get_gradient()
   */
  gradient_type
  get_gradient(const unsigned int q_point) const;

  /** @copydoc FEEvaluationBase<dim,1,Number,is_face>::get_normal_derivative()
   */
  value_type
  get_normal_derivative(const unsigned int q_point) const;

  /** @copydoc FEEvaluationBase<dim,1,Number,is_face>::submit_gradient()
   */
  void
  submit_gradient(const gradient_type grad_in, const unsigned int q_point);

  /** @copydoc FEEvaluationBase<dim,1,Number,is_face>::submit_normal_derivative()
   */
  void
  submit_normal_derivative(const value_type   grad_in,
                           const unsigned int q_point);

  /** @copydoc FEEvaluationBase<dim,1,Number,is_face>::get_hessian()
   */
  Tensor<2, dim, VectorizedArrayType>
  get_hessian(unsigned int q_point) const;

  /** @copydoc FEEvaluationBase<dim,1,Number,is_face>::get_hessian_diagonal()
   */
  gradient_type
  get_hessian_diagonal(const unsigned int q_point) const;

  /** @copydoc FEEvaluationBase<dim,1,Number,is_face>::get_laplacian()
   */
  value_type
  get_laplacian(const unsigned int q_point) const;

  /** @copydoc FEEvaluationBase<dim,1,Number,is_face>::integrate_value()
   */
  value_type
  integrate_value() const;

protected:
  /**
   * Constructor. Made protected to avoid initialization in user code. Takes
   * all data stored in MatrixFree. If applied to problems with more than one
   * finite element or more than one quadrature formula selected during
   * construction of @p matrix_free, @p first_selected_component and @p
   * quad_no allow to select the appropriate components.
   */
  FEEvaluationAccess(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    const unsigned int                                  dof_no,
    const unsigned int first_selected_component,
    const unsigned int quad_no,
    const unsigned int fe_degree,
    const unsigned int n_q_points,
    const bool         is_interior_face = true);

  /**
   * Constructor with reduced functionality for similar usage of FEEvaluation
   * as FEValues, including matrix assembly.
   */
  template <int n_components_other>
  FEEvaluationAccess(const Mapping<dim> &      mapping,
                     const FiniteElement<dim> &fe,
                     const Quadrature<1> &     quadrature,
                     const UpdateFlags         update_flags,
                     const unsigned int        first_selected_component,
                     const FEEvaluationBase<dim,
                                            n_components_other,
                                            Number,
                                            is_face,
                                            VectorizedArrayType> *other);

  /**
   * Copy constructor
   */
  FEEvaluationAccess(const FEEvaluationAccess &other);

  /**
   * Copy assignment operator
   */
  FEEvaluationAccess &
  operator=(const FEEvaluationAccess &other);
};



/**
 * This class provides access to the data fields of the FEEvaluation classes.
 * Partial specialization for fields with as many components as the underlying
 * space dimension, i.e., values are of type Tensor<1,dim> and gradients of
 * type Tensor<2,dim>. Provides some additional functions for access, like the
 * symmetric gradient and divergence.
 *
 * @ingroup matrixfree
 *
 * @author Katharina Kormann and Martin Kronbichler, 2010, 2011
 */
template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
class FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType>
  : public FEEvaluationBase<dim, dim, Number, is_face, VectorizedArrayType>
{
  static_assert(
    std::is_same<Number, typename VectorizedArrayType::value_type>::value,
    "Type of Number and of VectorizedArrayType do not match.");

public:
  using number_type                       = Number;
  using value_type                        = Tensor<1, dim, VectorizedArrayType>;
  using gradient_type                     = Tensor<2, dim, VectorizedArrayType>;
  static constexpr unsigned int dimension = dim;
  static constexpr unsigned int n_components = dim;
  using BaseClass =
    FEEvaluationBase<dim, dim, Number, is_face, VectorizedArrayType>;

  /** @copydoc FEEvaluationBase<dim,dim,Number,is_face>::get_gradient()
   */
  gradient_type
  get_gradient(const unsigned int q_point) const;

  /**
   * Return the divergence of a vector-valued finite element at quadrature
   * point number @p q_point after a call to @p evaluate(...,true,...).
   */
  VectorizedArrayType
  get_divergence(const unsigned int q_point) const;

  /**
   * Return the symmetric gradient of a vector-valued finite element at
   * quadrature point number @p q_point after a call to @p
   * evaluate(...,true,...). It corresponds to <tt>0.5
   * (grad+grad<sup>T</sup>)</tt>.
   */
  SymmetricTensor<2, dim, VectorizedArrayType>
  get_symmetric_gradient(const unsigned int q_point) const;

  /**
   * Return the curl of the vector field, $\nabla \times v$ after a call to @p
   * evaluate(...,true,...).
   */
  Tensor<1, (dim == 2 ? 1 : dim), VectorizedArrayType>
  get_curl(const unsigned int q_point) const;

  /** @copydoc FEEvaluationBase<dim,dim,Number,is_face>::get_hessian()
   */
  Tensor<3, dim, VectorizedArrayType>
  get_hessian(const unsigned int q_point) const;

  /** @copydoc FEEvaluationBase<dim,dim,Number,is_face>::get_hessian_diagonal()
   */
  gradient_type
  get_hessian_diagonal(const unsigned int q_point) const;

  /** @copydoc FEEvaluationBase<dim,dim,Number,is_face>::submit_gradient()
   */
  void
  submit_gradient(const gradient_type grad_in, const unsigned int q_point);

  /**
   * Write a contribution that is tested by the gradient to the field
   * containing the values on quadrature points with component @p q_point.
   * This function is an alternative to the other submit_gradient function
   * when using a system of fixed number of equations which happens to
   * coincide with the dimension for some dimensions, but not all. To allow
   * for dimension-independent programming, this function can be used instead.
   */
  void
  submit_gradient(
    const Tensor<1, dim, Tensor<1, dim, VectorizedArrayType>> grad_in,
    const unsigned int                                        q_point);

  /**
   * Write a contribution that is tested by the divergence to the field
   * containing the values on quadrature points with component @p q_point.
   * Access to the same field as through @p get_gradient. If applied before
   * the function @p integrate(...,true) is called, this specifies what is
   * tested by all basis function gradients on the current cell and integrated
   * over.
   */
  void
  submit_divergence(const VectorizedArrayType div_in,
                    const unsigned int        q_point);

  /**
   * Write a contribution that is tested by the symmetric gradient to the field
   * containing the values on quadrature points with component @p q_point.
   * Access to the same field as through @p get_symmetric_gradient. If applied before
   * the function @p integrate(...,true) is called, this specifies the
   * symmetric gradient which is tested by all basis function symmetric
   * gradients on the current cell and integrated over.
   */
  void
  submit_symmetric_gradient(
    const SymmetricTensor<2, dim, VectorizedArrayType> grad_in,
    const unsigned int                                 q_point);

  /**
   * Write the components of a curl containing the values on quadrature point
   * @p q_point. Access to the same data field as through @p get_gradient.
   */
  void
  submit_curl(const Tensor<1, dim == 2 ? 1 : dim, VectorizedArrayType> curl_in,
              const unsigned int                                       q_point);

protected:
  /**
   * Constructor. Made protected to avoid initialization in user code. Takes
   * all data stored in MatrixFree. If applied to problems with more than one
   * finite element or more than one quadrature formula selected during
   * construction of @p matrix_free, @p first_selected_component and @p
   * quad_no allow to select the appropriate components.
   */
  FEEvaluationAccess(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    const unsigned int                                  dof_no,
    const unsigned int first_selected_component,
    const unsigned int quad_no,
    const unsigned int dofs_per_cell,
    const unsigned int n_q_points,
    const bool         is_interior_face = true);

  /**
   * Constructor with reduced functionality for similar usage of FEEvaluation
   * as FEValues, including matrix assembly.
   */
  template <int n_components_other>
  FEEvaluationAccess(const Mapping<dim> &      mapping,
                     const FiniteElement<dim> &fe,
                     const Quadrature<1> &     quadrature,
                     const UpdateFlags         update_flags,
                     const unsigned int        first_selected_component,
                     const FEEvaluationBase<dim,
                                            n_components_other,
                                            Number,
                                            is_face,
                                            VectorizedArrayType> *other);

  /**
   * Copy constructor
   */
  FEEvaluationAccess(const FEEvaluationAccess &other);

  /**
   * Copy assignment operator
   */
  FEEvaluationAccess &
  operator=(const FEEvaluationAccess &other);
};


/**
 * This class provides access to the data fields of the FEEvaluation classes.
 * Partial specialization for scalar fields in 1d that defines access with
 * simple data fields, i.e., scalars for the values and Tensor<1,1> for the
 * gradients.
 *
 * @ingroup matrixfree
 *
 * @author Katharina Kormann and Martin Kronbichler, 2010, 2011, Shiva
 * Rudraraju, 2014
 */
template <typename Number, bool is_face, typename VectorizedArrayType>
class FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>
  : public FEEvaluationBase<1, 1, Number, is_face, VectorizedArrayType>
{
  static_assert(
    std::is_same<Number, typename VectorizedArrayType::value_type>::value,
    "Type of Number and of VectorizedArrayType do not match.");

public:
  using number_type                       = Number;
  using value_type                        = VectorizedArrayType;
  using gradient_type                     = Tensor<1, 1, VectorizedArrayType>;
  static constexpr unsigned int dimension = 1;
  using BaseClass =
    FEEvaluationBase<1, 1, Number, is_face, VectorizedArrayType>;

  /** @copydoc FEEvaluationBase<1,1,Number,is_face>::get_dof_value()
   */
  value_type
  get_dof_value(const unsigned int dof) const;

  /** @copydoc FEEvaluationBase<1,1,Number,is_face>::submit_dof_value()
   */
  void
  submit_dof_value(const value_type val_in, const unsigned int dof);

  /** @copydoc FEEvaluationBase<1,1,Number,is_face>::get_value()
   */
  value_type
  get_value(const unsigned int q_point) const;

  /** @copydoc FEEvaluationBase<1,1,Number,is_face>::submit_value()
   */
  void
  submit_value(const value_type val_in, const unsigned int q_point);

  /** @copydoc FEEvaluationBase<1,1,Number,is_face>::submit_value()
   */
  void
  submit_value(const gradient_type val_in, const unsigned int q_point);

  /** @copydoc FEEvaluationBase<1,1,Number,is_face>::get_gradient()
   */
  gradient_type
  get_gradient(const unsigned int q_point) const;

  /** @copydoc FEEvaluationBase<dim,1,Number,is_face>::get_normal_derivative()
   */
  value_type
  get_normal_derivative(const unsigned int q_point) const;

  /** @copydoc FEEvaluationBase<1,1,Number,is_face>::submit_gradient()
   */
  void
  submit_gradient(const gradient_type grad_in, const unsigned int q_point);

  /** @copydoc FEEvaluationBase<1,1,Number,is_face>::submit_gradient()
   */
  void
  submit_gradient(const value_type grad_in, const unsigned int q_point);

  /** @copydoc FEEvaluationBase<1,1,Number,is_face>::submit_normal_derivative()
   */
  void
  submit_normal_derivative(const value_type   grad_in,
                           const unsigned int q_point);

  /** @copydoc FEEvaluationBase<1,1,Number,is_face>::submit_normal_derivative()
   */
  void
  submit_normal_derivative(const gradient_type grad_in,
                           const unsigned int  q_point);

  /** @copydoc FEEvaluationBase<1,1,Number,is_face>::get_hessian()
   */
  Tensor<2, 1, VectorizedArrayType>
  get_hessian(unsigned int q_point) const;

  /** @copydoc FEEvaluationBase<1,1,Number,is_face>::get_hessian_diagonal()
   */
  gradient_type
  get_hessian_diagonal(const unsigned int q_point) const;

  /** @copydoc FEEvaluationBase<1,1,Number,is_face>::get_laplacian()
   */
  value_type
  get_laplacian(const unsigned int q_point) const;

  /** @copydoc FEEvaluationBase<1,1,Number,is_face>::integrate_value()
   */
  value_type
  integrate_value() const;

protected:
  /**
   * Constructor. Made protected to avoid initialization in user code. Takes
   * all data stored in MatrixFree. If applied to problems with more than one
   * finite element or more than one quadrature formula selected during
   * construction of @p matrix_free, @p first_selected_component and @p
   * quad_no allow to select the appropriate components.
   */
  FEEvaluationAccess(
    const MatrixFree<1, Number, VectorizedArrayType> &matrix_free,
    const unsigned int                                dof_no,
    const unsigned int                                first_selected_component,
    const unsigned int                                quad_no,
    const unsigned int                                fe_degree,
    const unsigned int                                n_q_points,
    const bool                                        is_interior_face = true);

  /**
   * Constructor with reduced functionality for similar usage of FEEvaluation
   * as FEValues, including matrix assembly.
   */
  template <int n_components_other>
  FEEvaluationAccess(const Mapping<1> &      mapping,
                     const FiniteElement<1> &fe,
                     const Quadrature<1> &   quadrature,
                     const UpdateFlags       update_flags,
                     const unsigned int      first_selected_component,
                     const FEEvaluationBase<1,
                                            n_components_other,
                                            Number,
                                            is_face,
                                            VectorizedArrayType> *other);

  /**
   * Copy constructor
   */
  FEEvaluationAccess(const FEEvaluationAccess &other);

  /**
   * Copy assignment operator
   */
  FEEvaluationAccess &
  operator=(const FEEvaluationAccess &other);
};



/**
 * The class that provides all functions necessary to evaluate functions at
 * quadrature points and cell integrations. In functionality, this class is
 * similar to FEValues, however, it includes a lot of specialized functions
 * that make it much faster (between 5 and 500, depending on the polynomial
 * degree). For evaluation of face terms in DG, see the class
 * FEFaceEvaluation.
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
 *         for (unsigned int v=0;
 *              v<VectorizedArray<double>::n_array_elements;
 *              ++v)
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
 * based on the AffineConstraints object specified at the MatrixFree::reinit()
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
 *           (dofs_per_cell - i) :
 *           VectorizedArray<double>::n_array_elements;
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
 * get_value -> Tensor<1,n_components,VectorizedArray<double>>
 * get_gradient -> Tensor<1,n_components,Tensor<1,dim,VectorizedArray<double>>
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
 * <h3>Handling several integration tasks and data storage in quadrature
 * points</h3>
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
 * different DoFHandler or AffineConstraints arguments. More precisely, even
 * though the layout is going to be the same in serial, there is no guarantee
 * about the ordering for different DoFHandler/AffineConstraints in the MPI
 * case. The reason is that the algorithm detects cells that need data
 * exchange with MPI and those can change for different elements &mdash; FE_Q
 * with hanging node constraints connects to more neighbors than a FE_DGQ
 * element, for instance, and cells which need data exchange are put in
 * different positions inside the cell loop. Of course, if the exact same
 * DoFHandler, AffineConstraints, and options (such as the setting for thread
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
 * @ingroup matrixfree
 *
 * @author Katharina Kormann and Martin Kronbichler, 2010, 2011
 */
template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
class FEEvaluation : public FEEvaluationAccess<dim,
                                               n_components_,
                                               Number,
                                               false,
                                               VectorizedArrayType>
{
  static_assert(
    std::is_same<Number, typename VectorizedArrayType::value_type>::value,
    "Type of Number and of VectorizedArrayType do not match.");

public:
  /**
   * An alias to the base class.
   */
  using BaseClass =
    FEEvaluationAccess<dim, n_components_, Number, false, VectorizedArrayType>;

  /**
   * An underlying number type specified as template argument.
   */
  using number_type = Number;

  /**
   * The type of function values, e.g. `VectorizedArrayType` for
   * `n_components=1` or `Tensor<1,dim,VectorizedArrayType >` for
   * `n_components=dim`.
   */
  using value_type = typename BaseClass::value_type;

  /**
   * The type of gradients, e.g. `Tensor<1,dim,VectorizedArrayType>` for
   * `n_components=1` or `Tensor<2,dim,VectorizedArrayType >` for
   * `n_components=dim`.
   */
  using gradient_type = typename BaseClass::gradient_type;

  /**
   * The dimension given as template argument.
   */
  static constexpr unsigned int dimension = dim;

  /**
   * The number of solution components of the evaluator given as template
   * argument.
   */
  static constexpr unsigned int n_components = n_components_;

  /**
   * The static number of quadrature points determined from the given template
   * argument `n_q_points_1d`. Note that the actual number of quadrature
   * points, `n_q_points`, can be different if `fe_degree=-1` is given and
   * run-time loop lengths are used rather than compile time ones.
   */
  static constexpr unsigned int static_n_q_points =
    Utilities::pow(n_q_points_1d, dim);

  /**
   * The static number of degrees of freedom of a scalar component determined
   * from the given template argument `fe_degree`. Note that the actual number
   * of degrees of freedom `dofs_per_component` can be different if
   * `fe_degree=-1` is given or if the underlying is of more complicated type
   * than the usual FE_Q or FE_DGQ ones, such as FE_DGP.
   */
  static constexpr unsigned int static_dofs_per_component =
    Utilities::pow(fe_degree + 1, dim);

  /**
   * The static number of degrees of freedom of all components determined from
   * the given template argument `fe_degree`. Note that the actual number of
   * degrees of freedom `dofs_per_cell` can be different if `fe_degree=-1` is
   * given or if the underlying is of more complicated type than the usual
   * FE_Q or FE_DGQ ones, such as FE_DGP.
   */
  static constexpr unsigned int tensor_dofs_per_cell =
    static_dofs_per_component * n_components;

  /**
   * The static number of degrees of freedom of all components determined from
   * the given template argument `fe_degree`. Note that the actual number of
   * degrees of freedom `dofs_per_cell` can be different if `fe_degree=-1` is
   * given or if the underlying is of more complicated type than the usual
   * FE_Q or FE_DGQ ones, such as FE_DGP.
   */
  static constexpr unsigned int static_dofs_per_cell =
    static_dofs_per_component * n_components;

  /**
   * Constructor. Takes all data stored in MatrixFree. If applied to problems
   * with more than one finite element or more than one quadrature formula
   * selected during construction of @p matrix_free, the appropriate component
   * can be selected by the optional arguments.
   *
   * @param matrix_free Data object that contains all data
   *
   * @param dof_no If matrix_free was set up with multiple DoFHandler
   * objects, this parameter selects to which DoFHandler/AffineConstraints pair
   * the given evaluator should be attached to.
   *
   * @param quad_no If matrix_free was set up with multiple Quadrature
   * objects, this parameter selects the appropriate number of the quadrature
   * formula.
   *
   * @param first_selected_component If the dof_handler selected by dof_no
   * uses an FESystem consisting of more than one component, this parameter
   * allows for selecting the component where the current evaluation routine
   * should start. Note that one evaluator does not support combining
   * different shape functions in different components. In other words, the
   * same base element of a FESystem needs to be set for the components
   * between @p first_selected_component and
   * <code>first_selected_component+n_components_</code>.
   */
  FEEvaluation(const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
               const unsigned int                                  dof_no  = 0,
               const unsigned int                                  quad_no = 0,
               const unsigned int first_selected_component                 = 0);

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
  FEEvaluation(const Mapping<dim> &      mapping,
               const FiniteElement<dim> &fe,
               const Quadrature<1> &     quadrature,
               const UpdateFlags         update_flags,
               const unsigned int        first_selected_component = 0);

  /**
   * Constructor for the reduced functionality. This constructor is equivalent
   * to the other one except that it makes the object use a $Q_1$ mapping
   * (i.e., an object of type MappingQGeneric(1)) implicitly.
   */
  FEEvaluation(const FiniteElement<dim> &fe,
               const Quadrature<1> &     quadrature,
               const UpdateFlags         update_flags,
               const unsigned int        first_selected_component = 0);

  /**
   * Constructor for the reduced functionality. Similar to the other
   * constructor with FiniteElement argument but using another
   * FEEvaluationBase object to provide information about the geometry. This
   * allows several FEEvaluation objects to share the geometry evaluation, i.e.,
   * the underlying mapping and quadrature points do only need to be evaluated
   * once. Make sure to not pass an optional object around when you intend to
   * use the FEEvaluation object in %parallel to the given one because
   * otherwise the intended sharing may create race conditions.
   */
  template <int n_components_other>
  FEEvaluation(const FiniteElement<dim> &                   fe,
               const FEEvaluationBase<dim,
                                      n_components_other,
                                      Number,
                                      false,
                                      VectorizedArrayType> &other,
               const unsigned int first_selected_component = 0);

  /**
   * Copy constructor. If FEEvaluationBase was constructed from a mapping, fe,
   * quadrature, and update flags, the underlying geometry evaluation based on
   * FEValues will be deep-copied in order to allow for using in parallel with
   * threads.
   */
  FEEvaluation(const FEEvaluation &other);

  /**
   * Copy assignment operator. If FEEvaluationBase was constructed from a
   * mapping, fe, quadrature, and update flags, the underlying geometry
   * evaluation based on FEValues will be deep-copied in order to allow for
   * using in parallel with threads.
   */
  FEEvaluation &
  operator=(const FEEvaluation &other);

  /**
   * Initialize the operation pointer to the current cell batch index. Unlike
   * the reinit functions taking a cell iterator as argument below and the
   * FEValues::reinit() methods, where the information related to a particular
   * cell is generated in the reinit call, this function is very cheap since
   * all data is pre-computed in @p matrix_free, and only a few indices have
   * to be set appropriately.
   */
  void
  reinit(const unsigned int cell_batch_index);

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
  void
  reinit(const TriaIterator<DoFCellAccessor<DoFHandlerType, level_dof_access>>
           &cell);

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
  void
  reinit(const typename Triangulation<dim>::cell_iterator &cell);

  /**
   * Evaluate the function values, the gradients, and the Hessians of the
   * polynomial interpolation from the DoF values in the input vector to the
   * quadrature points on the unit cell.  The function arguments specify which
   * parts shall actually be computed. This function has to be called first so
   * that the access functions @p get_value(), @p get_gradient() or @p
   * get_laplacian give useful information (unless these values have been set
   * manually).
   */
  void
  evaluate(const bool evaluate_values,
           const bool evaluate_gradients,
           const bool evaluate_hessians = false);

  /**
   * Evaluate the function values, the gradients, and the Hessians of the
   * polynomial interpolation from the DoF values in the input array @p
   * values_array to the quadrature points on the unit cell. If multiple
   * components are involved in the current FEEvaluation object, the sorting
   * in @p values_array is such that all degrees of freedom for the first
   * component come first, then all degrees of freedom for the second, and so
   * on. The function arguments specify which parts shall actually be
   * computed. This function has to be called first so that the access
   * functions @p get_value(), @p get_gradient() or @p get_laplacian give
   * useful information (unless these values have been set manually).
   */
  void
  evaluate(const VectorizedArrayType *values_array,
           const bool                 evaluate_values,
           const bool                 evaluate_gradients,
           const bool                 evaluate_hessians = false);

  /**
   * Read from the input vector and evaluates the function values, the
   * gradients, and the Hessians of the polynomial interpolation of the vector
   * entries from @p input_vector associated with the current cell to the
   * quadrature points on the unit cell. The function arguments specify which
   * parts shall actually be computed. This function has to be called first so
   * that the access functions @p get_value(), @p get_gradient() or @p
   * get_laplacian give useful information (unless these values have been set
   * manually).
   *
   * This call is equivalent to calling read_dof_values() followed by
   * evaluate(), but might internally use some additional optimizations.
   */
  template <typename VectorType>
  void
  gather_evaluate(const VectorType &input_vector,
                  const bool        evaluate_values,
                  const bool        evaluate_gradients,
                  const bool        evaluate_hessians = false);

  /**
   * This function takes the values and/or gradients that are stored on
   * quadrature points, tests them by all the basis functions/gradients on the
   * cell and performs the cell integration. The two function arguments
   * @p integrate_values and @p integrate_gradients are used to enable/disable
   * summation of the contributions submitted to the values or gradients slots,
   * respectively. The result is written into the internal data field
   * @p dof_values (that is usually written into the result vector by the
   * distribute_local_to_global() or set_dof_values() methods).
   */
  void
  integrate(const bool integrate_values, const bool integrate_gradients);

  /**
   * This function takes the values and/or gradients that are stored on
   * quadrature points, tests them by all the basis functions/gradients on the
   * cell and performs the cell integration. The two function arguments @p
   * integrate_values and @p integrate_gradients are used to enable/disable
   * summation of the contributions submitted to the values or gradients
   * slots, respectively. As opposed to the other integrate() method, this
   * call stores the result of the testing in the given array @p values_array,
   * whose previous results is overwritten, rather than writing it on the
   * internal data structures behind begin_dof_values().
   */
  void
  integrate(const bool           integrate_values,
            const bool           integrate_gradients,
            VectorizedArrayType *values_array);

  /**
   * This function takes the values and/or gradients that are stored on
   * quadrature points, tests them by all the basis functions/gradients on the
   * cell, performs the cell integration, and adds the result into the global
   * vector @p output_vector on the degrees of freedom associated with the
   * present cell index. The two function arguments @p integrate_values and
   * @p integrate_gradients are used to enable/disable summation of the
   * contributions submitted to the values or gradients slots, respectively.
   *
   * This call is equivalent to calling integrate() followed by
   * distribute_local_to_global(), but might internally use
   * some additional optimizations.
   */
  template <typename VectorType>
  void
  integrate_scatter(const bool  integrate_values,
                    const bool  integrate_gradients,
                    VectorType &output_vector);

  /**
   * Return the q-th quadrature point in real coordinates stored in
   * MappingInfo.
   */
  Point<dim, VectorizedArrayType>
  quadrature_point(const unsigned int q_point) const;

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
   * The number of quadrature points in use. If the number of quadrature
   * points in 1d is given as a template, this number is simply the
   * <tt>dim</tt>-th power of that value. If the element degree is set to -1
   * (dynamic selection of element degree), the static value of quadrature
   * points is inaccurate and this value must be used instead.
   */
  const unsigned int n_q_points;

private:
  /**
   * Checks if the template arguments regarding degree of the element
   * corresponds to the actual element used at initialization.
   */
  void
  check_template_arguments(const unsigned int fe_no,
                           const unsigned int first_selected_component);
};



/**
 * The class that provides all functions necessary to evaluate functions at
 * quadrature points and face integrations. The design of the class is similar
 * to FEEvaluation and most of the interfaces are shared with that class, in
 * particular most access functions that come from the common base classes
 * FEEvaluationAccess and FEEvaluationBase. Furthermore, the relation of this
 * class to FEEvaluation is similar to the relation between FEValues and
 * FEFaceValues.
 *
 * @tparam dim Dimension in which this class is to be used
 *
 * @tparam fe_degree Degree of the tensor product finite element with
 *                  fe_degree+1 degrees of freedom per coordinate
 *                  direction. If set to -1, the degree of the underlying
 *                  element will be used, which acts as a run time constant
 *                  rather than a compile time constant that slows down the
 *                  execution.
 *
 * @tparam n_q_points_1d Number of points in the quadrature formula in 1D,
 *                  usually chosen as fe_degree+1
 *
 * @tparam n_components Number of vector components when solving a system of
 *                  PDEs. If the same operation is applied to several
 *                  components of a PDE (e.g. a vector Laplace equation), they
 *                  can be applied simultaneously with one call (and often
 *                  more efficiently)
 *
 * @tparam Number Number format, usually @p double or @p float
 *
 * @tparam VectorizedArrayType Type of array to be woked on in a vectorized
 *                             fashion, defaults to VectorizedArray<Number>
 *
 * @note Currently only VectorizedArray<Number, width> is supported as
 *       VectorizedArrayType.
 *
 * @author Katharina Kormann and Martin Kronbichler, 2018
 */
template <int dim,
          int fe_degree,
          int n_q_points_1d            = fe_degree + 1,
          int n_components_            = 1,
          typename Number              = double,
          typename VectorizedArrayType = VectorizedArray<Number>>
class FEFaceEvaluation : public FEEvaluationAccess<dim,
                                                   n_components_,
                                                   Number,
                                                   true,
                                                   VectorizedArrayType>
{
  static_assert(
    std::is_same<Number, typename VectorizedArrayType::value_type>::value,
    "Type of Number and of VectorizedArrayType do not match.");

public:
  /**
   * An alias to the base class.
   */
  using BaseClass =
    FEEvaluationAccess<dim, n_components_, Number, true, VectorizedArrayType>;

  /**
   * A underlying number type specified as template argument.
   */
  using number_type = Number;

  /**
   * The type of function values, e.g. `VectorizedArrayType` for
   * `n_components=1` or `Tensor<1,dim,VectorizedArrayType >` for
   * `n_components=dim`.
   */
  using value_type = typename BaseClass::value_type;

  /**
   * The type of gradients, e.g. `Tensor<1,dim,VectorizedArrayType>` for
   * `n_components=1` or `Tensor<2,dim,VectorizedArrayType >` for
   * `n_components=dim`.
   */
  using gradient_type = typename BaseClass::gradient_type;

  /**
   * The dimension given as template argument.
   */
  static constexpr unsigned int dimension = dim;

  /**
   * The number of solution components of the evaluator given as template
   * argument.
   */
  static constexpr unsigned int n_components = n_components_;

  /**
   * The static number of quadrature points determined from the given template
   * argument `n_q_points_1d` taken to the power of dim-1. Note that the actual
   * number of quadrature points, `n_q_points`, can be different if
   * `fe_degree=-1` is given and run-time loop lengths are used rather than
   * compile time ones.
   */
  static constexpr unsigned int static_n_q_points =
    Utilities::pow(n_q_points_1d, dim - 1);

  /**
   * The static number of quadrature points on a cell with the same quadrature
   * formula. Note that this value is only present for simpler comparison with
   * the cell quadrature, as the actual number of points is given to a face by
   * the `static_n_q_points` variable.
   */
  static constexpr unsigned int static_n_q_points_cell =
    Utilities::pow(n_q_points_1d, dim);

  /**
   * The static number of degrees of freedom of a scalar component determined
   * from the given template argument `fe_degree`. Note that the actual number
   * of degrees of freedom `dofs_per_component` can be different if
   * `fe_degree=-1` is given.
   */
  static constexpr unsigned int static_dofs_per_component =
    Utilities::pow(fe_degree + 1, dim);

  /**
   * The static number of degrees of freedom of all components determined from
   * the given template argument `fe_degree`. Note that the actual number of
   * degrees of freedom `dofs_per_cell` can be different if `fe_degree=-1` is
   * given.
   */
  static constexpr unsigned int tensor_dofs_per_cell =
    static_dofs_per_component * n_components;

  /**
   * The static number of degrees of freedom of all components determined from
   * the given template argument `fe_degree`. Note that the actual number of
   * degrees of freedom `dofs_per_cell` can be different if `fe_degree=-1` is
   * given.
   */
  static constexpr unsigned int static_dofs_per_cell =
    static_dofs_per_component * n_components;

  /**
   * Constructor. Takes all data stored in MatrixFree. If applied to problems
   * with more than one finite element or more than one quadrature formula
   * selected during construction of @p matrix_free, the appropriate component
   * can be selected by the optional arguments.
   *
   * @param matrix_free Data object that contains all data
   *
   * @param is_interior_face This selects which of the two cells of an
   * internal face the current evaluator will be based upon. The interior face
   * is the main face along which the normal vectors are oriented. The
   * exterior face coming from the other side provides the same normal vector
   * as the interior side, so if the outer normal vector to that side is
   * desired, it must be multiplied by -1.
   *
   * @param dof_no If matrix_free was set up with multiple DoFHandler
   * objects, this parameter selects to which DoFHandler/AffineConstraints pair
   * the given evaluator should be attached to.
   *
   * @param quad_no If matrix_free was set up with multiple Quadrature
   * objects, this parameter selects the appropriate number of the quadrature
   * formula.
   *
   * @param first_selected_component If the dof_handler selected by dof_no
   * uses an FESystem consisting of more than one base element, this parameter
   * selects the number of the base element in FESystem. Note that this does
   * not directly relate to the component of the respective element due to the
   * possibility for a multiplicity in the element.
   */
  FEFaceEvaluation(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    const bool                                          is_interior_face = true,
    const unsigned int                                  dof_no           = 0,
    const unsigned int                                  quad_no          = 0,
    const unsigned int first_selected_component                          = 0);

  /**
   * Initializes the operation pointer to the current face. This method is the
   * default choice for face integration as the data stored in MappingInfo is
   * stored according to this numbering. Unlike the reinit functions taking a
   * cell iterator as argument below and the FEValues::reinit() methods, where
   * the information related to a particular cell is generated in the reinit
   * call, this function is very cheap since all data is pre-computed in
   * @p matrix_free, and only a few indices and pointers have to be set
   * appropriately.
   */
  void
  reinit(const unsigned int face_batch_number);

  /**
   * As opposed to the reinit() method from the base class, this reinit()
   * method initializes for a given number of cells and a face number. This
   * method is less efficient than the other reinit() method taking a
   * numbering of the faces because it needs to copy the data associated with
   * the faces to the cells in this call.
   */
  void
  reinit(const unsigned int cell_batch_number, const unsigned int face_number);

  /**
   * Evaluates the function values, the gradients, and the Laplacians of the
   * FE function given at the DoF values stored in the internal data field
   * `dof_values` (that is usually filled by the read_dof_values() method) at
   * the quadrature points on the unit cell.  The function arguments specify
   * which parts shall actually be computed. Needs to be called before the
   * functions get_value(), get_gradient() or get_normal_derivative() give
   * useful information (unless these values have been set manually by
   * accessing the internal data pointers).
   */
  void
  evaluate(const bool evaluate_values, const bool evaluate_gradients);

  /**
   * Evaluates the function values, the gradients, and the Laplacians of the
   * FE function given at the DoF values in the input array `values_array` at
   * the quadrature points on the unit cell. If multiple components are
   * involved in the current FEEvaluation object, the sorting in values_array
   * is such that all degrees of freedom for the first component come first,
   * then all degrees of freedom for the second, and so on. The function
   * arguments specify which parts shall actually be computed. Needs to be
   * called before the functions get_value(), get_gradient(), or
   * get_normal_derivative() give useful information (unless these values have
   * been set manually).
   */
  void
  evaluate(const VectorizedArrayType *values_array,
           const bool                 evaluate_values,
           const bool                 evaluate_gradients);

  /**
   * Reads from the input vector and evaluates the function values, the
   * gradients, and the Laplacians of the FE function at the quadrature points
   * on the unit cell. The function arguments specify which parts shall
   * actually be computed. Needs to be called before the functions
   * get_value(), get_gradient(), or get_normal_derivative() give useful
   * information.
   *
   * This call is equivalent to calling read_dof_values() followed by
   * evaluate(), but might internally use some additional optimizations.
   */
  template <typename VectorType>
  void
  gather_evaluate(const VectorType &input_vector,
                  const bool        evaluate_values,
                  const bool        evaluate_gradients);

  /**
   * This function takes the values and/or gradients that are stored on
   * quadrature points, tests them by all the basis functions/gradients on the
   * cell and performs the cell integration. The two function arguments
   * `integrate_val` and `integrate_grad` are used to enable/disable some of
   * values or gradients. The result is written into the internal data field
   * `dof_values` (that is usually written into the result vector by the
   * distribute_local_to_global() or set_dof_values() methods).
   */
  void
  integrate(const bool integrate_values, const bool integrate_gradients);

  /**
   * This function takes the values and/or gradients that are stored on
   * quadrature points, tests them by all the basis functions/gradients on the
   * cell and performs the cell integration. The two function arguments
   * `integrate_val` and `integrate_grad` are used to enable/disable some of
   * values or gradients. As opposed to the other integrate() method, this
   * call stores the result of the testing in the given array `values_array`.
   */
  void
  integrate(const bool           integrate_values,
            const bool           integrate_gradients,
            VectorizedArrayType *values_array);

  /**
   * This function takes the values and/or gradients that are stored on
   * quadrature points, tests them by all the basis functions/gradients on the
   * cell and performs the cell integration. The two function arguments
   * `integrate_val` and `integrate_grad` are used to enable/disable some of
   * values or gradients.
   *
   * This call is equivalent to calling integrate() followed by
   * distribute_local_to_global(), but might internally use some additional
   * optimizations.
   */
  template <typename VectorType>
  void
  integrate_scatter(const bool  integrate_values,
                    const bool  integrate_gradients,
                    VectorType &output_vector);

  /**
   * Returns the q-th quadrature point on the face in real coordinates stored
   * in MappingInfo.
   */
  Point<dim, VectorizedArrayType>
  quadrature_point(const unsigned int q_point) const;

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
   * The number of quadrature points in use. If the number of quadrature
   * points in 1d is given as a template, this number is simply the
   * <tt>dim-1</tt>-th power of that value. If the element degree is set to -1
   * (dynamic selection of element degree), the static value of quadrature
   * points is inaccurate and this value must be used instead.
   */
  const unsigned int n_q_points;
};



namespace internal
{
  namespace MatrixFreeFunctions
  {
    // a helper function to compute the number of DoFs of a DGP element at
    // compile time, depending on the degree
    template <int dim, int degree>
    struct DGP_dofs_per_component
    {
      // this division is always without remainder
      static constexpr unsigned int value =
        (DGP_dofs_per_component<dim - 1, degree>::value * (degree + dim)) / dim;
    };

    // base specialization: 1d elements have 'degree+1' degrees of freedom
    template <int degree>
    struct DGP_dofs_per_component<1, degree>
    {
      static constexpr unsigned int value = degree + 1;
    };
  } // namespace MatrixFreeFunctions
} // namespace internal


/*----------------------- Inline functions ----------------------------------*/

#ifndef DOXYGEN



/*----------------------- FEEvaluationBase ----------------------------------*/

template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline FEEvaluationBase<dim,
                        n_components_,
                        Number,
                        is_face,
                        VectorizedArrayType>::
  FEEvaluationBase(const MatrixFree<dim, Number, VectorizedArrayType> &data_in,
                   const unsigned int                                  dof_no,
                   const unsigned int first_selected_component,
                   const unsigned int quad_no_in,
                   const unsigned int fe_degree,
                   const unsigned int n_q_points,
                   const bool         is_interior_face)
  : scratch_data_array(data_in.acquire_scratch_data())
  , quad_no(quad_no_in)
  , n_fe_components(data_in.get_dof_info(dof_no).start_components.back())
  , active_fe_index(fe_degree != numbers::invalid_unsigned_int ?
                      data_in.get_dof_info(dof_no).fe_index_from_degree(
                        first_selected_component,
                        fe_degree) :
                      0)
  , active_quad_index(fe_degree != numbers::invalid_unsigned_int ?
                        (is_face ? data_in.get_mapping_info()
                                     .face_data[quad_no_in]
                                     .quad_index_from_n_q_points(n_q_points) :
                                   data_in.get_mapping_info()
                                     .cell_data[quad_no_in]
                                     .quad_index_from_n_q_points(n_q_points)) :
                        0)
  , n_quadrature_points(fe_degree != numbers::invalid_unsigned_int ?
                          n_q_points :
                          (is_face ? data_in
                                       .get_shape_info(dof_no,
                                                       quad_no_in,
                                                       active_fe_index,
                                                       active_quad_index)
                                       .n_q_points_face :
                                     data_in
                                       .get_shape_info(dof_no,
                                                       quad_no_in,
                                                       active_fe_index,
                                                       active_quad_index)
                                       .n_q_points))
  , matrix_info(&data_in)
  , dof_info(&data_in.get_dof_info(dof_no))
  , mapping_data(
      internal::MatrixFreeFunctions::
        MappingInfoCellsOrFaces<dim, Number, is_face, VectorizedArrayType>::get(
          data_in.get_mapping_info(),
          quad_no))
  , data(&data_in.get_shape_info(
      dof_no,
      quad_no_in,
      dof_info->component_to_base_index[first_selected_component],
      active_fe_index,
      active_quad_index))
  , jacobian(nullptr)
  , J_value(nullptr)
  , normal_vectors(nullptr)
  , normal_x_jacobian(nullptr)
  , quadrature_weights(
      mapping_data->descriptor[active_quad_index].quadrature_weights.begin())
  , cell(numbers::invalid_unsigned_int)
  , is_interior_face(is_interior_face)
  , dof_access_index(
      is_face ?
        (is_interior_face ?
           internal::MatrixFreeFunctions::DoFInfo::dof_access_face_interior :
           internal::MatrixFreeFunctions::DoFInfo::dof_access_face_exterior) :
        internal::MatrixFreeFunctions::DoFInfo::dof_access_cell)
  , cell_type(internal::MatrixFreeFunctions::general)
  , dof_values_initialized(false)
  , values_quad_initialized(false)
  , gradients_quad_initialized(false)
  , hessians_quad_initialized(false)
  , values_quad_submitted(false)
  , gradients_quad_submitted(false)
  , first_selected_component(first_selected_component)
{
  set_data_pointers();
  Assert(matrix_info->mapping_initialized() == true, ExcNotInitialized());
  AssertDimension(matrix_info->get_size_info().vectorization_length,
                  VectorizedArrayType::n_array_elements);
  AssertDimension((is_face ? data->n_q_points_face : data->n_q_points),
                  n_quadrature_points);
  AssertDimension(n_quadrature_points,
                  mapping_data->descriptor[active_quad_index].n_q_points);
  Assert(
    dof_info->start_components.back() == 1 ||
      static_cast<int>(n_components_) <=
        static_cast<int>(
          dof_info->start_components
            [dof_info->component_to_base_index[first_selected_component] + 1]) -
          first_selected_component,
    ExcMessage(
      "You tried to construct a vector-valued evaluator with " +
      Utilities::to_string(n_components) +
      " components. However, "
      "the current base element has only " +
      Utilities::to_string(
        dof_info->start_components
          [dof_info->component_to_base_index[first_selected_component] + 1] -
        first_selected_component) +
      " components left when starting from local element index " +
      Utilities::to_string(
        first_selected_component -
        dof_info->start_components
          [dof_info->component_to_base_index[first_selected_component]]) +
      " (global index " + Utilities::to_string(first_selected_component) +
      ")"));

  // do not check for correct dimensions of data fields here, should be done
  // in derived classes
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
template <int n_components_other>
inline FEEvaluationBase<dim,
                        n_components_,
                        Number,
                        is_face,
                        VectorizedArrayType>::
  FEEvaluationBase(const Mapping<dim> &      mapping,
                   const FiniteElement<dim> &fe,
                   const Quadrature<1> &     quadrature,
                   const UpdateFlags         update_flags,
                   const unsigned int        first_selected_component,
                   const FEEvaluationBase<dim,
                                          n_components_other,
                                          Number,
                                          is_face,
                                          VectorizedArrayType> *other)
  : scratch_data_array(new AlignedVector<VectorizedArrayType>())
  , quad_no(numbers::invalid_unsigned_int)
  , n_fe_components(n_components_)
  , active_fe_index(numbers::invalid_unsigned_int)
  , active_quad_index(numbers::invalid_unsigned_int)
  , n_quadrature_points(
      Utilities::fixed_power < is_face ? dim - 1 : dim > (quadrature.size()))
  , matrix_info(nullptr)
  , dof_info(nullptr)
  , mapping_data(nullptr)
  ,
  // select the correct base element from the given FE component
  data(new internal::MatrixFreeFunctions::ShapeInfo<VectorizedArrayType>(
    quadrature,
    fe,
    fe.component_to_base_index(first_selected_component).first))
  , jacobian(nullptr)
  , J_value(nullptr)
  , normal_vectors(nullptr)
  , normal_x_jacobian(nullptr)
  , quadrature_weights(nullptr)
  , cell(0)
  , cell_type(internal::MatrixFreeFunctions::general)
  , is_interior_face(true)
  , dof_access_index(internal::MatrixFreeFunctions::DoFInfo::dof_access_cell)
  , dof_values_initialized(false)
  , values_quad_initialized(false)
  , gradients_quad_initialized(false)
  , hessians_quad_initialized(false)
  , values_quad_submitted(false)
  , gradients_quad_submitted(false)
  ,
  // keep the number of the selected component within the current base element
  // for reading dof values
  first_selected_component(first_selected_component)
{
  set_data_pointers();

  Assert(other == nullptr || other->mapped_geometry.get() != nullptr,
         ExcInternalError());
  if (other != nullptr &&
      other->mapped_geometry->get_quadrature() == quadrature)
    mapped_geometry = other->mapped_geometry;
  else
    mapped_geometry =
      std::make_shared<internal::MatrixFreeFunctions::
                         MappingDataOnTheFly<dim, Number, VectorizedArrayType>>(
        mapping, quadrature, update_flags);
  cell = 0;

  mapping_data = &mapped_geometry->get_data_storage();
  jacobian     = mapped_geometry->get_data_storage().jacobians[0].begin();
  J_value      = mapped_geometry->get_data_storage().JxW_values.begin();

  const unsigned int base_element_number =
    fe.component_to_base_index(first_selected_component).first;
  Assert(fe.element_multiplicity(base_element_number) == 1 ||
           fe.element_multiplicity(base_element_number) -
               first_selected_component >=
             n_components_,
         ExcMessage("The underlying element must at least contain as many "
                    "components as requested by this class"));
  (void)base_element_number;
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline FEEvaluationBase<dim,
                        n_components_,
                        Number,
                        is_face,
                        VectorizedArrayType>::
  FEEvaluationBase(const FEEvaluationBase<dim,
                                          n_components_,
                                          Number,
                                          is_face,
                                          VectorizedArrayType> &other)
  : scratch_data_array(other.matrix_info == nullptr ?
                         new AlignedVector<VectorizedArrayType>() :
                         other.matrix_info->acquire_scratch_data())
  , quad_no(other.quad_no)
  , n_fe_components(other.n_fe_components)
  , active_fe_index(other.active_fe_index)
  , active_quad_index(other.active_quad_index)
  , n_quadrature_points(other.n_quadrature_points)
  , matrix_info(other.matrix_info)
  , dof_info(other.dof_info)
  , mapping_data(other.mapping_data)
  , data(other.matrix_info == nullptr ?
           new internal::MatrixFreeFunctions::ShapeInfo<VectorizedArrayType>(
             *other.data) :
           other.data)
  , jacobian(nullptr)
  , J_value(nullptr)
  , normal_vectors(nullptr)
  , normal_x_jacobian(nullptr)
  , quadrature_weights(
      other.matrix_info == nullptr ?
        nullptr :
        mapping_data->descriptor[active_quad_index].quadrature_weights.begin())
  , cell(numbers::invalid_unsigned_int)
  , cell_type(internal::MatrixFreeFunctions::general)
  , is_interior_face(other.is_interior_face)
  , dof_access_index(other.dof_access_index)
  , dof_values_initialized(false)
  , values_quad_initialized(false)
  , gradients_quad_initialized(false)
  , hessians_quad_initialized(false)
  , values_quad_submitted(false)
  , gradients_quad_submitted(false)
  , first_selected_component(other.first_selected_component)
{
  set_data_pointers();

  // Create deep copy of mapped geometry for use in parallel...
  if (other.mapped_geometry.get() != nullptr)
    {
      mapped_geometry.reset(
        new internal::MatrixFreeFunctions::
          MappingDataOnTheFly<dim, Number, VectorizedArrayType>(
            other.mapped_geometry->get_fe_values().get_mapping(),
            other.mapped_geometry->get_quadrature(),
            other.mapped_geometry->get_fe_values().get_update_flags()));
      mapping_data = &mapped_geometry->get_data_storage();
      cell         = 0;

      jacobian = mapped_geometry->get_data_storage().jacobians[0].begin();
      J_value  = mapped_geometry->get_data_storage().JxW_values.begin();
    }
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline FEEvaluationBase<dim,
                        n_components_,
                        Number,
                        is_face,
                        VectorizedArrayType> &
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
operator=(const FEEvaluationBase<dim,
                                 n_components_,
                                 Number,
                                 is_face,
                                 VectorizedArrayType> &other)
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

  matrix_info  = other.matrix_info;
  dof_info     = other.dof_info;
  mapping_data = other.mapping_data;
  if (other.matrix_info == nullptr)
    {
      data = new internal::MatrixFreeFunctions::ShapeInfo<VectorizedArrayType>(
        *other.data);
      scratch_data_array = new AlignedVector<VectorizedArrayType>();
    }
  else
    {
      data               = other.data;
      scratch_data_array = matrix_info->acquire_scratch_data();
    }
  set_data_pointers();

  quadrature_weights =
    (mapping_data != nullptr ?
       mapping_data->descriptor[active_quad_index].quadrature_weights.begin() :
       nullptr);
  cell             = numbers::invalid_unsigned_int;
  cell_type        = internal::MatrixFreeFunctions::general;
  is_interior_face = other.is_interior_face;
  dof_access_index = other.dof_access_index;

  // Create deep copy of mapped geometry for use in parallel...
  if (other.mapped_geometry.get() != nullptr)
    {
      mapped_geometry.reset(
        new internal::MatrixFreeFunctions::
          MappingDataOnTheFly<dim, Number, VectorizedArrayType>(
            other.mapped_geometry->get_fe_values().get_mapping(),
            other.mapped_geometry->get_quadrature(),
            other.mapped_geometry->get_fe_values().get_update_flags()));
      cell         = 0;
      mapping_data = &mapped_geometry->get_data_storage();
      jacobian     = mapped_geometry->get_data_storage().jacobians[0].begin();
      J_value      = mapped_geometry->get_data_storage().JxW_values.begin();
    }

  return *this;
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline FEEvaluationBase<dim,
                        n_components_,
                        Number,
                        is_face,
                        VectorizedArrayType>::~FEEvaluationBase()
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
      data = nullptr;
    }
  scratch_data_array = nullptr;
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline void
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  set_data_pointers()
{
  Assert(scratch_data_array != nullptr, ExcInternalError());

  const unsigned int tensor_dofs_per_component =
    Utilities::fixed_power<dim>(this->data->fe_degree + 1);
  const unsigned int dofs_per_component =
    this->data->dofs_per_component_on_cell;
  const unsigned int n_quadrature_points =
    is_face ? this->data->n_q_points_face : this->data->n_q_points;

  const unsigned int shift =
    std::max(tensor_dofs_per_component + 1, dofs_per_component) *
      n_components_ * 3 +
    2 * n_quadrature_points;
  const unsigned int allocated_size =
    shift + n_components_ * dofs_per_component +
    (n_components_ * (dim * dim + 2 * dim + 1) * n_quadrature_points);
  scratch_data_array->resize_fast(allocated_size);

  // set the pointers to the correct position in the data array
  for (unsigned int c = 0; c < n_components_; ++c)
    {
      this->values_dofs[c] =
        scratch_data_array->begin() + c * dofs_per_component;
      this->values_quad[c] = scratch_data_array->begin() +
                             n_components * dofs_per_component +
                             c * n_quadrature_points;
      for (unsigned int d = 0; d < dim; ++d)
        this->gradients_quad[c][d] =
          scratch_data_array->begin() +
          n_components * (dofs_per_component + n_quadrature_points) +
          (c * dim + d) * n_quadrature_points;
      for (unsigned int d = 0; d < (dim * dim + dim) / 2; ++d)
        this->hessians_quad[c][d] =
          scratch_data_array->begin() +
          n_components *
            ((dim + 1) * n_quadrature_points + dofs_per_component) +
          (c * (dim * dim + dim) + d) * n_quadrature_points;
    }
  scratch_data =
    scratch_data_array->begin() + n_components_ * dofs_per_component +
    (n_components_ * (dim * dim + 2 * dim + 1) * n_quadrature_points);
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline unsigned int
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  get_cell_data_number() const
{
  return get_mapping_data_index_offset();
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline unsigned int
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  get_mapping_data_index_offset() const
{
  if (matrix_info == nullptr)
    return 0;
  else
    {
      AssertIndexRange(cell, this->mapping_data->data_index_offsets.size());
      return this->mapping_data->data_index_offsets[cell];
    }
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline internal::MatrixFreeFunctions::GeometryType
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  get_cell_type() const
{
  Assert(cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  return cell_type;
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline const internal::MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> &
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  get_shape_info() const
{
  Assert(data != nullptr, ExcInternalError());
  return *data;
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline void
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  fill_JxW_values(AlignedVector<VectorizedArrayType> &JxW_values) const
{
  AssertDimension(JxW_values.size(), n_quadrature_points);
  Assert(J_value != nullptr, ExcNotInitialized());
  if (this->cell_type <= internal::MatrixFreeFunctions::affine)
    {
      VectorizedArrayType J = J_value[0];
      for (unsigned int q = 0; q < this->n_quadrature_points; ++q)
        JxW_values[q] = J * this->quadrature_weights[q];
    }
  else
    for (unsigned int q = 0; q < n_quadrature_points; ++q)
      JxW_values[q] = J_value[q];
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE Tensor<1, dim, VectorizedArrayType>
                             FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  get_normal_vector(const unsigned int q_index) const
{
  AssertIndexRange(q_index, n_quadrature_points);
  Assert(normal_vectors != nullptr, ExcMessage("Did not call reinit()!"));
  if (this->cell_type <= internal::MatrixFreeFunctions::flat_faces)
    return normal_vectors[0];
  else
    return normal_vectors[q_index];
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE VectorizedArrayType
                             FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::JxW(
  const unsigned int q_index) const
{
  AssertIndexRange(q_index, n_quadrature_points);
  Assert(J_value != nullptr, ExcNotInitialized());
  if (this->cell_type <= internal::MatrixFreeFunctions::affine)
    {
      Assert(this->quadrature_weights != nullptr, ExcInternalError());
      return J_value[0] * this->quadrature_weights[q_index];
    }
  else
    return J_value[q_index];
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline Tensor<2, dim, VectorizedArrayType>
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  inverse_jacobian(const unsigned int q_index) const
{
  AssertIndexRange(q_index, n_quadrature_points);
  Assert(this->jacobian != nullptr, ExcNotImplemented());
  if (this->cell_type <= internal::MatrixFreeFunctions::affine)
    return jacobian[0];
  else
    return jacobian[q_index];
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline VectorizedArrayType
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  read_cell_data(const AlignedVector<VectorizedArrayType> &array) const
{
  Assert(matrix_info != nullptr, ExcNotImplemented());
  AssertDimension(array.size(),
                  matrix_info->get_task_info().cell_partition_data.back());
  if (is_face)
    {
      VectorizedArrayType out = make_vectorized_array<Number>(Number(1.));
      const unsigned int *cells =
        is_interior_face ?
          &this->matrix_info->get_face_info(cell).cells_interior[0] :
          &this->matrix_info->get_face_info(cell).cells_exterior[0];
      for (unsigned int i = 0; i < VectorizedArrayType::n_array_elements; ++i)
        if (cells[i] != numbers::invalid_unsigned_int)
          out[i] = array[cells[i] / VectorizedArrayType::n_array_elements]
                        [cells[i] % VectorizedArrayType::n_array_elements];
      return out;
    }
  else
    return array[cell];
}



namespace internal
{
  // allows to select between block vectors and non-block vectors, which
  // allows to use a unified interface for extracting blocks on block vectors
  // and doing nothing on usual vectors
  template <typename VectorType, bool>
  struct BlockVectorSelector
  {};

  template <typename VectorType>
  struct BlockVectorSelector<VectorType, true>
  {
    using BaseVectorType = typename VectorType::BlockType;

    static BaseVectorType *
    get_vector_component(VectorType &vec, const unsigned int component)
    {
      AssertIndexRange(component, vec.n_blocks());
      return &vec.block(component);
    }
  };

  template <typename VectorType>
  struct BlockVectorSelector<VectorType, false>
  {
    using BaseVectorType = VectorType;

    static BaseVectorType *
    get_vector_component(VectorType &vec, const unsigned int component)
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
  struct BlockVectorSelector<std::vector<VectorType>, false>
  {
    using BaseVectorType = VectorType;

    static BaseVectorType *
    get_vector_component(std::vector<VectorType> &vec,
                         const unsigned int       component)
    {
      AssertIndexRange(component, vec.size());
      return &vec[component];
    }
  };

  template <typename VectorType>
  struct BlockVectorSelector<std::vector<VectorType *>, false>
  {
    using BaseVectorType = VectorType;

    static BaseVectorType *
    get_vector_component(std::vector<VectorType *> &vec,
                         const unsigned int         component)
    {
      AssertIndexRange(component, vec.size());
      return vec[component];
    }
  };
} // namespace internal



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
template <typename VectorType, typename VectorOperation>
inline void
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  read_write_operation(
    const VectorOperation &                                   operation,
    VectorType *                                              src[],
    const std::bitset<VectorizedArrayType::n_array_elements> &mask,
    const bool apply_constraints) const
{
  // Case 1: No MatrixFree object given, simple case because we do not need to
  // process constraints and need not care about vectorization -> go to
  // separate function
  if (matrix_info == nullptr)
    {
      read_write_operation_global(operation, src);
      return;
    }

  Assert(dof_info != nullptr, ExcNotInitialized());
  Assert(matrix_info->indices_initialized() == true, ExcNotInitialized());
  if (n_fe_components == 1)
    for (unsigned int comp = 0; comp < n_components; ++comp)
      internal::check_vector_compatibility(*src[comp], *dof_info);
  else
    {
      internal::check_vector_compatibility(*src[0], *dof_info);
    }

  // Case 2: contiguous indices which use reduced storage of indices and can
  // use vectorized load/store operations -> go to separate function
  AssertIndexRange(cell,
                   dof_info->index_storage_variants[dof_access_index].size());
  if (dof_info->index_storage_variants
        [is_face ? dof_access_index :
                   internal::MatrixFreeFunctions::DoFInfo::dof_access_cell]
        [cell] >=
      internal::MatrixFreeFunctions::DoFInfo::IndexStorageVariants::contiguous)
    {
      read_write_operation_contiguous(operation, src, mask);
      return;
    }

  // Case 3: standard operation with one index per degree of freedom -> go on
  // here
  constexpr unsigned int n_vectorization =
    VectorizedArrayType::n_array_elements;
  Assert(mask.count() == n_vectorization,
         ExcNotImplemented("Masking currently not implemented for "
                           "non-contiguous DoF storage"));

  std::integral_constant<bool,
                         internal::is_vectorizable<VectorType, Number>::value>
    vector_selector;

  const unsigned int dofs_per_component =
    this->data->dofs_per_component_on_cell;
  if (dof_info->index_storage_variants
        [is_face ? dof_access_index :
                   internal::MatrixFreeFunctions::DoFInfo::dof_access_cell]
        [cell] ==
      internal::MatrixFreeFunctions::DoFInfo::IndexStorageVariants::interleaved)
    {
      const unsigned int *dof_indices =
        dof_info->dof_indices_interleaved.data() +
        dof_info->row_starts[cell * n_fe_components * n_vectorization].first +
        dof_info->component_dof_indices_offset[active_fe_index]
                                              [first_selected_component] *
          n_vectorization;
      if (n_components == 1 || n_fe_components == 1)
        for (unsigned int i = 0; i < dofs_per_component;
             ++i, dof_indices += n_vectorization)
          for (unsigned int comp = 0; comp < n_components; ++comp)
            operation.process_dof_gather(dof_indices,
                                         *src[comp],
                                         0,
                                         values_dofs[comp][i],
                                         vector_selector);
      else
        for (unsigned int comp = 0; comp < n_components; ++comp)
          for (unsigned int i = 0; i < dofs_per_component;
               ++i, dof_indices += n_vectorization)
            operation.process_dof_gather(
              dof_indices, *src[0], 0, values_dofs[comp][i], vector_selector);
      return;
    }

  const unsigned int *  dof_indices[n_vectorization];
  VectorizedArrayType **values_dofs =
    const_cast<VectorizedArrayType **>(&this->values_dofs[0]);

  unsigned int        cells_copied[n_vectorization];
  const unsigned int *cells;
  unsigned int        n_vectorization_actual =
    dof_info->n_vectorization_lanes_filled[dof_access_index][cell];
  bool has_constraints = false;
  if (is_face)
    {
      if (dof_access_index ==
          internal::MatrixFreeFunctions::DoFInfo::dof_access_cell)
        for (unsigned int v = 0; v < n_vectorization_actual; ++v)
          cells_copied[v] = cell * VectorizedArrayType::n_array_elements + v;
      cells = dof_access_index ==
                  internal::MatrixFreeFunctions::DoFInfo::dof_access_cell ?
                &cells_copied[0] :
                (is_interior_face ?
                   &this->matrix_info->get_face_info(cell).cells_interior[0] :
                   &this->matrix_info->get_face_info(cell).cells_exterior[0]);
      for (unsigned int v = 0; v < n_vectorization_actual; ++v)
        {
          Assert(cells[v] < dof_info->row_starts.size() - 1,
                 ExcInternalError());
          has_constraints =
            has_constraints &&
            dof_info
                ->row_starts[cells[v] * n_fe_components +
                             first_selected_component + n_components]
                .second != dof_info
                             ->row_starts[cells[v] * n_fe_components +
                                          first_selected_component]
                             .second;
          dof_indices[v] = dof_info->dof_indices.data() +
                           dof_info
                             ->row_starts[cells[v] * n_fe_components +
                                          first_selected_component]
                             .first;
        }
      for (unsigned int v = n_vectorization_actual; v < n_vectorization; ++v)
        dof_indices[v] = nullptr;
    }
  else
    {
      AssertIndexRange((cell + 1) * n_vectorization * n_fe_components,
                       dof_info->row_starts.size());
      const unsigned int n_components_read =
        n_fe_components > 1 ? n_components : 1;
      for (unsigned int v = 0; v < n_vectorization_actual; ++v)
        {
          if (dof_info
                ->row_starts[(cell * n_vectorization + v) * n_fe_components +
                             first_selected_component + n_components_read]
                .second !=
              dof_info
                ->row_starts[(cell * n_vectorization + v) * n_fe_components +
                             first_selected_component]
                .second)
            has_constraints = true;
          Assert(
            dof_info
                  ->row_starts[(cell * n_vectorization + v) * n_fe_components +
                               first_selected_component + n_components_read]
                  .first ==
                dof_info
                  ->row_starts[(cell * n_vectorization + v) * n_fe_components +
                               first_selected_component]
                  .first ||
              dof_info
                  ->row_starts[(cell * n_vectorization + v) * n_fe_components +
                               first_selected_component]
                  .first < dof_info->dof_indices.size(),
            ExcIndexRange(
              0,
              dof_info
                ->row_starts[(cell * n_vectorization + v) * n_fe_components +
                             first_selected_component]
                .first,
              dof_info->dof_indices.size()));
          dof_indices[v] =
            dof_info->dof_indices.data() +
            dof_info
              ->row_starts[(cell * n_vectorization + v) * n_fe_components +
                           first_selected_component]
              .first;
        }
      for (unsigned int v = n_vectorization_actual; v < n_vectorization; ++v)
        dof_indices[v] = nullptr;
    }

  // Case where we have no constraints throughout the whole cell: Can go
  // through the list of DoFs directly
  if (!has_constraints)
    {
      if (n_vectorization_actual < n_vectorization)
        for (unsigned int comp = 0; comp < n_components; ++comp)
          for (unsigned int i = 0; i < dofs_per_component; ++i)
            operation.process_empty(values_dofs[comp][i]);
      if (n_components == 1 || n_fe_components == 1)
        {
          for (unsigned int v = 0; v < n_vectorization_actual; ++v)
            for (unsigned int i = 0; i < dofs_per_component; ++i)
              for (unsigned int comp = 0; comp < n_components; ++comp)
                operation.process_dof(dof_indices[v][i],
                                      *src[comp],
                                      values_dofs[comp][i][v]);
        }
      else
        {
          for (unsigned int comp = 0; comp < n_components; ++comp)
            for (unsigned int v = 0; v < n_vectorization_actual; ++v)
              for (unsigned int i = 0; i < dofs_per_component; ++i)
                operation.process_dof(
                  dof_indices[v][comp * dofs_per_component + i],
                  *src[0],
                  values_dofs[comp][i][v]);
        }
      return;
    }

  // In the case where there are some constraints to be resolved, loop over
  // all vector components that are filled and then over local dofs. ind_local
  // holds local number on cell, index iterates over the elements of
  // index_local_to_global and dof_indices points to the global indices stored
  // in index_local_to_global
  if (n_vectorization_actual < n_vectorization)
    for (unsigned int comp = 0; comp < n_components; ++comp)
      for (unsigned int i = 0; i < dofs_per_component; ++i)
        operation.process_empty(values_dofs[comp][i]);
  for (unsigned int v = 0; v < n_vectorization_actual; ++v)
    {
      unsigned int       index_indicators, next_index_indicators;
      const unsigned int n_components_read =
        n_fe_components > 1 ? n_components : 1;
      if (is_face)
        {
          index_indicators = dof_info
                               ->row_starts[cells[v] * n_fe_components +
                                            first_selected_component]
                               .second;
          next_index_indicators = dof_info
                                    ->row_starts[cells[v] * n_fe_components +
                                                 first_selected_component + 1]
                                    .second;
        }
      else
        {
          index_indicators =
            dof_info
              ->row_starts[(cell * n_vectorization + v) * n_fe_components +
                           first_selected_component]
              .second;
          next_index_indicators =
            dof_info
              ->row_starts[(cell * n_vectorization + v) * n_fe_components +
                           first_selected_component + 1]
              .second;
        }

      if (apply_constraints == false &&
          dof_info
              ->row_starts[(cell * n_vectorization + v) * n_fe_components +
                           first_selected_component]
              .second !=
            dof_info
              ->row_starts[(cell * n_vectorization + v) * n_fe_components +
                           first_selected_component + n_components_read]
              .second)
        {
          Assert(
            dof_info->row_starts_plain_indices[cell * n_vectorization + v] !=
              numbers::invalid_unsigned_int,
            ExcNotInitialized());
          dof_indices[v] =
            dof_info->plain_dof_indices.data() +
            dof_info->component_dof_indices_offset[active_fe_index]
                                                  [first_selected_component] +
            (is_face ?
               dof_info->row_starts_plain_indices[cells[v]] :
               dof_info->row_starts_plain_indices[cell * n_vectorization + v]);
          next_index_indicators = index_indicators;
        }

      if (n_components == 1 || n_fe_components == 1)
        {
          for (unsigned int c = 0; c < n_components; ++c)
            Assert(src[c] != nullptr,
                   ExcMessage(
                     "The finite element underlying this FEEvaluation "
                     "object is scalar, but you requested " +
                     std::to_string(n_components) +
                     " components via the template argument in "
                     "FEEvaluation. In that case, you must pass an "
                     "std::vector<VectorType> or a BlockVector to " +
                     "read_dof_values and distribute_local_to_global."));

          unsigned int ind_local = 0;
          for (; index_indicators != next_index_indicators; ++index_indicators)
            {
              const std::pair<unsigned short, unsigned short> indicator =
                dof_info->constraint_indicator[index_indicators];
              // run through values up to next constraint
              for (unsigned int j = 0; j < indicator.first; ++j)
                for (unsigned int comp = 0; comp < n_components; ++comp)
                  operation.process_dof(dof_indices[v][j],
                                        *src[comp],
                                        values_dofs[comp][ind_local + j][v]);

              ind_local += indicator.first;
              dof_indices[v] += indicator.first;

              // constrained case: build the local value as a linear
              // combination of the global value according to constraints
              Number value[n_components];
              for (unsigned int comp = 0; comp < n_components; ++comp)
                operation.pre_constraints(values_dofs[comp][ind_local][v],
                                          value[comp]);

              const Number *data_val =
                matrix_info->constraint_pool_begin(indicator.second);
              const Number *end_pool =
                matrix_info->constraint_pool_end(indicator.second);
              for (; data_val != end_pool; ++data_val, ++dof_indices[v])
                for (unsigned int comp = 0; comp < n_components; ++comp)
                  operation.process_constraint(*dof_indices[v],
                                               *data_val,
                                               *src[comp],
                                               value[comp]);

              for (unsigned int comp = 0; comp < n_components; ++comp)
                operation.post_constraints(value[comp],
                                           values_dofs[comp][ind_local][v]);
              ind_local++;
            }

          AssertIndexRange(ind_local, dofs_per_component + 1);

          for (; ind_local < dofs_per_component; ++dof_indices[v], ++ind_local)
            for (unsigned int comp = 0; comp < n_components; ++comp)
              operation.process_dof(*dof_indices[v],
                                    *src[comp],
                                    values_dofs[comp][ind_local][v]);
        }
      else
        {
          // case with vector-valued finite elements where all components are
          // included in one single vector. Assumption: first come all entries
          // to the first component, then all entries to the second one, and
          // so on. This is ensured by the way MatrixFree reads out the
          // indices.
          for (unsigned int comp = 0; comp < n_components; ++comp)
            {
              unsigned int ind_local = 0;

              // check whether there is any constraint on the current cell
              for (; index_indicators != next_index_indicators;
                   ++index_indicators)
                {
                  const std::pair<unsigned short, unsigned short> indicator =
                    dof_info->constraint_indicator[index_indicators];

                  // run through values up to next constraint
                  for (unsigned int j = 0; j < indicator.first; ++j)
                    operation.process_dof(dof_indices[v][j],
                                          *src[0],
                                          values_dofs[comp][ind_local + j][v]);
                  ind_local += indicator.first;
                  dof_indices[v] += indicator.first;

                  // constrained case: build the local value as a linear
                  // combination of the global value according to constraints
                  Number value;
                  operation.pre_constraints(values_dofs[comp][ind_local][v],
                                            value);

                  const Number *data_val =
                    matrix_info->constraint_pool_begin(indicator.second);
                  const Number *end_pool =
                    matrix_info->constraint_pool_end(indicator.second);

                  for (; data_val != end_pool; ++data_val, ++dof_indices[v])
                    operation.process_constraint(*dof_indices[v],
                                                 *data_val,
                                                 *src[0],
                                                 value);

                  operation.post_constraints(value,
                                             values_dofs[comp][ind_local][v]);
                  ind_local++;
                }

              AssertIndexRange(ind_local, dofs_per_component + 1);

              // get the dof values past the last constraint
              for (; ind_local < dofs_per_component;
                   ++dof_indices[v], ++ind_local)
                {
                  AssertIndexRange(*dof_indices[v], src[0]->size());
                  operation.process_dof(*dof_indices[v],
                                        *src[0],
                                        values_dofs[comp][ind_local][v]);
                }

              if (apply_constraints == true && comp + 1 < n_components)
                {
                  if (is_face)
                    next_index_indicators =
                      dof_info
                        ->row_starts[cells[v] * n_fe_components +
                                     first_selected_component + comp + 2]
                        .second;
                  else
                    next_index_indicators =
                      dof_info
                        ->row_starts[(cell * n_vectorization + v) *
                                       n_fe_components +
                                     first_selected_component + comp + 2]
                        .second;
                }
            }
        }
    }
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
template <typename VectorType, typename VectorOperation>
inline void
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  read_write_operation_global(const VectorOperation &operation,
                              VectorType *           src[]) const
{
  Assert(!local_dof_indices.empty(), ExcNotInitialized());

  unsigned int index =
    first_selected_component * data->dofs_per_component_on_cell;
  for (unsigned int comp = 0; comp < n_components; ++comp)
    {
      for (unsigned int i = 0; i < data->dofs_per_component_on_cell;
           ++i, ++index)
        {
          operation.process_empty(values_dofs[comp][i]);
          operation.process_dof_global(
            local_dof_indices[data->lexicographic_numbering[index]],
            *src[0],
            values_dofs[comp][i][0]);
        }
    }
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
template <typename VectorType, typename VectorOperation>
inline void
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  read_write_operation_contiguous(
    const VectorOperation &                                   operation,
    VectorType *                                              src[],
    const std::bitset<VectorizedArrayType::n_array_elements> &mask) const
{
  // This functions processes the functions read_dof_values,
  // distribute_local_to_global, and set_dof_values with the same code for
  // contiguous cell indices (DG case). The distinction between these three
  // cases is made by the input VectorOperation that either reads values from
  // a vector and puts the data into the local data field or write local data
  // into the vector. Certain operations are no-ops for the given use case.

  std::integral_constant<bool,
                         internal::is_vectorizable<VectorType, Number>::value>
                                                               vector_selector;
  const internal::MatrixFreeFunctions::DoFInfo::DoFAccessIndex ind =
    is_face ? dof_access_index :
              internal::MatrixFreeFunctions::DoFInfo::dof_access_cell;
  const unsigned int n_lanes = mask.count();

  const std::vector<unsigned int> &dof_indices_cont =
    dof_info->dof_indices_contiguous[ind];

  // Simple case: We have contiguous storage, so we can simply copy out the
  // data
  if (dof_info->index_storage_variants[ind][cell] ==
        internal::MatrixFreeFunctions::DoFInfo::IndexStorageVariants::
          interleaved_contiguous &&
      n_lanes == VectorizedArrayType::n_array_elements)
    {
      const unsigned int dof_index =
        dof_indices_cont[cell * VectorizedArrayType::n_array_elements] +
        dof_info->component_dof_indices_offset[active_fe_index]
                                              [first_selected_component] *
          VectorizedArrayType::n_array_elements;
      if (n_components == 1 || n_fe_components == 1)
        for (unsigned int comp = 0; comp < n_components; ++comp)
          operation.process_dofs_vectorized(data->dofs_per_component_on_cell,
                                            dof_index,
                                            *src[comp],
                                            values_dofs[comp],
                                            vector_selector);
      else
        operation.process_dofs_vectorized(data->dofs_per_component_on_cell *
                                            n_components,
                                          dof_index,
                                          *src[0],
                                          values_dofs[0],
                                          vector_selector);
      return;
    }

  // More general case: Must go through the components one by one and apply
  // some transformations
  const unsigned int vectorization_populated =
    dof_info->n_vectorization_lanes_filled[ind][this->cell];

  unsigned int dof_indices[VectorizedArrayType::n_array_elements];
  for (unsigned int v = 0; v < vectorization_populated; ++v)
    dof_indices[v] =
      dof_indices_cont[cell * VectorizedArrayType::n_array_elements + v] +
      dof_info->component_dof_indices_offset[active_fe_index]
                                            [first_selected_component] *
        dof_info->dof_indices_interleave_strides
          [ind][cell * VectorizedArrayType::n_array_elements + v];

  for (unsigned int v = vectorization_populated;
       v < VectorizedArrayType::n_array_elements;
       ++v)
    dof_indices[v] = numbers::invalid_unsigned_int;

  // In the case with contiguous cell indices, we know that there are no
  // constraints and that the indices within each element are contiguous
  if (vectorization_populated == VectorizedArrayType::n_array_elements &&
      n_lanes == VectorizedArrayType::n_array_elements)
    {
      if (dof_info->index_storage_variants[ind][cell] ==
          internal::MatrixFreeFunctions::DoFInfo::IndexStorageVariants::
            contiguous)
        {
          if (n_components == 1 || n_fe_components == 1)
            for (unsigned int comp = 0; comp < n_components; ++comp)
              operation.process_dofs_vectorized_transpose(
                data->dofs_per_component_on_cell,
                dof_indices,
                *src[comp],
                values_dofs[comp],
                vector_selector);
          else
            operation.process_dofs_vectorized_transpose(
              data->dofs_per_component_on_cell * n_components,
              dof_indices,
              *src[0],
              &values_dofs[0][0],
              vector_selector);
        }
      else if (dof_info->index_storage_variants[ind][cell] ==
               internal::MatrixFreeFunctions::DoFInfo::IndexStorageVariants::
                 interleaved_contiguous_strided)
        {
          if (n_components == 1 || n_fe_components == 1)
            for (unsigned int i = 0; i < data->dofs_per_component_on_cell; ++i)
              {
                for (unsigned int comp = 0; comp < n_components; ++comp)
                  operation.process_dof_gather(
                    dof_indices,
                    *src[comp],
                    i * VectorizedArrayType::n_array_elements,
                    values_dofs[comp][i],
                    vector_selector);
              }
          else
            for (unsigned int comp = 0; comp < n_components; ++comp)
              for (unsigned int i = 0; i < data->dofs_per_component_on_cell;
                   ++i)
                {
                  operation.process_dof_gather(
                    dof_indices,
                    *src[0],
                    (comp * data->dofs_per_component_on_cell + i) *
                      VectorizedArrayType::n_array_elements,
                    values_dofs[comp][i],
                    vector_selector);
                }
        }
      else
        {
          Assert(dof_info->index_storage_variants[ind][cell] ==
                   internal::MatrixFreeFunctions::DoFInfo::
                     IndexStorageVariants::interleaved_contiguous_mixed_strides,
                 ExcNotImplemented());
          const unsigned int *offsets =
            &dof_info->dof_indices_interleave_strides
               [ind][VectorizedArrayType::n_array_elements * cell];
          if (n_components == 1 || n_fe_components == 1)
            for (unsigned int i = 0; i < data->dofs_per_component_on_cell; ++i)
              {
                for (unsigned int comp = 0; comp < n_components; ++comp)
                  operation.process_dof_gather(dof_indices,
                                               *src[comp],
                                               0,
                                               values_dofs[comp][i],
                                               vector_selector);
                DEAL_II_OPENMP_SIMD_PRAGMA
                for (unsigned int v = 0;
                     v < VectorizedArrayType::n_array_elements;
                     ++v)
                  dof_indices[v] += offsets[v];
              }
          else
            for (unsigned int comp = 0; comp < n_components; ++comp)
              for (unsigned int i = 0; i < data->dofs_per_component_on_cell;
                   ++i)
                {
                  operation.process_dof_gather(dof_indices,
                                               *src[0],
                                               0,
                                               values_dofs[comp][i],
                                               vector_selector);
                  DEAL_II_OPENMP_SIMD_PRAGMA
                  for (unsigned int v = 0;
                       v < VectorizedArrayType::n_array_elements;
                       ++v)
                    dof_indices[v] += offsets[v];
                }
        }
    }
  else
    for (unsigned int comp = 0; comp < n_components; ++comp)
      {
        for (unsigned int i = 0; i < data->dofs_per_component_on_cell; ++i)
          operation.process_empty(values_dofs[comp][i]);
        if (dof_info->index_storage_variants[ind][cell] ==
            internal::MatrixFreeFunctions::DoFInfo::IndexStorageVariants::
              contiguous)
          {
            if (n_components == 1 || n_fe_components == 1)
              {
                for (unsigned int v = 0; v < vectorization_populated; ++v)
                  if (mask[v] == true)
                    for (unsigned int i = 0;
                         i < data->dofs_per_component_on_cell;
                         ++i)
                      operation.process_dof(dof_indices[v] + i,
                                            *src[comp],
                                            values_dofs[comp][i][v]);
              }
            else
              {
                for (unsigned int v = 0; v < vectorization_populated; ++v)
                  if (mask[v] == true)
                    for (unsigned int i = 0;
                         i < data->dofs_per_component_on_cell;
                         ++i)
                      operation.process_dof(
                        dof_indices[v] + i +
                          comp * data->dofs_per_component_on_cell,
                        *src[0],
                        values_dofs[comp][i][v]);
              }
          }
        else
          {
            const unsigned int *offsets =
              &dof_info->dof_indices_interleave_strides
                 [ind][VectorizedArrayType::n_array_elements * cell];
            for (unsigned int v = 0; v < vectorization_populated; ++v)
              AssertIndexRange(offsets[v],
                               VectorizedArrayType::n_array_elements + 1);
            if (n_components == 1 || n_fe_components == 1)
              for (unsigned int v = 0; v < vectorization_populated; ++v)
                {
                  if (mask[v] == true)
                    for (unsigned int i = 0;
                         i < data->dofs_per_component_on_cell;
                         ++i)
                      operation.process_dof(dof_indices[v] + i * offsets[v],
                                            *src[comp],
                                            values_dofs[comp][i][v]);
                }
            else
              {
                for (unsigned int v = 0; v < vectorization_populated; ++v)
                  if (mask[v] == true)
                    for (unsigned int i = 0;
                         i < data->dofs_per_component_on_cell;
                         ++i)
                      operation.process_dof(
                        dof_indices[v] +
                          (i + comp * data->dofs_per_component_on_cell) *
                            offsets[v],
                        *src[0],
                        values_dofs[comp][i][v]);
              }
          }
      }
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
template <typename VectorType>
inline void
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  read_dof_values(const VectorType &src, const unsigned int first_index)
{
  // select between block vectors and non-block vectors. Note that the number
  // of components is checked in the internal data
  typename internal::BlockVectorSelector<
    VectorType,
    IsBlockVector<VectorType>::value>::BaseVectorType *src_data[n_components];
  for (unsigned int d = 0; d < n_components; ++d)
    src_data[d] =
      internal::BlockVectorSelector<VectorType,
                                    IsBlockVector<VectorType>::value>::
        get_vector_component(const_cast<VectorType &>(src), d + first_index);

  internal::VectorReader<Number, VectorizedArrayType> reader;
  read_write_operation(
    reader,
    src_data,
    std::bitset<VectorizedArrayType::n_array_elements>().flip(),
    true);

#  ifdef DEBUG
  dof_values_initialized = true;
#  endif
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
template <typename VectorType>
inline void
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  read_dof_values_plain(const VectorType &src, const unsigned int first_index)
{
  // select between block vectors and non-block vectors. Note that the number
  // of components is checked in the internal data
  typename internal::BlockVectorSelector<
    VectorType,
    IsBlockVector<VectorType>::value>::BaseVectorType *src_data[n_components];
  for (unsigned int d = 0; d < n_components; ++d)
    src_data[d] =
      internal::BlockVectorSelector<VectorType,
                                    IsBlockVector<VectorType>::value>::
        get_vector_component(const_cast<VectorType &>(src), d + first_index);

  internal::VectorReader<Number, VectorizedArrayType> reader;
  read_write_operation(
    reader,
    src_data,
    std::bitset<VectorizedArrayType::n_array_elements>().flip(),
    false);

#  ifdef DEBUG
  dof_values_initialized = true;
#  endif
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
template <typename VectorType>
inline void
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  distribute_local_to_global(
    VectorType &                                              dst,
    const unsigned int                                        first_index,
    const std::bitset<VectorizedArrayType::n_array_elements> &mask) const
{
  Assert(dof_values_initialized == true,
         internal::ExcAccessToUninitializedField());

  // select between block vectors and non-block vectors. Note that the number
  // of components is checked in the internal data
  typename internal::BlockVectorSelector<
    VectorType,
    IsBlockVector<VectorType>::value>::BaseVectorType *dst_data[n_components];
  for (unsigned int d = 0; d < n_components; ++d)
    dst_data[d] = internal::BlockVectorSelector<
      VectorType,
      IsBlockVector<VectorType>::value>::get_vector_component(dst,
                                                              d + first_index);

  internal::VectorDistributorLocalToGlobal<Number, VectorizedArrayType>
    distributor;
  read_write_operation(distributor, dst_data, mask);
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
template <typename VectorType>
inline void
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  set_dof_values(
    VectorType &                                              dst,
    const unsigned int                                        first_index,
    const std::bitset<VectorizedArrayType::n_array_elements> &mask) const
{
  Assert(dof_values_initialized == true,
         internal::ExcAccessToUninitializedField());

  // select between block vectors and non-block vectors. Note that the number
  // of components is checked in the internal data
  typename internal::BlockVectorSelector<
    VectorType,
    IsBlockVector<VectorType>::value>::BaseVectorType *dst_data[n_components];
  for (unsigned int d = 0; d < n_components; ++d)
    dst_data[d] = internal::BlockVectorSelector<
      VectorType,
      IsBlockVector<VectorType>::value>::get_vector_component(dst,
                                                              d + first_index);

  internal::VectorSetter<Number, VectorizedArrayType> setter;
  read_write_operation(setter, dst_data, mask);
}



/*------------------------------ access to data fields ----------------------*/

template <int dim,
          int n_components,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline const std::vector<unsigned int> &
FEEvaluationBase<dim, n_components, Number, is_face, VectorizedArrayType>::
  get_internal_dof_numbering() const
{
  return data->lexicographic_numbering;
}



template <int dim,
          int n_components,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline ArrayView<VectorizedArrayType>
FEEvaluationBase<dim, n_components, Number, is_face, VectorizedArrayType>::
  get_scratch_data() const
{
  return ArrayView<VectorizedArrayType>(
    const_cast<VectorizedArrayType *>(scratch_data),
    scratch_data_array->end() - scratch_data);
}



template <int dim,
          int n_components,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline const VectorizedArrayType *
FEEvaluationBase<dim, n_components, Number, is_face, VectorizedArrayType>::
  begin_dof_values() const
{
  return &values_dofs[0][0];
}



template <int dim,
          int n_components,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline VectorizedArrayType *
FEEvaluationBase<dim, n_components, Number, is_face, VectorizedArrayType>::
  begin_dof_values()
{
#  ifdef DEBUG
  dof_values_initialized = true;
#  endif
  return &values_dofs[0][0];
}



template <int dim,
          int n_components,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline const VectorizedArrayType *
FEEvaluationBase<dim, n_components, Number, is_face, VectorizedArrayType>::
  begin_values() const
{
  Assert(values_quad_initialized || values_quad_submitted, ExcNotInitialized());
  return &values_quad[0][0];
}



template <int dim,
          int n_components,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline VectorizedArrayType *
FEEvaluationBase<dim, n_components, Number, is_face, VectorizedArrayType>::
  begin_values()
{
#  ifdef DEBUG
  values_quad_initialized = true;
  values_quad_submitted   = true;
#  endif
  return &values_quad[0][0];
}



template <int dim,
          int n_components,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline const VectorizedArrayType *
FEEvaluationBase<dim, n_components, Number, is_face, VectorizedArrayType>::
  begin_gradients() const
{
  Assert(gradients_quad_initialized || gradients_quad_submitted,
         ExcNotInitialized());
  return &gradients_quad[0][0][0];
}



template <int dim,
          int n_components,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline VectorizedArrayType *
FEEvaluationBase<dim, n_components, Number, is_face, VectorizedArrayType>::
  begin_gradients()
{
#  ifdef DEBUG
  gradients_quad_submitted   = true;
  gradients_quad_initialized = true;
#  endif
  return &gradients_quad[0][0][0];
}



template <int dim,
          int n_components,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline const VectorizedArrayType *
FEEvaluationBase<dim, n_components, Number, is_face, VectorizedArrayType>::
  begin_hessians() const
{
  Assert(hessians_quad_initialized, ExcNotInitialized());
  return &hessians_quad[0][0][0];
}



template <int dim,
          int n_components,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline VectorizedArrayType *
FEEvaluationBase<dim, n_components, Number, is_face, VectorizedArrayType>::
  begin_hessians()
{
#  ifdef DEBUG
  hessians_quad_initialized = true;
#  endif
  return &hessians_quad[0][0][0];
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE Tensor<1, n_components_, VectorizedArrayType>
                             FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  get_dof_value(const unsigned int dof) const
{
  AssertIndexRange(dof, this->data->dofs_per_component_on_cell);
  Tensor<1, n_components_, VectorizedArrayType> return_value;
  for (unsigned int comp = 0; comp < n_components; comp++)
    return_value[comp] = this->values_dofs[comp][dof];
  return return_value;
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE Tensor<1, n_components_, VectorizedArrayType>
                             FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  get_value(const unsigned int q_point) const
{
  Assert(this->values_quad_initialized == true,
         internal::ExcAccessToUninitializedField());
  AssertIndexRange(q_point, this->n_quadrature_points);
  Tensor<1, n_components_, VectorizedArrayType> return_value;
  for (unsigned int comp = 0; comp < n_components; comp++)
    return_value[comp] = this->values_quad[comp][q_point];
  return return_value;
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE
  Tensor<1, n_components_, Tensor<1, dim, VectorizedArrayType>>
  FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
    get_gradient(const unsigned int q_point) const
{
  Assert(this->gradients_quad_initialized == true,
         internal::ExcAccessToUninitializedField());
  AssertIndexRange(q_point, this->n_quadrature_points);

  Assert(jacobian != nullptr, ExcNotInitialized());

  Tensor<1, n_components_, Tensor<1, dim, VectorizedArrayType>> grad_out;

  // Cartesian cell
  if (!is_face && this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      for (unsigned int comp = 0; comp < n_components; comp++)
        for (unsigned int d = 0; d < dim; ++d)
          grad_out[comp][d] =
            (this->gradients_quad[comp][d][q_point] * jacobian[0][d][d]);
    }
  // cell with general/affine Jacobian
  else
    {
      const Tensor<2, dim, VectorizedArrayType> &jac =
        jacobian[this->cell_type > internal::MatrixFreeFunctions::affine ?
                   q_point :
                   0];
      for (unsigned int comp = 0; comp < n_components; comp++)
        for (unsigned int d = 0; d < dim; ++d)
          {
            grad_out[comp][d] =
              jac[d][0] * this->gradients_quad[comp][0][q_point];
            for (unsigned int e = 1; e < dim; ++e)
              grad_out[comp][d] +=
                jac[d][e] * this->gradients_quad[comp][e][q_point];
          }
    }
  return grad_out;
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE Tensor<1, n_components_, VectorizedArrayType>
                             FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  get_normal_derivative(const unsigned int q_point) const
{
  AssertIndexRange(q_point, this->n_quadrature_points);
  Assert(this->gradients_quad_initialized == true,
         internal::ExcAccessToUninitializedField());

  Assert(normal_x_jacobian != nullptr, ExcNotInitialized());

  Tensor<1, n_components, VectorizedArrayType> grad_out;
  if (this->cell_type == internal::MatrixFreeFunctions::cartesian)
    for (unsigned int comp = 0; comp < n_components; comp++)
      grad_out[comp] = this->gradients_quad[comp][dim - 1][q_point] *
                       (this->normal_x_jacobian[0][dim - 1]);
  else
    {
      const unsigned int index =
        this->cell_type <= internal::MatrixFreeFunctions::affine ? 0 : q_point;
      for (unsigned int comp = 0; comp < n_components; comp++)
        {
          grad_out[comp] = this->gradients_quad[comp][0][q_point] *
                           this->normal_x_jacobian[index][0];
          for (unsigned int d = 1; d < dim; ++d)
            grad_out[comp] += this->gradients_quad[comp][d][q_point] *
                              this->normal_x_jacobian[index][d];
        }
    }
  return grad_out;
}



namespace internal
{
  // compute tmp = hess_unit(u) * J^T. do this manually because we do not
  // store the lower diagonal because of symmetry
  template <typename VectorizedArrayType>
  inline void
  hessian_unit_times_jac(const Tensor<2, 1, VectorizedArrayType> &jac,
                         const VectorizedArrayType *const hessians_quad[1],
                         const unsigned int               q_point,
                         VectorizedArrayType (&tmp)[1][1])
  {
    tmp[0][0] = jac[0][0] * hessians_quad[0][q_point];
  }

  template <typename VectorizedArrayType>
  inline void
  hessian_unit_times_jac(const Tensor<2, 2, VectorizedArrayType> &jac,
                         const VectorizedArrayType *const hessians_quad[3],
                         const unsigned int               q_point,
                         VectorizedArrayType (&tmp)[2][2])
  {
    for (unsigned int d = 0; d < 2; ++d)
      {
        tmp[0][d] = (jac[d][0] * hessians_quad[0][q_point] +
                     jac[d][1] * hessians_quad[2][q_point]);
        tmp[1][d] = (jac[d][0] * hessians_quad[2][q_point] +
                     jac[d][1] * hessians_quad[1][q_point]);
      }
  }

  template <typename VectorizedArrayType>
  inline void
  hessian_unit_times_jac(const Tensor<2, 3, VectorizedArrayType> &jac,
                         const VectorizedArrayType *const hessians_quad[6],
                         const unsigned int               q_point,
                         VectorizedArrayType (&tmp)[3][3])
  {
    for (unsigned int d = 0; d < 3; ++d)
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
} // namespace internal



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline Tensor<1, n_components_, Tensor<2, dim, VectorizedArrayType>>
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  get_hessian(const unsigned int q_point) const
{
  Assert(!is_face, ExcNotImplemented());
  Assert(this->hessians_quad_initialized == true,
         internal::ExcAccessToUninitializedField());
  AssertIndexRange(q_point, this->n_quadrature_points);

  Assert(jacobian != nullptr, ExcNotImplemented());
  const Tensor<2, dim, VectorizedArrayType> &jac =
    jacobian[this->cell_type <= internal::MatrixFreeFunctions::affine ?
               0 :
               q_point];

  Tensor<2, dim, VectorizedArrayType> hessian_out[n_components];

  // Cartesian cell
  if (this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      for (unsigned int comp = 0; comp < n_components; comp++)
        for (unsigned int d = 0; d < dim; ++d)
          {
            hessian_out[comp][d][d] =
              (this->hessians_quad[comp][d][q_point] * jac[d][d] * jac[d][d]);
            switch (dim)
              {
                case 1:
                  break;
                case 2:
                  hessian_out[comp][0][1] =
                    (this->hessians_quad[comp][2][q_point] * jac[0][0] *
                     jac[1][1]);
                  break;
                case 3:
                  hessian_out[comp][0][1] =
                    (this->hessians_quad[comp][3][q_point] * jac[0][0] *
                     jac[1][1]);
                  hessian_out[comp][0][2] =
                    (this->hessians_quad[comp][4][q_point] * jac[0][0] *
                     jac[2][2]);
                  hessian_out[comp][1][2] =
                    (this->hessians_quad[comp][5][q_point] * jac[1][1] *
                     jac[2][2]);
                  break;
                default:
                  Assert(false, ExcNotImplemented());
              }
            for (unsigned int e = d + 1; e < dim; ++e)
              hessian_out[comp][e][d] = hessian_out[comp][d][e];
          }
    }
  // cell with general Jacobian, but constant within the cell
  else if (this->cell_type == internal::MatrixFreeFunctions::affine)
    {
      for (unsigned int comp = 0; comp < n_components; comp++)
        {
          // compute laplacian before the gradient because it needs to access
          // unscaled gradient data
          VectorizedArrayType tmp[dim][dim];
          internal::hessian_unit_times_jac(jac,
                                           this->hessians_quad[comp],
                                           q_point,
                                           tmp);

          // compute first part of hessian, J * tmp = J * hess_unit(u) * J^T
          for (unsigned int d = 0; d < dim; ++d)
            for (unsigned int e = d; e < dim; ++e)
              {
                hessian_out[comp][d][e] = jac[d][0] * tmp[0][e];
                for (unsigned int f = 1; f < dim; ++f)
                  hessian_out[comp][d][e] += jac[d][f] * tmp[f][e];
              }

          // no J' * grad(u) part here because the Jacobian is constant
          // throughout the cell and hence, its derivative is zero

          // take symmetric part
          for (unsigned int d = 0; d < dim; ++d)
            for (unsigned int e = d + 1; e < dim; ++e)
              hessian_out[comp][e][d] = hessian_out[comp][d][e];
        }
    }
  // cell with general Jacobian
  else
    {
      const Tensor<1, dim *(dim + 1) / 2, Tensor<1, dim, VectorizedArrayType>>
        &jac_grad =
          mapping_data->jacobian_gradients
            [1 - this->is_interior_face]
            [this->mapping_data->data_index_offsets[this->cell] + q_point];
      for (unsigned int comp = 0; comp < n_components; comp++)
        {
          // compute laplacian before the gradient because it needs to access
          // unscaled gradient data
          VectorizedArrayType tmp[dim][dim];
          internal::hessian_unit_times_jac(jac,
                                           this->hessians_quad[comp],
                                           q_point,
                                           tmp);

          // compute first part of hessian, J * tmp = J * hess_unit(u) * J^T
          for (unsigned int d = 0; d < dim; ++d)
            for (unsigned int e = d; e < dim; ++e)
              {
                hessian_out[comp][d][e] = jac[d][0] * tmp[0][e];
                for (unsigned int f = 1; f < dim; ++f)
                  hessian_out[comp][d][e] += jac[d][f] * tmp[f][e];
              }

          // add diagonal part of J' * grad(u)
          for (unsigned int d = 0; d < dim; ++d)
            for (unsigned int e = 0; e < dim; ++e)
              hessian_out[comp][d][d] +=
                (jac_grad[d][e] * this->gradients_quad[comp][e][q_point]);

          // add off-diagonal part of J' * grad(u)
          for (unsigned int d = 0, count = dim; d < dim; ++d)
            for (unsigned int e = d + 1; e < dim; ++e, ++count)
              for (unsigned int f = 0; f < dim; ++f)
                hessian_out[comp][d][e] +=
                  (jac_grad[count][f] * this->gradients_quad[comp][f][q_point]);

          // take symmetric part
          for (unsigned int d = 0; d < dim; ++d)
            for (unsigned int e = d + 1; e < dim; ++e)
              hessian_out[comp][e][d] = hessian_out[comp][d][e];
        }
    }
  return Tensor<1, n_components_, Tensor<2, dim, VectorizedArrayType>>(
    hessian_out);
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline Tensor<1, n_components_, Tensor<1, dim, VectorizedArrayType>>
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  get_hessian_diagonal(const unsigned int q_point) const
{
  Assert(!is_face, ExcNotImplemented());
  Assert(this->hessians_quad_initialized == true,
         internal::ExcAccessToUninitializedField());
  AssertIndexRange(q_point, this->n_quadrature_points);

  Assert(jacobian != nullptr, ExcNotImplemented());
  const Tensor<2, dim, VectorizedArrayType> &jac =
    jacobian[this->cell_type <= internal::MatrixFreeFunctions::affine ?
               0 :
               q_point];

  Tensor<1, n_components_, Tensor<1, dim, VectorizedArrayType>> hessian_out;

  // Cartesian cell
  if (this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      for (unsigned int comp = 0; comp < n_components; comp++)
        for (unsigned int d = 0; d < dim; ++d)
          hessian_out[comp][d] =
            (this->hessians_quad[comp][d][q_point] * jac[d][d] * jac[d][d]);
    }
  // cell with general Jacobian, but constant within the cell
  else if (this->cell_type == internal::MatrixFreeFunctions::affine)
    {
      for (unsigned int comp = 0; comp < n_components; comp++)
        {
          // compute laplacian before the gradient because it needs to access
          // unscaled gradient data
          VectorizedArrayType tmp[dim][dim];
          internal::hessian_unit_times_jac(jac,
                                           this->hessians_quad[comp],
                                           q_point,
                                           tmp);

          // compute only the trace part of hessian, J * tmp = J *
          // hess_unit(u) * J^T
          for (unsigned int d = 0; d < dim; ++d)
            {
              hessian_out[comp][d] = jac[d][0] * tmp[0][d];
              for (unsigned int f = 1; f < dim; ++f)
                hessian_out[comp][d] += jac[d][f] * tmp[f][d];
            }
        }
    }
  // cell with general Jacobian
  else
    {
      const Tensor<1, dim *(dim + 1) / 2, Tensor<1, dim, VectorizedArrayType>>
        &jac_grad =
          mapping_data->jacobian_gradients
            [0][this->mapping_data->data_index_offsets[this->cell] + q_point];
      for (unsigned int comp = 0; comp < n_components; comp++)
        {
          // compute laplacian before the gradient because it needs to access
          // unscaled gradient data
          VectorizedArrayType tmp[dim][dim];
          internal::hessian_unit_times_jac(jac,
                                           this->hessians_quad[comp],
                                           q_point,
                                           tmp);

          // compute only the trace part of hessian, J * tmp = J *
          // hess_unit(u) * J^T
          for (unsigned int d = 0; d < dim; ++d)
            {
              hessian_out[comp][d] = jac[d][0] * tmp[0][d];
              for (unsigned int f = 1; f < dim; ++f)
                hessian_out[comp][d] += jac[d][f] * tmp[f][d];
            }

          for (unsigned int d = 0; d < dim; ++d)
            for (unsigned int e = 0; e < dim; ++e)
              hessian_out[comp][d] +=
                (jac_grad[d][e] * this->gradients_quad[comp][e][q_point]);
        }
    }
  return hessian_out;
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline Tensor<1, n_components_, VectorizedArrayType>
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  get_laplacian(const unsigned int q_point) const
{
  Assert(is_face == false, ExcNotImplemented());
  Assert(this->hessians_quad_initialized == true,
         internal::ExcAccessToUninitializedField());
  AssertIndexRange(q_point, this->n_quadrature_points);

  Tensor<1, n_components_, VectorizedArrayType> laplacian_out;
  const Tensor<1, n_components_, Tensor<1, dim, VectorizedArrayType>>
    hess_diag = get_hessian_diagonal(q_point);
  for (unsigned int comp = 0; comp < n_components; ++comp)
    {
      laplacian_out[comp] = hess_diag[comp][0];
      for (unsigned int d = 1; d < dim; ++d)
        laplacian_out[comp] += hess_diag[comp][d];
    }
  return laplacian_out;
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  submit_dof_value(const Tensor<1, n_components_, VectorizedArrayType> val_in,
                   const unsigned int                                  dof)
{
#  ifdef DEBUG
  this->dof_values_initialized = true;
#  endif
  AssertIndexRange(dof, this->data->dofs_per_component_on_cell);
  for (unsigned int comp = 0; comp < n_components; comp++)
    this->values_dofs[comp][dof] = val_in[comp];
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  submit_value(const Tensor<1, n_components_, VectorizedArrayType> val_in,
               const unsigned int                                  q_point)
{
#  ifdef DEBUG
  Assert(this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  AssertIndexRange(q_point, this->n_quadrature_points);
  Assert(this->J_value != nullptr, ExcNotInitialized());
  this->values_quad_submitted = true;
#  endif

  if (this->cell_type <= internal::MatrixFreeFunctions::affine)
    {
      const VectorizedArrayType JxW = J_value[0] * quadrature_weights[q_point];
      for (unsigned int comp = 0; comp < n_components; ++comp)
        this->values_quad[comp][q_point] = val_in[comp] * JxW;
    }
  else
    {
      const VectorizedArrayType JxW = J_value[q_point];
      for (unsigned int comp = 0; comp < n_components; ++comp)
        this->values_quad[comp][q_point] = val_in[comp] * JxW;
    }
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  submit_gradient(
    const Tensor<1, n_components_, Tensor<1, dim, VectorizedArrayType>> grad_in,
    const unsigned int                                                  q_point)
{
#  ifdef DEBUG
  Assert(this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  AssertIndexRange(q_point, this->n_quadrature_points);
  this->gradients_quad_submitted = true;
  Assert(this->J_value != nullptr, ExcNotInitialized());
  Assert(this->jacobian != nullptr, ExcNotInitialized());
#  endif

  if (!is_face && this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      const VectorizedArrayType JxW = J_value[0] * quadrature_weights[q_point];
      for (unsigned int comp = 0; comp < n_components; comp++)
        for (unsigned int d = 0; d < dim; ++d)
          this->gradients_quad[comp][d][q_point] =
            (grad_in[comp][d] * jacobian[0][d][d] * JxW);
    }
  else
    {
      const Tensor<2, dim, VectorizedArrayType> &jac =
        this->cell_type > internal::MatrixFreeFunctions::affine ?
          jacobian[q_point] :
          jacobian[0];
      const VectorizedArrayType JxW =
        this->cell_type > internal::MatrixFreeFunctions::affine ?
          J_value[q_point] :
          J_value[0] * quadrature_weights[q_point];
      for (unsigned int comp = 0; comp < n_components; ++comp)
        for (unsigned int d = 0; d < dim; ++d)
          {
            VectorizedArrayType new_val = jac[0][d] * grad_in[comp][0];
            for (unsigned int e = 1; e < dim; ++e)
              new_val += (jac[e][d] * grad_in[comp][e]);
            this->gradients_quad[comp][d][q_point] = new_val * JxW;
          }
    }
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  submit_normal_derivative(
    const Tensor<1, n_components_, VectorizedArrayType> grad_in,
    const unsigned int                                  q_point)
{
#  ifdef DEBUG
  AssertIndexRange(q_point, this->n_quadrature_points);
  this->gradients_quad_submitted = true;
  Assert(this->normal_x_jacobian != nullptr, ExcNotInitialized());
#  endif

  if (this->cell_type == internal::MatrixFreeFunctions::cartesian)
    for (unsigned int comp = 0; comp < n_components; comp++)
      {
        for (unsigned int d = 0; d < dim - 1; ++d)
          this->gradients_quad[comp][d][q_point] = VectorizedArrayType();
        this->gradients_quad[comp][dim - 1][q_point] =
          grad_in[comp] *
          (this->normal_x_jacobian[0][dim - 1] * this->J_value[0] *
           this->quadrature_weights[q_point]);
      }
  else
    {
      const unsigned int index =
        this->cell_type <= internal::MatrixFreeFunctions::affine ? 0 : q_point;
      for (unsigned int comp = 0; comp < n_components; comp++)
        {
          VectorizedArrayType factor = grad_in[comp] * this->J_value[index];
          if (this->cell_type <= internal::MatrixFreeFunctions::affine)
            factor = factor * this->quadrature_weights[q_point];
          for (unsigned int d = 0; d < dim; ++d)
            this->gradients_quad[comp][d][q_point] =
              factor * this->normal_x_jacobian[index][d];
        }
    }
}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline Tensor<1, n_components_, VectorizedArrayType>
FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>::
  integrate_value() const
{
#  ifdef DEBUG
  Assert(this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  Assert(this->values_quad_submitted == true,
         internal::ExcAccessToUninitializedField());
#  endif
  Tensor<1, n_components_, VectorizedArrayType> return_value;
  for (unsigned int comp = 0; comp < n_components; ++comp)
    return_value[comp] = this->values_quad[comp][0];
  const unsigned int n_q_points = this->n_quadrature_points;
  for (unsigned int q = 1; q < n_q_points; ++q)
    for (unsigned int comp = 0; comp < n_components; ++comp)
      return_value[comp] += this->values_quad[comp][q];
  return (return_value);
}



/*----------------------- FEEvaluationAccess --------------------------------*/


template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline FEEvaluationAccess<dim,
                          n_components_,
                          Number,
                          is_face,
                          VectorizedArrayType>::
  FEEvaluationAccess(
    const MatrixFree<dim, Number, VectorizedArrayType> &data_in,
    const unsigned int                                  dof_no,
    const unsigned int first_selected_component,
    const unsigned int quad_no_in,
    const unsigned int fe_degree,
    const unsigned int n_q_points,
    const bool         is_interior_face)
  : FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>(
      data_in,
      dof_no,
      first_selected_component,
      quad_no_in,
      fe_degree,
      n_q_points,
      is_interior_face)
{}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
template <int n_components_other>
inline FEEvaluationAccess<dim,
                          n_components_,
                          Number,
                          is_face,
                          VectorizedArrayType>::
  FEEvaluationAccess(const Mapping<dim> &      mapping,
                     const FiniteElement<dim> &fe,
                     const Quadrature<1> &     quadrature,
                     const UpdateFlags         update_flags,
                     const unsigned int        first_selected_component,
                     const FEEvaluationBase<dim,
                                            n_components_other,
                                            Number,
                                            is_face,
                                            VectorizedArrayType> *other)
  : FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>(
      mapping,
      fe,
      quadrature,
      update_flags,
      first_selected_component,
      other)
{}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline FEEvaluationAccess<dim,
                          n_components_,
                          Number,
                          is_face,
                          VectorizedArrayType>::
  FEEvaluationAccess(const FEEvaluationAccess<dim,
                                              n_components_,
                                              Number,
                                              is_face,
                                              VectorizedArrayType> &other)
  : FEEvaluationBase<dim, n_components_, Number, is_face, VectorizedArrayType>(
      other)
{}



template <int dim,
          int n_components_,
          typename Number,
          bool is_face,
          typename VectorizedArrayType>
inline FEEvaluationAccess<dim,
                          n_components_,
                          Number,
                          is_face,
                          VectorizedArrayType> &
FEEvaluationAccess<dim, n_components_, Number, is_face, VectorizedArrayType>::
operator=(const FEEvaluationAccess<dim,
                                   n_components_,
                                   Number,
                                   is_face,
                                   VectorizedArrayType> &other)
{
  this->FEEvaluationBase<dim,
                         n_components_,
                         Number,
                         is_face,
                         VectorizedArrayType>::operator=(other);
  return *this;
}



/*-------------------- FEEvaluationAccess scalar ----------------------------*/


template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>::
  FEEvaluationAccess(
    const MatrixFree<dim, Number, VectorizedArrayType> &data_in,
    const unsigned int                                  dof_no,
    const unsigned int first_selected_component,
    const unsigned int quad_no_in,
    const unsigned int fe_degree,
    const unsigned int n_q_points,
    const bool         is_interior_face)
  : FEEvaluationBase<dim, 1, Number, is_face, VectorizedArrayType>(
      data_in,
      dof_no,
      first_selected_component,
      quad_no_in,
      fe_degree,
      n_q_points,
      is_interior_face)
{}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
template <int n_components_other>
inline FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>::
  FEEvaluationAccess(const Mapping<dim> &      mapping,
                     const FiniteElement<dim> &fe,
                     const Quadrature<1> &     quadrature,
                     const UpdateFlags         update_flags,
                     const unsigned int        first_selected_component,
                     const FEEvaluationBase<dim,
                                            n_components_other,
                                            Number,
                                            is_face,
                                            VectorizedArrayType> *other)
  : FEEvaluationBase<dim, 1, Number, is_face, VectorizedArrayType>(
      mapping,
      fe,
      quadrature,
      update_flags,
      first_selected_component,
      other)
{}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>::
  FEEvaluationAccess(
    const FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>
      &other)
  : FEEvaluationBase<dim, 1, Number, is_face, VectorizedArrayType>(other)
{}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType> &
FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>::operator=(
  const FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType> &other)
{
  this->FEEvaluationBase<dim, 1, Number, is_face, VectorizedArrayType>::
  operator=(other);
  return *this;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE VectorizedArrayType
                             FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>::get_dof_value(
  const unsigned int dof) const
{
  AssertIndexRange(dof, this->data->dofs_per_component_on_cell);
  return this->values_dofs[0][dof];
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE VectorizedArrayType
                             FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>::get_value(
  const unsigned int q_point) const
{
  Assert(this->values_quad_initialized == true,
         internal::ExcAccessToUninitializedField());
  AssertIndexRange(q_point, this->n_quadrature_points);
  return this->values_quad[0][q_point];
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE VectorizedArrayType
                             FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>::
  get_normal_derivative(const unsigned int q_point) const
{
  return BaseClass::get_normal_derivative(q_point)[0];
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE Tensor<1, dim, VectorizedArrayType>
                             FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>::get_gradient(
  const unsigned int q_point) const
{
  // could use the base class gradient, but that involves too many expensive
  // initialization operations on tensors

  Assert(this->gradients_quad_initialized == true,
         internal::ExcAccessToUninitializedField());
  AssertIndexRange(q_point, this->n_quadrature_points);

  Assert(this->jacobian != nullptr, ExcNotInitialized());

  Tensor<1, dim, VectorizedArrayType> grad_out;

  if (!is_face && this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      for (unsigned int d = 0; d < dim; ++d)
        grad_out[d] =
          (this->gradients_quad[0][d][q_point] * this->jacobian[0][d][d]);
    }
  // cell with general/affine Jacobian
  else
    {
      const Tensor<2, dim, VectorizedArrayType> &jac =
        this->jacobian[this->cell_type > internal::MatrixFreeFunctions::affine ?
                         q_point :
                         0];
      for (unsigned int d = 0; d < dim; ++d)
        {
          grad_out[d] = jac[d][0] * this->gradients_quad[0][0][q_point];
          for (unsigned int e = 1; e < dim; ++e)
            grad_out[d] += jac[d][e] * this->gradients_quad[0][e][q_point];
        }
    }
  return grad_out;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline Tensor<2, dim, VectorizedArrayType>
FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>::get_hessian(
  const unsigned int q_point) const
{
  return BaseClass::get_hessian(q_point)[0];
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline Tensor<1, dim, VectorizedArrayType>
FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>::
  get_hessian_diagonal(const unsigned int q_point) const
{
  return BaseClass::get_hessian_diagonal(q_point)[0];
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline VectorizedArrayType
FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>::get_laplacian(
  const unsigned int q_point) const
{
  return BaseClass::get_laplacian(q_point)[0];
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline void DEAL_II_ALWAYS_INLINE
            FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>::
  submit_dof_value(const VectorizedArrayType val_in, const unsigned int dof)
{
#  ifdef DEBUG
  this->dof_values_initialized = true;
  AssertIndexRange(dof, this->data->dofs_per_component_on_cell);
#  endif
  this->values_dofs[0][dof] = val_in;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline void DEAL_II_ALWAYS_INLINE
            FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>::submit_value(
  const VectorizedArrayType val_in,
  const unsigned int        q_index)
{
#  ifdef DEBUG
  Assert(this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  AssertIndexRange(q_index, this->n_quadrature_points);
  Assert(this->J_value != nullptr, ExcNotInitialized());
  this->values_quad_submitted = true;
#  endif
  if (this->cell_type <= internal::MatrixFreeFunctions::affine)
    {
      const VectorizedArrayType JxW =
        this->J_value[0] * this->quadrature_weights[q_index];
      this->values_quad[0][q_index] = val_in * JxW;
    }
  else // if (this->cell_type < internal::MatrixFreeFunctions::general)
    {
      this->values_quad[0][q_index] = val_in * this->J_value[q_index];
    }
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>::submit_value(
  const Tensor<1, 1, VectorizedArrayType> val_in,
  const unsigned int                      q_point)
{
  submit_value(val_in[0], q_point);
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>::
  submit_normal_derivative(const VectorizedArrayType grad_in,
                           const unsigned int        q_point)
{
  Tensor<1, 1, VectorizedArrayType> grad;
  grad[0] = grad_in;
  BaseClass::submit_normal_derivative(grad, q_point);
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>::
  submit_gradient(const Tensor<1, dim, VectorizedArrayType> grad_in,
                  const unsigned int                        q_index)
{
#  ifdef DEBUG
  Assert(this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  AssertIndexRange(q_index, this->n_quadrature_points);
  this->gradients_quad_submitted = true;
  Assert(this->J_value != nullptr, ExcNotInitialized());
  Assert(this->jacobian != nullptr, ExcNotInitialized());
#  endif

  if (!is_face && this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      const VectorizedArrayType JxW =
        this->J_value[0] * this->quadrature_weights[q_index];
      for (unsigned int d = 0; d < dim; ++d)
        this->gradients_quad[0][d][q_index] =
          (grad_in[d] * this->jacobian[0][d][d] * JxW);
    }
  // general/affine cell type
  else
    {
      const Tensor<2, dim, VectorizedArrayType> &jac =
        this->cell_type > internal::MatrixFreeFunctions::affine ?
          this->jacobian[q_index] :
          this->jacobian[0];
      const VectorizedArrayType JxW =
        this->cell_type > internal::MatrixFreeFunctions::affine ?
          this->J_value[q_index] :
          this->J_value[0] * this->quadrature_weights[q_index];
      for (unsigned int d = 0; d < dim; ++d)
        {
          VectorizedArrayType new_val = jac[0][d] * grad_in[0];
          for (unsigned int e = 1; e < dim; ++e)
            new_val += jac[e][d] * grad_in[e];
          this->gradients_quad[0][d][q_index] = new_val * JxW;
        }
    }
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline VectorizedArrayType
FEEvaluationAccess<dim, 1, Number, is_face, VectorizedArrayType>::
  integrate_value() const
{
  return BaseClass::integrate_value()[0];
}



/*----------------- FEEvaluationAccess vector-valued ------------------------*/


template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType>::
  FEEvaluationAccess(
    const MatrixFree<dim, Number, VectorizedArrayType> &data_in,
    const unsigned int                                  dof_no,
    const unsigned int first_selected_component,
    const unsigned int quad_no_in,
    const unsigned int fe_degree,
    const unsigned int n_q_points,
    const bool         is_interior_face)
  : FEEvaluationBase<dim, dim, Number, is_face, VectorizedArrayType>(
      data_in,
      dof_no,
      first_selected_component,
      quad_no_in,
      fe_degree,
      n_q_points,
      is_interior_face)
{}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
template <int n_components_other>
inline FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType>::
  FEEvaluationAccess(const Mapping<dim> &      mapping,
                     const FiniteElement<dim> &fe,
                     const Quadrature<1> &     quadrature,
                     const UpdateFlags         update_flags,
                     const unsigned int        first_selected_component,
                     const FEEvaluationBase<dim,
                                            n_components_other,
                                            Number,
                                            is_face,
                                            VectorizedArrayType> *other)
  : FEEvaluationBase<dim, dim, Number, is_face, VectorizedArrayType>(
      mapping,
      fe,
      quadrature,
      update_flags,
      first_selected_component,
      other)
{}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType>::
  FEEvaluationAccess(
    const FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType>
      &other)
  : FEEvaluationBase<dim, dim, Number, is_face, VectorizedArrayType>(other)
{}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType> &
FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType>::operator=(
  const FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType>
    &other)
{
  this->FEEvaluationBase<dim, dim, Number, is_face, VectorizedArrayType>::
  operator=(other);
  return *this;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE Tensor<2, dim, VectorizedArrayType>
                             FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType>::
  get_gradient(const unsigned int q_point) const
{
  return BaseClass::get_gradient(q_point);
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE VectorizedArrayType
                             FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType>::
  get_divergence(const unsigned int q_point) const
{
  Assert(this->gradients_quad_initialized == true,
         internal::ExcAccessToUninitializedField());
  AssertIndexRange(q_point, this->n_quadrature_points);
  Assert(this->jacobian != nullptr, ExcNotInitialized());

  VectorizedArrayType divergence;

  // Cartesian cell
  if (!is_face && this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      divergence =
        (this->gradients_quad[0][0][q_point] * this->jacobian[0][0][0]);
      for (unsigned int d = 1; d < dim; ++d)
        divergence +=
          (this->gradients_quad[d][d][q_point] * this->jacobian[0][d][d]);
    }
  // cell with general/constant Jacobian
  else
    {
      const Tensor<2, dim, VectorizedArrayType> &jac =
        this->cell_type == internal::MatrixFreeFunctions::general ?
          this->jacobian[q_point] :
          this->jacobian[0];
      divergence = (jac[0][0] * this->gradients_quad[0][0][q_point]);
      for (unsigned int e = 1; e < dim; ++e)
        divergence += (jac[0][e] * this->gradients_quad[0][e][q_point]);
      for (unsigned int d = 1; d < dim; ++d)
        for (unsigned int e = 0; e < dim; ++e)
          divergence += (jac[d][e] * this->gradients_quad[d][e][q_point]);
    }
  return divergence;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE SymmetricTensor<2, dim, VectorizedArrayType>
                             FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType>::
  get_symmetric_gradient(const unsigned int q_point) const
{
  // copy from generic function into dim-specialization function
  const Tensor<2, dim, VectorizedArrayType> grad = get_gradient(q_point);
  VectorizedArrayType                       symmetrized[(dim * dim + dim) / 2];
  VectorizedArrayType half = make_vectorized_array<Number>(0.5);
  for (unsigned int d = 0; d < dim; ++d)
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
        Assert(false, ExcNotImplemented());
    }
  return SymmetricTensor<2, dim, VectorizedArrayType>(symmetrized);
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE
  Tensor<1, (dim == 2 ? 1 : dim), VectorizedArrayType>
  FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType>::get_curl(
    const unsigned int q_point) const
{
  // copy from generic function into dim-specialization function
  const Tensor<2, dim, VectorizedArrayType> grad = get_gradient(q_point);
  Tensor<1, (dim == 2 ? 1 : dim), VectorizedArrayType> curl;
  switch (dim)
    {
      case 1:
        Assert(false,
               ExcMessage(
                 "Computing the curl in 1d is not a useful operation"));
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
        Assert(false, ExcNotImplemented());
    }
  return curl;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE Tensor<2, dim, VectorizedArrayType>
                             FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType>::
  get_hessian_diagonal(const unsigned int q_point) const
{
  return BaseClass::get_hessian_diagonal(q_point);
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE Tensor<3, dim, VectorizedArrayType>
                             FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType>::get_hessian(
  const unsigned int q_point) const
{
  Assert(this->hessians_quad_initialized == true,
         internal::ExcAccessToUninitializedField());
  AssertIndexRange(q_point, this->n_quadrature_points);
  return BaseClass::get_hessian(q_point);
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType>::
  submit_gradient(const Tensor<2, dim, VectorizedArrayType> grad_in,
                  const unsigned int                        q_point)
{
  BaseClass::submit_gradient(grad_in, q_point);
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType>::
  submit_gradient(
    const Tensor<1, dim, Tensor<1, dim, VectorizedArrayType>> grad_in,
    const unsigned int                                        q_point)
{
  BaseClass::submit_gradient(grad_in, q_point);
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType>::
  submit_divergence(const VectorizedArrayType div_in,
                    const unsigned int        q_point)
{
#  ifdef DEBUG
  Assert(this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  AssertIndexRange(q_point, this->n_quadrature_points);
  this->gradients_quad_submitted = true;
  Assert(this->J_value != nullptr, ExcNotInitialized());
  Assert(this->jacobian != nullptr, ExcNotInitialized());
#  endif

  if (!is_face && this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      const VectorizedArrayType fac =
        this->J_value[0] * this->quadrature_weights[q_point] * div_in;
      for (unsigned int d = 0; d < dim; ++d)
        {
          this->gradients_quad[d][d][q_point] = (fac * this->jacobian[0][d][d]);
          for (unsigned int e = d + 1; e < dim; ++e)
            {
              this->gradients_quad[d][e][q_point] = VectorizedArrayType();
              this->gradients_quad[e][d][q_point] = VectorizedArrayType();
            }
        }
    }
  else
    {
      const Tensor<2, dim, VectorizedArrayType> &jac =
        this->cell_type == internal::MatrixFreeFunctions::general ?
          this->jacobian[q_point] :
          this->jacobian[0];
      const VectorizedArrayType fac =
        (this->cell_type == internal::MatrixFreeFunctions::general ?
           this->J_value[q_point] :
           this->J_value[0] * this->quadrature_weights[q_point]) *
        div_in;
      for (unsigned int d = 0; d < dim; ++d)
        {
          for (unsigned int e = 0; e < dim; ++e)
            this->gradients_quad[d][e][q_point] = jac[d][e] * fac;
        }
    }
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType>::
  submit_symmetric_gradient(
    const SymmetricTensor<2, dim, VectorizedArrayType> sym_grad,
    const unsigned int                                 q_point)
{
  // could have used base class operator, but that involves some overhead
  // which is inefficient. it is nice to have the symmetric tensor because
  // that saves some operations
#  ifdef DEBUG
  Assert(this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  AssertIndexRange(q_point, this->n_quadrature_points);
  this->gradients_quad_submitted = true;
  Assert(this->J_value != nullptr, ExcNotInitialized());
  Assert(this->jacobian != nullptr, ExcNotInitialized());
#  endif

  if (!is_face && this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      const VectorizedArrayType JxW =
        this->J_value[0] * this->quadrature_weights[q_point];
      for (unsigned int d = 0; d < dim; ++d)
        this->gradients_quad[d][d][q_point] =
          (sym_grad.access_raw_entry(d) * JxW * this->jacobian[0][d][d]);
      for (unsigned int e = 0, counter = dim; e < dim; ++e)
        for (unsigned int d = e + 1; d < dim; ++d, ++counter)
          {
            const VectorizedArrayType value =
              sym_grad.access_raw_entry(counter) * JxW;
            this->gradients_quad[e][d][q_point] =
              (value * this->jacobian[0][d][d]);
            this->gradients_quad[d][e][q_point] =
              (value * this->jacobian[0][e][e]);
          }
    }
  // general/affine cell type
  else
    {
      const VectorizedArrayType JxW =
        this->cell_type == internal::MatrixFreeFunctions::general ?
          this->J_value[q_point] :
          this->J_value[0] * this->quadrature_weights[q_point];
      const Tensor<2, dim, VectorizedArrayType> &jac =
        this->cell_type == internal::MatrixFreeFunctions::general ?
          this->jacobian[q_point] :
          this->jacobian[0];
      VectorizedArrayType weighted[dim][dim];
      for (unsigned int i = 0; i < dim; ++i)
        weighted[i][i] = sym_grad.access_raw_entry(i) * JxW;
      for (unsigned int i = 0, counter = dim; i < dim; ++i)
        for (unsigned int j = i + 1; j < dim; ++j, ++counter)
          {
            const VectorizedArrayType value =
              sym_grad.access_raw_entry(counter) * JxW;
            weighted[i][j] = value;
            weighted[j][i] = value;
          }
      for (unsigned int comp = 0; comp < dim; ++comp)
        for (unsigned int d = 0; d < dim; ++d)
          {
            VectorizedArrayType new_val = jac[0][d] * weighted[comp][0];
            for (unsigned int e = 1; e < dim; ++e)
              new_val += jac[e][d] * weighted[comp][e];
            this->gradients_quad[comp][d][q_point] = new_val;
          }
    }
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationAccess<dim, dim, Number, is_face, VectorizedArrayType>::submit_curl(
  const Tensor<1, dim == 2 ? 1 : dim, VectorizedArrayType> curl,
  const unsigned int                                       q_point)
{
  Tensor<2, dim, VectorizedArrayType> grad;
  switch (dim)
    {
      case 1:
        Assert(false,
               ExcMessage(
                 "Testing by the curl in 1d is not a useful operation"));
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
        Assert(false, ExcNotImplemented());
    }
  submit_gradient(grad, q_point);
}


/*-------------------- FEEvaluationAccess scalar for 1d ---------------------*/


template <typename Number, bool is_face, typename VectorizedArrayType>
inline FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::
  FEEvaluationAccess(const MatrixFree<1, Number, VectorizedArrayType> &data_in,
                     const unsigned int                                dof_no,
                     const unsigned int first_selected_component,
                     const unsigned int quad_no_in,
                     const unsigned int fe_degree,
                     const unsigned int n_q_points,
                     const bool         is_interior_face)
  : FEEvaluationBase<1, 1, Number, is_face, VectorizedArrayType>(
      data_in,
      dof_no,
      first_selected_component,
      quad_no_in,
      fe_degree,
      n_q_points,
      is_interior_face)
{}



template <typename Number, bool is_face, typename VectorizedArrayType>
template <int n_components_other>
inline FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::
  FEEvaluationAccess(const Mapping<1> &      mapping,
                     const FiniteElement<1> &fe,
                     const Quadrature<1> &   quadrature,
                     const UpdateFlags       update_flags,
                     const unsigned int      first_selected_component,
                     const FEEvaluationBase<1,
                                            n_components_other,
                                            Number,
                                            is_face,
                                            VectorizedArrayType> *other)
  : FEEvaluationBase<1, 1, Number, is_face, VectorizedArrayType>(
      mapping,
      fe,
      quadrature,
      update_flags,
      first_selected_component,
      other)
{}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::
  FEEvaluationAccess(
    const FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType> &other)
  : FEEvaluationBase<1, 1, Number, is_face, VectorizedArrayType>(other)
{}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType> &
FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::operator=(
  const FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType> &other)
{
  this->FEEvaluationBase<1, 1, Number, is_face, VectorizedArrayType>::operator=(
    other);
  return *this;
}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE VectorizedArrayType
                             FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::get_dof_value(
  const unsigned int dof) const
{
  AssertIndexRange(dof, this->data->dofs_per_component_on_cell);
  return this->values_dofs[0][dof];
}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE VectorizedArrayType
                             FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::get_value(
  const unsigned int q_point) const
{
  Assert(this->values_quad_initialized == true,
         internal::ExcAccessToUninitializedField());
  AssertIndexRange(q_point, this->n_quadrature_points);
  return this->values_quad[0][q_point];
}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE Tensor<1, 1, VectorizedArrayType>
                             FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::get_gradient(
  const unsigned int q_point) const
{
  // could use the base class gradient, but that involves too many inefficient
  // initialization operations on tensors

  Assert(this->gradients_quad_initialized == true,
         internal::ExcAccessToUninitializedField());
  AssertIndexRange(q_point, this->n_quadrature_points);

  const Tensor<2, 1, VectorizedArrayType> &jac =
    this->cell_type == internal::MatrixFreeFunctions::general ?
      this->jacobian[q_point] :
      this->jacobian[0];

  Tensor<1, 1, VectorizedArrayType> grad_out;
  grad_out[0] = jac[0][0] * this->gradients_quad[0][0][q_point];

  return grad_out;
}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE VectorizedArrayType
                             FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::
  get_normal_derivative(const unsigned int q_point) const
{
  return BaseClass::get_normal_derivative(q_point)[0];
}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE Tensor<2, 1, VectorizedArrayType>
                             FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::get_hessian(
  const unsigned int q_point) const
{
  return BaseClass::get_hessian(q_point)[0];
}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE Tensor<1, 1, VectorizedArrayType>
                             FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::
  get_hessian_diagonal(const unsigned int q_point) const
{
  return BaseClass::get_hessian_diagonal(q_point)[0];
}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE VectorizedArrayType
                             FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::get_laplacian(
  const unsigned int q_point) const
{
  return BaseClass::get_laplacian(q_point)[0];
}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void DEAL_II_ALWAYS_INLINE
                                  FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::
  submit_dof_value(const VectorizedArrayType val_in, const unsigned int dof)
{
#  ifdef DEBUG
  this->dof_values_initialized = true;
  AssertIndexRange(dof, this->data->dofs_per_component_on_cell);
#  endif
  this->values_dofs[0][dof] = val_in;
}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::submit_value(
  const VectorizedArrayType val_in,
  const unsigned int        q_point)
{
#  ifdef DEBUG
  Assert(this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  AssertIndexRange(q_point, this->n_quadrature_points);
  this->values_quad_submitted = true;
#  endif
  if (this->cell_type == internal::MatrixFreeFunctions::general)
    {
      const VectorizedArrayType JxW = this->J_value[q_point];
      this->values_quad[0][q_point] = val_in * JxW;
    }
  else // if (this->cell_type == internal::MatrixFreeFunctions::general)
    {
      const VectorizedArrayType JxW =
        this->J_value[0] * this->quadrature_weights[q_point];
      this->values_quad[0][q_point] = val_in * JxW;
    }
}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::submit_value(
  const Tensor<1, 1, VectorizedArrayType> val_in,
  const unsigned int                      q_point)
{
  submit_value(val_in[0], q_point);
}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::submit_gradient(
  const Tensor<1, 1, VectorizedArrayType> grad_in,
  const unsigned int                      q_point)
{
  submit_gradient(grad_in[0], q_point);
}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::submit_gradient(
  const VectorizedArrayType grad_in,
  const unsigned int        q_point)
{
#  ifdef DEBUG
  Assert(this->cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  AssertIndexRange(q_point, this->n_quadrature_points);
  this->gradients_quad_submitted = true;
#  endif

  const Tensor<2, 1, VectorizedArrayType> &jac =
    this->cell_type == internal::MatrixFreeFunctions::general ?
      this->jacobian[q_point] :
      this->jacobian[0];
  const VectorizedArrayType JxW =
    this->cell_type == internal::MatrixFreeFunctions::general ?
      this->J_value[q_point] :
      this->J_value[0] * this->quadrature_weights[q_point];

  this->gradients_quad[0][0][q_point] = jac[0][0] * grad_in * JxW;
}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::
  submit_normal_derivative(const VectorizedArrayType grad_in,
                           const unsigned int        q_point)
{
  Tensor<1, 1, VectorizedArrayType> grad;
  grad[0] = grad_in;
  BaseClass::submit_normal_derivative(grad, q_point);
}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE void
FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::
  submit_normal_derivative(const Tensor<1, 1, VectorizedArrayType> grad_in,
                           const unsigned int                      q_point)
{
  BaseClass::submit_normal_derivative(grad_in, q_point);
}



template <typename Number, bool is_face, typename VectorizedArrayType>
inline VectorizedArrayType
FEEvaluationAccess<1, 1, Number, is_face, VectorizedArrayType>::
  integrate_value() const
{
  return BaseClass::integrate_value()[0];
}



/*-------------------------- FEEvaluation -----------------------------------*/


template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline FEEvaluation<dim,
                    fe_degree,
                    n_q_points_1d,
                    n_components_,
                    Number,
                    VectorizedArrayType>::
  FEEvaluation(const MatrixFree<dim, Number, VectorizedArrayType> &data_in,
               const unsigned int                                  fe_no,
               const unsigned int                                  quad_no,
               const unsigned int first_selected_component)
  : BaseClass(data_in,
              fe_no,
              first_selected_component,
              quad_no,
              fe_degree,
              static_n_q_points)
  , dofs_per_component(this->data->dofs_per_component_on_cell)
  , dofs_per_cell(this->data->dofs_per_component_on_cell * n_components_)
  , n_q_points(this->data->n_q_points)
{
  check_template_arguments(fe_no, 0);
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline FEEvaluation<dim,
                    fe_degree,
                    n_q_points_1d,
                    n_components_,
                    Number,
                    VectorizedArrayType>::
  FEEvaluation(const Mapping<dim> &      mapping,
               const FiniteElement<dim> &fe,
               const Quadrature<1> &     quadrature,
               const UpdateFlags         update_flags,
               const unsigned int        first_selected_component)
  : BaseClass(mapping,
              fe,
              quadrature,
              update_flags,
              first_selected_component,
              static_cast<
                FEEvaluationBase<dim, 1, Number, false, VectorizedArrayType> *>(
                nullptr))
  , dofs_per_component(this->data->dofs_per_component_on_cell)
  , dofs_per_cell(this->data->dofs_per_component_on_cell * n_components_)
  , n_q_points(this->data->n_q_points)
{
  check_template_arguments(numbers::invalid_unsigned_int, 0);
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline FEEvaluation<dim,
                    fe_degree,
                    n_q_points_1d,
                    n_components_,
                    Number,
                    VectorizedArrayType>::
  FEEvaluation(const FiniteElement<dim> &fe,
               const Quadrature<1> &     quadrature,
               const UpdateFlags         update_flags,
               const unsigned int        first_selected_component)
  : BaseClass(StaticMappingQ1<dim>::mapping,
              fe,
              quadrature,
              update_flags,
              first_selected_component,
              static_cast<
                FEEvaluationBase<dim, 1, Number, false, VectorizedArrayType> *>(
                nullptr))
  , dofs_per_component(this->data->dofs_per_component_on_cell)
  , dofs_per_cell(this->data->dofs_per_component_on_cell * n_components_)
  , n_q_points(this->data->n_q_points)
{
  check_template_arguments(numbers::invalid_unsigned_int, 0);
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
template <int n_components_other>
inline FEEvaluation<dim,
                    fe_degree,
                    n_q_points_1d,
                    n_components_,
                    Number,
                    VectorizedArrayType>::
  FEEvaluation(const FiniteElement<dim> &                   fe,
               const FEEvaluationBase<dim,
                                      n_components_other,
                                      Number,
                                      false,
                                      VectorizedArrayType> &other,
               const unsigned int first_selected_component)
  : BaseClass(other.mapped_geometry->get_fe_values().get_mapping(),
              fe,
              other.mapped_geometry->get_quadrature(),
              other.mapped_geometry->get_fe_values().get_update_flags(),
              first_selected_component,
              &other)
  , dofs_per_component(this->data->dofs_per_component_on_cell)
  , dofs_per_cell(this->data->dofs_per_component_on_cell * n_components_)
  , n_q_points(this->data->n_q_points)
{
  check_template_arguments(numbers::invalid_unsigned_int, 0);
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline FEEvaluation<dim,
                    fe_degree,
                    n_q_points_1d,
                    n_components_,
                    Number,
                    VectorizedArrayType>::FEEvaluation(const FEEvaluation
                                                         &other)
  : BaseClass(other)
  , dofs_per_component(this->data->dofs_per_component_on_cell)
  , dofs_per_cell(this->data->dofs_per_component_on_cell * n_components_)
  , n_q_points(this->data->n_q_points)
{
  check_template_arguments(numbers::invalid_unsigned_int, 0);
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline FEEvaluation<dim,
                    fe_degree,
                    n_q_points_1d,
                    n_components_,
                    Number,
                    VectorizedArrayType> &
FEEvaluation<dim,
             fe_degree,
             n_q_points_1d,
             n_components_,
             Number,
             VectorizedArrayType>::operator=(const FEEvaluation &other)
{
  BaseClass::operator=(other);
  check_template_arguments(numbers::invalid_unsigned_int, 0);
  return *this;
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline void
FEEvaluation<dim,
             fe_degree,
             n_q_points_1d,
             n_components_,
             Number,
             VectorizedArrayType>::
  check_template_arguments(const unsigned int dof_no,
                           const unsigned int first_selected_component)
{
  (void)dof_no;
  (void)first_selected_component;

#  ifdef DEBUG
  // print error message when the dimensions do not match. Propose a possible
  // fix
  if ((static_cast<unsigned int>(fe_degree) != numbers::invalid_unsigned_int &&
       static_cast<unsigned int>(fe_degree) != this->data->fe_degree) ||
      n_q_points != this->n_quadrature_points)
    {
      std::string message =
        "-------------------------------------------------------\n";
      message += "Illegal arguments in constructor/wrong template arguments!\n";
      message += "    Called -->   FEEvaluation<dim,";
      message += Utilities::int_to_string(fe_degree) + ",";
      message += Utilities::int_to_string(n_q_points_1d);
      message += "," + Utilities::int_to_string(n_components);
      message += ",Number>(data";
      if (first_selected_component != numbers::invalid_unsigned_int)
        {
          message += ", " + Utilities::int_to_string(dof_no) + ", ";
          message += Utilities::int_to_string(this->quad_no) + ", ";
          message += Utilities::int_to_string(first_selected_component);
        }
      message += ")\n";

      // check whether some other vector component has the correct number of
      // points
      unsigned int proposed_dof_comp  = numbers::invalid_unsigned_int,
                   proposed_fe_comp   = numbers::invalid_unsigned_int,
                   proposed_quad_comp = numbers::invalid_unsigned_int;
      if (dof_no != numbers::invalid_unsigned_int)
        {
          if (static_cast<unsigned int>(fe_degree) == this->data->fe_degree)
            {
              proposed_dof_comp = dof_no;
              proposed_fe_comp  = first_selected_component;
            }
          else
            for (unsigned int no = 0; no < this->matrix_info->n_components();
                 ++no)
              for (unsigned int nf = 0;
                   nf < this->matrix_info->n_base_elements(no);
                   ++nf)
                if (this->matrix_info
                      ->get_shape_info(no, 0, nf, this->active_fe_index, 0)
                      .fe_degree == static_cast<unsigned int>(fe_degree))
                  {
                    proposed_dof_comp = no;
                    proposed_fe_comp  = nf;
                    break;
                  }
          if (n_q_points ==
              this->mapping_data->descriptor[this->active_quad_index]
                .n_q_points)
            proposed_quad_comp = this->quad_no;
          else
            for (unsigned int no = 0;
                 no < this->matrix_info->get_mapping_info().cell_data.size();
                 ++no)
              if (this->matrix_info->get_mapping_info()
                    .cell_data[no]
                    .descriptor[this->active_quad_index]
                    .n_q_points == n_q_points)
                {
                  proposed_quad_comp = no;
                  break;
                }
        }
      if (proposed_dof_comp != numbers::invalid_unsigned_int &&
          proposed_quad_comp != numbers::invalid_unsigned_int)
        {
          if (proposed_dof_comp != first_selected_component)
            message += "Wrong vector component selection:\n";
          else
            message += "Wrong quadrature formula selection:\n";
          message += "    Did you mean FEEvaluation<dim,";
          message += Utilities::int_to_string(fe_degree) + ",";
          message += Utilities::int_to_string(n_q_points_1d);
          message += "," + Utilities::int_to_string(n_components);
          message += ",Number>(data";
          if (dof_no != numbers::invalid_unsigned_int)
            {
              message +=
                ", " + Utilities::int_to_string(proposed_dof_comp) + ", ";
              message += Utilities::int_to_string(proposed_quad_comp) + ", ";
              message += Utilities::int_to_string(proposed_fe_comp);
            }
          message += ")?\n";
          std::string correct_pos;
          if (proposed_dof_comp != dof_no)
            correct_pos = " ^ ";
          else
            correct_pos = "   ";
          if (proposed_quad_comp != this->quad_no)
            correct_pos += " ^ ";
          else
            correct_pos += "   ";
          if (proposed_fe_comp != first_selected_component)
            correct_pos += " ^\n";
          else
            correct_pos += "  \n";
          message += "                                                     " +
                     correct_pos;
        }
      // ok, did not find the numbers specified by the template arguments in
      // the given list. Suggest correct template arguments
      const unsigned int proposed_n_q_points_1d = static_cast<unsigned int>(
        std::pow(1.001 * this->n_quadrature_points, 1. / dim));
      message += "Wrong template arguments:\n";
      message += "    Did you mean FEEvaluation<dim,";
      message += Utilities::int_to_string(this->data->fe_degree) + ",";
      message += Utilities::int_to_string(proposed_n_q_points_1d);
      message += "," + Utilities::int_to_string(n_components);
      message += ",Number>(data";
      if (dof_no != numbers::invalid_unsigned_int)
        {
          message += ", " + Utilities::int_to_string(dof_no) + ", ";
          message += Utilities::int_to_string(this->quad_no);
          message += ", " + Utilities::int_to_string(first_selected_component);
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

      Assert(static_cast<unsigned int>(fe_degree) == this->data->fe_degree &&
               n_q_points == this->n_quadrature_points,
             ExcMessage(message));
    }
  if (dof_no != numbers::invalid_unsigned_int)
    AssertDimension(
      n_q_points,
      this->mapping_data->descriptor[this->active_quad_index].n_q_points);
#  endif
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline void
FEEvaluation<dim,
             fe_degree,
             n_q_points_1d,
             n_components_,
             Number,
             VectorizedArrayType>::reinit(const unsigned int cell_index)
{
  Assert(this->mapped_geometry == nullptr,
         ExcMessage("FEEvaluation was initialized without a matrix-free object."
                    " Integer indexing is not possible"));
  if (this->mapped_geometry != nullptr)
    return;

  Assert(this->dof_info != nullptr, ExcNotInitialized());
  Assert(this->mapping_data != nullptr, ExcNotInitialized());
  this->cell = cell_index;
  this->cell_type =
    this->matrix_info->get_mapping_info().get_cell_type(cell_index);

  const unsigned int offsets =
    this->mapping_data->data_index_offsets[cell_index];
  this->jacobian = &this->mapping_data->jacobians[0][offsets];
  this->J_value  = &this->mapping_data->JxW_values[offsets];

#  ifdef DEBUG
  this->dof_values_initialized     = false;
  this->values_quad_initialized    = false;
  this->gradients_quad_initialized = false;
  this->hessians_quad_initialized  = false;
#  endif
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
template <typename DoFHandlerType, bool level_dof_access>
inline void
FEEvaluation<dim,
             fe_degree,
             n_q_points_1d,
             n_components_,
             Number,
             VectorizedArrayType>::
  reinit(
    const TriaIterator<DoFCellAccessor<DoFHandlerType, level_dof_access>> &cell)
{
  Assert(this->matrix_info == nullptr,
         ExcMessage("Cannot use initialization from cell iterator if "
                    "initialized from MatrixFree object. Use variant for "
                    "on the fly computation with arguments as for FEValues "
                    "instead"));
  Assert(this->mapped_geometry.get() != nullptr, ExcNotInitialized());
  this->mapped_geometry->reinit(
    static_cast<typename Triangulation<dim>::cell_iterator>(cell));
  this->local_dof_indices.resize(cell->get_fe().dofs_per_cell);
  if (level_dof_access)
    cell->get_mg_dof_indices(this->local_dof_indices);
  else
    cell->get_dof_indices(this->local_dof_indices);
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline void
FEEvaluation<dim,
             fe_degree,
             n_q_points_1d,
             n_components_,
             Number,
             VectorizedArrayType>::
  reinit(const typename Triangulation<dim>::cell_iterator &cell)
{
  Assert(this->matrix_info == 0,
         ExcMessage("Cannot use initialization from cell iterator if "
                    "initialized from MatrixFree object. Use variant for "
                    "on the fly computation with arguments as for FEValues "
                    "instead"));
  Assert(this->mapped_geometry.get() != 0, ExcNotInitialized());
  this->mapped_geometry->reinit(cell);
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline Point<dim, VectorizedArrayType>
FEEvaluation<dim,
             fe_degree,
             n_q_points_1d,
             n_components_,
             Number,
             VectorizedArrayType>::quadrature_point(const unsigned int q) const
{
  if (this->matrix_info == nullptr)
    {
      Assert((this->mapped_geometry->get_fe_values().get_update_flags() |
              update_quadrature_points),
             ExcNotInitialized());
    }
  else
    {
      Assert(this->mapping_data->quadrature_point_offsets.empty() == false,
             ExcNotInitialized());
    }

  AssertIndexRange(q, n_q_points);

  const unsigned int n_q_points_1d_actual =
    fe_degree == -1 ? this->data->n_q_points_1d : n_q_points_1d;

  // Cartesian mesh: not all quadrature points are stored, only the
  // diagonal. Hence, need to find the tensor product index and retrieve the
  // value from that
  const Point<dim, VectorizedArrayType> *quadrature_points =
    &this->mapping_data->quadrature_points
       [this->mapping_data->quadrature_point_offsets[this->cell]];
  if (this->cell_type == internal::MatrixFreeFunctions::cartesian)
    {
      Point<dim, VectorizedArrayType> point;
      switch (dim)
        {
          case 1:
            return quadrature_points[q];
          case 2:
            point[0] = quadrature_points[q % n_q_points_1d_actual][0];
            point[1] = quadrature_points[q / n_q_points_1d_actual][1];
            return point;
          case 3:
            point[0] = quadrature_points[q % n_q_points_1d_actual][0];
            point[1] = quadrature_points[(q / n_q_points_1d_actual) %
                                         n_q_points_1d_actual][1];
            point[2] = quadrature_points[q / (n_q_points_1d_actual *
                                              n_q_points_1d_actual)][2];
            return point;
          default:
            Assert(false, ExcNotImplemented());
            return point;
        }
    }
  // all other cases: just return the respective data as it is fully stored
  else
    return quadrature_points[q];
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline void
FEEvaluation<dim,
             fe_degree,
             n_q_points_1d,
             n_components_,
             Number,
             VectorizedArrayType>::evaluate(const bool evaluate_values,
                                            const bool evaluate_gradients,
                                            const bool evaluate_hessians)
{
  Assert(this->dof_values_initialized == true,
         internal::ExcAccessToUninitializedField());
  evaluate(this->values_dofs[0],
           evaluate_values,
           evaluate_gradients,
           evaluate_hessians);
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline void
FEEvaluation<dim,
             fe_degree,
             n_q_points_1d,
             n_components_,
             Number,
             VectorizedArrayType>::evaluate(const VectorizedArrayType
                                              *        values_array,
                                            const bool evaluate_values,
                                            const bool evaluate_gradients,
                                            const bool evaluate_hessians)
{
  SelectEvaluator<
    dim,
    fe_degree,
    n_q_points_1d,
    n_components,
    VectorizedArrayType>::evaluate(*this->data,
                                   const_cast<VectorizedArrayType *>(
                                     values_array),
                                   this->values_quad[0],
                                   this->gradients_quad[0][0],
                                   this->hessians_quad[0][0],
                                   this->scratch_data,
                                   evaluate_values,
                                   evaluate_gradients,
                                   evaluate_hessians);

#  ifdef DEBUG
  if (evaluate_values == true)
    this->values_quad_initialized = true;
  if (evaluate_gradients == true)
    this->gradients_quad_initialized = true;
  if (evaluate_hessians == true)
    this->hessians_quad_initialized = true;
#  endif
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
template <typename VectorType>
inline void
FEEvaluation<
  dim,
  fe_degree,
  n_q_points_1d,
  n_components_,
  Number,
  VectorizedArrayType>::gather_evaluate(const VectorType &input_vector,
                                        const bool        evaluate_values,
                                        const bool        evaluate_gradients,
                                        const bool        evaluate_hessians)
{
  // If the index storage is interleaved and contiguous and the vector storage
  // has the correct alignment, we can directly pass the pointer into the
  // vector to the evaluate() call, without reading the vector entries into a
  // separate data field. This saves some operations.
  if (std::is_same<typename VectorType::value_type, Number>::value &&
      this->dof_info->index_storage_variants
          [internal::MatrixFreeFunctions::DoFInfo::dof_access_cell]
          [this->cell] == internal::MatrixFreeFunctions::DoFInfo::
                            IndexStorageVariants::interleaved_contiguous &&
      reinterpret_cast<std::size_t>(
        input_vector.begin() +
        this->dof_info->dof_indices_contiguous
          [internal::MatrixFreeFunctions::DoFInfo::dof_access_cell]
          [this->cell * VectorizedArrayType::n_array_elements]) %
          sizeof(VectorizedArrayType) ==
        0)
    {
      const VectorizedArrayType *vec_values =
        reinterpret_cast<const VectorizedArrayType *>(
          input_vector.begin() +
          this->dof_info->dof_indices_contiguous
            [internal::MatrixFreeFunctions::DoFInfo::dof_access_cell]
            [this->cell * VectorizedArrayType::n_array_elements] +
          this->dof_info
              ->component_dof_indices_offset[this->active_fe_index]
                                            [this->first_selected_component] *
            VectorizedArrayType::n_array_elements);

      evaluate(vec_values,
               evaluate_values,
               evaluate_gradients,
               evaluate_hessians);
    }
  else
    {
      this->read_dof_values(input_vector);
      evaluate(this->begin_dof_values(),
               evaluate_values,
               evaluate_gradients,
               evaluate_hessians);
    }
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline void
FEEvaluation<dim,
             fe_degree,
             n_q_points_1d,
             n_components_,
             Number,
             VectorizedArrayType>::integrate(const bool integrate_values,
                                             const bool integrate_gradients)
{
  integrate(integrate_values, integrate_gradients, this->values_dofs[0]);

#  ifdef DEBUG
  this->dof_values_initialized = true;
#  endif
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline void
FEEvaluation<dim,
             fe_degree,
             n_q_points_1d,
             n_components_,
             Number,
             VectorizedArrayType>::integrate(const bool integrate_values,
                                             const bool integrate_gradients,
                                             VectorizedArrayType *values_array)
{
  if (integrate_values == true)
    Assert(this->values_quad_submitted == true,
           internal::ExcAccessToUninitializedField());
  if (integrate_gradients == true)
    Assert(this->gradients_quad_submitted == true,
           internal::ExcAccessToUninitializedField());
  Assert(this->matrix_info != nullptr ||
           this->mapped_geometry->is_initialized(),
         ExcNotInitialized());

  SelectEvaluator<dim,
                  fe_degree,
                  n_q_points_1d,
                  n_components,
                  VectorizedArrayType>::integrate(*this->data,
                                                  values_array,
                                                  this->values_quad[0],
                                                  this->gradients_quad[0][0],
                                                  this->scratch_data,
                                                  integrate_values,
                                                  integrate_gradients,
                                                  false);

#  ifdef DEBUG
  this->dof_values_initialized = true;
#  endif
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
template <typename VectorType>
inline void
FEEvaluation<
  dim,
  fe_degree,
  n_q_points_1d,
  n_components_,
  Number,
  VectorizedArrayType>::integrate_scatter(const bool  integrate_values,
                                          const bool  integrate_gradients,
                                          VectorType &destination)
{
  // If the index storage is interleaved and contiguous and the vector storage
  // has the correct alignment, we can directly pass the pointer into the
  // vector to the integrate() call, without writing temporary results into a
  // separate data field that will later be added into the vector. This saves
  // some operations.
  if (std::is_same<typename VectorType::value_type, Number>::value &&
      this->dof_info->index_storage_variants
          [internal::MatrixFreeFunctions::DoFInfo::dof_access_cell]
          [this->cell] == internal::MatrixFreeFunctions::DoFInfo::
                            IndexStorageVariants::interleaved_contiguous &&
      reinterpret_cast<std::size_t>(
        destination.begin() +
        this->dof_info->dof_indices_contiguous
          [internal::MatrixFreeFunctions::DoFInfo::dof_access_cell]
          [this->cell * VectorizedArrayType::n_array_elements]) %
          sizeof(VectorizedArrayType) ==
        0)
    {
      VectorizedArrayType *vec_values = reinterpret_cast<VectorizedArrayType *>(
        destination.begin() +
        this->dof_info->dof_indices_contiguous
          [internal::MatrixFreeFunctions::DoFInfo::dof_access_cell]
          [this->cell * VectorizedArrayType::n_array_elements] +
        this->dof_info
            ->component_dof_indices_offset[this->active_fe_index]
                                          [this->first_selected_component] *
          VectorizedArrayType::n_array_elements);
      SelectEvaluator<dim,
                      fe_degree,
                      n_q_points_1d,
                      n_components,
                      VectorizedArrayType>::integrate(*this->data,
                                                      vec_values,
                                                      this->values_quad[0],
                                                      this
                                                        ->gradients_quad[0][0],
                                                      this->scratch_data,
                                                      integrate_values,
                                                      integrate_gradients,
                                                      true);
    }
  else
    {
      integrate(integrate_values,
                integrate_gradients,
                this->begin_dof_values());
      this->distribute_local_to_global(destination);
    }
}



/*-------------------------- FEFaceEvaluation ---------------------------*/



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline FEFaceEvaluation<dim,
                        fe_degree,
                        n_q_points_1d,
                        n_components_,
                        Number,
                        VectorizedArrayType>::
  FEFaceEvaluation(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    const bool                                          is_interior_face,
    const unsigned int                                  dof_no,
    const unsigned int                                  quad_no,
    const unsigned int first_selected_component)
  : BaseClass(matrix_free,
              dof_no,
              first_selected_component,
              quad_no,
              fe_degree,
              static_n_q_points,
              is_interior_face)
  , dofs_per_component(this->data->dofs_per_component_on_cell)
  , dofs_per_cell(this->data->dofs_per_component_on_cell * n_components_)
  , n_q_points(this->data->n_q_points_face)
{}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline void
FEFaceEvaluation<dim,
                 fe_degree,
                 n_q_points_1d,
                 n_components_,
                 Number,
                 VectorizedArrayType>::reinit(const unsigned int face_index)
{
  Assert(this->mapped_geometry == nullptr,
         ExcMessage("FEEvaluation was initialized without a matrix-free object."
                    " Integer indexing is not possible"));
  if (this->mapped_geometry != nullptr)
    return;

  this->cell = face_index;
  this->dof_access_index =
    this->is_interior_face ?
      internal::MatrixFreeFunctions::DoFInfo::dof_access_face_interior :
      internal::MatrixFreeFunctions::DoFInfo::dof_access_face_exterior;
  Assert(this->mapping_data != nullptr, ExcNotInitialized());
  const unsigned int n_vectors = VectorizedArrayType::n_array_elements;
  const internal::MatrixFreeFunctions::FaceToCellTopology<n_vectors> &faces =
    this->matrix_info->get_face_info(face_index);
  if (face_index >=
        this->matrix_info->get_task_info().face_partition_data.back() &&
      face_index <
        this->matrix_info->get_task_info().boundary_partition_data.back())
    Assert(this->is_interior_face,
           ExcMessage("Boundary faces do not have a neighbor"));

  this->face_no =
    (this->is_interior_face ? faces.interior_face_no : faces.exterior_face_no);
  this->subface_index = faces.subface_index;
  if (this->is_interior_face == true)
    {
      this->subface_index = GeometryInfo<dim>::max_children_per_cell;
      if (faces.face_orientation > 8)
        this->face_orientation = faces.face_orientation - 8;
      else
        this->face_orientation = 0;
    }
  else
    {
      if (faces.face_orientation < 8)
        this->face_orientation = faces.face_orientation;
      else
        this->face_orientation = 0;
    }

  this->values_quad_submitted = false;

  this->cell_type = this->matrix_info->get_mapping_info().face_type[face_index];
  const unsigned int offsets =
    this->mapping_data->data_index_offsets[face_index];
  this->J_value        = &this->mapping_data->JxW_values[offsets];
  this->normal_vectors = &this->mapping_data->normal_vectors[offsets];
  this->jacobian =
    &this->mapping_data->jacobians[!this->is_interior_face][offsets];
  this->normal_x_jacobian =
    &this->mapping_data
       ->normals_times_jacobians[!this->is_interior_face][offsets];

#  ifdef DEBUG
  this->dof_values_initialized     = false;
  this->values_quad_initialized    = false;
  this->gradients_quad_initialized = false;
  this->hessians_quad_initialized  = false;
#  endif
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline void
FEFaceEvaluation<dim,
                 fe_degree,
                 n_q_points_1d,
                 n_components_,
                 Number,
                 VectorizedArrayType>::reinit(const unsigned int cell_index,
                                              const unsigned int face_number)
{
  Assert(
    this->quad_no <
      this->matrix_info->get_mapping_info().face_data_by_cells.size(),
    ExcMessage(
      "You must set MatrixFree::AdditionalData::mapping_update_flags_faces_by_cells to use the present reinit method."));
  AssertIndexRange(face_number, GeometryInfo<dim>::faces_per_cell);
  AssertIndexRange(cell_index,
                   this->matrix_info->get_mapping_info().cell_type.size());
  Assert(this->mapped_geometry == nullptr,
         ExcMessage("FEEvaluation was initialized without a matrix-free object."
                    " Integer indexing is not possible"));
  Assert(this->is_interior_face == true,
         ExcMessage(
           "Cell-based FEFaceEvaluation::reinit only possible for the "
           "interior face with second argument to constructor as true"));
  if (this->mapped_geometry != nullptr)
    return;
  Assert(this->matrix_info != nullptr, ExcNotInitialized());

  this->cell_type = this->matrix_info->get_mapping_info().cell_type[cell_index];
  this->cell      = cell_index;
  this->face_orientation = 0;
  this->subface_index    = GeometryInfo<dim>::max_children_per_cell;
  this->face_no          = face_number;
  this->dof_access_index =
    internal::MatrixFreeFunctions::DoFInfo::dof_access_cell;

  const unsigned int offsets =
    this->matrix_info->get_mapping_info()
      .face_data_by_cells[this->quad_no]
      .data_index_offsets[cell_index * GeometryInfo<dim>::faces_per_cell +
                          face_number];
  AssertIndexRange(offsets,
                   this->matrix_info->get_mapping_info()
                     .face_data_by_cells[this->quad_no]
                     .JxW_values.size());
  this->J_value = &this->matrix_info->get_mapping_info()
                     .face_data_by_cells[this->quad_no]
                     .JxW_values[offsets];
  this->normal_vectors = &this->matrix_info->get_mapping_info()
                            .face_data_by_cells[this->quad_no]
                            .normal_vectors[offsets];
  this->jacobian = &this->matrix_info->get_mapping_info()
                      .face_data_by_cells[this->quad_no]
                      .jacobians[0][offsets];
  this->normal_x_jacobian = &this->matrix_info->get_mapping_info()
                               .face_data_by_cells[this->quad_no]
                               .normals_times_jacobians[0][offsets];

#  ifdef DEBUG
  this->dof_values_initialized     = false;
  this->values_quad_initialized    = false;
  this->gradients_quad_initialized = false;
  this->hessians_quad_initialized  = false;
#  endif
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components,
          typename Number,
          typename VectorizedArrayType>
inline void
FEFaceEvaluation<dim,
                 fe_degree,
                 n_q_points_1d,
                 n_components,
                 Number,
                 VectorizedArrayType>::evaluate(const bool evaluate_values,
                                                const bool evaluate_gradients)
{
  Assert(this->dof_values_initialized, ExcNotInitialized());

  evaluate(this->values_dofs[0], evaluate_values, evaluate_gradients);
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components,
          typename Number,
          typename VectorizedArrayType>
inline void
FEFaceEvaluation<dim,
                 fe_degree,
                 n_q_points_1d,
                 n_components,
                 Number,
                 VectorizedArrayType>::evaluate(const VectorizedArrayType
                                                  *        values_array,
                                                const bool evaluate_values,
                                                const bool evaluate_gradients)
{
  if (!(evaluate_values + evaluate_gradients))
    return;

  internal::FEFaceEvaluationSelector<
    dim,
    fe_degree,
    n_q_points_1d,
    n_components,
    Number,
    VectorizedArrayType>::evaluate(*this->data,
                                   values_array,
                                   this->begin_values(),
                                   this->begin_gradients(),
                                   this->scratch_data,
                                   evaluate_values,
                                   evaluate_gradients,
                                   this->face_no,
                                   this->subface_index,
                                   this->face_orientation,
                                   this->mapping_data
                                     ->descriptor[this->active_fe_index]
                                     .face_orientations);

#  ifdef DEBUG
  if (evaluate_values == true)
    this->values_quad_initialized = true;
  if (evaluate_gradients == true)
    this->gradients_quad_initialized = true;
#  endif
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components,
          typename Number,
          typename VectorizedArrayType>
inline void
FEFaceEvaluation<dim,
                 fe_degree,
                 n_q_points_1d,
                 n_components,
                 Number,
                 VectorizedArrayType>::integrate(const bool integrate_values,
                                                 const bool integrate_gradients)
{
  integrate(integrate_values, integrate_gradients, this->values_dofs[0]);

#  ifdef DEBUG
  this->dof_values_initialized = true;
#  endif
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components,
          typename Number,
          typename VectorizedArrayType>
inline void
FEFaceEvaluation<dim,
                 fe_degree,
                 n_q_points_1d,
                 n_components,
                 Number,
                 VectorizedArrayType>::integrate(const bool integrate_values,
                                                 const bool integrate_gradients,
                                                 VectorizedArrayType
                                                   *values_array)
{
  if (!(integrate_values + integrate_gradients))
    return;

  internal::FEFaceEvaluationSelector<
    dim,
    fe_degree,
    n_q_points_1d,
    n_components,
    Number,
    VectorizedArrayType>::integrate(*this->data,
                                    values_array,
                                    this->begin_values(),
                                    this->begin_gradients(),
                                    this->scratch_data,
                                    integrate_values,
                                    integrate_gradients,
                                    this->face_no,
                                    this->subface_index,
                                    this->face_orientation,
                                    this->mapping_data
                                      ->descriptor[this->active_fe_index]
                                      .face_orientations);
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
template <typename VectorType>
inline void
FEFaceEvaluation<
  dim,
  fe_degree,
  n_q_points_1d,
  n_components_,
  Number,
  VectorizedArrayType>::gather_evaluate(const VectorType &input_vector,
                                        const bool        evaluate_values,
                                        const bool        evaluate_gradients)
{
  static_assert(internal::has_begin<VectorType>::value &&
                  (std::is_same<decltype(std::declval<VectorType>().begin()),
                                double *>::value ||
                   std::is_same<decltype(std::declval<VectorType>().begin()),
                                float *>::value),
                "This function requires a vector type with begin() function "
                "evaluating to a pointer to basic number (float,double). "
                "Use read_dof_values() followed by evaluate() instead.");

  if (!internal::FEFaceEvaluationSelector<dim,
                                          fe_degree,
                                          n_q_points_1d,
                                          n_components,
                                          Number,
                                          VectorizedArrayType>::
        gather_evaluate(input_vector.begin(),
                        *this->data,
                        *this->dof_info,
                        this->begin_values(),
                        this->begin_gradients(),
                        this->scratch_data,
                        evaluate_values,
                        evaluate_gradients,
                        this->active_fe_index,
                        this->first_selected_component,
                        this->cell,
                        this->face_no,
                        this->subface_index,
                        this->dof_access_index,
                        this->face_orientation,
                        this->mapping_data->descriptor[this->active_fe_index]
                          .face_orientations))
    {
      this->read_dof_values(input_vector);
      this->evaluate(evaluate_values, evaluate_gradients);
    }

#  ifdef DEBUG
  if (evaluate_values == true)
    this->values_quad_initialized = true;
  if (evaluate_gradients == true)
    this->gradients_quad_initialized = true;
#  endif
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
template <typename VectorType>
inline void
FEFaceEvaluation<
  dim,
  fe_degree,
  n_q_points_1d,
  n_components_,
  Number,
  VectorizedArrayType>::integrate_scatter(const bool  integrate_values,
                                          const bool  integrate_gradients,
                                          VectorType &destination)
{
  static_assert(internal::has_begin<VectorType>::value &&
                  (std::is_same<decltype(std::declval<VectorType>().begin()),
                                double *>::value ||
                   std::is_same<decltype(std::declval<VectorType>().begin()),
                                float *>::value),
                "This function requires a vector type with begin() function "
                "evaluating to a pointer to basic number (float,double). "
                "Use integrate() followed by distribute_local_to_global() "
                "instead.");

  if (!internal::FEFaceEvaluationSelector<dim,
                                          fe_degree,
                                          n_q_points_1d,
                                          n_components,
                                          Number,
                                          VectorizedArrayType>::
        integrate_scatter(destination.begin(),
                          *this->data,
                          *this->dof_info,
                          this->begin_dof_values(),
                          this->begin_values(),
                          this->begin_gradients(),
                          this->scratch_data,
                          integrate_values,
                          integrate_gradients,
                          this->active_fe_index,
                          this->first_selected_component,
                          this->cell,
                          this->face_no,
                          this->subface_index,
                          this->dof_access_index,
                          this->face_orientation,
                          this->mapping_data->descriptor[this->active_fe_index]
                            .face_orientations))
    {
      // if we arrive here, writing into the destination vector did not succeed
      // because some of the assumptions in integrate_scatter were not
      // fulfilled (e.g. an element or degree that does not support direct
      // writing), so we must do it here
      this->distribute_local_to_global(destination);
    }
}



template <int dim,
          int fe_degree,
          int n_q_points_1d,
          int n_components_,
          typename Number,
          typename VectorizedArrayType>
inline Point<dim, VectorizedArrayType>
FEFaceEvaluation<dim,
                 fe_degree,
                 n_q_points_1d,
                 n_components_,
                 Number,
                 VectorizedArrayType>::quadrature_point(const unsigned int q)
  const
{
  AssertIndexRange(q, n_q_points);
  if (this->dof_access_index < 2)
    {
      Assert(this->mapping_data->quadrature_point_offsets.empty() == false,
             ExcNotImplemented());
      AssertIndexRange(this->cell,
                       this->mapping_data->quadrature_point_offsets.size());
      return this->mapping_data->quadrature_points
        [this->mapping_data->quadrature_point_offsets[this->cell] + q];
    }
  else
    {
      Assert(this->matrix_info->get_mapping_info()
                 .face_data_by_cells[this->quad_no]
                 .quadrature_point_offsets.empty() == false,
             ExcNotImplemented());
      const unsigned int index =
        this->cell * GeometryInfo<dim>::faces_per_cell + this->face_no;
      AssertIndexRange(index,
                       this->matrix_info->get_mapping_info()
                         .face_data_by_cells[this->quad_no]
                         .quadrature_point_offsets.size());
      return this->matrix_info->get_mapping_info()
        .face_data_by_cells[this->quad_no]
        .quadrature_points[this->matrix_info->get_mapping_info()
                             .face_data_by_cells[this->quad_no]
                             .quadrature_point_offsets[index] +
                           q];
    }
}



/*------------------------- end FEFaceEvaluation ------------------------- */


#endif // ifndef DOXYGEN


DEAL_II_NAMESPACE_CLOSE

#endif
