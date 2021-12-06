// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2021 by the deal.II authors
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


#ifndef dealii_matrix_free_fe_evaluation_base_data_h
#define dealii_matrix_free_fe_evaluation_base_data_h


#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/matrix_free/dof_info.h>
#include <deal.II/matrix_free/mapping_info_storage.h>
#include <deal.II/matrix_free/shape_info.h>
#include <deal.II/matrix_free/type_traits.h>


DEAL_II_NAMESPACE_OPEN



namespace internal
{
  DeclException0(ExcAccessToUninitializedField);

  DeclException1(
    ExcMatrixFreeAccessToUninitializedMappingField,
    std::string,
    << "You are requesting information from an FEEvaluation/FEFaceEvaluation "
    << "object for which this kind of information has not been computed. What "
    << "information these objects compute is determined by the update_* flags you "
    << "pass to MatrixFree::reinit() via MatrixFree::AdditionalData. Here, "
    << "the operation you are attempting requires the <" << arg1
    << "> flag to be set, but it was apparently not specified "
    << "upon initialization.");
} // namespace internal

// forward declarations
template <int dim,
          int fe_degree,
          int n_q_points_1d            = fe_degree + 1,
          int n_components_            = 1,
          typename Number              = double,
          typename VectorizedArrayType = VectorizedArray<Number>>
class FEEvaluation;

template <int dim,
          int n_components_,
          typename Number,
          bool     = false,
          typename = VectorizedArray<Number>>
class FEEvaluationBase;

namespace internal
{
  namespace MatrixFreeFunctions
  {
    template <int, typename, typename>
    class MappingDataOnTheFly;
  }
} // namespace internal


/**
 * This base class of the FEEvaluation and FEFaceEvaluation classes handles
 * mapping-related information independent of the degrees of freedom and
 * finite element in use. This class provides access functionality for user
 * code but is otherwise invisible without any public constructor. The usage
 * is through the class FEEvaluation instead.
 *
 * This class has four template arguments:
 *
 * @tparam dim Dimension in which this class is to be used
 *
 * @tparam Number Number format, usually @p double or @p float
 *
 * @tparam is_face Whether the class is used for a cell integrator (with
 * quadrature dimension the same as the space dimension) or for a face
 * integrator (with quadrature dimension one less)
 *
 * @tparam VectorizedArrayType Type of array to be woked on in a vectorized
 *                             fashion, defaults to VectorizedArray<Number>
 *
 * @note Currently only VectorizedArray<Number, width> is supported as
 *       VectorizedArrayType.
 *
 *
 * @ingroup matrixfree
 */
template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
class FEEvaluationBaseData
{
  static_assert(
    std::is_same<Number, typename VectorizedArrayType::value_type>::value,
    "Type of Number and of VectorizedArrayType do not match.");

  using ShapeInfoType =
    internal::MatrixFreeFunctions::ShapeInfo<VectorizedArrayType>;
  using MappingInfoStorageType =
    internal::MatrixFreeFunctions::MappingInfoStorage<(is_face ? dim - 1 : dim),
                                                      dim,
                                                      Number,
                                                      VectorizedArrayType>;
  using DoFInfo = internal::MatrixFreeFunctions::DoFInfo;

public:
  static constexpr unsigned int dimension = dim;

  /**
   * Copy constructor.
   */
  FEEvaluationBaseData(const FEEvaluationBaseData &other) = default;

  /**
   * Copy assignment operator.
   */
  FEEvaluationBaseData &
  operator=(const FEEvaluationBaseData &other);

  /**
   * @name 1: Access to geometry data at quadrature points
   */
  //@{

  /**
   * Return an ArrayView to internal memory for temporary use. Note that some
   * of this memory is overwritten during evaluate() and integrate() calls so
   * do not assume it to be stable over those calls. The maximum size you can
   * write into is 3*dofs_per_cell+2*n_q_points.
   */
  ArrayView<VectorizedArrayType>
  get_scratch_data() const;

  /**
   * Return the determinant of the Jacobian from the unit to the real cell
   * times the quadrature weight.
   */
  VectorizedArrayType
  JxW(const unsigned int q_point) const;

  /**
   * Return the inverse and transposed version $J^{-\mathrm T}$ of the
   * Jacobian of the mapping between the unit to the real cell defined as
   * $J_{ij} = d x_i / d\hat x_j$. The $(i,j)$ entry of the returned tensor
   * contains $d\hat x_j/dx_i$, i.e., columns refer to reference space
   * coordinates and rows to real cell coordinates. Thus, the returned tensor
   * represents a covariant transformation, which is used in the
   * FEEvaluationBase::get_gradient() function to transform the unit cell
   * gradients to gradients on the real cell by a multiplication $J^{-\mathrm
   * T} \hat{\nabla} u_h$.
   */
  Tensor<2, dim, VectorizedArrayType>
  inverse_jacobian(const unsigned int q_point) const;

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

  //@}

  /**
   * @name 2: Access to internal data arrays
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
   * first component on all quadrature points, then the yy-component,
   * zz-component in (3D), then the xy-component, and so on. Next comes the xx-
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
   * first component on all quadrature points, then the yy-component,
   * zz-component in (3D), then the xy-component, and so on. Next comes the
   * xx-component of the second component, and so on. This is related to the
   * internal data structures used in this class. The raw data after a call to
   * @p evaluate only contains unit cell operations, so possible
   * transformations, quadrature weights etc. must be applied manually. In
   * general, it is safer to use the get_laplacian() or get_hessian()
   * functions instead, which does all the transformation internally.
   */
  VectorizedArrayType *
  begin_hessians();

  //@}

  /**
   * @name 3: Information about the current cell this class operates on
   */
  //@{

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
  const ShapeInfoType &
  get_shape_info() const;

  /**
   * Return a reference to the DoFInfo object currently in use.
   */
  const internal::MatrixFreeFunctions::DoFInfo &
  get_dof_info() const;

  /**
   * Return the numbering of local degrees of freedom within the evaluation
   * routines of FEEvaluation in terms of the standard numbering on finite
   * elements.
   */
  const std::vector<unsigned int> &
  get_internal_dof_numbering() const;

  /**
   * Return the number of the quadrature formula of the present cell.
   */
  unsigned int
  get_quadrature_index() const;

  /**
   * Return index of the current cell or face.
   */
  unsigned int
  get_current_cell_index() const;

  /**
   * Return the active FE index for this class for efficient indexing in the
   * hp-case.
   */
  unsigned int
  get_active_fe_index() const;

  /**
   * Return the active quadrature index for this class for efficient indexing in
   * the hp-case.
   */
  unsigned int
  get_active_quadrature_index() const;

  //@}

protected:
  /**
   * Constructor. Made protected to prevent users from directly using this
   * class. Takes all data stored in MatrixFree. If applied to problems with
   * more than one quadrature formula selected during construction of
   * `matrix_free`, `quad_no` allows to select the appropriate formula.
   */
  FEEvaluationBaseData(const std::tuple<const ShapeInfoType *,
                                        const DoFInfo *,
                                        unsigned int,
                                        unsigned int> &info,
                       const MappingInfoStorageType &  mapping_data,
                       const unsigned int              quad_no,
                       const bool                      is_interior_face,
                       const unsigned int              face_type);

  /**
   * Constructor that comes with reduced functionality and works similar as
   * FEValues.
   */
  FEEvaluationBaseData(
    const std::shared_ptr<
      internal::MatrixFreeFunctions::
        MappingDataOnTheFly<dim, Number, VectorizedArrayType>> &mapping_data);

  /**
   * Sets the pointers for values, gradients, hessians to the central
   * scratch_data_array inside the given scratch array, for a given number of
   * components as provided by one of the derived classes.
   */
  void
  set_data_pointers(AlignedVector<VectorizedArrayType> *scratch_data,
                    const unsigned int                  n_components);

  /**
   * A pointer to the unit cell shape data, i.e., values, gradients and
   * Hessians in 1D at the quadrature points that constitute the tensor
   * product. Also contained in matrix_info, but it simplifies code if we
   * store a reference to it.
   */
  const ShapeInfoType *data;

  /**
   * A pointer to the underlying DoF indices and constraint description
   * for the component specified at construction. Also contained in
   * matrix_info, but it simplifies code if we store a reference to it.
   */
  const DoFInfo *dof_info;

  /**
   * A pointer to the underlying transformation data from unit to real cells
   * for the given quadrature formula specified at construction. Also
   * contained in matrix_info, but it simplifies code if we store a reference
   * to it.
   */
  const MappingInfoStorageType *mapping_data;

  /**
   * The number of the quadrature formula of the present cell among all
   * quadrature formulas available in the MatrixFree objects pass to derived
   * classes.
   */
  const unsigned int quad_no;

  /**
   * The active FE index for this class for efficient indexing in the hp-case.
   */
  const unsigned int active_fe_index;

  /**
   * The active quadrature index for this class for efficient indexing in the
   * hp-case.
   */
  const unsigned int active_quad_index;

  /**
   * A pointer to the underlying quadrature formula specified at construction.
   * Also contained in mapping_data, but it simplifies code if we store a
   * reference to it.
   */
  const typename MappingInfoStorageType::QuadratureDescriptor *descriptor;

  /**
   * The number of quadrature points in the current evaluation context.
   */
  const unsigned int n_quadrature_points;

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
   * This is the user-visible part of FEEvaluationBase::scratch_data_array,
   * only showing the part that can be consumed by various users. It is set
   * during the allocation of the internal data structures in
   * FEEvaluationBase.
   */
  mutable ArrayView<VectorizedArrayType> scratch_data;

  /**
   * This field stores the values for local degrees of freedom (e.g. after
   * reading out from a vector but before applying unit cell transformations
   * or before distributing them into a result vector). The methods
   * get_dof_value() and submit_dof_value() read from or write to this field.
   *
   * The values of this array are stored in the start section of
   * @p scratch_data_array. Due to its access as a thread local memory, the
   * memory can get reused between different calls.
   */
  VectorizedArrayType *values_dofs;

  /**
   * This field stores the values of the finite element function on quadrature
   * points after applying unit cell transformations or before integrating.
   * The methods get_value() and submit_value() access this field.
   *
   * The values of this array are stored in the start section of
   * @p scratch_data_array. Due to its access as a thread local memory, the
   * memory can get reused between different calls.
   */
  VectorizedArrayType *values_quad;

  /**
   * This field stores the gradients of the finite element function on
   * quadrature points after applying unit cell transformations or before
   * integrating. The methods get_gradient() and submit_gradient() (as well as
   * some specializations like get_symmetric_gradient() or get_divergence())
   * access this field.
   *
   * The values of this array are stored in the start section of
   * @p scratch_data_array. Due to its access as a thread local memory, the
   * memory can get reused between different calls.
   */
  VectorizedArrayType *gradients_quad;

  /**
   * This field stores the gradients of the finite element function on
   * quadrature points after applying unit cell transformations or before
   * integrating. The methods get_hessian() and submit_hessian() (as well as
   * some specializations like get_hessian_diagonal() or get_laplacian())
   * access this field for general cell/face types.
   *
   * The values of this array are stored in the start section of
   * @p scratch_data_array. Due to its access as a thread local memory, the
   * memory can get reused between different calls.
   */
  VectorizedArrayType *gradients_from_hessians_quad;

  /**
   * This field stores the Hessians of the finite element function on
   * quadrature points after applying unit cell transformations. The methods
   * get_hessian(), get_laplacian(), get_hessian_diagonal() access this field.
   *
   * The values of this array are stored in the start section of
   * @p scratch_data_array. Due to its access as a thread local memory, the
   * memory can get reused between different calls.
   */
  VectorizedArrayType *hessians_quad;

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
   * Debug information to track whether hessians on quadrature points have
   * been submitted for integration before the integration is actually stared.
   * Used to control exceptions when uninitialized data is used.
   */
  bool hessians_quad_submitted;

  /**
   * After a call to reinit(), stores the number of the cell we are currently
   * working with.
   */
  unsigned int cell;

  /**
   * Flag holding information whether a face is an interior or exterior face
   * according to the defined direction of the normal. For cells it defines if
   * the dof values should be read from the actual cell corresponding to the
   * interior face or the neighboring cell corresponding to the exterior face.
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
   * Geometry data that can be generated FEValues on the fly with the
   * respective constructor, as an alternative to the entry point with
   * MatrixFree.
   */
  std::shared_ptr<internal::MatrixFreeFunctions::
                    MappingDataOnTheFly<dim, Number, VectorizedArrayType>>
    mapped_geometry;

  // Make FEEvaluation and FEEvaluationBase objects friends for access to
  // protected member mapped_geometry.
  template <int, int, typename, bool, typename>
  friend class FEEvaluationBase;

  template <int, int, int, int, typename, typename>
  friend class FEEvaluation;
};



/*----------------------- Inline functions ----------------------------------*/

#ifndef DOXYGEN

template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::
  FEEvaluationBaseData(const std::tuple<const ShapeInfoType *,
                                        const DoFInfo *,
                                        unsigned int,
                                        unsigned int> &shape_dof_info,
                       const MappingInfoStorageType &  mapping_data,
                       const unsigned int              quad_no,
                       const bool                      is_interior_face,
                       const unsigned int              face_type)
  : data(std::get<0>(shape_dof_info))
  , dof_info(std::get<1>(shape_dof_info))
  , mapping_data(&mapping_data)
  , quad_no(quad_no)
  , active_fe_index(std::get<2>(shape_dof_info))
  , active_quad_index(std::get<3>(shape_dof_info))
  , descriptor(
      &mapping_data.descriptor
         [is_face ?
            (active_quad_index * std::max<unsigned int>(1, dim - 1) +
             (face_type == numbers::invalid_unsigned_int ? 0 : face_type)) :
            active_quad_index])
  , n_quadrature_points(descriptor->n_q_points)
  , jacobian(nullptr)
  , J_value(nullptr)
  , normal_vectors(nullptr)
  , normal_x_jacobian(nullptr)
  , quadrature_weights(descriptor->quadrature_weights.begin())
  , dof_values_initialized(false)
  , values_quad_initialized(false)
  , gradients_quad_initialized(false)
  , hessians_quad_initialized(false)
  , values_quad_submitted(false)
  , gradients_quad_submitted(false)
  , cell(numbers::invalid_unsigned_int)
  , is_interior_face(is_interior_face)
  , dof_access_index(
      is_face ?
        (is_interior_face ?
           internal::MatrixFreeFunctions::DoFInfo::dof_access_face_interior :
           internal::MatrixFreeFunctions::DoFInfo::dof_access_face_exterior) :
        internal::MatrixFreeFunctions::DoFInfo::dof_access_cell)
  , face_no(0)
  , face_orientation(0)
  , subface_index(0)
  , cell_type(internal::MatrixFreeFunctions::general)
{}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::
  FEEvaluationBaseData(
    const std::shared_ptr<
      internal::MatrixFreeFunctions::
        MappingDataOnTheFly<dim, Number, VectorizedArrayType>> &mapped_geometry)
  : data(nullptr)
  , dof_info(nullptr)
  , mapping_data(nullptr)
  , quad_no(numbers::invalid_unsigned_int)
  , active_fe_index(numbers::invalid_unsigned_int)
  , active_quad_index(numbers::invalid_unsigned_int)
  , descriptor(nullptr)
  , n_quadrature_points(
      mapped_geometry->get_data_storage().descriptor[0].n_q_points)
  , jacobian(nullptr)
  , J_value(nullptr)
  , normal_vectors(nullptr)
  , normal_x_jacobian(nullptr)
  , quadrature_weights(nullptr)
  , cell(0)
  , cell_type(internal::MatrixFreeFunctions::general)
  , is_interior_face(true)
  , dof_access_index(internal::MatrixFreeFunctions::DoFInfo::dof_access_cell)
  , mapped_geometry(mapped_geometry)
{
  mapping_data = &mapped_geometry->get_data_storage();
  jacobian     = mapped_geometry->get_data_storage().jacobians[0].begin();
  J_value      = mapped_geometry->get_data_storage().JxW_values.begin();
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType> &
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::operator=(
  const FEEvaluationBaseData &other)
{
  AssertDimension(quad_no, other.quad_no);
  AssertDimension(active_fe_index, other.active_fe_index);
  AssertDimension(active_quad_index, other.active_quad_index);
  AssertDimension(n_quadrature_points, descriptor->n_q_points);

  data               = other.data;
  dof_info           = other.dof_info;
  mapping_data       = other.mapping_data;
  descriptor         = other.descriptor;
  jacobian           = nullptr;
  J_value            = nullptr;
  normal_vectors     = nullptr;
  normal_x_jacobian  = nullptr;
  quadrature_weights = other.quadrature_weights;

  dof_values_initialized     = false;
  values_quad_initialized    = false;
  gradients_quad_initialized = false;
  hessians_quad_initialized  = false;
  values_quad_submitted      = false;
  gradients_quad_submitted   = false;

  cell             = numbers::invalid_unsigned_int;
  is_interior_face = other.is_interior_face;
  dof_access_index =
    is_face ?
      (is_interior_face ?
         internal::MatrixFreeFunctions::DoFInfo::dof_access_face_interior :
         internal::MatrixFreeFunctions::DoFInfo::dof_access_face_exterior) :
      internal::MatrixFreeFunctions::DoFInfo::dof_access_cell;
  face_no          = 0;
  face_orientation = 0;
  subface_index    = 0;
  cell_type        = internal::MatrixFreeFunctions::general;

  return *this;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE Tensor<1, dim, VectorizedArrayType>
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::
  get_normal_vector(const unsigned int q_point) const
{
  AssertIndexRange(q_point, n_quadrature_points);
  Assert(normal_vectors != nullptr,
         internal::ExcMatrixFreeAccessToUninitializedMappingField(
           "update_normal_vectors"));
  if (cell_type <= internal::MatrixFreeFunctions::flat_faces)
    return normal_vectors[0];
  else
    return normal_vectors[q_point];
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE VectorizedArrayType
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::JxW(
  const unsigned int q_point) const
{
  AssertIndexRange(q_point, n_quadrature_points);
  Assert(J_value != nullptr,
         internal::ExcMatrixFreeAccessToUninitializedMappingField(
           "update_values|update_gradients"));
  if (cell_type <= internal::MatrixFreeFunctions::affine)
    {
      Assert(quadrature_weights != nullptr, ExcInternalError());
      return J_value[0] * quadrature_weights[q_point];
    }
  else
    return J_value[q_point];
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline DEAL_II_ALWAYS_INLINE Tensor<2, dim, VectorizedArrayType>
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::
  inverse_jacobian(const unsigned int q_point) const
{
  AssertIndexRange(q_point, n_quadrature_points);
  Assert(jacobian != nullptr,
         internal::ExcMatrixFreeAccessToUninitializedMappingField(
           "update_gradients"));
  if (cell_type <= internal::MatrixFreeFunctions::affine)
    return jacobian[0];
  else
    return jacobian[q_point];
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline const VectorizedArrayType *
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::
  begin_dof_values() const
{
  return values_dofs;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline VectorizedArrayType *
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::
  begin_dof_values()
{
#  ifdef DEBUG
  dof_values_initialized = true;
#  endif
  return values_dofs;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline const VectorizedArrayType *
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::begin_values()
  const
{
#  ifdef DEBUG
  Assert(values_quad_initialized || values_quad_submitted, ExcNotInitialized());
#  endif
  return values_quad;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline VectorizedArrayType *
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::begin_values()
{
#  ifdef DEBUG
  values_quad_initialized = true;
  values_quad_submitted   = true;
#  endif
  return values_quad;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline const VectorizedArrayType *
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::
  begin_gradients() const
{
#  ifdef DEBUG
  Assert(gradients_quad_initialized || gradients_quad_submitted,
         ExcNotInitialized());
#  endif
  return gradients_quad;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline VectorizedArrayType *
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::
  begin_gradients()
{
#  ifdef DEBUG
  gradients_quad_submitted   = true;
  gradients_quad_initialized = true;
#  endif
  return gradients_quad;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline const VectorizedArrayType *
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::
  begin_hessians() const
{
#  ifdef DEBUG
  Assert(hessians_quad_initialized, ExcNotInitialized());
#  endif
  return hessians_quad;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline VectorizedArrayType *
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::
  begin_hessians()
{
#  ifdef DEBUG
  hessians_quad_initialized = true;
#  endif
  return hessians_quad;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline unsigned int
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::
  get_mapping_data_index_offset() const
{
  Assert(mapping_data != nullptr, ExcInternalError());

  if (dof_info == nullptr)
    return 0;
  else
    {
      AssertIndexRange(cell, mapping_data->data_index_offsets.size());
      return mapping_data->data_index_offsets[cell];
    }
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline internal::MatrixFreeFunctions::GeometryType
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::get_cell_type()
  const
{
  Assert(cell != numbers::invalid_unsigned_int, ExcNotInitialized());
  return cell_type;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline const internal::MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> &
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::
  get_shape_info() const
{
  Assert(data != nullptr, ExcInternalError());
  return *data;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline const internal::MatrixFreeFunctions::DoFInfo &
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::get_dof_info()
  const
{
  Assert(dof_info != nullptr,
         ExcMessage(
           "FEEvaluation was not initialized with a MatrixFree object!"));
  return *dof_info;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline const std::vector<unsigned int> &
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::
  get_internal_dof_numbering() const
{
  return data->lexicographic_numbering;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline unsigned int
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::
  get_quadrature_index() const
{
  return quad_no;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline unsigned int
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::
  get_current_cell_index() const
{
  if (is_face && dof_access_index ==
                   internal::MatrixFreeFunctions::DoFInfo::dof_access_cell)
    return cell * GeometryInfo<dim>::faces_per_cell + face_no;
  else
    return cell;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline unsigned int
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::
  get_active_fe_index() const
{
  return active_fe_index;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline unsigned int
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::
  get_active_quadrature_index() const
{
  return active_quad_index;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline ArrayView<VectorizedArrayType>
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::
  get_scratch_data() const
{
  return scratch_data;
}



template <int dim, typename Number, bool is_face, typename VectorizedArrayType>
inline void
FEEvaluationBaseData<dim, Number, is_face, VectorizedArrayType>::
  set_data_pointers(AlignedVector<VectorizedArrayType> *scratch_data_array,
                    const unsigned int                  n_components)
{
  Assert(scratch_data_array != nullptr, ExcInternalError());

  const unsigned int tensor_dofs_per_component =
    Utilities::fixed_power<dim>(data->data.front().fe_degree + 1);
  const unsigned int dofs_per_component = data->dofs_per_component_on_cell;

  const unsigned int size_scratch_data =
    std::max(tensor_dofs_per_component + 1, dofs_per_component) * n_components *
      3 +
    2 * n_quadrature_points;
  const unsigned int size_data_arrays =
    n_components * dofs_per_component +
    (n_components * ((dim * (dim + 1)) / 2 + 2 * dim + 1) *
     n_quadrature_points);

  const unsigned int allocated_size = size_scratch_data + size_data_arrays;
  scratch_data_array->resize_fast(allocated_size);
  scratch_data.reinit(scratch_data_array->begin() + size_data_arrays,
                      size_scratch_data);

  // set the pointers to the correct position in the data array
  values_dofs = scratch_data_array->begin();
  values_quad = scratch_data_array->begin() + n_components * dofs_per_component;
  gradients_quad = scratch_data_array->begin() +
                   n_components * (dofs_per_component + n_quadrature_points);
  gradients_from_hessians_quad =
    scratch_data_array->begin() +
    n_components * (dofs_per_component + (dim + 1) * n_quadrature_points);
  hessians_quad =
    scratch_data_array->begin() +
    n_components * (dofs_per_component + (2 * dim + 1) * n_quadrature_points);
}



#endif // ifndef DOXYGEN


DEAL_II_NAMESPACE_CLOSE

#endif
