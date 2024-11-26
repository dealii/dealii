// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_matrix_free_fe_remote_evaluation_h
#define dealii_matrix_free_fe_remote_evaluation_h

#include <deal.II/base/mpi_remote_point_evaluation.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/fe_point_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include <algorithm>
#include <variant>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  /**
   * A class that stores values and/or gradients, e.g. at quadrature points that
   * are accessed inside matrix-free loops.
   */
  template <int dim, int n_components, typename value_type_>
  struct PrecomputedEvaluationData
  {
    using value_type = typename internal::FEPointEvaluation::
      EvaluatorTypeTraits<dim, dim, n_components, value_type_>::value_type;

    using gradient_type = typename internal::FEPointEvaluation::
      EvaluatorTypeTraits<dim, dim, n_components, value_type_>::
        real_gradient_type;

    /**
     * Values at quadrature points.
     */
    AlignedVector<value_type> values;

    /**
     * Gradients at quadrature points.
     */
    AlignedVector<gradient_type> gradients;
  };

  /**
   * A class that stores a CRS like structure to access
   * PrecomputedEvaluationData. In most cases a one-level CRS structure
   * is enough. In this case @c ptrs has to be constructed and the shift can be
   * obtained with `get_shift(index)`. The field @c ptrs_ptrs stays empty.
   * It is only filled if a two level structure is needed. In this case
   * `get_shift(cell_index, face_number)` return the correct shift.
   */
  struct PrecomputedEvaluationDataView
  {
    /**
     * Get a pointer to data at @p index.
     */
    unsigned int
    get_shift(const unsigned int index) const;

    /**
     * Get a pointer to data at (@p cell_index, @p face_number).
     */
    unsigned int
    get_shift(const unsigned int cell_index,
              const unsigned int face_number) const;

    /**
     * Get the number of stored values.
     */
    unsigned int
    size() const;

    /**
     * This parameter can be used if indices do not start with 0.
     */
    unsigned int start = 0;

    /**
     * Pointers to @c ptrs.
     */
    std::vector<unsigned int> ptrs_ptrs;

    /**
     * Pointers to data.
     */
    std::vector<unsigned int> ptrs;
  };

  /**
   * A class helps to access PecomputedEvaluationData in a thread-safe
   * manner. Each thread has to create this wrapper class on its own to
   * avoid race-conditions during reinit().
   */
  template <int dim, int n_components, typename value_type_>
  class PrecomputedEvaluationDataAccessor
  {
    using value_type =
      typename PrecomputedEvaluationData<dim, n_components, value_type_>::
        value_type;
    using gradient_type =
      typename PrecomputedEvaluationData<dim, n_components, value_type_>::
        gradient_type;

  public:
    /**
     * Constructor.
     */
    PrecomputedEvaluationDataAccessor(
      const PrecomputedEvaluationData<dim, n_components, value_type_> &data,
      const PrecomputedEvaluationDataView                             &view);

    /**
     * Get the value at quadrature point @p q. The entity on which the values
     * are defined is set via `reinit()`.
     *
     * @param[in] q Quadrature point at which the value is queried.
     */
    const value_type
    get_value(const unsigned int q) const;

    /**
     * Get the gradients at quadrature point @p q. The entity on which the
     * gradients are defined is set via `reinit()`.
     *
     * @param[in] q Quadrature point at which the gradient is queried.
     */
    const gradient_type
    get_gradient(const unsigned int q) const;

    /**
     * This function has to be called before `get_value()` and/or
     * `get_gradient()`.
     *
     * @param[in] index Entity index at which quadrature points are accessed.
     * This can be, e.g., a cell index, a cell batch index, or a face batch
     * index.
     */
    void
    reinit(const unsigned int index);

    /**
     * Set cell and face number at which quadrature points are accessed.
     *
     * @param[in] index_0 cell index.
     * @param[in] index_1 cell-local face number.
     */
    void
    reinit(const unsigned int index_0, const unsigned int index_1);

  private:
    /**
     * PrecomputedEvaluationDataView provides information where values are
     * located.
     */
    const PrecomputedEvaluationDataView &view;

    /**
     * PrecomputedEvaluationData stores the actual values.
     */
    const PrecomputedEvaluationData<dim, n_components, value_type_> &data;

    /**
     * Offset to data after last call of `reinit()`.
     */
    unsigned int data_offset;
  };

} // namespace internal



/**
 * Communication objects know about the communication pattern.
 * In case of (matrix-free) cells batches or faces batches
 * a RemotePointEvaluation object stores the location of the
 * remote points. @c batch_id_n_entities relates these points
 * to the corresponding quadrature points of entity batches.
 * For this the field stores batch IDs and the number of entities
 * in the batch.
 */
template <int dim>
struct FERemoteCommunicationObjectEntityBatches
{
  /**
   * Object that is reinitialized with the remote points we want to access.
   */
  std::shared_ptr<Utilities::MPI::RemotePointEvaluation<dim>> rpe;

  /**
   * A vector that stores the batch IDs and the number of active entries in the
   * batch related to the the remote points.
   */
  std::vector<std::pair<unsigned int, unsigned int>> batch_id_n_entities;

  /**
   * Function that gives access to @c batch_id_n_entities. The function is only
   * used for a unified access to the vector that stores additional information
   * on the points that are processed via @c rpe.
   */
  std::vector<std::pair<unsigned int, unsigned int>>
  get_communication_object_pntrs() const;
};

/**
 * Similar as @c FERemoteCommunicationObjectEntityBatches.
 * To relate the points from @c RemotePointEvaluation to
 * quadrature points on corresponding entities, entity indices
 * have to be stored.
 */
template <int dim>
struct FERemoteCommunicationObject
{
  /**
   * Object that is reinitialized with the remote points we want to access.
   */
  std::shared_ptr<Utilities::MPI::RemotePointEvaluation<dim>> rpe;

  /**
   * A vector that stores the entity indices related to the the remote points.
   */
  std::vector<unsigned int> indices;

  /**
   * Function that gives access to @c indices. The function is only used
   * for a unified access to the vector that stores additional information on
   * the points that are processed via @c rpe.
   */
  std::vector<unsigned int>
  get_communication_object_pntrs() const;
};

/**
 * Similar as @c FERemoteCommunicationObject.
 * To relate the points from @c RemotePointEvaluation to
 * quadrature points on corresponding faces, cell iterators
 * and face numbers have to be stored.
 */
template <int dim>
struct FERemoteCommunicationObjectTwoLevel
{
  /**
   * Object that is reinitialized with the remote points we want to access.
   */
  std::shared_ptr<Utilities::MPI::RemotePointEvaluation<dim>> rpe;

  /**
   * A vector that stores the cell iterators and the face number related to the
   * the remote points.
   */
  std::vector<
    std::pair<typename Triangulation<dim>::cell_iterator, unsigned int>>
    cell_face_nos;

  /**
   * Function that gives access to @c cell_face_nos. The function is only used
   * for a unified access to the vector that stores additional information on
   * the points that are processed via @c rpe.
   */
  std::vector<
    std::pair<typename Triangulation<dim>::cell_iterator, unsigned int>>
  get_communication_object_pntrs() const;
};

/**
 * A class to fill the fields in PrecomputedEvaluationData.
 */
template <int dim>
class FERemoteEvaluationCommunicator : public EnableObserverPointer
{
public:
  /**
   * This function stores given communication objects
   * and constructs a @c PrecomputedEvaluationDataView
   * object if remote points are related to matrix-free face batches.
   */
  void
  reinit_faces(const std::vector<FERemoteCommunicationObjectEntityBatches<dim>>
                                                           &comm_objects,
               const std::pair<unsigned int, unsigned int> &face_batch_range,
               const std::vector<unsigned int>             &quadrature_sizes);

  /**
   * This function stores given communication objects
   * and constructs a @c PrecomputedEvaluationDataView
   * object if remote points are related to faces.
   */
  void
  reinit_faces(
    const std::vector<FERemoteCommunicationObject<dim>> &comm_objects,
    const std::pair<unsigned int, unsigned int>         &face_range,
    const std::vector<unsigned int>                     &quadrature_sizes);

  /**
   * This function stores given communication objects
   * and constructs a @c PrecomputedEvaluationDataView
   * object if remote points are related to faces of given
   * cells.
   */
  template <typename Iterator>
  void
  reinit_faces(
    const std::vector<FERemoteCommunicationObjectTwoLevel<dim>> &comm_objects,
    const IteratorRange<Iterator>                &cell_iterator_range,
    const std::vector<std::vector<unsigned int>> &quadrature_sizes);



  /**
   * Fill the fields stored in PrecomputedEvaluationData.
   */
  template <int n_components,
            typename PrecomputedEvaluationDataType,
            typename MeshType,
            typename VectorType>
  void
  update_ghost_values(
    PrecomputedEvaluationDataType         &dst,
    const MeshType                        &mesh,
    const VectorType                      &src,
    const EvaluationFlags::EvaluationFlags eval_flags,
    const unsigned int                     first_selected_component,
    const VectorTools::EvaluationFlags::EvaluationFlags vec_flags) const;

  /**
   * Provide access to @c PrecomputedEvaluationDataView.
   */
  const internal::PrecomputedEvaluationDataView &
  get_view() const;

private:
  /**
   * CRS like data structure that describes the data positions at given
   * indices.
   */
  internal::PrecomputedEvaluationDataView view;

  /**
   * A variant for all possible communication objects.
   */
  std::vector<std::variant<FERemoteCommunicationObjectEntityBatches<dim>,
                           FERemoteCommunicationObject<dim>,
                           FERemoteCommunicationObjectTwoLevel<dim>>>
    communication_objects;

  /**
   * Functions in this class only deals with copying data to @c
   * PrecomputedEvaluationData.
   */
  class CopyInstructions
  {
  public:
    /**
     * Copy data from @c src to @c dst. Overload for @c
     * FERemoteCommunicationObject.
     */
    template <typename T1, typename T2>
    static void
    copy_data(const internal::PrecomputedEvaluationDataView &view,
              AlignedVector<T1>                             &dst,
              const std::vector<T2>                         &src,
              const std::vector<unsigned int>               &indices);

    /**
     * Copy data from @c src to @c dst. Overload for @c
     * FERemoteCommunicationObjectTwoLevel.
     */
    template <typename T1, typename T2>
    static void
    copy_data(
      const internal::PrecomputedEvaluationDataView &view,
      AlignedVector<T1>                             &dst,
      const std::vector<T2>                         &src,
      const std::vector<std::pair<typename Triangulation<dim>::cell_iterator,
                                  unsigned int>>    &cell_face_nos);

    /**
     * Copy data from @c src to @c dst. Overload for @c
     * FERemoteCommunicationObjectEntityBatches.
     */
    template <typename T1, typename T2>
    static void
    copy_data(const internal::PrecomputedEvaluationDataView &view,
              AlignedVector<T1>                             &dst,
              const std::vector<T2>                         &src,
              const std::vector<std::pair<unsigned int, unsigned int>>
                &batch_id_n_entities);

  private:
    /**
     * Copy data to the correct position in a @c VectorizedArray.
     */
    template <typename T1, std::size_t n_lanes>
    static void
    copy_data_entries(VectorizedArray<T1, n_lanes> &dst,
                      const unsigned int            v,
                      const T1                     &src);

    /**
     * Similar as @c copy_data_entries() above.
     */
    template <typename T1, int rank_, std::size_t n_lanes, int dim_>
    static void
    copy_data_entries(Tensor<rank_, dim_, VectorizedArray<T1, n_lanes>> &dst,
                      const unsigned int                                 v,
                      const Tensor<rank_, dim_, T1>                     &src);

    /**
     * Similar as @c copy_data_entries() above.
     */
    template <typename T1,
              int         rank_,
              std::size_t n_lanes,
              int         n_components_,
              int         dim_>
    static void
    copy_data_entries(
      Tensor<rank_,
             n_components_,
             Tensor<rank_, dim_, VectorizedArray<T1, n_lanes>>>   &dst,
      const unsigned int                                           v,
      const Tensor<rank_, n_components_, Tensor<rank_, dim_, T1>> &src);

    /**
     * Throw a runtime exception if @c copy_data_entries() has not been
     * implemented for a given type.
     */
    template <typename T1, typename T2>
    static void
    copy_data_entries(T1 &, const unsigned int, const T2 &);
  };
};



/**
 * The namespace collects convenience functions to set up
 * FERemoteEvaluationCommunicator for typical use cases.
 */
namespace Utilities
{
  /**
   * A factory function for the FERemoteEvaluationCommunicator in the case of
   * point-to-point interpolation.
   *
   * @param[in] matrix_free MatrixFree object that is used to distribute
   * quadrature points at non-matching faces. In case of point-to-point
   * interpolation standard quadrature rules are used on faces that are
   * connected to non-matching interfaces.
   * @param[in] non_matching_faces_marked_vertices A vector of boundary face IDs
   * that relate to non-matching interfaces. Each boundary face ID is
   * accompanied by a lambda function that marks the vertices of cells to be
   * considered during the search of remote points (quadrature points of faces
   * with given boundary face ID).
   * @param[in] quad_no Quadrature number in @p matrix_free.
   * @param[in] dof_no DoFHandler number in @p matrix_free.
   * @param[in] tolerance Tolerance to find remote points.
   */
  template <int dim,
            typename Number,
            typename VectorizedArrayType = VectorizedArray<Number>>
  FERemoteEvaluationCommunicator<dim>
  compute_remote_communicator_faces_point_to_point_interpolation(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    const std::vector<
      std::pair<types::boundary_id, std::function<std::vector<bool>()>>>
                      &non_matching_faces_marked_vertices,
    const unsigned int quad_no   = 0,
    const unsigned int dof_no    = 0,
    const double       tolerance = 1e-9);



  /**
   * A factory function for the FERemoteEvaluationCommunicator in the case of
   * Nitsche-type mortaring.
   *
   * @param[in] matrix_free MatrixFree object that holds the DoFHandler,
   * associated with the non-matching grid.
   * @param[in] non_matching_faces_marked_vertices A vector of boundary face IDs
   * that relate to non-matching interfaces. Each boundary face ID is
   * accompanied by a lambda function that marks the vertices of cells to be
   * considered during the process of computing intersections.
   * @param[in] n_q_pnts_1D The number of 1D quadrature points per intersection.
   * @param[in] dof_no DoFHandler number in @p matrix_free.
   * distributed to each intersection (given for one coordinate direction).
   * @param[in] nm_mapping_info In case nm_mapping_info is not a `nullptr` it is
   * set up for cell based loops, such that it can be used in combination with
   * FERemoteEvaluation.
   * @param[in] tolerance Tolerance to find intersections.
   */
  template <int dim,
            typename Number,
            typename VectorizedArrayType = VectorizedArray<Number>>
  FERemoteEvaluationCommunicator<dim>
  compute_remote_communicator_faces_nitsche_type_mortaring(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    const std::vector<
      std::pair<types::boundary_id, std::function<std::vector<bool>()>>>
                      &non_matching_faces_marked_vertices,
    const unsigned int n_q_pnts_1D,
    const unsigned int dof_no                                   = 0,
    NonMatching::MappingInfo<dim, dim, Number> *nm_mapping_info = nullptr,
    const double                                tolerance       = 1e-9);
} // namespace Utilities



/**
 * Class to access data in a matrix-free loop for non-matching discretizations.
 * Interfaces are named with FEEvaluation in mind.
 * The main difference is, that `gather_evaluate()` updates and caches
 * all values at once. Therefore, it has to be called on one thread before a
 * matrix-free loop.
 *
 * To access values and gradients in a thread safe way, @c get_data_accessor()
 * has to be called on every thread. It provides the functions `get_value()` and
 * `get_gradient()`.
 */
template <int dim, int n_components, typename value_type>
class FERemoteEvaluation
{
public:
  /**
   * The constructor needs a corresponding FERemoteEvaluationCommunicator
   * which has to be set up outside of this class. This design choice is
   * motivated since the same FERemoteEvaluationCommunicator can be used
   * for different MeshTypes and number of components.
   *
   * @param[in] comm FERemoteEvaluationCommunicator.
   * @param[in] mesh Triangulation or DoFHandler.
   * @param[in] evaluation_flags Specify treatment of values at points which are
   * found on multiple cells.
   * @param[in] first_selected_component Select first component of evaluation in
   * DoFHandlers with multiple components.
   */
  template <typename MeshType>
  FERemoteEvaluation(const FERemoteEvaluationCommunicator<dim> &comm,
                     const MeshType                            &mesh,
                     const unsigned int first_selected_component = 0,
                     const VectorTools::EvaluationFlags::EvaluationFlags
                       evaluation_flags = VectorTools::EvaluationFlags::avg);

  /**
   * Update the data which can be accessed via `get_value()` and
   * `get_gradient()`.
   *
   * @param[in] src Solution vector used to update data.
   * @param[in] flags Evaluation flags. Currently supported are
   * EvaluationFlags::values and EvaluationFlags::gradients.
   */
  template <typename VectorType>
  void
  gather_evaluate(const VectorType                      &src,
                  const EvaluationFlags::EvaluationFlags flags);

  /**
   * @c FERemoteEvaluation does not provide the functions `get_value()` and
   * `get_gradient()`. To access values and/or gradients call @c
   * get_data_accessor() on every thread, e.g., `auto remote_evaluator =
   * get_data_accessor();` The returned object can be used as follows.
   * @code
   * for (unsigned int entity = range.first; entity < range.second; ++entity)
   * {
   *    remote_evaluator.reinit(entity);
   *    for(unsigned int q : quadrature_point_indices())
   *      remote_evaluator.get_value(q)
   * }
   * @endcode
   */
  internal::PrecomputedEvaluationDataAccessor<dim, n_components, value_type>
  get_data_accessor() const;

private:
  /**
   * Use Triangulation as MeshType.
   */
  void
  set_mesh(const Triangulation<dim> &tria);

  /**
   * Use DoFHandler as MeshType.
   */
  void
  set_mesh(const DoFHandler<dim> &dof_handler);

  /**
   * Precomputed values and/or gradients at remote locations.
   */
  internal::PrecomputedEvaluationData<dim, n_components, value_type> data;

  /**
   * Underlying communicator which handles update of the ghost values.
   */
  ObserverPointer<const FERemoteEvaluationCommunicator<dim>> comm;

  /**
   * Pointer to MeshType if used with Triangulation.
   */
  ObserverPointer<const Triangulation<dim>> tria;

  /**
   * Pointer to MeshType if used with DoFHandler.
   */
  ObserverPointer<const DoFHandler<dim>> dof_handler;

  /**
   * First selected component.
   */
  const unsigned int first_selected_component;

  /**
   * Flags that indicate which ghost values are updated.
   */
  const VectorTools::EvaluationFlags::EvaluationFlags evaluation_flags;
};



namespace internal
{
  unsigned int
  PrecomputedEvaluationDataView::get_shift(const unsigned int index) const
  {
    Assert(ptrs_ptrs.size() == 0, ExcMessage("Two level CRS set up"));

    Assert(index != numbers::invalid_unsigned_int,
           ExcMessage("Index has to be valid!"));

    Assert(start <= index, ExcInternalError());
    AssertIndexRange(index - start, ptrs.size());
    return ptrs[index - start];
  }

  unsigned int
  PrecomputedEvaluationDataView::get_shift(const unsigned int cell_index,
                                           const unsigned int face_number) const
  {
    Assert(ptrs_ptrs.size() > 0, ExcMessage("No two level CRS set up"));

    Assert(cell_index != numbers::invalid_unsigned_int,
           ExcMessage("Cell index has to be valid!"));
    Assert(face_number != numbers::invalid_unsigned_int,
           ExcMessage("Face number has to be valid!"));

    Assert(start <= cell_index, ExcInternalError());

    AssertIndexRange(cell_index - start, ptrs_ptrs.size());
    const unsigned int face_index = ptrs_ptrs[cell_index - start] + face_number;
    AssertIndexRange(face_index, ptrs.size());
    return ptrs[face_index];
  }

  unsigned int
  PrecomputedEvaluationDataView::size() const
  {
    Assert(ptrs.size() > 0, ExcInternalError());
    return ptrs.back();
  }

  template <int dim, int n_components, typename value_type_>
  PrecomputedEvaluationDataAccessor<dim, n_components, value_type_>::
    PrecomputedEvaluationDataAccessor(
      const PrecomputedEvaluationData<dim, n_components, value_type_> &data,
      const PrecomputedEvaluationDataView                             &view)
    : view(view)
    , data(data)
    , data_offset(numbers::invalid_unsigned_int)
  {}

  template <int dim, int n_components, typename value_type_>
  const typename PrecomputedEvaluationData<dim,
                                           n_components,
                                           value_type_>::value_type
  PrecomputedEvaluationDataAccessor<dim, n_components, value_type_>::get_value(
    const unsigned int q) const
  {
    Assert(data_offset != numbers::invalid_unsigned_int,
           ExcMessage("reinit() not called."));
    AssertIndexRange(data_offset + q, data.values.size());
    return data.values[data_offset + q];
  }

  template <int dim, int n_components, typename value_type_>
  const typename PrecomputedEvaluationData<dim, n_components, value_type_>::
    gradient_type
    PrecomputedEvaluationDataAccessor<dim, n_components, value_type_>::
      get_gradient(const unsigned int q) const
  {
    Assert(data_offset != numbers::invalid_unsigned_int,
           ExcMessage("reinit() not called."));
    AssertIndexRange(data_offset + q, data.gradients.size());
    return data.gradients[data_offset + q];
  }

  template <int dim, int n_components, typename value_type_>
  void
  PrecomputedEvaluationDataAccessor<dim, n_components, value_type_>::reinit(
    const unsigned int index)
  {
    data_offset = view.get_shift(index);
  }

  template <int dim, int n_components, typename value_type_>
  void
  PrecomputedEvaluationDataAccessor<dim, n_components, value_type_>::reinit(
    const unsigned int index_0,
    const unsigned int index_1)
  {
    data_offset = view.get_shift(index_0, index_1);
  }

} // namespace internal



template <int dim>
std::vector<std::pair<unsigned int, unsigned int>>
FERemoteCommunicationObjectEntityBatches<dim>::get_communication_object_pntrs()
  const
{
  return batch_id_n_entities;
}

template <int dim>
std::vector<unsigned int>
FERemoteCommunicationObject<dim>::get_communication_object_pntrs() const
{
  return indices;
}

template <int dim>
std::vector<std::pair<typename Triangulation<dim>::cell_iterator, unsigned int>>
FERemoteCommunicationObjectTwoLevel<dim>::get_communication_object_pntrs() const
{
  return cell_face_nos;
}



template <int dim>
void
FERemoteEvaluationCommunicator<dim>::reinit_faces(
  const std::vector<FERemoteCommunicationObjectEntityBatches<dim>>
                                              &comm_objects,
  const std::pair<unsigned int, unsigned int> &face_batch_range,
  const std::vector<unsigned int>             &quadrature_sizes)
{
  // erase type by converting to the base object
  communication_objects.clear();
  for (const auto &co : comm_objects)
    communication_objects.push_back(co);

  // fetch points and update communication patterns
  const unsigned int n_cells = quadrature_sizes.size();
  AssertDimension(n_cells, face_batch_range.second - face_batch_range.first);

  // construct view:
  view.start = face_batch_range.first;

  view.ptrs.resize(n_cells + 1);

  view.ptrs[0] = 0;
  for (unsigned int face = 0; face < n_cells; ++face)
    {
      view.ptrs[face + 1] = view.ptrs[face] + quadrature_sizes[face];
    }
}

template <int dim>
void
FERemoteEvaluationCommunicator<dim>::reinit_faces(
  const std::vector<FERemoteCommunicationObject<dim>> &comm_objects,
  const std::pair<unsigned int, unsigned int>         &face_range,
  const std::vector<unsigned int>                     &quadrature_sizes)
{
  // erase type
  communication_objects.clear();
  for (const auto &co : comm_objects)
    communication_objects.push_back(co);

  const unsigned int n_faces = quadrature_sizes.size();
  AssertDimension(n_faces, face_range.second - face_range.first);

  // construct view:
  view.start = face_range.first;

  view.ptrs.resize(n_faces + 1);

  view.ptrs[0] = 0;
  for (unsigned int face = 0; face < n_faces; ++face)
    view.ptrs[face + 1] = view.ptrs[face] + quadrature_sizes[face];
}

template <int dim>
template <typename Iterator>
void
FERemoteEvaluationCommunicator<dim>::reinit_faces(
  const std::vector<FERemoteCommunicationObjectTwoLevel<dim>> &comm_objects,
  const IteratorRange<Iterator>                &cell_iterator_range,
  const std::vector<std::vector<unsigned int>> &quadrature_sizes)
{
  // erase type
  communication_objects.clear();
  for (const auto &co : comm_objects)
    communication_objects.push_back(co);

  const unsigned int n_cells = quadrature_sizes.size();
  AssertDimension(n_cells,
                  std::distance(cell_iterator_range.begin(),
                                cell_iterator_range.end()));

  // construct view:
  auto &cell_ptrs = view.ptrs_ptrs;
  auto &face_ptrs = view.ptrs;

  view.start = 0;
  cell_ptrs.resize(n_cells);
  unsigned int n_faces = 0;
  for (const auto &cell : cell_iterator_range)
    {
      cell_ptrs[cell->active_cell_index()] = n_faces;
      n_faces += cell->n_faces();
    }

  face_ptrs.resize(n_faces + 1);
  face_ptrs[0] = 0;
  for (const auto &cell : cell_iterator_range)
    {
      for (const auto &f : cell->face_indices())
        {
          const unsigned int face_index =
            cell_ptrs[cell->active_cell_index()] + f;

          face_ptrs[face_index + 1] =
            face_ptrs[face_index] +
            quadrature_sizes[cell->active_cell_index()][f];
        }
    }
}

template <int dim>
template <int n_components,
          typename PrecomputedEvaluationDataType,
          typename MeshType,
          typename VectorType>
void
FERemoteEvaluationCommunicator<dim>::update_ghost_values(
  PrecomputedEvaluationDataType                      &dst,
  const MeshType                                     &mesh,
  const VectorType                                   &src,
  const EvaluationFlags::EvaluationFlags              eval_flags,
  const unsigned int                                  first_selected_component,
  const VectorTools::EvaluationFlags::EvaluationFlags vec_flags) const
{
  const bool has_ghost_elements = src.has_ghost_elements();

  if (has_ghost_elements == false)
    src.update_ghost_values();


  for (const auto &communication_object : communication_objects)
    {
      if (eval_flags & EvaluationFlags::values)
        {
          std::visit(
            [&](const auto &obj) {
              CopyInstructions::copy_data(
                view,
                dst.values,
                VectorTools::point_values<n_components>(
                  *obj.rpe, mesh, src, vec_flags, first_selected_component),
                obj.get_communication_object_pntrs());
            },
            communication_object);
        }

      if (eval_flags & EvaluationFlags::gradients)
        {
          std::visit(
            [&](const auto &obj) {
              CopyInstructions::copy_data(
                view,
                dst.gradients,
                VectorTools::point_gradients<n_components>(
                  *obj.rpe, mesh, src, vec_flags, first_selected_component),
                obj.get_communication_object_pntrs());
            },
            communication_object);
        }

      Assert(!(eval_flags & EvaluationFlags::hessians), ExcNotImplemented());
    }

  if (has_ghost_elements == false)
    src.zero_out_ghost_values();
}

template <int dim>
const internal::PrecomputedEvaluationDataView &
FERemoteEvaluationCommunicator<dim>::get_view() const
{
  return view;
}

template <int dim>
template <typename T1, typename T2>
void
FERemoteEvaluationCommunicator<dim>::CopyInstructions::copy_data(
  const internal::PrecomputedEvaluationDataView &view,
  AlignedVector<T1>                             &dst,
  const std::vector<T2>                         &src,
  const std::vector<unsigned int>               &indices)
{
  dst.resize(view.size());

  unsigned int c = 0;
  for (const auto idx : indices)
    {
      for (unsigned int j = view.get_shift(idx); j < view.get_shift(idx + 1);
           ++j, ++c)
        {
          AssertIndexRange(j, dst.size());
          AssertIndexRange(c, src.size());
          dst[j] = src[c];
        }
    }
}

template <int dim>
template <typename T1, typename T2>
void
FERemoteEvaluationCommunicator<dim>::CopyInstructions::copy_data(
  const internal::PrecomputedEvaluationDataView &view,
  AlignedVector<T1>                             &dst,
  const std::vector<T2>                         &src,
  const std::vector<std::pair<typename Triangulation<dim>::cell_iterator,
                              unsigned int>>    &cell_face_nos)
{
  dst.resize(view.size());

  unsigned int c = 0;
  for (const auto &[cell, f] : cell_face_nos)
    {
      for (unsigned int j = view.get_shift(cell->active_cell_index(), f);
           j < view.get_shift(cell->active_cell_index(), f + 1);
           ++j, ++c)
        {
          AssertIndexRange(j, dst.size());
          AssertIndexRange(c, src.size());

          dst[j] = src[c];
        }
    }
}

template <int dim>
template <typename T1, typename T2>
void
FERemoteEvaluationCommunicator<dim>::CopyInstructions::copy_data(
  const internal::PrecomputedEvaluationDataView            &view,
  AlignedVector<T1>                                        &dst,
  const std::vector<T2>                                    &src,
  const std::vector<std::pair<unsigned int, unsigned int>> &batch_id_n_entities)
{
  dst.resize(view.size());

  unsigned int c = 0;
  for (const auto &[batch_id, n_entries] : batch_id_n_entities)
    {
      for (unsigned int v = 0; v < n_entries; ++v)
        for (unsigned int j = view.get_shift(batch_id);
             j < view.get_shift(batch_id + 1);
             ++j, ++c)
          {
            AssertIndexRange(j, dst.size());
            AssertIndexRange(c, src.size());

            copy_data_entries(dst[j], v, src[c]);
          }
    }
}

template <int dim>
template <typename T1, std::size_t n_lanes>
void
FERemoteEvaluationCommunicator<dim>::CopyInstructions::copy_data_entries(
  VectorizedArray<T1, n_lanes> &dst,
  const unsigned int            v,
  const T1                     &src)
{
  AssertIndexRange(v, n_lanes);

  dst[v] = src;
}

template <int dim>
template <typename T1, int rank_, std::size_t n_lanes, int dim_>
void
FERemoteEvaluationCommunicator<dim>::CopyInstructions::copy_data_entries(
  Tensor<rank_, dim_, VectorizedArray<T1, n_lanes>> &dst,
  const unsigned int                                 v,
  const Tensor<rank_, dim_, T1>                     &src)
{
  AssertIndexRange(v, n_lanes);

  if constexpr (rank_ == 1)
    {
      for (unsigned int i = 0; i < dim_; ++i)
        dst[i][v] = src[i];
    }
  else
    {
      for (unsigned int i = 0; i < rank_; ++i)
        for (unsigned int j = 0; j < dim_; ++j)
          dst[i][j][v] = src[i][j];
    }
}

template <int dim>
template <typename T1,
          int         rank_,
          std::size_t n_lanes,
          int         n_components_,
          int         dim_>
void
FERemoteEvaluationCommunicator<dim>::CopyInstructions::copy_data_entries(
  Tensor<rank_,
         n_components_,
         Tensor<rank_, dim_, VectorizedArray<T1, n_lanes>>>   &dst,
  const unsigned int                                           v,
  const Tensor<rank_, n_components_, Tensor<rank_, dim_, T1>> &src)
{
  if constexpr (rank_ == 1)
    {
      for (unsigned int i = 0; i < n_components_; ++i)
        copy_data(dst[i], v, src[i]);
    }
  else
    {
      for (unsigned int i = 0; i < rank_; ++i)
        for (unsigned int j = 0; j < n_components_; ++j)
          dst[i][j][v] = src[i][j];
    }
}

template <int dim>
template <typename T1, typename T2>
void
FERemoteEvaluationCommunicator<dim>::CopyInstructions::copy_data_entries(
  T1 &,
  const unsigned int,
  const T2 &)
{
  Assert(false,
         ExcMessage(
           "copy_data_entries() not implemented for given arguments."));
}



namespace Utilities
{
  template <int dim, typename Number, typename VectorizedArrayType>
  FERemoteEvaluationCommunicator<dim>
  compute_remote_communicator_faces_point_to_point_interpolation(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    const std::vector<
      std::pair<types::boundary_id, std::function<std::vector<bool>()>>>
                      &non_matching_faces_marked_vertices,
    const unsigned int quad_no,
    const unsigned int dof_no,
    const double       tolerance)
  {
    const auto &dof_handler = matrix_free.get_dof_handler(dof_no);
    const auto &tria        = dof_handler.get_triangulation();
    const auto &mapping     = *matrix_free.get_mapping_info().mapping;

    // Communication objects know about the communication pattern. I.e.,
    // they know about the cells and quadrature points that have to be
    // evaluated at remote faces. This information is given via
    // RemotePointEvaluation. Additionally, the communication objects
    // have to be able to match the quadrature points of the remote
    // points (that provide exterior information) to the quadrature points
    // defined at the interior cell. In case of point-to-point interpolation
    // a vector of pairs with face batch Ids and the number of faces in the
    // batch is needed. @c FERemoteCommunicationObjectEntityBatches
    // is a container to store this information.
    //
    // We need multiple communication objects (one for each non-matching face
    // ID).
    std::vector<FERemoteCommunicationObjectEntityBatches<dim>> comm_objects;

    // Additionally to the communication objects we need a vector
    // that stores quadrature rule sizes for every face batch.
    // The quadrature can have size zero in case of non non-matching faces,
    // i.e. boundary faces. Internally this information is needed to correctly
    // access values over multiple communication objects.
    std::vector<unsigned int> global_quadrature_sizes(
      matrix_free.n_boundary_face_batches(), numbers::invalid_unsigned_int);

    // Get the range of face batches we have to look at during construction of
    // the communication objects. We only have to look at boundary faces.
    const auto face_batch_range =
      std::make_pair(matrix_free.n_inner_face_batches(),
                     matrix_free.n_inner_face_batches() +
                       matrix_free.n_boundary_face_batches());

    // Iterate over all non-matching face IDs.
    for (const auto &[nm_face, marked_vertices] :
         non_matching_faces_marked_vertices)
      {
        // Construct the communication object for every face ID:
        // 1) RemotePointEvaluation with user specified function for marked
        // vertices.
        auto rpe = std::make_shared<Utilities::MPI::RemotePointEvaluation<dim>>(
          tolerance, false, 0, marked_vertices);

        // 2) Face batch IDs and number of faces in batch.
        std::vector<std::pair<unsigned int, unsigned int>>
          face_batch_id_n_faces;

        // Points that are searched by rpe.
        std::vector<Point<dim>> points;

        // Temporarily set up FEFaceEvaluation to access the quadrature points
        // at the faces on the non-matching interface.
        FEFaceEvaluation<dim, -1, 0, 1, Number> phi(matrix_free,
                                                    true,
                                                    dof_no,
                                                    quad_no);

        // Iterate over the boundary faces.
        for (unsigned int bface = 0;
             bface < face_batch_range.second - face_batch_range.first;
             ++bface)
          {
            const unsigned int face = face_batch_range.first + bface;

            if (matrix_free.get_boundary_id(face) == nm_face)
              {
                phi.reinit(face);

                // If @c face is on the current side of the non-matching
                // interface. Add the face batch ID and the number of faces in
                // the batch to the corresponding data structure.
                const unsigned int n_faces =
                  matrix_free.n_active_entries_per_face_batch(face);
                face_batch_id_n_faces.emplace_back(face, n_faces);

                // Append the quadrature points to the points we need to search
                // for.
                for (unsigned int v = 0; v < n_faces; ++v)
                  {
                    for (unsigned int q : phi.quadrature_point_indices())
                      {
                        const auto point = phi.quadrature_point(q);
                        Point<dim> temp;
                        for (unsigned int i = 0; i < dim; ++i)
                          temp[i] = point[i][v];

                        points.push_back(temp);
                      }
                  }

                // Insert the quadrature size into the global vector.
                // First check that each face is only considered once.
                Assert(global_quadrature_sizes[bface] ==
                         numbers::invalid_unsigned_int,
                       ExcMessage(
                         "Quadrature for given face already provided."));

                global_quadrature_sizes[bface] = phi.n_q_points;
              }
          }

        // Reinit RPE and ensure all points are found.
        rpe->reinit(points, tria, mapping);
        Assert(rpe->all_points_found(),
               ExcMessage("Not all remote points found."));

        // Add communication object to the list of objects.
        FERemoteCommunicationObjectEntityBatches<dim> co;
        co.batch_id_n_entities = face_batch_id_n_faces;
        co.rpe                 = rpe;
        comm_objects.push_back(co);
      }

    // Reinit the communicator `FERemoteEvaluationCommunicator`
    // with the communication objects.
    FERemoteEvaluationCommunicator<dim> remote_communicator;

    // if no quadrature size is set, an empty quadrature is considered
    std::replace(global_quadrature_sizes.begin(),
                 global_quadrature_sizes.end(),
                 numbers::invalid_unsigned_int,
                 0u);

    remote_communicator.reinit_faces(comm_objects,
                                     face_batch_range,
                                     global_quadrature_sizes);

    return remote_communicator;
  }



  template <int dim, typename Number, typename VectorizedArrayType>
  FERemoteEvaluationCommunicator<dim>
  compute_remote_communicator_faces_nitsche_type_mortaring(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    const std::vector<
      std::pair<types::boundary_id, std::function<std::vector<bool>()>>>
                      &non_matching_faces_marked_vertices,
    const unsigned int n_q_pnts_1D,
    const unsigned int dof_no,
    NonMatching::MappingInfo<dim, dim, Number> *nm_mapping_info,
    const double                                tolerance)
  {
    const auto &dof_handler = matrix_free.get_dof_handler(dof_no);
    const auto &tria        = dof_handler.get_triangulation();
    const auto &mapping     = *matrix_free.get_mapping_info().mapping;

    constexpr unsigned int n_lanes = VectorizedArray<Number>::size();

    std::pair<unsigned int, unsigned int> face_range =
      std::make_pair(matrix_free.n_inner_face_batches(),
                     matrix_free.n_inner_face_batches() +
                       matrix_free.n_boundary_face_batches());

    std::vector<Quadrature<dim - 1>> global_quadrature_vector(
      (matrix_free.n_inner_face_batches() +
       matrix_free.n_boundary_face_batches()) *
      n_lanes);

    // In case of Nitsche-type mortaring a vector of face indices is
    // needed as communication object.
    // @c FERemoteCommunicationObjectFaces is a container to store this
    // information.
    //
    // We need multiple communication objects (one for each non-matching face
    // ID).
    std::vector<FERemoteCommunicationObject<dim>> comm_objects;

    // Create bounding boxes and GridTools::Cache which is needed in
    // the following loop.
    std::vector<BoundingBox<dim>> local_boxes;
    for (const auto &cell : tria.active_cell_iterators())
      if (cell->is_locally_owned())
        local_boxes.emplace_back(mapping.get_bounding_box(cell));

    // Create r-tree of bounding boxes
    const auto local_tree = pack_rtree(local_boxes);

    // Compress r-tree to a minimal set of bounding boxes
    std::vector<std::vector<BoundingBox<dim>>> global_bboxes(1);
    global_bboxes[0] = extract_rtree_level(local_tree, 0);

    const GridTools::Cache<dim, dim> cache(tria, mapping);

    // Iterate over all sides of the non-matching interface.
    for (const auto &[nm_face, marked_vertices] :
         non_matching_faces_marked_vertices)
      {
        // 1) compute cell face pairs
        std::vector<
          std::pair<typename Triangulation<dim>::cell_iterator, unsigned int>>
          cell_face_pairs;

        std::vector<unsigned int> indices;

        for (unsigned int face = face_range.first; face < face_range.second;
             ++face)
          {
            if (matrix_free.get_boundary_id(face) == nm_face)
              {
                for (unsigned int v = 0;
                     v < matrix_free.n_active_entries_per_face_batch(face);
                     ++v)
                  {
                    const auto &[c, f] = matrix_free.get_face_iterator(face, v);

                    cell_face_pairs.emplace_back(std::make_pair(c, f));
                    indices.push_back(face * n_lanes + v);
                  }
              }
          }

        // 2) Create RPE.
        // In the Nitsche-type case we do not collect points for the setup
        // of RemotePointEvaluation. Instead we compute intersections between
        // the faces and set up RemotePointEvaluation with the computed
        // intersections.

        // Build intersection requests. Intersection requests
        // correspond to vertices at faces.
        std::vector<std::vector<Point<dim>>> intersection_requests;
        for (const auto &[cell, f] : cell_face_pairs)
          {
            std::vector<Point<dim>> vertices(cell->face(f)->n_vertices());
            std::copy_n(mapping.get_vertices(cell, f).begin(),
                        cell->face(f)->n_vertices(),
                        vertices.begin());
            intersection_requests.emplace_back(vertices);
          }

        // Compute intersection data with user specified function for marked
        // vertices.
        auto intersection_data =
          GridTools::internal::distributed_compute_intersection_locations<
            dim - 1>(cache,
                     intersection_requests,
                     global_bboxes,
                     marked_vertices(),
                     tolerance);

        // Convert to RPE.
        std::vector<Quadrature<dim>> mapped_quadratures_recv_comp;

        auto rpe =
          std::make_shared<Utilities::MPI::RemotePointEvaluation<dim>>();
        rpe->reinit(
          intersection_data
            .template convert_to_distributed_compute_point_locations_internal<
              dim>(n_q_pnts_1D, tria, mapping, &mapped_quadratures_recv_comp),
          tria,
          mapping);

        // 3) Fill global quadrature vector.
        for (unsigned int i = 0; i < intersection_requests.size(); ++i)
          {
            const auto idx = indices[i];

            // We do not use a structural binding here, since with
            // C++17 capturing structural bindings in lambdas leads
            // to an ill formed program.
            const auto &cell = std::get<0>(cell_face_pairs[i]);
            const auto &f    = std::get<1>(cell_face_pairs[i]);

            std::vector<Point<dim - 1>> q_points;
            std::vector<double>         weights;
            for (unsigned int ptr = intersection_data.recv_ptrs[i];
                 ptr < intersection_data.recv_ptrs[i + 1];
                 ++ptr)
              {
                const auto &quad = mapped_quadratures_recv_comp[ptr];

                const auto &ps = quad.get_points();
                std::transform(
                  ps.begin(),
                  ps.end(),
                  std::back_inserter(q_points),
                  [&](const Point<dim> &p) {
                    return mapping.project_real_point_to_unit_point_on_face(
                      cell, f, p);
                  });

                const auto &ws = quad.get_weights();
                weights.insert(weights.end(), ws.begin(), ws.end());
              }
            Quadrature<dim - 1> quad(q_points, weights);

            Assert(global_quadrature_vector[idx].size() == 0,
                   ExcMessage("Quadrature for given face already provided."));

            global_quadrature_vector[idx] = quad;
          }

        // Add communication object.
        FERemoteCommunicationObject<dim> co;
        co.indices = indices;
        co.rpe     = rpe;
        comm_objects.push_back(co);
      }

    // Reinit the communicator with the communication objects.
    FERemoteEvaluationCommunicator<dim> remote_communicator;

    std::vector<unsigned int> global_quadrature_sizes(
      global_quadrature_vector.size());
    std::transform(global_quadrature_vector.cbegin(),
                   global_quadrature_vector.cend(),
                   global_quadrature_sizes.begin(),
                   [](const auto &q) { return q.size(); });

    remote_communicator.reinit_faces(
      comm_objects,
      std::make_pair(0, global_quadrature_vector.size()),
      global_quadrature_sizes);

    if (nm_mapping_info != nullptr)
      {
        std::vector<
          std::pair<typename DoFHandler<dim>::cell_iterator, unsigned int>>
          vector_face_accessors;
        vector_face_accessors.reserve((matrix_free.n_inner_face_batches() +
                                       matrix_free.n_boundary_face_batches()) *
                                      n_lanes);

        // fill container for inner face batches
        unsigned int face_batch = 0;
        for (; face_batch < matrix_free.n_inner_face_batches(); ++face_batch)
          {
            for (unsigned int v = 0; v < n_lanes; ++v)
              {
                if (v < matrix_free.n_active_entries_per_face_batch(face_batch))
                  vector_face_accessors.push_back(
                    matrix_free.get_face_iterator(face_batch, v));
                else
                  vector_face_accessors.push_back(
                    matrix_free.get_face_iterator(face_batch, 0));
              }
          }
        // and boundary face batches
        for (; face_batch < (matrix_free.n_inner_face_batches() +
                             matrix_free.n_boundary_face_batches());
             ++face_batch)
          {
            for (unsigned int v = 0; v < n_lanes; ++v)
              {
                if (v < matrix_free.n_active_entries_per_face_batch(face_batch))
                  vector_face_accessors.push_back(
                    matrix_free.get_face_iterator(face_batch, v));
                else
                  vector_face_accessors.push_back(
                    matrix_free.get_face_iterator(face_batch, 0));
              }
          }

        nm_mapping_info->reinit_faces(vector_face_accessors,
                                      global_quadrature_vector);
      }

    return remote_communicator;
  }
} // namespace Utilities



template <int dim, int n_components, typename value_type>
template <typename MeshType>
FERemoteEvaluation<dim, n_components, value_type>::FERemoteEvaluation(
  const FERemoteEvaluationCommunicator<dim>          &comm,
  const MeshType                                     &mesh,
  const unsigned int                                  first_selected_component,
  const VectorTools::EvaluationFlags::EvaluationFlags evaluation_flags)
  : comm(&comm)
  , first_selected_component(first_selected_component)
  , evaluation_flags(evaluation_flags)
{
  set_mesh(mesh);
}

template <int dim, int n_components, typename value_type>
template <typename VectorType>
void
FERemoteEvaluation<dim, n_components, value_type>::gather_evaluate(
  const VectorType                      &src,
  const EvaluationFlags::EvaluationFlags flags)
{
  if (tria)
    {
      Assert(n_components == 1, ExcNotImplemented());
      comm->template update_ghost_values<n_components>(this->data,
                                                       *tria,
                                                       src,
                                                       flags,
                                                       first_selected_component,
                                                       evaluation_flags);
    }
  else if (dof_handler)
    {
      comm->template update_ghost_values<n_components>(this->data,
                                                       *dof_handler,
                                                       src,
                                                       flags,
                                                       first_selected_component,
                                                       evaluation_flags);
    }
  else
    DEAL_II_NOT_IMPLEMENTED();
}

template <int dim, int n_components, typename value_type>
internal::PrecomputedEvaluationDataAccessor<dim, n_components, value_type>
FERemoteEvaluation<dim, n_components, value_type>::get_data_accessor() const
{
  internal::PrecomputedEvaluationDataAccessor data_accessor(data,
                                                            comm->get_view());
  return data_accessor;
}

template <int dim, int n_components, typename value_type>
void
FERemoteEvaluation<dim, n_components, value_type>::set_mesh(
  const Triangulation<dim> &tria)
{
  this->tria = &tria;
}

template <int dim, int n_components, typename value_type>
void
FERemoteEvaluation<dim, n_components, value_type>::set_mesh(
  const DoFHandler<dim> &dof_handler)
{
  this->dof_handler = &dof_handler;
}

DEAL_II_NAMESPACE_CLOSE

#endif
