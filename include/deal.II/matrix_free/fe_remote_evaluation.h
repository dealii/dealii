// ---------------------------------------------------------------------
//
// Copyright (C) 2023 by the deal.II authors
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


#ifndef dealii_matrix_free_fe_remote_evaluation_h
#define dealii_matrix_free_fe_remote_evaluation_h

#include <deal.II/base/mpi_remote_point_evaluation.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/fe_point_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/vector_tools.h>


DEAL_II_NAMESPACE_OPEN

namespace internal
{
  /**
   * Check if given object is FEEvaluation object.
   */
  template <typename T>
  struct is_FEEvaluation : std::false_type
  {};

  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename Number,
            typename VectorizedArrayType>
  struct is_FEEvaluation<FEEvaluation<dim,
                                      fe_degree,
                                      n_q_points_1d,
                                      n_components,
                                      Number,
                                      VectorizedArrayType>> : std::true_type
  {};



  /**
   * Check if given object is FEFaceEvaluation object.
   */
  template <typename T>
  struct is_FEFaceEvaluation : std::false_type
  {};

  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components,
            typename Number,
            typename VectorizedArrayType>
  struct is_FEFaceEvaluation<FEFaceEvaluation<dim,
                                              fe_degree,
                                              n_q_points_1d,
                                              n_components,
                                              Number,
                                              VectorizedArrayType>>
    : std::true_type
  {};



  /**
   * Check if given object is FEPointEvaluation object.
   */
  template <typename T>
  struct is_FEPointEvaluation : std::false_type
  {};

  template <int n_components, int dim, int spacedim, typename Number>
  struct is_FEPointEvaluation<
    dealii::FEPointEvaluation<n_components, dim, spacedim, Number>>
    : std::true_type
  {};



  /**
   * Type traits for supported FEEvaluationTypes. Different FEEvaluationTypes
   * need different communication objects and different access to data at
   * quadrature points. Each specialization defines its @p CommunicationObjectType
   * and whether a two level CRS structure is needed to access the data by the
   * memeber `cell_face_pairs`. The same type with a different number of
   * components can be obtained with @p FEEvaluationTypeComponents.
   */
  template <typename FEEvaluationType, bool is_face>
  struct FEEvaluationTypeTraits
  {};



  /**
   * Specialization for FEEvaluation.
   */
  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components_,
            typename Number,
            typename VectorizedArrayType>
  struct FEEvaluationTypeTraits<FEEvaluation<dim,
                                             fe_degree,
                                             n_q_points_1d,
                                             n_components_,
                                             Number,
                                             VectorizedArrayType>,
                                false>
  {
    static constexpr bool cell_face_pairs = false;

    template <int n_components>
    using FEEvaluationTypeComponents = FEEvaluation<dim,
                                                    fe_degree,
                                                    n_q_points_1d,
                                                    n_components,
                                                    Number,
                                                    VectorizedArrayType>;

    using CommunicationObjectType = std::vector<
      std::pair<std::shared_ptr<Utilities::MPI::RemotePointEvaluation<dim>>,
                std::vector<std::pair<unsigned int, unsigned int>>>>;
  };



  /**
   * Specialization for FEFaceEvaluation.
   */
  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components_,
            typename Number,
            typename VectorizedArrayType>
  struct FEEvaluationTypeTraits<FEFaceEvaluation<dim,
                                                 fe_degree,
                                                 n_q_points_1d,
                                                 n_components_,
                                                 Number,
                                                 VectorizedArrayType>,
                                true>
  {
    static constexpr bool cell_face_pairs = false;

    template <int n_components>
    using FEEvaluationTypeComponents = FEEvaluation<dim,
                                                    fe_degree,
                                                    n_q_points_1d,
                                                    n_components,
                                                    Number,
                                                    VectorizedArrayType>;

    using CommunicationObjectType = std::vector<
      std::pair<std::shared_ptr<Utilities::MPI::RemotePointEvaluation<dim>>,
                std::vector<std::pair<unsigned int, unsigned int>>>>;
  };



  /**
   * Specialization for FEPointEvaluation.
   */
  template <int n_components_, int dim, int spacedim, typename Number>
  struct FEEvaluationTypeTraits<
    dealii::FEPointEvaluation<n_components_, dim, spacedim, Number>,
    false>
  {
    static constexpr bool cell_face_pairs = false;

    template <int n_components>
    using FEEvaluationTypeComponents =
      dealii::FEPointEvaluation<n_components, dim, spacedim, Number>;

    using CommunicationObjectType = std::vector<
      std::pair<std::shared_ptr<Utilities::MPI::RemotePointEvaluation<dim>>,
                std::vector<typename Triangulation<dim>::cell_iterator>>>;
  };



  /**
   * Specialization for FEPointEvaluation for faces.
   */
  template <int n_components_, int dim, int spacedim, typename Number>
  struct FEEvaluationTypeTraits<
    dealii::FEPointEvaluation<n_components_, dim, spacedim, Number>,
    true>
  {
    static constexpr bool cell_face_pairs = true;

    template <int n_components>
    using FEEvaluationTypeComponents =
      dealii::FEPointEvaluation<n_components, dim, spacedim, Number>;

    using CommunicationObjectType = std::vector<std::pair<
      std::shared_ptr<Utilities::MPI::RemotePointEvaluation<dim>>,
      std::vector<
        std::pair<typename Triangulation<dim>::cell_iterator, unsigned int>>>>;
  };



  /**
   * A class that stores values and/or gradients at quadrature points
   * corresponding to a FEEvaluationType (FEEvaluation, FEFaceEvaluation,
   * FEPointEvaluation).
   */
  template <typename FEEvaluationType>
  struct FERemoteEvaluationData
  {
    static constexpr unsigned int n_components = FEEvaluationType::n_components;

    using value_type    = typename FEEvaluationType::value_type;
    using gradient_type = typename FEEvaluationType::gradient_type;

    /**
     * values at quadrature points.
     */
    std::vector<value_type> values;

    /**
     * gradients at quadrature points.
     */
    std::vector<gradient_type> gradients;
  };



  /**
   * A class that stores a CRS like structure to access
   * FERemoteEvaluationData. If cell_face_pairs=false a CRS like structure
   * is created and the offset to the data can be obtained by
   * `get_shift(index)`. This case is used if quadrature points are only related
   * to a unique cell ID or face ID. If quadrature points are related to, e.g.,
   * a face on a given cell, cell_face_pairs=true and a two level CRS structure
   * is created. The offset to the data can be obtained by
   * `get_shift(cell_index, face_number)`.
   */
  template <bool cell_face_pairs>
  class FERemoteEvaluationDataView
  {};



  /**
   * Specialization for `cell_face_pairs=false`.
   */
  template <>
  class FERemoteEvaluationDataView<false>
  {
  public:
    /**
     * Get a pointer to data at index.
     */
    unsigned int
    get_shift(const unsigned int index) const
    {
      Assert(index != numbers::invalid_unsigned_int,
             ExcMessage("Index has to be valid!"));

      Assert(start <= index, ExcInternalError());
      AssertIndexRange(index - start, ptrs.size());
      return ptrs[index - start];
    }

    /**
     * Get the number of stored values.
     */
    unsigned int
    size() const
    {
      Assert(ptrs.size() > 0, ExcInternalError());
      return ptrs.back();
    }

    /**
     * Fill class from outside.
     */
    void
    reinit(const std::vector<unsigned int> &ptrs_in,
           const unsigned int               start_in = 0)
    {
      ptrs  = ptrs_in;
      start = start_in;
    }

  private:
    /**
     * This parameter can be used if indices do not start with 0.
     */
    unsigned int start = 0;
    /**
     * Pointers to data at index.
     */
    std::vector<unsigned int> ptrs;
  };



  /**
   * Specialization for `cell_face_pairs=true`.
   */
  template <>
  class FERemoteEvaluationDataView<true>
  {
  public:
    /**
     * Get a pointer to data at (cell_index, face_number).
     */
    unsigned int
    get_shift(const unsigned int cell_index,
              const unsigned int face_number) const
    {
      Assert(cell_index != numbers::invalid_unsigned_int,
             ExcMessage("Cell index has to be valid!"));
      Assert(face_number != numbers::invalid_unsigned_int,
             ExcMessage("Face number has to be valid!"));

      Assert(cell_start <= cell_index, ExcInternalError());

      AssertIndexRange(cell_index - cell_start, cell_ptrs.size());
      const unsigned int face_index =
        cell_ptrs[cell_index - cell_start] + face_number;
      AssertIndexRange(face_index, face_ptrs.size());
      return face_ptrs[face_index];
    }

    /**
     * Get the number of stored values.
     */
    unsigned int
    size() const
    {
      Assert(face_ptrs.size() > 0, ExcInternalError());
      return face_ptrs.back();
    }

    /**
     * Fill class from outside.
     */
    void
    reinit(const std::vector<unsigned int> &cell_ptrs_in,
           const std::vector<unsigned int> &face_ptrs_in,
           const unsigned int               cell_start_in = 0)
    {
      cell_ptrs  = cell_ptrs_in;
      face_ptrs  = face_ptrs_in;
      cell_start = cell_start_in;
    }

  private:
    /**
     * This parameter can be used if cell_indices do not start with 0.
     */
    unsigned int cell_start = 0;

    /**
     * Pointers to first face of given cell_index.
     */
    std::vector<unsigned int> cell_ptrs;

    /**
     * Pointers to data at (cell_index, face_number).
     */
    std::vector<unsigned int> face_ptrs;
  };

} // namespace internal



/**
 * A class to fill the fields in FERemoteEvaluationData.
 * FERemoteEvaluation is thought to be used with another @p FEEvaluationType
 * (FEEvaluation, FEFaceEvaluation, or FEPointEvaluation). @p is_face specifies
 * if @p FEEvaluationType works on faces.
 */
template <typename FEEvaluationType_, bool is_face_>
class FERemoteEvaluationCommunicator : public Subscriptor
{
  using FEETT =
    typename internal::FEEvaluationTypeTraits<FEEvaluationType_, is_face_>;

public:
  using FEEvaluationType = FEEvaluationType_;
  using FERemoteEvaluationDataViewType =
    internal::FERemoteEvaluationDataView<FEETT::cell_face_pairs>;
  static constexpr unsigned int dim     = FEEvaluationType_::dimension;
  static constexpr bool         is_face = is_face_;

  /**
   * Initialize communication patterns between (cell) quadrature points of
   * FEEvaluation and another triangulation.
   *
   * @param[in] this_fe_eval FEEvaluation object that describes the position of
   * the quadrature points.
   * @param[in] other_tria Remote triangulation object from which `values` and
   * `gradients` can be accessed in quadrature points of @p this_fe_eval.
   * @param[in] other_mapping Corresponding mapping for @p other_tria.
   * @param[in] tol Tolerance used to find quadrature points at other
   * triangulation.
   */
  template <typename T = FEEvaluationType_, bool F = is_face_>
  typename std::enable_if<true == internal::is_FEEvaluation<T>::value &&
                            false == F,
                          void>::type
  initialize_cells(FEEvaluationType         &this_fe_eval,
                   const Triangulation<dim> &other_tria,
                   const Mapping<dim>       &other_mapping,
                   const double              tol = 1e-9)
  {
    std::vector<unsigned int> cell_ptrs(
      this_fe_eval.get_matrix_free().n_cell_batches() + 1);
    cell_ptrs[0] = 0;

    std::vector<std::pair<unsigned int, unsigned int>> cells;
    std::vector<Point<dim>>                            points;
    Point<dim>                                         point;

    for (unsigned int cell = 0;
         cell < this_fe_eval.get_matrix_free().n_cell_batches();
         ++cell)
      {
        this_fe_eval.reinit(cell);
        cell_ptrs[cell + 1] = cell_ptrs[cell] + this_fe_eval.n_q_points;
        cells.push_back(
          cell,
          this_fe_eval.get_matrix_free().n_active_entries_per_cell_batch(cell));

        for (unsigned int v = 0;
             v < this_fe_eval.get_matrix_free().n_active_entries_per_cell_batch(
                   cell);
             ++v)
          {
            for (unsigned int q = 0; q < this_fe_eval.n_q_points; ++q)
              {
                for (unsigned int i = 0; i < dim; ++i)
                  point[i] = this_fe_eval.quadrature_point(q)[i][v];
                points.push_back(point);
              }
          }
      }

    auto rpe =
      std::make_shared<Utilities::MPI::RemotePointEvaluation<dim>>(tol,
                                                                   false,
                                                                   0);
    rpe->reinit(points, other_tria, other_mapping);

    communication_objects.push_back(std::make_pair(rpe, cells));

    view.reinit(cell_ptrs);
  }

  /**
   * Initialize communication patterns between (face) quadrature points of
   * FEFaceEvaluation and another triangulation.
   *
   * @param[in] faces_marked_vertices Boundary ID and marked vertices which are
   * possible connected. The ID corresponds to the faces which provide the
   * quadrature points. Quadrature points are searched in cells with at least
   * one marked vortex.
   * @param[in] this_fe_eval FEFaceEvaluation object that describes the position
   * of the quadrature points.
   * @param[in] other_tria Remote triangulation object from which `values` and
   * `gradients` can be accessed in quadrature points of @p this_fe_eval.
   * @param[in] other_mapping Corresponding mapping for @p other_tria.
   * @param[in] tol Tolerance used to find quadrature points at other
   * triangulation.
   */
  template <typename T = FEEvaluationType_, bool F = is_face_>
  typename std::enable_if<true == internal::is_FEFaceEvaluation<T>::value &&
                            true == F,
                          void>::type
  initialize_boundary_faces(
    const std::vector<
      std::pair<types::boundary_id, const std::function<std::vector<bool>()>>>
                             &faces_marked_vertices,
    FEEvaluationType         &this_fe_eval,
    const Triangulation<dim> &other_tria,
    const Mapping<dim>       &other_mapping,
    const double              tol = 1e-9)
  {
    const auto face_start =
      this_fe_eval.get_matrix_free().n_inner_face_batches();

    // data view and communication is setup in seperate steps because each face
    // pair needs its own communication object to be able to correctly mark
    // vertices for RPE.

    // setup data view
    {
      std::set<types::boundary_id> faces;
      for (const auto &face_marked_v : faces_marked_vertices)
        faces.insert(face_marked_v.first);

      std::vector<unsigned int> face_ptrs;
      face_ptrs.resize(
        this_fe_eval.get_matrix_free().n_boundary_face_batches() + 1);
      face_ptrs[0] = 0;

      for (unsigned int bface = 0;
           bface < this_fe_eval.get_matrix_free().n_boundary_face_batches();
           ++bface)
        {
          const unsigned int face = bface + face_start;

          if (faces.find(this_fe_eval.get_matrix_free().get_boundary_id(
                face)) != faces.end())
            {
              this_fe_eval.reinit(face);
              face_ptrs[bface + 1] = face_ptrs[bface] + this_fe_eval.n_q_points;
            }
          else
            {
              face_ptrs[bface + 1] = face_ptrs[bface];
            }
        }

      view.reinit(face_ptrs, face_start);
    }

    // setup communication
    {
      for (const auto &face_marked_v : faces_marked_vertices)
        {
          std::vector<std::pair<unsigned int, unsigned int>> faces;
          std::vector<Point<dim>>                            points;
          Point<dim>                                         point;

          for (unsigned int bface = 0;
               bface < this_fe_eval.get_matrix_free().n_boundary_face_batches();
               ++bface)
            {
              const unsigned int face = bface + face_start;

              if (this_fe_eval.get_matrix_free().get_boundary_id(face) ==
                  face_marked_v.first)
                {
                  this_fe_eval.reinit(face);
                  faces.emplace_back(face,
                                     this_fe_eval.get_matrix_free()
                                       .n_active_entries_per_face_batch(face));

                  for (unsigned int v = 0;
                       v < this_fe_eval.get_matrix_free()
                             .n_active_entries_per_face_batch(face);
                       ++v)
                    for (unsigned int q = 0; q < this_fe_eval.n_q_points; ++q)
                      {
                        for (unsigned int i = 0; i < dim; ++i)
                          point[i] = this_fe_eval.quadrature_point(q)[i][v];

                        points.push_back(point);
                      }
                }
            }


          auto rpe =
            std::make_shared<Utilities::MPI::RemotePointEvaluation<dim>>(
              tol, false, 0, face_marked_v.second);

          rpe->reinit(points, other_tria, other_mapping);

          Assert(rpe->all_points_found(),
                 ExcMessage("Not all remote points found."));

          communication_objects.push_back(std::make_pair(rpe, faces));
        }
    }
  }

  /**
   * Same as above, assuming the first boundary ID is connected to faces of the
   * second boundary ID.
   */
  template <typename T = FEEvaluationType_, bool F = is_face_>
  typename std::enable_if<true == internal::is_FEFaceEvaluation<T>::value &&
                            true == F,
                          void>::type
  initialize_face_pairs(
    const std::vector<std::pair<types::boundary_id, types::boundary_id>>
                             &face_pairs,
    FEEvaluationType         &this_fe_eval,
    const Triangulation<dim> &other_tria,
    const Mapping<dim>       &other_mapping,
    const double              tol = 1e-9)
  {
    std::vector<
      std::pair<types::boundary_id, const std::function<std::vector<bool>()>>>
      faces_marked_vertices;

    for (const auto &face_pair : face_pairs)
      {
        const auto id_dst    = face_pair.first;
        const auto id_remote = face_pair.second;
        faces_marked_vertices.push_back(std::make_pair(id_dst, [&]() {
          // only search points at cells on remote faces
          std::vector<bool> mask(other_tria.n_vertices(), false);

          for (const auto &face : other_tria.active_face_iterators())
            if (face->at_boundary() && face->boundary_id() == id_remote)
              for (const auto v : face->vertex_indices())
                mask[face->vertex_index(v)] = true;

          return mask;
        }));
      }

    initialize_boundary_faces(
      faces_marked_vertices, this_fe_eval, other_tria, other_mapping, tol);
  }

  /**
   * Same as above, assuming all face pairs are located in the same
   * triangulation.
   */
  template <typename T = FEEvaluationType_, bool F = is_face_>
  typename std::enable_if<true == internal::is_FEFaceEvaluation<T>::value &&
                            true == F,
                          void>::type
  initialize_face_pairs(
    const std::vector<std::pair<types::boundary_id, types::boundary_id>>
                     &face_pairs,
    FEEvaluationType &this_fe_eval,
    const double      tol = 1e-9)
  {
    initialize_face_pairs(
      face_pairs,
      this_fe_eval,
      this_fe_eval.get_matrix_free().get_dof_handler().get_triangulation(),
      *this_fe_eval.get_matrix_free().get_mapping_info().mapping,
      tol);
  }

  /**
   * Fill the fields stored in FERemoteEvaluationData.
   */
  template <typename FERemoteEvaluationDataType,
            typename MeshType,
            typename VectorType>
  void
  update_ghost_values(
    FERemoteEvaluationDataType                         &dst,
    const MeshType                                     &mesh,
    const VectorType                                   &src,
    const EvaluationFlags::EvaluationFlags              eval_flags,
    const VectorTools::EvaluationFlags::EvaluationFlags vec_flags,
    const unsigned int first_selected_component) const
  {
    constexpr unsigned int n_components =
      FERemoteEvaluationDataType::n_components;

    const bool has_ghost_elements = src.has_ghost_elements();

    if (has_ghost_elements == false)
      src.update_ghost_values();


    if (eval_flags & EvaluationFlags::values)
      {
        for (const auto &communication_object : communication_objects)
          copy_data(
            dst.values,
            VectorTools::point_values<n_components>(*communication_object.first,
                                                    mesh,
                                                    src,
                                                    vec_flags,
                                                    first_selected_component),
            communication_object.second);
      }

    if (eval_flags & EvaluationFlags::gradients)
      {
        for (const auto &communication_object : communication_objects)
          copy_data(dst.gradients,
                    VectorTools::point_gradients<n_components>(
                      *communication_object.first,
                      mesh,
                      src,
                      vec_flags,
                      first_selected_component),
                    communication_object.second);
      }

    Assert(!(eval_flags & EvaluationFlags::hessians), ExcNotImplemented());


    if (has_ghost_elements == false)
      src.zero_out_ghost_values();
  }

  /**
   * Get a pointer to data at index.
   */
  template <typename T = FEEvaluationType, bool F = is_face>
  typename std::enable_if<
    false == internal::FEEvaluationTypeTraits<T, F>::cell_face_pairs,
    unsigned int>::type
  get_shift(const unsigned int index) const
  {
    return view.get_shift(index);
  }

  /**
   * Get a pointer to data at (cell_index, face_number).
   */
  template <typename T = FEEvaluationType, bool F = is_face>
  typename std::enable_if<
    true == internal::FEEvaluationTypeTraits<T, F>::cell_face_pairs,
    unsigned int>::type
  get_shift(const unsigned int cell_index, const unsigned int face_number) const
  {
    return view.get_shift(cell_index, face_number);
  }

private:
  /**
   * Copy data obtained with RemotePointEvaluation to corresponding field
   * FERemoteEvaluationData.
   */
  template <typename T1, typename T2>
  void
  copy_data(std::vector<T1>       &dst,
            const std::vector<T2> &src,
            const std::vector<typename Triangulation<dim>::cell_iterator>
              &data_ptrs) const
  {
    dst.resize(view.size());

    if (data_ptrs.size() == 0)
      return;


    unsigned int c = 0;
    for (const auto &data_ptr : data_ptrs)
      {
        for (unsigned int j = get_shift(data_ptr.first->active_cell_index());
             j < get_shift(data_ptr->active_cell_index() + 1);
             ++j, ++c)
          {
            AssertIndexRange(j, dst.size());
            AssertIndexRange(c, src.size());

            dst[j] = src[c];
          }
      }
  }

  /**
   * Copy data obtained with RemotePointEvaluation to corresponding field
   * FERemoteEvaluationData.
   */
  template <typename T1, typename T2>
  void
  copy_data(
    std::vector<T1>                            &dst,
    const std::vector<T2>                      &src,
    const std::vector<std::pair<typename Triangulation<dim>::cell_iterator,
                                unsigned int>> &data_ptrs) const
  {
    dst.resize(view.size());

    if (data_ptrs.size() == 0)
      return;

    unsigned int c = 0;
    for (const auto &data_ptr : data_ptrs)
      {
        for (unsigned int j =
               get_shift(data_ptr.first->active_cell_index(), data_ptr.second);
             j < get_shift(data_ptr.first->active_cell_index(),
                           data_ptr.second + 1);
             ++j, ++c)
          {
            AssertIndexRange(j, dst.size());
            AssertIndexRange(c, src.size());

            dst[j] = src[c];
          }
      }
  }

  /**
   * Copy data obtained with RemotePointEvaluation to corresponding field
   * FERemoteEvaluationData.
   */
  template <typename T1, typename T2>
  void
  copy_data(
    std::vector<T1>                                          &dst,
    const std::vector<T2>                                    &src,
    const std::vector<std::pair<unsigned int, unsigned int>> &data_ptrs) const
  {
    dst.resize(view.size());

    unsigned int c = 0;
    for (const auto &data_ptr : data_ptrs)
      {
        const unsigned int bface     = data_ptr.first;
        const unsigned int n_entries = data_ptr.second;

        for (unsigned int v = 0; v < n_entries; ++v)
          for (unsigned int j = get_shift(bface); j < get_shift(bface + 1);
               ++j, ++c)
            {
              AssertIndexRange(j, dst.size());
              AssertIndexRange(c, src.size());

              copy_data(dst[j], v, src[c]);
            }
      }
  }

  /**
   * Copy data between different data layouts.
   */
  template <typename T1, std::size_t n_lanes>
  void
  copy_data(VectorizedArray<T1, n_lanes> &dst,
            const unsigned int            v,
            const T1                     &src) const
  {
    AssertIndexRange(v, n_lanes);

    dst[v] = src;
  }

  /**
   * Copy data between different data layouts.
   */
  template <typename T1, int rank_, std::size_t n_lanes, int dim_>
  void
  copy_data(Tensor<rank_, dim_, VectorizedArray<T1, n_lanes>> &dst,
            const unsigned int                                 v,
            const Tensor<rank_, dim_, T1>                     &src) const
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

  /**
   * Copy data between different data layouts.
   */
  template <typename T1,
            int         rank_,
            std::size_t n_lanes,
            int         n_components_,
            int         dim_>
  void
  copy_data(
    Tensor<rank_,
           n_components_,
           Tensor<rank_, dim_, VectorizedArray<T1, n_lanes>>>   &dst,
    const unsigned int                                           v,
    const Tensor<rank_, n_components_, Tensor<rank_, dim_, T1>> &src) const
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

  /**
   * CRS like data structure that describes the data positions at given
   * indices.
   */
  FERemoteEvaluationDataViewType view;
  /**
   * RemotePointEvaluation objects and indices to points used in
   * RemotePointEvaluation.
   */
  typename FEETT::CommunicationObjectType communication_objects;
};

/**
 * Class to access data in matrix-free loops for non-matching discretizations.
 * Interfaces are named with FEEvaluation, FEFaceEvaluation or FEPointEvaluation
 * in mind. The main difference is, that `gather_evaluate()` updates and caches
 * all values at once. Therefore, it has to be called only once before a
 * matrix-free loop.
 *
 * FERemoteEvaluation is thought to be used with another @p FEEvaluationType
 * (FEEvaluation, FEFaceEvaluation, or FEPointEvaluation).
 * FERemoteEvaluationCommunicator knows the type. However,
 * FERemoteEvaluationCommunicator is independent of @p n_components.
 */
template <typename FERemoteEvaluationCommunicatorType,
          unsigned int n_components_ = 1>
class FERemoteEvaluation
{
  static constexpr bool is_face = FERemoteEvaluationCommunicatorType::is_face;

  using FEETT = typename internal::FEEvaluationTypeTraits<
    typename FERemoteEvaluationCommunicatorType::FEEvaluationType,
    is_face>;

  using FEEvaluationType =
    typename FEETT::FEEvaluationTypeComponents<n_components_>;

  using FERemoteEvaluationDataType =
    internal::FERemoteEvaluationData<FEEvaluationType>;

public:
  static constexpr unsigned int n_components = n_components_;
  static constexpr unsigned int dimension    = FEEvaluationType::dimension;

  using Number        = typename FEEvaluationType::number_type;
  using value_type    = typename FERemoteEvaluationDataType::value_type;
  using gradient_type = typename FERemoteEvaluationDataType::gradient_type;

  /**
   * The constructor needs a corresponding FERemoteEvaluationCommunicator
   * which has to be setup outside of this class. This design choice is
   * motivated since the same FERemoteEvaluationCommunicator can be used
   * for different MeshTypes and number of components.
   *
   * @param[in] comm FERemoteEvaluationCommunicator.
   * @param[in] mesh Triangulation or DoFHandler.
   * @param[in] vt_flags Specify treatment of values at points which are found
   * on multiple cells.
   * @param[in] first_selected_component Select first component of evaluation in
   * DoFHandlers with multiple components.
   */
  template <typename MeshType>
  FERemoteEvaluation(const FERemoteEvaluationCommunicatorType &comm,
                     const MeshType                           &mesh,
                     const VectorTools::EvaluationFlags::EvaluationFlags
                       vt_flags = VectorTools::EvaluationFlags::avg,
                     const unsigned int first_selected_component = 0)
    : comm(&comm)
    , vt_flags(vt_flags)
    , first_selected_component(first_selected_component)
    , data_offset(numbers::invalid_unsigned_int)

  {
    set_mesh(mesh);
  }

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
                  const EvaluationFlags::EvaluationFlags flags)
  {
    if (tria)
      {
        AssertThrow(n_components == 1, ExcNotImplemented());
        comm->update_ghost_values(
          this->data, *tria, src, flags, vt_flags, first_selected_component);
      }
    else if (dof_handler)
      {
        comm->update_ghost_values(this->data,
                                  *dof_handler,
                                  src,
                                  flags,
                                  vt_flags,
                                  first_selected_component);
      }
    else
      AssertThrow(false, ExcNotImplemented());
  }

  /**
   * Set entity index at which quadrature points are accessed. This can, e.g.,
   * a cell index, a cell batch index, or a face batch index.
   */
  template <typename T = FERemoteEvaluationCommunicatorType>
  typename std::enable_if<
    std::is_same<typename T::FERemoteEvaluationDataViewType,
                 internal::FERemoteEvaluationDataView<false>>::value,
    void>::type
  reinit(const unsigned int index)
  {
    data_offset = comm->get_shift(index);
  }

  /**
   * Set cell and face_number at which quadrature points are accessed.
   */
  template <typename T = FERemoteEvaluationCommunicatorType>
  typename std::enable_if<
    std::is_same<typename T::FERemoteEvaluationDataViewType,
                 internal::FERemoteEvaluationDataView<true>>::value,
    void>::type
  reinit(const unsigned int cell_index, const unsigned int face_number)
  {
    data_offset = comm->get_shift(cell_index, face_number);
  }

  /**
   * Get the value at quadrature point @p q. The entity on which the values
   * are defined is set via `reinit()`.
   */
  const value_type
  get_value(const unsigned int q) const
  {
    Assert(data_offset != numbers::invalid_unsigned_int,
           ExcMessage("reinit() not called."));
    AssertIndexRange(data_offset + q, data.values.size());
    return data.values[data_offset + q];
  }

  /**
   * Get the gradients at quadrature pointt @p q. The entity on which the
   * gradients are defined is set via `reinit()`.
   */
  const gradient_type
  get_gradient(const unsigned int q) const
  {
    Assert(data_offset != numbers::invalid_unsigned_int,
           ExcMessage("reinit() not called."));
    AssertIndexRange(data_offset + q, data.gradients.size());
    return data.gradients[data_offset + q];
  }

private:
  /**
   * Use Triangulation as MeshType.
   */
  void
  set_mesh(const Triangulation<dimension> &tria)
  {
    this->tria = &tria;
  }

  /**
   * Use DoFHandler as MeshType.
   */
  void
  set_mesh(const DoFHandler<dimension> &dof_handler)
  {
    this->dof_handler = &dof_handler;
  }

  /**
   * Data that is accessed by `get_value()` and `get_gradient()`.
   */
  FERemoteEvaluationDataType data;

  /**
   * Underlying communicator which handles update of the ghost values and
   * gives position of values and gradients stored in
   * FERemoteEvaluationData.
   */
  SmartPointer<const FERemoteEvaluationCommunicatorType> comm;

  /**
   * Pointer to MeshType if used with Triangulation.
   */
  SmartPointer<const Triangulation<dimension>> tria;
  /**
   * Pointer to MeshType if used with DoFHandler.
   */
  SmartPointer<const DoFHandler<dimension>> dof_handler;

  /**
   * Flags that indicate which ghost values are updated.
   */
  const VectorTools::EvaluationFlags::EvaluationFlags vt_flags;

  /**
   * First selected component.
   */
  const unsigned int first_selected_component;

  /**
   * Offset to data after last call of `reinit()`.
   */
  unsigned int data_offset;
};


DEAL_II_NAMESPACE_CLOSE

#endif
