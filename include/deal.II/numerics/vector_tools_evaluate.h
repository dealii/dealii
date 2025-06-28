// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_vector_tools_evaluation_h
#define dealii_vector_tools_evaluation_h

#include <deal.II/base/config.h>

#include <deal.II/base/mpi_remote_point_evaluation.h>

#include <deal.II/distributed/tria_base.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_point_evaluation.h>

DEAL_II_NAMESPACE_OPEN

namespace VectorTools
{
  /**
   * Namespace for the flags for point_values() and point_gradients().
   */
  namespace EvaluationFlags
  {
    /**
     * Flags for point_values() and point_gradients().
     */
    enum EvaluationFlags
    {
      /**
       * Compute average.
       */
      avg = 0,
      /**
       * Compute maximum.
       *
       * @note Only available for scalar values.
       */
      max = 1,
      /**
       * Compute minimum.
       *
       * @note Only available for scalar values.
       */
      min = 2,
      /**
       * Take any value.
       */
      insert = 3
    };
  } // namespace EvaluationFlags

  /**
   * Given a (distributed) solution vector @p vector, reinitialize @p cache and
   * evaluate the values at the (arbitrary and even remote) points specified by
   * @p evaluation_points.
   *
   * The following code snippet shows the usage of this function. Given
   * a Mapping object, a DoFHandler object, and solution vector as well as a
   * vector filled with points at which the vector should be evaluated, this
   * function returns a vector with values at those points. Furthermore, the
   * function initializes the communication pattern within the cache
   * Utilities::MPI::RemotePointEvaluation, which can be efficiently used in
   * further function calls (see also the function below).
   *
   * @code
   * Utilities::MPI::RemotePointEvaluation<dim, spacedim> cache;
   *
   * // Set up the cache and perform evaluation. This overload always
   * // reinitializes the RemotePointEvaluation cache object so, typically, it
   * // should only be called once.
   * const auto result_1 = VectorTools::point_values(
   *   mapping, dof_handler_1, vector_1, evaluation_points, cache);
   *
   * // This overload does not initialize or reinitialize the cache.
   * const auto result_2 = VectorTools::point_values(
   *   cache, dof_handler_2, vector_2);
   * @endcode
   *
   * Note that different DoFHandler objects can be passed to different
   * calls of this function. However, the underlying Triangulation object
   * needs to be the same if the cache should be reused.
   *
   * Alternatively, the user can set up the cache via
   * Utilities::MPI::RemotePointEvaluation::reinit() manually:
   *
   * @code
   * // set up cache manually
   * Utilities::MPI::RemotePointEvaluation<dim, spacedim> cache;
   * cache.reinit(evaluation_points, triangulation, mapping);
   *
   * // use the cache
   * const auto result_1 = VectorTools::point_values(
   *   cache, dof_handler_1, vector_1);
   *
   * const auto result_2 = VectorTools::point_values(
   *   cache, dof_handler_2, vector_2);
   * @endcode
   *
   * The function also works with FiniteElement objects with multiple
   * components. If one is interested only in a range of components, one can
   * select these by the parameters @p first_selected_component and
   * @p n_components. For further details on supported FiniteElement
   * objects, see the documentation of FEPointEvaluation.
   *
   * The function can also be used to evaluate cell-data vectors. For this
   * purpose, one passes in a Triangulation instead of a DoFHandler and a
   * vector of size Triangulation::n_active_cells() or a vector, which
   * has been initialized with the partitioner returned by
   * parallel::TriangulationBase::global_active_cell_index_partitioner().
   *
   * @note If a point cannot be found, the result for these points will
   *   be undefined (most probably 0). If you want to be sure that all
   *   points received a valid result, call `cache.all_points_found()`
   *   after this function call.
   *
   * @warning This is a collective call that needs to be executed by all
   *   processors in the communicator.
   *
   * @dealiiConceptRequires{(concepts::is_dealii_vector_type<VectorType> &&
   *    concepts::is_triangulation_or_dof_handler<MeshType<dim, spacedim>>)}
   */
  template <int n_components,
            template <int, int>
            class MeshType,
            int dim,
            int spacedim,
            typename VectorType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_vector_type<VectorType> &&
     concepts::is_triangulation_or_dof_handler<MeshType<dim, spacedim>>))
  std::vector<
    typename FEPointEvaluation<n_components,
                               dim,
                               spacedim,
                               typename VectorType::value_type>::
      value_type> point_values(const Mapping<dim>            &mapping,
                               const MeshType<dim, spacedim> &mesh,
                               const VectorType              &vector,
                               const std::vector<Point<spacedim>>
                                 &evaluation_points,
                               Utilities::MPI::RemotePointEvaluation<dim,
                                                                     spacedim>
                                                                     &cache,
                               const EvaluationFlags::EvaluationFlags flags =
                                 EvaluationFlags::avg,
                               const unsigned int first_selected_component = 0);

  /**
   * Given a (distributed) solution vector @p vector, evaluate the values at
   * the points specified by @p cache which might have been set up by the
   * above function.
   *
   * @note Refinement/coarsening/repartitioning leads to the invalidation of the
   *   cache so that the above function has to be called again. See also the
   *   example above.
   *
   * @warning This is a collective call that needs to be executed by all
   *   processors in the communicator.
   *
   * @dealiiConceptRequires{(concepts::is_dealii_vector_type<VectorType> &&
   *    concepts::is_triangulation_or_dof_handler<MeshType<dim, spacedim>>)}
   */
  template <int n_components,
            template <int, int>
            class MeshType,
            int dim,
            int spacedim,
            typename VectorType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_vector_type<VectorType> &&
     concepts::is_triangulation_or_dof_handler<MeshType<dim, spacedim>>))
  std::vector<
    typename FEPointEvaluation<n_components,
                               dim,
                               spacedim,
                               typename VectorType::value_type>::
      value_type> point_values(const Utilities::MPI::
                                 RemotePointEvaluation<dim, spacedim> &cache,
                               const MeshType<dim, spacedim>          &mesh,
                               const VectorType                       &vector,
                               const EvaluationFlags::EvaluationFlags  flags =
                                 EvaluationFlags::avg,
                               const unsigned int first_selected_component = 0);

  /**
   * Given a (distributed) solution vector @p vector, evaluate the gradients at
   * the (arbitrary and even remote) points specified by @p evaluation_points.
   *
   * @note The same comments as in the case of point_values() are true for this
   *   function.
   *
   * @warning This is a collective call that needs to be executed by all
   *   processors in the communicator.
   *
   * @dealiiConceptRequires{(concepts::is_dealii_vector_type<VectorType> &&
   *    concepts::is_triangulation_or_dof_handler<MeshType<dim, spacedim>>)}
   */
  template <int n_components,
            template <int, int>
            class MeshType,
            int dim,
            int spacedim,
            typename VectorType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_vector_type<VectorType> &&
     concepts::is_triangulation_or_dof_handler<MeshType<dim, spacedim>>))
  std::vector<
    typename FEPointEvaluation<n_components,
                               dim,
                               spacedim,
                               typename VectorType::value_type>::
      gradient_type> point_gradients(const Mapping<dim>            &mapping,
                                     const MeshType<dim, spacedim> &mesh,
                                     const VectorType              &vector,
                                     const std::vector<Point<spacedim>>
                                       &evaluation_points,
                                     Utilities::MPI::RemotePointEvaluation<
                                       dim,
                                       spacedim> &cache,
                                     const EvaluationFlags::EvaluationFlags
                                       flags = EvaluationFlags::avg,
                                     const unsigned int
                                       first_selected_component = 0);

  /**
   * Given a (distributed) solution vector @p vector, evaluate the gradients at
   * the points specified by @p cache which might have been set up by the
   * above function.
   *
   * @note The same comments as in the case of point_values() are true for this
   *   function.
   *
   * @warning This is a collective call that needs to be executed by all
   *   processors in the communicator.
   *
   * @dealiiConceptRequires{(concepts::is_dealii_vector_type<VectorType> &&
   *    concepts::is_triangulation_or_dof_handler<MeshType<dim, spacedim>>)}
   */
  template <int n_components,
            template <int, int>
            class MeshType,
            int dim,
            int spacedim,
            typename VectorType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_vector_type<VectorType> &&
     concepts::is_triangulation_or_dof_handler<MeshType<dim, spacedim>>))
  std::vector<
    typename FEPointEvaluation<n_components,
                               dim,
                               spacedim,
                               typename VectorType::value_type>::
      gradient_type> point_gradients(const Utilities::MPI::
                                       RemotePointEvaluation<dim, spacedim>
                                                                   &cache,
                                     const MeshType<dim, spacedim> &mesh,
                                     const VectorType              &vector,
                                     const EvaluationFlags::EvaluationFlags
                                       flags = EvaluationFlags::avg,
                                     const unsigned int
                                       first_selected_component = 0);



  // inlined functions


#ifndef DOXYGEN
  template <int n_components,
            template <int, int>
            class MeshType,
            int dim,
            int spacedim,
            typename VectorType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_vector_type<VectorType> &&
     concepts::is_triangulation_or_dof_handler<MeshType<dim, spacedim>>))
  inline std::vector<
    typename FEPointEvaluation<n_components,
                               dim,
                               spacedim,
                               typename VectorType::value_type>::
      value_type> point_values(const Mapping<dim>            &mapping,
                               const MeshType<dim, spacedim> &mesh,
                               const VectorType              &vector,
                               const std::vector<Point<spacedim>>
                                 &evaluation_points,
                               Utilities::MPI::RemotePointEvaluation<dim,
                                                                     spacedim>
                                                                     &cache,
                               const EvaluationFlags::EvaluationFlags flags,
                               const unsigned int first_selected_component)
  {
    cache.reinit(evaluation_points, mesh.get_triangulation(), mapping);

    return point_values<n_components>(
      cache, mesh, vector, flags, first_selected_component);
  }



  template <int n_components,
            template <int, int>
            class MeshType,
            int dim,
            int spacedim,
            typename VectorType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_vector_type<VectorType> &&
     concepts::is_triangulation_or_dof_handler<MeshType<dim, spacedim>>))
  inline std::vector<
    typename FEPointEvaluation<n_components,
                               dim,
                               spacedim,
                               typename VectorType::value_type>::
      gradient_type> point_gradients(const Mapping<dim>            &mapping,
                                     const MeshType<dim, spacedim> &mesh,
                                     const VectorType              &vector,
                                     const std::vector<Point<spacedim>>
                                       &evaluation_points,
                                     Utilities::MPI::RemotePointEvaluation<
                                       dim,
                                       spacedim> &cache,
                                     const EvaluationFlags::EvaluationFlags
                                       flags,
                                     const unsigned int
                                       first_selected_component)
  {
    cache.reinit(evaluation_points, mesh.get_triangulation(), mapping);

    return point_gradients<n_components>(
      cache, mesh, vector, flags, first_selected_component);
  }



  namespace internal
  {
    /**
     * Perform reduction for scalars.
     */
    template <typename T>
    T
    reduce(const EvaluationFlags::EvaluationFlags &flags,
           const ArrayView<const T>               &values)
    {
      switch (flags)
        {
          case EvaluationFlags::avg:
            {
              return std::accumulate(values.begin(), values.end(), T{}) /
                     (T(1.0) * values.size());
            }
          case EvaluationFlags::max:
            return *std::max_element(values.begin(), values.end());
          case EvaluationFlags::min:
            return *std::min_element(values.begin(), values.end());
          case EvaluationFlags::insert:
            return values[0];
          default:
            DEAL_II_NOT_IMPLEMENTED();
            return values[0];
        }
    }



    /**
     * Perform reduction for tensors.
     */
    template <int rank, int dim, typename Number>
    Tensor<rank, dim, Number>
    reduce(const EvaluationFlags::EvaluationFlags           &flags,
           const ArrayView<const Tensor<rank, dim, Number>> &values)
    {
      switch (flags)
        {
          case EvaluationFlags::avg:
            {
              return std::accumulate(values.begin(),
                                     values.end(),
                                     Tensor<rank, dim, Number>{}) /
                     (Number(1.0) * values.size());
            }
          case EvaluationFlags::insert:
            return values[0];
          default:
            DEAL_II_NOT_IMPLEMENTED();
            return values[0];
        }
    }



    /**
     * Perform reduction for tensors of tensors (e.g., gradient of
     * vectorial quantities).
     */
    template <int n_components, int rank, int dim, typename Number>
    Tensor<1, n_components, Tensor<rank, dim, Number>>
    reduce(
      const EvaluationFlags::EvaluationFlags &flags,
      const ArrayView<const Tensor<1, n_components, Tensor<rank, dim, Number>>>
        &values)
    {
      switch (flags)
        {
          case EvaluationFlags::avg:
            {
              Tensor<1, n_components, Tensor<rank, dim, Number>> temp;

              for (unsigned int j = 0; j < values.size(); ++j)
                for (unsigned int i = 0; i < n_components; ++i)
                  temp[i] = temp[i] + values[j][i];

              for (unsigned int i = 0; i < n_components; ++i)
                temp[i] /= Number(values.size());

              return temp;
            }
          case EvaluationFlags::insert:
            return values[0];
          default:
            DEAL_II_NOT_IMPLEMENTED();
            return values[0];
        }
    }



    template <int n_components,
              int dim,
              int spacedim,
              typename VectorType,
              typename value_type>
    void
    process_cell(
      const unsigned int i,
      const typename Utilities::MPI::RemotePointEvaluation<dim,
                                                           spacedim>::CellData
                                                                 &cell_data,
      const Utilities::MPI::RemotePointEvaluation<dim, spacedim> &cache,
      const DoFHandler<dim, spacedim>                            &dof_handler,
      const VectorType                                           &vector,
      const UpdateFlags                                           update_flags,
      const dealii::EvaluationFlags::EvaluationFlags evaluation_flags,
      const unsigned int                             first_selected_component,
      const std::function<
        value_type(const FEPointEvaluation<n_components,
                                           dim,
                                           spacedim,
                                           typename VectorType::value_type> &,
                   const unsigned int &)>           process_quadrature_point,
      const ArrayView<value_type>                  &values,
      std::vector<typename VectorType::value_type> &solution_values,
      std::vector<
        std::unique_ptr<FEPointEvaluation<n_components,
                                          dim,
                                          spacedim,
                                          typename VectorType::value_type>>>
        &evaluators)
    {
      if (evaluators.empty())
        evaluators.resize(dof_handler.get_fe_collection().size());

      typename DoFHandler<dim>::active_cell_iterator cell = {
        &cache.get_triangulation(),
        cell_data.cells[i].first,
        cell_data.cells[i].second,
        &dof_handler};

      const ArrayView<const Point<dim>> unit_points(
        cell_data.reference_point_values.data() +
          cell_data.reference_point_ptrs[i],
        cell_data.reference_point_ptrs[i + 1] -
          cell_data.reference_point_ptrs[i]);

      solution_values.resize(
        dof_handler.get_fe(cell->active_fe_index()).n_dofs_per_cell());
      cell->get_dof_values(vector,
                           solution_values.begin(),
                           solution_values.end());

      if (evaluators[cell->active_fe_index()] == nullptr)
        evaluators[cell->active_fe_index()] =
          std::make_unique<FEPointEvaluation<n_components,
                                             dim,
                                             spacedim,
                                             typename VectorType::value_type>>(
            cache.get_mapping(),
            cell->get_fe(),
            update_flags,
            first_selected_component);
      auto &evaluator = *evaluators[cell->active_fe_index()];

      evaluator.reinit(cell, unit_points);
      evaluator.evaluate(solution_values, evaluation_flags);

      for (unsigned int q = 0; q < unit_points.size(); ++q)
        values[q + cell_data.reference_point_ptrs[i]] =
          process_quadrature_point(evaluator, q);
    }



    template <int dim, int spacedim, typename Number>
    Number
    get_value(
      const Triangulation<dim, spacedim>                                &tria,
      const Vector<Number>                                              &vector,
      const typename Triangulation<dim, spacedim>::active_cell_iterator &cell)
    {
      AssertDimension(tria.n_active_cells(), vector.size());
      return vector[cell->active_cell_index()];
    }



    template <int dim, int spacedim, typename Number>
    Number
    get_value(
      const Triangulation<dim, spacedim>                                &tria,
      const LinearAlgebra::distributed::Vector<Number>                  &vector,
      const typename Triangulation<dim, spacedim>::active_cell_iterator &cell)
    {
      const auto distributed_tria =
        dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(&tria);

      const bool use_distributed_path =
        (distributed_tria == nullptr) ?
          false :
          (vector.get_partitioner().get() ==
           distributed_tria->global_active_cell_index_partitioner()
             .lock()
             .get());

      if (use_distributed_path)
        {
          return vector[cell->global_active_cell_index()];
        }
      else
        {
          AssertDimension(tria.n_active_cells(), vector.locally_owned_size());
          return vector[cell->active_cell_index()];
        }
    }



    template <typename Number, typename Number2>
    void
    set_value(Number &dst, const Number2 &src)
    {
      dst = src;
    }



    template <typename Number, int rank, int dim, typename Number2>
    void
    set_value(Tensor<rank, dim, Number> &, const Number2 &)
    {
      Assert(false,
             ExcMessage(
               "A cell-data vector can only have a single component."));
    }



    template <int n_components,
              int dim,
              int spacedim,
              typename VectorType,
              typename value_type>
    void
    process_cell(
      const unsigned int i,
      const typename Utilities::MPI::RemotePointEvaluation<dim,
                                                           spacedim>::CellData
        &cell_data,
      const Utilities::MPI::RemotePointEvaluation<dim, spacedim> &,
      const Triangulation<dim, spacedim> &triangulation,
      const VectorType                   &vector,
      const UpdateFlags,
      const dealii::EvaluationFlags::EvaluationFlags evaluation_flags,
      const unsigned int                             first_selected_component,
      const std::function<
        value_type(const FEPointEvaluation<n_components,
                                           dim,
                                           spacedim,
                                           typename VectorType::value_type> &,
                   const unsigned int &)>,
      const ArrayView<value_type> &values,
      std::vector<typename VectorType::value_type> &,
      std::vector<
        std::unique_ptr<FEPointEvaluation<n_components,
                                          dim,
                                          spacedim,
                                          typename VectorType::value_type>>> &)
    {
      Assert(n_components == 1 && first_selected_component == 0,
             ExcMessage(
               "A cell-data vector can only have a single component."));

      Assert(evaluation_flags ==
               dealii::EvaluationFlags::EvaluationFlags::values,
             ExcMessage("For cell-data vectors, only values can be queried."));

      typename Triangulation<dim>::active_cell_iterator cell = {
        &triangulation, cell_data.cells[i].first, cell_data.cells[i].second};

      const auto value = get_value(triangulation, vector, cell);

      for (unsigned int q = cell_data.reference_point_ptrs[i];
           q < cell_data.reference_point_ptrs[i + 1];
           ++q)
        set_value(values[q], value);
    }



    template <int n_components,
              int dim,
              int spacedim,
              typename MeshType,
              typename VectorType,
              typename value_type>
    DEAL_II_CXX20_REQUIRES(
      concepts::is_dealii_vector_type<VectorType>
        &&concepts::is_triangulation_or_dof_handler<MeshType>)
    inline std::vector<value_type> evaluate_at_points(
      const Utilities::MPI::RemotePointEvaluation<dim, spacedim> &cache,
      const MeshType                                             &mesh,
      const VectorType                                           &vector,
      const EvaluationFlags::EvaluationFlags                      flags,
      const unsigned int                             first_selected_component,
      const UpdateFlags                              update_flags,
      const dealii::EvaluationFlags::EvaluationFlags evaluation_flags,
      const std::function<
        value_type(const FEPointEvaluation<n_components,
                                           dim,
                                           spacedim,
                                           typename VectorType::value_type> &,
                   const unsigned int &)> process_quadrature_point)
    {
      Assert(cache.is_ready(),
             ExcMessage(
               "Utilities::MPI::RemotePointEvaluation is not ready yet! "
               "Please call Utilities::MPI::RemotePointEvaluation::reinit() "
               "yourself or another function that does this for you."));

      Assert(
        &mesh.get_triangulation() == &cache.get_triangulation(),
        ExcMessage(
          "The provided Utilities::MPI::RemotePointEvaluation and DoFHandler "
          "object have been set up with different Triangulation objects, "
          "a scenario not supported!"));

      // evaluate values at points if possible
      const auto evaluation_point_results = [&]() {
        // helper function for accessing the global vector and interpolating
        // the results onto the points
        const auto evaluation_function = [&](auto       &values,
                                             const auto &cell_data) {
          std::vector<typename VectorType::value_type> solution_values;

          std::vector<
            std::unique_ptr<FEPointEvaluation<n_components,
                                              dim,
                                              spacedim,
                                              typename VectorType::value_type>>>
            evaluators;

          for (unsigned int i = 0; i < cell_data.cells.size(); ++i)
            process_cell<n_components, dim, spacedim, VectorType, value_type>(
              i,
              cell_data,
              cache,
              mesh,
              vector,
              update_flags,
              evaluation_flags,
              first_selected_component,
              process_quadrature_point,
              values,
              solution_values,
              evaluators);
        };

        std::vector<value_type> evaluation_point_results;
        std::vector<value_type> buffer;

        cache.template evaluate_and_process<value_type>(
          evaluation_point_results, buffer, evaluation_function);

        return evaluation_point_results;
      }();

      if (cache.is_map_unique())
        {
          // each point has exactly one result (unique map)
          return evaluation_point_results;
        }
      else
        {
          // map is not unique (multiple or no results): postprocessing is
          // needed
          std::vector<value_type> unique_evaluation_point_results(
            cache.get_point_ptrs().size() - 1);

          const auto &ptr = cache.get_point_ptrs();

          for (unsigned int i = 0; i < ptr.size() - 1; ++i)
            {
              const auto n_entries = ptr[i + 1] - ptr[i];
              if (n_entries == 0)
                continue;

              unique_evaluation_point_results[i] =
                reduce(flags,
                       ArrayView<const value_type>(
                         evaluation_point_results.data() + ptr[i], n_entries));
            }

          return unique_evaluation_point_results;
        }
    }
  } // namespace internal

  template <int n_components,
            template <int, int>
            class MeshType,
            int dim,
            int spacedim,
            typename VectorType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_vector_type<VectorType> &&
     concepts::is_triangulation_or_dof_handler<MeshType<dim, spacedim>>))
  inline std::vector<
    typename FEPointEvaluation<n_components,
                               dim,
                               spacedim,
                               typename VectorType::value_type>::
      value_type> point_values(const Utilities::MPI::
                                 RemotePointEvaluation<dim, spacedim> &cache,
                               const MeshType<dim, spacedim>          &mesh,
                               const VectorType                       &vector,
                               const EvaluationFlags::EvaluationFlags  flags,
                               const unsigned int first_selected_component)
  {
    return internal::evaluate_at_points<
      n_components,
      dim,
      spacedim,
      MeshType<dim, spacedim>,
      VectorType,
      typename FEPointEvaluation<n_components,
                                 dim,
                                 spacedim,
                                 typename VectorType::value_type>::value_type>(
      cache,
      mesh,
      vector,
      flags,
      first_selected_component,
      update_values,
      dealii::EvaluationFlags::values,
      [](const auto &evaluator, const auto &q) {
        return evaluator.get_value(q);
      });
  }



  template <int n_components,
            template <int, int>
            class MeshType,
            int dim,
            int spacedim,
            typename VectorType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_dealii_vector_type<VectorType> &&
     concepts::is_triangulation_or_dof_handler<MeshType<dim, spacedim>>))
  inline std::vector<
    typename FEPointEvaluation<n_components,
                               dim,
                               spacedim,
                               typename VectorType::value_type>::
      gradient_type> point_gradients(const Utilities::MPI::
                                       RemotePointEvaluation<dim, spacedim>
                                                                   &cache,
                                     const MeshType<dim, spacedim> &mesh,
                                     const VectorType              &vector,
                                     const EvaluationFlags::EvaluationFlags
                                       flags,
                                     const unsigned int
                                       first_selected_component)
  {
    return internal::evaluate_at_points<
      n_components,
      dim,
      spacedim,
      MeshType<dim, spacedim>,
      VectorType,
      typename FEPointEvaluation<
        n_components,
        dim,
        spacedim,
        typename VectorType::value_type>::gradient_type>(
      cache,
      mesh,
      vector,
      flags,
      first_selected_component,
      update_gradients,
      dealii::EvaluationFlags::gradients,
      [](const auto &evaluator, const unsigned &q) {
        return evaluator.get_gradient(q);
      });
  }

#endif
} // namespace VectorTools

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_vector_tools_boundary_h
