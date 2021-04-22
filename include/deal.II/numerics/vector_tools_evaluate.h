// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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


#ifndef dealii_vector_tools_evaluation_h
#define dealii_vector_tools_evaluation_h

#include <deal.II/base/config.h>

#include <deal.II/base/mpi_remote_point_evaluation.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_point_evaluation.h>

#include <map>

DEAL_II_NAMESPACE_OPEN

namespace VectorTools
{
  /**
   * Namespace for the flags for evaluate_at_points().
   */
  namespace EvaluationFlags
  {
    /**
     * Flags for evaluate_at_points().
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
   * Given a (distributed) solution vector @p vector, evaluate the values at
   * the (arbitrary and even remote) points specified by @p evaluation_points.
   */
  template <int n_components, int dim, int spacedim, typename VectorType>
  std::vector<typename FEPointEvaluation<n_components, dim>::value_type>
  evaluate_at_points(
    const Mapping<dim> &                                  mapping,
    const DoFHandler<dim, spacedim> &                     dof_handler,
    const VectorType &                                    vector,
    const std::vector<Point<spacedim>> &                  evaluation_points,
    Utilities::MPI::RemotePointEvaluation<dim, spacedim> &cache,
    const EvaluationFlags::EvaluationFlags flags = EvaluationFlags::avg);

  /**
   * Given a (distributed) solution vector @p vector, evaluate the values at
   * the points specified by @p cache which might have been set up by the
   * above function.
   *
   * @note Refinement/coarsening/repartitioning leads to the invalidation of the
   *   cache so that the above function has to be called again.
   */
  template <int n_components, int dim, int spacedim, typename VectorType>
  std::vector<typename FEPointEvaluation<n_components, dim>::value_type>
  evaluate_at_points(
    const Utilities::MPI::RemotePointEvaluation<dim, spacedim> &cache,
    const DoFHandler<dim, spacedim> &                           dof_handler,
    const VectorType &                                          vector,
    const EvaluationFlags::EvaluationFlags flags = EvaluationFlags::avg);



  // inlined functions


#ifndef DOXYGEN
  template <int n_components, int dim, int spacedim, typename VectorType>
  inline std::vector<typename FEPointEvaluation<n_components, dim>::value_type>
  evaluate_at_points(
    const Mapping<dim> &                                  mapping,
    const DoFHandler<dim, spacedim> &                     dof_handler,
    const VectorType &                                    vector,
    const std::vector<Point<spacedim>> &                  evaluation_points,
    Utilities::MPI::RemotePointEvaluation<dim, spacedim> &cache,
    const EvaluationFlags::EvaluationFlags                flags)
  {
    cache.reinit(evaluation_points, dof_handler.get_triangulation(), mapping);

    return evaluate_at_points<n_components>(cache, dof_handler, vector, flags);
  }



  namespace internal
  {
    /**
     * Perform reduction for scalars.
     */
    template <typename T>
    T
    reduce(const EvaluationFlags::EvaluationFlags &flags,
           const ArrayView<const T> &              values)
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
            Assert(false, ExcNotImplemented());
            return values[0];
        }
    }



    /**
     * Perform reduction for tensors.
     */
    template <int rank, int dim, typename Number>
    Tensor<rank, dim, Number>
    reduce(const EvaluationFlags::EvaluationFlags &          flags,
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
            Assert(false, ExcNotImplemented());
            return values[0];
        }
    }
  } // namespace internal



  template <int n_components, int dim, int spacedim, typename VectorType>
  inline std::vector<typename FEPointEvaluation<n_components, dim>::value_type>
  evaluate_at_points(
    const Utilities::MPI::RemotePointEvaluation<dim, spacedim> &cache,
    const DoFHandler<dim, spacedim> &                           dof_handler,
    const VectorType &                                          vector,
    const EvaluationFlags::EvaluationFlags                      flags)
  {
    using value_type =
      typename FEPointEvaluation<n_components, dim>::value_type;

    Assert(cache.is_ready(),
           ExcMessage(
             "Utilties::MPI::RemotePointEvaluation is not ready yet! "
             "Please call Utilties::MPI::RemotePointEvaluation::reinit() "
             "yourself or the other evaluate_at_points(), which does this for"
             "you."));

    Assert(&dof_handler.get_triangulation() == &cache.get_triangulation(),
           ExcMessage(
             "The provided Utilties::MPI::RemotePointEvaluation and DoFHandler "
             "object have been set up with different Triangulation objects, "
             "a scenario not supported!"));

    // evaluate values at points if possible
    const auto evaluation_point_results = [&]() {
      // helper function for accessing the global vector and interpolating
      // the results onto the points
      const auto evaluation_function = [&](auto &      values,
                                           const auto &cell_data) {
        std::vector<typename VectorType::value_type> solution_values;

        std::vector<std::unique_ptr<FEPointEvaluation<n_components, dim>>>
          evaluators(dof_handler.get_fe_collection().size());

        const auto get_evaluator = [&](const unsigned int active_fe_index)
          -> FEPointEvaluation<n_components, dim> & {
          if (evaluators[active_fe_index] == nullptr)
            evaluators[active_fe_index] =
              std::make_unique<FEPointEvaluation<n_components, dim>>(
                cache.get_mapping(), dof_handler.get_fe(active_fe_index));

          return *evaluators[active_fe_index];
        };

        for (unsigned int i = 0; i < cell_data.cells.size(); ++i)
          {
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

            auto &evaluator = get_evaluator(cell->active_fe_index());

            evaluator.evaluate(cell,
                               unit_points,
                               solution_values,
                               dealii::EvaluationFlags::values);

            for (unsigned int q = 0; q < unit_points.size(); ++q)
              values[q + cell_data.reference_point_ptrs[i]] =
                evaluator.get_value(q);
          }
      };

      std::vector<value_type> evaluation_point_results;
      std::vector<value_type> buffer;

      cache.template evaluate_and_process<value_type>(evaluation_point_results,
                                                      buffer,
                                                      evaluation_function);

      return evaluation_point_results;
    }();

    if (cache.is_map_unique())
      {
        // each point has exactly one result (unique map)
        return evaluation_point_results;
      }
    else
      {
        // map is not unique (multiple or no results): postprocessing is needed
        std::vector<value_type> unique_evaluation_point_results(
          cache.get_point_ptrs().size() - 1);

        const auto &ptr = cache.get_point_ptrs();

        for (unsigned int i = 0; i < ptr.size() - 1; ++i)
          {
            const auto n_entries = ptr[i + 1] - ptr[i];
            if (n_entries == 0)
              continue;

            unique_evaluation_point_results[i] =
              internal::reduce(flags,
                               ArrayView<const value_type>(
                                 evaluation_point_results.data() + ptr[i],
                                 n_entries));
          }

        return unique_evaluation_point_results;
      }
  }
#endif
} // namespace VectorTools

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_vector_tools_boundary_h
