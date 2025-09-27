// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/config.h>

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/tria_base.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/hp/refinement.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/vector.h>

#include <limits>

DEAL_II_NAMESPACE_OPEN

namespace hp
{
  namespace Refinement
  {
    /**
     * Setting p-adaptivity flags
     */
    template <int dim, int spacedim>
    void
    full_p_adaptivity(const DoFHandler<dim, spacedim> &dof_handler)
    {
      if (dof_handler.get_fe_collection().empty())
        // nothing to do
        return;

      Assert(dof_handler.has_hp_capabilities(),
             (typename DoFHandler<dim, spacedim>::ExcOnlyAvailableWithHP()));

      std::vector<bool> p_flags(
        dof_handler.get_triangulation().n_active_cells(), true);

      p_adaptivity_from_flags(dof_handler, p_flags);
    }



    template <int dim, int spacedim>
    void
    p_adaptivity_from_flags(const DoFHandler<dim, spacedim> &dof_handler,
                            const std::vector<bool>         &p_flags)
    {
      if (dof_handler.get_fe_collection().empty())
        // nothing to do
        return;

      Assert(dof_handler.has_hp_capabilities(),
             (typename DoFHandler<dim, spacedim>::ExcOnlyAvailableWithHP()));
      AssertDimension(dof_handler.get_triangulation().n_active_cells(),
                      p_flags.size());

      for (const auto &cell : dof_handler.active_cell_iterators())
        if (cell->is_locally_owned() && p_flags[cell->active_cell_index()])
          {
            if (cell->refine_flag_set())
              {
                const unsigned int super_fe_index =
                  dof_handler.get_fe_collection().next_in_hierarchy(
                    cell->active_fe_index());

                // Reject update if already most superordinate element.
                if (super_fe_index != cell->active_fe_index())
                  cell->set_future_fe_index(super_fe_index);
              }
            else if (cell->coarsen_flag_set())
              {
                const unsigned int sub_fe_index =
                  dof_handler.get_fe_collection().previous_in_hierarchy(
                    cell->active_fe_index());

                // Reject update if already least subordinate element.
                if (sub_fe_index != cell->active_fe_index())
                  cell->set_future_fe_index(sub_fe_index);
              }
          }
    }



    template <int dim, typename Number, int spacedim>
    void
    p_adaptivity_from_absolute_threshold(
      const DoFHandler<dim, spacedim> &dof_handler,
      const Vector<Number>            &criteria,
      const Number                     p_refine_threshold,
      const Number                     p_coarsen_threshold,
      const ComparisonFunction<std_cxx20::type_identity_t<Number>>
        &compare_refine,
      const ComparisonFunction<std_cxx20::type_identity_t<Number>>
        &compare_coarsen)
    {
      if (dof_handler.get_fe_collection().empty())
        // nothing to do
        return;

      Assert(dof_handler.has_hp_capabilities(),
             (typename DoFHandler<dim, spacedim>::ExcOnlyAvailableWithHP()));
      AssertDimension(dof_handler.get_triangulation().n_active_cells(),
                      criteria.size());

      std::vector<bool> p_flags(
        dof_handler.get_triangulation().n_active_cells(), false);

      for (const auto &cell : dof_handler.active_cell_iterators())
        if (cell->is_locally_owned() &&
            ((cell->refine_flag_set() &&
              compare_refine(criteria[cell->active_cell_index()],
                             p_refine_threshold)) ||
             (cell->coarsen_flag_set() &&
              compare_coarsen(criteria[cell->active_cell_index()],
                              p_coarsen_threshold))))
          p_flags[cell->active_cell_index()] = true;

      p_adaptivity_from_flags(dof_handler, p_flags);
    }



    template <int dim, typename Number, int spacedim>
    void
    p_adaptivity_from_relative_threshold(
      const DoFHandler<dim, spacedim> &dof_handler,
      const Vector<Number>            &criteria,
      const double                     p_refine_fraction,
      const double                     p_coarsen_fraction,
      const ComparisonFunction<std_cxx20::type_identity_t<Number>>
        &compare_refine,
      const ComparisonFunction<std_cxx20::type_identity_t<Number>>
        &compare_coarsen)
    {
      if (dof_handler.get_fe_collection().empty())
        // nothing to do
        return;

      Assert(dof_handler.has_hp_capabilities(),
             (typename DoFHandler<dim, spacedim>::ExcOnlyAvailableWithHP()));
      AssertDimension(dof_handler.get_triangulation().n_active_cells(),
                      criteria.size());
      Assert((p_refine_fraction >= 0) && (p_refine_fraction <= 1),
             GridRefinement::ExcInvalidParameterValue());
      Assert((p_coarsen_fraction >= 0) && (p_coarsen_fraction <= 1),
             GridRefinement::ExcInvalidParameterValue());

      // We first have to determine the maximal and minimal values of the
      // criteria of all flagged cells.
      Number max_criterion_refine  = std::numeric_limits<Number>::lowest(),
             min_criterion_refine  = std::numeric_limits<Number>::max();
      Number max_criterion_coarsen = max_criterion_refine,
             min_criterion_coarsen = min_criterion_refine;

      for (const auto &cell : dof_handler.active_cell_iterators())
        if (cell->is_locally_owned())
          {
            if (cell->refine_flag_set())
              {
                max_criterion_refine =
                  std::max(max_criterion_refine,
                           criteria(cell->active_cell_index()));
                min_criterion_refine =
                  std::min(min_criterion_refine,
                           criteria(cell->active_cell_index()));
              }
            else if (cell->coarsen_flag_set())
              {
                max_criterion_coarsen =
                  std::max(max_criterion_coarsen,
                           criteria(cell->active_cell_index()));
                min_criterion_coarsen =
                  std::min(min_criterion_coarsen,
                           criteria(cell->active_cell_index()));
              }
          }

      const parallel::TriangulationBase<dim, spacedim> *parallel_tria =
        dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
          &dof_handler.get_triangulation());
      if (parallel_tria != nullptr &&
          dynamic_cast<const parallel::shared::Triangulation<dim, spacedim> *>(
            &dof_handler.get_triangulation()) == nullptr)
        {
          max_criterion_refine =
            Utilities::MPI::max(max_criterion_refine,
                                parallel_tria->get_mpi_communicator());
          min_criterion_refine =
            Utilities::MPI::min(min_criterion_refine,
                                parallel_tria->get_mpi_communicator());
          max_criterion_coarsen =
            Utilities::MPI::max(max_criterion_coarsen,
                                parallel_tria->get_mpi_communicator());
          min_criterion_coarsen =
            Utilities::MPI::min(min_criterion_coarsen,
                                parallel_tria->get_mpi_communicator());
        }

      // Absent any better strategies, we will set the threshold by linear
      // interpolation for both classes of cells individually.
      const Number threshold_refine =
                     min_criterion_refine +
                     p_refine_fraction *
                       (max_criterion_refine - min_criterion_refine),
                   threshold_coarsen =
                     min_criterion_coarsen +
                     p_coarsen_fraction *
                       (max_criterion_coarsen - min_criterion_coarsen);

      p_adaptivity_from_absolute_threshold(dof_handler,
                                           criteria,
                                           threshold_refine,
                                           threshold_coarsen,
                                           compare_refine,
                                           compare_coarsen);
    }



    template <int dim, typename Number, int spacedim>
    void
    p_adaptivity_fixed_number(
      const DoFHandler<dim, spacedim> &dof_handler,
      const Vector<Number>            &criteria,
      const double                     p_refine_fraction,
      const double                     p_coarsen_fraction,
      const ComparisonFunction<std_cxx20::type_identity_t<Number>>
        &compare_refine,
      const ComparisonFunction<std_cxx20::type_identity_t<Number>>
        &compare_coarsen)
    {
      if (dof_handler.get_fe_collection().empty())
        // nothing to do
        return;

      Assert(dof_handler.has_hp_capabilities(),
             (typename DoFHandler<dim, spacedim>::ExcOnlyAvailableWithHP()));
      AssertDimension(dof_handler.get_triangulation().n_active_cells(),
                      criteria.size());
      Assert((p_refine_fraction >= 0) && (p_refine_fraction <= 1),
             GridRefinement::ExcInvalidParameterValue());
      Assert((p_coarsen_fraction >= 0) && (p_coarsen_fraction <= 1),
             GridRefinement::ExcInvalidParameterValue());

      // ComparisonFunction returning 'true' or 'false' for any set of
      // parameters. These will be used to overwrite user-provided comparison
      // functions whenever no actual comparison is required in the decision
      // process, i.e. when no or all cells will be refined or coarsened.
      const ComparisonFunction<Number> compare_false =
        [](const Number &, const Number &) { return false; };
      const ComparisonFunction<Number> compare_true =
        [](const Number &, const Number &) { return true; };

      // 1.) First extract from the vector of indicators the ones that
      //     correspond to cells that we locally own.
      unsigned int   n_flags_refinement = 0;
      unsigned int   n_flags_coarsening = 0;
      Vector<Number> indicators_refinement(
        dof_handler.get_triangulation().n_active_cells());
      Vector<Number> indicators_coarsening(
        dof_handler.get_triangulation().n_active_cells());
      for (const auto &cell :
           dof_handler.get_triangulation().active_cell_iterators())
        if (!cell->is_artificial() && cell->is_locally_owned())
          {
            if (cell->refine_flag_set())
              indicators_refinement(n_flags_refinement++) =
                criteria(cell->active_cell_index());
            else if (cell->coarsen_flag_set())
              indicators_coarsening(n_flags_coarsening++) =
                criteria(cell->active_cell_index());
          }
      indicators_refinement.grow_or_shrink(n_flags_refinement);
      indicators_coarsening.grow_or_shrink(n_flags_coarsening);

      // 2.) Determine the number of cells for p-refinement and p-coarsening on
      //     basis of the flagged cells.
      //
      // 3.) Find thresholds for p-refinement and p-coarsening on only those
      //     cells flagged for adaptation.
      //
      //     For cases in which no or all cells flagged for refinement and/or
      //     coarsening are subject to p-adaptation, we usually pick thresholds
      //     that apply to all or none of the cells at once. However here, we
      //     do not know which threshold would suffice for this task because the
      //     user could provide any comparison function. Thus if necessary, we
      //     overwrite the user's choice with suitable functions simply
      //     returning 'true' and 'false' for any cell with reference wrappers.
      //     Thus, no function object copies are stored.
      //
      // 4.) Perform p-adaptation with absolute thresholds.
      Number threshold_refinement      = 0.;
      Number threshold_coarsening      = 0.;
      auto   reference_compare_refine  = std::cref(compare_refine);
      auto   reference_compare_coarsen = std::cref(compare_coarsen);

      const parallel::TriangulationBase<dim, spacedim> *parallel_tria =
        dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
          &dof_handler.get_triangulation());
      if (parallel_tria != nullptr &&
          dynamic_cast<const parallel::shared::Triangulation<dim, spacedim> *>(
            &dof_handler.get_triangulation()) == nullptr)
        {
#ifndef DEAL_II_WITH_P4EST
          DEAL_II_ASSERT_UNREACHABLE();
#else
          //
          // parallel implementation with distributed memory
          //

          MPI_Comm mpi_communicator = parallel_tria->get_mpi_communicator();

          // 2.) Communicate the number of cells scheduled for p-adaptation
          //     globally.
          const unsigned int n_global_flags_refinement =
            Utilities::MPI::sum(n_flags_refinement, mpi_communicator);
          const unsigned int n_global_flags_coarsening =
            Utilities::MPI::sum(n_flags_coarsening, mpi_communicator);

          const unsigned int target_index_refinement =
            static_cast<unsigned int>(
              std::floor(p_refine_fraction * n_global_flags_refinement));
          const unsigned int target_index_coarsening =
            static_cast<unsigned int>(
              std::ceil((1 - p_coarsen_fraction) * n_global_flags_coarsening));

          // 3.) Figure out the global max and min of the criteria. We don't
          //     need it here, but it's a collective communication call.
          const std::pair<Number, Number> global_min_max_refinement =
            internal::parallel::distributed::GridRefinement::
              compute_global_min_and_max_at_root(indicators_refinement,
                                                 mpi_communicator);

          const std::pair<Number, Number> global_min_max_coarsening =
            internal::parallel::distributed::GridRefinement::
              compute_global_min_and_max_at_root(indicators_coarsening,
                                                 mpi_communicator);

          // 3.) Compute thresholds if necessary.
          if (target_index_refinement == 0)
            reference_compare_refine = std::cref(compare_false);
          else if (target_index_refinement == n_global_flags_refinement)
            reference_compare_refine = std::cref(compare_true);
          else
            threshold_refinement = internal::parallel::distributed::
              GridRefinement::RefineAndCoarsenFixedNumber::compute_threshold(
                indicators_refinement,
                global_min_max_refinement,
                target_index_refinement,
                mpi_communicator);

          if (target_index_coarsening == n_global_flags_coarsening)
            reference_compare_coarsen = std::cref(compare_false);
          else if (target_index_coarsening == 0)
            reference_compare_coarsen = std::cref(compare_true);
          else
            threshold_coarsening = internal::parallel::distributed::
              GridRefinement::RefineAndCoarsenFixedNumber::compute_threshold(
                indicators_coarsening,
                global_min_max_coarsening,
                target_index_coarsening,
                mpi_communicator);
#endif
        }
      else
        {
          //
          // serial implementation (and parallel::shared implementation)
          //

          // 2.) Determine the number of cells scheduled for p-adaptation.
          const unsigned int n_p_refine_cells = static_cast<unsigned int>(
            std::floor(p_refine_fraction * n_flags_refinement));
          const unsigned int n_p_coarsen_cells = static_cast<unsigned int>(
            std::floor(p_coarsen_fraction * n_flags_coarsening));

          // 3.) Compute thresholds if necessary.
          if (n_p_refine_cells == 0)
            reference_compare_refine = std::cref(compare_false);
          else if (n_p_refine_cells == n_flags_refinement)
            reference_compare_refine = std::cref(compare_true);
          else
            {
              std::nth_element(indicators_refinement.begin(),
                               indicators_refinement.begin() +
                                 n_p_refine_cells - 1,
                               indicators_refinement.end(),
                               std::greater<Number>());
              threshold_refinement =
                *(indicators_refinement.begin() + n_p_refine_cells - 1);
            }

          if (n_p_coarsen_cells == 0)
            reference_compare_coarsen = std::cref(compare_false);
          else if (n_p_coarsen_cells == n_flags_coarsening)
            reference_compare_coarsen = std::cref(compare_true);
          else
            {
              std::nth_element(indicators_coarsening.begin(),
                               indicators_coarsening.begin() +
                                 n_p_coarsen_cells - 1,
                               indicators_coarsening.end(),
                               std::less<Number>());
              threshold_coarsening =
                *(indicators_coarsening.begin() + n_p_coarsen_cells - 1);
            }
        }

      // 4.) Finally perform adaptation.
      p_adaptivity_from_absolute_threshold(dof_handler,
                                           criteria,
                                           threshold_refinement,
                                           threshold_coarsening,
                                           std::cref(reference_compare_refine),
                                           std::cref(
                                             reference_compare_coarsen));
    }



    template <int dim, typename Number, int spacedim>
    void
    p_adaptivity_from_regularity(const DoFHandler<dim, spacedim> &dof_handler,
                                 const Vector<Number> &sobolev_indices)
    {
      if (dof_handler.get_fe_collection().empty())
        // nothing to do
        return;

      Assert(dof_handler.has_hp_capabilities(),
             (typename DoFHandler<dim, spacedim>::ExcOnlyAvailableWithHP()));
      AssertDimension(dof_handler.get_triangulation().n_active_cells(),
                      sobolev_indices.size());

      for (const auto &cell : dof_handler.active_cell_iterators())
        if (cell->is_locally_owned())
          {
            if (cell->refine_flag_set())
              {
                const unsigned int super_fe_index =
                  dof_handler.get_fe_collection().next_in_hierarchy(
                    cell->active_fe_index());

                // Reject update if already most superordinate element.
                if (super_fe_index != cell->active_fe_index())
                  {
                    const unsigned int super_fe_degree =
                      dof_handler.get_fe_collection()[super_fe_index].degree;

                    if (sobolev_indices[cell->active_cell_index()] >
                        super_fe_degree)
                      cell->set_future_fe_index(super_fe_index);
                  }
              }
            else if (cell->coarsen_flag_set())
              {
                const unsigned int sub_fe_index =
                  dof_handler.get_fe_collection().previous_in_hierarchy(
                    cell->active_fe_index());

                // Reject update if already least subordinate element.
                if (sub_fe_index != cell->active_fe_index())
                  {
                    const unsigned int sub_fe_degree =
                      dof_handler.get_fe_collection()[sub_fe_index].degree;

                    if (sobolev_indices[cell->active_cell_index()] <
                        sub_fe_degree)
                      cell->set_future_fe_index(sub_fe_index);
                  }
              }
          }
    }



    template <int dim, typename Number, int spacedim>
    void
    p_adaptivity_from_reference(
      const DoFHandler<dim, spacedim> &dof_handler,
      const Vector<Number>            &criteria,
      const Vector<Number>            &references,
      const ComparisonFunction<std_cxx20::type_identity_t<Number>>
        &compare_refine,
      const ComparisonFunction<std_cxx20::type_identity_t<Number>>
        &compare_coarsen)
    {
      if (dof_handler.get_fe_collection().empty())
        // nothing to do
        return;

      Assert(dof_handler.has_hp_capabilities(),
             (typename DoFHandler<dim, spacedim>::ExcOnlyAvailableWithHP()));
      AssertDimension(dof_handler.get_triangulation().n_active_cells(),
                      criteria.size());
      AssertDimension(dof_handler.get_triangulation().n_active_cells(),
                      references.size());

      std::vector<bool> p_flags(
        dof_handler.get_triangulation().n_active_cells(), false);

      for (const auto &cell : dof_handler.active_cell_iterators())
        if (cell->is_locally_owned() &&
            ((cell->refine_flag_set() &&
              compare_refine(criteria[cell->active_cell_index()],
                             references[cell->active_cell_index()])) ||
             (cell->coarsen_flag_set() &&
              compare_coarsen(criteria[cell->active_cell_index()],
                              references[cell->active_cell_index()]))))
          p_flags[cell->active_cell_index()] = true;

      p_adaptivity_from_flags(dof_handler, p_flags);
    }



    /**
     * Error prediction
     */
    template <int dim, typename Number, int spacedim>
    void
    predict_error(const DoFHandler<dim, spacedim> &dof_handler,
                  const Vector<Number>            &error_indicators,
                  Vector<Number>                  &predicted_errors,
                  const double                     gamma_p,
                  const double                     gamma_h,
                  const double                     gamma_n)
    {
      if (dof_handler.get_fe_collection().empty())
        // nothing to do
        return;

      AssertDimension(dof_handler.get_triangulation().n_active_cells(),
                      error_indicators.size());
      AssertDimension(dof_handler.get_triangulation().n_active_cells(),
                      predicted_errors.size());
      Assert(0 < gamma_p && gamma_p < 1,
             GridRefinement::ExcInvalidParameterValue());
      Assert(0 < gamma_h, GridRefinement::ExcInvalidParameterValue());
      Assert(0 < gamma_n, GridRefinement::ExcInvalidParameterValue());

      // auxiliary variables
      unsigned int future_fe_degree       = numbers::invalid_unsigned_int;
      unsigned int parent_future_fe_index = numbers::invalid_unsigned_int;
      // store all determined future finite element indices on parent cells for
      // coarsening
      std::map<typename DoFHandler<dim, spacedim>::cell_iterator, unsigned int>
        future_fe_indices_on_coarsened_cells;

      // deep copy error indicators
      predicted_errors = error_indicators;

      for (const auto &cell : dof_handler.active_cell_iterators() |
                                IteratorFilters::LocallyOwnedCell())
        {
          // current cell will not be adapted
          if (!(cell->future_fe_index_set()) && !(cell->refine_flag_set()) &&
              !(cell->coarsen_flag_set()))
            {
              predicted_errors[cell->active_cell_index()] *= gamma_n;
              continue;
            }

          // current cell will be adapted
          // determine degree of its future finite element
          if (cell->coarsen_flag_set())
            {
              Assert(cell->level() > 0,
                     ExcMessage("A coarse cell is flagged for coarsening. "
                                "Please read the note in the documentation "
                                "of predict_error()."));

              // cell will be coarsened, thus determine future finite element
              // on parent cell
              const auto &parent = cell->parent();
              if (future_fe_indices_on_coarsened_cells.find(parent) ==
                  future_fe_indices_on_coarsened_cells.end())
                {
                  if constexpr (running_in_debug_mode())
                    {
                      for (const auto &child : parent->child_iterators())
                        Assert(child->is_active() && child->coarsen_flag_set(),
                               typename Triangulation<
                                 dim>::ExcInconsistentCoarseningFlags());
                    }

                  parent_future_fe_index =
                    internal::hp::DoFHandlerImplementation::
                      dominated_future_fe_on_children<dim, spacedim>(parent);

                  future_fe_indices_on_coarsened_cells.insert(
                    {parent, parent_future_fe_index});
                }
              else
                {
                  parent_future_fe_index =
                    future_fe_indices_on_coarsened_cells[parent];
                }

              future_fe_degree =
                dof_handler.get_fe_collection()[parent_future_fe_index].degree;
            }
          else
            {
              // future finite element on current cell is already set
              future_fe_degree =
                dof_handler.get_fe_collection()[cell->future_fe_index()].degree;
            }

          // step 1: exponential decay with p-adaptation
          if (cell->future_fe_index_set())
            {
              if (future_fe_degree > cell->get_fe().degree)
                predicted_errors[cell->active_cell_index()] *=
                  Utilities::pow(gamma_p,
                                 future_fe_degree - cell->get_fe().degree);
              else if (future_fe_degree < cell->get_fe().degree)
                predicted_errors[cell->active_cell_index()] /=
                  Utilities::pow(gamma_p,
                                 cell->get_fe().degree - future_fe_degree);
              else
                {
                  // The two degrees are the same; we do not need to
                  // adapt the predicted error
                }
            }

          // step 2: algebraic decay with h-adaptation
          if (cell->refine_flag_set())
            {
              predicted_errors[cell->active_cell_index()] *=
                (gamma_h * Utilities::pow(.5, future_fe_degree));

              // predicted error will be split on children cells
              // after adaptation via CellDataTransfer
            }
          else if (cell->coarsen_flag_set())
            {
              predicted_errors[cell->active_cell_index()] /=
                (gamma_h * Utilities::pow(.5, future_fe_degree));

              // predicted error will be summed up on parent cell
              // after adaptation via CellDataTransfer
            }
        }
    }



    /**
     * Decide between h- and p-adaptivity
     */
    template <int dim, int spacedim>
    void
    force_p_over_h(const DoFHandler<dim, spacedim> &dof_handler)
    {
      if (dof_handler.get_fe_collection().empty())
        // nothing to do
        return;

      Assert(dof_handler.has_hp_capabilities(),
             (typename DoFHandler<dim, spacedim>::ExcOnlyAvailableWithHP()));

      for (const auto &cell : dof_handler.active_cell_iterators())
        if (cell->is_locally_owned() && cell->future_fe_index_set())
          {
            cell->clear_refine_flag();
            cell->clear_coarsen_flag();
          }
    }



    template <int dim, int spacedim>
    void
    choose_p_over_h(const DoFHandler<dim, spacedim> &dof_handler)
    {
      if (dof_handler.get_fe_collection().empty())
        // nothing to do
        return;

      Assert(dof_handler.has_hp_capabilities(),
             (typename DoFHandler<dim, spacedim>::ExcOnlyAvailableWithHP()));

      // Ghost siblings might occur on parallel Triangulation objects.
      // We need information about refinement flags and future FE indices
      // on all locally relevant cells here, and thus communicate them.
      if (dealii::parallel::distributed::Triangulation<dim, spacedim> *tria =
            dynamic_cast<
              dealii::parallel::distributed::Triangulation<dim, spacedim> *>(
              const_cast<dealii::Triangulation<dim, spacedim> *>(
                &dof_handler.get_triangulation())))
        {
          dealii::internal::parallel::distributed::TriangulationImplementation::
            exchange_refinement_flags(*tria);
        }

      internal::hp::DoFHandlerImplementation::communicate_future_fe_indices(
        const_cast<DoFHandler<dim, spacedim> &>(dof_handler));

      // Now: choose p-adaptation over h-adaptation.
      for (const auto &cell : dof_handler.active_cell_iterators())
        if (cell->is_locally_owned() && cell->future_fe_index_set())
          {
            // This cell is flagged for p-adaptation.

            // Remove any h-refinement flags.
            cell->clear_refine_flag();

            // A cell will only be coarsened into its parent if all of its
            // siblings are flagged for h-coarsening as well. We must take this
            // into account for our decision whether we would like to impose h-
            // or p-adaptivity.
            if (cell->coarsen_flag_set())
              {
                if (cell->level() == 0)
                  {
                    // This cell is a coarse cell and has neither parent nor
                    // siblings, thus it cannot be h-coarsened.
                    // Clear the flag and move on to the next cell.
                    cell->clear_coarsen_flag();
                    continue;
                  }

                const auto        &parent     = cell->parent();
                const unsigned int n_children = parent->n_children();

                unsigned int h_flagged_children = 0, p_flagged_children = 0;
                for (const auto &child : parent->child_iterators())
                  {
                    if (child->is_active())
                      {
                        Assert(child->is_artificial() == false,
                               ExcInternalError());

                        if (child->coarsen_flag_set())
                          ++h_flagged_children;

                        // The public interface does not allow to access
                        // future FE indices on ghost cells. However, we
                        // need this information here and thus call the
                        // internal function that does not check for cell
                        // ownership.
                        if (internal::DoFCellAccessorImplementation::
                              Implementation::
                                future_fe_index_set<dim, spacedim, false>(
                                  *child))
                          ++p_flagged_children;
                      }
                  }

                if (h_flagged_children == n_children &&
                    p_flagged_children != n_children)
                  {
                    // Perform pure h-coarsening and
                    // drop all p-adaptation flags.
                    for (const auto &child : parent->child_iterators())
                      {
                        // h_flagged_children == n_children implies
                        // that all children are active
                        Assert(child->is_active(), ExcInternalError());
                        if (child->is_locally_owned())
                          child->clear_future_fe_index();
                      }
                  }
                else
                  {
                    // Perform p-adaptation (if scheduled) and
                    // drop all h-coarsening flags.
                    for (const auto &child : parent->child_iterators())
                      {
                        if (child->is_active() && child->is_locally_owned())
                          child->clear_coarsen_flag();
                      }
                  }
              }
          }
    }



    /**
     * Optimize p-level distribution
     */
    template <int dim, int spacedim>
    bool
    limit_p_level_difference(const DoFHandler<dim, spacedim> &dof_handler,
                             const unsigned int               max_difference,
                             const unsigned int               contains_fe_index)
    {
      if (dof_handler.get_fe_collection().empty())
        // nothing to do
        return false;

      Assert(dof_handler.has_hp_capabilities(),
             (typename DoFHandler<dim, spacedim>::ExcOnlyAvailableWithHP()));
      Assert(
        max_difference > 0,
        ExcMessage(
          "This function does not serve any purpose for max_difference = 0."));
      AssertIndexRange(contains_fe_index,
                       dof_handler.get_fe_collection().size());

      //
      // establish hierarchy
      //
      // - create bimap between hierarchy levels and FE indices

      // there can be as many levels in the hierarchy as active FE indices are
      // possible
      using level_type         = types::fe_index;
      const auto invalid_level = static_cast<level_type>(-1);

      // map from FE index to level in hierarchy
      // FE indices that are not covered in the hierarchy are not in the map
      const std::vector<unsigned int> fe_index_for_hierarchy_level =
        dof_handler.get_fe_collection().get_hierarchy_sequence(
          contains_fe_index);

      // map from level in hierarchy to FE index
      // FE indices that are not covered in the hierarchy will be mapped to
      // invalid_level
      std::vector<unsigned int> hierarchy_level_for_fe_index(
        dof_handler.get_fe_collection().size(), invalid_level);
      for (unsigned int l = 0; l < fe_index_for_hierarchy_level.size(); ++l)
        hierarchy_level_for_fe_index[fe_index_for_hierarchy_level[l]] = l;


      //
      // parallelization
      //
      // - create distributed vector of level indices
      // - update ghost values in each iteration (see later)
      // - no need to compress, since the owning processor will have the correct
      //   level index

      // HOTFIX: dealii::Vector does not accept integral types
      LinearAlgebra::distributed::Vector<float> future_levels;
      if (const auto parallel_tria =
            dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
              &(dof_handler.get_triangulation())))
        {
          future_levels.reinit(
            parallel_tria->global_active_cell_index_partitioner().lock());
        }
      else
        {
          future_levels.reinit(
            dof_handler.get_triangulation().n_active_cells());
        }

      for (const auto &cell : dof_handler.active_cell_iterators() |
                                IteratorFilters::LocallyOwnedCell())
        future_levels[cell->global_active_cell_index()] =
          hierarchy_level_for_fe_index[cell->future_fe_index()];


      //
      // limit level difference of neighboring cells
      //
      // - go over all locally relevant cells, and adjust the level indices of
      //   locally owned neighbors to match the level difference (as a
      //   consequence, indices on ghost cells will be updated only on the
      //   owning processor)
      // - always raise levels to match criterion, never lower them
      // - exchange level indices on ghost cells

      // Function that updates the level of neighbor to fulfill difference
      // criterion, and returns whether it was changed.
      const auto update_neighbor_level =
        [&future_levels, max_difference, invalid_level](
          const auto &neighbor, const level_type cell_level) -> bool {
        Assert(neighbor->is_active(), ExcInternalError());
        // We only care about locally owned neighbors. If neighbor is a ghost
        // cell, its future FE index will be updated on the owning process and
        // communicated at the next loop iteration.
        if (neighbor->is_locally_owned())
          {
            const level_type neighbor_level = static_cast<level_type>(
              future_levels[neighbor->global_active_cell_index()]);

            // ignore neighbors that are not part of the hierarchy
            if (neighbor_level == invalid_level)
              return false;

            if ((cell_level - max_difference) > neighbor_level)
              {
                future_levels[neighbor->global_active_cell_index()] =
                  cell_level - max_difference;

                return true;
              }
          }

        return false;
      };

      // For cells to be h-coarsened, we need to determine a future FE for the
      // parent cell, which will be the dominated FE among all children
      // However, if we want to enforce the max_difference criterion on all
      // cells on the updated mesh, we will need to simulate the updated mesh on
      // the current mesh.
      //
      // As we are working on p-levels, we will set all siblings that will be
      // coarsened to the highest p-level among them. The parent cell will be
      // assigned exactly this level in form of the corresponding FE index in
      // the adaptation process in
      // Triangulation::execute_coarsening_and_refinement().
      //
      // This function takes a cell and sets all its siblings to the highest
      // p-level among them. Returns whether any future levels have been
      // changed.
      const auto prepare_level_for_parent = [&](const auto &neighbor) -> bool {
        Assert(neighbor->is_active(), ExcInternalError());
        if (neighbor->coarsen_flag_set() && neighbor->is_locally_owned())
          {
            const auto parent = neighbor->parent();

            std::vector<unsigned int> future_levels_children;
            future_levels_children.reserve(parent->n_children());
            for (const auto &child : parent->child_iterators())
              {
                Assert(child->is_active() && child->coarsen_flag_set(),
                       (typename Triangulation<dim, spacedim>::
                          ExcInconsistentCoarseningFlags()));

                const level_type child_level = static_cast<level_type>(
                  future_levels[child->global_active_cell_index()]);
                Assert(child_level != invalid_level,
                       ExcMessage(
                         "The FiniteElement on one of the siblings of "
                         "a cell you are trying to coarsen is not part "
                         "of the registered p-adaptation hierarchy."));
                future_levels_children.push_back(child_level);
              }
            Assert(!future_levels_children.empty(), ExcInternalError());

            const unsigned int max_level_children =
              *std::max_element(future_levels_children.begin(),
                                future_levels_children.end());

            bool children_changed = false;
            for (const auto &child : parent->child_iterators())
              // We only care about locally owned children. If child is a ghost
              // cell, its future FE index will be updated on the owning process
              // and communicated at the next loop iteration.
              if (child->is_locally_owned() &&
                  future_levels[child->global_active_cell_index()] !=
                    max_level_children)
                {
                  future_levels[child->global_active_cell_index()] =
                    max_level_children;
                  children_changed = true;
                }
            return children_changed;
          }

        return false;
      };

      bool levels_changed = false;
      bool levels_changed_in_cycle;
      do
        {
          levels_changed_in_cycle = false;

          future_levels.update_ghost_values();

          for (const auto &cell : dof_handler.active_cell_iterators())
            if (!cell->is_artificial())
              {
                const level_type cell_level = static_cast<level_type>(
                  future_levels[cell->global_active_cell_index()]);

                // ignore cells that are not part of the hierarchy
                if (cell_level == invalid_level)
                  continue;

                // ignore lowest levels of the hierarchy that always fulfill the
                // max_difference criterion
                if (cell_level <= max_difference)
                  continue;

                for (unsigned int f = 0; f < cell->n_faces(); ++f)
                  if (cell->face(f)->at_boundary() == false)
                    {
                      if (cell->face(f)->has_children())
                        {
                          for (unsigned int sf = 0;
                               sf < cell->face(f)->n_children();
                               ++sf)
                            {
                              const auto neighbor =
                                cell->neighbor_child_on_subface(f, sf);

                              levels_changed_in_cycle |=
                                update_neighbor_level(neighbor, cell_level);

                              levels_changed_in_cycle |=
                                prepare_level_for_parent(neighbor);
                            }
                        }
                      else
                        {
                          const auto neighbor = cell->neighbor(f);

                          levels_changed_in_cycle |=
                            update_neighbor_level(neighbor, cell_level);

                          levels_changed_in_cycle |=
                            prepare_level_for_parent(neighbor);
                        }
                    }
              }

          levels_changed_in_cycle =
            Utilities::MPI::logical_or(levels_changed_in_cycle,
                                       dof_handler.get_mpi_communicator());
          levels_changed |= levels_changed_in_cycle;
        }
      while (levels_changed_in_cycle);

      // update future FE indices on locally owned cells
      for (const auto &cell : dof_handler.active_cell_iterators() |
                                IteratorFilters::LocallyOwnedCell())
        {
          const level_type cell_level = static_cast<level_type>(
            future_levels[cell->global_active_cell_index()]);

          if (cell_level != invalid_level)
            {
              const unsigned int fe_index =
                fe_index_for_hierarchy_level[cell_level];

              if (fe_index != cell->active_fe_index())
                cell->set_future_fe_index(fe_index);
              else
                cell->clear_future_fe_index();
            }
        }

      return levels_changed;
    }
  } // namespace Refinement
} // namespace hp


// explicit instantiations
#include "hp/refinement.inst"

DEAL_II_NAMESPACE_CLOSE
