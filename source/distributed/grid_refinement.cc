// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/utilities.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#ifdef DEAL_II_WITH_P4EST

#  include <deal.II/distributed/grid_refinement.h>

#  include <deal.II/grid/filtered_iterator.h>
#  include <deal.II/grid/grid_refinement.h>
#  include <deal.II/grid/tria.h>
#  include <deal.II/grid/tria_accessor.h>
#  include <deal.II/grid/tria_iterator.h>

#  include <algorithm>
#  include <functional>
#  include <limits>
#  include <numeric>


DEAL_II_NAMESPACE_OPEN


namespace
{
  template <int dim, int spacedim>
  unsigned int
  n_locally_owned_active_cells(const Triangulation<dim, spacedim> &tria)
  {
    if (const auto parallel_tria =
          dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
            &tria))
      return parallel_tria->n_locally_owned_active_cells();
    else
      return tria.n_active_cells();
  }

  template <typename number>
  inline number
  max_element(const dealii::Vector<number> &criteria)
  {
    return (criteria.size() > 0) ?
             (*std::max_element(criteria.begin(), criteria.end())) :
             std::numeric_limits<number>::min();
  }



  template <typename number>
  inline number
  min_element(const dealii::Vector<number> &criteria)
  {
    return (criteria.size() > 0) ?
             (*std::min_element(criteria.begin(), criteria.end())) :
             std::numeric_limits<number>::max();
  }



  /**
   * Compute the global sum over the elements of the vectors passed to this
   * function on all processors. This number is returned only on the processor
   * with rank zero, all others get zero.
   */
  template <typename number>
  double
  compute_global_sum(const dealii::Vector<number> &criteria,
                     const MPI_Comm                mpi_communicator)
  {
    double my_sum =
      std::accumulate(criteria.begin(),
                      criteria.end(),
                      /* do accumulation in the correct data type: */
                      number());

    double result = 0;
    // compute the minimum on processor zero
    const int ierr =
      MPI_Reduce(&my_sum, &result, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_communicator);
    AssertThrowMPI(ierr);

    // make sure only processor zero got something
    if (Utilities::MPI::this_mpi_process(mpi_communicator) != 0)
      Assert(result == 0, ExcInternalError());

    return result;
  }



  /**
   * Given a vector of refinement criteria for all cells of a mesh (locally
   * owned or not), extract those that pertain to locally owned cells.
   */
  template <int dim, int spacedim, typename Number>
  void
  get_locally_owned_indicators(const dealii::Triangulation<dim, spacedim> &tria,
                               const dealii::Vector<Number> &criteria,
                               Vector<Number> &locally_owned_indicators)
  {
    Assert(locally_owned_indicators.size() ==
             n_locally_owned_active_cells(tria),
           ExcInternalError());

    unsigned int owned_index = 0;
    for (const auto &cell :
         tria.active_cell_iterators() | IteratorFilters::LocallyOwnedCell())
      {
        locally_owned_indicators(owned_index) =
          criteria(cell->active_cell_index());
        ++owned_index;
      }
    Assert(owned_index == n_locally_owned_active_cells(tria),
           ExcInternalError());
  }


  // we compute refinement thresholds by bisection of the interval spanned by
  // the smallest and largest error indicator. this leads to a small problem:
  // if, for example, we want to coarsen zero per cent of the cells, then we
  // need to pick a threshold equal to the smallest indicator, but of course
  // the bisection algorithm can never find a threshold equal to one of the
  // end points of the interval. So we slightly increase the interval before
  // we even start
  void
  adjust_interesting_range(double (&interesting_range)[2])
  {
    Assert(interesting_range[0] <= interesting_range[1], ExcInternalError());

    if (interesting_range[0] > 0)
      {
        // In this case, we calculate the first interval split point `m` in the
        // `compute_threshold` functions in the optimized way: We exploit that
        // the logarithms of all criteria are more uniformly distributed than
        // their actual values, i.e. m=sqrt(b*e).
        //
        // Both factors will modify the split point only slightly by a factor of
        // sqrt(1.01*0.99) = sqrt(0.9999) ~ 0.9950.
        interesting_range[0] *= 0.99;
        interesting_range[1] *= 1.01;
      }
    else
      {
        // In all other cases, we begin with an the arithmetic mean as the
        // standard interval split point, i.e. m=(b+e)/2.
        //
        // Both increments will add up to zero when calculating the initial
        // split point in the `compute_threshold` functions.
        const double difference =
          std::abs(interesting_range[1] - interesting_range[0]);
        interesting_range[0] -= 0.01 * difference;
        interesting_range[1] += 0.01 * difference;
      }
  }



  /**
   * Given a vector of criteria and bottom and top thresholds for coarsening and
   * refinement, mark all those cells that we locally own as appropriate for
   * coarsening or refinement.
   */
  template <int dim, int spacedim, typename Number>
  void
  mark_cells(dealii::Triangulation<dim, spacedim> &tria,
             const dealii::Vector<Number>         &criteria,
             const double                          top_threshold,
             const double                          bottom_threshold)
  {
    dealii::GridRefinement::refine(tria, criteria, top_threshold);
    dealii::GridRefinement::coarsen(tria, criteria, bottom_threshold);

    // as a final good measure, delete all flags again from cells that we don't
    // locally own
    for (const auto &cell : tria.active_cell_iterators())
      if (cell->subdomain_id() != tria.locally_owned_subdomain())
        {
          cell->clear_refine_flag();
          cell->clear_coarsen_flag();
        }
  }



  /**
   * Fixed fraction algorithm without a specified vector norm.
   *
   * Entries of the criteria vector and fractions are taken as is, so this
   * function basically evaluates norms on the vector or its subsets as
   * l1-norms.
   */
  template <int dim, int spacedim, typename Number>
  void
  refine_and_coarsen_fixed_fraction_via_l1_norm(
    dealii::Triangulation<dim, spacedim> &tria,
    const dealii::Vector<Number>         &criteria,
    const double                          top_fraction_of_error,
    const double                          bottom_fraction_of_error)
  {
    // first extract from the vector of indicators the ones that correspond
    // to cells that we locally own
    Vector<Number> locally_owned_indicators(n_locally_owned_active_cells(tria));
    get_locally_owned_indicators(tria, criteria, locally_owned_indicators);

    MPI_Comm mpi_communicator = tria.get_mpi_communicator();

    // figure out the global max and min of the indicators. we don't need it
    // here, but it's a collective communication call
    const std::pair<double, double> global_min_and_max =
      dealii::internal::parallel::distributed::GridRefinement::
        compute_global_min_and_max_at_root(locally_owned_indicators,
                                           mpi_communicator);

    const double total_error =
      compute_global_sum(locally_owned_indicators, mpi_communicator);

    double top_target_error    = top_fraction_of_error * total_error,
           bottom_target_error = (1. - bottom_fraction_of_error) * total_error;

    double top_threshold, bottom_threshold;
    top_threshold = dealii::internal::parallel::distributed::GridRefinement::
      RefineAndCoarsenFixedFraction::compute_threshold(locally_owned_indicators,
                                                       global_min_and_max,
                                                       top_target_error,
                                                       mpi_communicator);

    // compute bottom threshold only if necessary. otherwise use the lowest
    // threshold possible
    if (bottom_fraction_of_error > 0)
      bottom_threshold = dealii::internal::parallel::distributed::
        GridRefinement::RefineAndCoarsenFixedFraction::compute_threshold(
          locally_owned_indicators,
          global_min_and_max,
          bottom_target_error,
          mpi_communicator);
    else
      bottom_threshold = std::numeric_limits<Number>::lowest();

    // now refine the mesh
    mark_cells(tria, criteria, top_threshold, bottom_threshold);
  }
} // namespace



namespace internal
{
  namespace parallel
  {
    namespace distributed
    {
      namespace GridRefinement
      {
        template <typename number>
        std::pair<number, number>
        compute_global_min_and_max_at_root(
          const dealii::Vector<number> &criteria,
          const MPI_Comm                mpi_communicator)
        {
          // we'd like to compute the global max and min from the local ones in
          // one MPI communication. we can do that by taking the elementwise
          // minimum of the local min and the negative maximum over all
          // processors

          const double local_min = min_element(criteria),
                       local_max = max_element(criteria);
          double comp[2]         = {local_min, -local_max};
          double result[2]       = {0, 0};

          // compute the minimum on processor zero
          const int ierr = MPI_Reduce(
            comp, result, 2, MPI_DOUBLE, MPI_MIN, 0, mpi_communicator);
          AssertThrowMPI(ierr);

          // make sure only processor zero got something
          if (Utilities::MPI::this_mpi_process(mpi_communicator) != 0)
            Assert((result[0] == 0) && (result[1] == 0), ExcInternalError());

          return std::make_pair(result[0], -result[1]);
        }



        namespace RefineAndCoarsenFixedNumber
        {
          template <typename number>
          number
          compute_threshold(const dealii::Vector<number>    &criteria,
                            const std::pair<double, double> &global_min_and_max,
                            const types::global_cell_index   n_target_cells,
                            const MPI_Comm                   mpi_communicator)
          {
            double interesting_range[2] = {global_min_and_max.first,
                                           global_min_and_max.second};
            adjust_interesting_range(interesting_range);

            const unsigned int root_mpi_rank = 0;
            unsigned int       iteration     = 0;

            do
              {
                int ierr = MPI_Bcast(interesting_range,
                                     2,
                                     MPI_DOUBLE,
                                     root_mpi_rank,
                                     mpi_communicator);
                AssertThrowMPI(ierr);

                if (interesting_range[0] == interesting_range[1])
                  return interesting_range[0];

                const double test_threshold =
                  (interesting_range[0] > 0 ?
                     std::sqrt(interesting_range[0] * interesting_range[1]) :
                     (interesting_range[0] + interesting_range[1]) / 2);

                // Count how many of our own elements would be above this
                // threshold. Use a 64bit result type if we are compiling with
                // 64bit indices to avoid an overflow when computing the sum
                // below.
                const types::global_cell_index my_count =
                  std::count_if(criteria.begin(),
                                criteria.end(),
                                [test_threshold](const double c) {
                                  return c > test_threshold;
                                });
                const types::global_cell_index total_count =
                  Utilities::MPI::sum(my_count, mpi_communicator);

                // now adjust the range. if we have too many cells, we take the
                // upper half of the previous range, otherwise the lower half.
                // if we have hit the right number, then set the range to the
                // exact value. non-root nodes also update their own
                // interesting_range, however their results are not significant
                // since the values will be overwritten by MPI_Bcast from the
                // root node in next loop.
                if (total_count > n_target_cells)
                  interesting_range[0] = test_threshold;
                else if (total_count < n_target_cells)
                  interesting_range[1] = test_threshold;
                else
                  interesting_range[0] = interesting_range[1] = test_threshold;

                // terminate the iteration after 25 go-arounds. this is
                // necessary because oftentimes error indicators on cells have
                // exactly the same value, and so there may not be a particular
                // value that cuts the indicators in such a way that we can
                // achieve the desired number of cells. using a maximal number
                // of iterations means that we terminate the iteration after a
                // fixed number N of steps if the indicators were perfectly
                // badly distributed, and we make at most a mistake of 1/2^N in
                // the number of cells flagged if indicators are perfectly
                // equidistributed
                ++iteration;
                if (iteration == 25)
                  interesting_range[0] = interesting_range[1] = test_threshold;
              }
            while (true);

            DEAL_II_ASSERT_UNREACHABLE();
            return -1;
          }
        } // namespace RefineAndCoarsenFixedNumber



        namespace RefineAndCoarsenFixedFraction
        {
          template <typename number>
          number
          compute_threshold(const dealii::Vector<number>    &criteria,
                            const std::pair<double, double> &global_min_and_max,
                            const double                     target_error,
                            const MPI_Comm                   mpi_communicator)
          {
            double interesting_range[2] = {global_min_and_max.first,
                                           global_min_and_max.second};
            adjust_interesting_range(interesting_range);

            const unsigned int root_mpi_rank = 0;
            unsigned int       iteration     = 0;

            do
              {
                int ierr = MPI_Bcast(interesting_range,
                                     2,
                                     MPI_DOUBLE,
                                     root_mpi_rank,
                                     mpi_communicator);
                AssertThrowMPI(ierr);

                if (interesting_range[0] == interesting_range[1])
                  {
                    // so we have found our threshold. since we adjust the range
                    // at the top of the function to be slightly larger than the
                    // actual extremes of the refinement criteria values, we can
                    // end up in a situation where the threshold is in fact
                    // larger than the maximal refinement indicator. in such
                    // cases, we get no refinement at all. thus, cap the
                    // threshold by the actual largest value
                    double final_threshold =
                      std::min(interesting_range[0], global_min_and_max.second);
                    ierr = MPI_Bcast(&final_threshold,
                                     1,
                                     MPI_DOUBLE,
                                     root_mpi_rank,
                                     mpi_communicator);
                    AssertThrowMPI(ierr);

                    return final_threshold;
                  }

                const double test_threshold =
                  (interesting_range[0] > 0 ?
                     std::sqrt(interesting_range[0] * interesting_range[1]) :
                     (interesting_range[0] + interesting_range[1]) / 2);

                // accumulate the error of those our own elements above this
                // threshold and then add to it the number for all the others
                double my_error = 0;
                for (unsigned int i = 0; i < criteria.size(); ++i)
                  if (criteria(i) > test_threshold)
                    my_error += criteria(i);

                double total_error = 0.;

                ierr = MPI_Reduce(&my_error,
                                  &total_error,
                                  1,
                                  MPI_DOUBLE,
                                  MPI_SUM,
                                  root_mpi_rank,
                                  mpi_communicator);
                AssertThrowMPI(ierr);

                // now adjust the range. if we have too many cells, we take the
                // upper half of the previous range, otherwise the lower half.
                // if we have hit the right number, then set the range to the
                // exact value. non-root nodes also update their own
                // interesting_range, however their results are not significant
                // since the values will be overwritten by MPI_Bcast from the
                // root node in next loop.
                if (total_error > target_error)
                  interesting_range[0] = test_threshold;
                else if (total_error < target_error)
                  interesting_range[1] = test_threshold;
                else
                  interesting_range[0] = interesting_range[1] = test_threshold;

                // terminate the iteration after 25 go-arounds. this is
                // necessary because oftentimes error indicators on cells
                // have exactly the same value, and so there may not be a
                // particular value that cuts the indicators in such a way
                // that we can achieve the desired number of cells. using a
                // max of 25 iterations means that we terminate the
                // iteration after 25 steps if the indicators were perfectly
                // badly distributed, and we make at most a mistake of
                // 1/2^25 in the number of cells flagged if indicators are
                // perfectly equidistributed
                ++iteration;
                if (iteration == 25)
                  interesting_range[0] = interesting_range[1] = test_threshold;
              }
            while (true);

            DEAL_II_ASSERT_UNREACHABLE();
            return -1;
          }
        } // namespace RefineAndCoarsenFixedFraction
      }   // namespace GridRefinement
    }     // namespace distributed
  }       // namespace parallel
} // namespace internal



namespace parallel
{
  namespace distributed
  {
    namespace GridRefinement
    {
      template <int dim, typename Number, int spacedim>
      void
      refine_and_coarsen_fixed_number(
        dealii::Triangulation<dim, spacedim> &tria,
        const dealii::Vector<Number>         &criteria,
        const double                          top_fraction_of_cells,
        const double                          bottom_fraction_of_cells,
        const types::global_cell_index        max_n_cells)
      {
        Assert(criteria.size() == tria.n_active_cells(),
               ExcDimensionMismatch(criteria.size(), tria.n_active_cells()));
        Assert((top_fraction_of_cells >= 0) && (top_fraction_of_cells <= 1),
               dealii::GridRefinement::ExcInvalidParameterValue());
        Assert((bottom_fraction_of_cells >= 0) &&
                 (bottom_fraction_of_cells <= 1),
               dealii::GridRefinement::ExcInvalidParameterValue());
        Assert(top_fraction_of_cells + bottom_fraction_of_cells <= 1,
               dealii::GridRefinement::ExcInvalidParameterValue());
        Assert(criteria.is_non_negative(),
               dealii::GridRefinement::ExcNegativeCriteria());

        const std::pair<double, double> adjusted_fractions =
          dealii::GridRefinement::adjust_refine_and_coarsen_number_fraction<
            dim>(tria.n_global_active_cells(),
                 max_n_cells,
                 top_fraction_of_cells,
                 bottom_fraction_of_cells);

        // first extract from the vector of indicators the ones that correspond
        // to cells that we locally own
        Vector<Number> locally_owned_indicators(
          n_locally_owned_active_cells(tria));
        get_locally_owned_indicators(tria, criteria, locally_owned_indicators);

        MPI_Comm mpi_communicator = tria.get_mpi_communicator();

        // figure out the global max and min of the indicators. we don't need it
        // here, but it's a collective communication call
        const std::pair<Number, Number> global_min_and_max =
          dealii::internal::parallel::distributed::GridRefinement::
            compute_global_min_and_max_at_root(locally_owned_indicators,
                                               mpi_communicator);


        double top_threshold, bottom_threshold;
        top_threshold = dealii::internal::parallel::distributed::
          GridRefinement::RefineAndCoarsenFixedNumber::compute_threshold(
            locally_owned_indicators,
            global_min_and_max,
            static_cast<types::global_cell_index>(adjusted_fractions.first *
                                                  tria.n_global_active_cells()),
            mpi_communicator);

        // compute bottom threshold only if necessary. otherwise use the lowest
        // threshold possible
        if (adjusted_fractions.second > 0)
          bottom_threshold = dealii::internal::parallel::distributed::
            GridRefinement::RefineAndCoarsenFixedNumber::compute_threshold(
              locally_owned_indicators,
              global_min_and_max,
              static_cast<types::global_cell_index>(
                std::ceil((1. - adjusted_fractions.second) *
                          tria.n_global_active_cells())),
              mpi_communicator);
        else
          bottom_threshold = std::numeric_limits<Number>::lowest();

        // now refine the mesh
        mark_cells(tria, criteria, top_threshold, bottom_threshold);
      }



      template <int dim, typename Number, int spacedim>
      void
      refine_and_coarsen_fixed_fraction(
        dealii::Triangulation<dim, spacedim> &tria,
        const dealii::Vector<Number>         &criteria,
        const double                          top_fraction_of_error,
        const double                          bottom_fraction_of_error,
        const VectorTools::NormType           norm_type)
      {
        Assert(criteria.size() == tria.n_active_cells(),
               ExcDimensionMismatch(criteria.size(), tria.n_active_cells()));
        Assert((top_fraction_of_error >= 0) && (top_fraction_of_error <= 1),
               dealii::GridRefinement::ExcInvalidParameterValue());
        Assert((bottom_fraction_of_error >= 0) &&
                 (bottom_fraction_of_error <= 1),
               dealii::GridRefinement::ExcInvalidParameterValue());
        Assert(top_fraction_of_error + bottom_fraction_of_error <= 1,
               dealii::GridRefinement::ExcInvalidParameterValue());
        Assert(criteria.is_non_negative(),
               dealii::GridRefinement::ExcNegativeCriteria());

        switch (norm_type)
          {
            case VectorTools::L1_norm:
              // evaluate norms on subsets and compare them as
              //   c_0 + c_1 + ... < fraction * l1-norm(c)
              refine_and_coarsen_fixed_fraction_via_l1_norm(
                tria,
                criteria,
                top_fraction_of_error,
                bottom_fraction_of_error);
              break;

            case VectorTools::L2_norm:
              {
                // we do not want to evaluate norms on subsets as:
                //   sqrt(c_0^2 + c_1^2 + ...) < fraction * l2-norm(c)
                // instead take the square of both sides of the equation
                // and evaluate:
                //   c_0^2 + c_1^2 + ... < fraction^2 * l1-norm(c.c)
                // we adjust all parameters accordingly
                Vector<Number> criteria_squared(criteria.size());
                std::transform(criteria.begin(),
                               criteria.end(),
                               criteria_squared.begin(),
                               [](Number c) { return c * c; });

                refine_and_coarsen_fixed_fraction_via_l1_norm(
                  tria,
                  criteria_squared,
                  top_fraction_of_error * top_fraction_of_error,
                  bottom_fraction_of_error * bottom_fraction_of_error);
              }
              break;

            default:
              DEAL_II_NOT_IMPLEMENTED();
              break;
          }
      }
    } // namespace GridRefinement
  }   // namespace distributed
} // namespace parallel


// explicit instantiations
#  include "distributed/grid_refinement.inst"

DEAL_II_NAMESPACE_CLOSE

#endif
