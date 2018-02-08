// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2017 by the deal.II authors
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


#include <deal.II/base/utilities.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>

#ifdef DEAL_II_WITH_P4EST

#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/distributed/grid_refinement.h>

#include <numeric>
#include <algorithm>
#include <limits>
#include <functional>


DEAL_II_NAMESPACE_OPEN


namespace
{
  template <typename number>
  inline
  number
  max_element (const dealii::Vector<number> &criteria)
  {
    return (criteria.size()>0)
           ?
           (*std::max_element(criteria.begin(), criteria.end()))
           :
           std::numeric_limits<number>::min();
  }



  template <typename number>
  inline
  number
  min_element (const dealii::Vector<number> &criteria)
  {
    return (criteria.size()>0)
           ?
           (*std::min_element(criteria.begin(), criteria.end()))
           :
           std::numeric_limits<number>::max();
  }


  /**
   * Compute the global max and min of the criteria vector. These are returned
   * only on the processor with rank zero, all others get a pair of zeros.
   */
  template <typename number>
  std::pair<number,number>
  compute_global_min_and_max_at_root (const dealii::Vector<number> &criteria,
                                      MPI_Comm                      mpi_communicator)
  {
    // we'd like to compute the global max and min from the local ones in one
    // MPI communication. we can do that by taking the elementwise minimum of
    // the local min and the negative maximum over all processors

    const double local_min = min_element (criteria),
                 local_max = max_element (criteria);
    double comp[2] = { local_min, -local_max };
    double result[2] = { 0, 0 };

    // compute the minimum on processor zero
    const int ierr = MPI_Reduce (comp, result, 2, MPI_DOUBLE,
                                 MPI_MIN, 0, mpi_communicator);
    AssertThrowMPI(ierr);

    // make sure only processor zero got something
    if (Utilities::MPI::this_mpi_process (mpi_communicator) != 0)
      Assert ((result[0] == 0) && (result[1] == 0),
              ExcInternalError());

    return std::make_pair (result[0], -result[1]);
  }



  /**
   * Compute the global sum over the elements of the vectors passed to this
   * function on all processors. This number is returned only on the processor
   * with rank zero, all others get zero.
   */
  template <typename number>
  double
  compute_global_sum (const dealii::Vector<number> &criteria,
                      MPI_Comm                      mpi_communicator)
  {
    double my_sum = std::accumulate (criteria.begin(),
                                     criteria.end(),
                                     /* do accumulation in the correct data type: */
                                     number());

    double result = 0;
    // compute the minimum on processor zero
    const int ierr = MPI_Reduce (&my_sum, &result, 1, MPI_DOUBLE,
                                 MPI_SUM, 0, mpi_communicator);
    AssertThrowMPI(ierr);

    // make sure only processor zero got something
    if (Utilities::MPI::this_mpi_process (mpi_communicator) != 0)
      Assert (result == 0, ExcInternalError());

    return result;
  }



  /**
   * Given a vector of refinement criteria for all cells of a mesh (locally
   * owned or not), extract those that pertain to locally owned cells.
   */
  template <int dim, int spacedim, typename Number>
  void
  get_locally_owned_indicators (const parallel::distributed::Triangulation<dim,spacedim> &tria,
                                const dealii::Vector<Number>                             &criteria,
                                Vector<Number>                  &locally_owned_indicators)
  {
    Assert (locally_owned_indicators.size() == tria.n_locally_owned_active_cells(),
            ExcInternalError());

    unsigned int owned_index = 0;
    for (typename Triangulation<dim,spacedim>::active_cell_iterator
         cell = tria.begin_active();
         cell != tria.end(); ++cell)
      if (cell->subdomain_id() == tria.locally_owned_subdomain())
        {
          locally_owned_indicators(owned_index)
            = criteria(cell->active_cell_index());
          ++owned_index;
        }
    Assert (owned_index == tria.n_locally_owned_active_cells(),
            ExcInternalError());
  }


  // we compute refinement thresholds by bisection of the interval spanned by
  // the smallest and largest error indicator. this leads to a small problem:
  // if, for example, we want to coarsen zero per cent of the cells, then we
  // need to pick a threshold equal to the smallest indicator, but of course
  // the bisection algorithm can never find a threshold equal to one of the
  // end points of the interval. So we slightly increase the interval before
  // we even start
  void adjust_interesting_range (double (&interesting_range)[2])
  {
    Assert (interesting_range[0] <= interesting_range[1],
            ExcInternalError());

    Assert (interesting_range[0] >= 0,
            ExcInternalError());

    // adjust the lower bound only if the end point is not equal to zero,
    // otherwise it could happen, that the result becomes negative
    if (interesting_range[0] > 0)
      interesting_range[0] *= 0.99;

    if (interesting_range[1] > 0)
      interesting_range[1] *= 1.01;
    else
      interesting_range[1]
      += 0.01 * (interesting_range[1] - interesting_range[0]);
  }



  /**
   * Given a vector of criteria and bottom and top thresholds for coarsening and
   * refinement, mark all those cells that we locally own as appropriate for
   * coarsening or refinement.
   */
  template <int dim, int spacedim, typename Number>
  void
  mark_cells (parallel::distributed::Triangulation<dim,spacedim> &tria,
              const dealii::Vector<Number>                       &criteria,
              const double                                        top_threshold,
              const double                                        bottom_threshold)
  {
    dealii::GridRefinement::refine (tria, criteria, top_threshold);
    dealii::GridRefinement::coarsen (tria, criteria, bottom_threshold);

    // as a final good measure, delete all flags again from cells that we don't
    // locally own
    for (typename Triangulation<dim,spacedim>::active_cell_iterator
         cell = tria.begin_active();
         cell != tria.end(); ++cell)
      if (cell->subdomain_id() != tria.locally_owned_subdomain())
        {
          cell->clear_refine_flag ();
          cell->clear_coarsen_flag ();
        }
  }




  namespace RefineAndCoarsenFixedNumber
  {
    /**
     * Compute a threshold value so that exactly n_target_cells have a value
     * that is larger.
     */
    template <typename number>
    number
    compute_threshold (const dealii::Vector<number>   &criteria,
                       const std::pair<double,double>  global_min_and_max,
                       const unsigned int              n_target_cells,
                       MPI_Comm                        mpi_communicator)
    {
      double interesting_range[2] = { global_min_and_max.first,
                                      global_min_and_max.second
                                    };
      adjust_interesting_range (interesting_range);

      const unsigned int master_mpi_rank = 0;
      unsigned int iteration = 0;

      do
        {
          int ierr = MPI_Bcast (interesting_range, 2, MPI_DOUBLE,
                                master_mpi_rank, mpi_communicator);
          AssertThrowMPI(ierr);

          if (interesting_range[0] == interesting_range[1])
            return interesting_range[0];

          const double test_threshold
            = (interesting_range[0] > 0
               ?
               std::sqrt(interesting_range[0] * interesting_range[1])
               :
               (interesting_range[0] + interesting_range[1]) / 2);

          // count how many of our own elements would be above this threshold
          // and then add to it the number for all the others
          unsigned int
          my_count = std::count_if (criteria.begin(),
                                    criteria.end(),
                                    [test_threshold](const double c)
          {
            return c>test_threshold;
          });

          unsigned int total_count;
          ierr = MPI_Reduce (&my_count, &total_count, 1, MPI_UNSIGNED,
                             MPI_SUM, master_mpi_rank, mpi_communicator);
          AssertThrowMPI(ierr);

          // now adjust the range. if we have to many cells, we take the upper
          // half of the previous range, otherwise the lower half. if we have
          // hit the right number, then set the range to the exact value.
          // slave nodes also update their own interesting_range, however their
          // results are not significant since the values will be overwritten by
          // MPI_Bcast from the master node in next loop.
          if (total_count > n_target_cells)
            interesting_range[0] = test_threshold;
          else if (total_count < n_target_cells)
            interesting_range[1] = test_threshold;
          else
            interesting_range[0] = interesting_range[1] = test_threshold;

          // terminate the iteration after 25 go-arounds. this is necessary
          // because oftentimes error indicators on cells have exactly the
          // same value, and so there may not be a particular value that cuts
          // the indicators in such a way that we can achieve the desired
          // number of cells. using a maximal number of iterations means that
          // we terminate the iteration after a fixed number N of steps if the
          // indicators were perfectly badly distributed, and we make at most
          // a mistake of 1/2^N in the number of cells flagged if indicators
          // are perfectly equidistributed
          ++iteration;
          if (iteration == 25)
            interesting_range[0] = interesting_range[1] = test_threshold;
        }
      while (true);

      Assert (false, ExcInternalError());
      return -1;
    }
  }




  namespace RefineAndCoarsenFixedFraction
  {
    /**
     * Compute a threshold value so that the error
     * accumulated over all criteria[i] so that
     *     criteria[i] > threshold
     * is larger than target_error.
     */
    template <typename number>
    number
    compute_threshold (const dealii::Vector<number>   &criteria,
                       const std::pair<double,double>  global_min_and_max,
                       const double                    target_error,
                       MPI_Comm                        mpi_communicator)
    {
      double interesting_range[2] = { global_min_and_max.first,
                                      global_min_and_max.second
                                    };
      adjust_interesting_range (interesting_range);

      const unsigned int master_mpi_rank = 0;
      unsigned int iteration = 0;

      do
        {
          int ierr = MPI_Bcast (interesting_range, 2, MPI_DOUBLE,
                                master_mpi_rank, mpi_communicator);
          AssertThrowMPI(ierr);

          if (interesting_range[0] == interesting_range[1])
            {
              // so we have found our threshold. since we adjust the range
              // at the top of the function to be slightly larger than the
              // actual extremes of the refinement criteria values, we can end
              // up in a situation where the threshold is in fact larger than
              // the maximal refinement indicator. in such cases, we get no
              // refinement at all. thus, cap the threshold by the actual
              // largest value
              double final_threshold =  std::min (interesting_range[0],
                                                  global_min_and_max.second);
              ierr = MPI_Bcast (&final_threshold, 1, MPI_DOUBLE,
                                master_mpi_rank, mpi_communicator);
              AssertThrowMPI(ierr);

              return final_threshold;
            }

          const double test_threshold
            = (interesting_range[0] > 0
               ?
               std::sqrt(interesting_range[0] * interesting_range[1])
               :
               (interesting_range[0] + interesting_range[1]) / 2);

          // accumulate the error of those our own elements above this threshold
          // and then add to it the number for all the others
          double my_error = 0;
          for (unsigned int i=0; i<criteria.size(); ++i)
            if (criteria(i) > test_threshold)
              my_error += criteria(i);

          double total_error;
          ierr = MPI_Reduce (&my_error, &total_error, 1, MPI_DOUBLE,
                             MPI_SUM, master_mpi_rank, mpi_communicator);
          AssertThrowMPI(ierr);

          // now adjust the range. if we have to many cells, we take the upper
          // half of the previous range, otherwise the lower half. if we have
          // hit the right number, then set the range to the exact value.
          // slave nodes also update their own interesting_range, however their
          // results are not significant since the values will be overwritten by
          // MPI_Bcast from the master node in next loop.
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

      Assert (false, ExcInternalError());
      return -1;
    }
  }
}



namespace parallel
{
  namespace distributed
  {
    namespace GridRefinement
    {
      template <int dim, typename Number, int spacedim>
      void
      refine_and_coarsen_fixed_number
      (parallel::distributed::Triangulation<dim,spacedim> &tria,
       const dealii::Vector<Number>                       &criteria,
       const double                                        top_fraction_of_cells,
       const double                                        bottom_fraction_of_cells,
       const unsigned int                                  max_n_cells)
      {
        Assert (criteria.size() == tria.n_active_cells(),
                ExcDimensionMismatch (criteria.size(), tria.n_active_cells()));
        Assert ((top_fraction_of_cells>=0) && (top_fraction_of_cells<=1),
                dealii::GridRefinement::ExcInvalidParameterValue());
        Assert ((bottom_fraction_of_cells>=0) && (bottom_fraction_of_cells<=1),
                dealii::GridRefinement::ExcInvalidParameterValue());
        Assert (top_fraction_of_cells+bottom_fraction_of_cells <= 1,
                dealii::GridRefinement::ExcInvalidParameterValue());
        Assert (criteria.is_non_negative (),
                dealii::GridRefinement::ExcNegativeCriteria());

        const std::pair<double, double> adjusted_fractions =
          dealii::GridRefinement::adjust_refine_and_coarsen_number_fraction<dim> (
            tria.n_global_active_cells(),
            max_n_cells,
            top_fraction_of_cells,
            bottom_fraction_of_cells);

        // first extract from the vector of indicators the ones that correspond
        // to cells that we locally own
        Vector<Number>
        locally_owned_indicators (tria.n_locally_owned_active_cells());
        get_locally_owned_indicators (tria,
                                      criteria,
                                      locally_owned_indicators);

        MPI_Comm mpi_communicator = tria.get_communicator ();

        // figure out the global max and min of the indicators. we don't need it
        // here, but it's a collective communication call
        const std::pair<Number,Number> global_min_and_max
          = compute_global_min_and_max_at_root (locally_owned_indicators,
                                                mpi_communicator);


        double top_threshold, bottom_threshold;
        top_threshold =
          RefineAndCoarsenFixedNumber::
          compute_threshold (locally_owned_indicators,
                             global_min_and_max,
                             static_cast<unsigned int>
                             (adjusted_fractions.first *
                              tria.n_global_active_cells()),
                             mpi_communicator);

        // compute bottom threshold only if necessary. otherwise use a threshold
        // lower than the smallest value we have locally
        if (adjusted_fractions.second > 0)
          bottom_threshold =
            RefineAndCoarsenFixedNumber::
            compute_threshold (locally_owned_indicators,
                               global_min_and_max,
                               static_cast<unsigned int>
                               ((1-adjusted_fractions.second) *
                                tria.n_global_active_cells()),
                               mpi_communicator);
        else
          {
            bottom_threshold = *std::min_element (criteria.begin(),
                                                  criteria.end());
            bottom_threshold -= std::fabs(bottom_threshold);
          }

        // now refine the mesh
        mark_cells (tria, criteria, top_threshold, bottom_threshold);
      }


      template <int dim, typename Number, int spacedim>
      void
      refine_and_coarsen_fixed_fraction
      (parallel::distributed::Triangulation<dim,spacedim> &tria,
       const dealii::Vector<Number>                               &criteria,
       const double                                        top_fraction_of_error,
       const double                                        bottom_fraction_of_error)
      {
        Assert (criteria.size() == tria.n_active_cells(),
                ExcDimensionMismatch (criteria.size(), tria.n_active_cells()));
        Assert ((top_fraction_of_error>=0) && (top_fraction_of_error<=1),
                dealii::GridRefinement::ExcInvalidParameterValue());
        Assert ((bottom_fraction_of_error>=0) && (bottom_fraction_of_error<=1),
                dealii::GridRefinement::ExcInvalidParameterValue());
        Assert (top_fraction_of_error+bottom_fraction_of_error <= 1,
                dealii::GridRefinement::ExcInvalidParameterValue());
        Assert (criteria.is_non_negative (),
                dealii::GridRefinement::ExcNegativeCriteria());

        // first extract from the vector of indicators the ones that correspond
        // to cells that we locally own
        Vector<Number> locally_owned_indicators (tria.n_locally_owned_active_cells());
        get_locally_owned_indicators (tria,
                                      criteria,
                                      locally_owned_indicators);

        MPI_Comm mpi_communicator = tria.get_communicator ();

        // figure out the global max and min of the indicators. we don't need it
        // here, but it's a collective communication call
        const std::pair<double,double> global_min_and_max
          = compute_global_min_and_max_at_root (locally_owned_indicators,
                                                mpi_communicator);

        const double total_error
          = compute_global_sum (locally_owned_indicators,
                                mpi_communicator);
        double top_threshold, bottom_threshold;
        top_threshold =
          RefineAndCoarsenFixedFraction::
          compute_threshold (locally_owned_indicators,
                             global_min_and_max,
                             top_fraction_of_error *
                             total_error,
                             mpi_communicator);
        // compute bottom threshold only if necessary. otherwise use a threshold
        // lower than the smallest value we have locally
        if (bottom_fraction_of_error > 0)
          bottom_threshold =
            RefineAndCoarsenFixedFraction::
            compute_threshold (locally_owned_indicators,
                               global_min_and_max,
                               (1-bottom_fraction_of_error) *
                               total_error,
                               mpi_communicator);
        else
          {
            bottom_threshold = *std::min_element (criteria.begin(),
                                                  criteria.end());
            bottom_threshold -= std::fabs(bottom_threshold);
          }

        // now refine the mesh
        mark_cells (tria, criteria, top_threshold, bottom_threshold);
      }
    }
  }
}


// explicit instantiations
#include "grid_refinement.inst"

DEAL_II_NAMESPACE_CLOSE

#endif
