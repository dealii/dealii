// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_distributed_grid_refinement_h
#define dealii_distributed_grid_refinement_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/numerics/vector_tools_common.h>

#include <limits>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace parallel
  {
    namespace distributed
    {
      namespace GridRefinement
      {
        /**
         * Compute the global max and min of the criteria vector. These are
         * returned only on the processor with rank zero, all others get a pair
         * of zeros.
         */
        template <typename number>
        std::pair<number, number>
        compute_global_min_and_max_at_root(
          const dealii::Vector<number> &criteria,
          const MPI_Comm                mpi_communicator);

        namespace RefineAndCoarsenFixedNumber
        {
          /**
           * Compute a threshold value so that exactly n_target_cells have a
           * value that is larger.
           */
          template <typename number>
          number
          compute_threshold(const dealii::Vector<number>    &criteria,
                            const std::pair<double, double> &global_min_and_max,
                            const types::global_cell_index   n_target_cells,
                            const MPI_Comm                   mpi_communicator);
        } // namespace RefineAndCoarsenFixedNumber

        namespace RefineAndCoarsenFixedFraction
        {
          /**
           * Compute a threshold value so that the error accumulated over all
           * criteria[i] so that
           * @code
           *     criteria[i] > threshold
           * @endcode
           * is larger than target_error.
           */
          template <typename number>
          number
          compute_threshold(const dealii::Vector<number>    &criteria,
                            const std::pair<double, double> &global_min_and_max,
                            const double                     target_error,
                            const MPI_Comm                   mpi_communicator);
        } // namespace RefineAndCoarsenFixedFraction
      }   // namespace GridRefinement
    }     // namespace distributed
  }       // namespace parallel
} // namespace internal



namespace parallel
{
  namespace distributed
  {
    /**
     * This namespace provides a collection of functions that aid in refinement
     * and coarsening of triangulations. Despite the name of the namespace, the
     * functions do not actually <i>refine</i> the triangulation, but only
     * <i>mark cells for refinement or coarsening</i>. In other words, they
     * perform the "mark" part of the typical "solve-estimate-mark-refine"
     * cycle of the adaptive finite element loop.
     *
     * In contrast to the functions in namespace dealii::GridRefinement,
     * the functions in the current namespace are intended for parallel
     * computations, i.e., computations, e.g., on objects of type
     * parallel::distributed::Triangulation.
     *
     * @ingroup grid
     */
    namespace GridRefinement
    {
      /**
       * Like dealii::GridRefinement::refine_and_coarsen_fixed_number, but
       * designed for parallel computations, where each process has only
       * information about locally owned cells.
       *
       * The vector of criteria needs to be a vector of refinement criteria
       * for all cells active on the current triangulation, i.e.,
       * it needs to be of length <code>tria.n_active_cells()</code> (and not
       * <code>tria.n_locally_owned_active_cells()</code>). In other words,
       * the vector needs to include entries for ghost and artificial
       * cells. However, the current
       * function will only look at the indicators that correspond to those
       * cells that are actually locally owned, and ignore the indicators for
       * all other cells. The function will then coordinate among all
       * processors that store part of the triangulation so that at the end
       * a fraction @p top_fraction_of_cells of all Triangulation::n_global_active_cells()
       * active cells are refined, rather than a fraction of the
       * Triangulation::n_locally_active_cells on each processor individually.
       * In other words, it may be that on some processors, no cells are
       * refined at all.
       *
       * The same is true for the fraction of cells that is coarsened.
       *
       * @param[in,out] tria The triangulation whose cells this function is
       * supposed to mark for coarsening and refinement.
       *
       * @param[in] criteria The refinement criterion for each mesh cell active
       * on the current triangulation. Entries may not be negative.
       *
       * @param[in] top_fraction_of_cells The fraction of cells to be refined.
       * If this number is zero, no cells will be refined. If it equals one,
       * the result will be flagging for global refinement.
       *
       * @param[in] bottom_fraction_of_cells The fraction of cells to be
       * coarsened. If this number is zero, no cells will be coarsened.
       *
       * @param[in] max_n_cells This argument can be used to specify a maximal
       * number of cells. If this number is going to be exceeded upon
       * refinement, then refinement and coarsening fractions are going to be
       * adjusted in an attempt to reach the maximum number of cells. Be aware
       * though that through proliferation of refinement due to
       * Triangulation::MeshSmoothing, this number is only an indicator. The
       * default value of this argument is to impose no limit on the number of
       * cells.
       */
      template <int dim, typename Number, int spacedim>
      void
      refine_and_coarsen_fixed_number(
        dealii::Triangulation<dim, spacedim> &tria,
        const dealii::Vector<Number>         &criteria,
        const double                          top_fraction_of_cells,
        const double                          bottom_fraction_of_cells,
        const types::global_cell_index        max_n_cells =
          std::numeric_limits<types::global_cell_index>::max());

      /**
       * Like dealii::GridRefinement::refine_and_coarsen_fixed_fraction, but
       * designed for parallel computations, where each process only has
       * information about locally owned cells.
       *
       * The vector of criteria needs to be a vector of refinement criteria
       * for all cells active on the current triangulation, i.e.,
       * it needs to be of length <code>tria.n_active_cells()</code> (and not
       * <code>tria.n_locally_owned_active_cells()</code>). In other words,
       * the vector needs to include entries for ghost and artificial
       * cells. However, the current
       * function will only look at the indicators that correspond to those
       * cells that are actually locally owned, and ignore the indicators for
       * all other cells. The function will then coordinate among all
       * processors that store part of the triangulation so that at the end
       * the smallest fraction of Triangulation::n_global_active_cells (not
       * Triangulation::n_locally_owned_active_cells() on each processor
       * individually)
       * is refined that together make up a total of @p top_fraction_of_error
       * of the total error. In other words, it may be that on some
       * processors, no cells are refined at all.
       *
       * The same is true for the fraction of cells that is coarsened.
       *
       * @param[in,out] tria The triangulation whose cells this function is
       * supposed to mark for coarsening and refinement.
       *
       * @param[in] criteria The refinement criterion computed on each mesh cell
       * active on the current triangulation. Entries may not be negative.
       *
       * @param[in] top_fraction_of_error The fraction of the total estimate
       * which should be refined. If this number is zero, no cells will be
       * refined. If it equals one, the result will be flagging for global
       * refinement.
       *
       * @param[in] bottom_fraction_of_error The fraction of the estimate
       * coarsened. If this number is zero, no cells will be coarsened.
       *
       * @param[in] norm_type To determine thresholds, combined errors on
       * subsets of cells are calculated as norms of the criteria on these
       * cells. Different types of norms can be used for this purpose, from
       * which VectorTools::L1_norm and
       * VectorTools::L2_norm are currently supported.
       */
      template <int dim, typename Number, int spacedim>
      void
      refine_and_coarsen_fixed_fraction(
        dealii::Triangulation<dim, spacedim> &tria,
        const dealii::Vector<Number>         &criteria,
        const double                          top_fraction_of_error,
        const double                          bottom_fraction_of_error,
        const VectorTools::NormType           norm_type = VectorTools::L1_norm);
    } // namespace GridRefinement
  }   // namespace distributed
} // namespace parallel


DEAL_II_NAMESPACE_CLOSE

#endif // dealii_distributed_grid_refinement_h
