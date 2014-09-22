// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
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

#ifndef __deal2__grid_refinement_h
#define __deal2__grid_refinement_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/grid/tria.h>

#include <vector>
#include <limits>

DEAL_II_NAMESPACE_OPEN

// forward declarations
template <int dim, int spacedim> class Triangulation;
template <class T> class Vector;


/**
 * Collection of functions controlling refinement and coarsening of
 * Triangulation objects.
 *
 * The functions in this namespace are in two classes. There are the
 * auxiliary functions refine() and coarsen(). More important for
 * users are the other functions, which implement refinement
 * strategies, as being found in the literature on adaptive finite
 * element methods. For mathematical discussion of these methods,
 * consider works by D&ouml;rfler, Morin, Nochetto, Rannacher,
 * Stevenson and many more.
 *
 * @ingroup grid
 * @author Wolfgang Bangerth, Thomas Richter, Guido Kanschat 1998, 2000, 2009
 */
namespace GridRefinement
{
  /**
   * This function provides a refinement strategy with predictable growth of the mesh.
   *
   * The function takes a vector of refinement @p criteria and two values
   * between zero and one denoting the fractions of cells to be refined
   * and coarsened. It flags cells for further processing by
   * Triangulation::execute_coarsening_and_refinement() according to the
   * following greedy algorithm:
   *
   * <ol>
   *
   *  <li> Sort the cells according to descending values of @p criteria.
   *
   *  <li> Set the refinement threshold to be the criterion belonging to
   *  the cell at position @p top_fraction_of_cells times
   *  Triangulation::n_active_cells().
   *
   * <li> Set the coarsening threshold accordingly using the cell @p
   *  bottom_fraction_of_cells times Triangulation::n_active_cells()
   *  from the end of the sorted list.
   *
   *  <li> Use these two thresholds in calls to refine() and coarsen(),
   *  respectively.
   *
   * </ol>
   *
   * As an example, with no coarsening, setting @p top_fraction_of_cells
   * to 1/3 will result in approximately doubling the number of cells in
   * two dimensions. The same effect in three dimensions is achieved by
   * refining 1/7th of the cells. These values are good initial guesses,
   * but should be adjusted depending on the singularity of approximated
   * function.
   *
   * The sorting of criteria is not done actually, since we only need
   * the threshold values in order to call refine() and coarsen(). The
   * order of cells with higher and of those with lower criteria is
   * irrelevant. Getting this value is accomplished by the @p
   * nth_element function of the <tt>C++</tt> standard library, which
   * takes only linear time in the number of elements, rather than <tt>N
   * log N</tt> for sorting all values.
   *
   * @warning This function only sets the coarsening and refinement
   * flags. The mesh is not changed until you call
   * Triangulation::execute_coarsening_and_refinement().
   *
   * @arg @p criteria: the refinement criterion computed on each mesh
   * cell. Entries may not be negative.
   *
   * @arg @p top_fraction_of_cells is the fraction of cells to be
   * refined. If this number is zero, no cells will be refined. If it
   * equals one, the result will be flagging for global refinement.
   *
   * @arg @p bottom_fraction_of_cells is the fraction of cells to be
   * coarsened. If this number is zero, no cells will be coarsened.
   *
   * @arg @p max_n_cells can be used to specify a maximal number of
   * cells. If this number is going to be exceeded upon refinement, then
   * refinement and coarsening fractions are going to be adjusted in an
   * attempt to reach the maximum number of cells. Be aware though that
   * through proliferation of refinement due to
   * Triangulation::MeshSmoothing, this number is only an indicator. The
   * default value of this argument is to impose no limit on the number
   * of cells.
   */
  template <int dim, class Vector, int spacedim>
  void
  refine_and_coarsen_fixed_number (
    Triangulation<dim,spacedim> &tria,
    const Vector                &criteria,
    const double                top_fraction_of_cells,
    const double                bottom_fraction_of_cells,
    const unsigned int          max_n_cells = std::numeric_limits<unsigned int>::max());

  /**
   * This function provides a refinement strategy controlling the
   * reduction of the error estimate.
   *
   * Also known as the <b>bulk criterion</b>, this function computes the
   * thresholds for refinement and coarsening such that the @p criteria
   * of cells getting flagged for refinement make up for a certain
   * fraction of the total error. We explain its operation for
   * refinement, coarsening works analogously.
   *
   * Let <i>c<sub>K</sub></i> be the criterion of cell <i>K</i>. Then
   * the total error estimate is computed by the formula
   * @f[
   * E = \sum_{K\in \cal T} c_K.
   * @f]
   *
   * If <i> 0 &lt; a &lt; 1</i> is @p top_fraction, then we refine the
   * smallest subset $\cal M$ of the Triangulation $\cal T$ such that
   * @f[
   * a E \le \sum_{K\in \cal M} c_K
   * @f]
   *
   * The algorithm is performed by the greedy algorithm described in
   * refine_and_coarsen_fixed_number().
   *
   * @note The often used formula with squares on the left and right is
   * recovered by actually storing the square of <i>c<sub>K</sub></i> in
   * the vector @p criteria.
   *
   * From the point of view of implementation, this time we really need
   * to sort the array of criteria.  Just like the other strategy
   * described above, this function only computes the threshold values
   * and then passes over to refine() and coarsen().
   *
   * @arg @p criteria: the refinement criterion computed on each mesh
   * cell. Entries may not be negative.
   *
   * @arg @p top_fraction is the fraction of the total estimate which
   * should be refined. If this number is zero, no cells will be refined. If it
   * equals one, the result will be flagging for global refinement.
   *
   * @arg @p bottom_fraction is the fraction of the estimate
   * coarsened. If this number is zero, no cells will be coarsened.
   *
   * @arg @p max_n_cells can be used to specify a maximal number of
   * cells. If this number is going to be exceeded upon refinement, then
   * refinement and coarsening fractions are going to be adjusted in an
   * attempt to reach the maximum number of cells. Be aware though that
   * through proliferation of refinement due to
   * Triangulation::MeshSmoothing, this number is only an indicator. The
   * default value of this argument is to impose no limit on the number
   * of cells.
   */
  template <int dim, class Vector, int spacedim>
  void
  refine_and_coarsen_fixed_fraction (
    Triangulation<dim,spacedim> &tria,
    const Vector                &criteria,
    const double                top_fraction,
    const double                bottom_fraction,
    const unsigned int          max_n_cells = std::numeric_limits<unsigned int>::max());



  /**
   * Refine the triangulation by flagging
   * certain cells to reach an optimal
   * grid: We try to minimize the error
   * multiplied with the number of cells in
   * the new grid. All cells with large
   * error indicator are refined to
   * generate an optimal grid in the above
   * sense.  We assume that the error in
   * one cell is reduced to 1-2^{-order}
   * after refinement, if 'order' is the
   * expected order of convergence. This
   * expected order of convergence must be
   * passed as an argument but is defaulted
   * to 2.  The new triangulation has
   * ($2^d-1$) new cells for every flagged
   * cell (the original cell is replaced by
   * $2^d$ cells but it then made
   * inactive).
   *
   * Refer to the general doc of
   * this class for more
   * information.
   */
  template <int dim, class Vector, int spacedim>
  void
  refine_and_coarsen_optimize (Triangulation<dim,spacedim> &tria,
                               const Vector                &criteria,
                               const unsigned int           order=2);

  /**
   * Flag all mesh cells for which the value in @p criteria exceeds @p
   * threshold for refinement, but only flag up to @p max_to_mark cells.
   *
   * The vector @p criteria contains a nonnegative value for each active
   * cell, ordered in the canonical order of of
   * Triangulation::active_cell_iterator.
   *
   * The cells are only flagged for refinement, they are not actually
   * refined. To do so, you have to call
   * Triangulation::execute_coarsening_and_refinement().
   *
   * This function does not implement a refinement strategy, it is more
   * a helper function for the actual strategies.
   */
  template <int dim, class Vector, int spacedim>
  void refine (Triangulation<dim,spacedim> &tria,
               const Vector                &criteria,
               const double                threshold,
               const unsigned int max_to_mark = numbers::invalid_unsigned_int);

  /**
   * Flag all mesh cells for which the value in @p criteria
   * is less than @p threshold for coarsening.
   *
   * The vector @p criteria contains a nonnegative value for each active cell,
   * ordered in the canonical order of of
   * Triangulation::active_cell_iterator.
   *
   * The cells are only flagged for coarsening, they are not actually
   * coarsened. To do so, you have to call
   * Triangulation::execute_coarsening_and_refinement().
   *
   * This function does not implement a refinement strategy, it is more
   * a helper function for the actual strategies.
   */
  template <int dim, class Vector, int spacedim>
  void coarsen (Triangulation<dim,spacedim> &tria,
                const Vector                &criteria,
                const double                threshold);

  /**
   * An exception thrown if the
   * vector with cell criteria contains
   * negative values
   */
  DeclException0(ExcNegativeCriteria);

  /**
   * One of the threshold parameters
   * causes trouble. Or the
   * refinement and coarsening
   * thresholds overlap.
   */
  DeclException0 (ExcInvalidParameterValue);
}



DEAL_II_NAMESPACE_CLOSE

#endif //__deal2__grid_refinement_h
