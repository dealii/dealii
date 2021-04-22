// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2020 by the deal.II authors
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

#ifndef dealii_grid_refinement_h
#define dealii_grid_refinement_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <deal.II/numerics/vector_tools_common.h>

#include <limits>

DEAL_II_NAMESPACE_OPEN

// forward declarations
#ifndef DOXYGEN
template <int dim, int spacedim>
class Triangulation;
template <typename Number>
class Vector;
#endif

/**
 * This namespace provides a collection of functions that aid in refinement
 * and coarsening of triangulations. Despite the name of the namespace, the
 * functions do not actually <i>refine</i> the triangulation, but only
 * <i>mark cells for refinement or coarsening</i>. In other words, they
 * perform the "mark" part of the typical "solve-estimate-mark-refine"
 * cycle of the adaptive finite element loop.
 *
 * The functions in this namespace form two categories. There are the
 * auxiliary functions refine() and coarsen(). More important for users are
 * the other functions, which implement refinement strategies, as being found
 * in the literature on adaptive finite element methods. For mathematical
 * discussion of these methods, consider works by D&ouml;rfler, Morin,
 * Nochetto, Rannacher, Stevenson, and others.
 *
 * @ingroup grid
 */
namespace GridRefinement
{
  /**
   * Return a pair of double values of which the first is adjusted refinement
   * fraction of cells and the second is adjusted coarsening fraction of
   * cells.
   *
   *
   * @param[in] current_n_cells The current cell number.
   *
   * @param[in] max_n_cells The maximal number of cells. If current cell
   * number @p current_n_cells is already exceeded maximal cell number @p
   * max_n_cells, refinement fraction of cells will be set to zero and
   * coarsening fraction of cells will be adjusted to reduce cell number to @
   * max_n_cells. If cell number is going to be exceeded only upon refinement,
   * then refinement and coarsening fractions are going to be adjusted with a
   * same ratio in an attempt to reach the maximum number of cells. Be aware
   * though that through proliferation of refinement due to
   * Triangulation::MeshSmoothing, this number is only an indicator. The
   * default value of this argument is to impose no limit on the number of
   * cells.
   *
   * @param[in] top_fraction_of_cells The requested fraction of active
   * cells to be refined.
   *
   * @param[in] bottom_fraction_of_cells The requested fraction of
   * active cells to be coarsened.
   *
   * @note Usually you do not need to call this function explicitly. Pass @p
   * max_n_cells to function refine_and_coarsen_fixed_number() or function
   * refine_and_coarsen_fixed_fraction() and they will call this function if
   * necessary.
   */
  template <int dim>
  std::pair<double, double>
  adjust_refine_and_coarsen_number_fraction(
    const types::global_cell_index current_n_cells,
    const types::global_cell_index max_n_cells,
    const double                   top_fraction_of_cells,
    const double                   bottom_fraction_of_cells);

  /**
   * This function provides a strategy to mark cells for refinement and
   * coarsening with the goal of providing predictable growth in
   * the size of the mesh by refining a given fraction of all cells.
   *
   * The function takes a vector of refinement @p criteria and two values
   * between zero and one denoting the fractions of cells to be refined and
   * coarsened. It flags cells for further processing by
   * Triangulation::execute_coarsening_and_refinement() according to the
   * following greedy algorithm:
   *
   * <ol>
   *
   * <li> Sort the cells according to descending values of @p criteria.
   *
   * <li> Mark the @p top_fraction_of_cells times
   * Triangulation::n_active_cells() active cells with the largest
   * refinement criteria for refinement.
   *
   * <li> Mark the @p bottom_fraction_of_cells times
   * Triangulation::n_active_cells() active cells with the smallest
   * refinement criteria for coarsening.
   *
   * </ol>
   *
   * As an example, with no coarsening, setting @p top_fraction_of_cells to
   * 1/3 will result in approximately doubling the number of cells in two
   * dimensions. That is because each of these 1/3 of cells will be replaced by
   * its four children, resulting in $4\times \frac 13 N$ cells, whereas the
   * remaining 2/3 of cells remains untouched -- thus yielding a total of
   * $4\times \frac 13 N + \frac 23 N = 2N$ cells.
   * The same effect in three dimensions is achieved by refining
   * 1/7th of the cells. These values are therefore frequently used because
   * they ensure that the cost of computations on subsequent meshes become
   * expensive sufficiently quickly that the fraction of time spent on
   * the coarse meshes is not too large. On the other hand, the fractions
   * are small enough that mesh adaptation does not refine too many cells
   * in each step.
   *
   * @note This function only sets the coarsening and refinement flags. The
   * mesh is not changed until you call
   * Triangulation::execute_coarsening_and_refinement().
   *
   * @param[in,out] triangulation The triangulation whose cells this function
   * is supposed to mark for coarsening and refinement.
   *
   * @param[in] criteria The refinement criterion for each mesh cell. Entries
   * may not be negative.
   *
   * @param[in] top_fraction_of_cells The fraction of cells to be refined. If
   * this number is zero, no cells will be refined. If it equals one, the
   * result will be flagging for global refinement.
   *
   * @param[in] bottom_fraction_of_cells The fraction of cells to be
   * coarsened. If this number is zero, no cells will be coarsened.
   *
   * @param[in] max_n_cells This argument can be used to specify a maximal
   * number of cells. If this number is going to be exceeded upon refinement,
   * then refinement and coarsening fractions are going to be adjusted in an
   * attempt to reach the maximum number of cells. Be aware though that
   * through proliferation of refinement due to Triangulation::MeshSmoothing,
   * this number is only an indicator. The default value of this argument is
   * to impose no limit on the number of cells.
   */
  template <int dim, typename Number, int spacedim>
  void
  refine_and_coarsen_fixed_number(
    Triangulation<dim, spacedim> &triangulation,
    const Vector<Number> &        criteria,
    const double                  top_fraction_of_cells,
    const double                  bottom_fraction_of_cells,
    const unsigned int max_n_cells = std::numeric_limits<unsigned int>::max());

  /**
   * This function provides a strategy to mark cells for refinement and
   * coarsening with the goal of controlling the reduction of
   * the error estimate.
   *
   * Also known as the <b>bulk criterion</b> or D&ouml;rfler marking,
   * this function computes the thresholds for refinement and coarsening
   * such that the @p criteria of cells getting flagged for refinement make
   * up for a certain fraction of the total error. We explain its operation
   * for refinement, coarsening works analogously.
   *
   * Let <i>c<sub>K</sub></i> be the criterion of cell <i>K</i>. Then the
   * total error estimate is computed by the formula
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
   * recovered by actually storing the square of <i>c<sub>K</sub></i> in the
   * vector @p criteria.
   *
   * From the point of view of implementation, this time we really need to
   * sort the array of criteria.  Just like the other strategy described
   * above, this function only computes the threshold values and then passes
   * over to refine() and coarsen().
   *
   * @param[in,out] tria The triangulation whose cells this function is
   * supposed to mark for coarsening and refinement.
   *
   * @param[in] criteria The refinement criterion computed on each mesh cell.
   * Entries may not be negative.
   *
   * @param[in] top_fraction The fraction of the total estimate which should
   * be refined. If this number is zero, no cells will be refined. If it
   * equals one, the result will be flagging for global refinement.
   *
   * @param[in] bottom_fraction The fraction of the estimate coarsened. If
   * this number is zero, no cells will be coarsened.
   *
   * @param[in] max_n_cells This argument can be used to specify a maximal
   * number of cells. If this number is going to be exceeded upon refinement,
   * then refinement and coarsening fractions are going to be adjusted in an
   * attempt to reach the maximum number of cells. Be aware though that
   * through proliferation of refinement due to Triangulation::MeshSmoothing,
   * this number is only an indicator. The default value of this argument is
   * to impose no limit on the number of cells.
   *
   * @param[in] norm_type To determine thresholds, combined errors on
   * subsets of cells are calculated as norms of the criteria on these
   * cells. Different types of norms can be used for this purpose, from
   * which VectorTools::NormType::L1_norm and
   * VectorTools::NormType::L2_norm are currently supported.
   */
  template <int dim, typename Number, int spacedim>
  void
  refine_and_coarsen_fixed_fraction(
    Triangulation<dim, spacedim> &tria,
    const Vector<Number> &        criteria,
    const double                  top_fraction,
    const double                  bottom_fraction,
    const unsigned int max_n_cells = std::numeric_limits<unsigned int>::max(),
    const VectorTools::NormType norm_type = VectorTools::NormType::L1_norm);



  /**
   * This function flags cells of a triangulation for refinement with the
   * aim to reach a grid that
   * is optimal with respect to an objective function that tries to balance
   * reducing the error and increasing the numerical cost when the mesh is
   * refined. Specifically, this function makes the assumption that if you
   * refine a cell $K$ with error indicator $\eta_K$ provided by the second
   * argument to this function, then the error on the children (for all
   * children together) will only be $2^{-\text{order}}\eta_K$ where
   * <code>order</code> is the third argument of this function. This makes the
   * assumption that the error is only a local property on a mesh and can be
   * reduced by local refinement -- an assumption that is true for the
   * interpolation operator, but not for the usual Galerkin projection,
   * although it is approximately true for elliptic problems where the Greens
   * function decays quickly and the error here is not too much affected by a
   * too coarse mesh somewhere else.
   *
   * With this, we can define the objective function this function tries to
   * optimize. Let us assume that the mesh currently has $N_0$ cells. Then, if
   * we refine the $m$ cells with the largest errors, we expect to get (in $d$
   * space dimensions)
   * @f[
   *   N(m) = (N_0-m) + 2^d m = N_0 + (2^d-1)m
   * @f]
   * cells ($N_0-m$ are not refined, and each of the $m$ cells we refine yield
   * $2^d$ child cells. On the other hand, with refining $m$ cells, and using
   * the assumptions above, we expect that the error will be
   * @f[
   *   \eta^\text{exp}(m)
   *   =
   *   \sum_{K, K\; \text{will not be refined}} \eta_K
   *   +
   *   \sum_{K, K\; \text{will be refined}} 2^{-\text{order}}\eta_K
   * @f]
   * where the first sum extends over $N_0-m$ cells and the second over the
   * $m$ cells that will be refined. Note that $N(m)$ is an increasing
   * function of $m$ whereas $\eta^\text{exp}(m)$ is a decreasing function.
   *
   * This function then tries to find that number $m$ of cells to mark for
   * refinement for which the objective function
   * @f[
   *   J(m) = N(m)^{\text{order}/d} \eta^\text{exp}(m)
   * @f]
   * is minimal.
   *
   * The rationale for this function is two-fold. First, compared to the
   * refine_and_coarsen_fixed_fraction() and refine_and_coarsen_fixed_number()
   * functions, this function has the property that if all refinement
   * indicators are the same (i.e., we have achieved a mesh where the error
   * per cell is equilibrated), then the entire mesh is refined. This is based
   * on the observation that a mesh with equilibrated error indicators is the
   * optimal mesh (i.e., has the least overall error) among all meshes with
   * the same number of cells. (For proofs of this, see R. Becker, M. Braack,
   * R. Rannacher: "Numerical simulation of laminar flames at low Mach number
   * with adaptive finite elements", Combustion Theory and Modelling, Vol. 3,
   * Nr. 3, p. 503-534 1999; and W. Bangerth, R. Rannacher: "Adaptive Finite
   * Element Methods for Differential Equations", Birkhauser, 2003.)
   *
   * Second, the function uses the observation that ideally, the error behaves
   * like $e \approx c N^{-\alpha}$ with some constant $\alpha$ that depends
   * on the dimension and the finite element degree. It should - given optimal
   * mesh refinement - not depend so much on the regularity of the solution,
   * as it is based on the idea, that all singularities can be resolved by
   * refinement. Mesh refinement is then based on the idea that we want to
   * make $c=e N^\alpha$ small. This corresponds to the functional $J(m)$
   * above.
   *
   * @note This function was originally implemented by Thomas Richter. It
   * follows a strategy described in T. Richter, "Parallel Multigrid Method
   * for Adaptive Finite Elements with Application to 3D Flow Problems", PhD
   * thesis, University of Heidelberg, 2005. See in particular Section 4.3,
   * pp. 42-43.
   */
  template <int dim, typename Number, int spacedim>
  void
  refine_and_coarsen_optimize(Triangulation<dim, spacedim> &tria,
                              const Vector<Number> &        criteria,
                              const unsigned int            order = 2);

  /**
   * Mark all mesh cells for which the value in @p criteria exceeds @p
   * threshold for refinement, but only flag up to @p max_to_mark cells.
   *
   * The vector @p criteria contains a nonnegative value for each active cell,
   * ordered in the canonical order of Triangulation::active_cell_iterator.
   *
   * The cells are only flagged for refinement, they are not actually refined.
   * To do so, you have to call
   * Triangulation::execute_coarsening_and_refinement().
   *
   * This function does not implement a refinement strategy, it is more a
   * helper function for the actual strategies.
   */
  template <int dim, typename Number, int spacedim>
  void
  refine(Triangulation<dim, spacedim> &tria,
         const Vector<Number> &        criteria,
         const double                  threshold,
         const unsigned int max_to_mark = numbers::invalid_unsigned_int);

  /**
   * Mark all mesh cells for which the value in @p criteria is less than @p
   * threshold for coarsening.
   *
   * The vector @p criteria contains a nonnegative value for each active cell,
   * ordered in the canonical order of Triangulation::active_cell_iterator.
   *
   * The cells are only flagged for coarsening, they are not actually
   * coarsened. To do so, you have to call
   * Triangulation::execute_coarsening_and_refinement().
   *
   * This function does not implement a refinement strategy, it is more a
   * helper function for the actual strategies.
   */
  template <int dim, typename Number, int spacedim>
  void
  coarsen(Triangulation<dim, spacedim> &tria,
          const Vector<Number> &        criteria,
          const double                  threshold);

  /**
   * An exception thrown if the vector with cell criteria contains negative
   * values
   */
  DeclException0(ExcNegativeCriteria);

  /**
   * One of the threshold parameters causes trouble. Or the refinement and
   * coarsening thresholds overlap.
   */
  DeclException0(ExcInvalidParameterValue);
} // namespace GridRefinement



DEAL_II_NAMESPACE_CLOSE

#endif // dealii_grid_refinement_h
