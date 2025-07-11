// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2000 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/template_constraints.h>

#include <deal.II/distributed/tria_base.h>

#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <limits>
#include <numeric>

DEAL_II_NAMESPACE_OPEN

namespace
{
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
    Triangulation<dim, spacedim> &tria,
    const Vector<Number>         &criteria,
    const double                  top_fraction,
    const double                  bottom_fraction,
    const unsigned int            max_n_cells)
  {
    // sort the criteria in descending order in an auxiliary vector, which we
    // have to sum up and compare with @p{fraction_of_error*total_error}
    Vector<Number> criteria_sorted = criteria;
    std::sort(criteria_sorted.begin(),
              criteria_sorted.end(),
              std::greater<double>());

    const double total_error = criteria_sorted.l1_norm();

    // compute thresholds
    typename Vector<Number>::const_iterator pp = criteria_sorted.begin();
    for (double sum = 0;
         (sum < top_fraction * total_error) && (pp != criteria_sorted.end());
         ++pp)
      sum += *pp;
    double top_threshold =
      (pp != criteria_sorted.begin() ? (*pp + *(pp - 1)) / 2 : *pp);

    typename Vector<Number>::const_iterator qq = criteria_sorted.end() - 1;
    for (double sum = 0; (sum < bottom_fraction * total_error) &&
                         (qq != criteria_sorted.begin() - 1);
         --qq)
      sum += *qq;
    double bottom_threshold =
      ((qq != criteria_sorted.end() - 1) ? (*qq + *(qq + 1)) / 2 : 0.);

    // we now have an idea how many cells we
    // are going to refine and coarsen. we use
    // this information to see whether we are
    // over the limit and if so use a function
    // that knows how to deal with this
    // situation

    // note, that at this point, we have no
    // information about anisotropically refined
    // cells, thus use the situation of purely
    // isotropic refinement as guess for a mixed
    // refinemnt as well.
    const unsigned int refine_cells  = pp - criteria_sorted.begin(),
                       coarsen_cells = criteria_sorted.end() - 1 - qq;

    if (static_cast<unsigned int>(
          tria.n_active_cells() +
          refine_cells * (GeometryInfo<dim>::max_children_per_cell - 1) -
          (coarsen_cells * (GeometryInfo<dim>::max_children_per_cell - 1) /
           GeometryInfo<dim>::max_children_per_cell)) > max_n_cells)
      {
        GridRefinement::refine_and_coarsen_fixed_number(tria,
                                                        criteria,
                                                        1. * refine_cells /
                                                          criteria.size(),
                                                        1. * coarsen_cells /
                                                          criteria.size(),
                                                        max_n_cells);
        return;
      }

    // in some rare cases it may happen that
    // both thresholds are the same (e.g. if
    // there are many cells with the same
    // error indicator). That would mean that
    // all cells will be flagged for
    // refinement or coarsening, but some will
    // be flagged for both, namely those for
    // which the indicator equals the
    // thresholds. This is forbidden, however.
    //
    // In some rare cases with very few cells
    // we also could get integer round off
    // errors and get problems with
    // the top and bottom fractions.
    //
    // In these case we arbitrarily reduce the
    // bottom threshold by one permille below
    // the top threshold
    //
    // Finally, in some cases
    // (especially involving symmetric
    // solutions) there are many cells
    // with the same error indicator
    // values. if there are many with
    // indicator equal to the top
    // threshold, no refinement will
    // take place below; to avoid this
    // case, we also lower the top
    // threshold if it equals the
    // largest indicator and the
    // top_fraction!=1
    const double max_criterion = *(criteria_sorted.begin()),
                 min_criterion = *(criteria_sorted.end() - 1);

    if ((top_threshold == max_criterion) && (top_fraction != 1))
      top_threshold *= 0.999;

    if (bottom_threshold >= top_threshold)
      bottom_threshold = 0.999 * top_threshold;

    // actually flag cells
    if (top_threshold < max_criterion)
      GridRefinement::refine(tria, criteria, top_threshold, refine_cells);

    if (bottom_threshold > min_criterion)
      GridRefinement::coarsen(tria, criteria, bottom_threshold);
  }
} // namespace



template <int dim, typename Number, int spacedim>
void
GridRefinement::refine(Triangulation<dim, spacedim> &tria,
                       const Vector<Number>         &criteria,
                       const double                  threshold,
                       const unsigned int            max_to_mark)
{
  Assert(criteria.size() == tria.n_active_cells(),
         ExcDimensionMismatch(criteria.size(), tria.n_active_cells()));
  Assert(criteria.is_non_negative(), ExcNegativeCriteria());

  // when all indicators are zero we
  // do not need to refine but only
  // to coarsen
  if (criteria.all_zero())
    return;

  const unsigned int n_cells = criteria.size();

  // TODO: This is undocumented, looks fishy and seems unnecessary

  double new_threshold = threshold;
  // when threshold==0 find the
  // smallest value in criteria
  // greater 0
  if (new_threshold == 0)
    {
      new_threshold = criteria(0);
      for (unsigned int index = 1; index < n_cells; ++index)
        if (criteria(index) > 0 && (criteria(index) < new_threshold))
          new_threshold = criteria(index);
    }

  unsigned int marked = 0;
  for (const auto &cell : tria.active_cell_iterators())
    if ((dynamic_cast<parallel::DistributedTriangulationBase<dim, spacedim> *>(
           &tria) == nullptr ||
         cell->is_locally_owned()) &&
        std::fabs(criteria(cell->active_cell_index())) >= new_threshold)
      {
        if (max_to_mark != numbers::invalid_unsigned_int &&
            marked >= max_to_mark)
          break;
        ++marked;
        cell->set_refine_flag();
      }
}



template <int dim, typename Number, int spacedim>
void
GridRefinement::coarsen(Triangulation<dim, spacedim> &tria,
                        const Vector<Number>         &criteria,
                        const double                  threshold)
{
  Assert(criteria.size() == tria.n_active_cells(),
         ExcDimensionMismatch(criteria.size(), tria.n_active_cells()));
  Assert(criteria.is_non_negative(), ExcNegativeCriteria());

  for (const auto &cell : tria.active_cell_iterators())
    if ((dynamic_cast<parallel::DistributedTriangulationBase<dim, spacedim> *>(
           &tria) == nullptr ||
         cell->is_locally_owned()) &&
        std::fabs(criteria(cell->active_cell_index())) <= threshold)
      if (!cell->refine_flag_set())
        cell->set_coarsen_flag();
}



template <int dim>
std::pair<double, double>
GridRefinement::adjust_refine_and_coarsen_number_fraction(
  const types::global_cell_index current_n_cells,
  const types::global_cell_index max_n_cells,
  const double                   top_fraction,
  const double                   bottom_fraction)
{
  Assert(top_fraction >= 0, ExcInvalidParameterValue());
  Assert(top_fraction <= 1, ExcInvalidParameterValue());
  Assert(bottom_fraction >= 0, ExcInvalidParameterValue());
  Assert(bottom_fraction <= 1, ExcInvalidParameterValue());
  Assert(top_fraction + bottom_fraction <=
           1 + 10 * std::numeric_limits<double>::epsilon(),
         ExcInvalidParameterValue());

  double refine_cells  = current_n_cells * top_fraction;
  double coarsen_cells = current_n_cells * bottom_fraction;

  const double cell_increase_on_refine =
    GeometryInfo<dim>::max_children_per_cell - 1.0;
  const double cell_decrease_on_coarsen =
    1.0 - 1.0 / GeometryInfo<dim>::max_children_per_cell;

  std::pair<double, double> adjusted_fractions(top_fraction, bottom_fraction);
  // first we have to see whether we
  // currently already exceed the target
  // number of cells
  if (current_n_cells >= max_n_cells)
    {
      // if yes, then we need to stop
      // refining cells and instead try to
      // only coarsen as many as it would
      // take to get to the target

      // as we have no information on cells
      // being refined isotropically or
      // anisotropically, assume isotropic
      // refinement here, though that may
      // result in a worse approximation
      adjusted_fractions.first = 0;
      coarsen_cells =
        (current_n_cells - max_n_cells) / cell_decrease_on_coarsen;
      adjusted_fractions.second =
        std::min(coarsen_cells / current_n_cells, 1.0);
    }
  // otherwise, see if we would exceed the
  // maximum desired number of cells with the
  // number of cells that are likely going to
  // result from refinement. here, each cell
  // to be refined is replaced by
  // C=GeometryInfo<dim>::max_children_per_cell
  // new cells, i.e. there will be C-1 more
  // cells than before. similarly, C cells
  // will be replaced by 1

  // again, this is true for isotropically
  // refined cells. we take this as an
  // approximation of a mixed refinement.
  else if (static_cast<types::global_cell_index>(
             current_n_cells + refine_cells * cell_increase_on_refine -
             coarsen_cells * cell_decrease_on_coarsen) > max_n_cells)
    {
      // we have to adjust the
      // fractions. assume we want
      // alpha*refine_fraction and
      // alpha*coarsen_fraction as new
      // fractions and the resulting number
      // of cells to be equal to
      // max_n_cells. this leads to the
      // following equation for alpha
      const double alpha = 1. * (max_n_cells - current_n_cells) /
                           (refine_cells * cell_increase_on_refine -
                            coarsen_cells * cell_decrease_on_coarsen);

      adjusted_fractions.first  = alpha * top_fraction;
      adjusted_fractions.second = alpha * bottom_fraction;
    }
  return (adjusted_fractions);
}



template <int dim, typename Number, int spacedim>
void
GridRefinement::refine_and_coarsen_fixed_number(
  Triangulation<dim, spacedim> &tria,
  const Vector<Number>         &criteria,
  const double                  top_fraction,
  const double                  bottom_fraction,
  const unsigned int            max_n_cells)
{
  // correct number of cells is
  // checked in @p{refine}
  Assert((top_fraction >= 0) && (top_fraction <= 1),
         ExcInvalidParameterValue());
  Assert((bottom_fraction >= 0) && (bottom_fraction <= 1),
         ExcInvalidParameterValue());
  Assert(top_fraction + bottom_fraction <=
           1 + 10 * std::numeric_limits<double>::epsilon(),
         ExcInvalidParameterValue());
  Assert(criteria.is_non_negative(), ExcNegativeCriteria());

  const std::pair<double, double> adjusted_fractions =
    adjust_refine_and_coarsen_number_fraction<dim>(criteria.size(),
                                                   max_n_cells,
                                                   top_fraction,
                                                   bottom_fraction);

  const int refine_cells =
    static_cast<int>(adjusted_fractions.first * criteria.size());
  const int coarsen_cells =
    static_cast<int>(adjusted_fractions.second * criteria.size());

  if (refine_cells || coarsen_cells)
    {
      Vector<Number> tmp(criteria);
      if (refine_cells)
        {
          if (static_cast<std::size_t>(refine_cells) == criteria.size())
            refine(tria, criteria, std::numeric_limits<double>::lowest());
          else
            {
              std::nth_element(tmp.begin(),
                               tmp.begin() + refine_cells - 1,
                               tmp.end(),
                               std::greater<double>());
              refine(tria, criteria, *(tmp.begin() + refine_cells - 1));
            }
        }

      if (coarsen_cells)
        {
          if (static_cast<std::size_t>(coarsen_cells) == criteria.size())
            coarsen(tria, criteria, std::numeric_limits<double>::max());
          else
            {
              std::nth_element(tmp.begin(),
                               tmp.begin() + tmp.size() - coarsen_cells,
                               tmp.end(),
                               std::greater<double>());
              coarsen(tria,
                      criteria,
                      *(tmp.begin() + tmp.size() - coarsen_cells));
            }
        }
    }
}



template <int dim, typename Number, int spacedim>
void
GridRefinement::refine_and_coarsen_fixed_fraction(
  Triangulation<dim, spacedim> &tria,
  const Vector<Number>         &criteria,
  const double                  top_fraction,
  const double                  bottom_fraction,
  const unsigned int            max_n_cells,
  const VectorTools::NormType   norm_type)
{
  // correct number of cells is checked in @p{refine}
  Assert((top_fraction >= 0) && (top_fraction <= 1),
         ExcInvalidParameterValue());
  Assert((bottom_fraction >= 0) && (bottom_fraction <= 1),
         ExcInvalidParameterValue());
  Assert(top_fraction + bottom_fraction <=
           1 + 10 * std::numeric_limits<double>::epsilon(),
         ExcInvalidParameterValue());
  Assert(criteria.is_non_negative(), ExcNegativeCriteria());

  switch (norm_type)
    {
      case VectorTools::L1_norm:
        // evaluate norms on subsets and compare them as
        //   c_0 + c_1 + ... < fraction * l1-norm(c)
        refine_and_coarsen_fixed_fraction_via_l1_norm(
          tria, criteria, top_fraction, bottom_fraction, max_n_cells);
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

          refine_and_coarsen_fixed_fraction_via_l1_norm(tria,
                                                        criteria_squared,
                                                        top_fraction *
                                                          top_fraction,
                                                        bottom_fraction *
                                                          bottom_fraction,
                                                        max_n_cells);
        }
        break;

      default:
        DEAL_II_NOT_IMPLEMENTED();
        break;
    }
}



template <int dim, typename Number, int spacedim>
void
GridRefinement::refine_and_coarsen_optimize(Triangulation<dim, spacedim> &tria,
                                            const Vector<Number> &criteria,
                                            const unsigned int    order)
{
  Assert(criteria.size() == tria.n_active_cells(),
         ExcDimensionMismatch(criteria.size(), tria.n_active_cells()));
  Assert(criteria.is_non_negative(), ExcNegativeCriteria());

  // get a decreasing order on the error indicator
  std::vector<unsigned int> cell_indices(criteria.size());
  std::iota(cell_indices.begin(), cell_indices.end(), 0u);

  std::sort(cell_indices.begin(),
            cell_indices.end(),
            [&criteria](const unsigned int left, const unsigned int right) {
              return criteria[left] > criteria[right];
            });

  double       expected_error_reduction = 0;
  const double original_error           = criteria.l1_norm();

  const std::size_t N = criteria.size();

  // minimize the cost functional discussed in the documentation
  double      min_cost = std::numeric_limits<double>::max();
  std::size_t min_arg  = 0;

  const double reduction_factor = (1. - std::pow(2., -1. * order));
  for (std::size_t M = 0; M < criteria.size(); ++M)
    {
      expected_error_reduction += reduction_factor * criteria(cell_indices[M]);

      const double cost =
        std::pow(((Utilities::fixed_power<dim>(2) - 1) * (1 + M) + N),
                 static_cast<double>(order) / dim) *
        (original_error - expected_error_reduction);
      if (cost <= min_cost)
        {
          min_cost = cost;
          min_arg  = M;
        }
    }

  refine(tria, criteria, criteria(cell_indices[min_arg]));
}


// explicit instantiations
#include "grid/grid_refinement.inst"

DEAL_II_NAMESPACE_CLOSE
