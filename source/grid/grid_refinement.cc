// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2015 by the deal.II authors
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


#include <deal.II/base/template_constraints.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector_base.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>

#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria.h>

#include <numeric>
#include <algorithm>
#include <cmath>
#include <functional>
#include <fstream>

DEAL_II_NAMESPACE_OPEN


namespace
{
  namespace internal
  {
    template <typename number>
    inline
    number
    max_element (const Vector<number> &criteria)
    {
      return *std::max_element(criteria.begin(), criteria.end());
    }


    template <typename number>
    inline
    number
    min_element (const Vector<number> &criteria)
    {
      return *std::min_element(criteria.begin(), criteria.end());
    }

    // Silence a (bogus) warning in clang-3.6 about the following four
    // functions being unused:
    DEAL_II_DISABLE_EXTRA_DIAGNOSTICS

#ifdef DEAL_II_WITH_PETSC
    inline
    PetscScalar
    max_element (const PETScWrappers::Vector &criteria)
    {
      // this is horribly slow (since we have
      // to get the array of values from PETSc
      // in every iteration), but works
      PetscScalar m = 0;
      for (unsigned int i=0; i<criteria.size(); ++i)
        m = std::max (m, criteria(i));
      return m;
    }


    inline
    PetscScalar
    min_element (const PETScWrappers::Vector &criteria)
    {
      // this is horribly slow (since we have
      // to get the array of values from PETSc
      // in every iteration), but works
      PetscScalar m = criteria(0);
      for (unsigned int i=1; i<criteria.size(); ++i)
        m = std::min (m, criteria(i));
      return m;
    }
#endif


#ifdef DEAL_II_WITH_TRILINOS
    inline
    TrilinosScalar
    max_element (const TrilinosWrappers::Vector &criteria)
    {
      TrilinosScalar m = 0;
      criteria.trilinos_vector().MaxValue(&m);
      return m;
    }


    inline
    TrilinosScalar
    min_element (const TrilinosWrappers::Vector &criteria)
    {
      TrilinosScalar m = 0;
      criteria.trilinos_vector().MinValue(&m);
      return m;
    }
#endif

    DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

  } /* namespace internal */


  template <typename Vector>
  typename constraint_and_return_value<!IsBlockVector<Vector>::value,
           typename Vector::value_type>::type
           min_element (const Vector &criteria)
  {
    return internal::min_element (criteria);
  }


  template <typename Vector>
  typename constraint_and_return_value<!IsBlockVector<Vector>::value,
           typename Vector::value_type>::type
           max_element (const Vector &criteria)
  {
    return internal::max_element (criteria);
  }


  template <typename Vector>
  typename constraint_and_return_value<IsBlockVector<Vector>::value,
           typename Vector::value_type>::type
           min_element (const Vector &criteria)
  {
    typename Vector::value_type t = internal::min_element(criteria.block(0));
    for (unsigned int b=1; b<criteria.n_blocks(); ++b)
      t = std::min (t, internal::min_element(criteria.block(b)));

    return t;
  }


  template <typename Vector>
  typename constraint_and_return_value<IsBlockVector<Vector>::value,
           typename Vector::value_type>::type
           max_element (const Vector &criteria)
  {
    typename Vector::value_type t = internal::max_element(criteria.block(0));
    for (unsigned int b=1; b<criteria.n_blocks(); ++b)
      t = std::max (t, internal::max_element(criteria.block(b)));

    return t;
  }

}


namespace
{
  /**
   * Sorts the vector @p ind as an index vector of @p a in increasing order.
   * This implementation of quicksort seems to be faster than the standard
   * library version and is needed in @p refine_and_coarsen_optimize.
   */

  template <class Vector>
  void qsort_index (const Vector              &a,
                    std::vector<unsigned int> &ind,
                    int                        l,
                    int                        r)
  {
    int i,j;
    typename Vector::value_type v;

    if (r<=l)
      return;

    v = a(ind[r]);
    i = l-1;
    j = r;
    do
      {
        do
          {
            ++i;
          }
        while ((a(ind[i])>v) && (i<r));
        do
          {
            --j;
          }
        while ((a(ind[j])<v) && (j>0));

        if (i<j)
          std::swap (ind[i], ind[j]);
        else
          std::swap (ind[i], ind[r]);
      }
    while (i<j);
    qsort_index(a,ind,l,i-1);
    qsort_index(a,ind,i+1,r);
  }
}




template <int dim, class Vector, int spacedim>
void GridRefinement::refine (Triangulation<dim,spacedim> &tria,
                             const Vector       &criteria,
                             const double        threshold,
                             const unsigned int max_to_mark)
{
  Assert (criteria.size() == tria.n_active_cells(),
          ExcDimensionMismatch(criteria.size(), tria.n_active_cells()));
  Assert (criteria.is_non_negative (), ExcNegativeCriteria());

  // when all indicators are zero we
  // do not need to refine but only
  // to coarsen
  if (criteria.all_zero())
    return;

  const unsigned int n_cells = criteria.size();

//TODO: This is undocumented, looks fishy and seems unnecessary

  double new_threshold=threshold;
  // when threshold==0 find the
  // smallest value in criteria
  // greater 0
  if (new_threshold==0)
    {
      new_threshold = criteria(0);
      for (unsigned int index=1; index<n_cells; ++index)
        if (criteria(index)>0
            && (criteria(index)<new_threshold))
          new_threshold=criteria(index);
    }

  unsigned int marked=0;
  for (typename Triangulation<dim,spacedim>::active_cell_iterator cell = tria.begin_active();
       cell != tria.end(); ++cell)
    if (std::fabs(criteria(cell->active_cell_index())) >= new_threshold)
      {
        if (max_to_mark!=numbers::invalid_unsigned_int && marked>=max_to_mark)
          break;
        marked++;
        cell->set_refine_flag();
      }
}



template <int dim, class Vector, int spacedim>
void GridRefinement::coarsen (Triangulation<dim,spacedim> &tria,
                              const Vector       &criteria,
                              const double        threshold)
{
  Assert (criteria.size() == tria.n_active_cells(),
          ExcDimensionMismatch(criteria.size(), tria.n_active_cells()));
  Assert (criteria.is_non_negative (), ExcNegativeCriteria());

  for (typename Triangulation<dim,spacedim>::active_cell_iterator cell = tria.begin_active();
       cell != tria.end(); ++cell)
    if (std::fabs(criteria(cell->active_cell_index())) <= threshold)
      if (!cell->refine_flag_set())
        cell->set_coarsen_flag();
}



template <int dim, class Vector, int spacedim>
void
GridRefinement::refine_and_coarsen_fixed_number (Triangulation<dim,spacedim> &tria,
                                                 const Vector       &criteria,
                                                 const double        top_fraction,
                                                 const double        bottom_fraction,
                                                 const unsigned int  max_n_cells,
                                                 const double        refine_priority)
{
  // correct number of cells is
  // checked in @p{refine}
  Assert ((top_fraction>=0) && (top_fraction<=1), ExcInvalidParameterValue());
  Assert ((bottom_fraction>=0) && (bottom_fraction<=1), ExcInvalidParameterValue());
  Assert (top_fraction+bottom_fraction <= 1, ExcInvalidParameterValue());
  Assert (criteria.is_non_negative (), ExcNegativeCriteria());
  Assert(refine_priority >= 0.0, ExcInvalidParameterValue());
  Assert(refine_priority <= 1.0, ExcInvalidParameterValue());

  // Convert variables to double to avoid annoying type casting and round error.
  const double current_cell_number = static_cast<double>(criteria.size());
  const double max_cell_number = static_cast<double>(max_n_cells);

  // as we have no information on cells being refined isotropically or
  // anisotropically, assume isotropic refinement here, though that may
  // result in a worse approximation
  const double cell_increase_on_refine  = GeometryInfo<dim>::max_children_per_cell - 1.0;
  const double cell_decrease_on_coarsen = 1.0 - 1.0/GeometryInfo<dim>::max_children_per_cell;

  // first we estimate the cell number after refinement and coarsening along the
  // the requested fractions.
  const double refine_cells_requested  = top_fraction * current_cell_number;
  const double coarsen_cells_requested = bottom_fraction * current_cell_number;
  const double n_cell_after_adaptation = current_cell_number +
                                         refine_cells_requested * cell_increase_on_refine -
                                         coarsen_cells_requested * cell_decrease_on_coarsen;

  int refine_cells  = static_cast<int> (refine_cells_requested);
  int coarsen_cells = static_cast<int> (coarsen_cells_requested);

  if (n_cell_after_adaptation > max_cell_number)
    {
      // Let's clarify the simple math behind this situation.
      // Define symbols:
      // N     : @var current_cell_number
      // N_n   : @var n_cell_after_adaptation
      // N_l   : @var max_cell_number
      // N_r   : number of cell to refine
      // N_c   : number of cell to coarsen
      // N_r_r : @var refine_cells_requested
      // N_c_r : @var coarsen_cells_requested
      // C_r   : @var cell_increase_on_refine
      // C_c   : @var cell_decrease_on_coarsen
      //
      // Now, the problem is:
      //       N_n = N + N_r_r*C_r - N_c_r*C_c > N_l   ...(0)
      // And we want to bring the number of cell after adaptation down to its
      // upper limit:
      //       N_n = N + N_r*C_r - N_c*C_c = N_l       ...(1)
      //
      // It is easy to see that generally solution to equation (1) is not unique,
      // so there is a feasible region for N_r and N_c.
      // Note that equation (1) is linear on both N_r and N_c and larger N_r
      // corresponds to larger N_c. So we just have to find out the two
      // extreme points of the feasible region, and interpolate between them
      // using parameter @refine_priority to get the final decision.
      //
      // In order to find out the extreme values, let's check the constraints
      // that N_r and N_c shoud subject to:
      // 1) if we don't want to break both of the input frations, we need
      //         N_r <= N_r_r,
      //       N_c_r <= N_c,
      // 2) naturally
      //                 0 <= N_r,
      //               N_c <= N
      //         N_r + N_c <= N.
      // Thus finally
      //          0 <= N_r <= N_r_r,
      //      N_c_r <= N_c <= N - N_r.

      // Try lower extreme of N_r at minimum N_c, and limit the result from below
      // with its lower constraint. Limiting the result from above is not necessary
      // at this point since we are already inside the if block.
      const double refine_cells_min = std::max(
                                        (max_cell_number - current_cell_number + coarsen_cells_requested*
                                         cell_decrease_on_coarsen)/cell_increase_on_refine
                                        ,
                                        0.0);
      // Lower extreme of N_c is explicit.
      const double coarsen_cells_min = coarsen_cells_requested;


      // Try upper extreme of N_r at maximum N_c = N - N_r, and limit the result
      // from both above and below.
      const double refine_cells_max = std::max (
                                        std::min(
                                          (max_cell_number - (1.0-cell_decrease_on_coarsen) *
                                           current_cell_number)/(cell_increase_on_refine+cell_decrease_on_coarsen)
                                          ,
                                          refine_cells_requested)
                                        ,
                                        0.0);
      // Solve equation (1) with N_r = max(N_r) for upper extreme of N_c, and limit the result
      // from above.
      const double coarsen_cells_max = std::min(
                                         (current_cell_number + refine_cells_max * cell_increase_on_refine -
                                          max_cell_number)/cell_decrease_on_coarsen
                                         ,
                                         current_cell_number - refine_cells_max);

      // Interpolate between the two extremes
      const double refine_cells_actual =
        (1.0 - refine_priority) * refine_cells_min
        +      refine_priority  * refine_cells_max;
      const double coarsen_cells_actual =
        (1.0 - refine_priority) * coarsen_cells_min
        +      refine_priority  * coarsen_cells_max;

      refine_cells  = static_cast<int> (refine_cells_actual);
      coarsen_cells = static_cast<int> (coarsen_cells_actual);
    }

  if (refine_cells || coarsen_cells)
    {
      dealii::Vector<typename Vector::value_type> tmp (criteria);
      if (refine_cells)
        {
          std::nth_element (tmp.begin(), tmp.begin()+refine_cells,
                            tmp.end(),
                            std::greater<double>());
          refine (tria, criteria, *(tmp.begin() + refine_cells));
        }

      if (coarsen_cells)
        {
          std::nth_element (tmp.begin(), tmp.begin()+tmp.size()-coarsen_cells,
                            tmp.end(),
                            std::greater<double>());
          coarsen (tria, criteria,
                   *(tmp.begin() + tmp.size() - coarsen_cells));
        }
    }
}



template <int dim, class Vector, int spacedim>
void
GridRefinement::refine_and_coarsen_fixed_fraction (Triangulation<dim,spacedim> &tria,
                                                   const Vector       &criteria,
                                                   const double        top_fraction,
                                                   const double        bottom_fraction,
                                                   const unsigned int  max_n_cells)
{
  // correct number of cells is
  // checked in @p{refine}
  Assert ((top_fraction>=0) && (top_fraction<=1), ExcInvalidParameterValue());
  Assert ((bottom_fraction>=0) && (bottom_fraction<=1), ExcInvalidParameterValue());
  Assert (top_fraction+bottom_fraction <= 1, ExcInvalidParameterValue());
  Assert (criteria.is_non_negative (), ExcNegativeCriteria());

  // let tmp be the cellwise square of the
  // error, which is what we have to sum
  // up and compare with
  // @p{fraction_of_error*total_error}.
  dealii::Vector<typename Vector::value_type> tmp;
  tmp = criteria;
  const double total_error = tmp.l1_norm();

  // sort the largest criteria to the
  // beginning of the vector
  std::sort (tmp.begin(), tmp.end(), std::greater<double>());

  // compute thresholds
  typename dealii::Vector<typename Vector::value_type>::const_iterator
  pp=tmp.begin();
  for (double sum=0;
       (sum<top_fraction*total_error) && (pp!=(tmp.end()-1));
       ++pp)
    sum += *pp;
  double top_threshold = ( pp != tmp.begin () ?
                           (*pp+*(pp-1))/2 :
                           *pp );
  typename dealii::Vector<typename Vector::value_type>::const_iterator
  qq=(tmp.end()-1);
  for (double sum=0;
       (sum<bottom_fraction*total_error) && (qq!=tmp.begin());
       --qq)
    sum += *qq;
  double bottom_threshold = ( qq != (tmp.end()-1) ?
                              (*qq + *(qq+1))/2 :
                              0);

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
  {
    const unsigned int refine_cells  = pp - tmp.begin(),
                       coarsen_cells = tmp.end() - qq;

    if (static_cast<unsigned int>
        (tria.n_active_cells()
         + refine_cells * (GeometryInfo<dim>::max_children_per_cell - 1)
         - (coarsen_cells *
            (GeometryInfo<dim>::max_children_per_cell - 1) /
            GeometryInfo<dim>::max_children_per_cell))
        >
        max_n_cells)
      {
        refine_and_coarsen_fixed_number (tria,
                                         criteria,
                                         1.*refine_cells/criteria.size(),
                                         1.*coarsen_cells/criteria.size(),
                                         max_n_cells);
        return;
      }
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
  if ((top_threshold == max_element(criteria)) &&
      (top_fraction != 1))
    top_threshold *= 0.999;

  if (bottom_threshold>=top_threshold)
    bottom_threshold = 0.999*top_threshold;

  // actually flag cells
  if (top_threshold < max_element(criteria))
    refine (tria, criteria, top_threshold, pp - tmp.begin());

  if (bottom_threshold > min_element(criteria))
    coarsen (tria, criteria, bottom_threshold);
}



template <int dim, class Vector, int spacedim>
void
GridRefinement::refine_and_coarsen_optimize (Triangulation<dim,spacedim> &tria,
                                             const Vector       &criteria,
                                             const unsigned int  order)
{
  Assert (criteria.size() == tria.n_active_cells(),
          ExcDimensionMismatch(criteria.size(), tria.n_active_cells()));
  Assert (criteria.is_non_negative (), ExcNegativeCriteria());

  // get an increasing order on
  // the error indicator
  std::vector<unsigned int> tmp(criteria.size());
  for (unsigned int i=0; i<criteria.size(); ++i)
    tmp[i] = i;

  qsort_index (criteria, tmp, 0, criteria.size()-1);

  double expected_error_reduction = 0;
  const double original_error     = criteria.l1_norm();

  const unsigned int N = criteria.size();

  // minimize the cost functional discussed in the documentation
  double min_cost = std::numeric_limits<double>::max();
  unsigned int min_arg = 0;

  for (unsigned int M = 0; M<criteria.size(); ++M)
    {
      expected_error_reduction += (1-std::pow(2.,-1.*order)) * criteria(tmp[M]);

      const double cost = std::pow(((std::pow(2.,dim)-1)*(1+M)+N),
                                   (double)order/dim) *
                          (original_error-expected_error_reduction);
      if (cost <= min_cost)
        {
          min_cost = cost;
          min_arg = M;
        }
    }

  refine (tria, criteria, criteria(tmp[min_arg]));
}


// explicit instantiations
#include "grid_refinement.inst"

DEAL_II_NAMESPACE_CLOSE

