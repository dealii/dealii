// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_adaptation_strategies_h
#define dealii_adaptation_strategies_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <deal.II/grid/tria_accessor.h>

#include <algorithm>
#include <numeric>
#include <typeinfo>
#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * When data is transferred during adaptation, it is not trivial to decide
 * how to process data from former cells on the old mesh that have been changed
 * into current cells on the new mesh. Or in other words, how data should be
 * stored in the cells on the adapted mesh.
 *
 * In this namespace, we offer a few strategies that cope with this
 * problem. Such strategies can be passed to the CellDataTransfer and
 * parallel::distributed::CellDataTransfer constructors.
 */
namespace AdaptationStrategies
{
  /**
   * For refinement, all strategies take the parent cell and its associated
   * data. They return a vector containing data for each individual child that
   * the parent cell will be refined to.
   *
   * The ordering of values in the vector for children data corresponds to the
   * index when calling TriaAccessor::child_index.
   */
  namespace Refinement
  {
    /**
     * Return a vector containing copies of data of the parent cell for each
     * child.
     *
     * @f[
     *   d_{K_c} = d_{K_p}
     *   \qquad
     *   \forall K_c \text{ children of } K_p
     * @f]
     */
    template <int dim, int spacedim, typename value_type>
    std::vector<value_type>
    preserve(const typename dealii::Triangulation<dim, spacedim>::cell_iterator
                             &parent,
             const value_type parent_value);

    /**
     * Return a vector which contains data of the parent cell being equally
     * divided among all children.
     *
     * @f[
     *   d_{K_c} = d_{K_p} / n_\text{children}
     *   \qquad
     *   \forall K_c \text{ children of } K_p
     * @f]
     *
     * This strategy preserves the $l_1$-norm of the corresponding global data
     * Vector before and after adaptation.
     */
    template <int dim, int spacedim, typename value_type>
    std::vector<value_type>
    split(const typename dealii::Triangulation<dim, spacedim>::cell_iterator
                          &parent,
          const value_type parent_value);

    /**
     * Return a vector which contains squared data of the parent cell being
     * equally divided among the squares of all children.
     *
     * @f[
     *   d_{K_c}^2 = d_{K_p}^2 / n_\text{children}
     *   \qquad
     *   \forall K_c \text{ children of } K_p
     * @f]
     *
     * This strategy preserves the $l_2$-norm of the corresponding global data
     * Vector before and after adaptation.
     */
    template <int dim, int spacedim, typename value_type>
    std::vector<value_type>
    l2_norm(const typename dealii::Triangulation<dim, spacedim>::cell_iterator
                            &parent,
            const value_type parent_value);
  } // namespace Refinement

  /**
   * For coarsening, all strategies take the parent cell and a vector of data
   * that belonged to its former children. They return the value that will be
   * assigned to the parent cell.
   *
   * The ordering of values in the vector for children data corresponds to the
   * index when calling TriaAccessor::child_index.
   */
  namespace Coarsening
  {
    /**
     * Check if data on all children match, and return value of the first child.
     *
     * @f[
     *   d_{K_p} = d_{K_c}
     *   \qquad
     *   \forall K_c \text{ children of } K_p
     * @f]
     */
    template <int dim, int spacedim, typename value_type>
    value_type
    check_equality(
      const typename dealii::Triangulation<dim, spacedim>::cell_iterator
                                    &parent,
      const std::vector<value_type> &children_values);

    /**
     * Return sum of data on all children.
     *
     * @f[
     *   d_{K_p} = \sum d_{K_c}
     *   \qquad
     *   \forall K_c \text{ children of } K_p
     * @f]
     *
     * This strategy preserves the $l_1$-norm of the corresponding global data
     * vector before and after adaptation.
     */
    template <int dim, int spacedim, typename value_type>
    value_type
    sum(const typename dealii::Triangulation<dim, spacedim>::cell_iterator
                                      &parent,
        const std::vector<value_type> &children_values);

    /**
     * Return $ l_2 $-norm of data on all children.
     *
     * @f[
     *   d_{K_p}^2 = \sum d_{K_c}^2
     *   \qquad
     *   \forall K_c \text{ children of } K_p
     * @f]
     *
     * This strategy preserves the $l_2$-norm of the corresponding global data
     * vector before and after adaptation.
     */
    template <int dim, int spacedim, typename value_type>
    value_type
    l2_norm(const typename dealii::Triangulation<dim, spacedim>::cell_iterator
                                          &parent,
            const std::vector<value_type> &children_values);

    /**
     * Return mean value of data on all children.
     *
     * @f[
     *   d_{K_p} = \sum d_{K_c} / n_\text{children}
     *   \qquad
     *   \forall K_c \text{ children of } K_p
     * @f]
     */
    template <int dim, int spacedim, typename value_type>
    value_type
    mean(const typename dealii::Triangulation<dim, spacedim>::cell_iterator
                                       &parent,
         const std::vector<value_type> &children_values);

    /**
     * Return maximum value of data on all children.
     *
     * @f[
     *   d_{K_p} = \max \left( d_{K_c} \right)
     *   \qquad
     *   \forall K_c \text{ children of } K_p
     * @f]
     */
    template <int dim, int spacedim, typename value_type>
    value_type
    max(const typename dealii::Triangulation<dim, spacedim>::cell_iterator
                                      &parent,
        const std::vector<value_type> &children_values);
  } // namespace Coarsening
} // namespace AdaptationStrategies



/* ---------------- template functions ---------------- */

#ifndef DOXYGEN

namespace AdaptationStrategies
{
  namespace Refinement
  {
    template <int dim, int spacedim, typename value_type>
    std::vector<value_type>
    preserve(const typename dealii::Triangulation<dim, spacedim>::cell_iterator
                             &parent,
             const value_type parent_value)
    {
      Assert(parent->n_children() > 0, ExcInternalError());
      return std::vector<value_type>(parent->n_children(), parent_value);
    }



    template <int dim, int spacedim, typename value_type>
    std::vector<value_type>
    split(const typename dealii::Triangulation<dim, spacedim>::cell_iterator
                          &parent,
          const value_type parent_value)
    {
      static_assert(std::is_arithmetic_v<value_type> &&
                      !std::is_same_v<value_type, bool>,
                    "The provided value_type may not meet the requirements "
                    "of this function.");

      Assert(parent->n_children() > 0, ExcInternalError());
      return std::vector<value_type>(parent->n_children(),
                                     parent_value / parent->n_children());
    }



    template <int dim, int spacedim, typename value_type>
    std::vector<value_type>
    l2_norm(const typename dealii::Triangulation<dim, spacedim>::cell_iterator
                            &parent,
            const value_type parent_value)
    {
      static_assert(std::is_arithmetic_v<value_type> &&
                      !std::is_same_v<value_type, bool>,
                    "The provided value_type may not meet the requirements "
                    "of this function.");

      Assert(parent->n_children() > 0, ExcInternalError());
      return std::vector<value_type>(parent->n_children(),
                                     parent_value /
                                       std::sqrt(parent->n_children()));
    }
  } // namespace Refinement



  namespace Coarsening
  {
    template <int dim, int spacedim, typename value_type>
    value_type
    check_equality(
      const typename dealii::Triangulation<dim, spacedim>::cell_iterator &,
      const std::vector<value_type> &children_values)
    {
      Assert(!children_values.empty(), ExcInternalError());

      const auto first_child = children_values.cbegin();
      for (auto other_child = first_child + 1;
           other_child != children_values.cend();
           ++other_child)
        Assert(*first_child == *other_child,
               ExcMessage(
                 "Values on cells that will be coarsened are not equal!"));

      return *first_child;
    }



    template <int dim, int spacedim, typename value_type>
    value_type
    sum(const typename dealii::Triangulation<dim, spacedim>::cell_iterator &,
        const std::vector<value_type> &children_values)
    {
      static_assert(std::is_arithmetic_v<value_type> &&
                      !std::is_same_v<value_type, bool>,
                    "The provided value_type may not meet the requirements "
                    "of this function.");

      Assert(!children_values.empty(), ExcInternalError());
      return std::accumulate(children_values.cbegin(),
                             children_values.cend(),
                             static_cast<value_type>(0));
    }



    template <int dim, int spacedim, typename value_type>
    value_type
    l2_norm(
      const typename dealii::Triangulation<dim, spacedim>::cell_iterator &,
      const std::vector<value_type> &children_values)
    {
      static_assert(std::is_arithmetic_v<value_type> &&
                      !std::is_same_v<value_type, bool>,
                    "The provided value_type may not meet the requirements "
                    "of this function.");

      Assert(!children_values.empty(), ExcInternalError());
      return std::sqrt(std::inner_product(children_values.cbegin(),
                                          children_values.cend(),
                                          children_values.cbegin(),
                                          static_cast<value_type>(0)));
    }



    template <int dim, int spacedim, typename value_type>
    value_type
    mean(const typename dealii::Triangulation<dim, spacedim>::cell_iterator
                                       &parent,
         const std::vector<value_type> &children_values)
    {
      return sum<dim, spacedim, value_type>(parent, children_values) /
             children_values.size();
    }



    template <int dim, int spacedim, typename value_type>
    value_type
    max(const typename dealii::Triangulation<dim, spacedim>::cell_iterator &,
        const std::vector<value_type> &children_values)
    {
      Assert(!children_values.empty(), ExcInternalError());
      return *std::max_element(children_values.cbegin(),
                               children_values.cend());
    }
  } // namespace Coarsening
} // namespace AdaptationStrategies

#endif

DEAL_II_NAMESPACE_CLOSE

#endif /* dealii_adaptation_strategies_h */
