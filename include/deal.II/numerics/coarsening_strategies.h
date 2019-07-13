// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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

#ifndef dealii_coarsening_strategies_h
#define dealii_coarsening_strategies_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <algorithm>
#include <numeric>
#include <typeinfo>
#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * When data is transferred during coarsening, it is not trivial to decide
 * how to handle data of child cells which will be coarsened. Or in other
 * words, which data should be stored in the corresponding parent cell.
 *
 * In this namespace, we offer a few strategies that cope with this
 * problem. Such strategies can be passed to the CellDataTransfer and
 * parallel::distributed::CellDataTransfer constructors.
 *
 * @author Marc Fehling, 2019
 */
namespace CoarseningStrategies
{
  /**
   * Check if data on all children match, and return that value.
   */
  template <typename value_type>
  value_type
  check_equality(const std::vector<value_type> &children_values);

  /**
   * Return sum of data on all children.
   */
  template <typename value_type>
  value_type
  sum(const std::vector<value_type> &children_values);

  /**
   * Return mean value of data on all children.
   */
  template <typename value_type>
  value_type
  mean(const std::vector<value_type> &children_values);

  /**
   * Return maximum value of data on all children.
   */
  template <typename value_type>
  value_type
  max(const std::vector<value_type> &children_values);
} // namespace CoarseningStrategies



/* ---------------- template functions ---------------- */

#ifndef DOXYGEN

namespace CoarseningStrategies
{
  template <typename value_type>
  value_type
  check_equality(const std::vector<value_type> &children_values)
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



  template <typename value_type>
  value_type
  sum(const std::vector<value_type> &children_values)
  {
    static_assert(std::is_arithmetic<value_type>::value &&
                    !std::is_same<value_type, bool>::value,
                  "The provided value_type may not meet the requirements "
                  "of this function.");

    Assert(!children_values.empty(), ExcInternalError());
    return std::accumulate(children_values.cbegin(),
                           children_values.cend(),
                           static_cast<value_type>(0));
  }



  template <typename value_type>
  value_type
  mean(const std::vector<value_type> &children_values)
  {
    return sum<value_type>(children_values) / children_values.size();
  }



  template <typename value_type>
  value_type
  max(const std::vector<value_type> &children_values)
  {
    Assert(!children_values.empty(), ExcInternalError());
    return *std::max_element(children_values.cbegin(), children_values.cend());
  }
} // namespace CoarseningStrategies

#endif

DEAL_II_NAMESPACE_CLOSE

#endif /* dealii_coarsening_strategies_h */
