// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/lac/sparsity_pattern_base.h>

#include <boost/container/small_vector.hpp>

#include <algorithm>
#include <utility>

DEAL_II_NAMESPACE_OPEN

void
SparsityPatternBase::add_entries(
  const ArrayView<const std::pair<size_type, size_type>> &inputs)
{
  // We always want to sort so that we can use the optimized add_row_entries()
  // function
  boost::container::small_vector<std::pair<size_type, size_type>, 128> entries(
    inputs.begin(), inputs.end());
  std::sort(entries.begin(), entries.end());
  boost::container::small_vector<size_type, 128> columns;

  auto entry = entries.begin();
  while (entry != entries.end())
    {
      const auto row     = entry->first;
      auto       row_end = entry;
      while (row_end != entries.end() && row_end->first == row)
        ++row_end;

      columns.resize(0);
      columns.reserve(row_end - entry);
      // Simultaneously transform and uniquify
      columns.push_back(entry->second);
      ++entry;
      while (entry != row_end)
        {
          if (columns.back() != entry->second)
            columns.push_back(entry->second);
          ++entry;
        }

      add_row_entries(row,
                      make_array_view(columns.begin(), columns.end()),
                      true);
    }
}

DEAL_II_NAMESPACE_CLOSE
