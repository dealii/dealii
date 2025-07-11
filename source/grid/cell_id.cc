// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/grid/cell_id.h>
#include <deal.II/grid/tria.h>

#include <limits>
#include <sstream>

DEAL_II_NAMESPACE_OPEN


CellId::CellId()
  : coarse_cell_id(numbers::invalid_coarse_cell_id)
  , n_child_indices(numbers::invalid_unsigned_int)
{
  // initialize the child indices to invalid values
  // (the only allowed values are between zero and
  // GeometryInfo<dim>::max_children_per_cell)
  std::fill(child_indices.begin(),
            child_indices.end(),
            std::numeric_limits<char>::max());
}



CellId::CellId(const types::coarse_cell_id      coarse_cell_id,
               const std::vector<std::uint8_t> &id)
  : coarse_cell_id(coarse_cell_id)
  , n_child_indices(id.size())
{
  Assert(n_child_indices < child_indices.size(), ExcInternalError());
  std::copy(id.begin(), id.end(), child_indices.begin());
}



CellId::CellId(const types::coarse_cell_id coarse_cell_id,
               const unsigned int          n_child_indices,
               const std::uint8_t         *id)
  : coarse_cell_id(coarse_cell_id)
  , n_child_indices(n_child_indices)
{
  Assert(n_child_indices < child_indices.size(), ExcInternalError());
  std::memcpy(child_indices.data(), id, n_child_indices);
}



CellId::CellId(const CellId::binary_type &binary_representation)
{
  // The first entry stores the coarse cell id
  coarse_cell_id = binary_representation[0];

  // The rightmost two bits of the second entry store the dimension,
  // the rest stores the number of child indices.
  const unsigned int two_bit_mask = (1 << 2) - 1;
  const unsigned int dim          = binary_representation[1] & two_bit_mask;
  n_child_indices                 = (binary_representation[1] >> 2);

  Assert(n_child_indices < child_indices.size(), ExcInternalError());

  // Each child requires 'dim' bits to store its index
  const unsigned int children_per_value =
    sizeof(binary_type::value_type) * 8 / dim;
  const unsigned int child_mask = (1 << dim) - 1;

  // Loop until all child indices have been read
  unsigned int child_level  = 0;
  unsigned int binary_entry = 2;
  while (child_level < n_child_indices)
    {
      for (unsigned int j = 0; j < children_per_value; ++j)
        {
          // Read the current child index by shifting to the current
          // index's position and doing a bitwise-and with the child_mask.
          child_indices[child_level] =
            (binary_representation[binary_entry] >> (dim * j)) & child_mask;
          ++child_level;
          if (child_level == n_child_indices)
            break;
        }
      ++binary_entry;
    }
}



CellId::CellId(const std::string &string_representation)
{
  std::istringstream ss(string_representation);
  ss >> *this;
}



template <int dim>
CellId::binary_type
CellId::to_binary() const
{
  CellId::binary_type binary_representation;
  binary_representation.fill(0);

  Assert(n_child_indices < child_indices.size(), ExcInternalError());

  // The first entry stores the coarse cell id
  binary_representation[0] = coarse_cell_id;

  // The rightmost two bits of the second entry store the dimension,
  // the rest stores the number of child indices.
  binary_representation[1] = (n_child_indices << 2);
  binary_representation[1] |= dim;

  // Each child requires 'dim' bits to store its index
  const unsigned int children_per_value =
    sizeof(binary_type::value_type) * 8 / dim;
  unsigned int child_level  = 0;
  unsigned int binary_entry = 2;

  // Loop until all child indices have been written
  while (child_level < n_child_indices)
    {
      Assert(binary_entry < binary_representation.size(), ExcInternalError());

      for (unsigned int j = 0; j < children_per_value; ++j)
        {
          const unsigned int child_index =
            static_cast<unsigned int>(child_indices[child_level]);
          // Shift the child index to its position in the unsigned int and store
          // it
          binary_representation[binary_entry] |= (child_index << (j * dim));
          ++child_level;
          if (child_level == n_child_indices)
            break;
        }
      ++binary_entry;
    }

  return binary_representation;
}



std::string
CellId::to_string() const
{
  std::ostringstream ss;
  ss << *this;
  return ss.str();
}


// explicit instantiations
#include "grid/cell_id.inst"

DEAL_II_NAMESPACE_CLOSE
