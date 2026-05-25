// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2015 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#include <deal.II/grid/cell_id.h>
#include <deal.II/grid/tria.h>

#include <limits>
#include <sstream>

DEAL_II_NAMESPACE_OPEN


CellId::CellId()
  : coarse_cell_id(numbers::invalid_coarse_cell_id)
  , n_child_indices(std::numeric_limits<std::uint8_t>::max())
{
  // initialize the child indices to invalid values
  std::fill(child_indices.begin(),
            child_indices.end(),
            std::numeric_limits<std::uint8_t>::max());
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
  : coarse_cell_id(binary_representation[0])
  , n_child_indices(static_cast<std::uint8_t>(binary_representation[2] & 0xFF))
{
  Assert(n_child_indices <= child_indices.size(), ExcInternalError());
  child_indices.fill(std::numeric_limits<std::uint8_t>::max());

  // This is the reverse of to_binary()
  unsigned int binary_entry_index        = 1;
  int          index_within_binary_entry = 7;
  for (std::uint8_t child_no = 0; child_no < n_child_indices; ++child_no)
    {
      binary_type::value_type mask = 0;
      if (child_no % 2 == 0)
        mask |= (0x0F << 4);
      else
        mask |= 0x0F;

      AssertIndexRange(binary_entry_index, binary_representation.size());
      auto value = (binary_representation[binary_entry_index] >>
                    (8 * index_within_binary_entry)) &
                   mask;
      if (child_no % 2 == 0)
        value >>= 4;

      // We don't have dim available here, so make sure everything fits in the
      // 3d limit instead
      AssertIndexRange(value, ReferenceCells::max_n_children<3>());
      child_indices[child_no] = static_cast<std::uint8_t>(value);

      if (child_no % 2 == 1)
        --index_within_binary_entry;
      if (index_within_binary_entry == -1)
        {
          ++binary_entry_index;
          index_within_binary_entry = 7;
        }
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
  binary_type binary_representation;
  binary_representation.fill(0);

  binary_representation[0]               = coarse_cell_id;
  unsigned int binary_entry_index        = 1;
  int          index_within_binary_entry = 7;
  // This encoding only makes sense if we can use 7 as an index of sorts (see
  // below) into that binary type, so future-proof a little by asserting that
  // binary_type::value_type is big enough for this to work:
  Assert(index_within_binary_entry < int(sizeof(binary_type::value_type)),
         ExcInternalError());
  // Encode each child as a distinct power of 2 inside the binary
  // representation. Rather than memcpy(), this representation is constructed
  // from high values to low values so that we have a unique numerical value in
  // the integers of binary_representation for any cell: that way we get results
  // independent of the endianness of the machine (i.e., text_oarchive will have
  // the same results across platforms). Even though we don't need it, this
  // particular choice preserves lexical ordering (i.e., if cell_id_1 <
  // cell_id_2, then cell_id_1.to_binary() < cell_id_2.to_binary()).
  for (std::uint8_t child_no = 0; child_no < n_child_indices; ++child_no)
    {
      Assert(child_indices[child_no] < ReferenceCells::max_n_children<dim>(),
             ExcInternalError());
      Assert(ReferenceCells::max_n_children<dim>() <= 0x0F, ExcInternalError());

      binary_type::value_type entry = 0;

      if (child_no % 2 == 0)
        entry |= (child_indices[child_no] & 0x0F) << 4;
      else
        entry |= child_indices[child_no] & 0x0F;

      AssertIndexRange(binary_entry_index, binary_representation.size());
      binary_representation[binary_entry_index] |=
        (entry << 8 * index_within_binary_entry);

      if (child_no % 2 == 1)
        --index_within_binary_entry;
      if (index_within_binary_entry == -1)
        {
          ++binary_entry_index;
          index_within_binary_entry = 7;
        }
    }

  if (binary_entry_index == 2)
    Assert(index_within_binary_entry > 0, ExcInternalError());
  // finally store the number of children at the end:
  binary_representation[2] |= n_child_indices;

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
