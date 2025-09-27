// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/memory_consumption.h>

#include <deal.II/fe/block_mask.h>

#include <iostream>


DEAL_II_NAMESPACE_OPEN

std::ostream &
operator<<(std::ostream &out, const BlockMask &mask)
{
  if (mask.block_mask.empty())
    out << "[all blocks selected]";
  else
    {
      out << '[';
      for (unsigned int i = 0; i < mask.block_mask.size(); ++i)
        {
          out << (mask.block_mask[i] ? "true" : "false");
          if (i != mask.block_mask.size() - 1)
            out << ',';
        }
      out << ']';
    }

  return out;
}



std::size_t
BlockMask::memory_consumption() const
{
  return sizeof(*this) + MemoryConsumption::memory_consumption(block_mask);
}


DEAL_II_NAMESPACE_CLOSE
