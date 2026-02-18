// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 1998 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#include <deal.II/algorithms/any_data.h>
#include <deal.II/algorithms/named_selection.h>

DEAL_II_NAMESPACE_OPEN

void
NamedSelection::initialize(const AnyData &data)
{
  indices.resize(names.size());
  for (unsigned int i = 0; i < names.size(); ++i)
    indices[i] = data.find(names[i]);
}

DEAL_II_NAMESPACE_CLOSE
