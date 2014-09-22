// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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

#include <deal.II/base/named_data.h>
#include <deal.II/algorithms/any_data.h>

DEAL_II_NAMESPACE_OPEN

void
NamedSelection::initialize(const AnyData &data)
{
  indices.resize(names.size());
  for (unsigned int i=0; i<names.size(); ++i)
    indices[i] = data.find(names[i]);
}

DEAL_II_NAMESPACE_CLOSE

