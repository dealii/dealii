// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2015 by the deal.II authors
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

#include <deal.II/numerics/matrix_creator.templates.h>


DEAL_II_NAMESPACE_OPEN

// explicit instantiations
#define SPLIT_INSTANTIATIONS_COUNT 3
#ifndef SPLIT_INSTANTIATIONS_INDEX
#define SPLIT_INSTANTIATIONS_INDEX 0
#endif
#include "matrix_creator.inst"

DEAL_II_NAMESPACE_CLOSE
