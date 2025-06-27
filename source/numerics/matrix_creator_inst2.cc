// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/numerics/matrix_creator.templates.h>


DEAL_II_NAMESPACE_OPEN

// explicit instantiations
#define SPLIT_INSTANTIATIONS_COUNT 3
#define SPLIT_INSTANTIATIONS_INDEX 1

#include "numerics/matrix_creator.inst"

DEAL_II_NAMESPACE_CLOSE
