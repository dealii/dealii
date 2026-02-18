// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2023 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#ifndef dealii_cxx17_algorithm_h
#define dealii_cxx17_algorithm_h

#include <deal.II/base/config.h>

#include <algorithm>

#ifndef DEAL_II_BUILDING_CXX20_MODULE
DEAL_II_WARNING("This file is deprecated."
                "Use the corresponding C++17 header algorithm instead.")
#endif

DEAL_II_NAMESPACE_OPEN
namespace std_cxx17
{
  using std::clamp;
}
DEAL_II_NAMESPACE_CLOSE

#endif
