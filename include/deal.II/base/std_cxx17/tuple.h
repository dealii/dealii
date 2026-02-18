// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2019 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#ifndef dealii_cxx17_tuple_h
#define dealii_cxx17_tuple_h

#include <deal.II/base/config.h>

#include <tuple>

#ifndef DEAL_II_BUILDING_CXX20_MODULE
DEAL_II_WARNING("This file is deprecated."
                "Use the corresponding C++17 header tuple instead.")
#endif

DEAL_II_NAMESPACE_OPEN
namespace std_cxx17
{
  using std::apply;
}
DEAL_II_NAMESPACE_CLOSE

#endif
