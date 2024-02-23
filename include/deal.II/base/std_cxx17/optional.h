// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------
#ifndef dealii_cxx17_optional_h
#define dealii_cxx17_optional_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_EARLY_DEPRECATIONS
DEAL_II_WARNING(
  "This file is deprecated. Simply use the corresponding C++17 header <optional>.")
#endif

#include <optional>

DEAL_II_NAMESPACE_OPEN
namespace std_cxx17
{
  using std::optional;
} // namespace std_cxx17
DEAL_II_NAMESPACE_CLOSE

#endif // dealii_cxx17_optional_h
