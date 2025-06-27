// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_cxx17_variant_h
#define dealii_cxx17_variant_h

#include <deal.II/base/config.h>

#include <variant>

#ifdef DEAL_II_EARLY_DEPRECATIONS
DEAL_II_WARNING("This file is deprecated."
                "Use the corresponding C++17 header variant instead.")
#endif

DEAL_II_NAMESPACE_OPEN
namespace std_cxx17
{
  using std::get;
  using std::holds_alternative;
  using std::variant;
} // namespace std_cxx17
DEAL_II_NAMESPACE_CLOSE

#endif
