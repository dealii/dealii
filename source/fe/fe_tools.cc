// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2000 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/fe/fe_tools.templates.h>

DEAL_II_NAMESPACE_OPEN


namespace FETools
{
  namespace internal
  {
    namespace FEToolsAddFENameHelper
    {
      std::shared_mutex fe_name_map_lock;
    }
  } // namespace internal
} // namespace FETools


/*-------------- Explicit Instantiations -------------------------------*/
#include "fe/fe_tools.inst"

// these do not fit into the templates of the dimension in the inst file
namespace FETools
{
  template std::vector<unsigned int>
  hierarchic_to_lexicographic_numbering<0>(unsigned int);

  template std::vector<unsigned int>
  lexicographic_to_hierarchic_numbering<0>(unsigned int);
} // namespace FETools

DEAL_II_NAMESPACE_CLOSE
