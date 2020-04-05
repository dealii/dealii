// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


#include <deal.II/fe/fe_tools.templates.h>

DEAL_II_NAMESPACE_OPEN

/*-------------- Explicit Instantiations -------------------------------*/
#include "fe_tools.inst"

// these do not fit into the templates of the dimension in the inst file
namespace FETools
{
  template void
  hierarchic_to_lexicographic_numbering<0>(unsigned int,
                                           std::vector<unsigned int> &);

  template std::vector<unsigned int>
  hierarchic_to_lexicographic_numbering<0>(unsigned int);

  template std::vector<unsigned int>
  lexicographic_to_hierarchic_numbering<0>(unsigned int);
} // namespace FETools

DEAL_II_NAMESPACE_CLOSE
