// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2020 by the deal.II authors
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


#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include "../tests.h"

#include "dof_tools_common.h"
#include "dof_tools_common_fake_hp.h"

// check
//   DoFTools::
//   make_sparsity_pattern (const DoFHandler<dim>     &,
//                          Table<2,Coupling> &,
//                      DynamicSparsityPattern &);



template <typename DoFHandlerType>
void
check_this(const DoFHandlerType &dof_handler)
{
  // set up X-shape mask
  const unsigned int n_components = dof_handler.get_fe().n_components();
  Table<2, DoFTools::Coupling> mask(n_components, n_components);
  for (unsigned int i = 0; i < n_components; ++i)
    for (unsigned int j = 0; j < n_components; ++j)
      mask(i, j) = DoFTools::none;
  for (unsigned int i = 0; i < n_components; ++i)
    mask[i][i] = mask[i][n_components - i - 1] = DoFTools::always;

  // create sparsity pattern
  DynamicSparsityPattern sp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, mask, sp);
  sp.compress();

  // write out 20 lines of this
  // pattern (if we write out the
  // whole pattern, the output file
  // would be in the range of 40 MB)
  for (unsigned int l = 0; l < 20; ++l)
    {
      const unsigned int line = l * (sp.n_rows() / 20);
      for (unsigned int c = 0; c < sp.row_length(line); ++c)
        deallog << sp.column_number(line, c) << " ";
      deallog << std::endl;
    }

  // write out some other indicators
  deallog << sp.bandwidth() << std::endl
          << sp.max_entries_per_row() << std::endl
          << sp.n_nonzero_elements() << std::endl;

  unsigned int hash = 0;
  for (unsigned int l = 0; l < sp.n_rows(); ++l)
    hash += l * sp.row_length(l);
  deallog << hash << std::endl;
}
