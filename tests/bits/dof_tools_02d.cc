// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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


#include "../tests.h"
#include "dof_tools_common.h"
#include <deal.II/lac/block_sparsity_pattern.h>

// check
//   DoFTools::
//   make_sparsity_pattern (const DoFHandler<dim> &,
//                          std::vector<std::vector<bool> > &,
//                      BlockCompressedSparsityPattern  &);

std::string output_file_name = "output";


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler)
{
  // set up X-shape mask
  const unsigned int n_components = dof_handler.get_fe().n_components();
  std::vector<std::vector<bool> > mask (n_components,
                                        std::vector<bool>(n_components,false));
  for (unsigned int i=0; i<n_components; ++i)
    mask[i][i] = mask[i][n_components-i-1] = true;

  // we split up the matrix into
  // blocks according to the number
  // of dofs in each component. this
  // fails if the element is not
  // primitive, so skip this test for
  // such elements
  if (dof_handler.get_fe().is_primitive() != true)
    return;

  // create sparsity pattern
  BlockCompressedSparsityPattern sp (n_components,
                                     n_components);
  std::vector<types::global_dof_index> dofs_per_component(n_components);
  DoFTools::count_dofs_per_component (dof_handler,
                                      dofs_per_component);
  for (unsigned int i=0; i<n_components; ++i)
    for (unsigned int j=0; j<n_components; ++j)
      sp.block(i,j).reinit(dofs_per_component[i],
                           dofs_per_component[j]);
  sp.collect_sizes ();

  DoFTools::make_sparsity_pattern (dof_handler, mask, sp);
  sp.compress ();

  // write out 20 lines of this
  // pattern (if we write out the
  // whole pattern, the output file
  // would be in the range of 40 MB)
  for (unsigned int l=0; l<20; ++l)
    {
      const unsigned int line = l*(sp.n_rows()/20);
      std::pair<unsigned int,unsigned int>
      block_row = sp.get_row_indices().global_to_local(line);
      for (unsigned int col=0; col<n_components; ++col)
        {
          for (unsigned int c=0;
               c<sp.block(block_row.first,col).row_length(block_row.second);
               ++c)
            deallog << sp.block(block_row.first,col).column_number(block_row.second,c)
                    << " ";
          deallog << std::endl;
        }
    }

  // write out some other indicators
  for (unsigned int r=0; r<n_components; ++r)
    for (unsigned int c=0; c<n_components; ++c)
      {
        const CompressedSparsityPattern &x = sp.block(r,c);
        deallog << x.bandwidth () << std::endl
                << x.max_entries_per_row () << std::endl
                << x.n_nonzero_elements () << std::endl;

        unsigned int hash = 0;
        for (unsigned int l=0; l<x.n_rows(); ++l)
          hash += l*x.row_length(l);
        deallog << hash << std::endl;
      }
}

