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
//   make_boundary_sparsity_pattern (const DoFHandler<dim> &,
//                                   const std::vector<unsigned int> &
//                               BlockCompressedSparsityPattern  &);

std::string output_file_name = "output";


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler)
{
  // test doesn't make much sense if
  // no boundary dofs exist
  if (dof_handler.get_fe().dofs_per_face == 0)
    return;

  std::vector<types::global_dof_index> map (dof_handler.n_dofs());
  DoFTools::map_dof_to_boundary_indices (dof_handler, map);

  const unsigned int n_blocks = std::min (
                                  static_cast<types::global_dof_index>(dof_handler.get_fe().n_components()),
                                  dof_handler.n_boundary_dofs());
  BlockCompressedSparsityPattern sp (n_blocks,
                                     n_blocks);
  // split dofs almost arbitrarily to
  // blocks
  std::vector<types::global_dof_index> dofs_per_block(n_blocks);
  for (unsigned int i=0; i<n_blocks-1; ++i)
    dofs_per_block[i] = dof_handler.n_boundary_dofs()/n_blocks;
  dofs_per_block.back() = (dof_handler.n_boundary_dofs() -
                           (dof_handler.n_boundary_dofs()/n_blocks)*(n_blocks-1));

  for (unsigned int i=0; i<n_blocks; ++i)
    for (unsigned int j=0; j<n_blocks; ++j)
      sp.block(i,j).reinit(dofs_per_block[i],
                           dofs_per_block[j]);
  sp.collect_sizes ();

  DoFTools::make_boundary_sparsity_pattern (dof_handler, map, sp);
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
      for (unsigned int col=0; col<n_blocks; ++col)
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
  for (unsigned int r=0; r<n_blocks; ++r)
    for (unsigned int c=0; c<n_blocks; ++c)
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

