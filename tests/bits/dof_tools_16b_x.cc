// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2014 by the deal.II authors
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
#include <deal.II/lac/compressed_set_sparsity_pattern.h>

// check
//   DoFTools::
//   make_boundary_sparsity_pattern (const DoFHandler<dim>     &,
//                                   const typename FunctionMap<dim>::type &
//                                   const std::vector<unsigned int> &
//                               CompressedSetSparsityPattern &);

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
  std::set<types::boundary_id> set;
  set.insert (0);
  DoFTools::map_dof_to_boundary_indices (dof_handler, set, map);

  // create sparsity pattern
  typename FunctionMap<dim>::type boundary_ids;
  boundary_ids[0] = 0;
  CompressedSetSparsityPattern sp (dof_handler.n_boundary_dofs(boundary_ids));
  DoFTools::make_boundary_sparsity_pattern (dof_handler, boundary_ids,
                                            map, sp);
  sp.compress ();

  // write out 20 lines of this
  // pattern (if we write out the
  // whole pattern, the output file
  // would be in the range of 40 MB)
  for (unsigned int l=0; l<20; ++l)
    {
      const unsigned int line = l*(sp.n_rows()/20);
      for (CompressedSetSparsityPattern::row_iterator
           c = sp.row_begin(line); c!=sp.row_end(line); ++c)
        deallog << *c << " ";
      deallog << std::endl;
    }

  // write out some other indicators
  deallog << sp.bandwidth () << std::endl
          << sp.max_entries_per_row () << std::endl
          << sp.n_nonzero_elements () << std::endl;

  unsigned int hash = 0;
  for (unsigned int l=0; l<sp.n_rows(); ++l)
    hash += l*sp.row_length(l);
  deallog << hash << std::endl;
}
