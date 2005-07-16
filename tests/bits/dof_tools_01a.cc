//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------

#include "../tests.h"
#include "dof_tools_common.cc"
#include <lac/sparsity_pattern.h>

// check
//   DoFTools::make_sparsity_pattern (
//     const DoFHandler<dim> &,
//     SparsityPattern       &);
// and
//   DoFTools::compute_row_length_vector
//    const DoFHandler<dim>&,
//    std::vector<unsigned int>&,
//    const Coupling)


std::string output_file_name = "dof_tools_01a.output";


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler)
{
                                   // create sparsity pattern
  std::vector<unsigned int> row_lengths(dof_handler.n_dofs());
  DoFTools::compute_row_length_vector(dof_handler, row_lengths);
  unsigned int total_length=0;
  for (unsigned int i=0;i<row_lengths.size();++i)
    total_length += row_lengths[i];
  
  SparsityPattern sp (dof_handler.n_dofs(), row_lengths);
  DoFTools::make_sparsity_pattern (dof_handler, sp);
  sp.compress ();
  
                                   // write out 20 lines of this
                                   // pattern (if we write out the
                                   // whole pattern, the output file
                                   // would be in the range of 40 MB)
  for (unsigned int l=0; l<20; ++l)
    {
      const unsigned int line = l*(sp.n_rows()/20);
      for (unsigned int c=0; c<sp.row_length(line); ++c)
        deallog << sp.column_number(line,c) << " ";
      deallog << std::endl;
    }

                                   // write out some other indicators
  deallog << "Bandwidth " << sp.bandwidth ()
          << "  Max per row " << sp.max_entries_per_row ()
          << " Nonzero " << sp.n_nonzero_elements ()
	  << " of "  << total_length << std::endl;

  unsigned int hash = 0;
  for (unsigned int l=0; l<sp.n_rows(); ++l)
    hash += l*(sp.row_length(l) +
               sp.get_rowstart_indices()[l] +
               sp.get_column_numbers()[sp.get_rowstart_indices()[l]
                                       +
                                       (sp.row_length(l)>1 ? 1 : 0)]);
  deallog << hash << std::endl;
}
