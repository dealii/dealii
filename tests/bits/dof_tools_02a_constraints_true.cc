//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2003, 2004, 2007, 2009 by the deal.II authors
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
//   DoFTools::
//   make_sparsity_pattern (const DoFHandler<dim> &,
//                          std::vector<std::vector<bool> > &,
//	                    SparsityPattern       &,
//                          ConstraintMatrix,
//                          true);


std::string output_file_name = "dof_tools_02a_constraints_true/output";


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler)
{
                                   // set up X-shape mask
  const unsigned int n_components = dof_handler.get_fe().n_components();
  Table<2,DoFTools::Coupling> mask (n_components,n_components);
  for (unsigned int i=0; i<n_components; ++i)
    for (unsigned int j=0; j<n_components; ++j)
      mask[i][j] = DoFTools::none;
  for (unsigned int i=0; i<n_components; ++i)
    mask[i][i] = mask[i][n_components-i-1] = DoFTools::always;

  ConstraintMatrix cm;
  DoFTools::make_hanging_node_constraints (dof_handler, cm);
  cm.close ();
  
                                   // create sparsity pattern
  SparsityPattern sp (dof_handler.n_dofs(),
                      dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, mask, sp, cm, true);
  sp.compress ();
  
                                   // write out 10 lines of this
                                   // pattern (if we write out the
                                   // whole pattern, the output file
                                   // would be in the range of 40 MB)
  for (unsigned int l=0; l<10; ++l)
    {
      const unsigned int line = l*(sp.n_rows()/10);
      for (unsigned int c=0; c<sp.row_length(line); ++c)
        deallog << sp.column_number(line,c) << " ";
      deallog << std::endl;
    }

                                   // write out some other indicators
  deallog << sp.bandwidth () << std::endl
          << sp.max_entries_per_row () << std::endl
          << sp.n_nonzero_elements () << std::endl;

  unsigned int hash = 0;
  for (unsigned int l=0; l<sp.n_rows(); ++l)
    hash += l*(sp.row_length(l) +
               sp.get_rowstart_indices()[l] +
               sp.get_column_numbers()[sp.get_rowstart_indices()[l]
                                       +
                                       (sp.row_length(l)>1 ? 1 : 0)]);
  deallog << hash << std::endl;
}
