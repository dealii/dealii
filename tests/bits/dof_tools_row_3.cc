//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2003, 2004, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------

#include "../tests.h"
#include "dof_tools_common.cc"
#include <lac/block_sparsity_pattern.h>

#include <iomanip>

// check
//   DoFTools::compute_row_length_vector(...


std::string output_file_name = "dof_tools_row_3/output";


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler)
{
  if (dim==1)
    return;
                                   // create sparsity pattern
  const unsigned int n_components = dof_handler.get_fe().n_components();
  const unsigned int n_blocks = dof_handler.get_fe().n_blocks();
  
  std::vector<unsigned int> dofs_per_block(n_blocks);
  DoFTools::count_dofs_per_block (dof_handler,
				  dofs_per_block);
  
  std::vector<std::vector<unsigned int> >
    row_length(n_blocks, std::vector<unsigned int>(dof_handler.n_dofs()));
  Table<2,DoFTools::Coupling> couplings(n_components, n_components);
  couplings.fill(DoFTools::always);
  DoFTools::compute_row_length_vector(dof_handler, row_length,
				      couplings, couplings);
  
  for (unsigned int i=0;i<dof_handler.n_dofs();++i)
    {
      deallog << std::setw(5) << i;
      for (unsigned int j=0;j<n_blocks;++j)
	deallog << std::setw(5) << row_length[j][i];
      deallog << std::endl;
    }
}

