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

// check
//   DoFTools::compute_row_length_vector(...


std::string output_file_name = "dof_tools_row_1/output";


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler)
{
  std::vector<unsigned int>
    row_length(dof_handler.n_dofs());
  DoFTools::compute_row_length_vector(dof_handler, row_length, DoFTools::none);
  
  for (unsigned int i=0;i<dof_handler.n_dofs();++i)
    {
      deallog << "  " << row_length[i];
      deallog << std::endl;
    }
}

