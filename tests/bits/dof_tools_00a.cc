//----------------------------  dof_tools_1c.cc  ---------------------------
//    $Id$
//    Version: 
//
//    Copyright (C) 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  dof_tools_1c.cc  ---------------------------

#include "../tests.h"
#include "dof_tools_common.cc"
#include <lac/block_sparsity_pattern.h>

// check
//   DoFTools::
//   count_dofs_per_component (...);

std::string output_file_name = "dof_tools_00a.output";


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler)
{
  const unsigned int n_components = dof_handler.get_fe().n_components();

  deallog << "n_dofs:" << dof_handler.n_dofs() << std::endl;
  
  std::vector<unsigned int> dofs_per_component(n_components);
  DoFTools::count_dofs_per_component (dof_handler,
				      dofs_per_component);

  for (unsigned int i=0;i<n_components;++i)
    deallog << ' ' << dofs_per_component[i];
  deallog << std::endl;

  if (n_components>1)
    {
      std::vector<unsigned int> target_component(n_components,0U);

      DoFTools::count_dofs_per_component (dof_handler,
					  dofs_per_component,
					  target_component);
      for (unsigned int i=0;i<n_components;++i)
	deallog << ' ' << dofs_per_component[i];
      deallog << std::endl;

      for (unsigned int i=n_components/2;i<n_components;++i)
	target_component[i] = 1;
      
      DoFTools::count_dofs_per_component (dof_handler,
					  dofs_per_component,
					  target_component);
      for (unsigned int i=0;i<n_components;++i)
	deallog << ' ' << dofs_per_component[i];
      deallog << std::endl;
    }
  
}
