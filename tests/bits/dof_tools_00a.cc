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

// check
//   DoFTools::
//   count_dofs_per_component (...);

std::string output_file_name = "output";


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler)
{
  const unsigned int n_components = dof_handler.get_fe().n_components();

  deallog << "n_dofs:" << dof_handler.n_dofs() << std::endl;

  std::vector<types::global_dof_index> dofs_per_component(n_components);
  DoFTools::count_dofs_per_component (dof_handler,
                                      dofs_per_component);

  for (unsigned int i=0; i<n_components; ++i)
    deallog << ' ' << dofs_per_component[i];
  deallog << std::endl;

  DoFTools::count_dofs_per_component (dof_handler,
                                      dofs_per_component, true);

  for (unsigned int i=0; i<n_components; ++i)
    deallog << ' ' << dofs_per_component[i];
  deallog << std::endl;

  if (n_components>1)
    {
      std::vector<unsigned int> target_component(n_components,0U);
      dofs_per_component.resize (1);
      DoFTools::count_dofs_per_component (dof_handler,
                                          dofs_per_component,
                                          false,
                                          target_component);
      for (unsigned int i=0; i<std::min(n_components, (unsigned int)dofs_per_component.size()); ++i)
        deallog << ' ' << dofs_per_component[i];
      deallog << std::endl;

      DoFTools::count_dofs_per_component (dof_handler,
                                          dofs_per_component,
                                          true,
                                          target_component);
      for (unsigned int i=0; i<std::min(n_components, (unsigned int)dofs_per_component.size()); ++i)
        deallog << ' ' << dofs_per_component[i];
      deallog << std::endl;

      for (unsigned int i=n_components/2; i<n_components; ++i)
        target_component[i] = 1;
      dofs_per_component.resize (2);

      DoFTools::count_dofs_per_component (dof_handler,
                                          dofs_per_component,
                                          false,
                                          target_component);
      for (unsigned int i=0; i<std::min(n_components, (unsigned int)dofs_per_component.size()); ++i)
        deallog << ' ' << dofs_per_component[i];
      deallog << std::endl;
      DoFTools::count_dofs_per_component (dof_handler,
                                          dofs_per_component,
                                          true,
                                          target_component);
      for (unsigned int i=0; i<std::min(n_components, (unsigned int)dofs_per_component.size()); ++i)
        deallog << ' ' << dofs_per_component[i];
      deallog << std::endl;
    }

}
