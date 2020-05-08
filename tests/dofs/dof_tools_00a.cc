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


#include "../tests.h"

#include "dof_tools_common.h"
#include "dof_tools_common_fake_hp.h"

// check
//   DoFTools::
//   count_dofs_per_component (...);



template <typename DoFHandlerType>
void
check_this(const DoFHandlerType &dof_handler)
{
  const unsigned int n_components = dof_handler.get_fe().n_components();

  deallog << "n_dofs:" << dof_handler.n_dofs() << std::endl;

  std::vector<types::global_dof_index> dofs_per_component =
    DoFTools::count_dofs_per_fe_component(dof_handler);

  for (unsigned int i = 0; i < n_components; ++i)
    deallog << ' ' << dofs_per_component[i];
  deallog << std::endl;

  dofs_per_component = DoFTools::count_dofs_per_fe_component(dof_handler, true);

  for (unsigned int i = 0; i < n_components; ++i)
    deallog << ' ' << dofs_per_component[i];
  deallog << std::endl;

  if (n_components > 1)
    {
      std::vector<unsigned int> target_component(n_components, 0U);
      dofs_per_component =
        DoFTools::count_dofs_per_fe_component(dof_handler,
                                              false,
                                              target_component);
      for (unsigned int i = 0;
           i < std::min(n_components, (unsigned int)dofs_per_component.size());
           ++i)
        deallog << ' ' << dofs_per_component[i];
      deallog << std::endl;

      dofs_per_component =
        DoFTools::count_dofs_per_fe_component(dof_handler,
                                              true,
                                              target_component);
      for (unsigned int i = 0;
           i < std::min(n_components, (unsigned int)dofs_per_component.size());
           ++i)
        deallog << ' ' << dofs_per_component[i];
      deallog << std::endl;

      for (unsigned int i = n_components / 2; i < n_components; ++i)
        target_component[i] = 1;
      dofs_per_component.resize(2);

      dofs_per_component =
        DoFTools::count_dofs_per_fe_component(dof_handler,
                                              false,
                                              target_component);
      for (unsigned int i = 0;
           i < std::min(n_components, (unsigned int)dofs_per_component.size());
           ++i)
        deallog << ' ' << dofs_per_component[i];
      deallog << std::endl;
      dofs_per_component =
        DoFTools::count_dofs_per_fe_component(dof_handler,
                                              true,
                                              target_component);
      for (unsigned int i = 0;
           i < std::min(n_components, (unsigned int)dofs_per_component.size());
           ++i)
        deallog << ' ' << dofs_per_component[i];
      deallog << std::endl;
    }
}
