// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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

#ifndef dealii_fully_distributed_tests_h
#define dealii_fully_distributed_tests_h

#include "../tests.h"

using namespace dealii;

template <int dim, int spacedim>
void
print_statistics(
  const parallel::fullydistributed::Triangulation<dim, spacedim> &tria,
  bool                                                            do_mg = false)
{
  deallog << "n_levels:                  " << tria.n_levels() << std::endl;
  deallog << "n_cells:                   " << tria.n_cells() << std::endl;
  deallog << "n_active_cells:            " << tria.n_active_cells()
          << std::endl;

  if (do_mg)
    {
      for (auto level = 0u; level < tria.n_levels(); level++)
        deallog << "n_cells on level=" << level << ":        "
                << tria.n_cells(level) << std::endl;

      for (auto level = 0u; level < tria.n_levels(); level++)
        deallog << "n_active_cells on level=" << level << ": "
                << tria.n_active_cells(level) << std::endl;
    }

  deallog << std::endl;
}

template <int dim, int spacedim>
void
print_statistics(const DoFHandler<dim, spacedim> &dof_handler,
                 bool                             do_mg = false)
{
  deallog << "n_dofs:                             " << dof_handler.n_dofs()
          << std::endl;
  deallog << "n_locally_owned_dofs:               "
          << dof_handler.n_locally_owned_dofs() << std::endl;

  const auto n_levels = dof_handler.get_triangulation().n_levels();

  if (do_mg)
    {
      for (auto level = 0u; level < n_levels; level++)
        deallog << "n_dofs on level=" << level << ":                  "
                << dof_handler.n_dofs(level) << std::endl;

      for (auto level = 0u; level < n_levels; level++)
        deallog << "n_locally_owned_mg_dofs on level=" << level << ": "
                << dof_handler.locally_owned_mg_dofs(level).n_elements()
                << std::endl;
    }

  deallog << std::endl;
}

#endif
