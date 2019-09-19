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

namespace dealii
{
  namespace parallel
  {
    namespace fullydistributed
    {
      template <int dim>
      bool
      operator==(const CellData<dim> &t1, const CellData<dim> &t2)
      {
        if (t1.id != t2.id)
          return false;
        if (t1.subdomain_id != t2.subdomain_id)
          return false;
        if (t1.level_subdomain_id != t2.level_subdomain_id)
          return false;
        if (t1.manifold_id != t2.manifold_id)
          return false;
        if (dim >= 2 && t1.manifold_line_ids != t2.manifold_line_ids)
          return false;
        if (dim >= 3 && t1.manifold_quad_ids != t2.manifold_quad_ids)
          return false;
        if (t1.boundary_ids != t2.boundary_ids)
          return false;

        return true;
      }

      template <int dim, int spacedim>
      bool
      operator==(const ConstructionData<dim, spacedim> &t1,
                 const ConstructionData<dim, spacedim> &t2)
      {
        if (t1.coarse_cells != t2.coarse_cells)
          return false;
        if (t1.coarse_cell_vertices != t2.coarse_cell_vertices)
          return false;
        if (t1.coarse_cell_index_to_coarse_cell_id !=
            t2.coarse_cell_index_to_coarse_cell_id)
          return false;
        if (t1.cell_infos != t2.cell_infos)
          return false;
        if (t1.settings != t2.settings)
          return false;

        return true;
      }
    } // namespace fullydistributed
  }   // namespace parallel
} // namespace dealii

#endif
