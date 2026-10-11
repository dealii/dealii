// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// Test map_boundary_to_bulk_dof_iterators

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"



template <int dim>
void
test()
{
  Triangulation<dim>          triangulation;
  Triangulation<dim - 1, dim> surface_triangulation;

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof_handler(triangulation);

  FE_Q<dim - 1, dim>       surface_fe(1);
  DoFHandler<dim - 1, dim> surface_dof_handler(surface_triangulation);

  GridGenerator::half_hyper_ball(triangulation);
  triangulation.refine_global(4 - dim);

  surface_triangulation.set_manifold(0, SphericalManifold<dim - 1, dim>());
  const auto surface_to_bulk_map =
    GridGenerator::extract_boundary_mesh(triangulation,
                                         surface_triangulation,
                                         {0});

  deallog << "Bulk mesh active cells:" << triangulation.n_active_cells()
          << std::endl
          << "Surface mesh active cells:"
          << surface_triangulation.n_active_cells() << std::endl;

  dof_handler.distribute_dofs(fe);
  surface_dof_handler.distribute_dofs(surface_fe);

  // Log degrees of freedom:
  deallog << "Bulk mesh degrees of freedom:" << dof_handler.n_dofs()
          << std::endl
          << "Surface mesh degrees of freedom:" << surface_dof_handler.n_dofs()
          << std::endl;

  // Extract the mapping between surface and bulk degrees of freedom:
  const auto surface_to_bulk_dof_iterator_map =
    DoFTools::map_boundary_to_bulk_dof_iterators(surface_to_bulk_map,
                                                 dof_handler,
                                                 surface_dof_handler);

  // Loop over the map, and print some information:
  for (const auto &p : surface_to_bulk_dof_iterator_map)
    {
      const auto &surface_cell = p.first;
      const auto &bulk_cell    = p.second.first;
      const auto &bulk_face    = p.second.second;
      deallog << "Surface cell " << surface_cell << " coincides with face "
              << bulk_face << " of bulk cell " << bulk_cell << std::endl;
    }
}



int
main()
{
  initlog();

  test<2>();
  test<3>();
}
