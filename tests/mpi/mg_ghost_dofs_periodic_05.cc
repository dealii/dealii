// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// this tests distribute_mg_dofs on a mesh with periodic boundary
// conditions. we used to forget to enforce the 2:1 level balance over
// periodic boundaries in artificial cells, which lead to the case where the
// coarser MG level did not provide the correct ghost layer (leading to a
// deadlock in the distribute_dofs function for 6 procs in 2D and 28 procs in
// 3D)

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> tria(
    MPI_COMM_WORLD,
    dealii::Triangulation<dim>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);
  GridGenerator::subdivided_hyper_cube(tria, 1);
  // set periodic boundary conditions in x and z directions
  for (typename Triangulation<dim>::cell_iterator cell = tria.begin();
       cell != tria.end();
       ++cell)
    for (const unsigned int f : GeometryInfo<dim>::face_indices())
      if (f / 2 != 1 && cell->at_boundary(f))
        cell->face(f)->set_all_boundary_ids(f + 10);

  std::vector<
    GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator>>
    periodic_faces;
  GridTools::collect_periodic_faces(tria, 0 + 10, 1 + 10, 0, periodic_faces);
  if (dim == 3)
    GridTools::collect_periodic_faces(tria, 4 + 10, 5 + 10, 2, periodic_faces);
  tria.add_periodicity(periodic_faces);
  tria.refine_global(5);

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);
  deallog << "Number of cells: " << tria.n_global_active_cells() << std::endl
          << "Number of DoFs: " << dof_handler.n_dofs() << std::endl;

  dof_handler.distribute_mg_dofs();
  deallog << "Number of DoFs per level: ";
  for (unsigned int level = 0; level < tria.n_global_levels(); ++level)
    deallog << dof_handler.n_dofs(level) << ' ';
  deallog << std::endl;
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  mpi_initlog();

  test<2>();
  test<3>();

  return 0;
}
