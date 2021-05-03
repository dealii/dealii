// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2019 by the deal.II authors
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


// Test map_dofs_to_support_points on an hyper cube
// with two components. The points obtained for component 0
// of the dofhandler should be identical (and in the same order)
// as the points obtained for component 1


#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);

  FESystem<dim>   fe(FE_Q<dim>(1), 2);
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);

  ComponentMask mask(fe.n_components(), true);
  mask.set(1, false);

  MappingQ<dim>                                 mapping(1);
  std::map<types::global_dof_index, Point<dim>> support_points_0;
  DoFTools::map_dofs_to_support_points(mapping, dof, support_points_0, mask);

  mask.set(0, false);
  mask.set(1, true);

  std::map<types::global_dof_index, Point<dim>> support_points_1;
  DoFTools::map_dofs_to_support_points(mapping, dof, support_points_1, mask);

  auto it_0 = support_points_0.begin();
  auto it_1 = support_points_1.begin();

  for (; it_0 != support_points_0.end(); ++it_0, ++it_1)
    {
      deallog << it_0->second << " " << it_1->second << std::endl;
    }
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll init;

  deallog.push("2d");
  test<2>();
  deallog.pop();
  deallog.push("3d");
  test<3>();
  deallog.pop();
}
