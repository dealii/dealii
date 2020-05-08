// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2020 by the deal.II authors
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



// create a parallel hp::DoFHandler on a single CPU
//
// like the test without the hp_ prefix, but for hp::DoFHandler

#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_collection.h>

#include <fstream>

#include "../tests.h"

#include "coarse_grid_common.h"


template <int dim>
void
test()
{
  deallog << "hyper_cube" << std::endl;

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr);
  tr.refine_global(1);

  hp::FECollection<dim> fe;
  fe.push_back(FE_Q<dim>(2));
  hp::DoFHandler<dim> dofh(tr);
  dofh.distribute_dofs(fe);

  typename hp::DoFHandler<dim>::active_cell_iterator cell = dofh.begin_active();

  const unsigned int dofs_per_cell = dofh.get_fe(0).dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);



  for (; cell != dofh.end(); ++cell)
    {
      cell->get_dof_indices(local_dof_indices);

      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        deallog << local_dof_indices[i] << " ";

      deallog << std::endl;
    }
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();

  deallog.push("2d");
  test<2>();
  deallog.pop();
}
