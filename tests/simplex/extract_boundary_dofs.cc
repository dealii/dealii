// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Like the test mpi/extract_boundary_dofs but for simplex meshes.

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_wedge_p.h>

#include <deal.II/grid/grid_generator.h>

#include "../tests.h"



template <int dim>
void
test()
{
  parallel::shared::Triangulation<dim> tr(
    MPI_COMM_WORLD,
    Triangulation<dim>::MeshSmoothing::none,
    true,
    parallel::shared::Triangulation<dim>::Settings::partition_zorder);

  GridGenerator::subdivided_hyper_cube_with_simplices(tr, 4);

  const FE_SimplexP<dim> fe(2);
  DoFHandler<dim>        dofh(tr);
  dofh.distribute_dofs(fe);

  IndexSet boundary_dofs =
    DoFTools::extract_boundary_dofs(dofh, ComponentMask(1, true));
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    boundary_dofs.write(deallog.get_file_stream());

  // the result of extract_boundary_dofs is supposed to be a subset of the
  // locally relevant dofs, so test this
  const IndexSet relevant_set = DoFTools::extract_locally_relevant_dofs(dofh);
  boundary_dofs.subtract_set(relevant_set);
  AssertThrow(boundary_dofs.n_elements() == 0, ExcInternalError());
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      initlog();

      deallog.push("2d");
      test<2>();
      deallog.pop();

      deallog.push("3d");
      test<3>();
      deallog.pop();
    }
  else
    {
      deallog.push("2d");
      test<2>();
      deallog.pop();

      deallog.push("3d");
      test<3>();
      deallog.pop();
    }
}
