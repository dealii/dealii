// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2021 by the deal.II authors
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
    DoFTools::extract_boundary_dofs(dofh, std::vector<bool>(1, true));
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    boundary_dofs.write(deallog.get_file_stream());

  // the result of extract_boundary_dofs is supposed to be a subset of the
  // locally relevant dofs, so test this
  IndexSet relevant_set;
  DoFTools::extract_locally_relevant_dofs(dofh, relevant_set);
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
