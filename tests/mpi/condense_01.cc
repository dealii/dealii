// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test AffineConstraints<double>::condense(in, out)

#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/petsc_vector.h>

#include "../tests.h"



template <int dim>
void
test()
{
  const unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(triangulation, -0.5, 0.5);

  triangulation.refine_global(1);

  /*
  for (typename Triangulation<dim>::active_cell_iterator cell =
  triangulation.begin_active (); cell != triangulation.end (); ++cell) if
  (cell->is_locally_owned ())
    {
      if(cell->center().square() < 0.3*0.3)
        cell->set_refine_flag ();
    }
  */
  typename Triangulation<dim>::active_cell_iterator cell =
    triangulation.begin_active();
  if (cell->is_locally_owned())
    cell->set_refine_flag();

  triangulation.execute_coarsening_and_refinement();
  if (myid == 0)
    deallog << "#cells = " << triangulation.n_global_active_cells()
            << std::endl;

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  AffineConstraints<PetscScalar> constraints;
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();

  const IndexSet &locally_owned_dofs = dof_handler.locally_owned_dofs();
  const IndexSet  locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof_handler);

  PETScWrappers::MPI::Vector vec(locally_owned_dofs, MPI_COMM_WORLD);
  PETScWrappers::MPI::Vector vec_ghosted(locally_owned_dofs,
                                         locally_relevant_dofs,
                                         MPI_COMM_WORLD);
  vec         = 1;
  vec_ghosted = vec;
  vec         = -1;

  constraints.condense(vec_ghosted, vec);

  vec.print(deallog.get_file_stream());

  if (myid == 0)
    deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  test<2>();
}
