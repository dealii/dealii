// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check that AffineConstraints<double>.distribute() is not doing anything in a
// distributed computation for a vector that already has the entries set
// correctly

#include <deal.II/base/function.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr);

  tr.refine_global(1);
  tr.begin_active()->set_refine_flag();
  tr.execute_coarsening_and_refinement();

  DoFHandler<dim> dofh(tr);

  static FE_Q<dim> fe(1);

  dofh.distribute_dofs(fe);

  const IndexSet &owned_set    = dofh.locally_owned_dofs();
  const IndexSet  relevant_set = DoFTools::extract_locally_relevant_dofs(dofh);

  TrilinosWrappers::MPI::Vector x_ref;
  x_ref.reinit(owned_set, MPI_COMM_WORLD);
  VectorTools::interpolate(dofh, Functions::ConstantFunction<dim>(1.), x_ref);

  TrilinosWrappers::MPI::Vector x1(x_ref);

  // we have interpolated values, so
  // AffineConstraints<double>::distribute should not do
  // anything
  x1 = x_ref;
  AffineConstraints<double> cm(owned_set, relevant_set);
  DoFTools::make_hanging_node_constraints(dofh, cm);
  cm.close();
  cm.distribute(x1);

  x1 -= x_ref;
  double err = x1.linfty_norm();
  if (err > 1.0e-12)
    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
      deallog << "err:" << err << std::endl;

  // now test the same thing with a fresh vector
  // that we manually fill with ones, not by a
  // function in interpolate
  TrilinosWrappers::MPI::Vector x2(owned_set, MPI_COMM_WORLD);
  x2 = 1;
  cm.distribute(x2);
  x2 -= x_ref;
  err = x2.linfty_norm();
  if (err > 1.0e-12)
    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
      deallog << "err:" << err << std::endl;

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "OK" << std::endl;
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
    }
  else
    test<2>();
}
