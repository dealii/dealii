// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check FETools::back_interpolate on parallel vector

#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_tools.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>
#include <vector>

#include "../tests.h"


void
test()
{
  const unsigned int dim = 2;

  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(3);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  FE_Q<dim>       fe1(1), fe2(2);
  DoFHandler<dim> dof1(tria), dof2(tria);
  dof1.distribute_dofs(fe1);
  dof2.distribute_dofs(fe2);
  AffineConstraints<double> c1, c2;
  DoFTools::make_hanging_node_constraints(dof1, c1);
  c1.close();
  DoFTools::make_hanging_node_constraints(dof2, c2);
  c2.close();

  const IndexSet locally_relevant_dofs2 =
    DoFTools::extract_locally_relevant_dofs(dof2);

  LinearAlgebra::distributed::Vector<double> v2(dof2.locally_owned_dofs(),
                                                locally_relevant_dofs2,
                                                MPI_COMM_WORLD),
    v2_interpolated(v2);

  // set first vector to 1
  VectorTools::interpolate(dof2, Functions::ConstantFunction<dim>(1.), v2);
  for (unsigned int i = 0; i < v2.locally_owned_size(); ++i)
    Assert(v2.local_element(i) == 1., ExcInternalError());

  v2.update_ghost_values();
  FETools::back_interpolate(dof2, c2, v2, dof1, c1, v2_interpolated);
  for (unsigned int i = 0; i < v2_interpolated.locally_owned_size(); ++i)
    Assert(v2_interpolated.local_element(i) == 1., ExcInternalError());
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      initlog();
      deallog << std::setprecision(4);

      test();
    }
  else
    test();
}
