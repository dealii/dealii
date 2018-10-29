// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2018 by the deal.II authors
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



// we got an email saying that VectorTools::interpolate had trouble
// when one involved processor had no locally owned cells and would
// deadlock. but it turns out that it works. there's never anything
// wrong with having too many tests, though.

#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <sstream>

#include "../tests.h"


template <int dim>
void
test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  // create a mesh so that all but one
  // processor are empty
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tr);
  FE_Q<dim>       fe(1);
  DoFHandler<dim> dofh(tr);
  dofh.distribute_dofs(fe);

  IndexSet                      owned_set = dofh.locally_owned_dofs();
  TrilinosWrappers::MPI::Vector x;

  x.reinit(owned_set, MPI_COMM_WORLD);

  VectorTools::interpolate(dofh, Functions::ConstantFunction<dim>(1), x);
  const double norm = x.l2_norm();
  if (myid == 0)
    deallog << dofh.n_locally_owned_dofs() << ' ' << dofh.n_dofs() << std::endl
            << norm << std::endl;
}


int
main(int argc, char *argv[])
{
#ifdef DEAL_II_WITH_MPI
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());
#else
  (void)argc;
  (void)argv;
  compile_time_error;
#endif

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
      test<2>();
      test<3>();
    }
}
