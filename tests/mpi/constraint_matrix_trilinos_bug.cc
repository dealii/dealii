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



// AffineConstraints<double>.distribute() produces a different result when using
// a Trilinos::Vector with ghost elements (e.g. owned vs. active), which is a
// bug. Now distribute() throws an Exception when called with a Vector with
// ghost elements. check that this happens.

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

#include "../tests.h"



template <int dim>
void
test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr);

  tr.refine_global(1);
  for (unsigned int step = 0; step < 5; ++step)
    {
      typename Triangulation<dim>::active_cell_iterator cell =
                                                          tr.begin_active(),
                                                        endc = tr.end();

      for (; cell != endc; ++cell)
        if (Testing::rand() % 42 == 1)
          cell->set_refine_flag();

      tr.execute_coarsening_and_refinement();
    }
  DoFHandler<dim> dofh(tr);

  static FESystem<dim> fe(FE_Q<dim>(1 + 1), dim, FE_Q<dim>(1), 1);

  dofh.distribute_dofs(fe);

  const IndexSet &owned_set = dofh.locally_owned_dofs();

  const IndexSet relevant_set = DoFTools::extract_locally_relevant_dofs(dofh);

  TrilinosWrappers::MPI::Vector x;
  x.reinit(owned_set, MPI_COMM_WORLD);
  x = 2.0;

  TrilinosWrappers::MPI::Vector x_rel;
  x_rel.reinit(relevant_set, MPI_COMM_WORLD);

  AffineConstraints<double> cm(owned_set, relevant_set);
  DoFTools::make_hanging_node_constraints(dofh, cm);
  ComponentMask velocity_mask(dim + 1, true);

  velocity_mask.set(dim, false);

  VectorTools::interpolate_boundary_values(
    dofh, 0, Functions::ZeroFunction<dim>(dim + 1), cm, velocity_mask);

  cm.close();

  TrilinosWrappers::MPI::Vector x_test;
  x_test.reinit(x_rel);

  x_test = x;

  bool throwing = false;
  deal_II_exceptions::disable_abort_on_exception();
  try
    {
      cm.distribute(x_test);
    }
  catch (const ExceptionBase &e)
    {
      if (myid == 0)
        deallog << "Exception: " << e.get_exc_name() << std::endl;
      throwing = true;
    }
  Assert(throwing, ExcInternalError());

  cm.distribute(x);
  x_rel = x;

  // l2_norm() not possible for ghosted vectors...
  // double a=0;//x_test.l2_norm();
  // double b=0;//x_rel.l2_norm();

  /*    if (myid==0)
        deallog << a << " vs " << b << std::endl;
  */
  /*    Assert (x_test.l2_norm() == x_rel.l2_norm(),
        ExcInternalError());
  */
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
