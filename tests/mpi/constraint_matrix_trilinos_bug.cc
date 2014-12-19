// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// ConstraintMatrix.distribute() produces a different result when using a
// Trilinos::Vector with ghost elements (e.g. owned vs. active), which is a
// bug. Now distribute() throws an Exception when called with a Vector with
// ghost elements. check that this happens.

#include "../tests.h"
#include "coarse_grid_common.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/grid/filtered_iterator.h>

#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <fstream>
#include <sstream>



template<int dim>
void test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr);

  tr.refine_global (1);
  for (unsigned int step=0; step<5; ++step)
    {
      typename Triangulation<dim>::active_cell_iterator
      cell = tr.begin_active(),
      endc = tr.end();

      for (; cell!=endc; ++cell)
        if (Testing::rand()%42==1)
          cell->set_refine_flag ();

      tr.execute_coarsening_and_refinement ();
    }
  DoFHandler<dim> dofh(tr);

  static FESystem<dim> fe (FE_Q<dim>(1+1), dim,
                           FE_Q<dim>(1), 1);

  dofh.distribute_dofs (fe);

  IndexSet owned_set = dofh.locally_owned_dofs();

  IndexSet relevant_set;
  DoFTools::extract_locally_relevant_dofs (dofh, relevant_set);

  TrilinosWrappers::MPI::Vector x;
  x.reinit(owned_set, MPI_COMM_WORLD);
  x=2.0;

  TrilinosWrappers::MPI::Vector x_rel;
  x_rel.reinit(relevant_set, MPI_COMM_WORLD);

  ConstraintMatrix cm(relevant_set);
  DoFTools::make_hanging_node_constraints (dofh, cm);
  std::vector<bool> velocity_mask (dim+1, true);

  velocity_mask[dim] = false;

  VectorTools::interpolate_boundary_values (dofh,
                                            0,
                                            ZeroFunction<dim>(dim+1),
                                            cm,
                                            velocity_mask);

  cm.close ();

  TrilinosWrappers::MPI::Vector x_test;
  x_test.reinit(x_rel);

  x_test=x;

  bool throwing=false;
  deal_II_exceptions::disable_abort_on_exception();
  try
    {
      cm.distribute(x_test);
    }
  catch (const ExceptionBase &e)
    {
      if (myid==0)
        deallog << "Exception: " << e.get_exc_name() << std::endl;
      throwing=true;
    }
  Assert(throwing, ExcInternalError());

  cm.distribute(x);
  x_rel = x;

  //l2_norm() not possible for ghosted vectors...
  //double a=0;//x_test.l2_norm();
  //double b=0;//x_rel.l2_norm();

  /*    if (myid==0)
        deallog << a << " vs " << b << std::endl;
  */
  /*    Assert (x_test.l2_norm() == x_rel.l2_norm(),
        ExcInternalError());
  */
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile("output");
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      deallog.push("2d");
      test<2>();
      deallog.pop();
    }
  else
    test<2>();

}
