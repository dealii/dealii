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
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// This is a modified copy from interpolate_01.cc.
// VectorTools::interpolate used to accidentally read from non-owned
// dofs when called with a component mask that excluded some components.
// This triggered an assertion. Test that the solution,
// namely skipping all non-local dofs that are not selected to be interpolated, works.

#include "../tests.h"
#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>
#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/distributed/tria.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <sstream>

template <int dim>
void
test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  // create a distributed mesh
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);

  // Create a system with more than one component
  FESystem<dim>   fe(FE_Q<dim>(1), 2);
  DoFHandler<dim> dofh(tr);
  dofh.distribute_dofs(fe);

  IndexSet                      owned_set = dofh.locally_owned_dofs();
  TrilinosWrappers::MPI::Vector x;

  x.reinit(owned_set, MPI_COMM_WORLD);

  // Select the first of the two components
  std::vector<bool> components(2);
  components[0] = true;

  // This failed before this test was introduced.
  // Interpolate onto the first component only.
  VectorTools::interpolate(
    dofh, Functions::ConstantFunction<dim>(1, 2), x, ComponentMask(components));

  // Integrate the difference in the first component, if everything went
  // well, this should be zero.
  IndexSet relevant_set;
  DoFTools::extract_locally_relevant_dofs(dofh, relevant_set);
  TrilinosWrappers::MPI::Vector x_rel(relevant_set, MPI_COMM_WORLD);
  x_rel = x;
  Vector<double>               error(tr.n_active_cells());
  ComponentSelectFunction<dim> right_component_select(0, 2);
  ComponentSelectFunction<dim> wrong_component_select(1, 2);

  VectorTools::integrate_difference(dofh,
                                    x_rel,
                                    Functions::ConstantFunction<dim>(1, 2),
                                    error,
                                    QMidpoint<dim>(),
                                    VectorTools::L2_norm,
                                    &right_component_select);

  double norm
    = VectorTools::compute_global_error(tr, error, VectorTools::L2_norm);

  if(myid == 0)
    deallog << dofh.n_locally_owned_dofs() << ' ' << dofh.n_dofs() << std::endl
            << "Error of interpolated component: " << norm << std::endl;

  // Integrate the difference in the second component. Since we did not interpolate
  // the function into the finite element space, this should be equal to the
  // integral of the ConstantFunction (=1) over the domain (unit square/cube).
  // Thus the integral should be one.
  VectorTools::integrate_difference(dofh,
                                    x_rel,
                                    Functions::ConstantFunction<dim>(1, 2),
                                    error,
                                    QMidpoint<dim>(),
                                    VectorTools::L2_norm,
                                    &wrong_component_select);

  norm = VectorTools::compute_global_error(tr, error, VectorTools::L2_norm);

  if(myid == 0)
    deallog << "Error of not interpolated component: " << norm << std::endl;
}

int
main(int argc, char* argv[])
{
#ifdef DEAL_II_WITH_MPI
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());
#else
  (void) argc;
  (void) argv;
  compile_time_error;
#endif

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  deallog.push(Utilities::int_to_string(myid));

  if(myid == 0)
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
