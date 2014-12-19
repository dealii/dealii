// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
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



// this tests the correctness of matrix free matrix-vector products when
// threads are additionally used to MPI by comparing the result with a
// non-parallel version. Otherwise same problem as matrix_vector_10.cc

#include "../tests.h"

#include "matrix_vector_mf.h"

#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/parallel_vector.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>




template <int dim, int fe_degree>
void test ()
{
  typedef double number;

  parallel::distributed::Triangulation<dim> tria (MPI_COMM_WORLD);
  GridGenerator::hyper_cube (tria);
  tria.refine_global(1);
  typename Triangulation<dim>::active_cell_iterator
  cell = tria.begin_active (),
  endc = tria.end();
  cell = tria.begin_active ();
  for (; cell!=endc; ++cell)
    if (cell->is_locally_owned())
      if (cell->center().norm()<0.2)
        cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  if (fe_degree < 2)
    tria.refine_global(2);
  else
    tria.refine_global(1);
  if (tria.begin(tria.n_levels()-1)->is_locally_owned())
    tria.begin(tria.n_levels()-1)->set_refine_flag();
  if (tria.last()->is_locally_owned())
    tria.last()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  cell = tria.begin_active ();
  for (unsigned int i=0; i<11-3*dim; ++i)
    {
      cell = tria.begin_active ();
      unsigned int counter = 0;
      for (; cell!=endc; ++cell, ++counter)
        if (cell->is_locally_owned())
          if (counter % (7-i) == 0)
            cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }

  FE_Q<dim> fe (fe_degree);
  DoFHandler<dim> dof (tria);
  dof.distribute_dofs(fe);

  IndexSet owned_set = dof.locally_owned_dofs();
  IndexSet relevant_set;
  DoFTools::extract_locally_relevant_dofs (dof, relevant_set);

  ConstraintMatrix constraints (relevant_set);
  DoFTools::make_hanging_node_constraints(dof, constraints);
  VectorTools::interpolate_boundary_values (dof, 0, ZeroFunction<dim>(),
                                            constraints);
  constraints.close();

  deallog << "Testing " << dof.get_fe().get_name() << std::endl;
  //std::cout << "Number of cells: " << tria.n_global_active_cells() << std::endl;
  //std::cout << "Number of degrees of freedom: " << dof.n_dofs() << std::endl;
  //std::cout << "Number of constraints on 0: " << constraints.n_constraints() << std::endl;

  MatrixFree<dim,number> mf_data;
  {
    const QGauss<1> quad (fe_degree+1);
    typename MatrixFree<dim,number>::AdditionalData data;
    data.mpi_communicator = MPI_COMM_WORLD;
    data.tasks_parallel_scheme =
      MatrixFree<dim,number>::AdditionalData::none;
    data.tasks_block_size = 7;
    mf_data.reinit (dof, constraints, quad, data);
  }

  MatrixFreeTest<dim,fe_degree,number,parallel::distributed::Vector<number> > mf (mf_data);
  parallel::distributed::Vector<number> in, out, ref;
  mf_data.initialize_dof_vector (in);
  out.reinit (in);
  ref.reinit (in);

  for (unsigned int i=0; i<in.local_size(); ++i)
    {
      const unsigned int glob_index =
        owned_set.nth_index_in_set (i);
      if (constraints.is_constrained(glob_index))
        continue;
      in.local_element(i) = (double)Testing::rand()/RAND_MAX;
    }

  mf.vmult (ref, in);

  for (unsigned int parallel_option = 0; parallel_option < 3; ++parallel_option)
    {
      const QGauss<1> quad (fe_degree+1);
      typename MatrixFree<dim,number>::AdditionalData data;
      data.mpi_communicator = MPI_COMM_WORLD;
      if (parallel_option == 0)
        {
          data.tasks_parallel_scheme =
            MatrixFree<dim,number>::AdditionalData::partition_partition;
          deallog << "Parallel option: partition partition" << std::endl;
        }
      else if (parallel_option == 1)
        {
          data.tasks_parallel_scheme =
            MatrixFree<dim,number>::AdditionalData::partition_color;
          deallog << "Parallel option: partition color" << std::endl;
        }
      else if (parallel_option == 2)
        {
          data.tasks_parallel_scheme =
            MatrixFree<dim,number>::AdditionalData::color;
          deallog << "Parallel option: color" << std::endl;
        }

      data.tasks_block_size = 3;
      mf_data.reinit (dof, constraints, quad, data);
      MatrixFreeTest<dim, fe_degree, number,parallel::distributed::Vector<number> > mf (mf_data);
      deallog << "Norm of difference:";

      // run 10 times to make a possible error more
      // likely to show up
      for (unsigned int run=0; run<10; ++run)
        {
          mf.vmult (out, in);
          out -= ref;
          const double diff_norm = out.linfty_norm();
          deallog << " " << diff_norm;
        }
      deallog << std::endl;
    }
  deallog << std::endl;
}


int main (int argc, char **argv)
{
  Utilities::System::MPI_InitFinalize mpi_initialization(argc, argv,
                                                         numbers::invalid_unsigned_int);

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile("output");
      deallog.attach(logfile);
      deallog << std::setprecision(4);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      deallog.push("2d");
      test<2,1>();
      test<2,2>();
      deallog.pop();

      deallog.push("3d");
      test<3,1>();
      test<3,2>();
      deallog.pop();
    }
  else
    {
      deallog.depth_console(0);
      test<2,1>();
      test<2,2>();
      test<3,1>();
      test<3,2>();
    }
}

