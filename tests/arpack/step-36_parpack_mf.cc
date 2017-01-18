/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2016 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * This file tests the PARPACK interface for a symmetric operator taken from step-36
 * using matrix-free operators.
 *
 * We test that the computed vectors are eigenvectors and mass-orthonormal, i.e.
 * a) (A*x_i-\lambda*B*x_i).L2() == 0
 * b) x_j*B*x_i = \delta_{i,j}
 *
 */

#include "../tests.h"

#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/index_set.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/lac/iterative_inverse.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/numerics/vector_tools.h>


#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/lac/vector.h>

#include <deal.II/lac/parpack_solver.h>

#include <fstream>
#include <iostream>


const unsigned int dim = 2;

using namespace dealii;

const double eps = 1e-10;

const unsigned int fe_degree = 1;

void test ()
{
  const unsigned int global_mesh_refinement_steps = 5;
  const unsigned int number_of_eigenvalues        = 5;

  MPI_Comm mpi_communicator = MPI_COMM_WORLD;
  const unsigned int n_mpi_processes = Utilities::MPI::n_mpi_processes(mpi_communicator);
  const unsigned int this_mpi_process = Utilities::MPI::this_mpi_process(mpi_communicator);

  parallel::distributed::Triangulation<dim> triangulation (mpi_communicator);
  GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global (global_mesh_refinement_steps);


  DoFHandler<dim> dof_handler(triangulation);
  FE_Q<dim> fe(fe_degree);
  dof_handler.distribute_dofs (fe);


  IndexSet locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs (dof_handler,
                                           locally_relevant_dofs);
  ConstraintMatrix constraints;
  constraints.reinit (locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints  (dof_handler, constraints);
  VectorTools::interpolate_boundary_values (dof_handler,
                                            0,
                                            ZeroFunction<dim> (),
                                            constraints);
  constraints.close ();

  MatrixFree<dim,double> mf_data;
  {
    const QGauss<1> quad (fe_degree+1);
    typename MatrixFree<dim,double>::AdditionalData data;
    data.tasks_parallel_scheme =
      MatrixFree<dim,double>::AdditionalData::partition_color;
    data.mapping_update_flags = update_values | update_gradients | update_JxW_values;
    mf_data.reinit (dof_handler, constraints, quad, data);
  }


  std::vector<LinearAlgebra::distributed::Vector<double> > eigenfunctions;
  std::vector<double>                                      eigenvalues;
  MatrixFreeOperators::MassOperator<dim,fe_degree, fe_degree+1, 1, double> mass;
  MatrixFreeOperators::LaplaceOperator<dim,fe_degree, fe_degree+1, 1, double> laplace;
  mass.initialize(mf_data);
  laplace.initialize(mf_data);

  eigenfunctions.resize (number_of_eigenvalues);
  eigenvalues.resize (number_of_eigenvalues);
  for (unsigned int i=0; i<eigenfunctions.size (); ++i)
    mf_data.initialize_dof_vector (eigenfunctions[i]);

  // test PArpack with matrix-free
  {
    std::vector<std::complex<double> > lambda(number_of_eigenvalues);

    SolverControl solver_control     (dof_handler.n_dofs(), 1e-9,/*log_history*/false,/*log_results*/false);
    SolverControl solver_control_lin (dof_handler.n_dofs(), 1e-10,/*log_history*/false,/*log_results*/false);


    PreconditionIdentity preconditioner;
    IterativeInverse<LinearAlgebra::distributed::Vector<double> > shift_and_invert;
    shift_and_invert.initialize(laplace,preconditioner);
    shift_and_invert.solver.select("cg");
    static ReductionControl inner_control_c(dof_handler.n_dofs(), 0.0, 1.e-13);
    shift_and_invert.solver.set_control(inner_control_c);

    const unsigned int num_arnoldi_vectors = 2*eigenvalues.size() + 2;

    PArpackSolver<LinearAlgebra::distributed::Vector<double> >::AdditionalData
    additional_data(num_arnoldi_vectors,
                    PArpackSolver<LinearAlgebra::distributed::Vector<double> >::largest_magnitude,
                    true);

    PArpackSolver<LinearAlgebra::distributed::Vector<double> > eigensolver
    (solver_control,
     mpi_communicator,
     additional_data);
    eigensolver.reinit(eigenfunctions[0]);
    // make sure initial vector is orthogonal to the space due to constraints
    {
      LinearAlgebra::distributed::Vector<double> init_vector;
      mf_data.initialize_dof_vector(init_vector);
      init_vector = 1.;
      constraints.set_zero(init_vector);
      eigensolver.set_initial_vector(init_vector);
    }
    // avoid output of iterative solver:
    const unsigned int previous_depth = deallog.depth_file(0);
    eigensolver.solve (laplace,
                       mass,
                       shift_and_invert,
                       lambda,
                       eigenfunctions,
                       eigenvalues.size());
    deallog.depth_file(previous_depth);

    for (unsigned int i = 0; i < lambda.size(); i++)
      eigenvalues[i] = lambda[i].real();

    for (unsigned int i=0; i < eigenvalues.size(); i++)
      deallog << eigenvalues[i] << std::endl;

    // make sure that we have eigenvectors and they are mass-orthonormal:
    // a) (A*x_i-\lambda*B*x_i).L2() == 0
    // b) x_j*B*x_i=\delta_{ij}
    {
      const double precision = 1e-7;
      LinearAlgebra::distributed::Vector<double> Ax(eigenfunctions[0]), Bx(eigenfunctions[0]);
      for (unsigned int i=0; i < eigenfunctions.size(); ++i)
        {
          mass.vmult(Bx,eigenfunctions[i]);

          for (unsigned int j=0; j < eigenfunctions.size(); j++)
            Assert( std::abs( eigenfunctions[j] * Bx - (i==j))< precision,
                    ExcMessage("Eigenvectors " +
                               Utilities::int_to_string(i) +
                               " and " +
                               Utilities::int_to_string(j) +
                               " are not orthonormal!"));

          laplace.vmult(Ax,eigenfunctions[i]);
          Ax.add(-1.0*eigenvalues[i],Bx);
          Assert (Ax.l2_norm() < precision,
                  ExcMessage("Returned vector " +
                             Utilities::int_to_string(i) +
                             " is not an eigenvector!"));
        }
    }
  }


  dof_handler.clear ();
  deallog << "Ok"<<std::endl;
}


int main (int argc,char **argv)
{
  std::ofstream logfile("output");
  deallog.attach(logfile,/*do not print job id*/false);
  deallog.threshold_double(eps);

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      {
        test ();
      }

    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    };
}
