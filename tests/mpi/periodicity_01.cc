// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2015 by the deal.II authors
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

//
// check solution for periodicity. The used test case is similar to the one
// in step-45. The right hand side is not periodic and hence the solution will
// also be not periodic unless the periodicity constraints on distributed
// triangulations are correctly implemented.
// In the 2D case we require periodicity in y-direction and in 3D both in y-
// and z-direction.
// In both cases we refine two times adaptively using the Kelly error estimator.
//

#include "../tests.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/compressed_simple_sparsity_pattern.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>

#include <fstream>
#include <iostream>

namespace Step40
{
  using namespace dealii;

  template <int dim>
  class LaplaceProblem
  {
  public:
    LaplaceProblem ();
    ~LaplaceProblem ();

    void run ();

  private:
    void setup_system ();
    void assemble_system ();
    void solve ();
    void refine_grid ();
    void get_point_value (const Point<dim> point, const int proc,
                          Vector<PetscScalar> &value) const;
    void check_periodicity (const unsigned int cycle) const;
    void output_results (const unsigned int cycle) const;

    MPI_Comm mpi_communicator;

    parallel::distributed::Triangulation<dim>   triangulation;

    DoFHandler<dim>      dof_handler;
    FE_Q<dim>            fe;

    IndexSet             locally_owned_dofs;
    IndexSet             locally_relevant_dofs;

    ConstraintMatrix     constraints;

    PETScWrappers::MPI::SparseMatrix system_matrix;
    PETScWrappers::MPI::Vector locally_relevant_solution;
    PETScWrappers::MPI::Vector system_rhs;

    ConditionalOStream                pcout;
  };

  template <int dim>
  LaplaceProblem<dim>::LaplaceProblem ()
    :
    mpi_communicator (MPI_COMM_WORLD),
    triangulation (mpi_communicator),
    dof_handler (triangulation),
    fe (2),
    pcout (Utilities::MPI::this_mpi_process(mpi_communicator)
           == 0
           ?
           deallog.get_file_stream()
           :
           std::cout,
           (Utilities::MPI::this_mpi_process(mpi_communicator)
            == 0))
  {}



  template <int dim>
  LaplaceProblem<dim>::~LaplaceProblem ()
  {
    dof_handler.clear ();
  }

  template <int dim>
  void LaplaceProblem<dim>::setup_system ()
  {
    dof_handler.distribute_dofs (fe);

    locally_owned_dofs = dof_handler.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (dof_handler,
                                             locally_relevant_dofs);
    locally_relevant_solution.reinit (locally_owned_dofs,
                                      locally_relevant_dofs,
                                      mpi_communicator);
    system_rhs.reinit (mpi_communicator,
                       dof_handler.n_dofs(),
                       dof_handler.n_locally_owned_dofs());
    system_rhs = 0;

    //Periodic Conditions
    constraints.clear ();
    constraints.reinit (locally_relevant_dofs);
    for (int i=1; i<dim; ++i)
      DoFTools::make_periodicity_constraints(dof_handler,
                                             /*b_id1*/ 2*i,
                                             /*b_id2*/ 2*i+1,
                                             /*direction*/ i,
                                             constraints);

    VectorTools::interpolate_boundary_values (dof_handler,
                                              0,
                                              ZeroFunction<dim>(),
                                              constraints);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    constraints.close ();

    CompressedSimpleSparsityPattern csp (dof_handler.n_dofs(),
                                         dof_handler.n_dofs(),
                                         locally_relevant_dofs);
    DoFTools::make_sparsity_pattern (dof_handler,
                                     csp,
                                     constraints, false);
    SparsityTools::distribute_sparsity_pattern (csp,
                                                dof_handler.n_locally_owned_dofs_per_processor(),
                                                mpi_communicator,
                                                locally_relevant_dofs);
    system_matrix.reinit (mpi_communicator,
                          csp,
                          dof_handler.n_locally_owned_dofs_per_processor(),
                          dof_handler.n_locally_owned_dofs_per_processor(),
                          Utilities::MPI::this_mpi_process(mpi_communicator));
  }

  template <int dim>
  void LaplaceProblem<dim>::assemble_system ()
  {
    const QGauss<dim>  quadrature_formula(3);

    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values    |  update_gradients |
                             update_quadrature_points |
                             update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    FullMatrix<PetscScalar>   cell_matrix (dofs_per_cell, dofs_per_cell);
    Vector<PetscScalar>       cell_rhs (dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          cell_matrix = PetscScalar();
          cell_rhs = PetscScalar();

          fe_values.reinit (cell);

          for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
            {
              PetscScalar rhs_value
                = (std::cos(2*numbers::PI*fe_values.quadrature_point(q_point)[0]) *
                   std::exp(-1*fe_values.quadrature_point(q_point)[0]) *
                   std::cos(2*numbers::PI*fe_values.quadrature_point(q_point)[1]) *
                   std::exp(-2*fe_values.quadrature_point(q_point)[1]));

              if (dim==3)
                rhs_value*=
                  std::cos(2*numbers::PI*fe_values.quadrature_point(q_point)[2]) *
                  std::exp (- 3 * fe_values.quadrature_point(q_point)[2]);

              for (unsigned int i=0; i<dofs_per_cell; ++i)
                {
                  for (unsigned int j=0; j<dofs_per_cell; ++j)
                    cell_matrix(i,j) += (fe_values.shape_grad(i,q_point) *
                                         fe_values.shape_grad(j,q_point) *
                                         fe_values.JxW(q_point));

                  cell_rhs(i) += (rhs_value *
                                  fe_values.shape_value(i,q_point) *
                                  fe_values.JxW(q_point));
                }
            }

          cell->get_dof_indices (local_dof_indices);
          constraints.distribute_local_to_global (cell_matrix,
                                                  cell_rhs,
                                                  local_dof_indices,
                                                  system_matrix,
                                                  system_rhs);
        }

    system_matrix.compress (VectorOperation::add);
    system_rhs.compress (VectorOperation::add);
  }


  template <int dim>
  void LaplaceProblem<dim>::solve ()
  {
    PETScWrappers::MPI::Vector
    completely_distributed_solution (mpi_communicator,
                                     dof_handler.n_dofs(),
                                     dof_handler.n_locally_owned_dofs());

    SolverControl solver_control (dof_handler.n_dofs(), 1e-12, false, false);

    PETScWrappers::SolverCG solver(solver_control, mpi_communicator);

#ifndef PETSC_USE_COMPLEX
    // Ask for a symmetric preconditioner by setting the first parameter in
    // AdditionalData to true.
    PETScWrappers::PreconditionBoomerAMG
    preconditioner(system_matrix,
                   PETScWrappers::PreconditionBoomerAMG::AdditionalData(true));

    solver.solve (system_matrix, completely_distributed_solution, system_rhs,
                  preconditioner);
#else
    solver.solve (system_matrix, completely_distributed_solution, system_rhs,
                  PETScWrappers::PreconditionJacobi(system_matrix));
#endif

    constraints.distribute (completely_distributed_solution);

    locally_relevant_solution = completely_distributed_solution;
  }

  template <int dim>
  void LaplaceProblem<dim>::refine_grid ()
  {
    Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
    KellyErrorEstimator<dim>::estimate (dof_handler,
                                        QGauss<dim-1>(3),
                                        typename FunctionMap<dim>::type(),
                                        locally_relevant_solution,
                                        estimated_error_per_cell);
    parallel::distributed::GridRefinement::
    refine_and_coarsen_fixed_number (triangulation,
                                     estimated_error_per_cell,
                                     0.3, 0.03);

    triangulation.execute_coarsening_and_refinement ();
  }

  template <int dim>
  void LaplaceProblem<dim>::get_point_value
  (const Point<dim> point, const int proc, Vector<PetscScalar> &value) const
  {
    typename DoFHandler<dim>::active_cell_iterator cell
      = GridTools::find_active_cell_around_point (dof_handler, point);

    if (cell->is_locally_owned())
      VectorTools::point_value (dof_handler, locally_relevant_solution,
                                point, value);

    std::vector<double> tmp (value.size());
    std::vector<double> tmp2 (value.size());
    for (unsigned int i=0; i<value.size(); ++i)
      tmp[i]=get_real_assert_zero_imag(value[i]);

    MPI_Reduce(&(tmp[0]), &(tmp2[0]), value.size(), MPI_DOUBLE,
               MPI_SUM, proc, mpi_communicator);

    for (unsigned int i=0; i<value.size(); ++i)
      value[i]=tmp2[i];
  }

  template <int dim>
  void LaplaceProblem<dim>::check_periodicity
  (const unsigned int cycle) const
  {}

  template <>
  void LaplaceProblem<2>::check_periodicity(const unsigned int cycle) const
  {
    unsigned int n_points = 2;
    for (unsigned int i = 0; i<cycle; i++)
      n_points*=2;

    //don't test exactly at the support points, since point_value is not stable there
    const double eps = 1./(16.*n_points);

    for (unsigned int i=1; i< n_points; i++)
      {
        Vector<PetscScalar> value1(1);
        Vector<PetscScalar> value2(1);

        Point <2> point1;
        point1(0)=1.*i/n_points+eps;
        point1(1)=0.;
        Point <2> point2;
        point2(0)=1.*i/n_points+eps;
        point2(1)=1.;

        get_point_value (point1, 0, value1);
        get_point_value (point2, 0, value2);

        if (Utilities::MPI::this_mpi_process(mpi_communicator)==0)
          {
            if (std::abs(value2[0]-value1[0])>1e-8)
              {
                pcout << point1 << "\t" << "fail" << std::endl;
                std::cout<<point1<< "\t" << value1[0] << std::endl;
                std::cout<<point2<< "\t" << value2[0] << std::endl;
                Assert(false, ExcInternalError());
              }
            else
              {
                pcout << point1 << "\t" << "pass" << std::endl;
              }
          }
      }
  }

  template <>
  void LaplaceProblem<3>::check_periodicity(const unsigned int cycle) const
  {
    unsigned int n_points = 2;
    for (unsigned int i = 0; i<cycle; i++)
      n_points*=2;

    //don't test exactly at the support points, since point_value is not stable there
    const double eps = 1./(16.*n_points);

    for (unsigned int i=1; i< n_points; i++)
      for (unsigned int j=1; j< n_points; j++)
        {
          Vector<PetscScalar> value1(1);
          Vector<PetscScalar> value2(1);
          Vector<PetscScalar> value3(1);
          Vector<PetscScalar> value4(1);

          Point <3> point1;
          point1(0)=1.*i/n_points+eps;
          point1(1)=1.*j/n_points+eps;
          point1(2)=0;
          Point <3> point2;
          point2(0)=1.*i/n_points+eps;
          point2(1)=1.*j/n_points+eps;
          point2(2)=1.;
          Point <3> point3;
          point3(0)=1.*i/n_points+eps;
          point3(1)=0.;
          point3(2)=1.*j/n_points+eps;
          Point <3> point4;
          point4(0)=1.*i/n_points+eps;
          point4(1)=1.;
          point4(2)=1.*j/n_points+eps;

          get_point_value (point1, 0, value1);
          get_point_value (point2, 0, value2);
          get_point_value (point3, 0, value3);
          get_point_value (point4, 0, value4);

          if (Utilities::MPI::this_mpi_process(mpi_communicator)==0)
            {
              if (std::abs(value2[0]-value1[0])>1e-8)
                {
                  pcout << point1 << "\t fail check 0" << std::endl;
                  std::cout<<point1<< "\t" << value1[0] << std::endl;
                  std::cout<<point2<< "\t" << value2[0] << std::endl;
                  Assert(false, ExcInternalError());
                }
              else
                pcout << point1 << "\t pass check 0" << std::endl;

              if (std::abs(value4[0]-value3[0])>1e-8)
                {
                  pcout << point3 << "\t fail check 1" << std::endl;
                  std::cout<<point3<< "\t" << value3[0] << std::endl;
                  std::cout<<point4<< "\t" << value4[0] << std::endl;
                  Assert(false, ExcInternalError());
                }
              else
                pcout << point3 << "\t pass check 1" << std::endl;
            }
        }
  }

  //only needed for graphical output
  template <int dim>
  void LaplaceProblem<dim>::output_results (const unsigned int cycle) const
  {
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (locally_relevant_solution, "u");

    Vector<float> subdomain (triangulation.n_active_cells());
    for (unsigned int i=0; i<subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();
    data_out.add_data_vector (subdomain, "subdomain");

    data_out.build_patches (3);

    const std::string filename = ("solution-" +
                                  Utilities::int_to_string (cycle, 2) +
                                  "." +
                                  Utilities::int_to_string
                                  (triangulation.locally_owned_subdomain(), 4));
    std::ofstream output ((filename + ".vtu").c_str());
    data_out.write_vtu (output);

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
      {
        std::vector<std::string> filenames;
        for (unsigned int i=0;
             i<Utilities::MPI::n_mpi_processes(mpi_communicator);
             ++i)
          filenames.push_back ("solution-" +
                               Utilities::int_to_string (cycle, 2) +
                               "." +
                               Utilities::int_to_string (i, 4) +
                               ".vtu");

        std::ofstream master_output ((filename + ".pvtu").c_str());
        data_out.write_pvtu_record (master_output, filenames);
      }
  }

  template <int dim>
  void LaplaceProblem<dim>::run ()
  {
    pcout << std::endl<< "Testing for dim="<<dim<<std::endl;

    const unsigned int n_cycles = 3;
    for (unsigned int cycle=0; cycle<n_cycles; ++cycle)
      {
        pcout << std::endl << "Cycle " << cycle << ':' << std::endl;

        if (cycle == 0)
          {
            std::vector<unsigned int> reps;
            reps.push_back(2);
            reps.push_back(2);
            if (dim==3)
              reps.push_back(2);

            Point<dim> p1;
            Point<dim> p2;
            for (unsigned int i=0; i<dim; ++i)
              p2(i)=1.0;

            GridGenerator::subdivided_hyper_rectangle
            (triangulation,reps,p1,p2,true);


            std::vector<GridTools::PeriodicFacePair<typename parallel::distributed::Triangulation<dim>::cell_iterator> >
            periodicity_vector;

            for (int i=1; i<dim; ++i)
              GridTools::collect_periodic_faces
              ( triangulation, /*b_id1*/ 2*i, /*b_id2*/ 2*i+1,
                /*direction*/ i, periodicity_vector);

            triangulation.add_periodicity(periodicity_vector);
            triangulation.refine_global (1);
          }
        else
          {
            refine_grid ();
          }

        setup_system ();
        assemble_system ();
        solve ();
        //output_results (cycle);
        deallog.push(Utilities::int_to_string(dof_handler.n_dofs(),5));
        check_periodicity(cycle);
        deallog.pop();
      }
  }
}

int main(int argc, char *argv[])
{
  try
    {
      using namespace dealii;
      using namespace Step40;

      Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);

      if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD)==0)
        {
          std::ofstream logfile("output");
          deallog.attach(logfile, false);
          deallog.threshold_double(1.e-10);
          {
            LaplaceProblem<2> laplace_problem;
            laplace_problem.run ();
          }
          {
            LaplaceProblem<3> laplace_problem;
            laplace_problem.run ();
          }
        }
      else
        {
          {
            LaplaceProblem<2> laplace_problem;
            laplace_problem.run ();
          }
          {
            LaplaceProblem<3> laplace_problem;
            laplace_problem.run ();
          }
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
    }

  return 0;
}
