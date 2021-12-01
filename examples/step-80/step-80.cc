/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2021 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
 *
 * Authors: Hongfeng Ma, 2021.
 */


// @sect3{Include files}

// All the include files have already been discussed in previous tutorials.
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

// the following files are for petsc, if you want to use
// pesct instead of trilinos, you can uncomment the following

/*
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>
*/

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_precondition.h>

#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/generic_linear_algebra.h>

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/index_set.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>


#include <deal.II/base/logstream.h>
// this library is deprecated. use affine_constraints.h instead
//#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/base/timer.h>

// grid refine
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>

// remember to use the namespace dealii before
// defining a class that is inheritaged from Function<dim>
using namespace dealii;

// In general, it is more clear to separate boundary and initial condictions,
// as well as the right hand side function from a file holding all the things.
// To do so, in this work, the globalPara.h file defines physical constants, 
// laser parameters, and heat characteristics of materials involved. Boundary
// and initial conditions are defined in boundaryInit.h. The rightHandSide.h
// defines the heat source, which in this work is a moving Gaussian beam.

#ifndef GLOBAL_PARA
#define GLOBAL_PARA
#include "./globalPara.h"
#include "./boundaryInit.h"
#include "./rightHandSide.h"
#endif

// Now the main class is defined as following
template <int dim>
class LaserHeating
{
public:
  LaserHeating ();
  ~LaserHeating ();
  void run ();

private:
  void make_grid ();
  void setup_system();

  void assemble_system_matrix_init (double time_step);
  void dynamic_assemble_rhs_T  (double time, double time_step);

  void solve_T ();

  void refine_mesh();
  void output_results (int output_num) const;
  
  MPI_Comm                                    mpi_communicator;

  parallel::distributed::Triangulation<dim>   triangulation;
  FE_Q<dim>            fe;
  DoFHandler<dim>      dof_handler;

  AffineConstraints<double>     constraints_T;


  // system_matrix
  TrilinosWrappers::SparseMatrix system_matrix_T;

  // for storing left matrix
  TrilinosWrappers::SparseMatrix left_system_matrix_T;

  // for storing right matrix
  TrilinosWrappers::SparseMatrix right_system_matrix_T;

  // System_rhs, only locally owned cells
  TrilinosWrappers::MPI::Vector       system_rhs_T;

  // Solutions
  // Old Solutions with ghost cells, for output
  TrilinosWrappers::MPI::Vector       old_solution_T; 

  // Old Solutions only with locally owned cells     
  TrilinosWrappers::MPI::Vector       old_solution_T_cal; 

  // New Solutions only with locally owned cells      
  TrilinosWrappers::MPI::Vector       new_solution_T; 

  // Dynamic assembling of the righthandside terms
  TrilinosWrappers::MPI::Vector       dynamic_rhs_T;


  IndexSet              locally_owned_dofs;
  IndexSet              locally_relevant_dofs;

  ConditionalOStream    pcout;
  TimerOutput           computing_timer;

  double theta;



};

// the constructor
template <int dim>
LaserHeating<dim>::LaserHeating ()
  :
  mpi_communicator (MPI_COMM_WORLD),
  triangulation (mpi_communicator),
  fe (1),
  dof_handler (triangulation),
  pcout (std::cout,Utilities::MPI::this_mpi_process(mpi_communicator) == 0),
  computing_timer (mpi_communicator, pcout, TimerOutput::summary, TimerOutput::wall_times),
  theta(0.5)
{}

// the destructor
template <int dim>
LaserHeating<dim>::~LaserHeating ()
{
    dof_handler.clear();
}

// @sect4{LaserHeating::make_grid}
// make grid by importing msh file, and rescale
template <int dim>
void LaserHeating<dim>::make_grid ()
{
  TimerOutput::Scope t(computing_timer,"make_grid()");

  GridIn<dim> grid_in;
  grid_in.attach_triangulation (triangulation);
  std::ifstream input_file ("geometry.msh");
  grid_in.read_msh (input_file);
  GridTools::scale (1e-6,triangulation);

  pcout << "   Number of active cells: "
            << triangulation.n_global_active_cells()
            << std::endl
            << "   Total number of cells: "
            << triangulation.n_cells()
            << std::endl;
}

// @sect4{LaserHeating::setup_system}
// initialization
template <int dim>
void LaserHeating<dim>::setup_system ()
{

  TimerOutput::Scope t(computing_timer,"setup_system()");

  dof_handler.distribute_dofs (fe);

  pcout << "   Number of degrees of freedom: "
            << dof_handler.n_dofs()
            << std::endl;

  locally_owned_dofs = dof_handler.locally_owned_dofs();
  DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);

  // we want to output solution, so here should have ghost cells
  old_solution_T.reinit     (locally_owned_dofs,locally_relevant_dofs,mpi_communicator);

  // locally owned cells
  old_solution_T_cal.reinit (locally_owned_dofs,mpi_communicator);
  new_solution_T.reinit     (locally_owned_dofs,mpi_communicator);
  dynamic_rhs_T.reinit      (locally_owned_dofs,mpi_communicator);
  system_rhs_T.reinit       (locally_owned_dofs,mpi_communicator);

  constraints_T.clear();
  constraints_T.reinit (locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints (dof_handler, constraints_T);
/*
  VectorTools::interpolate_boundary_values (dof_handler,
                                            BOUNDARY_NUM,
                                            BoundaryValues<dim>(),
                                            constraints_T);                                      
*/
  constraints_T.close();

  DynamicSparsityPattern dsp_T(locally_relevant_dofs);
  DoFTools::make_sparsity_pattern (dof_handler, dsp_T,constraints_T,false);
  SparsityTools::distribute_sparsity_pattern (dsp_T,
                                              dof_handler.n_locally_owned_dofs_per_processor(),
                                              mpi_communicator,
                                              locally_relevant_dofs);

  left_system_matrix_T.reinit  (locally_owned_dofs,locally_owned_dofs,dsp_T,mpi_communicator);
  right_system_matrix_T.reinit (locally_owned_dofs,locally_owned_dofs,dsp_T,mpi_communicator);
  system_matrix_T.reinit       (locally_owned_dofs,locally_owned_dofs,dsp_T,mpi_communicator);

}


// @sect4{LaserHeating::assemble_system_matrix}

template <int dim>
void LaserHeating<dim>::assemble_system_matrix_init (double time_step)
{

  TimerOutput::Scope t(computing_timer,"assemble_system_matrix_init()");

  QGauss<dim>  quadrature_formula(2);


  const InitialValues<dim>   initial_value_func_T;


  const RhoC<dim>       rho_C_fun_T;
  const K_T<dim>        k_fun_T;
 

  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_values   | update_gradients |
                           update_quadrature_points | update_JxW_values);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();


  FullMatrix<double>   local_init_matrix    (dofs_per_cell, dofs_per_cell);
  Vector<double>       local_init_T_rhs     (dofs_per_cell);

  FullMatrix<double>   local_rho_c_T_matrix (dofs_per_cell, dofs_per_cell);
  FullMatrix<double>   local_k_T_matrix     (dofs_per_cell, dofs_per_cell);


  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  system_matrix_T = 0;
  left_system_matrix_T = 0;
  right_system_matrix_T = 0;

  system_rhs_T = 0;

  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();

  for (; cell!=endc; ++cell)
      if(cell->is_locally_owned())
    {

      fe_values.reinit (cell);
      local_init_matrix = 0;

      local_rho_c_T_matrix = 0;
      local_k_T_matrix = 0;

      local_init_T_rhs = 0;

      for (unsigned int q=0; q<n_q_points; ++q)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            const Tensor<1,dim>     div_phi_i_u = fe_values.shape_grad (i,q);
            const double            phi_i_u = fe_values.shape_value (i,q);

            for (unsigned int j=0; j<dofs_per_cell; ++j)
            {

              const Tensor<1,dim>       div_phi_j_u = fe_values.shape_grad(j,q);
              const double              phi_j_u = fe_values.shape_value (j,q);

              local_init_matrix(i,j) +=    (phi_i_u *
                                            phi_j_u *
                                            fe_values.JxW (q));

              local_rho_c_T_matrix(i,j) += (rho_C_fun_T.value(fe_values.quadrature_point(q)) *
                                            phi_i_u *
                                            phi_j_u *
                                            fe_values.JxW (q))
                                            +
                                            time_step * (theta) *  
                                            (k_fun_T.value(fe_values.quadrature_point(q)) *
                                            div_phi_i_u *
                                            div_phi_j_u *
                                            fe_values.JxW (q));

              local_k_T_matrix(i,j) += (rho_C_fun_T.value(fe_values.quadrature_point(q)) *
                                            phi_i_u *
                                            phi_j_u *
                                            fe_values.JxW (q))
                                            -
                                            time_step * (1.0-theta) * 
                                            (k_fun_T.value(fe_values.quadrature_point(q)) *
                                            div_phi_i_u *
                                            div_phi_j_u *
                                            fe_values.JxW (q));
             
            }

            local_init_T_rhs(i) += (phi_i_u *
                                    initial_value_func_T.value (fe_values.quadrature_point (q)) *
                                    fe_values.JxW (q));

          }

      cell->get_dof_indices (local_dof_indices);

      // copy to system_matrix_T and system_rhs_T for projecting initial values
      constraints_T.distribute_local_to_global(local_init_matrix,
                                             local_init_T_rhs,
                                             local_dof_indices,
                                             system_matrix_T,
                                             system_rhs_T);


      // store M + dt*theta*A as the left_system_matrix

      constraints_T.distribute_local_to_global(local_rho_c_T_matrix,
                                             local_dof_indices,
                                             left_system_matrix_T);

      // store M - dt*(1-theta)*A as the right_system_matrix
      constraints_T.distribute_local_to_global(local_k_T_matrix,
                                             local_dof_indices,
                                             right_system_matrix_T);


    }

  system_matrix_T.compress(VectorOperation::add);
  left_system_matrix_T.compress(VectorOperation::add);
  right_system_matrix_T.compress(VectorOperation::add);
  system_rhs_T.compress(VectorOperation::add);

}


// @sect4{LaserHeating::dynamic_assemble_rhs_T}
// The right hand side is assembled each time during running, which is necessary as
// the laser source is moving. To separate the heat source and the right hand 
// side assembling, the right hand side function is defined as RightHandside<dim>.
template <int dim>
void LaserHeating<dim>::dynamic_assemble_rhs_T (double time, double time_step)
{

  TimerOutput::Scope t(computing_timer,"assemble_rhs_T()");

  QGauss<dim>  quadrature_formula(2);

  RightHandside<dim> rhs_func_T_1;
  rhs_func_T_1.set_time(time);

  RightHandside<dim> rhs_func_T_2;
  rhs_func_T_2.set_time(time-time_step);

  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_values   |
                           update_quadrature_points | update_JxW_values);


  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  Vector<double>       local_rhs_vector_T (dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  std::vector<double>   Rnp_cal_assemble (n_q_points);


  dynamic_rhs_T = 0 ;

  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();

  for (; cell!=endc; ++cell)
      if(cell->is_locally_owned())
    {
      fe_values.reinit (cell);
      local_rhs_vector_T = 0;

      for (unsigned int q=0; q<n_q_points; ++q)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            const double            phi_i_u = fe_values.shape_value (i,q);


            local_rhs_vector_T(i) +=     time_step * theta *
                                         (phi_i_u *
                                         rhs_func_T_1.value_v2 (fe_values.quadrature_point (q)) *
                                         fe_values.JxW (q))
                                         +
                                         time_step * (1.0 - theta) *
                                         (phi_i_u *
                                         rhs_func_T_2.value_v2 (fe_values.quadrature_point (q)) *
                                         fe_values.JxW (q));
          }

      cell->get_dof_indices (local_dof_indices);


      constraints_T.distribute_local_to_global(local_rhs_vector_T,
                                             local_dof_indices,
                                             dynamic_rhs_T);

    }

  dynamic_rhs_T.compress(VectorOperation::add);


}


// @sect4{LaserHeating::solve}
// Solving the equation is direct. Recall that we have defined several matrices and vectors,
// to avoid ambiguous, here, we only use system_matrix_T as the system matrix, system_rhs_T
// as the system right hand side. The vector completely_distributed_solution is used to store
// the obtained solution.
template <int dim>
void LaserHeating<dim>::solve_T ()
{
  TimerOutput::Scope t(computing_timer,"solve_T()");

  TrilinosWrappers::MPI::Vector   completely_distributed_solution (locally_owned_dofs,mpi_communicator);
  SolverControl     solver_control (1*system_rhs_T.size(),1e-12*system_rhs_T.l2_norm(),true);

  TrilinosWrappers::SolverCG      solver (solver_control);

  TrilinosWrappers::PreconditionAMG  preconditioner;
  TrilinosWrappers::PreconditionAMG::AdditionalData  data;

  preconditioner.initialize(system_matrix_T,data);

  solver.solve (system_matrix_T,completely_distributed_solution,system_rhs_T,preconditioner);


  // Print the number of iterations by hand.

  pcout     << "   " << solver_control.last_step()
            << " CG iterations needed to obtain convergence." << std::endl
            << "\t initial convergence value = " << solver_control.initial_value() << std::endl
            << "\t final convergence value = " << solver_control.last_value() << std::endl
            << std::endl;

  constraints_T.distribute (completely_distributed_solution);
  new_solution_T = completely_distributed_solution;

}


// @sect4{LaserHeating::refine_mesh}

template <int dim>
void LaserHeating<dim>::refine_mesh()
{

  TimerOutput::Scope t(computing_timer,"refine_mesh_at_beginning");


  QGauss<dim>  quadrature_formula(2);

  FEValues<dim> fe_values (fe, quadrature_formula,update_quadrature_points);


// only refine mesh 5um above the TiO2 and glass interface 

  for (typename Triangulation<dim>::active_cell_iterator
          cell = triangulation.begin_active();
          cell != triangulation.end(); ++cell)
      if(cell->is_locally_owned())
      {
          fe_values.reinit(cell);
          if(std::abs(fe_values.quadrature_point(0)[1]) <= global_film_thickness+5e-6)
          {
              cell->set_refine_flag();
          }
          else
          {
              cell->clear_refine_flag();
          }
      }
  triangulation.execute_coarsening_and_refinement();
  
}


// @sect4{LaserHeatingh::output_results}

template <int dim>
void LaserHeating<dim>::output_results (int output_num) const
 {

   DataOut<dim> data_out;

   data_out.attach_dof_handler (dof_handler);
   data_out.add_data_vector (old_solution_T,         "T");

   // the length of output numbering

   int step_N = 7;

   data_out.build_patches ();

   const std::string filename = ("solution-" +
                                Utilities::int_to_string (output_num, step_N) +
                                "." +
                                Utilities::int_to_string (triangulation.locally_owned_subdomain(),4) +
                                ".vtu");
   std::ofstream output (filename.c_str());
   data_out.write_vtu (output);


   // output the overall solution

   if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
   {
       std::vector<std::string> filenames;
       for (unsigned int i=0;i<Utilities::MPI::n_mpi_processes(mpi_communicator);++i)
           filenames.push_back ("solution-" +
                                Utilities::int_to_string (output_num, step_N) +
                                "." +
                                Utilities::int_to_string (i,4) +
                                ".vtu");
       std::ofstream master_output (("solution-" + 
                                    Utilities::int_to_string (output_num,step_N)+
                                    ".pvtu").c_str());
       data_out.write_pvtu_record (master_output,filenames);

   }

 }




// @sect4{LaserHeating::run}

// This is the function which has the top-level control over everything. Apart
// from one line of additional output, it is the same as for the previous
// example.
template <int dim>
void LaserHeating<dim>::run ()
{
    pcout << "Solving problem in " << dim << " space dimensions." << std::endl;

    
    make_grid();
    refine_mesh();
    setup_system ();
    assemble_system_matrix_init (global_simulation_time_step);

    // projection of initial conditions by solving.
    // solution stored in new_solution_T; 

    solve_T ();

    old_solution_T        = new_solution_T;
    old_solution_T_cal    = new_solution_T;

    // reinitialization

    system_matrix_T = 0;
    system_rhs_T = 0;


    double time_step = global_simulation_time_step;
    double time = 0;
    int timestep_number = 0;

    // output initial values; need ghost cells 
    output_results (0);


      while(time < global_simulation_end_time)
      {

        time += time_step;
        timestep_number ++;

        pcout << "Time step " << timestep_number
                  << " at t=" << time
                  << " time_step = " << time_step
                  << std::endl;

        // the dynamic solving part
        {

            right_system_matrix_T.vmult(system_rhs_T,old_solution_T_cal);

            dynamic_assemble_rhs_T (time,time_step);
            system_rhs_T.add(1,dynamic_rhs_T);

            system_matrix_T.copy_from (left_system_matrix_T);

            {
                BoundaryValues<dim> boundary_values_function;
                std::map<types::global_dof_index,double> boundary_values;

                VectorTools::interpolate_boundary_values (dof_handler,
                                                            BOUNDARY_NUM,
                                                            boundary_values_function,
                                                            boundary_values);

                MatrixTools::apply_boundary_values (boundary_values,
                                                        system_matrix_T,
                                                        new_solution_T,
                                                        system_rhs_T,
                                                        false);
            }


            solve_T ();

            // old_solution_T is used for output, holding ghost cells
            // old_solution_T_cal is used for calculation, holding only
            // locally owned cells.
            old_solution_T          = new_solution_T;
            old_solution_T_cal      = new_solution_T;

            if (Utilities::MPI::n_mpi_processes(mpi_communicator) <= 96 && (timestep_number % 50  == 0 ))
            {
                TimerOutput::Scope t(computing_timer,"output");
                output_results (timestep_number);
            }

            computing_timer.print_summary ();
            computing_timer.reset();
            
            pcout << std::endl;

        }
      }

}


// @sect3{The <code>main</code> function}

int main (int argc, char *argv[])
{
    try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      LaserHeating<2> laserHeating_2d;
      laserHeating_2d.run ();


    }
    catch (std::exception &exc)
    {
        std::cerr << std::endl << std::endl
            << "--------------------------------------------"
            <<std::endl;
    }


  return 0;
}
