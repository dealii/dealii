// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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



// Tests passing of settings in TrilinosWrappers::SolverBase::AdditionalData

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"

namespace LA = dealii::LinearAlgebraTrilinos;


class Test_Solver_Output
{
public:
  Test_Solver_Output();

  ~Test_Solver_Output();

  void
  run();

private:
  void
  make_grid();

  void
  setup_system();

  void
  assemble_system();

  void
  solve_base();
  void
  solve_cg();
  void
  solve_cgs();
  void
  solve_gmres();
  void
  solve_bicgstab();
  void
  solve_tfqmr();

  void
  refine_grid();

  void
  output(unsigned int cycle);

  MPI_Comm           mpi_comm;
  const unsigned int n_mpi_proc;
  const unsigned int this_mpi_proc;

  parallel::distributed::Triangulation<2> triangulation;

  DoFHandler<2> dof_handler;
  FE_Q<2>       fe;

  IndexSet locally_owned_dofs;
  IndexSet locally_relevant_dofs;

  AffineConstraints<double> constraints;

  LA::MPI::SparseMatrix system_matrix;
  LA::MPI::Vector       locally_relevant_solution;
  LA::MPI::Vector       system_rhs;

  ConditionalOStream pcout;
  TimerOutput        timer;
};

Test_Solver_Output::Test_Solver_Output()
  : mpi_comm(MPI_COMM_WORLD)
  , n_mpi_proc(Utilities::MPI::n_mpi_processes(mpi_comm))
  , this_mpi_proc(Utilities::MPI::this_mpi_process(mpi_comm))
  , triangulation(mpi_comm,
                  typename Triangulation<2>::MeshSmoothing(
                    Triangulation<2>::smoothing_on_refinement |
                    Triangulation<2>::smoothing_on_coarsening))
  , dof_handler(triangulation)
  , fe(1)
  , pcout(std::cout, (Utilities::MPI::this_mpi_process(mpi_comm) == 0))
  ,
  // pcout(deallog.get_file_stream(),
  //       (Utilities::MPI::this_mpi_process(mpi_comm) == 0)),
  timer(mpi_comm, pcout, TimerOutput::never, TimerOutput::wall_times)
{}

Test_Solver_Output::~Test_Solver_Output()
{
  dof_handler.clear();
}

void
Test_Solver_Output::run()
{
  const unsigned int n_cycles = 2;
  for (unsigned int cycle = 0; cycle < n_cycles; ++cycle)
    {
      pcout << "   Cycle: " << cycle << std::endl;
      if (cycle == 0)
        {
          make_grid();
        }
      else
        {
          refine_grid();
        }

      setup_system();

      pcout << "   Number of active cells:       "
            << triangulation.n_global_active_cells() << std::endl
            << "   Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;

      assemble_system();
      solve_base();
      solve_cg();
      solve_cgs();
      solve_gmres();
      solve_bicgstab();
      solve_tfqmr();

      // {
      //   TimerOutput::Scope t(timer, "output");
      //   output(cycle);
      // }
      //
      // timer.print_summary();
      // timer.reset();

      pcout << std::endl << std::endl << std::endl;
    }
}

void
Test_Solver_Output::make_grid()
{
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(5);
}

void
Test_Solver_Output::setup_system()
{
  TimerOutput::Scope t(timer, "setup");

  dof_handler.distribute_dofs(fe);
  locally_owned_dofs = dof_handler.locally_owned_dofs();
  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

  locally_relevant_solution.reinit(locally_owned_dofs,
                                   locally_relevant_dofs,
                                   mpi_comm);
  system_rhs.reinit(locally_owned_dofs, mpi_comm);

  constraints.clear();
  constraints.reinit(locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           Functions::ZeroFunction<2>(),
                                           constraints);
  constraints.close();

  DynamicSparsityPattern dsp(locally_relevant_dofs);

  DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);

  SparsityTools::distribute_sparsity_pattern(dsp,
                                             locally_owned_dofs,
                                             mpi_comm,
                                             locally_relevant_dofs);

  system_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp, mpi_comm);
}

void
Test_Solver_Output::assemble_system()
{
  TimerOutput::Scope t(timer, "assembly");

  const QGauss<2> quadrature_formula(3);

  FEValues<2> fe_values(fe,
                        quadrature_formula,
                        update_values | update_gradients |
                          update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  for (auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          cell_matrix = 0;
          cell_rhs    = 0;

          fe_values.reinit(cell);

          for (unsigned int qp = 0; qp < n_q_points; ++qp)
            {
              const double rhs_value =
                (fe_values.quadrature_point(qp)[1] >
                     0.5 + 0.25 * std::sin(4.0 * numbers::PI *
                                           fe_values.quadrature_point(qp)[0]) ?
                   1 :
                   -1);
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      cell_matrix(i, j) +=
                        (fe_values.shape_grad(i, qp) *
                         fe_values.shape_grad(j, qp) * fe_values.JxW(qp));
                    }
                  cell_rhs(i) += (rhs_value * fe_values.shape_value(i, qp) *
                                  fe_values.JxW(qp));
                }
            }
          cell->get_dof_indices(local_dof_indices);
          constraints.distribute_local_to_global(cell_matrix,
                                                 cell_rhs,
                                                 local_dof_indices,
                                                 system_matrix,
                                                 system_rhs);
        }
    }
  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);
}

void
Test_Solver_Output::solve_base()
{
  TimerOutput::Scope t(timer, "solve_base");
  pcout << "Solving using SolverBase" << std::endl;

  LA::MPI::Vector completely_distributed_solution(locally_owned_dofs, mpi_comm);

  LA::MPI::PreconditionAMG                 prec;
  LA::MPI::PreconditionAMG::AdditionalData amgdata;
  prec.initialize(system_matrix, amgdata);

  SolverControl                                solver_control(100, 1e-12);
  TrilinosWrappers::SolverBase::AdditionalData solver_data(true);
  TrilinosWrappers::SolverBase solver(TrilinosWrappers::SolverBase::cg,
                                      solver_control,
                                      solver_data);
  solver.solve(system_matrix,
               completely_distributed_solution,
               system_rhs,
               prec);

  pcout << "   Solved in " << solver_control.last_step() << " iterations."
        << std::endl;

  constraints.distribute(completely_distributed_solution);
  locally_relevant_solution = completely_distributed_solution;
}

void
Test_Solver_Output::solve_cg()
{
  TimerOutput::Scope t(timer, "solve_cg");
  pcout << "Solving using SolverCG" << std::endl;

  LA::MPI::Vector completely_distributed_solution(locally_owned_dofs, mpi_comm);

  LA::MPI::PreconditionAMG                 prec;
  LA::MPI::PreconditionAMG::AdditionalData amgdata;
  prec.initialize(system_matrix, amgdata);

  SolverControl                solver_control(100, 1e-12);
  LA::SolverCG::AdditionalData solver_data(true);
  LA::SolverCG                 solver(solver_control, solver_data);
  solver.solve(system_matrix,
               completely_distributed_solution,
               system_rhs,
               prec);

  pcout << "   Solved in " << solver_control.last_step() << " iterations."
        << std::endl;

  constraints.distribute(completely_distributed_solution);
  locally_relevant_solution = completely_distributed_solution;
}

void
Test_Solver_Output::solve_cgs()
{
  TimerOutput::Scope t(timer, "solve_cgs");
  pcout << "Solving using SolverCGS" << std::endl;

  LA::MPI::Vector completely_distributed_solution(locally_owned_dofs, mpi_comm);

  LA::MPI::PreconditionAMG                 prec;
  LA::MPI::PreconditionAMG::AdditionalData amgdata;
  prec.initialize(system_matrix, amgdata);

  SolverControl                               solver_control(100, 1e-12);
  TrilinosWrappers::SolverCGS::AdditionalData solver_data(true);
  TrilinosWrappers::SolverCGS solver(solver_control, solver_data);
  solver.solve(system_matrix,
               completely_distributed_solution,
               system_rhs,
               prec);

  pcout << "   Solved in " << solver_control.last_step() << " iterations."
        << std::endl;

  constraints.distribute(completely_distributed_solution);
  locally_relevant_solution = completely_distributed_solution;
}

void
Test_Solver_Output::solve_gmres()
{
  TimerOutput::Scope t(timer, "solve_gmres");
  pcout << "Solving using SolverGMRES" << std::endl;

  LA::MPI::Vector completely_distributed_solution(locally_owned_dofs, mpi_comm);

  LA::MPI::PreconditionAMG                 prec;
  LA::MPI::PreconditionAMG::AdditionalData amgdata;
  prec.initialize(system_matrix, amgdata);

  SolverControl                   solver_control(100, 1e-12);
  LA::SolverGMRES::AdditionalData solver_data(true, 25);
  LA::SolverGMRES                 solver(solver_control, solver_data);
  solver.solve(system_matrix,
               completely_distributed_solution,
               system_rhs,
               prec);

  pcout << "   Solved in " << solver_control.last_step() << " iterations."
        << std::endl;

  constraints.distribute(completely_distributed_solution);
  locally_relevant_solution = completely_distributed_solution;
}

void
Test_Solver_Output::solve_bicgstab()
{
  TimerOutput::Scope t(timer, "solve_bicgstab");
  pcout << "Solving using SolverBicgstab" << std::endl;

  LA::MPI::Vector completely_distributed_solution(locally_owned_dofs, mpi_comm);

  LA::MPI::PreconditionAMG                 prec;
  LA::MPI::PreconditionAMG::AdditionalData amgdata;
  prec.initialize(system_matrix, amgdata);

  SolverControl                                    solver_control(100, 1e-12);
  TrilinosWrappers::SolverBicgstab::AdditionalData solver_data(true);
  TrilinosWrappers::SolverBicgstab solver(solver_control, solver_data);
  solver.solve(system_matrix,
               completely_distributed_solution,
               system_rhs,
               prec);

  pcout << "   Solved in " << solver_control.last_step() << " iterations."
        << std::endl;

  constraints.distribute(completely_distributed_solution);
  locally_relevant_solution = completely_distributed_solution;
}

void
Test_Solver_Output::solve_tfqmr()
{
  TimerOutput::Scope t(timer, "solve_tfqmr");
  pcout << "Solving using SolverTFQMR" << std::endl;

  LA::MPI::Vector completely_distributed_solution(locally_owned_dofs, mpi_comm);

  LA::MPI::PreconditionAMG                 prec;
  LA::MPI::PreconditionAMG::AdditionalData amgdata;
  prec.initialize(system_matrix, amgdata);

  SolverControl                                 solver_control(100, 1e-12);
  TrilinosWrappers::SolverTFQMR::AdditionalData solver_data(true);
  TrilinosWrappers::SolverTFQMR solver(solver_control, solver_data);
  solver.solve(system_matrix,
               completely_distributed_solution,
               system_rhs,
               prec);

  pcout << "   Solved in " << solver_control.last_step() << " iterations."
        << std::endl;

  constraints.distribute(completely_distributed_solution);
  locally_relevant_solution = completely_distributed_solution;
}

void
Test_Solver_Output::output(unsigned int cycle)
{
  DataOut<2> data_out;

  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(locally_relevant_solution, "u");

  Vector<float> subdomain(triangulation.n_active_cells());
  for (unsigned int i = 0; i < subdomain.size(); ++i)
    subdomain(i) = triangulation.locally_owned_subdomain();
  data_out.add_data_vector(subdomain, "subdomain");
  data_out.build_patches();

  const std::string filename =
    ("solution-" + Utilities::int_to_string(cycle, 2) + "." +
     Utilities::int_to_string(triangulation.locally_owned_subdomain(), 4));
  std::ofstream output((filename + ".vtu").c_str());
  data_out.write_vtu(output);

  if (Utilities::MPI::this_mpi_process(mpi_comm) == 0)
    {
      std::vector<std::string> filenames;
      for (unsigned int i = 0; i < Utilities::MPI::n_mpi_processes(mpi_comm);
           ++i)
        filenames.push_back("solution-" + Utilities::int_to_string(cycle, 2) +
                            "." + Utilities::int_to_string(i, 4) + ".vtu");
      std::ofstream master_output(
        ("solution-" + Utilities::int_to_string(cycle, 2) + ".pvtu").c_str());
      data_out.write_pvtu_record(master_output, filenames);
    }
}

void
Test_Solver_Output::refine_grid()
{
  TimerOutput::Scope t(timer, "refine");

  Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
  KellyErrorEstimator<2>::estimate(
    dof_handler,
    QGauss<2 - 1>(3),
    std::map<types::boundary_id, const Function<2> *>(),
    locally_relevant_solution,
    estimated_error_per_cell);
  parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number(
    triangulation, estimated_error_per_cell, 0.3, 0.03);
  triangulation.execute_coarsening_and_refinement();
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());
  mpi_initlog();
  deallog.depth_console(0);

  Test_Solver_Output problem;
  problem.run();

  // Trilinos dumps the output into std::cout
  // We catch this output and it is written to the stdout logfile
  // Since we're interested in this output we read it back in and
  // write parts of it to the logstream
  std::ifstream inputfile;
  inputfile.open("stdout");
  Assert(inputfile.good() && inputfile.is_open(), ExcIO());
  std::string       line;
  const std::string key = "*****";
  while (std::getline(inputfile, line))
    {
      if (line.find(key) != std::string::npos)
        deallog << line << std::endl;
    }
  inputfile.close();

  deallog << "OK" << std::endl;

  return 0;
}
