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



// solves a 2D Poisson equation for linear FE_DGP elements (SIP
// discretization) with AMG preconditioner

#include "../tests.h"
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/base/function.h>
#include <deal.II/grid/tria.h>

#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/assembler.h>
#include <deal.II/meshworker/loop.h>
#include <deal.II/integrators/laplace.h>


#include <fstream>
#include <iomanip>


template <int dim>
class MatrixIntegrator : public MeshWorker::LocalIntegrator<dim>
{
public:
  void cell(MeshWorker::DoFInfo<dim> &dinfo,
            typename MeshWorker::IntegrationInfo<dim> &info) const;
  void boundary(MeshWorker::DoFInfo<dim> &dinfo,
                typename MeshWorker::IntegrationInfo<dim> &info) const;
  void face(MeshWorker::DoFInfo<dim> &dinfo1,
            MeshWorker::DoFInfo<dim> &dinfo2,
            typename MeshWorker::IntegrationInfo<dim> &info1,
            typename MeshWorker::IntegrationInfo<dim> &info2) const;
};

template <int dim>
void MatrixIntegrator<dim>
::cell(MeshWorker::DoFInfo<dim> &dinfo,
       typename MeshWorker::IntegrationInfo<dim> &info) const
{
  LocalIntegrators::Laplace::cell_matrix(dinfo.matrix(0,false).matrix,
                                         info.fe_values());
}

template <int dim>
void MatrixIntegrator<dim>
::boundary(MeshWorker::DoFInfo<dim> &dinfo,
           typename MeshWorker::IntegrationInfo<dim> &info) const
{
  const unsigned int deg = info.fe_values(0).get_fe().degree;
  LocalIntegrators::Laplace
    ::nitsche_matrix(dinfo.matrix(0,false).matrix, info.fe_values(0),
                     LocalIntegrators::Laplace::
                     compute_penalty(dinfo, dinfo, deg, deg));
}

template <int dim>
void MatrixIntegrator<dim>
::face(MeshWorker::DoFInfo<dim> &dinfo1,
       MeshWorker::DoFInfo<dim> &dinfo2,
       typename MeshWorker::IntegrationInfo<dim> &info1,
       typename MeshWorker::IntegrationInfo<dim> &info2) const
{
  const unsigned int deg = info1.fe_values(0).get_fe().degree;
  LocalIntegrators::Laplace
    ::ip_matrix(dinfo1.matrix(0,false).matrix, dinfo1.matrix(0,true).matrix,
                dinfo2.matrix(0,true).matrix, dinfo2.matrix(0,false).matrix,
                info1.fe_values(0), info2.fe_values(0),
                LocalIntegrators::Laplace::compute_penalty(dinfo1, dinfo2, deg, deg));
}


template <int dim>
class Step4
{
public:
  Step4 ();
  void run ();

private:
  void make_grid ();
  void setup_system();
  void solve ();

  Triangulation<dim>   triangulation;
  FE_DGP<dim>          fe;
  DoFHandler<dim>      dof_handler;

  TrilinosWrappers::SparseMatrix system_matrix;

  Vector<double>       solution;
  Vector<double>       system_rhs;
};




template <int dim>
Step4<dim>::Step4 ()
  :
  fe (1),
  dof_handler (triangulation)
{}


template <int dim>
void Step4<dim>::make_grid ()
{
  GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global (6);
}



template <int dim>
void Step4<dim>::setup_system ()
{
  dof_handler.distribute_dofs (fe);

  CompressedSparsityPattern c_sparsity(dof_handler.n_dofs());
  DoFTools::make_flux_sparsity_pattern (dof_handler, c_sparsity);
  system_matrix.reinit (c_sparsity);

  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());

  MappingQ1<dim> mapping;
  MeshWorker::IntegrationInfoBox<dim> info_box;
  UpdateFlags update_flags = update_values | update_gradients;
  info_box.add_update_flags_all(update_flags);
  info_box.initialize(fe, mapping);

  MeshWorker::DoFInfo<dim> dof_info(dof_handler);
  MeshWorker::Assembler::MatrixSimple<TrilinosWrappers::SparseMatrix> assembler;
  assembler.initialize(system_matrix);
  MatrixIntegrator<dim> integrator;
  MeshWorker::integration_loop<dim, dim>(dof_handler.begin_active(),
                                         dof_handler.end(),
                                         dof_info, info_box,
                                         integrator, assembler);
      
  system_matrix.compress(VectorOperation::add);

  for (unsigned int i=0; i<system_rhs.size(); ++i)
    system_rhs(i) = 0.01*i-0.000001*i*i;
}



template <int dim>
void Step4<dim>::solve ()
{

  deallog.push(Utilities::int_to_string(dof_handler.n_dofs(),5));
  TrilinosWrappers::PreconditionAMG preconditioner;
  TrilinosWrappers::PreconditionAMG::AdditionalData data;
  DoFTools::extract_constant_modes(dof_handler, std::vector<bool>(1,true),
                                   data.constant_modes);
  data.smoother_sweeps = 2;
  {
    solution = 0;
    SolverControl           solver_control (1000, 1e-10);
    SolverCG<>              solver (solver_control);
    preconditioner.initialize(system_matrix, data);
    solver.solve (system_matrix, solution, system_rhs,
                  preconditioner);
  }
  deallog.pop();
}



template <int dim>
void Step4<dim>::run()
{
  for (unsigned int cycle = 0; cycle < 2; ++cycle)
    {
      if (cycle == 0)
        make_grid();
      else
        triangulation.refine_global(1);

      setup_system();
      solve();
    }
}


int main (int argc, char **argv)
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);

  try
    {
      Step4<2> test;
      test.run();
    }
  catch (std::exception &exc)
    {
      deallog << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;

      return 1;
    }
  catch (...)
    {
      deallog << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}
