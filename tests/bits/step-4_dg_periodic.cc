// ---------------------------------------------------------------------
//
// Copyright (C) 2015 by the deal.II authors
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



// solves a 2D Poisson equation with FE_DGQ elements and periodic boundary conditions

#include "../tests.h"
#include <deal.II/base/function.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/integrators/laplace.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/assembler.h>
#include <deal.II/meshworker/loop.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>


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
  void output_results (const unsigned int cycle) const;
  void check_periodicity (const unsigned int cycle) const;


  Triangulation<dim>   triangulation;
  FE_DGQ<dim>          fe;
  DoFHandler<dim>      dof_handler;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

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
  GridGenerator::hyper_cube (triangulation, -1, 1, true);
  typedef typename dealii::Triangulation<dim>::cell_iterator CellIteratorTria;
  std::vector<dealii::GridTools::PeriodicFacePair<CellIteratorTria> > periodic_faces;
  const unsigned int b_id1 = 2;
  const unsigned int b_id2 = 3;
  const unsigned int direction = 1;

  dealii::GridTools::collect_periodic_faces (triangulation, b_id1, b_id2,
                                             direction, periodic_faces, dealii::Tensor<1,dim>());
  triangulation.add_periodicity(periodic_faces);
  triangulation.refine_global (1);
}



template <int dim>
void Step4<dim>::setup_system ()
{
  dof_handler.distribute_dofs (fe);

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_flux_sparsity_pattern (dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);
  system_matrix.reinit (sparsity_pattern);

  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());

  MappingQGeneric<dim> mapping(1);
  MeshWorker::IntegrationInfoBox<dim> info_box;
  UpdateFlags update_flags = update_values | update_gradients;
  info_box.add_update_flags_all(update_flags);
  info_box.initialize(fe, mapping);

  MeshWorker::DoFInfo<dim> dof_info(dof_handler);
  MeshWorker::Assembler::MatrixSimple<SparseMatrix<double> > assembler;
  assembler.initialize(system_matrix);
  MatrixIntegrator<dim> integrator;
  MeshWorker::integration_loop<dim, dim>(dof_handler.begin_active(),
                                         dof_handler.end(),
                                         dof_info, info_box,
                                         integrator, assembler);

  for (unsigned int i=0; i<system_rhs.size(); ++i)
    system_rhs(i) = 0.01*i-0.000001*i*i;
}



template <int dim>
void Step4<dim>::solve ()
{
  solution = 0;
  SolverControl           solver_control (1000, 1e-10, false, false);
  SolverCG<>              solver (solver_control);

  PreconditionJacobi<SparseMatrix<double> > preconditioner;
  preconditioner.initialize(system_matrix);

  solver.solve (system_matrix, solution, system_rhs,
                preconditioner);
}



template <int dim>
void Step4<dim>::output_results (const unsigned int cycle) const
{
  std::string filename = "solution-"+dealii::Utilities::int_to_string(cycle,2);

  dealii::DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "u");

  data_out.build_patches (fe.degree);
  std::ofstream output ((filename + ".vtk").c_str());
  data_out.write_vtk (output);
}



template <int dim>
void Step4<dim>::check_periodicity (const unsigned int cycle) const
{}

template <>
void Step4<2>::check_periodicity (const unsigned int cycle) const
{
  unsigned int n_points = 4;
  for (unsigned int i = 0; i<cycle; i++)
    n_points*=2;

  //don't test exactly at the support points, since point_value is not stable there
  const double eps = 1./(16.*n_points);

  bool all_passed = true;

  for (unsigned int i=1; i< n_points; i++)
    {
      Vector<double> value1(1);
      Vector<double> value2(1);

      Point<2> point1;
      point1(0)=2*(1.*i/n_points+eps)-1;
      point1(1)=-1.;
      Point<2> point2;
      point2(0)=2*(1.*i/n_points+eps)-1;
      point2(1)=1.;

      VectorTools::point_value (dof_handler, solution, point1, value1);
      VectorTools::point_value (dof_handler, solution, point2, value2);

      const double rel_error = std::abs((value2[0]-value1[0])/value1[0]);
      const double rel_tol = 1./std::pow(2., cycle);

      if (rel_error < rel_tol)
        deallog << point1 << "\t pass" << std::endl;
      else
        {
          deallog << point1 << "\t fail" << std::endl;
          deallog << point1 << "\t" << value1[0] << "\t"
                  << value2[0] << "\t"
                  << rel_error << std::endl;
          all_passed = false;
        }
    }
  AssertThrow(all_passed, ExcInternalError());
}



template <int dim>
void Step4<dim>::run()
{
  for (unsigned int cycle = 0; cycle < 4; ++cycle)
    {
      if (cycle == 0)
        make_grid();
      else
        triangulation.refine_global(1);

      setup_system();
      solve();
      //output_results(cycle);
      deallog.push(Utilities::int_to_string(dof_handler.n_dofs(),5));
      check_periodicity(cycle);
      deallog.pop();
    }
}


int main (int argc, char **argv)
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

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
