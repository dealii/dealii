// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
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



/* Author: Wolfgang Bangerth, University of Heidelberg, 1999 */
/* adapted from step-4. */


#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_out_rotation.h>
#include <deal.II/numerics/data_out_faces.h>
#include <fstream>

#include <deal.II/base/logstream.h>


std::ofstream logfile("output");


template <int dim>
class LaplaceProblem
{
public:
  LaplaceProblem ();
  void run ();

private:
  void make_grid_and_dofs ();
  void solve ();
  void output_results () const;

  Triangulation<dim>   triangulation;
  FE_Q<dim>            fe;
  DoFHandler<dim>      dof_handler;

  Vector<double>       solution;
};


template <int dim>
LaplaceProblem<dim>::LaplaceProblem () :
  fe (1), dof_handler (triangulation)
{}



template <int dim>
void LaplaceProblem<dim>::make_grid_and_dofs ()
{
  GridGenerator::hyper_cube (triangulation, 0, 1);
  triangulation.refine_global (1);
  for (unsigned int i=0; i<2; ++i)
    {
      triangulation.begin_active()->set_refine_flag ();
      triangulation.execute_coarsening_and_refinement ();
    };


  deallog << "   Number of active cells: "
          << triangulation.n_active_cells()
          << std::endl
          << "   Total number of cells: "
          << triangulation.n_cells()
          << std::endl;

  dof_handler.distribute_dofs (fe);

  deallog << "   Number of degrees of freedom: "
          << dof_handler.n_dofs()
          << std::endl;

  solution.reinit (dof_handler.n_dofs());
}




template <int dim>
void LaplaceProblem<dim>::solve ()
{
  // dummy solve. just insert some
  // arbitrary values
  for (unsigned int i=0; i<solution.size(); ++i)
    solution(i) = i;
}



template <>
void LaplaceProblem<2>::output_results () const
{
  const unsigned int dim = 2;

  // test regular output in 2d
  if (true)
    {
      DataOut<dim> data_out;
      data_out.attach_dof_handler (dof_handler);
      data_out.add_data_vector (solution, "solution");
      data_out.build_patches ();
      data_out.write_dx (logfile);
      data_out.write_gmv (logfile);
      data_out.write_gnuplot (logfile);
      data_out.set_flags (DataOutBase::UcdFlags(true));
      data_out.write_ucd (logfile);
      data_out.write_povray (logfile);
      data_out.write_eps (logfile);

      const unsigned int number_of_time_steps = 3;
      std::vector<std::vector<std::string > > piece_names(number_of_time_steps);
      piece_names[0].push_back("subdomain-01.time_step_0.vtk");
      piece_names[0].push_back("subdomain-02.time_step_0.vtk");
      piece_names[1].push_back("subdomain-01.time_step_1.vtk");
      piece_names[1].push_back("subdomain-02.time_step_1.vtk");
      piece_names[2].push_back("subdomain-01.time_step_2.vtk");
      piece_names[2].push_back("subdomain-02.time_step_2.vtk");
      data_out.write_visit_record(logfile, piece_names);
      data_out.write_visit_record(logfile, piece_names[01]);
    };

  // test DataOutRotation in 2d
  if (true)
    {
      DataOutRotation<dim> data_out;
      data_out.attach_dof_handler (dof_handler);
      data_out.add_data_vector (solution, "solution");
      data_out.build_patches (3);
      data_out.write_dx (logfile);
      data_out.write_gmv (logfile);
      data_out.write_gnuplot (logfile);
      data_out.set_flags (DataOutBase::UcdFlags(true));
      data_out.write_ucd (logfile);
    };
}



template <>
void LaplaceProblem<3>::output_results () const
{
  const unsigned int dim = 3;

  // test regular output in 3d
  if (true)
    {
      DataOut<dim> data_out;
      data_out.attach_dof_handler (dof_handler);
      data_out.add_data_vector (solution, "solution");
      data_out.build_patches ();
      data_out.write_dx (logfile);
      data_out.write_gmv (logfile);
      data_out.write_gnuplot (logfile);
      data_out.set_flags (DataOutBase::UcdFlags(true));
      data_out.write_ucd (logfile);
    };

  // test DataOutFaces in 3d. note:
  // not all output formats support
  // this
  if (true)
    {
      DataOutFaces<dim> data_out;
      data_out.attach_dof_handler (dof_handler);
      data_out.add_data_vector (solution, "solution");
      data_out.build_patches (3);
      data_out.write_dx (logfile);
      data_out.write_gmv (logfile);
      data_out.write_gnuplot (logfile);
      data_out.set_flags (DataOutBase::UcdFlags(true));
      data_out.write_ucd (logfile);
    };
}



template <int dim>
void LaplaceProblem<dim>::run ()
{
  make_grid_and_dofs();
  solve ();
  output_results ();
}



int main ()
{
  deallog.depth_console (0);
  logfile << std::setprecision(2);
  deallog << std::setprecision(2);

  LaplaceProblem<2> laplace_problem_2d;
  laplace_problem_2d.run ();

  LaplaceProblem<3> laplace_problem_3d;
  laplace_problem_3d.run ();

  return 0;
}
