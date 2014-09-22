// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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



// a hp-ified version of step-11


#include "../tests.h"

#include <deal.II/base/logstream.h>
#include <fstream>
std::ofstream logfile("output");

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/table_handler.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include <deal.II/lac/compressed_sparsity_pattern.h>

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cmath>


template <int dim>
class LaplaceProblem
{
public:
  LaplaceProblem (const unsigned int mapping_degree);
  void run ();

private:
  void setup_system ();
  void assemble_and_solve ();
  void solve ();

  Triangulation<dim>   triangulation;
  hp::FECollection<dim>            fe;
  hp::DoFHandler<dim>      dof_handler;
  hp::MappingCollection<dim>        mapping;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;
  ConstraintMatrix     mean_value_constraints;

  Vector<double>       solution;
  Vector<double>       system_rhs;

  TableHandler         output_table;
};



template <int dim>
LaplaceProblem<dim>::LaplaceProblem (const unsigned int mapping_degree) :
  fe (FE_Q<dim>(1)),
  dof_handler (triangulation),
  mapping (MappingQ<dim>(mapping_degree))
{
  deallog << "Using mapping with degree " << mapping_degree << ":"
          << std::endl
          << "============================"
          << std::endl;
}



template <int dim>
void LaplaceProblem<dim>::setup_system ()
{
  dof_handler.distribute_dofs (fe);
  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());

  std::vector<bool> boundary_dofs (dof_handler.n_dofs(), false);
  DoFTools::extract_boundary_dofs (dof_handler, std::vector<bool>(1,true),
                                   boundary_dofs);

  const unsigned int first_boundary_dof
    = std::distance (boundary_dofs.begin(),
                     std::find (boundary_dofs.begin(),
                                boundary_dofs.end(),
                                true));

  mean_value_constraints.clear ();
  mean_value_constraints.add_line (first_boundary_dof);
  for (unsigned int i=first_boundary_dof+1; i<dof_handler.n_dofs(); ++i)
    if (boundary_dofs[i] == true)
      mean_value_constraints.add_entry (first_boundary_dof,
                                        i, -1);
  mean_value_constraints.close ();

  CompressedSparsityPattern csp (dof_handler.n_dofs(),
                                 dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, csp);
  mean_value_constraints.condense (csp);

  sparsity_pattern.copy_from (csp);
  system_matrix.reinit (sparsity_pattern);
}



template <int dim>
void LaplaceProblem<dim>::assemble_and_solve ()
{

  const unsigned int gauss_degree
    = std::max (static_cast<unsigned int>(std::ceil(1.*(static_cast<const MappingQ<dim>&>(mapping[0]).get_degree()+1)/2)),
                2U);
  MatrixTools::create_laplace_matrix (mapping, dof_handler,
                                      hp::QCollection<dim>(QGauss<dim>(gauss_degree)),
                                      system_matrix);
  VectorTools::create_right_hand_side (mapping, dof_handler,
                                       hp::QCollection<dim>(QGauss<dim>(gauss_degree)),
                                       ConstantFunction<dim>(-2),
                                       system_rhs);
  Vector<double> tmp (system_rhs.size());
  VectorTools::create_boundary_right_hand_side (mapping, dof_handler,
                                                hp::QCollection<dim-1>(QGauss<dim-1>(gauss_degree)),
                                                ConstantFunction<dim>(1),
                                                tmp);
  system_rhs += tmp;

  mean_value_constraints.condense (system_matrix);
  mean_value_constraints.condense (system_rhs);

  solve ();
  mean_value_constraints.distribute (solution);

  Vector<float> norm_per_cell (triangulation.n_active_cells());
  VectorTools::integrate_difference (mapping, dof_handler,
                                     solution,
                                     ZeroFunction<dim>(),
                                     norm_per_cell,
                                     hp::QCollection<dim>(QGauss<dim>(gauss_degree+1)),
                                     VectorTools::H1_seminorm);
  const double norm = norm_per_cell.l2_norm();

  output_table.add_value ("cells", triangulation.n_active_cells());
  output_table.add_value ("|u|_1", norm);
  output_table.add_value ("error", std::fabs(norm-std::sqrt(3.14159265358/2)));
}



template <int dim>
void LaplaceProblem<dim>::solve ()
{
  SolverControl           solver_control (1000, 1e-12);
  SolverCG<>              cg (solver_control);

  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(system_matrix, 1.2);

  cg.solve (system_matrix, solution, system_rhs,
            preconditioner);
}



template <int dim>
void LaplaceProblem<dim>::run ()
{
  GridGenerator::hyper_ball (triangulation);
  static const HyperBallBoundary<dim> boundary;
  triangulation.set_boundary (0, boundary);

  for (unsigned int cycle=0; cycle<6; ++cycle, triangulation.refine_global(1))
    {
      setup_system ();
      assemble_and_solve ();
    };

  output_table.set_precision("|u|_1", 6);
  output_table.set_precision("error", 6);
  output_table.write_text (deallog.get_file_stream());
  deallog << std::endl;
}



int main ()
{
  try
    {
      logfile.precision(2);
      deallog << std::setprecision(2);

      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      for (unsigned int mapping_degree=1; mapping_degree<=3; ++mapping_degree)
        LaplaceProblem<2>(mapping_degree).run ();
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

  return 0;
}
