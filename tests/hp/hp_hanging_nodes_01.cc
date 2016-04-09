// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2015 by the deal.II authors
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


#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>


// This a test for the hp capable version of the make_hanging_node_constraints
// method. It uses a triangulation with one refined element beside an
// unrefined element to create the constraints for this configuration.

template <int dim>
int generate_grid (Triangulation<dim> &tria)
{
  Point<dim> p1,
        p2;
  std::vector<unsigned int> sub_div;

  // Define a rectangular shape
  for (unsigned int d=0; d < dim; ++d)
    {
      p1(d) = 0;
      p2(d) = (d == 0) ? 2.0 : 1.0;
      sub_div.push_back ( (d == 0) ? 2 : 1);
    }
  GridGenerator::subdivided_hyper_rectangle (tria, sub_div, p1, p2, true);

  // Refine the first cell.
  tria.begin_active ()->set_refine_flag ();
  tria.execute_coarsening_and_refinement ();

  return (0);
}



template <int dim>
void test_constraints (hp::FECollection<dim> &fe_coll)
{
  Triangulation<dim> tria;

  // Setup a rectangular domain
  // where one cell is h-refined,
  // while the other cell is
  // unrefined. Furthermore every cell
  // gets a different active_fe_index.
  // This should serve as a testcase
  // for the hanging node constraints.
  generate_grid (tria);

  // Now assign increasing
  // active_fe_indices to
  // the different cells.
  hp::DoFHandler<dim> dof_handler (tria);
  typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active (),
                                                     endc = dof_handler.end ();
  unsigned int fe_indx = 0;
  for (; cell != endc; ++cell)
    {
      cell->set_active_fe_index (fe_indx);
      ++fe_indx;
    }

  // Distribute DoFs;
  dof_handler.distribute_dofs (fe_coll);
  deallog << "DoFs: " << dof_handler.n_dofs () << std::endl;

  // Create the constraints.
  ConstraintMatrix constraint_matrix;

  DoFTools::make_hanging_node_constraints (dof_handler,
                                           constraint_matrix);

  // Output the constraints
  constraint_matrix.print (deallog.get_file_stream ());
}


template <int dim>
void test_constraints_old (FiniteElement<dim> &fe)
{
  Triangulation<dim> tria;

  // Setup a rectangular domain
  // where one cell is h-refined,
  // while the other cell is
  // unrefined. Furthermore every cell
  // gets a different active_fe_index.
  // This should serve as a testcase
  // for the hanging node constraints.
  generate_grid (tria);

  // Now assign increasing
  // active_fe_indices to
  // the different cells.
  DoFHandler<dim> dof_handler (tria);

  // Distribute DoFs;
  dof_handler.distribute_dofs (fe);
  deallog << "DoFs: " << dof_handler.n_dofs () << std::endl;

  // Create the constraints.
  ConstraintMatrix constraint_matrix;

  DoFTools::make_hanging_node_constraints (dof_handler,
                                           constraint_matrix);

  // Output the constraints
  constraint_matrix.print (deallog.get_file_stream ());
}

int main ()
{
  std::ofstream logfile("output");
  logfile.precision(2);

  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  FE_Q<2> fe_1 (1);
  FE_Q<2> fe_2 (2);
  FE_Q<2> fe_3 (QIterated<1>(QTrapez<1>(),3));

  hp::FECollection<2> fe_coll2;
  fe_coll2.push_back (fe_3);
  fe_coll2.push_back (fe_2);
  fe_coll2.push_back (fe_2);

  fe_coll2.push_back (fe_2);
  fe_coll2.push_back (fe_3);

  test_constraints<2> (fe_coll2);

  test_constraints_old<2> (fe_1);

  return 0;
}
