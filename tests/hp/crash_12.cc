// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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



// check the complex case described in the hp paper by playing through all
// sorts of arrangements of finite elements on one coarse and one refined cell
//
// this code in particular tests some compensating code in
// dof_tools.cc, where we have to make sure that we select a suitable
// set of master dofs. this is mostly trivial in 2d and for most fe
// combinations in 3d as well. the exceptions are that it doesn't work
// as easily in 3d for the combinations Q4/Q3, Q5/Q3, and
// Q5/Q4. Higher order finite elements in 3d will probably only
// exacerbate the problem, but the code there appears to be robust.

char logname[] = "output";


#include "../tests.h"
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>

#include <deal.II/grid/tria.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/fe/fe_q.h>

#include <fstream>
#include <vector>



template <int dim>
void test ()
{
  // create a mesh like this (viewed
  // from top, if in 3d):
  // *---*---*
  // | 0 | 1 |
  // *---*---*
  Triangulation<dim>     triangulation;
  std::vector<unsigned int> subdivisions (dim, 1);
  subdivisions[0] = 2;
  GridGenerator::subdivided_hyper_rectangle (triangulation, subdivisions,
                                             Point<dim>(),
                                             (dim == 3 ?
                                              Point<dim>(2,1,1) :
                                              Point<dim>(2,1)));
  (++triangulation.begin_active())->set_refine_flag ();
  triangulation.execute_coarsening_and_refinement ();

  hp::FECollection<dim> fe;
  fe.push_back (FE_Q<dim>(1));
  fe.push_back (FE_Q<dim>(2));
  fe.push_back (FE_Q<dim>(3));
  fe.push_back (FE_Q<dim>(4));
  fe.push_back (FE_Q<dim>(5));

  hp::DoFHandler<dim>        dof_handler(triangulation);

  for (unsigned int i=0; i<fe.size(); ++i)
    for (unsigned int j=0; j<fe.size(); ++j)
      {
        deallog << "Testing " << fe[i].get_name()
                << " vs. " << fe[j].get_name()
                << std::endl;

        // set fe on coarse cell to 'i', on
        // all fine cells to 'j'
        typename hp::DoFHandler<dim>::active_cell_iterator
        cell = dof_handler.begin_active();
        cell->set_active_fe_index (i);
        ++cell;

        for (; cell != dof_handler.end(); ++cell)
          cell->set_active_fe_index (j);

        dof_handler.distribute_dofs (fe);

        ConstraintMatrix constraints;
        DoFTools::make_hanging_node_constraints (dof_handler,
                                                 constraints);
        constraints.close ();

        constraints.print (deallog.get_file_stream());
      }
}



int main ()
{
  std::ofstream logfile(logname);
  logfile.precision (7);

  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);


  test<2> ();
  test<3> ();
}

