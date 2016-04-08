// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2015 by the deal.II authors
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



// take hp/random and generate the graph of constraints


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/base/quadrature_lib.h>

#include <fstream>


template <int dim>
void test ()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global (1);
  tria.begin_active()->set_refine_flag ();
  tria.execute_coarsening_and_refinement ();
  tria.refine_global (1);

  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(FE_Q<dim> (1));
  fe_collection.push_back(FE_Q<dim> (2));
  fe_collection.push_back(FE_Q<dim> (QIterated<1>(QTrapez<1>(),3)));
  fe_collection.push_back(FE_Q<dim> (QIterated<1>(QTrapez<1>(),4)));

  hp::DoFHandler<dim> dof_handler(tria);

  for (typename hp::DoFHandler<dim>::active_cell_iterator
       cell = dof_handler.begin_active();
       cell != dof_handler.end(); ++cell)
    cell->set_active_fe_index (Testing::rand() % fe_collection.size());

  dof_handler.distribute_dofs(fe_collection);

  ConstraintMatrix cm;
  DoFTools::make_hanging_node_constraints (dof_handler, cm);
  cm.write_dot (deallog.get_file_stream());
}


int main ()
{
  initlog();
  deallogfile << std::setprecision(8);
  deallog.threshold_double(1.e-10);

  test<3> ();

  deallog << "OK" << std::endl;
}
