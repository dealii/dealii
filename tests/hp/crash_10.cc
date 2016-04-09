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



// a version of hp_hanging_node_02 that crashed at the time
// of writing the time

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_tools.h>

#include <fstream>

std::ofstream logfile("output");


template <int dim>
void test ()
{
  Triangulation<dim>     triangulation;
  hp::FECollection<dim>              fe;
  hp::DoFHandler<dim>        dof_handler(triangulation);
  ConstraintMatrix     hanging_node_constraints;

  FE_Q<dim> fe_1 (1), fe_2 (2), fe_3 (QIterated<1>(QTrapez<1>(),3)), fe_4 (QIterated<1>(QTrapez<1>(),4));

  fe.push_back (fe_1);
  fe.push_back (fe_2);
  fe.push_back (fe_3);
  fe.push_back (fe_4);



  GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global (2);
  deallog << "Number of active cells: "
          << triangulation.n_active_cells()
          << std::endl;
  deallog << "Total number of cells: "
          << triangulation.n_cells()
          << std::endl;

  // Now to the p-Method. Assign
  // random active_fe_indices to the
  // different cells.
  typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active (),
                                                     endc = dof_handler.end ();
  for (; cell != endc; ++cell)
    cell->set_active_fe_index (Testing::rand() % fe.size());

  dof_handler.distribute_dofs (fe);

  cell = dof_handler.begin_active ();
  for (; cell != endc; ++cell)
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
      {
        deallog << "cell=" << cell << ", face=" << f
                << ", fe_index=" << cell->active_fe_index()
                << ", dofs=";
        std::vector<types::global_dof_index> dofs (fe[cell->active_fe_index()].dofs_per_face);
        cell->face(f)->get_dof_indices (dofs,
                                        cell->active_fe_index());
        for (unsigned int i=0; i<dofs.size(); ++i)
          deallog << dofs[i] << ' ';
        deallog << std::endl;
      }

  DoFTools::make_hanging_node_constraints (dof_handler,
                                           hanging_node_constraints);

  hanging_node_constraints.print (deallog.get_file_stream ());

  hanging_node_constraints.close ();
}


int main ()
{
  logfile.precision(2);

  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  // this test depends on the right
  // starting value of the random
  // number generator. set it to the
  // same value it has in
  // hp_hanging_nodes_02
  for (unsigned int i=0; i<64; ++i)
    Testing::rand();
  test<3> ();

  return 0;
}
