// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/dofs/dof_tools.h>

// check
//   DoFTools::
//   count_dofs_per_component (...);
// for the hp case
//
// in contrast to the _01 test also check that this works if the element
// collection has more than one element


using namespace std;

template <int dim>
void test ()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(2);

  //define DoFhandler and FEs
  FE_Q<dim> u1(2);
  FE_Q<dim> p1(1);
  FESystem<dim> fe_system1(u1, 2, p1, 1);
  FE_Q<dim> u2(3);
  FE_Q<dim> p2(2);
  FESystem<dim> fe_system2(u2, 2, p2, 1);

  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(fe_system1);
  fe_collection.push_back(fe_system2);

  hp::DoFHandler<dim> hp_dof_handler(triangulation);
  hp_dof_handler.begin_active()->set_active_fe_index(1);

  //distribute dofs
  hp_dof_handler.distribute_dofs(fe_collection);

  //count dofs per component and show them on the screen
  std::vector<types::global_dof_index> dofs_per_component_hp(3,0);
  DoFTools::count_dofs_per_component(hp_dof_handler, dofs_per_component_hp);

  for (unsigned int i=0; i<3; i++)
    {
      deallog <<"DoFs in the " <<i<<". component for hp FE: "<< dofs_per_component_hp.at(i)<< std::endl;
    }

}

int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1> ();
  test<2> ();
  test<3> ();
  return 0;
}
