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
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/dofs/dof_tools.h>

// check
//   DoFTools::
//   count_dofs_per_component (...);
// for the hp case


using namespace std;

template <int dim>
void test ()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(3);

  //define DoFhandler and FEs
  FE_Q<dim> u(2);
  FE_Q<dim> p(1);

  FESystem<dim> fe_system(u, 2, p, 1);

  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(fe_system);

  hp::DoFHandler<dim> hp_dof_handler(triangulation);
  DoFHandler<dim> dof_handler(triangulation);

  //distribute dofs
  hp_dof_handler.distribute_dofs(fe_collection);
  dof_handler.distribute_dofs(fe_system);

  //count dofs per component and show them on the screen
  std::vector<types::global_dof_index> dofs_per_component(3,0);
  std::vector<types::global_dof_index> dofs_per_component_hp(3,0);
  DoFTools::count_dofs_per_component(dof_handler, dofs_per_component);
  DoFTools::count_dofs_per_component(hp_dof_handler, dofs_per_component_hp);

  for (unsigned int i=0; i<3; i++)
    {
      deallog <<"DoFs in the " <<i<<". component for classical FE: "<< dofs_per_component.at(i)<< std::endl;
      deallog <<"DoFs in the " <<i<<". component for hp FE: "<< dofs_per_component_hp.at(i)<< std::endl;

      Assert (dofs_per_component.at(i) == dofs_per_component_hp.at(i),
              ExcInternalError());
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
