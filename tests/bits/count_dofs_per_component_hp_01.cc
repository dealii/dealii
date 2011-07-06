//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2003, 2004, 2005, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------

#include "../tests.h"
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <hp/fe_collection.h>
#include <hp/dof_handler.h>
#include <dofs/dof_handler.h>
#include <fe/fe_q.h>
#include <fe/fe_system.h>
#include <dofs/dof_tools.h>

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
  std::vector<unsigned int> dofs_per_component(3,0);
  std::vector<unsigned int> dofs_per_component_hp(3,0);
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
  std::ofstream logfile("count_dofs_per_component_hp_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1> ();
  test<2> ();
  test<3> ();
  return 0;
}
