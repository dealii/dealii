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
#include <grid/tria_iterator.h>
#include <hp/fe_collection.h>
#include <hp/dof_handler.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <fe/fe_q.h>
#include <fe/fe_system.h>
#include <dofs/dof_tools.h>

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
  std::vector<unsigned int> dofs_per_component_hp(3,0);
  DoFTools::count_dofs_per_component(hp_dof_handler, dofs_per_component_hp);

  for (unsigned int i=0; i<3; i++)
    {
      deallog <<"DoFs in the " <<i<<". component for hp FE: "<< dofs_per_component_hp.at(i)<< std::endl;
    }

}

int main () 
{
  std::ofstream logfile("count_dofs_per_component_hp_02/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1> ();
  test<2> ();
  test<3> ();
  return 0;
}
