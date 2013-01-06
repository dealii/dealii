//----------------------------  extract_dofs_by_component_01_mg.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2003, 2004, 2007, 2010, 2012, 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  extract_dofs_by_component_01_mg.cc  ---------------------------


// test internal::extract_dofs_by_component for some corner cases that
// I was unsure about when refactoring some code in there
//
// this is a test for the MG version of this test


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/multigrid/mg_dof_accessor.h>

#include <fstream>




template <int dim>
void
check ()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr, -1,1);
  tr.refine_global (1);
  tr.begin_active()->set_refine_flag();
  tr.execute_coarsening_and_refinement();

  FESystem<dim> element (FE_Q<dim>(2), 1,
			 FE_Nedelec<dim>(0), 1);
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(element);
  dof.distribute_mg_dofs(element);

				   // try all possible component
				   // masks, which we encode as bit
				   // strings
  for (unsigned int int_mask=0; int_mask<(1U<<element.n_components()); ++int_mask)
  {
    std::vector<bool> component_mask (element.n_components());
    for (unsigned int c=0; c<element.n_components(); ++c)
      component_mask[c] = (int_mask & (1<<c));

    for (unsigned int level=0; level<tr.n_levels(); ++level)
      {
	deallog << "level=" << level << std::endl;

	std::vector<bool> dofs (dof.n_dofs(level));
	DoFTools::extract_level_dofs (level, dof, ComponentMask(component_mask), dofs);

	for (unsigned int d=0; d<dofs.size(); ++d)
	  deallog << dofs[d];
	deallog << std::endl;
      }
  }
}


int main ()
{
  std::ofstream logfile ("extract_dofs_by_component_01_mg/output");
  deallog << std::setprecision (2);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console (0);

  deallog.push ("2d");
  check<2> ();
  deallog.pop ();
  deallog.push ("3d");
  check<3> ();
  deallog.pop ();
}
