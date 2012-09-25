//----------------------------  extract_dofs_by_component_01.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2003, 2004, 2007, 2010, 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  extract_dofs_by_component_01.cc  ---------------------------


// test internal::extract_dofs_by_component for some corner cases that
// I was unsure about when refactoring some code in there
//
// this particular test checks the call path to
// internal::extract_dofs_by_component from DoFTools::extract_dofs via
// the component_select flag


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/dofs/dof_tools.h>

#include <fstream>




template <int dim>
void
check ()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr, -1,1);
  tr.refine_global (1);

  FESystem<dim> element (FE_Q<dim>(1), 1,
			 FE_RaviartThomas<dim>(0), 1,
			 FE_Q<dim>(1), 1,
			 FE_Nedelec<dim>(0), 1);
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(element);

				   // try all possible component
				   // masks, which we encode as bit
				   // strings
  for (unsigned int int_mask=0; int_mask<(1U<<element.n_components()); ++int_mask)
  {
    std::vector<bool> component_mask (element.n_components());
    for (unsigned int c=0; c<element.n_components(); ++c)
      component_mask[c] = (int_mask & (1<<c));

    std::vector<bool> dofs (dof.n_dofs());
    DoFTools::extract_dofs (dof, ComponentMask(component_mask), dofs);

    for (unsigned int d=0; d<dof.n_dofs(); ++d)
      deallog << dofs[d];
    deallog << std::endl;
  }
}


int main ()
{
  std::ofstream logfile ("extract_dofs_by_component_01/output");
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
