// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
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



// test internal::extract_dofs_by_component for some corner cases that
// I was unsure about when refactoring some code in there
//
// this particular test checks the call path to
// internal::extract_dofs_by_component from DoFTools::extract_dofs via
// the component_select flag
//
// this is a variant of the _01x test that shows that we were doing something
// wrong all along even on mainline but that nobody apparently ever realized
// this. the bug is that internal::extract_dofs_by_component got things wrong
// when a non-primitive element was part of an FESystem but it was
// accidentally correct whenever there was only one such element; this test
// verifies the same with two non-primitive elements in one FESystem


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_dgq.h>
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

  FESystem<dim> element (FE_Nedelec<dim>(0), 2);
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(element);

  // use a mask that only has the one
  // component set and cycle over
  // which one that is
  for (unsigned int comp=0; comp<element.n_components(); ++comp)
    {
      std::vector<bool> component_mask (element.n_components(), false);
      component_mask[comp] = true;

      std::vector<bool> dofs (dof.n_dofs());
      DoFTools::extract_dofs (dof, ComponentMask(component_mask), dofs);

      for (unsigned int d=0; d<dof.n_dofs(); ++d)
        deallog << dofs[d];
      deallog << std::endl;
    }
}


int main ()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (2);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console (0);

  deallog.push ("2d");
  check<2> ();
  deallog.pop ();
}
