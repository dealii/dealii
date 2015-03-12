// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2014 by the deal.II authors
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
  std::ofstream logfile ("output");
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
