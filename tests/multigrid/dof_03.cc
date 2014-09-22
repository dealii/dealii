// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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


// check DoFHandler::has_level_dofs and DoFHandler::has_active_dofs

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/multigrid/mg_dof_handler.h>

#include <fstream>
#include <iomanip>
#include <iomanip>
#include <algorithm>



template <int dim>
void check()
{
  FE_Q<dim> fe(1);

  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(1);

  // check the two functions mentioned above
  // in their natural order
  {
    DoFHandler<dim> dof(tr);
    deallog << "check " << dim << " before distribute "
            << dof.has_active_dofs() << ' '
            << dof.has_level_dofs()
            << std::endl;

    dof.distribute_dofs(fe);
    deallog << "check " << dim << " after  distribute "
            << dof.has_active_dofs() << ' '
            << dof.has_level_dofs()
            << std::endl;


    dof.distribute_mg_dofs(fe);
    deallog << "check " << dim << " level  distribute "
            << dof.has_active_dofs() << ' '
            << dof.has_level_dofs()
            << std::endl;

    dof.clear();
    deallog << "check " << dim << " after  clear      "
            << dof.has_active_dofs() << ' '
            << dof.has_level_dofs()
            << std::endl;
  }

  // now check them the other way around
  {
    DoFHandler<dim> dof(tr);
    deallog << "check " << dim << " before distribute "
            << dof.has_active_dofs() << ' '
            << dof.has_level_dofs()
            << std::endl;

    dof.distribute_mg_dofs(fe);
    deallog << "check " << dim << " level  distribute "
            << dof.has_active_dofs() << ' '
            << dof.has_level_dofs()
            << std::endl;

    dof.distribute_dofs(fe);
    deallog << "check " << dim << " after  distribute "
            << dof.has_active_dofs() << ' '
            << dof.has_level_dofs()
            << std::endl;

    dof.clear();
    deallog << "check " << dim << " after  clear      "
            << dof.has_active_dofs() << ' '
            << dof.has_level_dofs()
            << std::endl;
  }
}

int main()
{
  initlog(__FILE__);
  check<1> ();
  check<2> ();
  check<3> ();
}
