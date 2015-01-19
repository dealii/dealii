// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2014 by the deal.II authors
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


// check DoFTools::count_dofs_with_subdomain_association


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/dofs/dof_tools.h>

#include <fstream>
#include <algorithm>
#include <cmath>


std::ofstream logfile("output");


template <int dim>
void test ()
{
  deallog << dim << 'D' << std::endl;
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, -1, 1);
  tria.refine_global (2);

  // we now have a number of cells,
  // flag them with some subdomain
  // ids based on their position, in
  // particular we take the quadrant
  // (octant)
  typename Triangulation<dim>::active_cell_iterator
  cell = tria.begin_active (),
  endc = tria.end ();
  for (; cell!=endc; ++cell)
    {
      unsigned int subdomain = 0;
      for (unsigned int d=0; d<dim; ++d)
        if (cell->center()(d) > 0)
          subdomain |= (1<<d);
      Assert (subdomain < (1<<dim), ExcInternalError());

      cell->set_subdomain_id (subdomain);
    };

  // distribute some degrees of freedom and
  // output some information on them
  FESystem<dim> fe(FE_Q<dim>(2),dim, FE_DGQ<dim>(1), 1);
  DoFHandler<dim> dof_handler (tria);
  dof_handler.distribute_dofs (fe);
  deallog << dof_handler.n_dofs() << std::endl;

  std::vector<unsigned int> subdomain_association (dof_handler.n_dofs());
  DoFTools::get_subdomain_association (dof_handler,
                                       subdomain_association);
  for (unsigned int subdomain=0; subdomain<(1<<dim); ++subdomain)
    // check that the number of dofs
    // associated is also what the respective
    // function returns
    Assert (static_cast<unsigned int>
            (std::count (subdomain_association.begin(),
                         subdomain_association.end(), subdomain))
            ==
            DoFTools::count_dofs_with_subdomain_association (dof_handler,
                                                             subdomain),
            ExcInternalError());

  deallog << "OK" << std::endl;
}


int main ()
{
  deallog << std::setprecision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1> ();
  test<2> ();
  test<3> ();

  return 0;
}
