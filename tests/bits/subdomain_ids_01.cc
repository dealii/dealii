//----------------------------  subdomain_ids_01.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2001, 2002, 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  subdomain_ids_01.cc  ---------------------------

// check DoFTools::extract_subdomain_dofs


#include "../tests.h"
#include <base/logstream.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_dgq.h>
#include <fe/fe_q.h>
#include <fe/fe_system.h>
#include <dofs/dof_tools.h>

#include <fstream>
#include <algorithm>
#include <cmath>


std::ofstream logfile("subdomain_ids_01.output");


DeclException2 (ExcNumberMismatch,
		int, int,
		<< "The numbers " << arg1 << " and " << arg2
		<< " should be equal, but are not.");


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
  
  std::vector<bool> selected_dofs (dof_handler.n_dofs());
  for (unsigned int subdomain=0; subdomain<(1<<dim); ++subdomain)
    {
                                       // count number on dofs on
                                       // subdomain. note that they add up to
                                       // more than n_dofs() since the ones on
                                       // the interfaces count for each
                                       // subdomain
      DoFTools::extract_subdomain_dofs (dof_handler, subdomain,
                                        selected_dofs);
      deallog << std::count (selected_dofs.begin(),
                             selected_dofs.end(), true)
              << std::endl;
    }
}


int main ()
{
  logfile.precision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<1> ();
  test<2> ();
  test<3> ();
  
  return 0;
}
