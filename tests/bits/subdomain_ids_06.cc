//----------------------------  subdomain_ids_06.cc  ---------------------------
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
//----------------------------  subdomain_ids_06.cc  ---------------------------

// check GridTools::count_cells_with_subdomain_association


#include "../tests.h"
#include <base/logstream.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <grid/grid_tools.h>
#include <fe/fe_dgq.h>
#include <fe/fe_q.h>
#include <fe/fe_system.h>
#include <dofs/dof_tools.h>

#include <fstream>
#include <algorithm>
#include <cmath>


std::ofstream logfile("subdomain_ids_06.output");


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
    }

  std::vector<unsigned int> subdomain_association (tria.n_active_cells());
  GridTools::get_subdomain_association (tria,
                                        subdomain_association);
  for (unsigned int subdomain=0; subdomain<(1<<dim); ++subdomain)
    {
                                       // check that the number of cells
                                       // associated is also what the respective
                                       // function returns
      Assert (static_cast<unsigned int>
              (std::count (subdomain_association.begin(),
                           subdomain_association.end(), subdomain))
              ==
              GridTools::count_cells_with_subdomain_association (tria,
                                                                 subdomain),
              ExcInternalError());

                                       // ...and that this is also the correct
                                       // number
      Assert (GridTools::count_cells_with_subdomain_association (tria,
                                                                 subdomain)
              == (tria.n_active_cells() / (1<<dim)),
              ExcInternalError());
    }
  
  deallog << "OK" << std::endl;
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
