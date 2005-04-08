//----------------------------  subdomain_ids_02.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2001, 2002, 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  subdomain_ids_02.cc  ---------------------------

// check DoFTools::get_subdomain_association


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


std::ofstream logfile("subdomain_ids_02.output");


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
    {
                                       // count number on dofs on
                                       // subdomain. this time it should add
                                       // up, since each dof is uniquely
                                       // associated
      deallog << std::count (subdomain_association.begin(),
                             subdomain_association.end(), subdomain)
              << std::endl;
    }

                                   // make sure that all subdomain_ids are
                                   // really in the allowed range. if this is
                                   // the case, then we have also proven that
                                   // the numbers really add up correctly,
                                   // since every dof is assigned a valid
                                   // subdomain id
  for (unsigned int i=0;i<dof_handler.n_dofs(); ++i)
    Assert (subdomain_association[i] < (1<<dim),
            ExcInternalError());
}


int main ()
{
  logfile.precision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1> ();
  test<2> ();
  test<3> ();
  
  return 0;
}
