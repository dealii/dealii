//----------------------------  subdomain_ids.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  subdomain_ids.cc  ---------------------------

// check some things about subdomain_ids


#include <base/logstream.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <fe/fe_dgq.h>
#include <fe/fe_q.h>
#include <dofs/dof_tools.h>

#include <fstream>
#include <algorithm>
#include <cmath>


std::ofstream logfile("subdomain_ids.output");


DeclException2 (ExcNumberMismatch,
		int, int,
		<< "The numbers " << arg1 << " and " << arg2
		<< " should be equal, but are not.");


template <int dim>
void test ()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, -1, 1);
  tria.refine_global (3);

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

				   // refine twice (for 3d only once)
				   // and count again the numbers of
				   // cells in each subdomain
  tria.refine_global ((dim != 3) ? 2 : 1);
  if (true)
    {
      cell = tria.begin_active();
      endc = tria.end();
      std::vector<unsigned int> subdomain_cells(1<<dim, 0);
      for (; cell!=endc; ++cell)
	{
	  Assert (cell->subdomain_id() < (1<<dim), ExcInternalError());
	  ++subdomain_cells[cell->subdomain_id()];
	};
      for (unsigned int i=0; i<(1<<dim); ++i)
	Assert (subdomain_cells[i] == tria.n_active_cells()/(1<<dim),
		ExcNumberMismatch(subdomain_cells[i],
				  tria.n_active_cells()/(1<<dim)));
      deallog << "Check 1 (dim=" << dim << ") ok" << std::endl;
    };
  
				   // coarsen once and check again
  if (true)
    {
      cell = tria.begin_active();
      endc = tria.end();
      for (; cell!=endc; ++cell)
	cell->set_coarsen_flag ();
      tria.execute_coarsening_and_refinement ();
      
      cell = tria.begin_active();
      endc = tria.end();
      std::vector<unsigned int> subdomain_cells(1<<dim, 0);
      for (; cell!=endc; ++cell)
	{
	  Assert (cell->subdomain_id() < (1<<dim), ExcInternalError());
	  ++subdomain_cells[cell->subdomain_id()];
	};
      for (unsigned int i=0; i<(1<<dim); ++i)
	Assert (subdomain_cells[i] == tria.n_active_cells()/(1<<dim),
		ExcNumberMismatch(subdomain_cells[i],
				  tria.n_active_cells()/(1<<dim)));
      deallog << "Check 2 (dim=" << dim << ") ok" << std::endl;
    };

				   // check 3: assign DQ2 elements to
				   // cells and count degrees of
				   // freedom associated with each
				   // subdomain
  if (true)
    {
      FE_DGQ<dim> fe(2);
      DoFHandler<dim> dof_handler (tria);
      dof_handler.distribute_dofs (fe);
      std::vector<bool> selected_dofs (dof_handler.n_dofs());
      for (unsigned int subdomain=0; subdomain<(1<<dim); ++subdomain)
	{
	  DoFTools::extract_subdomain_dofs (dof_handler, subdomain,
					    selected_dofs);
	  Assert (static_cast<unsigned int>(std::count (selected_dofs.begin(),
							selected_dofs.end(),
							true))
		  ==
		  dof_handler.n_dofs() / (1<<dim),
		  ExcNumberMismatch(std::count (selected_dofs.begin(),
						selected_dofs.end(),
						true),
				    dof_handler.n_dofs() / (1<<dim)));
	}
      deallog << "Check 3 (dim=" << dim << ") ok" << std::endl;      
    };


				   // check 4: check again for
				   // continuous elements. note that
				   // the number of dofs here is
				   // different, since dofs can be on
				   // several subdomain at once
  if (true)
    {
      FE_Q<dim> fe(1);
      DoFHandler<dim> dof_handler (tria);
      dof_handler.distribute_dofs (fe);

      const unsigned int cells_per_direction
	= static_cast<unsigned int>(rint(pow(tria.n_active_cells(), 1./dim)));
      
      std::vector<bool> selected_dofs (dof_handler.n_dofs());
      for (unsigned int subdomain=0; subdomain<(1<<dim); ++subdomain)
	{
	  DoFTools::extract_subdomain_dofs (dof_handler, subdomain,
					    selected_dofs);
	  Assert (static_cast<unsigned int>(std::count (selected_dofs.begin(),
							selected_dofs.end(),
							true))
		  ==
		  pow(cells_per_direction/2+1,dim),
		  ExcNumberMismatch(std::count (selected_dofs.begin(),
						selected_dofs.end(),
						true),
				    (int)pow(cells_per_direction/2+1,dim)));
	}
      deallog << "Check 4 (dim=" << dim << ") ok" << std::endl;      
    };
};


int main ()
{
  logfile.precision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<1> ();
  test<2> ();
  test<3> ();
  
  return 0;
};
