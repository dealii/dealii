//----------------------------  fe_nothing_04.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  fe_nothing_04.cc  ---------------------------


// test that FE_Nothing works as intended


#include "../tests.h"
#include <base/logstream.h>
#include <base/quadrature_lib.h>
#include <fe/fe_nothing.h>
#include <fe/fe_q.h>
#include <fe/fe_dgq.h>
#include <fe/fe_raviart_thomas.h>
#include <fe/fe_system.h>
#include <hp/fe_collection.h>
#include <hp/dof_handler.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_refinement.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <hp/dof_handler.h>
#include <hp/fe_values.h>


#include <fstream>


template <int dim>
void test ()
{
  Triangulation<dim>       triangulation;
  GridGenerator :: hyper_cube (triangulation, -0.5, 0.5);
  triangulation.refine_global(4);

  hp::FECollection<dim>    fe_collection;
  
          // test the behavior of the FE_Nothing element
          // in multicomponent systems (for, e.g., solving
          // Stokes' equations).
  
  fe_collection.push_back (FESystem<dim>(FE_RaviartThomas<dim>(0), 1,
                                         FE_DGQ<dim>(0), 1));
                                         
  fe_collection.push_back (FESystem<dim>(FE_Nothing<dim>(dim), 1,
                                         FE_Nothing<dim>(), 1));


  hp::DoFHandler<dim>      dof_handler (triangulation);

				   // loop over cells, and set cells
				   // within a circle to be of type
				   // FE_Nothing, while outside the
				   // circle to be of type RT(0)/DG(0)
       
  typename hp::DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

  for(; cell != endc; cell++)
  {
      Point<dim> center = cell->center();
      if(std::sqrt(center.square()) < 0.25 )
        cell->set_active_fe_index(1);
      else
        cell->set_active_fe_index(0);
  }

        // attempt to distribute dofs

  dof_handler.distribute_dofs (fe_collection);

  deallog << "   Number of active cells:       "
	  << triangulation.n_active_cells()
	  << std::endl
	  << "   Number of degrees of freedom: "
	  << dof_handler.n_dofs()
	  << std::endl;
}



int main ()
{
  std::ofstream logfile("fe_nothing_04/output");
  logfile.precision(2);

  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  //test<1> ();
  test<2> ();
  test<3> ();

  deallog << "OK" << std::endl;
}
