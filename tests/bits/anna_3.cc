//----------------------------  anna_3.cc  ---------------------------
//    $Id$
//    Version: 
//
//    Copyright (C) 2002, 2003, 2004 by the deal.II authors and Anna Schneebeli
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  anna_3.cc  ---------------------------


// check some things about Nedelec elements, here that
// DoFTools::component_select and DoFTools::count_dofs_per_component
// works
//
// this program is a modified version of one by Anna Schneebeli,
// University of Basel

#include "../tests.h"
#include <base/logstream.h>
#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_system.h>		
#include <fe/fe_q.h>
#include <fe/fe_nedelec.h>
#include <fe/fe_base.h>
#include <dofs/dof_tools.h>
#include <dofs/dof_renumbering.h>
#include <iostream>
#include <fstream>


template <int dim>
class SystemTest 
{
  public:
    SystemTest ();
    void run ();    
				    
  private:
    void make_grid_and_dofs ();
    void check ();

				    
    Triangulation<dim>     triangulation;
    FESystem<dim>          fe;
    DoFHandler<dim>        dof_handler;

				   
};

template <int dim>
SystemTest<dim>::SystemTest () :
                fe (FE_Nedelec<dim>(1), 2,
                    FE_Q<dim>(1), 1),
		dof_handler (triangulation)
{}


template <int dim>
void SystemTest<dim>::make_grid_and_dofs ()
{			  
  GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global (0);
  deallog << "Number of active cells: " << triangulation.n_active_cells()
          << std::endl;
  deallog << "Total number of cells: " << triangulation.n_cells()
          << std::endl;
				  
  dof_handler.distribute_dofs (fe);
  deallog << "Number of degrees of freedom: " << dof_handler.n_dofs()
          << std::endl;
				  
}


template <int dim>
void SystemTest<dim>::check () 
{
  for (unsigned int c=0; c<fe.n_components(); ++c)
    {
      deallog << "Checking for component " << c << std::endl;
      std::vector<bool> x(fe.n_components(), false);
      x[c] = true;
      std::vector<bool> sel(dof_handler.n_dofs());
      DoFTools::extract_dofs (dof_handler, x, sel);

      for (unsigned int i=0; i<sel.size(); ++i)
        if (sel[i])
          deallog << "  DoF " << i << std::endl;
    };

  std::vector<unsigned int> dofs_per_component (fe.n_components(), 0U);
  DoFTools::count_dofs_per_component (dof_handler, dofs_per_component);
  deallog << "DoFs per component: ";
  for (unsigned int i=0; i<fe.n_components(); ++i)
    deallog << dofs_per_component[i] << ' ';
  deallog << std::endl;
}


template <int dim>
void SystemTest<dim>::run () 
{
  deallog << "************* " << dim << "D *************" << std::endl;
  make_grid_and_dofs ();
  check();

                                   // renumber degrees of freedom and
                                   // try again
  deallog << std::endl << "*** Renumbering ***" << std::endl;
  DoFRenumbering::component_wise (dof_handler);
  check ();
}

    

int main () 
{
  std::ofstream logfile("anna_3.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  SystemTest<2>().run();
  SystemTest<3>().run();  
  return 0;
}
