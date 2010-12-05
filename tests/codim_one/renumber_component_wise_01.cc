
//----------------------------  template.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005, 2008, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  template.cc  ---------------------------


// test DoFRenumbering::component_wise for codim=1

#include "../tests.h"
#include <fstream>
#include <base/logstream.h>
#include <base/function_lib.h>
#include <grid/tria.h>
#include <grid/grid_in.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <fe/fe_q.h>
#include <fe/fe_system.h>
#include <dofs/dof_renumbering.h>

#include <string>

std::ofstream logfile("renumber_component_wise_01/output");

template <int dim, int spacedim>
void test(std::string filename)
{
  Triangulation<dim, spacedim> tria;
  GridIn<dim, spacedim> gi;
  gi.attach_triangulation (tria);
  std::ifstream in (filename.c_str());
  gi.read_ucd (in);

  deallog << tria.n_active_cells() << " active cells" << std::endl;

  FESystem<dim,spacedim> fe(FE_Q<dim,spacedim> (2), spacedim);
  DoFHandler<dim,spacedim> dof_handler (tria);
  dof_handler.distribute_dofs (fe);

  deallog << dof_handler.n_dofs() << " degrees of freedom" << std::endl;
  DoFRenumbering::component_wise (dof_handler);

  for (typename DoFHandler<dim,spacedim>::active_cell_iterator
	 cell = dof_handler.begin_active(); cell!=dof_handler.end(); ++cell)
    {
      std::vector<unsigned int> x (cell->get_fe().dofs_per_cell);
      cell->get_dof_indices (x);

      deallog << cell << std::endl;
      for (unsigned int i=0; i<x.size(); ++i)
	deallog << "  " << x[i] << std::endl;
    }
}



int main ()
{
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<2,3>("grids/square.inp");
  test<2,3>("grids/sphere_1.inp");

  return 0;
}

