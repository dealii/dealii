//----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000 - 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------------

// Until version 1.50 of mg_dof_handler.cc, the
// MGDoFHandler::distribute_dofs function in 1d and 3d could not
// handle coarsened grids (unused vertices). Also, the
// MGDoFHandler::renumbering function could not handle coarsened grids
// (unused vertices, unused faces). Check that all this works now.

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/multigrid/mg_dof_handler.h>
#include <deal.II/multigrid/mg_dof_accessor.h>

#include <fstream>
#include <iomanip>
#include <iomanip>
#include <algorithm>

using namespace std;

template <int dim>
void check()
{
  FE_Q<dim> fe(3);
  deallog << fe.get_name() << std::endl;
  
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);
  typename Triangulation<dim>::active_cell_iterator cell=tria.begin_active();
  for (unsigned int i=0; i<GeometryInfo<dim>::max_children_per_cell; ++i, ++cell)
    cell->set_coarsen_flag();
  tria.execute_coarsening_and_refinement ();
  
  MGDoFHandler<dim> mg_dof_handler(tria);
  mg_dof_handler.distribute_dofs(fe);
  for (unsigned int level=0; level<tria.n_levels(); ++level)
    {
      const unsigned int n_dofs=mg_dof_handler.n_dofs(level);
      vector<unsigned int> new_numbers(n_dofs);
      for (unsigned int i=0; i<n_dofs; ++i)
	new_numbers[i]=n_dofs-1-i;

      mg_dof_handler.renumber_dofs(level, new_numbers);
    }
}


int main()
{
  std::ofstream logfile("renumbering_03/output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  check<1> ();
  check<2> ();
  check<3> ();

  deallog << "OK" << endl;
}
