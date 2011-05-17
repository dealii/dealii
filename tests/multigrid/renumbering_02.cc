//----------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------------

// Until version 1.50 of mg_dof_handler.cc, the
// MGDoFHandler::renumbering function could not handle coarsened grids
// (unused cells). Check that this works now.

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
  FE_DGQ<dim> fe(1);
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
  Point<dim> a;
  a(0)=1;
  for (unsigned int level=0; level<tria.n_levels(); ++level)
    DoFRenumbering::downstream_dg(mg_dof_handler, level, a);
}


int main()
{
  std::ofstream logfile("renumbering_02/output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  check<1> ();
  check<2> ();
  check<3> ();

  deallog << "OK" << endl;
}
