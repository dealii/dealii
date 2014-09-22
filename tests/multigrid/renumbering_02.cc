// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
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
#include <deal.II/dofs/dof_accessor.h>

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
    DoFRenumbering::downstream(mg_dof_handler, level, a);
}


int main()
{
  initlog(__FILE__);
  check<1> ();
  check<2> ();
  check<3> ();

  deallog << "OK" << endl;
}
