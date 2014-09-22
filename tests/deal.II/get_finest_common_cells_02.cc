// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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



#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/multigrid/mg_dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>

#include <fstream>


template<int dim>
void test()
{
  // create 2 triangulations with the
  // same coarse grid, and refine
  // them differently
  Triangulation<dim> tria[2];

  GridGenerator::hyper_cube (tria[0]);
  GridGenerator::hyper_cube (tria[1]);

  tria[0].refine_global (2);
  tria[1].refine_global (2);

  tria[0].begin_active()->set_refine_flag();
  tria[0].execute_coarsening_and_refinement ();

  tria[1].last_active()->set_refine_flag();
  tria[1].execute_coarsening_and_refinement ();

  tria[1].last_active()->set_refine_flag();
  tria[1].execute_coarsening_and_refinement ();

  DoFHandler<dim> dh0 (tria[0]);
  DoFHandler<dim> dh1 (tria[1]);

  typedef
  std::list<std::pair<typename DoFHandler<dim>::cell_iterator,
      typename DoFHandler<dim>::cell_iterator> >
      CellList;

  const CellList cell_list
    = GridTools::get_finest_common_cells (dh0, dh1);
  for (typename CellList::const_iterator cell_pair = cell_list.begin();
       cell_pair != cell_list.end(); ++cell_pair)
    deallog << cell_pair->first << ' ' << cell_pair->second
            << std::endl;
}


int main()
{
  std::ofstream logfile ("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1>();
  test<2>();
  test<3>();
}

