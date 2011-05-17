//----------------------------  get_finest_common_cells_04.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  get_finest_common_cells_04.cc  ---------------------------
// check GridTools::get_finest_common_cells for MGDoFHandler arguments


#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/multigrid/mg_dof_handler.h>
#include <deal.II/multigrid/mg_dof_accessor.h>

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

  MGDoFHandler<dim> dh0 (tria[0]);
  MGDoFHandler<dim> dh1 (tria[1]);
  
  typedef
    std::list<std::pair<typename MGDoFHandler<dim>::cell_iterator,
                        typename MGDoFHandler<dim>::cell_iterator> >
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
  std::ofstream logfile ("get_finest_common_cells_04/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1>();
  test<2>();
  test<3>();
}

