// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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



// Check that graph coloring colors every cells.


#include "../tests.h"
#include <fstream>
#include <vector>

#include <deal.II/base/logstream.h>
#include <deal.II/base/graph_coloring.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/lac/vector.h>



template <int dim>
std::vector<types::global_dof_index> get_conflict_indices(
    const typename DoFHandler<dim>::active_cell_iterator &it)
{
  std::vector<types::global_dof_index> local_dof_indices(it->get_fe().dofs_per_cell);
  it->get_dof_indices(local_dof_indices);

  return local_dof_indices;
}


template <int dim>
void check ()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_shell (triangulation,
      Point<dim>(),
      1, 2,
      (dim==3) ? 96 : 12,
      true);

  triangulation.refine_global (3);

  for (unsigned int i=0; i<3; ++i)
  {
    Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
    for (unsigned int i=0; i< estimated_error_per_cell.size(); ++i)
      estimated_error_per_cell(i) = i;
    GridRefinement::
      refine_and_coarsen_fixed_fraction (triangulation,
          estimated_error_per_cell,
          0.3, 0.1);
    triangulation.execute_coarsening_and_refinement();
  }
  
  FE_Q<dim> fe(1);
  DoFHandler<dim> stokes_dof_handler (triangulation);
  stokes_dof_handler.distribute_dofs (fe);

  for (typename DoFHandler<dim>::active_cell_iterator cell=stokes_dof_handler.begin_active();
      cell != stokes_dof_handler.end(); ++cell)
    cell->clear_user_flag ();
  std::vector<std::vector<typename DoFHandler<dim>::active_cell_iterator> > coloring 
    = GraphColoring::make_graph_coloring(stokes_dof_handler.begin_active(),
        stokes_dof_handler.end(),
        std_cxx11::function<std::vector<types::global_dof_index> (
          const typename DoFHandler<dim>::active_cell_iterator &)>(&get_conflict_indices<dim>));

  for (unsigned int c=0; c<coloring.size(); ++c)
    for (unsigned int i=0; i<coloring[c].size(); ++i)
      coloring[c][i]->set_user_flag();

  for (typename DoFHandler<dim>::active_cell_iterator cell=stokes_dof_handler.begin_active();
      cell != stokes_dof_handler.end(); ++cell)
    Assert (cell->user_flag_set () == true, ExcInternalError());

  deallog<<"OK"<<std::endl;
}

int main (int argc, char *argv[])
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console (0);

  check<2> ();

  return 0;
}
