// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2017 by the deal.II authors
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



/* Author: Wolfgang Bangerth, University of Heidelberg, 2001 */



#include "../tests.h"
#include <deal.II/base/function_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_system.h>




template <int dim>
void
print_dofs (const DoFHandler<dim> &dof)
{
  std::vector<types::global_dof_index> v (dof.get_fe().dofs_per_cell);
  for (typename DoFHandler<dim>::active_cell_iterator cell=dof.begin_active();
       cell != dof.end(); ++cell)
    {
      deallog << "Cell " << cell << " -- ";
      cell->get_dof_indices (v);
      for (unsigned int i=0; i<v.size(); ++i)
        deallog << v[i] << ' ';
      deallog << std::endl;
    }
}



template <int dim>
void
print_dofs (const DoFHandler<dim> &dof, unsigned int level)
{
  std::vector<types::global_dof_index> v (dof.get_fe().dofs_per_cell);
  for (typename DoFHandler<dim>::cell_iterator cell=dof.begin(level);
       cell != dof.end(level); ++cell)
    {
      deallog << "Cell " << cell << " -- ";
      cell->get_mg_dof_indices (v);
      for (unsigned int i=0; i<v.size(); ++i)
        deallog << v[i] << ' ';
      deallog << std::endl;
    }
}


template <int dim>
void
check_renumbering(DoFHandler<dim> &mgdof, bool discontinuous)
{
  const FiniteElement<dim> &element = mgdof.get_fe();
  DoFHandler<dim> &dof = mgdof;

  // Prepare a reordering of
  // components for later use
  std::vector<unsigned int> order(element.n_components());
  for (unsigned int i=0; i<order.size(); ++i) order[i] = order.size()-i-1;

  Tensor<1,dim> direction;
  for (unsigned int i=0; i<dim; ++i)
    direction[i] = std::pow(10.,static_cast<double>(i));

  // Check global ordering
  print_dofs (dof);

  if (discontinuous)
    {
      DoFRenumbering::downstream(dof, direction);
    }
  else
    {
      DoFRenumbering::Cuthill_McKee (dof);
      print_dofs (dof);
      DoFRenumbering::Cuthill_McKee (dof, true);
      print_dofs (dof);
      DoFRenumbering::Cuthill_McKee (dof, true, true);
      print_dofs (dof);
    }


  DoFRenumbering::component_wise (dof, order);
  print_dofs (dof);

  std::vector<bool> selected_dofs (dof.n_dofs(), false);
  for (unsigned int i=0; i<dof.n_dofs(); ++i) if (i%2==0) selected_dofs[i] = true;
  DoFRenumbering::sort_selected_dofs_back (dof, selected_dofs);
  print_dofs (dof);

  // Check level ordering
  for (unsigned int level=0; level<dof.get_triangulation().n_levels(); ++level)
    {
      print_dofs (mgdof, level);

// Reinsert after fixing
//        DoFRenumbering::Cuthill_McKee (mgdof, level);
//        print_dofs (mgdof, level);
//        DoFRenumbering::Cuthill_McKee (mgdof, level, true);
//        print_dofs (mgdof, level);

      if (discontinuous)
        {
          DoFRenumbering::downstream(mgdof, level, direction);
        }

      // renumber the non-MG part of the
      // DoFHandler again. this is probably
      // not what was meant in the first
      // place, but is compatible with the
      // behavior of the
      // DoFRenumbering::component_wise set
      // of functions before December 2005
      DoFRenumbering::component_wise (static_cast<DoFHandler<dim>&>(mgdof),
                                      order);
      print_dofs (mgdof, level);
    }

}


template <int dim>
void
check ()
{
  Triangulation<dim> tr(Triangulation<dim>::limit_level_difference_at_vertices);
  if (dim==2)
    GridGenerator::hyper_ball(tr, Point<dim>(), 1);
  else
    GridGenerator::hyper_cube(tr, -1,1);
  tr.reset_all_manifolds();
  tr.refine_global (1);
  tr.begin_active()->set_refine_flag ();
  tr.execute_coarsening_and_refinement ();
  if (dim==1)
    tr.refine_global(2);

  DoFHandler<dim> mgdof(tr);

  FESystem<dim> e1 (FE_Q<dim>(2), 2, FE_DGQ<dim>(1), 1);
  mgdof.distribute_dofs(e1);
  mgdof.distribute_mg_dofs(e1);
  check_renumbering(mgdof, false);
  mgdof.clear();

  FESystem<dim> e2 (FE_DGP<dim>(2), 2, FE_DGQ<dim>(1), 1);
  mgdof.distribute_dofs(e2);
  mgdof.distribute_mg_dofs(e2);
  check_renumbering(mgdof, true);
  mgdof.clear();
}


int main ()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (2);
  deallog << std::fixed;
  deallog.attach(logfile);

  deallog.push ("1d");
  check<1> ();
  deallog.pop ();
  deallog.push ("2d");
  check<2> ();
  deallog.pop ();
  deallog.push ("3d");
  check<3> ();
  deallog.pop ();
}
