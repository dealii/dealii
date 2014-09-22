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


// Check downstream numbering


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/multigrid/mg_dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <fstream>



template <int dim>
void
print_dofs (const DoFHandler<dim> &dof)
{
  const FiniteElement<dim> &fe = dof.get_fe();
  std::vector<types::global_dof_index> v (fe.dofs_per_cell);
  std_cxx11::shared_ptr<FEValues<dim> > fevalues;

  if (fe.has_support_points())
    {
      Quadrature<dim> quad(fe.get_unit_support_points());
      fevalues = std_cxx11::shared_ptr<FEValues<dim> >(new FEValues<dim>(fe, quad, update_q_points));
    }

  for (typename DoFHandler<dim>::active_cell_iterator cell=dof.begin_active();
       cell != dof.end(); ++cell)
    {
      Point<dim> p = cell->center();
      if (fevalues.get() != 0)
        fevalues->reinit(cell);

      cell->get_dof_indices (v);
      for (unsigned int i=0; i<v.size(); ++i)
        if (fevalues.get() != 0)
          deallog << fevalues->quadrature_point(i) << '\t' << v[i] << std::endl;
        else
          deallog << p << '\t' << v[i] << std::endl;
      deallog << std::endl;
    }
}



template <int dim>
void
print_dofs (const MGDoFHandler<dim> &dof, unsigned int level)
{
  const FiniteElement<dim> &fe = dof.get_fe();
  std::vector<types::global_dof_index> v (fe.dofs_per_cell);
  std_cxx11::shared_ptr<FEValues<dim> > fevalues;

  if (fe.has_support_points())
    {
      Quadrature<dim> quad(fe.get_unit_support_points());
      fevalues = std_cxx11::shared_ptr<FEValues<dim> >(new FEValues<dim>(fe, quad, update_q_points));
    }

  for (typename MGDoFHandler<dim>::cell_iterator cell=dof.begin(level);
       cell != dof.end(level); ++cell)
    {
      Point<dim> p = cell->center();
      if (fevalues.get() != 0)
        fevalues->reinit(cell);

      cell->get_mg_dof_indices (v);
      for (unsigned int i=0; i<v.size(); ++i)
        if (fevalues.get() != 0)
          deallog << fevalues->quadrature_point(i) << '\t' << v[i] << std::endl;
        else
          deallog << p << '\t' << v[i] << std::endl;
      deallog << std::endl;
    }
}


template <int dim>
void
check_renumbering(MGDoFHandler<dim> &mgdof)
{
  const FiniteElement<dim> &element = mgdof.get_fe();
  DoFHandler<dim> &dof = mgdof;
  deallog << element.get_name() << std::endl;

  // Prepare a reordering of
  // components for later use
  std::vector<unsigned int> order(element.n_components());
  for (unsigned int i=0; i<order.size(); ++i) order[i] = order.size()-i-1;

  Point<dim> direction;
  for (unsigned int i=0; i<dim; ++i)
    direction(i) = -5.0001+i;

  deallog << std::endl << "Downstream numbering cell-wise" << std::endl;
  DoFRenumbering::downstream(dof, direction);
  print_dofs (dof);
  // Check level ordering
  for (unsigned int level=0; level<dof.get_tria().n_levels(); ++level)
    {
      deallog << "Level " << level << std::endl;
      DoFRenumbering::downstream(mgdof, level, direction);
      print_dofs (mgdof, level);
    }

  deallog << std::endl << "Downstream numbering dof-wise" << std::endl;
  DoFRenumbering::downstream(dof, direction, true);
  print_dofs (dof);
  // Check level ordering
  for (unsigned int level=0; level<dof.get_tria().n_levels(); ++level)
    {
      deallog << "Level " << level << std::endl;
      DoFRenumbering::downstream(mgdof, level, direction, true);
      print_dofs (mgdof, level);
    }

}


template <int dim>
void
check ()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr, -1., 1.);
  tr.refine_global (1);
  tr.begin_active()->set_refine_flag ();
  tr.execute_coarsening_and_refinement ();
  tr.refine_global(3-dim);

  MGDoFHandler<dim> mgdof(tr);

  FE_Q<dim> q2(2);
  FE_DGQ<dim> dgq1(1);
  FESystem<dim> system (q2, 2, dgq1, 1);

  mgdof.distribute_dofs(q2);
  check_renumbering(mgdof);
  mgdof.clear();

  mgdof.distribute_dofs(dgq1);
  check_renumbering(mgdof);
  mgdof.clear();

  mgdof.distribute_dofs(system);
  check_renumbering(mgdof);
  mgdof.clear();
}


int main ()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (2);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console (0);

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
