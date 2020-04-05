// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// Check downstream numbering


#include <deal.II/base/function_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"



template <int dim>
void
print_dofs(const DoFHandler<dim> &dof)
{
  const FiniteElement<dim> &           fe = dof.get_fe();
  std::vector<types::global_dof_index> v(fe.dofs_per_cell);
  std::shared_ptr<FEValues<dim>>       fevalues;

  if (fe.has_support_points())
    {
      Quadrature<dim> quad(fe.get_unit_support_points());
      fevalues = std::shared_ptr<FEValues<dim>>(
        new FEValues<dim>(fe, quad, update_quadrature_points));
    }

  for (typename DoFHandler<dim>::active_cell_iterator cell = dof.begin_active();
       cell != dof.end();
       ++cell)
    {
      Point<dim> p = cell->center();
      if (fevalues.get() != nullptr)
        fevalues->reinit(cell);

      cell->get_dof_indices(v);
      for (unsigned int i = 0; i < v.size(); ++i)
        if (fevalues.get() != nullptr)
          deallog << fevalues->quadrature_point(i) << '\t' << v[i] << std::endl;
        else
          deallog << p << '\t' << v[i] << std::endl;
      deallog << std::endl;
    }
}



template <int dim>
void
print_dofs(const DoFHandler<dim> &dof, unsigned int level)
{
  const FiniteElement<dim> &           fe = dof.get_fe();
  std::vector<types::global_dof_index> v(fe.dofs_per_cell);
  std::shared_ptr<FEValues<dim>>       fevalues;

  if (fe.has_support_points())
    {
      Quadrature<dim> quad(fe.get_unit_support_points());
      fevalues = std::shared_ptr<FEValues<dim>>(
        new FEValues<dim>(fe, quad, update_quadrature_points));
    }

  for (typename DoFHandler<dim>::cell_iterator cell = dof.begin(level);
       cell != dof.end(level);
       ++cell)
    {
      Point<dim> p = cell->center();
      if (fevalues.get() != nullptr)
        fevalues->reinit(cell);

      cell->get_mg_dof_indices(v);
      for (unsigned int i = 0; i < v.size(); ++i)
        if (fevalues.get() != nullptr)
          deallog << fevalues->quadrature_point(i) << '\t' << v[i] << std::endl;
        else
          deallog << p << '\t' << v[i] << std::endl;
      deallog << std::endl;
    }
}


template <int dim>
void
check_renumbering(DoFHandler<dim> &mgdof)
{
  const FiniteElement<dim> &element = mgdof.get_fe();
  DoFHandler<dim> &         dof     = mgdof;
  deallog << element.get_name() << std::endl;

  // Prepare a reordering of
  // components for later use
  std::vector<unsigned int> order(element.n_components());
  for (unsigned int i = 0; i < order.size(); ++i)
    order[i] = order.size() - i - 1;

  Tensor<1, dim> direction;
  for (unsigned int i = 0; i < dim; ++i)
    direction[i] = -5.0001 + 1.13 * i;

  deallog << std::endl << "Downstream numbering cell-wise" << std::endl;
  DoFRenumbering::downstream(dof, direction);
  print_dofs(dof);
  // Check level ordering
  for (unsigned int level = 0; level < dof.get_triangulation().n_levels();
       ++level)
    {
      deallog << "Level " << level << std::endl;
      DoFRenumbering::downstream(mgdof, level, direction);
      print_dofs(mgdof, level);
    }

  deallog << std::endl << "Downstream numbering dof-wise" << std::endl;
  DoFRenumbering::downstream(dof, direction, true);
  print_dofs(dof);
  // Check level ordering
  for (unsigned int level = 0; level < dof.get_triangulation().n_levels();
       ++level)
    {
      deallog << "Level " << level << std::endl;
      DoFRenumbering::downstream(mgdof, level, direction, true);
      print_dofs(mgdof, level);
    }
}


template <int dim>
void
check()
{
  Triangulation<dim> tr(Triangulation<dim>::limit_level_difference_at_vertices);
  GridGenerator::hyper_cube(tr, -1., 1.);
  tr.refine_global(1);
  tr.begin_active()->set_refine_flag();
  tr.execute_coarsening_and_refinement();
  tr.refine_global(3 - dim);

  DoFHandler<dim> mgdof(tr);

  FE_Q<dim>     q2(2);
  FE_DGQ<dim>   dgq1(1);
  FESystem<dim> system(q2, 2, dgq1, 1);

  mgdof.distribute_dofs(q2);
  mgdof.distribute_mg_dofs();
  check_renumbering(mgdof);
  mgdof.clear();

  mgdof.distribute_dofs(dgq1);
  mgdof.distribute_mg_dofs();
  check_renumbering(mgdof);
  mgdof.clear();

  mgdof.distribute_dofs(system);
  mgdof.distribute_mg_dofs();
  check_renumbering(mgdof);
  mgdof.clear();
}


int
main()
{
  initlog();
  deallog << std::setprecision(2) << std::fixed;

  deallog.push("1d");
  check<1>();
  deallog.pop();
  deallog.push("2d");
  check<2>();
  deallog.pop();
  deallog.push("3d");
  check<3>();
  deallog.pop();
}
