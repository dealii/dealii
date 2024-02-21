// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check DoFRenumbering::hierarchical changes nothing for a regular refined mesh

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

#include <sstream>

#include "../tests.h"


template <int dim, class stream>
void
print_dofs(const DoFHandler<dim> &dof, stream &out)
{
  const FiniteElement<dim>            &fe = dof.get_fe();
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
          out << fevalues->quadrature_point(i) << '\t' << v[i] << std::endl;
        else
          out << p << '\t' << v[i] << std::endl;
      out << std::endl;
    }
}



template <int dim>
void
check()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr, -1., 1.);
  tr.refine_global(1);

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);

  std::ostringstream o1, o2;

  print_dofs(dof, o1);
  deallog << "**" << std::endl;

  DoFRenumbering::hierarchical(dof);

  print_dofs(dof, o2);

  if (o1.str() == o2.str())
    deallog << "OK" << std::endl;
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
