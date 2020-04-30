// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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



// hp::FEValues objects can be copied, but that sometimes led to
// surprising results because it kept shared pointers underneath (this
// has been fixed).


#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nothing.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>

#include "../tests.h"



template <int dim>
void
test()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(1);

  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(FE_Nothing<dim>());

  hp::DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe_collection);

  deallog << "   Number of active cells:       "
          << triangulation.n_active_cells() << std::endl
          << "   Number of degrees of freedom: " << dof_handler.n_dofs()
          << std::endl;


  hp::QCollection<dim> quadrature_collection;
  quadrature_collection.push_back(QMidpoint<dim>());

  // Create one hp::FEValues object, and initialize it for the first
  // cell. When we later ask for the location of the quadrature point,
  // it better be the midpoint of that cell
  const auto        cell1 = dof_handler.begin_active();
  hp::FEValues<dim> hp_fe_values1(fe_collection,
                                  quadrature_collection,
                                  update_quadrature_points);
  hp_fe_values1.reinit(cell1);

  // Now make a copy of the object and initialize it for the second cell
  const auto        cell2 = ++dof_handler.begin_active();
  hp::FEValues<dim> hp_fe_values2(hp_fe_values1);
  hp_fe_values2.reinit(cell2);

  // Output the quadrature points computed by the two objects. They
  // ought to correspond to the cell centers of the first and second
  // cell.
  deallog << hp_fe_values1.get_present_fe_values().get_quadrature_points()[0]
          << " should be at " << cell1->center() << std::endl;

  deallog << hp_fe_values2.get_present_fe_values().get_quadrature_points()[0]
          << " should be at " << cell2->center() << std::endl;
}



int
main()
{
  initlog();

  test<1>();
  test<2>();
  test<3>();
}
