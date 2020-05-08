// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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

// test that MappingQEulerian knows how to compute centers of cells

#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q_eulerian.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>
#include <vector>

#include "../tests.h"



template <int dim>
class Displacement : public Function<dim>
{
public:
  Displacement()
    : Function<dim>(dim)
  {}

  double
  value(const Point<dim> &p, const unsigned int component) const
  {
    return p[component] * p[0] / 2;
  }
};


template <int dim>
void
test()
{
  deallog << "dim=" << dim << std::endl;

  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation, -1, 1);

  FESystem<dim>   fe(FE_Q<dim>(2), dim);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  Vector<double> displacements(dof_handler.n_dofs());

  VectorTools::interpolate(dof_handler, Displacement<dim>(), displacements);

  MappingQEulerian<dim> euler(2, dof_handler, displacements);

  // now the actual test
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      deallog << "Center: [" << cell->center() << "], with mapping ["
              << euler.get_center(cell, false) << "]" << std::endl;
      deallog << "Center: [" << cell->center() << "], with mapping + flag ["
              << euler.get_center(cell, true) << "]" << std::endl;
    }
}



int
main()
{
  initlog(1);

  test<1>();
  test<2>();
  test<3>();
}
