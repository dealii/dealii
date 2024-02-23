// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Check the CellSimilarity with MappingFEField - since the field can be
// arbitrarily deformed, we should not get any similarity.

#include <deal.II/base/multithread_info.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe_field.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <int dim>
void
test(const unsigned int degree)
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);

  FESystem<dim>   fe_grid(FE_Q<dim>(degree), dim);
  DoFHandler<dim> dof_grid(tria);
  dof_grid.distribute_dofs(fe_grid);
  Vector<double> euler;
  euler.reinit(dof_grid.n_dofs());
  const ComponentMask mask(dim, true);
  VectorTools::get_position_vector(dof_grid, euler, mask);

  FE_Q<dim> fe(1);
  deallog << "FE=" << fe.get_name() << std::endl;

  MappingFEField<dim> mapping(dof_grid, euler, mask);
  deallog << "Mapping=MappingFEField(FE_Q(" << degree << "))" << std::endl;

  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);

  std::vector<Point<dim>> points(2);
  for (unsigned int d = 0; d < dim; ++d)
    points[0][d] = 0.1;
  for (unsigned int d = 0; d < dim; ++d)
    points[1][d] = 0.85;

  const Quadrature<dim> quadrature(points);
  FEValues<dim>         fe_values(mapping,
                          fe,
                          quadrature,
                          update_gradients | update_jacobians);

  const auto check = [&]() {
    for (const auto &cell : dof.active_cell_iterators())
      {
        fe_values.reinit(cell);

        deallog << "Jacobians: ";
        for (const auto q : fe_values.quadrature_point_indices())
          {
            deallog << "[ ";
            for (unsigned int d = 0; d < dim; ++d)
              for (unsigned int e = 0; e < dim; ++e)
                deallog << fe_values.jacobian(q)[d][e] << ' ';
            deallog << " ] ";
          }
        deallog << std::endl;
        deallog << "Derivatives of shape function: ";
        for (const auto q : fe_values.quadrature_point_indices())
          {
            deallog << "[ ";
            for (unsigned int d = 0; d < dim; ++d)
              deallog << fe_values.shape_grad(fe.dofs_per_cell / 2, q)[d]
                      << ' ';
            deallog << " ] ";
          }
        deallog << std::endl;
      }
  };

  deallog << "Undeformed configuration" << std::endl;
  check();

  euler(euler.size() - 1) += 0.01;
  deallog << "Last Euler point deformed" << std::endl;
  check();

  euler(0) += 0.01;
  deallog << "First and last Euler points deformed" << std::endl;
  check();
}



int
main()
{
  initlog();
  deallog << std::setprecision(8);
  MultithreadInfo::set_thread_limit(1);

  test<2>(1);
  test<2>(2);
  test<3>(1);
}
