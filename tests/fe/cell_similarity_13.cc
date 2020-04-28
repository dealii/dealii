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



// Check the CellSimilarity with MappingManifold on a somewhat arbitrarily
// deformed shape - a cube with one refinement in flat coordinates where we
// afterwards apply a spherical manifold to deform things.

#include <deal.II/base/multithread_info.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_manifold.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>

#include "../tests.h"



template <int dim>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0.2, 1.);
  tria.refine_global(1);
  SphericalManifold<dim> manifold;
  tria.set_all_manifold_ids(0);
  tria.set_manifold(0, manifold);

  FE_Q<dim> fe(1);
  deallog << "FE=" << fe.get_name() << std::endl;

  MappingManifold<dim> mapping;
  deallog << "Mapping=MappingManifold" << std::endl;

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

  for (const auto &cell : dof.active_cell_iterators())
    {
      fe_values.reinit(cell);

      deallog << "Jacobians: ";
      for (const auto q : fe_values.quadrature_point_indices())
        {
          deallog << "[ ";
          for (unsigned int d = 0; d < dim; ++d)
            for (unsigned int e = 0; e < dim; ++e)
              deallog << fe_values.jacobian(q)[d][e] << " ";
          deallog << " ] ";
        }
      deallog << std::endl;
      deallog << "Derivatives of shape function: ";
      for (const auto q : fe_values.quadrature_point_indices())
        {
          deallog << "[ ";
          for (unsigned int d = 0; d < dim; ++d)
            deallog << fe_values.shape_grad(fe.dofs_per_cell / 2, q)[d] << " ";
          deallog << " ] ";
        }
      deallog << std::endl;
    }
}



int
main()
{
  initlog();
  deallog << std::setprecision(8);
  MultithreadInfo::set_thread_limit(1);

  test<2>();
  test<3>();
}
