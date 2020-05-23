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

// check MappingFEField when initialized on all multigrid levels

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe_field.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/multigrid/mg_transfer_matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

template <int dim, int spacedim>
void
test()
{
  Triangulation<dim, spacedim> tria(
    Triangulation<dim, spacedim>::limit_level_difference_at_vertices);
  GridGenerator::hyper_ball(tria);

  tria.refine_global(1);

  FESystem<dim, spacedim>   fe(FE_Q<dim, spacedim>(2), spacedim);
  DoFHandler<dim, spacedim> dh(tria);

  dh.distribute_dofs(fe);
  dh.distribute_mg_dofs();

  deallog << "dim, spacedim: " << dim << ", " << spacedim << std::endl
          << "cells: " << tria.n_active_cells() << ", dofs: " << dh.n_dofs()
          << std::endl;

  // Create a Mapping
  LinearAlgebra::distributed::Vector<double> map_vector(dh.n_dofs());
  VectorTools::get_position_vector(dh, map_vector);
  MGLevelObject<LinearAlgebra::distributed::Vector<double>> level_vectors(
    0, tria.n_global_levels() - 1);
  for (unsigned int level = 0; level < tria.n_global_levels(); ++level)
    level_vectors[level].reinit(dh.n_dofs(level));

  MGTransferMatrixFree<dim, double> transfer;
  transfer.build(dh);
  transfer.interpolate_to_mg(dh, level_vectors, map_vector);
  MappingFEField<dim,
                 spacedim,
                 LinearAlgebra::distributed::Vector<double>,
                 DoFHandler<dim>>
                       mapping(dh, level_vectors);
  MappingQGeneric<dim> mapping_ref(fe.degree);

  QGauss<dim>   quad(1);
  FEValues<dim> fe_values_ref(mapping_ref,
                              fe,
                              quad,
                              update_jacobians | update_quadrature_points);
  FEValues<dim> fe_values(mapping,
                          fe,
                          quad,
                          update_jacobians | update_quadrature_points);

  for (const auto &cell : tria.cell_iterators())
    {
      fe_values_ref.reinit(cell);
      fe_values.reinit(cell);

      if (fe_values_ref.quadrature_point(0).distance(
            fe_values.quadrature_point(0)) > 1e-12)
        deallog << "Mapped point should be "
                << fe_values_ref.quadrature_point(0) << " and is "
                << fe_values.quadrature_point(0) << std::endl;
      Tensor<2, dim> jac_ref = fe_values_ref.jacobian(0),
                     jac     = fe_values.jacobian(0);
      if ((jac_ref - jac).norm() > 1e-12)
        deallog << "Jacobian should be " << jac_ref << " and is " << jac
                << std::endl;
    }
  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();
  test<2, 2>();
  test<3, 3>();
}
