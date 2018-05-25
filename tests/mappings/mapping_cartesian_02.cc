// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2018 by the deal.II authors
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



// there used to be a bug in MappingCartesian, where the inverse jacobian
// was not computed correctly. Check that all works as intended.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_cartesian.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int dim>
void
check(const Triangulation<dim> &tria)
{
  MappingCartesian<dim> mapping;
  FE_Q<dim>             fe(1);
  DoFHandler<dim>       dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  QGauss<dim> quadrature(1);

  UpdateFlags update_flags =
    update_quadrature_points | update_JxW_values | update_jacobians |
    update_jacobian_grads | update_jacobian_pushed_forward_grads |
    update_jacobian_2nd_derivatives |
    update_jacobian_pushed_forward_2nd_derivatives |
    update_jacobian_3rd_derivatives |
    update_jacobian_pushed_forward_3rd_derivatives | update_inverse_jacobians |
    // Transformation dependence
    update_covariant_transformation | update_contravariant_transformation |
    update_transformation_values | update_transformation_gradients |
    // Volume data
    update_volume_elements;

  FEValues<dim> fe_values(mapping, fe, quadrature, update_flags);

  fe_values.reinit(dof_handler.begin_active());

  for (unsigned int i = 0; i < dim; ++i)
    deallog << fe_values.inverse_jacobian(0)[i] << std::endl;

  deallog << std::endl;

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  {
    Triangulation<2> coarse_grid;
    GridGenerator::hyper_cube(coarse_grid);
    check(coarse_grid);
  }

  {
    Triangulation<3> coarse_grid;
    GridGenerator::hyper_cube(coarse_grid);
    check(coarse_grid);
  }
}
