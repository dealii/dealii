// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2019 by the deal.II authors
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


// evaluate jump(), average(), shape_value() of FEInterfaceValues for
// at_boundary()

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>

#include <fstream>
#include <iostream>

#include "../tests.h"


template <int dim>
void
test(unsigned int fe_degree)
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);

  DoFHandler<dim> dofh(tria);
  FE_DGQ<dim>     fe(fe_degree);
  deallog << fe.get_name() << std::endl;
  dofh.distribute_dofs(fe);

  MappingQ<dim> mapping(1);
  UpdateFlags   update_flags = update_values | update_gradients |
                             update_quadrature_points | update_JxW_values;

  FEInterfaceValues<dim> fiv(mapping,
                             fe,
                             QGauss<dim - 1>(fe.degree + 1),
                             update_flags);


  auto               cell = dofh.begin();
  const unsigned int f    = 2;
  fiv.reinit(cell, f);

  const unsigned int n_dofs = fiv.n_current_interface_dofs();
  Vector<double>     cell_vector(n_dofs);

  const auto &q_points = fiv.get_quadrature_points();

  cell_vector = 0.0;
  for (unsigned int qpoint = 0; qpoint < q_points.size(); ++qpoint)
    for (unsigned int i = 0; i < n_dofs; ++i)
      cell_vector(i) +=
        fiv.shape_value(true, i, qpoint) * fiv.get_JxW_values()[qpoint];
  deallog << "shape_value(true): " << cell_vector;

  cell_vector = 0.0;
  for (unsigned int qpoint = 0; qpoint < q_points.size(); ++qpoint)
    for (unsigned int i = 0; i < n_dofs; ++i)
      cell_vector(i) +=
        fiv.shape_value(false, i, qpoint) * fiv.get_JxW_values()[qpoint];
  deallog << "shape_value(false): " << cell_vector;

  cell_vector = 0.0;
  for (unsigned int qpoint = 0; qpoint < q_points.size(); ++qpoint)
    for (unsigned int i = 0; i < n_dofs; ++i)
      cell_vector(i) += fiv.jump(i, qpoint) * fiv.get_JxW_values()[qpoint];
  deallog << "jump(): " << cell_vector;

  cell_vector = 0.0;
  for (unsigned int qpoint = 0; qpoint < q_points.size(); ++qpoint)
    for (unsigned int i = 0; i < n_dofs; ++i)
      cell_vector(i) += fiv.average(i, qpoint) * fiv.get_JxW_values()[qpoint];
  deallog << "average(): " << cell_vector;
}



int
main()
{
  initlog();
  test<2>(0);
  test<2>(1);
  test<3>(0);
  test<3>(1);
}
