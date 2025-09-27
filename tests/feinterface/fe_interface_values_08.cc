// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// evaluate average_gradient average_of_shape_hessians( jump_gradient of
// FEInterfaceValues

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>

#include <fstream>
#include <iostream>

#include "../tests.h"

#include "../test_grids.h"


template <int dim>
void
test(const FiniteElement<dim> &fe)
{
  Triangulation<dim> tria;
  TestGrids::hyper_line(tria, 2);

  DoFHandler<dim> dofh(tria);
  deallog << fe.get_name() << std::endl;
  dofh.distribute_dofs(fe);

  MappingQ<dim> mapping(1);
  UpdateFlags   update_flags = update_values | update_gradients |
                             update_hessians | update_quadrature_points |
                             update_JxW_values;

  FEInterfaceValues<dim> fiv(mapping,
                             fe,
                             QGauss<dim - 1>(fe.degree + 1),
                             update_flags);


  auto cell = dofh.begin();

  {
    unsigned int f = 0;
    Assert(cell->at_boundary(f), ExcInternalError());

    fiv.reinit(cell, f);

    const unsigned int n_dofs = fiv.n_current_interface_dofs();
    Vector<double>     cell_vector(n_dofs);

    const auto &q_points = fiv.get_quadrature_points();

    cell_vector = 0.0;
    for (unsigned int qpoint = 0; qpoint < q_points.size(); ++qpoint)
      for (unsigned int i = 0; i < n_dofs; ++i)
        cell_vector(i) += fiv.average_of_shape_gradients(i, qpoint).norm() *
                          fiv.get_JxW_values()[qpoint];
    deallog << "average_of_shape_gradients.norm(): " << cell_vector
            << std::endl;

    cell_vector = 0.0;
    for (unsigned int qpoint = 0; qpoint < q_points.size(); ++qpoint)
      for (unsigned int i = 0; i < n_dofs; ++i)
        cell_vector(i) += fiv.average_of_shape_hessians(i, qpoint).norm() *
                          fiv.get_JxW_values()[qpoint];
    deallog << "average_of_shape_hessians.norm(): " << cell_vector << std::endl;

    cell_vector = 0.0;
    for (unsigned int qpoint = 0; qpoint < q_points.size(); ++qpoint)
      for (unsigned int i = 0; i < n_dofs; ++i)
        cell_vector(i) += fiv.jump_in_shape_gradients(i, qpoint).norm() *
                          fiv.get_JxW_values()[qpoint];
    deallog << "jump_in_shape_gradients.norm(): " << cell_vector << std::endl;
  }


  for (const unsigned int f : GeometryInfo<dim>::face_indices())
    if (!cell->at_boundary(f))
      {
        fiv.reinit(cell,
                   f,
                   numbers::invalid_unsigned_int,
                   cell->neighbor(f),
                   cell->neighbor_of_neighbor(f),
                   numbers::invalid_unsigned_int);

        const unsigned int n_dofs = fiv.n_current_interface_dofs();
        Vector<double>     cell_vector(n_dofs);

        const auto &q_points = fiv.get_quadrature_points();

        cell_vector = 0.0;
        for (unsigned int qpoint = 0; qpoint < q_points.size(); ++qpoint)
          for (unsigned int i = 0; i < n_dofs; ++i)
            cell_vector(i) += fiv.average_of_shape_gradients(i, qpoint).norm() *
                              fiv.get_JxW_values()[qpoint];
        deallog << "average_of_shape_gradients.norm(): " << cell_vector
                << std::endl;

        cell_vector = 0.0;
        for (unsigned int qpoint = 0; qpoint < q_points.size(); ++qpoint)
          for (unsigned int i = 0; i < n_dofs; ++i)
            cell_vector(i) += fiv.average_of_shape_hessians(i, qpoint).norm() *
                              fiv.get_JxW_values()[qpoint];
        deallog << "average_of_shape_hessians.norm(): " << cell_vector
                << std::endl;

        cell_vector = 0.0;
        for (unsigned int qpoint = 0; qpoint < q_points.size(); ++qpoint)
          for (unsigned int i = 0; i < n_dofs; ++i)
            cell_vector(i) += fiv.jump_in_shape_gradients(i, qpoint).norm() *
                              fiv.get_JxW_values()[qpoint];
        deallog << "jump_in_shape_gradients.norm(): " << cell_vector
                << std::endl;
      }
}



int
main()
{
  initlog();
  test<2>(FE_Q<2>(1));
  test<2>(FE_DGQ<2>(0));
  test<2>(FE_DGQ<2>(1));
  test<2>(FE_DGQ<2>(2));
  test<3>(FE_Q<3>(1));
  test<3>(FE_DGQ<3>(0));
  test<3>(FE_DGQ<3>(1));
}
