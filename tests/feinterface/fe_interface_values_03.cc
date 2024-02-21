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


// evaluate jump_in_shape_values(), average_of_shape_values(), shape_value() of
// FEInterfaceValues on an adaptive mesh

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>

#include <fstream>
#include <iostream>

#include "../tests.h"

#include "../test_grids.h"


template <int dim>
void
test(unsigned int fe_degree)
{
  Triangulation<dim> tria;
  TestGrids::hyper_line(tria, 2);

  tria.begin()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

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

  auto cell = dofh.begin(1);
  ++cell;

  for (const unsigned int f : GeometryInfo<dim>::face_indices())
    if (!cell->at_boundary(f))
      {
        if (!cell->neighbor_is_coarser(f))
          continue;

        auto nn = cell->neighbor_of_coarser_neighbor(f);
        fiv.reinit(cell,
                   f,
                   numbers::invalid_unsigned_int,
                   cell->neighbor(f),
                   nn.first,
                   nn.second);

        const unsigned int n_dofs = fiv.n_current_interface_dofs();
        Vector<double>     cell_vector(n_dofs);

        const auto &q_points = fiv.get_quadrature_points();
        for (unsigned int qpoint = 0; qpoint < q_points.size(); ++qpoint)
          deallog << "qpoint " << qpoint << ": " << q_points[qpoint]
                  << std::endl;

        for (unsigned int idx = 0; idx < n_dofs; ++idx)
          {
            const auto pair = fiv.interface_dof_to_dof_indices(idx);
            deallog << "  idx: " << idx
                    << " global: " << fiv.get_interface_dof_indices()[idx]
                    << " dof indices: " << static_cast<int>(pair[0]) << " | "
                    << static_cast<int>(pair[1]) << std::endl;
          }

        cell_vector = 0.0;
        for (unsigned int qpoint = 0; qpoint < q_points.size(); ++qpoint)
          for (unsigned int i = 0; i < n_dofs; ++i)
            cell_vector(i) +=
              fiv.shape_value(true, i, qpoint) * fiv.get_JxW_values()[qpoint];
        deallog << "shape_value(true): " << cell_vector << std::endl;

        cell_vector = 0.0;
        for (unsigned int qpoint = 0; qpoint < q_points.size(); ++qpoint)
          for (unsigned int i = 0; i < n_dofs; ++i)
            cell_vector(i) +=
              fiv.shape_value(false, i, qpoint) * fiv.get_JxW_values()[qpoint];
        deallog << "shape_value(false): " << cell_vector << std::endl;

        cell_vector = 0.0;
        for (unsigned int qpoint = 0; qpoint < q_points.size(); ++qpoint)
          for (unsigned int i = 0; i < n_dofs; ++i)
            cell_vector(i) += fiv.jump_in_shape_values(i, qpoint) *
                              fiv.get_JxW_values()[qpoint];
        deallog << "jump_in_shape_values(): " << cell_vector << std::endl;

        cell_vector = 0.0;
        for (unsigned int qpoint = 0; qpoint < q_points.size(); ++qpoint)
          for (unsigned int i = 0; i < n_dofs; ++i)
            cell_vector(i) += fiv.average_of_shape_values(i, qpoint) *
                              fiv.get_JxW_values()[qpoint];
        deallog << "average_of_shape_values(): " << cell_vector << std::endl;
      }
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
