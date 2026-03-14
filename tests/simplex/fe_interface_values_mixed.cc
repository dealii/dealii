// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2026 by the deal.II authors
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
// FEInterfaceValues for a mixed mesh, adopted from
// /tests/feinterace/fe_interace_values_02

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/q_collection.h>

#include <fstream>
#include <iostream>

#include "../tests.h"

#include "simplex_grids.h"

template <int dim, typename IteratorType>
void
test(const IteratorType &cell, FEInterfaceValues<dim> &fiv)
{
  for (const unsigned int f : cell->face_indices())
    if (!cell->at_boundary(f))
      {
        typename FEInterfaceValues<dim>::InterfaceData data(
          f,
          numbers::invalid_unsigned_int,
          cell->neighbor_of_neighbor(f),
          numbers::invalid_unsigned_int);
        data.fe_index               = cell->active_fe_index();
        data.q_index                = cell->active_fe_index();
        data.mapping_index          = cell->active_fe_index();
        data.fe_index_neighbor      = cell->neighbor(f)->active_fe_index();
        data.q_index_neighbor       = cell->neighbor(f)->active_fe_index();
        data.mapping_index_neighbor = cell->neighbor(f)->active_fe_index();

        fiv.reinit(cell, cell->neighbor(f), data);

        const unsigned int n_dofs = fiv.n_current_interface_dofs();
        Vector<double>     cell_vector(n_dofs);

        const auto &q_points = fiv.get_quadrature_points();

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


template <int dim>
void
test(unsigned int fe_degree)
{
  deallog << "Degree: " << fe_degree << std::endl;
  Triangulation<dim> tria;
  GridGenerator::subdivided_hyper_cube_with_simplices_mix(tria, 1);

  DoFHandler<dim> dofh(tria);

  for (const auto &cell : dofh.active_cell_iterators())
    if (cell->reference_cell().is_hyper_cube())
      cell->set_active_fe_index(1);
    else if (cell->reference_cell().is_simplex())
      cell->set_active_fe_index(0);
    else
      DEAL_II_ASSERT_UNREACHABLE();

  hp::FECollection<dim> fe;
  fe.push_back(FE_SimplexDGP<dim>(fe_degree));
  fe.push_back(FE_DGQ<dim>(fe_degree));

  dofh.distribute_dofs(fe);

  const hp::MappingCollection<dim> mapping(
    MappingFE<dim>(FE_SimplexDGP<dim>(1)), MappingQ1<dim>());
  const hp::QCollection<dim - 1> quadrature(QGaussSimplex<dim - 1>(fe_degree +
                                                                   1),
                                            QGauss<dim - 1>(fe_degree + 1));

  UpdateFlags update_flags = update_values | update_gradients |
                             update_quadrature_points | update_JxW_values;

  FEInterfaceValues<dim> fiv(mapping, fe, quadrature, update_flags);

  for (const auto &cell : dofh.active_cell_iterators())

    test(cell, fiv);
}



int
main()
{
  initlog();
  test<2>(1);
  test<2>(2);
  test<2>(3);
}
