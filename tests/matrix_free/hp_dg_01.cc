// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check if initialization of FEFaceEvaluation works for hp + DG in 3D.

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include "../tests.h"

template <int dim>
void
test()
{
  using VectorType = Vector<double>;

  Triangulation<dim> tria;

  std::vector<unsigned int> repetitions(dim, 1);
  repetitions[0] = 2;

  Point<dim> p1;
  Point<dim> p2;

  p2[0] = 2.0;
  for (unsigned int i = 1; i < dim; ++i)
    p2[i] = 1.0;

  GridGenerator::subdivided_hyper_rectangle(tria, repetitions, p1, p2);


  hp::FECollection<dim> fe;
  fe.push_back(FE_DGQ<dim>(1));
  fe.push_back(FE_DGQ<dim>(1));

  hp::QCollection<dim> quad;
  quad.push_back(QGauss<dim>(2));

  MappingQ1<dim> mapping;

  DoFHandler<dim> dof_handler(tria);

  for (const auto &cell : dof_handler.active_cell_iterators())
    cell->set_active_fe_index(cell->active_cell_index());

  dof_handler.distribute_dofs(fe);

  AffineConstraints<double> constraints;

  typename MatrixFree<dim>::AdditionalData ad;

  ad.mapping_update_flags_inner_faces = update_gradients;

  MatrixFree<dim> matrix_free;
  matrix_free.reinit(mapping, dof_handler, constraints, quad, ad);

  VectorType dst(dof_handler.n_dofs());
  VectorType src(dof_handler.n_dofs());

  matrix_free.template loop<VectorType, VectorType>(
    [](const auto &, auto &, const auto &, const auto) {},
    [](const auto &matrix_free, auto &, const auto &, const auto range) {
      FEFaceEvaluation<dim, -1, 0, 1, double> eval(matrix_free, range);
    },
    [](const auto &, auto &, const auto &, const auto) {},
    dst,
    src);
}

int
main()
{
  initlog();

  test<3>();

  deallog << "OK" << std::endl;
}
