// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
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

// Check if intialization of FEFaceEvaluation works for hp + DG in 3D.

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
