// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test MatrixFree::get_matrix_free_cell_index() and
// MatrixFree::get_cell_iterator(), whether they indeed return
// inverse information.


#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>

#include <set>

#include "../tests.h"

template <int dim>
void
test()
{
  const unsigned int fe_degree = 1;

  using Number              = double;
  using VectorizedArrayType = VectorizedArray<Number>;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(3);

  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  MappingQ<dim> mapping(1);

  QGauss<1> quad(fe_degree + 1);

  AffineConstraints<Number> constraints;

  MatrixFree<dim, Number, VectorizedArrayType> matrix_free;
  matrix_free.reinit(mapping, dof_handler, constraints, quad);

  for (const auto &cell : tria.active_cell_iterators())
    {
      const auto mf_index = matrix_free.get_matrix_free_cell_index(cell);

      AssertThrow(cell == matrix_free.get_cell_iterator(
                            mf_index / VectorizedArrayType::size(),
                            mf_index % VectorizedArrayType::size()),
                  ExcInternalError());
    }

  deallog << "OK!" << std::endl;

  for (unsigned int i = 0; i < matrix_free.n_cell_batches(); ++i)
    for (unsigned int v = 0; v < matrix_free.n_active_entries_per_cell_batch(i);
         ++v)
      AssertDimension(matrix_free.get_matrix_free_cell_index(
                        matrix_free.get_cell_iterator(i, v)),
                      i * VectorizedArrayType::size() + v);

  deallog << "OK!" << std::endl;
}


int
main()
{
  initlog();

  test<2>();
}
