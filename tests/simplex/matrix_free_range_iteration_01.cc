// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test ShapeData for FE_SimplexP and QGaussSimplex

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include <ios>

#include "../tests.h"


template <int dim>
void
test()
{
  const unsigned int degree = 1;

  Triangulation<dim> tria;
  GridGenerator::subdivided_hyper_cube(tria, 4);

  hp::FECollection<dim> fe{FE_Q<2>(degree), FE_Q<2>(degree)};
  MappingQ<dim>         mapping(1);
  QGauss<dim>           quadrature(degree + 1);

  DoFHandler<dim> dof_handler(tria);

  unsigned int counter = 0;

  for (auto &cell : dof_handler.cell_iterators())
    if (counter++ < tria.n_cells() / 2)
      cell->set_active_fe_index(1);

  dof_handler.distribute_dofs(fe);


  typename MatrixFree<dim, double>::AdditionalData additional_data;
  additional_data.mapping_update_flags_boundary_faces =
    update_gradients | update_values;
  additional_data.mapping_update_flags_inner_faces =
    update_gradients | update_values;
  additional_data.initialize_mapping = false;
  additional_data.tasks_parallel_scheme =
    MatrixFree<dim, double>::AdditionalData::none;

  AffineConstraints<double> constraints;

  MatrixFree<dim, double> matrix_free;
  matrix_free.reinit(
    mapping, dof_handler, constraints, quadrature, additional_data);

  using VectorType = Vector<double>;

  VectorType src, dst;
  matrix_free.initialize_dof_vector(src);
  matrix_free.initialize_dof_vector(dst);

  std::vector<std::vector<CellId>>                          cells(2);
  std::vector<std::vector<std::pair<CellId, unsigned int>>> boundary_faces(2);

  std::vector<std::vector<std::vector<std::pair<CellId, CellId>>>> inner_faces(
    2);
  for (unsigned int i = 0; i < 2; ++i)
    inner_faces[i].resize(2);

  matrix_free.template loop<VectorType, VectorType>(
    [&](const auto &data, auto &dst, const auto &src, const auto range) {
      const auto i = data.get_cell_active_fe_index(range);

      for (unsigned int cell = range.first; cell < range.second; ++cell)
        for (unsigned int v = 0; v < data.n_active_entries_per_cell_batch(cell);
             ++v)
          cells[i].push_back(data.get_cell_iterator(cell, v)->id());
    },
    [&](const auto &data, auto &dst, const auto &src, const auto range) {
      const auto i = data.get_face_active_fe_index(range, true);
      const auto j = data.get_face_active_fe_index(range, false);

      for (unsigned int face = range.first; face < range.second; ++face)
        for (unsigned int v = 0; v < data.n_active_entries_per_face_batch(face);
             ++v)
          inner_faces[i][j].emplace_back(
            data.get_face_iterator(face, v, true).first->id(),
            data.get_face_iterator(face, v, false).first->id());
    },
    [&](const auto &data, auto &dst, const auto &src, const auto range) {
      const auto i = data.get_face_active_fe_index(range);

      for (unsigned int face = range.first; face < range.second; ++face)
        for (unsigned int v = 0; v < data.n_active_entries_per_face_batch(face);
             ++v)
          boundary_faces[i].emplace_back(
            data.get_face_iterator(face, v).first->id(),
            data.get_face_iterator(face, v).second);
    },
    dst,
    src);

  deallog << "CELL categories:" << std::endl;
  for (auto &i : cells)
    {
      std::sort(i.begin(), i.end());

      for (const auto &j : i)
        deallog << j << ' ';
      deallog << std::endl;
    }

  deallog << "BOUNDARY_FACE categories:" << std::endl;
  for (auto &i : boundary_faces)
    {
      std::sort(i.begin(), i.end());
      for (const auto &j : i)
        deallog << j.first << '@' << j.second << "   ";
      deallog << std::endl;
    }

  deallog << "INNER_FACE categories:" << std::endl;
  for (auto &k : inner_faces)
    {
      for (auto &i : k)
        {
          std::sort(i.begin(), i.end());
          for (const auto &j : i)
            deallog << j.first << '@' << j.second << "   ";
          deallog << std::endl;
        }
    }
}

int
main()
{
  initlog();

  test<2>();
}
