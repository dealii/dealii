// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// tests MatrixFree::get_cell_category() and
// MatrixFree::get_cell_range_category()

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include "../tests.h"

template <int dim>
void
test()
{
  using Number              = double;
  using VectorizedArrayType = VectorizedArray<Number>;

  Triangulation<dim> tria;

  GridGenerator::hyper_cube(tria);
  tria.refine_global(3);

  // caterorization - not strict
  {
    DoFHandler<dim> dof_handler(tria);
    dof_handler.distribute_dofs(FE_Q<dim>(1));

    typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData data;
    data.tasks_parallel_scheme = MatrixFree<dim, double>::AdditionalData::none;
    data.mapping_update_flags  = update_quadrature_points;

    std::vector<unsigned int> cell_vectorization_category(
      tria.n_active_cells());
    for (unsigned int i = 0; i < cell_vectorization_category.size(); ++i)
      cell_vectorization_category[i] = i / 7;

    data.cell_vectorization_category          = cell_vectorization_category;
    data.cell_vectorization_categories_strict = false;

    MatrixFree<dim, Number, VectorizedArrayType> matrix_free;
    matrix_free.reinit(MappingQ1<dim>(),
                       dof_handler,
                       AffineConstraints<Number>(),
                       QGauss<dim>(2),
                       data);

    double dummy = 0;

    matrix_free.template cell_loop<double, double>(
      [&](const auto &matrix_free, auto &, const auto &, const auto range) {
        unsigned int max_category_range = 0;

        for (unsigned int cell = range.first; cell < range.second; ++cell)
          {
            unsigned int max_category = 0;
            for (unsigned int v = 0;
                 v < matrix_free.n_active_entries_per_cell_batch(cell);
                 ++v)
              {
                const auto category = cell_vectorization_category
                  [matrix_free.get_cell_iterator(cell, v)->active_cell_index()];

                // note: here and below we are using std::cout to output
                // information that is useful for debugging; however, since the
                // output depends on the hardware and might change if we
                // internally modify the way we loop over cells, we don't write
                // the data to the log but instead only check that the
                // properties that are described in the documentation are
                // fulfilled
                std::cout << category << " ";
                max_category = std::max(max_category, category);
              }

            std::cout << std::endl;

            AssertDimension(max_category,
                            matrix_free.get_cell_category(cell)); // test

            max_category_range = std::max(max_category_range, max_category);
          }

        AssertDimension(max_category_range,
                        matrix_free.get_cell_range_category(range)); // test
        std::cout << std::endl;
      },
      dummy,
      dummy);
  }
  deallog << "OK" << std::endl;

  // caterorization - strict
  {
    DoFHandler<dim> dof_handler(tria);
    dof_handler.distribute_dofs(FE_Q<dim>(1));

    typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData data;
    data.tasks_parallel_scheme = MatrixFree<dim, double>::AdditionalData::none;
    data.mapping_update_flags  = update_quadrature_points;
    data.mapping_update_flags_boundary_faces = update_quadrature_points;
    data.mapping_update_flags_inner_faces    = update_quadrature_points;
    data.hold_all_faces_to_owned_cells       = true;

    std::vector<unsigned int> cell_vectorization_category(
      tria.n_active_cells());
    for (unsigned int i = 0; i < cell_vectorization_category.size(); ++i)
      cell_vectorization_category[i] = i / 7;

    data.cell_vectorization_category          = cell_vectorization_category;
    data.cell_vectorization_categories_strict = true;

    MatrixFree<dim, Number, VectorizedArrayType> matrix_free;
    matrix_free.reinit(MappingQ1<dim>(),
                       dof_handler,
                       AffineConstraints<Number>(),
                       QGauss<dim>(2),
                       data);

    double dummy = 0;

    matrix_free.template cell_loop<double, double>(
      [&](const auto &matrix_free, auto &, const auto &, const auto range) {
        unsigned int max_category_range = 0;

        for (unsigned int cell = range.first; cell < range.second; ++cell)
          {
            for (unsigned int v = 0;
                 v < matrix_free.n_active_entries_per_cell_batch(cell);
                 ++v)
              {
                const auto category = cell_vectorization_category
                  [matrix_free.get_cell_iterator(cell, v)->active_cell_index()];
                std::cout << category << " ";
                AssertDimension(category,
                                matrix_free.get_cell_category(cell)); // test
              }

            std::cout << std::endl;

            max_category_range =
              std::max(max_category_range, matrix_free.get_cell_category(cell));
          }

        AssertDimension(max_category_range,
                        matrix_free.get_cell_range_category(range)); // test
        std::cout << std::endl;
      },
      dummy,
      dummy);
  }
  deallog << "OK" << std::endl;


  // hp
  {
    DoFHandler<dim> dof_handler(tria);

    unsigned int counter = 0;
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        cell->set_active_fe_index((counter++) / 7);

    hp::FECollection<dim> fe_collection;

    for (unsigned int i = 0; i < (tria.n_cells() + 6) / 7; ++i)
      fe_collection.push_back(FE_Q<dim>(1));

    dof_handler.distribute_dofs(fe_collection);

    typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData data;
    data.tasks_parallel_scheme = MatrixFree<dim, double>::AdditionalData::none;
    data.mapping_update_flags  = update_quadrature_points;
    data.mapping_update_flags_boundary_faces = update_quadrature_points;
    data.mapping_update_flags_inner_faces    = update_quadrature_points;
    data.hold_all_faces_to_owned_cells       = true;

    MatrixFree<dim, Number, VectorizedArrayType> matrix_free;
    matrix_free.reinit(MappingQ1<dim>(),
                       dof_handler,
                       AffineConstraints<Number>(),
                       QGauss<dim>(2),
                       data);


    double dummy = 0;

    matrix_free.template cell_loop<double, double>(
      [&](const auto &matrix_free, auto &, const auto &, const auto range) {
        unsigned int max_category_range = 0;

        for (unsigned int cell = range.first; cell < range.second; ++cell)
          {
            for (unsigned int v = 0;
                 v < matrix_free.n_active_entries_per_cell_batch(cell);
                 ++v)
              {
                const auto category =
                  matrix_free.get_cell_iterator(cell, v)->active_fe_index();
                std::cout << category << " ";
                AssertDimension(category,
                                matrix_free.get_cell_category(cell)); // test
                AssertDimension(
                  category, matrix_free.get_cell_range_category(range)); // test
              }

            std::cout << std::endl;
          }
        std::cout << std::endl;
      },
      dummy,
      dummy);
  }
  deallog << "OK" << std::endl;
}

int
main(int argc, char **argv)
{
  initlog();

  test<2>();
}
