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



// tests categorization of ghost cells

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

  parallel::shared::Triangulation<dim> tria(
    MPI_COMM_WORLD,
    ::Triangulation<dim>::none,
    true,
    parallel::shared::Triangulation<dim>::partition_custom_signal);

  tria.signals.create.connect([&]() {
    for (const auto &cell : tria.active_cell_iterators())
      if (cell->center()[1] < 0.5)
        cell->set_subdomain_id(0);
      else
        cell->set_subdomain_id(1);
  });

  GridGenerator::subdivided_hyper_rectangle(tria,
                                            {10, 10},
                                            {0.0, 0.0},
                                            {1.0, 1.0});

  // caterorization - not strict
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
      cell_vectorization_category[i] = i % 10;

    data.cell_vectorization_category          = cell_vectorization_category;
    data.cell_vectorization_categories_strict = false;

    MatrixFree<dim, Number, VectorizedArrayType> matrix_free;
    matrix_free.reinit(MappingQ1<dim>(),
                       dof_handler,
                       AffineConstraints<Number>(),
                       QGauss<dim>(2),
                       data);

    for (unsigned int cell = matrix_free.n_cell_batches();
         cell <
         matrix_free.n_cell_batches() + matrix_free.n_ghost_cell_batches();
         ++cell)
      {
        unsigned int category = 0;
        for (unsigned int v = 0;
             v < matrix_free.n_active_entries_per_cell_batch(cell);
             ++v)
          category = std::max(
            category,
            cell_vectorization_category[matrix_free.get_cell_iterator(cell, v)
                                          ->active_cell_index()]);

        AssertDimension(category, matrix_free.get_cell_category(cell));
      }
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
      cell_vectorization_category[i] = i % 10;

    data.cell_vectorization_category          = cell_vectorization_category;
    data.cell_vectorization_categories_strict = true;

    MatrixFree<dim, Number, VectorizedArrayType> matrix_free;
    matrix_free.reinit(MappingQ1<dim>(),
                       dof_handler,
                       AffineConstraints<Number>(),
                       QGauss<dim>(2),
                       data);

    for (unsigned int cell = matrix_free.n_cell_batches();
         cell <
         matrix_free.n_cell_batches() + matrix_free.n_ghost_cell_batches();
         ++cell)
      {
        for (unsigned int v = 0;
             v < matrix_free.n_active_entries_per_cell_batch(cell);
             ++v)
          AssertDimension(
            cell_vectorization_category[matrix_free.get_cell_iterator(cell, v)
                                          ->active_cell_index()],
            matrix_free.get_cell_category(cell));
      }
  }
  deallog << "OK" << std::endl;


  // hp
  {
    DoFHandler<dim> dof_handler(tria);

    unsigned int counter = 0;
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        cell->set_active_fe_index((counter++) % 10);

    hp::FECollection<dim> fe_collection;

    for (unsigned int i = 0; i < 10; ++i)
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

    for (unsigned int cell = matrix_free.n_cell_batches();
         cell <
         matrix_free.n_cell_batches() + matrix_free.n_ghost_cell_batches();
         ++cell)
      {
        for (unsigned int v = 0;
             v < matrix_free.n_active_entries_per_cell_batch(cell);
             ++v)
          AssertDimension(
            matrix_free.get_cell_iterator(cell, v)->active_fe_index(),
            matrix_free.get_cell_category(cell));
      }
  }
  deallog << "OK" << std::endl;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll log;

  AssertDimension(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD), 2);

  test<2>();
}
