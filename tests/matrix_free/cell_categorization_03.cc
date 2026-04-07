// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// Tests the categorization of cells in the SIMD vectorized case, where we set
// a separate category for each cell, using a lexicographic order of cells


#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/matrix_free/matrix_free.h>

#include "../tests.h"


template <std::size_t dim>
void
extract_lexicographic_index_of_cells(
  const std::array<int, dim>                       &lexicographic_index,
  const typename Triangulation<dim>::cell_iterator &cell,
  std::vector<std::pair<std::array<int, dim>,
                        typename Triangulation<dim>::active_cell_iterator>>
    &list_of_cells)
{
  if (cell->is_active() && cell->is_locally_owned())
    list_of_cells.emplace_back(lexicographic_index, cell);
  else if (cell->n_children() == Utilities::pow(2, dim))
    for (unsigned int child = 0, d2 = 0; d2 < (dim > 2 ? 2 : 1); ++d2)
      for (unsigned int d1 = 0; d1 < (dim > 1 ? 2 : 1); ++d1)
        for (unsigned int d0 = 0; d0 < 2; ++d0, ++child)
          {
            // Note that we switch the order for simpler indexing, placing the
            // z direction at the first array position
            std::array<int, dim> child_lexicographic;
            child_lexicographic[dim - 1] =
              lexicographic_index[dim - 1] * 2 + d0;
            if constexpr (dim > 1)
              child_lexicographic[dim - 2] =
                lexicographic_index[dim - 2] * 2 + d1;
            if constexpr (dim > 2)
              child_lexicographic[dim - 3] =
                lexicographic_index[dim - 3] * 2 + d2;
            extract_lexicographic_index_of_cells(child_lexicographic,
                                                 cell->child(child),
                                                 list_of_cells);
          }
  else
    Assert(false, ExcNotImplemented());
}


template <int dim, typename Number>
void
test()
{
  deallog << "Testing in " << dim << "d with type "
          << (std::is_same_v<Number, float> ? "float" : "double") << std::endl;
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0., 1.);
  tria.refine_global(6 - dim);
  FE_DGQ<dim>     fe(0);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);

  std::vector<std::pair<std::array<int, dim>,
                        typename Triangulation<dim>::active_cell_iterator>>
    list_of_cells;

  // hardcode for hypercube mesh; all other mesh kinds should be similar so do
  // not run more complicated case
  extract_lexicographic_index_of_cells(std::array<int, dim>{},
                                       tria.begin(0),
                                       list_of_cells);
  std::sort(list_of_cells.begin(), list_of_cells.end());

  typename MatrixFree<dim, Number>::AdditionalData mf_data;
  mf_data.cell_vectorization_category.resize(tria.n_active_cells());
  unsigned int counter = 0;
  for (const auto entry : list_of_cells)
    mf_data.cell_vectorization_category[entry.second->active_cell_index()] =
      counter++;
  mf_data.cell_vectorization_categories_strict = false;

  MatrixFree<dim, Number> mf;
  mf.reinit(
    MappingQ1<dim>(), dof, AffineConstraints<Number>(), QGauss<1>(1), mf_data);

  deallog << "Check order of cells with mergeable categories:" << std::endl;
  deallog << "Number of batches multiplied by SIMD width: "
          << mf.n_cell_batches() * VectorizedArray<Number>::size();
  for (unsigned int c = 0; c < mf.n_cell_batches(); ++c)
    {
      if (c * VectorizedArray<Number>::size() % 16 == 0)
        deallog << std::endl;
      for (unsigned int v = 0; v < VectorizedArray<Number>::size(); ++v)
        deallog << mf.get_cell_iterator(c, v)->id() << " ";
    }
  deallog << std::endl;

  mf_data.cell_vectorization_categories_strict = true;
  mf.reinit(
    MappingQ1<dim>(), dof, AffineConstraints<Number>(), QGauss<1>(1), mf_data);

  AssertDimension(mf.n_cell_batches(), tria.n_active_cells());

  deallog << "Check order of cells with strict categorization strategy:"
          << std::endl;
  deallog << "Number of batches: " << mf.n_cell_batches();
  for (unsigned int c = 0; c < mf.n_cell_batches(); ++c)
    {
      if (c % 16 == 0)
        deallog << std::endl;
      AssertDimension(mf.n_active_entries_per_cell_batch(c), 1);
      deallog << mf.get_cell_iterator(c, 0)->id() << " ";
    }
  deallog << std::endl;
  deallog << std::endl;
}


int
main()
{
  initlog();
  test<2, double>();
  test<3, double>();
  test<2, float>();
  test<3, float>();
}
