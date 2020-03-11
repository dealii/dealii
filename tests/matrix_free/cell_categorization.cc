// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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



// Tests the categorization of cells for vectorization


#include <deal.II/base/function.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/matrix_free/matrix_free.h>

#include "../tests.h"

#include "create_mesh.h"


template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  create_mesh(tria);

  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  typename Triangulation<dim>::active_cell_iterator cell, endc;
  cell = tria.begin_active();
  endc = tria.end();
  for (; cell != endc; ++cell)
    if (cell->center().norm() < 0.5)
      cell->set_refine_flag();
  tria.refine_global(1);
  tria.execute_coarsening_and_refinement();
  tria.begin(tria.n_levels() - 1)->set_refine_flag();
  tria.last()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  cell = tria.begin_active();
  for (unsigned int i = 0; i < 9 - 3 * dim; ++i)
    {
      cell                 = tria.begin_active();
      endc                 = tria.end();
      unsigned int counter = 0;
      for (; cell != endc; ++cell, ++counter)
        if (counter % (7 - i) == 0)
          cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }

  FE_DGQ<dim>     fe(1);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  AffineConstraints<double> constraints;
  constraints.close();

  MatrixFree<dim>                          mf_data;
  typename MatrixFree<dim>::AdditionalData data;
  data.tasks_parallel_scheme = MatrixFree<dim>::AdditionalData::none;
  data.mapping_update_flags_inner_faces =
    (update_gradients | update_JxW_values);
  data.mapping_update_flags_boundary_faces =
    (update_gradients | update_JxW_values);

  data.cell_vectorization_category.resize(tria.n_active_cells());
  for (const auto &cell : tria.active_cell_iterators())
    if (cell->is_locally_owned())
      data.cell_vectorization_category[cell->active_cell_index()] =
        static_cast<unsigned int>(cell->center()[1] * 10.);

  data.cell_vectorization_categories_strict = false;
  mf_data.reinit(dof, constraints, QGauss<1>(2), data);

  deallog << "Number of cell batches: " << mf_data.n_macro_cells() << std::endl;
  for (unsigned int i = 0; i < mf_data.n_macro_cells(); ++i)
    for (unsigned int c = 0; c < mf_data.n_components_filled(i); ++c)
      deallog << mf_data.get_cell_iterator(i, c)->id() << " with "
              << mf_data.get_cell_category(i) << std::endl;
  deallog << std::endl;

  data.cell_vectorization_categories_strict = true;
  mf_data.reinit(dof, constraints, QGauss<1>(2), data);
  deallog << "Number of cell batches: " << mf_data.n_macro_cells() << std::endl;
  for (unsigned int i = 0; i < mf_data.n_macro_cells(); ++i)
    for (unsigned int c = 0; c < mf_data.n_components_filled(i); ++c)
      deallog << mf_data.get_cell_iterator(i, c)->id() << " with "
              << mf_data.get_cell_category(i) << std::endl;
  deallog << std::endl;

  // set all categories to a fixed large number to make sure the memory
  // allocation does not depend on the actual number
  std::fill(data.cell_vectorization_category.begin(),
            data.cell_vectorization_category.end(),
            100000000);

  data.cell_vectorization_categories_strict = false;
  mf_data.reinit(dof, constraints, QGauss<1>(2), data);
  deallog << "Number of cell batches: " << mf_data.n_macro_cells() << std::endl;
  for (unsigned int i = 0; i < mf_data.n_macro_cells(); ++i)
    for (unsigned int c = 0; c < mf_data.n_components_filled(i); ++c)
      deallog << mf_data.get_cell_iterator(i, c)->id() << " with "
              << mf_data.get_cell_category(i) << std::endl;
  deallog << std::endl;
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc,
                                            argv,
                                            testing_max_num_threads());
  MPILogInitAll                    log;
  test<2>();
  test<3>();
}
