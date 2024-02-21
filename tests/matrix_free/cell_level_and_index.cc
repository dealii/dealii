// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test MatrixFree::get_cell_level_and_index()

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"



template <int dim, typename number>
void
compare_indices(const MatrixFree<dim, number> *mf_data)
{
  const unsigned int n_batches = mf_data->n_cell_batches();
  for (unsigned int batch_no = 0; batch_no < n_batches; ++batch_no)
    {
      const unsigned int n_lanes_filled =
        mf_data->n_active_entries_per_cell_batch(batch_no);
      for (unsigned int lane = 0; lane < n_lanes_filled; ++lane)
        {
          const auto cell = mf_data->get_cell_iterator(batch_no, lane);
          const auto level_index_pair =
            mf_data->get_cell_level_and_index(batch_no, lane);
          std::cout << "(cell->level(), cell->index) = (" << cell->level()
                    << ", " << cell->index() << ") "
                    << "versus get_cell_level_and_index() = ("
                    << level_index_pair.first << ", " << level_index_pair.second
                    << ")\n";
          AssertThrow(cell->level() == level_index_pair.first,
                      ExcMessage("mismatching cell levels"));
          AssertThrow(cell->index() == level_index_pair.second,
                      ExcMessage("mismatching cell indices"));
        }
    }
}



template <int dim, int fe_degree, typename number = double>
void
test(const bool adaptive_ref = true)
{
  Triangulation<dim> tria(
    Triangulation<dim>::limit_level_difference_at_vertices);
  GridGenerator::hyper_cube(tria, 0, 1, true);
  tria.refine_global(1);
  if (adaptive_ref)
    {
      for (auto &cell : tria.active_cell_iterators())
        if (cell->is_locally_owned())
          if (cell->center().norm() < 0.5)
            cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();
      for (auto &cell : tria.active_cell_iterators())
        if (cell->is_locally_owned())
          if (cell->center()[0] < 0.2)
            cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }
  else
    {
      tria.refine_global(1);
    }

  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  dof.distribute_mg_dofs();

  AffineConstraints<double> constraints;
  constraints.close();
  const QGauss<1>                                  quad(fe_degree + 1);
  typename MatrixFree<dim, number>::AdditionalData additional_data;
  additional_data.tasks_parallel_scheme =
    MatrixFree<dim, number>::AdditionalData::none;
  auto mf_data = std::make_shared<MatrixFree<dim, number>>();

  {
    std::cout << "Compare active indices." << std::endl;
    mf_data->reinit(MappingQ1<dim>{}, dof, constraints, quad, additional_data);
    compare_indices(mf_data.get());
  }

  {
    std::cout << "Compare level indices." << std::endl;
    const unsigned int level = tria.n_global_levels() - 1;
    additional_data.mg_level = level;
    mf_data->reinit(MappingQ1<dim>{}, dof, constraints, quad, additional_data);
    compare_indices(mf_data.get());
  }
}



int
main(int argc, char **argv)
{
  test<2, 1>();
  test<2, 1, float>();
  test<2, 1>(false);
}
