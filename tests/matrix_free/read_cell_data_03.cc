// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// tests matrix-free read_cell_data with Table


#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>

#include "../tests.h"


template <int dim,
          typename Number              = double,
          typename VectorizedArrayType = VectorizedArray<Number>>
class Test
{
public:
  const constexpr static unsigned int n_lanes   = VectorizedArrayType().size();
  const constexpr static unsigned int fe_degree = 2;

  const constexpr static unsigned int size1 = 3;
  const constexpr static unsigned int size2 = 4;


  void
  run()
  {
    Triangulation<dim> tria(
      Triangulation<dim>::limit_level_difference_at_vertices);
    GridGenerator::hyper_cube(tria);

    for (unsigned int i = 0; i < 3; ++i)
      {
        tria.begin_active(i)->set_refine_flag();
        tria.execute_coarsening_and_refinement();
      }

    FE_Q<dim>       fe(fe_degree);
    DoFHandler<dim> dof_handler(tria);
    dof_handler.distribute_dofs(fe);
    dof_handler.distribute_mg_dofs();

    MappingQ<dim> mapping(1);

    QGauss<1> quad(fe_degree + 1);



    unsigned int level = numbers::invalid_unsigned_int;
    {
      AffineConstraints<Number> constraint;

      typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData
        additional_data;
      additional_data.mapping_update_flags                = update_values;
      additional_data.mapping_update_flags_inner_faces    = update_values;
      additional_data.mapping_update_flags_boundary_faces = update_values;
      additional_data.mg_level                            = level;

      MatrixFree<dim, Number, VectorizedArrayType> matrix_free;
      matrix_free.reinit(
        mapping, dof_handler, constraint, quad, additional_data);

      this->setup_vector(matrix_free);
      FEEvaluation<dim,
                   fe_degree,
                   fe_degree + 1,
                   1,
                   Number,
                   VectorizedArrayType>
        fe_eval(matrix_free);

      for (unsigned int cell_batch = 0;
           cell_batch < matrix_free.n_cell_batches();
           ++cell_batch)
        {
          fe_eval.reinit(cell_batch);

          std::array<unsigned int, n_lanes> cell_ids = fe_eval.get_cell_ids();


          Table<2, VectorizedArrayType> dst_2d;
          dst_2d.reinit(TableIndices<2>(size1, size2));
          fe_eval.read_cell_data(src_3d, dst_2d);

          // Check correctness of read data

          for (unsigned int i = 0; i < size1; ++i)
            for (unsigned int j = 0; j < size2; ++j)
              for (unsigned int lane = 0; lane < n_lanes; ++lane)
                {
                  if (cell_ids[lane] == numbers::invalid_unsigned_int)
                    continue;
                  const unsigned int cell_index = cell_ids[lane] / n_lanes;
                  const unsigned int array_lane = cell_ids[lane] % n_lanes;

                  const Number expected_value =
                    src_3d(cell_index, i, j)[array_lane];

                  if (dst_2d(i, j)[lane] != expected_value)
                    deallog << "Error on cell batch " << cell_batch
                            << " at dst_2d(" << i << "," << j << ") lane "
                            << lane << ": got " << dst_2d(i, j)[lane]
                            << ", expected " << expected_value << std::endl;
                }
        }

      deallog << std::endl;
    }
  }

private:
  void
  setup_vector(const MatrixFree<dim, Number, VectorizedArrayType> &data)
  {
    // Test case: 3D Table to 2D Table
    // src is a batch of 2D arrays, dst is a single 2D array
    unsigned int n_batches =
      data.n_cell_batches() + data.n_ghost_cell_batches();

    src_3d.reinit(TableIndices<3>(n_batches, size1, size2));

    // Fill src with test data
    for (unsigned int b = 0; b < n_batches; ++b)
      for (unsigned int i = 0; i < size1; ++i)
        for (unsigned int j = 0; j < size2; ++j)
          for (unsigned int lane = 0; lane < n_lanes; ++lane)
            src_3d(b, i, j)[lane] = b * 1000 + i * 100 + j * 10 + lane;
  }


  Table<3, VectorizedArrayType> src_3d;
};


int
main()
{
  initlog();
  {
    deallog.push("2d");
    Test<2> runner;
    runner.run();
    deallog.pop();
  }
  {
    deallog.push("3d");
    Test<3> runner;
    runner.run();
    deallog.pop();
  }
}
