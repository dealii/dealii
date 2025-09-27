// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test MatrixFreeTools::compute_diagonal() for differently refined
// triangulations

#include "compute_diagonal_util.h"

template <int dim,
          int fe_degree,
          int n_points                 = fe_degree + 1,
          int n_components             = 1,
          typename Number              = double,
          typename VectorizedArrayType = VectorizedArray<Number>>
void
test(const std::vector<std::shared_ptr<Triangulation<dim>>> &trias)
{
  for (const auto &tria : trias)
    {
      const FE_Q<dim>     fe_q(fe_degree);
      const FESystem<dim> fe(fe_q, n_components);

      // setup dof-handlers
      DoFHandler<dim> dof_handler(*tria);
      dof_handler.distribute_dofs(fe);

      AffineConstraints<Number> constraint;
      DoFTools::make_hanging_node_constraints(dof_handler, constraint);
      constraint.close();

      typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData
        additional_data;
      additional_data.mapping_update_flags = update_values | update_gradients;

      MappingQ<dim> mapping(1);
      QGauss<1>     quad(fe_degree + 1);

      MatrixFree<dim, Number, VectorizedArrayType> matrix_free;
      matrix_free.reinit(
        mapping, dof_handler, constraint, quad, additional_data);


      Test<dim, fe_degree, n_points, n_components, Number, VectorizedArrayType>
        test(matrix_free,
             constraint,
             [](FEEvaluation<dim,
                             fe_degree,
                             n_points,
                             n_components,
                             Number,
                             VectorizedArrayType> &phi) {
               phi.evaluate(EvaluationFlags::gradients);
               for (unsigned int q = 0; q < phi.n_q_points; ++q)
                 phi.submit_gradient(phi.get_gradient(q), q);
               phi.integrate(EvaluationFlags::gradients);
             });

      test.do_test();
    }
}

template <int dim>
std::vector<std::shared_ptr<Triangulation<dim>>>
create_triangulations()
{
  std::vector<std::shared_ptr<Triangulation<dim>>> result;

  // refine one cell of a 3^dim subdivided_hyper_cube
  for (unsigned int i = 0; i < Utilities::pow<unsigned int>(3, dim); ++i)
    {
      std::shared_ptr<Triangulation<dim>> tria =
        std::make_shared<Triangulation<dim>>();
      GridGenerator::subdivided_hyper_cube(*tria, 3);
      CellAccessor<dim, dim>(tria.get(), 0, i).set_refine_flag();
      tria->execute_coarsening_and_refinement();
      result.push_back(tria);
    }

  // refine one cell of a hyper_cross
  for (unsigned int i = 0; i < 4 * dim + 1; ++i)
    {
      std::shared_ptr<Triangulation<dim>> tria =
        std::make_shared<Triangulation<dim>>();
      std::vector<unsigned int> sizes(dim * 2, 2);
      GridGenerator::hyper_cross(*tria, sizes);
      CellAccessor<dim, dim>(tria.get(), 0, i).set_refine_flag();
      tria->execute_coarsening_and_refinement();
      result.push_back(tria);
    }

  // refine all cells of a 2^dim subdivided_hyper_cube but one
  for (unsigned int i = 0; i < Utilities::pow<unsigned int>(2, dim); ++i)
    {
      std::shared_ptr<Triangulation<dim>> tria =
        std::make_shared<Triangulation<dim>>();
      GridGenerator::subdivided_hyper_cube(*tria, 2);

      unsigned int counter = 0;

      for (const auto &cell : tria->active_cell_iterators())
        {
          if (i != counter)
            cell->set_refine_flag();

          ++counter;
        }

      result.push_back(tria);
    }

  return result;
}

int
main(int argc, char **argv)
{
  initlog();

  {
    const auto trias = create_triangulations<2>();

    // 2D - linear
    test<2, 1, 2, 1>(trias); // scalar
    test<2, 1, 2, 2>(trias); // vector

    // 2D - quadratic
    test<2, 2, 3, 1>(trias); // scalar
    test<2, 2, 3, 2>(trias); // vector
  }

  {
    const auto trias = create_triangulations<3>();

    // 3D - linear
    test<3, 1, 2, 1>(trias); // scalar
    test<3, 1, 2, 3>(trias); // vector

    // 3D - quadratic
    test<3, 2, 3, 1>(trias); // scalar
    test<3, 2, 3, 3>(trias); // vector
  }
}
