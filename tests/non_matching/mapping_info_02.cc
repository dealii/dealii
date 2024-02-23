// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

/*
 * Test the NonMatching::MappingInfo class together with FEPointEvaluation when
 * `reinit_cells` is called with a vector of vector of qpoints, rather than a
 * vector of quadratures.
 */

#include <deal.II/base/function_signed_distance.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/matrix_free/fe_point_evaluation.h>

#include <deal.II/non_matching/fe_values.h>
#include <deal.II/non_matching/mapping_info.h>
#include <deal.II/non_matching/mesh_classifier.h>
#include <deal.II/non_matching/quadrature_generator.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

void
test_vec_of_qpoints()
{
  constexpr unsigned int dim    = 2;
  constexpr unsigned int degree = 1;

  FE_Q<dim> fe_q(degree);

  Triangulation<dim> tria;

  MappingQ<dim> mapping(degree);

  GridGenerator::subdivided_hyper_cube(tria, 4);

  DoFHandler<dim> dof_handler(tria);

  dof_handler.distribute_dofs(fe_q);

  Functions::SignedDistance::Sphere<dim> level_set;

  Vector<double> level_set_vec(dof_handler.n_dofs());

  VectorTools::interpolate(dof_handler, level_set, level_set_vec);

  NonMatching::MeshClassifier<dim> mesh_classifier(dof_handler, level_set_vec);
  mesh_classifier.reclassify();

  hp::QCollection<1> q_collection((QGauss<1>(degree)));

  NonMatching::DiscreteQuadratureGenerator<dim> quadrature_generator(
    q_collection, dof_handler, level_set_vec);

  // FEPointEvaluation
  NonMatching::MappingInfo<dim> mapping_info_cell(
    mapping, update_values | update_gradients | update_JxW_values);

  std::vector<std::vector<Point<dim>>> unit_points_vector;
  for (const auto &cell : tria.active_cell_iterators())
    {
      quadrature_generator.generate(cell);
      unit_points_vector.emplace_back(
        quadrature_generator.get_inside_quadrature().get_points());
    }

  mapping_info_cell.reinit_cells(tria.active_cell_iterators(),
                                 unit_points_vector);

  Vector<double>      src(dof_handler.n_dofs());
  std::vector<double> solution_values_in(fe_q.dofs_per_cell);

  for (auto &v : src)
    v = random_value<double>();
  FEPointEvaluation<1, dim, dim, double> fe_point_cell(mapping_info_cell, fe_q);
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      cell->get_dof_values(src,
                           solution_values_in.begin(),
                           solution_values_in.end());

      fe_point_cell.reinit(cell->active_cell_index());
      fe_point_cell.evaluate(solution_values_in,
                             EvaluationFlags::values |
                               EvaluationFlags::gradients);
      for (const auto q : fe_point_cell.quadrature_point_indices())
        {
          deallog << fe_point_cell.get_value(q) << std::endl;
          deallog << fe_point_cell.get_gradient(q) << std::endl;
        }
    }
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  initlog();
  test_vec_of_qpoints();
}
