// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
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
 * Test NonMatching::MappingInfo if reinit with quadratures which contain JxW as
 * weights.
 */

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/matrix_free/fe_point_evaluation.h>

#include <deal.II/non_matching/mapping_info.h>

#include "../tests.h"

void
test()
{
  constexpr unsigned int dim          = 2;
  constexpr unsigned int spacedim     = 2;
  constexpr unsigned int degree       = 1;
  constexpr unsigned int n_components = 1;
  using Number                        = double;

  Triangulation<dim, spacedim> tria;
  GridGenerator::subdivided_hyper_cube(tria, 4);

  FE_Q<dim>     fe(degree);
  MappingQ<dim> mapping(degree);

  NonMatching::MappingInfo<dim>::AdditionalData additional_data;
  additional_data.use_global_weights = true;
  NonMatching::MappingInfo<dim> mapping_info(mapping,
                                             update_JxW_values,
                                             additional_data);

  deallog << "Check JxW faces..." << std::endl;
  {
    // 1) build vector of quadratures
    std::vector<std::vector<Quadrature<dim - 1>>> quad_vec;
    // prescribe JxW (there is no meaning in the actual values, they just have
    // to stay the same when fetched with FEPointEvaluation)
    double JxW = 1.0;
    for (const auto &cell : tria.active_cell_iterators())
      {
        std::vector<Quadrature<dim - 1>> quad;
        for (auto f : cell->face_indices())
          {
            dealii::QGauss<dim - 1> face_quadrature(degree + 1);
            std::vector<double> weights(face_quadrature.get_weights().size());
            for (auto &w : weights)
              {
                w = JxW;
                JxW += 1.0;
              }

            quad.emplace_back(
              Quadrature<dim - 1>(face_quadrature.get_points(), weights));
          }
        quad_vec.push_back(quad);
      }

    // 2) reinit mapping info
    mapping_info.reinit_faces(tria.active_cell_iterators(), quad_vec);

    FEFacePointEvaluation<n_components, dim, spacedim, Number> fe_point_eval(
      mapping_info, fe);

    // 3) print JxW
    for (const auto &cell : tria.active_cell_iterators())
      {
        for (auto f : cell->face_indices())
          {
            fe_point_eval.reinit(cell->active_cell_index(), f);
            for (const unsigned int q :
                 fe_point_eval.quadrature_point_indices())
              deallog << fe_point_eval.JxW(q) << std::endl;
          }
      }
  }

  deallog << "\n\nCheck JxW cells..." << std::endl;
  {
    // 1) build vector of quadratures
    std::vector<Quadrature<dim>> quad_vec;
    // prescribe JxW (there is no meaning in the actual values, they just have
    // to stay the same when fetched with FEPointEvaluation)
    double JxW = 1.0;
    for (const auto &cell : tria.active_cell_iterators())
      {
        dealii::QGauss<dim> cell_quadrature(degree + 1);
        std::vector<double> weights(cell_quadrature.get_weights().size());
        for (auto &w : weights)
          {
            w = JxW;
            JxW += 1.0;
          }

        quad_vec.emplace_back(
          Quadrature<dim>(cell_quadrature.get_points(), weights));
      }

    // 2) reinit mapping info
    mapping_info.reinit_cells(tria.active_cell_iterators(), quad_vec);
    FEPointEvaluation<n_components, dim, spacedim, Number> fe_point_eval(
      mapping_info, fe);

    // 3) print JxW
    for (const auto &cell : tria.active_cell_iterators())
      {
        fe_point_eval.reinit(cell->active_cell_index());
        for (const unsigned int q : fe_point_eval.quadrature_point_indices())
          deallog << fe_point_eval.JxW(q) << std::endl;
      }
  }


  deallog << std::endl;
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  initlog();
  test();
}
