// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test that the fast matrix-free hanging-node algorithm is also working on
// adaptively refined 2D mixed meshes.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/matrix_free/matrix_free.h>

#include "../tests.h"



int
main()
{
  initlog();

  const unsigned int dim    = 2;
  const unsigned int degree = 1;

  Triangulation<dim> tria_0, tria_1, tria_2, tria_3, tria;

  GridGenerator::subdivided_hyper_rectangle(tria_0,
                                            {1, 1},
                                            {0.0, 0.0},
                                            {0.5, 0.5});
  GridGenerator::subdivided_hyper_rectangle_with_simplices(tria_1,
                                                           {1, 1},
                                                           {0.5, 0.0},
                                                           {1.0, 0.5});
  GridGenerator::subdivided_hyper_rectangle_with_simplices(tria_2,
                                                           {1, 1},
                                                           {0.0, 0.5},
                                                           {0.5, 1.0});
  GridGenerator::subdivided_hyper_rectangle_with_simplices(tria_3,
                                                           {1, 1},
                                                           {0.5, 0.5},
                                                           {1.0, 1.0});

  GridGenerator::merge_triangulations({&tria_0, &tria_1, &tria_2, &tria_3},
                                      tria);

  auto cell = tria.begin();

  cell->set_refine_flag();
  cell++;
  cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  if (false)
    {
      GridOut       grid_out;
      std::ofstream out("mesh.vtk");
      grid_out.write_vtk(tria, out);
    }

  DoFHandler<dim> dof_handler(tria);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->reference_cell() == ReferenceCells::Triangle)
        cell->set_active_fe_index(0);
      else if (cell->reference_cell() == ReferenceCells::Quadrilateral)
        cell->set_active_fe_index(1);
      else
        DEAL_II_NOT_IMPLEMENTED();
    }

  const hp::MappingCollection<2> mapping(MappingFE<2>(FE_SimplexP<2>(1)),
                                         MappingFE<2>(FE_Q<2>(1)));
  const hp::FECollection<2>      fe(FE_SimplexP<2>{degree}, FE_Q<2>{degree});
  const hp::QCollection<2> quadrature_formula(QGaussSimplex<2>(degree + 1),
                                              QGauss<2>(degree + 1));

  dof_handler.distribute_dofs(fe);

  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();

  const auto print = [](const auto &label, const auto &matrix_free) {
    deallog << label << std::endl;

    for (unsigned int c = 0;
         c < matrix_free.get_dof_info(0).row_starts.size() - 1;
         ++c)
      {
        deallog
          << std::setw(3)
          << (matrix_free.get_dof_info(0)
                    .hanging_node_constraint_masks.size() == 0 ?
                0 :
                static_cast<unsigned int>(
                  matrix_free.get_dof_info(0).hanging_node_constraint_masks[c]))
          << " : ";

        for (unsigned int i = matrix_free.get_dof_info(0).row_starts[c].first;
             i < matrix_free.get_dof_info(0).row_starts[c + 1].first;
             ++i)
          deallog << std::setw(3) << matrix_free.get_dof_info(0).dof_indices[i]
                  << ' ';
        deallog << std::endl;
      }
    deallog << std::endl;
  };

  {
    typename MatrixFree<dim, double, VectorizedArray<double, 1>>::AdditionalData
      additional_data;
    additional_data.mapping_update_flags = update_gradients | update_values;

    MatrixFree<dim, double, VectorizedArray<double, 1>> matrix_free;
    matrix_free.reinit(
      mapping, dof_handler, constraints, quadrature_formula, additional_data);

    print("use_fast_hanging_node_algorithm = true", matrix_free);
  }
}
