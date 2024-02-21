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


// Compute the right-hand side vector due to a Neumann boundary condition on
// different meshes (pure simplex, wedge, pyramid mesh).

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_wedge_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparse_matrix.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

#include "./simplex_grids.h"



template <int dim>
void
test(const unsigned int v, const unsigned int degree)
{
  Triangulation<dim> tria;

  std::shared_ptr<FiniteElement<dim>> fe;
  std::shared_ptr<Quadrature<dim>>    quad;
  std::shared_ptr<FiniteElement<dim>> fe_mapping;
  hp::QCollection<dim - 1>            face_quad;

  if (v == 0)
    {
      GridGenerator::subdivided_hyper_cube_with_simplices(tria,
                                                          dim == 2 ? 16 : 8);
      fe         = std::make_shared<FE_SimplexP<dim>>(degree);
      quad       = std::make_shared<QGaussSimplex<dim>>(degree + 1);
      fe_mapping = std::make_shared<FE_SimplexP<dim>>(1);
      face_quad  = hp::QCollection<dim - 1>{QGaussSimplex<dim - 1>(degree + 1)};
    }
  else if (v == 1)
    {
      GridGenerator::subdivided_hyper_cube_with_wedges(tria, dim == 2 ? 16 : 8);
      fe         = std::make_shared<FE_WedgeP<dim>>(degree);
      quad       = std::make_shared<QGaussWedge<dim>>(degree + 1);
      fe_mapping = std::make_shared<FE_WedgeP<dim>>(1);
      face_quad  = hp::QCollection<dim - 1>{QGaussSimplex<dim - 1>(degree + 1),
                                            QGaussSimplex<dim - 1>(degree + 1),
                                            QGauss<dim - 1>(degree + 1),
                                            QGauss<dim - 1>(degree + 1),
                                            QGauss<dim - 1>(degree + 1)};
    }
  else if (v == 2)
    {
      GridGenerator::subdivided_hyper_cube_with_pyramids(tria,
                                                         dim == 2 ? 16 : 8);
      fe         = std::make_shared<FE_PyramidP<dim>>(degree);
      quad       = std::make_shared<QGaussPyramid<dim>>(degree + 1);
      fe_mapping = std::make_shared<FE_PyramidP<dim>>(1);
      face_quad  = hp::QCollection<dim - 1>{QGauss<dim - 1>(degree + 1),
                                            QGaussSimplex<dim - 1>(degree + 1),
                                            QGaussSimplex<dim - 1>(degree + 1),
                                            QGaussSimplex<dim - 1>(degree + 1),
                                            QGaussSimplex<dim - 1>(degree + 1)};
    }
  else
    DEAL_II_NOT_IMPLEMENTED();

  MappingFE<dim> mapping(*fe_mapping);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(*fe);

  AffineConstraints<double> constraints;

  using VectorType = LinearAlgebra::distributed::Vector<double>;
  using number     = double;

  VectorType vec_mf, vec_mb;

  vec_mf.reinit(dof_handler.n_dofs());
  vec_mb.reinit(dof_handler.n_dofs());

  // compute vector with FEFaceEvaluation
  {
    typename MatrixFree<dim, double>::AdditionalData additional_data;
    additional_data.mapping_update_flags = update_gradients | update_values;
    additional_data.mapping_update_flags_boundary_faces =
      update_gradients | update_values;

    MatrixFree<dim, double> matrix_free;
    matrix_free.reinit(
      mapping, dof_handler, constraints, *quad, additional_data);

    matrix_free.template loop<VectorType, VectorType>(
      [&](const auto &data, auto &dst, const auto &src, const auto cell_range) {
        (void)data;
        (void)dst;
        (void)src;
        (void)cell_range;
      },
      [&](const auto &data, auto &dst, const auto &src, const auto face_range) {
        (void)data;
        (void)dst;
        (void)src;
        (void)face_range;
      },
      [&](const auto &data, auto &dst, const auto &src, const auto face_range) {
        FEFaceEvaluation<dim, -1, 0, 1, number> fe_eval(data, face_range, true);
        for (unsigned int face = face_range.first; face < face_range.second;
             face++)
          {
            fe_eval.reinit(face);
            for (unsigned int q = 0; q < fe_eval.n_q_points; ++q)
              fe_eval.submit_value(1.0, q);

            fe_eval.integrate_scatter(true, false, dst);
          }
      },
      vec_mf,
      vec_mf);
  }

  // compute vector with FEFaceValues
  {
    const UpdateFlags flag = update_JxW_values | update_values |
                             update_gradients | update_quadrature_points;
    FEValues<dim> fe_values(mapping, *fe, *quad, flag);

    FEFaceValues<dim> fe_face_values(mapping, *fe, face_quad, flag);

    const unsigned int dofs_per_cell = fe->dofs_per_cell;

    std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
    Vector<double>                       cell_rhs(dofs_per_cell);

    for (const auto &cell : dof_handler.cell_iterators())
      {
        if (!cell->is_locally_owned())
          continue;

        fe_values.reinit(cell);
        cell_rhs = 0;

        for (const auto &face : cell->face_iterators())
          if (face->at_boundary())
            {
              fe_face_values.reinit(cell, face);
              for (const auto q : fe_face_values.quadrature_point_indices())
                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  cell_rhs(i) +=
                    (1.0 *                              // 1.0
                     fe_face_values.shape_value(i, q) * // phi_i(x_q)
                     fe_face_values.JxW(q));            // dx
            }

        cell->get_dof_indices(dof_indices);

        constraints.distribute_local_to_global(cell_rhs, dof_indices, vec_mb);
      }
  }

#if false
  vec_mf.print(deallog.get_file_stream());
  vec_mb.print(deallog.get_file_stream());
#endif

  for (unsigned int i = 0; i < vec_mf.size(); ++i)
    Assert(std::abs(vec_mf[i] - vec_mb[i]) < 1e-8, ExcNotImplemented());

  deallog << " dim=" << dim << " degree=" << degree << ": OK" << std::endl;
}


int
main(int argc, char **argv)
{
  initlog();

  deallog.depth_file(2);

  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);

  for (unsigned int i = 0; i <= 2; ++i)
    {
      if (i == 0)
        deallog.push("SIMPLEX");
      else if (i == 1)
        deallog.push("WEDGE  ");
      else if (i == 2)
        deallog.push("PYRAMID");
      else
        DEAL_II_NOT_IMPLEMENTED();

      if (i == 0) // 2D makes only sense for simplex
        {
          test<2>(i, /*degree=*/1);
          test<2>(i, /*degree=*/2);
        }

      test<3>(i, /*degree=*/1);

      if (i != 2) // for pyramids: no quadratic elements implemented
        test<3>(i, /*degree=*/2);

      deallog.pop();
    }
}
