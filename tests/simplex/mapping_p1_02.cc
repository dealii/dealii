// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_p1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools_topology.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


// Test MappingP1's face values by comparing with MappingFE.

using namespace dealii;

template <int dim, int spacedim = dim>
void
test()
{
  deallog << "MappingP1<" << dim << ", " << spacedim << ">" << std::endl;

  const auto reference_cell = ReferenceCells::get_simplex<dim>();

  Triangulation<dim, spacedim> tria;
  if constexpr (dim == 1)
    GridGenerator::subdivided_hyper_cube(tria, 1, 1.0, 2.0, true);
  else if constexpr (dim == 2 && spacedim == 3)
    {
      Triangulation<dim> tria_2d;
      GridGenerator::subdivided_hyper_cube_with_simplices(
        tria_2d, 1, 1.0, 2.0, true);
      std::vector<Point<dim>>    points_2d;
      std::vector<CellData<dim>> cell_data;
      SubCellData                subcell_data;
      std::tie(points_2d, cell_data, subcell_data) =
        GridTools::get_coarse_mesh_description(tria_2d);

      std::vector<Point<spacedim>> points_3d;
      for (const auto &p : points_2d)
        points_3d.emplace_back(p[0], p[1], 1.0 + p[0] * p[0] + p[1] * p[1]);

      tria.create_triangulation(points_3d, cell_data, subcell_data);
    }
  else
    GridGenerator::subdivided_hyper_cube_with_simplices(
      tria, 1, 1.0, 2.0, true);
  FE_SimplexP<dim, spacedim> fe(1);
  MappingP1<dim, spacedim>   mapping;
  MappingFE<dim, spacedim>   mapping_fe(fe);

  QGaussSimplex<dim - 1> quadrature(2);

  UpdateFlags flags =
    update_quadrature_points | update_JxW_values | update_jacobians |
    update_jacobian_grads | update_inverse_jacobians |
    update_jacobian_pushed_forward_grads | update_jacobian_2nd_derivatives |
    update_jacobian_pushed_forward_2nd_derivatives |
    update_jacobian_3rd_derivatives |
    update_jacobian_pushed_forward_3rd_derivatives | update_boundary_forms;
  if (dim == spacedim)
    flags |= update_normal_vectors;

  FEFaceValues<dim, spacedim> fe_values(mapping, fe, quadrature, flags);
  FEFaceValues<dim, spacedim> fe_values_2(mapping_fe, fe, quadrature, flags);

  for (const auto &cell : tria.active_cell_iterators())
    {
      deallog << std::endl;

      deallog << "cell vertices" << std::endl;
      for (const unsigned int vertex_no : cell->vertex_indices())
        deallog << "  " << cell->vertex(vertex_no) << std::endl;
      deallog << std::endl;

      for (unsigned int face_no : cell->face_indices())
        {
          fe_values.reinit(cell, face_no);
          fe_values_2.reinit(cell, face_no);

          deallog << "JxW" << std::endl;
          for (unsigned int qp_n = 0; qp_n < quadrature.size(); ++qp_n)
            {
              deallog << "  " << fe_values.JxW(qp_n) << std::endl;

              Assert(std::abs((fe_values.JxW(qp_n) - fe_values_2.JxW(qp_n))) <
                       1e-12,
                     ExcInternalError());
            }
          deallog << std::endl;

          deallog << "quadrature points" << std::endl;
          for (unsigned int qp_n = 0; qp_n < quadrature.size(); ++qp_n)
            {
              const Point<dim> p = mapping.transform_real_to_unit_cell(
                cell, fe_values.quadrature_point(qp_n));
              deallog << "  " << p << " -> " << fe_values.quadrature_point(qp_n)
                      << std::endl;
              Assert(reference_cell.contains_point(p, 1e-12),
                     ExcInternalError());

              Assert((fe_values.quadrature_point(qp_n) -
                      fe_values_2.quadrature_point(qp_n))
                         .norm() < 1e-12,
                     ExcInternalError());
            }
          deallog << std::endl;

          deallog << "boundary forms" << std::endl;
          for (unsigned int qp_n = 0; qp_n < quadrature.size(); ++qp_n)
            {
              deallog << fe_values.normal_vector(qp_n) << std::endl;

              Assert((fe_values.boundary_form(qp_n) -
                      fe_values_2.boundary_form(qp_n))
                         .norm() < 1e-12,
                     ExcInternalError());
            }
          deallog << std::endl;

          if constexpr (dim == spacedim)
            {
              deallog << "normal vectors" << std::endl;
              for (unsigned int qp_n = 0; qp_n < quadrature.size(); ++qp_n)
                {
                  deallog << fe_values.normal_vector(qp_n) << std::endl;

                  Assert((fe_values.normal_vector(qp_n) -
                          fe_values_2.normal_vector(qp_n))
                             .norm() < 1e-12,
                         ExcInternalError());
                }
              deallog << std::endl;
            }

          deallog << "Jacobians" << std::endl;
          for (unsigned int qp_n = 0; qp_n < quadrature.size(); ++qp_n)
            {
              deallog << "  " << fe_values.jacobian(qp_n) << std::endl;

              for (unsigned int d = 0; d < spacedim; ++d)
                Assert((fe_values.jacobian(qp_n)[d] -
                        fe_values_2.jacobian(qp_n)[d])
                           .norm() < 1e-12,
                       ExcInternalError());
            }
          deallog << std::endl;

          deallog << "Jacobian gradients" << std::endl;
          for (unsigned int qp_n = 0; qp_n < quadrature.size(); ++qp_n)
            {
              deallog << "  all zero = "
                      << (fe_values.jacobian_grad(qp_n).norm() == 0)
                      << std::endl;

              for (unsigned int d = 0; d < spacedim; ++d)
                Assert((fe_values.jacobian_grad(qp_n)[d] -
                        fe_values_2.jacobian_grad(qp_n)[d])
                           .norm() < 1e-12,
                       ExcInternalError());
            }
          deallog << std::endl;

          deallog << "inverse jacobians" << std::endl;
          for (unsigned int qp_n = 0; qp_n < quadrature.size(); ++qp_n)
            {
              deallog << "  " << fe_values.inverse_jacobian(qp_n) << std::endl;

              for (unsigned int d = 0; d < dim; ++d)
                Assert((fe_values.inverse_jacobian(qp_n)[d] -
                        fe_values_2.inverse_jacobian(qp_n)[d])
                           .norm() < 1e-12,
                       ExcInternalError());
            }
          deallog << std::endl;

          deallog << "Jacobian push forward gradients" << std::endl;
          for (unsigned int qp_n = 0; qp_n < quadrature.size(); ++qp_n)
            {
              deallog << "  all zero = "
                      << (fe_values.jacobian_pushed_forward_grad(qp_n).norm() ==
                          0)
                      << std::endl;
            }
          deallog << std::endl;

          deallog << "Jacobian 2nd derivatives" << std::endl;
          for (unsigned int qp_n = 0; qp_n < quadrature.size(); ++qp_n)
            {
              deallog << "  all zero = "
                      << (fe_values.jacobian_2nd_derivative(qp_n).norm() == 0)
                      << std::endl;
            }
          deallog << std::endl;

          deallog << "Jacobian push forward 2nd derivatives" << std::endl;
          for (unsigned int qp_n = 0; qp_n < quadrature.size(); ++qp_n)
            {
              deallog << "  all zero = "
                      << (fe_values.jacobian_pushed_forward_2nd_derivative(qp_n)
                            .norm() == 0)
                      << std::endl;
            }
          deallog << std::endl;

          deallog << "Jacobian 3rd derivatives" << std::endl;
          for (unsigned int qp_n = 0; qp_n < quadrature.size(); ++qp_n)
            {
              deallog << "  all zero = "
                      << (fe_values.jacobian_3rd_derivative(qp_n).norm() == 0)
                      << std::endl;
            }
          deallog << std::endl;

          deallog << "Jacobian push forward 3rd derivatives" << std::endl;
          for (unsigned int qp_n = 0; qp_n < quadrature.size(); ++qp_n)
            {
              deallog << "  all zero = "
                      << (fe_values.jacobian_pushed_forward_3rd_derivative(qp_n)
                            .norm() == 0)
                      << std::endl;
            }
          deallog << std::endl;
        }
    }

  deallog << std::endl;
}

int
main()
{
  initlog();

  test<1>();
  test<1, 2>();
  test<2>();
  test<2, 3>();
  test<3>();

  deallog << "OK" << std::endl;
}
