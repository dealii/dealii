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


// Test MappingP1's cell values for FE quantities by comparing with MappingFE.

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
  FE_SimplexP<dim, spacedim> fe(3);
  MappingP1<dim, spacedim>   mapping;
  MappingFE<dim, spacedim>   mapping_fe(fe);

  QGaussSimplex<dim> quadrature(2);

  UpdateFlags flags =
    update_values | update_gradients | update_hessians | update_3rd_derivatives;

  FEValues<dim, spacedim> fe_values(mapping, fe, quadrature, flags);
  FEValues<dim, spacedim> fe_values_2(mapping_fe, fe, quadrature, flags);

  for (const auto &cell : tria.active_cell_iterators())
    {
      deallog << std::endl;
      fe_values.reinit(cell);
      fe_values_2.reinit(cell);

      deallog << "cell vertices" << std::endl;
      for (const unsigned int vertex_no : cell->vertex_indices())
        deallog << "  " << cell->vertex(vertex_no) << std::endl;
      deallog << std::endl;

      deallog << "shape_values" << std::endl;
      for (unsigned int qp_n = 0; qp_n < quadrature.size(); ++qp_n)
        {
          for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
            {
              // here and below, only print out values for a single shape
              // function so that the output isn't tens of thousands of lines
              // long
              if (i == fe.dofs_per_cell / 2)
                deallog << "  " << fe_values.shape_value(i, qp_n) << std::endl;
              Assert(std::abs(fe_values.shape_value(i, qp_n) -
                              fe_values_2.shape_value(i, qp_n)) < 1e-12,
                     ExcInternalError());
            }
        }
      deallog << std::endl;

      deallog << "shape_gradients" << std::endl;
      for (unsigned int qp_n = 0; qp_n < quadrature.size(); ++qp_n)
        {
          for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
            {
              if (i == fe.dofs_per_cell / 2)
                deallog << "  " << fe_values.shape_grad(i, qp_n) << std::endl;
              Assert((fe_values.shape_grad(i, qp_n) -
                      fe_values_2.shape_grad(i, qp_n))
                         .norm() < 1e-12,
                     ExcInternalError());
            }
        }
      deallog << std::endl;

      deallog << "shape_hessians" << std::endl;
      for (unsigned int qp_n = 0; qp_n < quadrature.size(); ++qp_n)
        {
          for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
            {
              if (i == fe.dofs_per_cell / 2)
                deallog << "  " << fe_values.shape_hessian(i, qp_n)
                        << std::endl;
              Assert((fe_values.shape_hessian(i, qp_n) -
                      fe_values_2.shape_hessian(i, qp_n))
                         .norm() < 1e-12,
                     ExcInternalError());
            }
        }
      deallog << std::endl;

      deallog << "shape_3rd_derivatives" << std::endl;
      for (unsigned int qp_n = 0; qp_n < quadrature.size(); ++qp_n)
        {
          for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
            {
              if (i == fe.dofs_per_cell / 2)
                deallog << "  " << fe_values.shape_3rd_derivative(i, qp_n)
                        << std::endl;
              Assert((fe_values.shape_3rd_derivative(i, qp_n) -
                      fe_values_2.shape_3rd_derivative(i, qp_n))
                         .norm() < 1e-08,
                     ExcInternalError());
            }
        }
      deallog << std::endl;
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
