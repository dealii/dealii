// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// like hp_constraints_03 but on mixed meshes
// check that computation of hp-constraints works for FE_SimplexP, FE_PyramidP
// and FE_WedgeP elements by interpolating a function on a mixed mesh

#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_wedge_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/reference_cell.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <string>
#include <vector>

#include "../tests.h"

#include "simplex_grids.h"

template <int dim>
void
test_interpolation(const hp::FECollection<dim> &fe,
                   const Triangulation<dim>    &triangulation)
{
  Assert(fe.size() == 2, ExcInternalError());
  Assert(triangulation.get_reference_cells().size() == 2, ExcInternalError());

  deallog << "Testing interpolation " << fe[0].get_name() << " vs. "
          << fe[1].get_name() << std::endl;


  DoFHandler<dim> dof_handler(triangulation);
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->reference_cell() == fe[0].reference_cell())
        cell->set_active_fe_index(0);
      else if (cell->reference_cell() == fe[1].reference_cell())
        cell->set_active_fe_index(1);
      else
        DEAL_II_ASSERT_UNREACHABLE();
    }

  dof_handler.distribute_dofs(fe);

  Vector<double> interpolant(dof_handler.n_dofs());
  Vector<float>  error(triangulation.n_active_cells());

  const unsigned int min_degree = std::min(fe[0].degree, fe[1].degree);
  // interpolate monomials of the for x1^a*x2^b*x3^c
  // tets can interpolate a + b + c <= degree exactly
  // wedges can interpolate a + b <= degree and c <= degree exactly
  // pyramids can interpolate a + b + c <= degree exactly
  // hex can interpolate a, b, c <= degree exactly
  // so now figure out what to use
  // on a mixed mesh use the most restricting space
  std::vector<Tensor<1, dim>> exponents;
  for (unsigned int a = 0; a <= min_degree; ++a)
    for (unsigned int b = 0; b <= min_degree; ++b)
      for (unsigned int c = 0; c <= min_degree; ++c)
        {
          Tensor<1, dim> exponent;

          bool use_exponent = true;
          for (unsigned int i = 0; i < fe.size(); ++i)
            {
              if (fe[i].reference_cell().is_simplex() ||
                  fe[i].reference_cell() == ReferenceCells::Pyramid)
                {
                  if (a + b + c > fe[i].degree)
                    use_exponent = false;
                }
              else if (fe[i].reference_cell() == ReferenceCells::Wedge)
                {
                  if (a + b > fe[i].degree)
                    use_exponent = false;
                }
            }

          if (use_exponent)
            {
              exponent[0] = a;
              if constexpr (dim > 1)
                exponent[1] = b;
              if constexpr (dim > 2)
                exponent[2] = c;

              exponents.emplace_back(exponent);
            }
        }

  deallog << "  Relative interpolation errors:";
  for (const auto &exponent : exponents)
    {
      const Functions::Monomial<dim> test_function(exponent, fe.n_components());

      hp::MappingCollection<dim> mapping_collection;
      mapping_collection.push_back(
        fe[0].reference_cell().template get_default_linear_mapping<dim>());
      mapping_collection.push_back(
        fe[1].reference_cell().template get_default_linear_mapping<dim>());

      // interpolate the function
      VectorTools::interpolate(mapping_collection,
                               dof_handler,
                               test_function,
                               interpolant);

      // then compute the interpolation error
      hp::QCollection<dim> quadrature_collection;
      for (unsigned int i = 0; i < fe.size(); ++i)
        {
          unsigned int n_quadrature_points;
          if (fe[i].reference_cell().is_hyper_cube())
            n_quadrature_points = fe[i].degree + 3;
          else if (fe[i].reference_cell().is_simplex() ||
                   fe[i].reference_cell() == ReferenceCells::Wedge)
            n_quadrature_points = 4;
          else
            n_quadrature_points = 2;

          quadrature_collection.push_back(
            fe[i].reference_cell().template get_gauss_type_quadrature<dim>(
              n_quadrature_points));
        }

      VectorTools::integrate_difference(mapping_collection,
                                        dof_handler,
                                        interpolant,
                                        test_function,
                                        error,
                                        quadrature_collection,
                                        VectorTools::L2_norm);

      if (error.l2_norm() / interpolant.l2_norm() < 1e-12)
        deallog << " ok";

      Assert(error.l2_norm() < 1e-12 * interpolant.l2_norm(),
             ExcInternalError());
    }
  deallog << std::endl;
}

template <int dim>
void
test_mixed_mesh(const Triangulation<dim> &tria,
                const unsigned int        degree,
                const unsigned int        degree_neighbor,
                const bool                swap_fe_indices)
{
  hp::FECollection<dim> fe;
  const auto            reference_cells = tria.get_reference_cells();

  for (unsigned int i = 0; i < reference_cells.size(); ++i)
    {
      const unsigned int current_degree = (i == 0) ? degree : degree_neighbor;
      const unsigned int index =
        swap_fe_indices ? reference_cells.size() - 1 - i : i;
      const auto reference_cell = reference_cells[index];

      if (reference_cell.is_hyper_cube())
        fe.push_back(FE_Q<dim>(current_degree));
      else if (reference_cell.is_simplex())
        fe.push_back(FE_SimplexP<dim>(current_degree));
      else if (reference_cell == ReferenceCells::Pyramid)
        fe.push_back(FE_PyramidP<dim>(current_degree));
      else if (reference_cell == ReferenceCells::Wedge)
        fe.push_back(FE_WedgeP<dim>(current_degree));
    }

  test_interpolation(fe, tria);
  deallog << std::endl;
}


int
main()
{
  initlog();

  Triangulation<3>         tria;
  std::vector<Point<3>>    vertices;
  std::vector<CellData<3>> cells;
  {
    // tet + pyramid mesh
    tria.clear();
    vertices.clear();
    cells.clear();

    vertices.push_back(Point<3>(0., 0., 0.));
    vertices.push_back(Point<3>(1., 0., 0.));
    vertices.push_back(Point<3>(0., 1., 0.));
    vertices.push_back(Point<3>(1., 1., 0.));
    vertices.push_back(Point<3>(0., 0., 1.));
    vertices.push_back(Point<3>(1., 1., 1.));
    {
      CellData<3> pyramid;
      pyramid.vertices = {0, 1, 2, 3, 4};
      cells.push_back(pyramid);
    }
    {
      CellData<3> tet;
      tet.vertices = {3, 4, 1, 5};
      cells.push_back(tet);
    }
    tria.create_triangulation(vertices, cells, SubCellData());

    for (unsigned int i = 1; i < 4; ++i)
      {
        test_mixed_mesh<3>(tria, i, 1, false);
        test_mixed_mesh<3>(tria, 1, i, true);
      }
  }

  {
    // tet + wedge mesh
    tria.clear();
    vertices.clear();
    cells.clear();

    vertices.push_back(Point<3>(0., 0., 0.));
    vertices.push_back(Point<3>(1., 0., 0.));
    vertices.push_back(Point<3>(0., 1., 0.));
    vertices.push_back(Point<3>(0., 0., 1.));
    vertices.push_back(Point<3>(1., 0., 1.));
    vertices.push_back(Point<3>(0., 1., 1.));
    vertices.push_back(Point<3>(0., 0., 2.));
    {
      CellData<3> wedge;
      wedge.vertices = {0, 1, 2, 3, 4, 5};
      cells.push_back(wedge);
    }
    {
      CellData<3> tet;
      tet.vertices = {3, 4, 5, 6};
      cells.push_back(tet);
    }
    tria.create_triangulation(vertices, cells, SubCellData());
    tria.refine_global(2);

    for (unsigned int i = 1; i < 4; ++i)
      for (unsigned int j = 1; j < 3; ++j)
        {
          test_mixed_mesh<3>(tria, i, j, false);
          test_mixed_mesh<3>(tria, j, i, true);
        }
  }

  {
    // hex + pyramid mesh
    tria.clear();
    vertices.clear();
    cells.clear();

    vertices.push_back(Point<3>(0., 0., 1.));
    vertices.push_back(Point<3>(1., 0., 1.));
    vertices.push_back(Point<3>(0., 1., 1.));
    vertices.push_back(Point<3>(1., 1., 1.));
    vertices.push_back(Point<3>(0., 0., 2.));
    vertices.push_back(Point<3>(0., 0., 0.));
    vertices.push_back(Point<3>(1., 0., 0.));
    vertices.push_back(Point<3>(0., 1., 0.));
    vertices.push_back(Point<3>(1., 1., 0.));
    {
      CellData<3> pyramid;
      pyramid.vertices = {0, 1, 2, 3, 4};
      cells.push_back(pyramid);
    }
    {
      CellData<3> hex;
      hex.vertices = {5, 6, 7, 8, 0, 1, 2, 3};
      cells.push_back(hex);
    }
    tria.create_triangulation(vertices, cells, SubCellData());

    for (unsigned int i = 1; i < 4; ++i)
      {
        test_mixed_mesh<3>(tria, 1, i, false);
        test_mixed_mesh<3>(tria, i, 1, true);
      }
  }

  {
    // hex + wedge mesh
    tria.clear();
    vertices.clear();
    cells.clear();

    vertices.push_back(Point<3>(0., 0., 0.));
    vertices.push_back(Point<3>(1., 0., 0.));
    vertices.push_back(Point<3>(0., 1., 0.));
    vertices.push_back(Point<3>(1., 1., 0.));
    vertices.push_back(Point<3>(0., 0., 1.));
    vertices.push_back(Point<3>(1., 0., 1.));
    vertices.push_back(Point<3>(0., 1., 1.));
    vertices.push_back(Point<3>(1., 1., 1.));
    vertices.push_back(Point<3>(2., 0.5, 0.));
    vertices.push_back(Point<3>(2., 0.5, 1.));
    {
      CellData<3> wedge;
      wedge.vertices = {2, 8, 3, 6, 9, 7};
      cells.push_back(wedge);
    }
    {
      CellData<3> hex;
      hex.vertices = {0, 1, 2, 3, 4, 5, 6, 7};
      cells.push_back(hex);
    }
    tria.create_triangulation(vertices, cells, SubCellData());
    tria.refine_global(2);

    for (unsigned int i = 1; i < 4; ++i)
      for (unsigned int j = 1; j < 3; ++j)
        {
          test_mixed_mesh<3>(tria, j, i, false);
          test_mixed_mesh<3>(tria, i, j, true);
        }
  }

  {
    // pyramid + wedge mesh at quad face
    tria.clear();
    vertices.clear();
    cells.clear();

    vertices.push_back(Point<3>(0., 0., 0.));
    vertices.push_back(Point<3>(1., 0., 0.));
    vertices.push_back(Point<3>(0., 1., 0.));
    vertices.push_back(Point<3>(1., 1., 0.));
    vertices.push_back(Point<3>(0., 0.5, 0.5));
    vertices.push_back(Point<3>(1., 0.5, 0.5));
    vertices.push_back(Point<3>(0., 0., -1.));
    {
      CellData<3> pyramid;
      pyramid.vertices = {2, 3, 0, 1, 6};
      cells.push_back(pyramid);
    }
    {
      CellData<3> wedge;
      wedge.vertices = {0, 2, 4, 1, 3, 5};
      cells.push_back(wedge);
    }
    tria.create_triangulation(vertices, cells, SubCellData());

    for (unsigned int i = 1; i < 3; ++i)
      {
        test_mixed_mesh<3>(tria, 1, i, false);
        test_mixed_mesh<3>(tria, i, 1, true);
      }
  }

  {
    // pyramid + wedge mesh at triangular face
    tria.clear();
    vertices.clear();
    cells.clear();

    vertices.push_back(Point<3>(0., 0., 0.));
    vertices.push_back(Point<3>(1., 0., 0.));
    vertices.push_back(Point<3>(0., 1., 0.));
    vertices.push_back(Point<3>(1., 1., 0.));
    vertices.push_back(Point<3>(0., 0., 1.));
    vertices.push_back(Point<3>(-1., 0., 0.));
    vertices.push_back(Point<3>(-1., 1., 0.));
    vertices.push_back(Point<3>(-1., 0., 1.));
    {
      CellData<3> pyramid;
      pyramid.vertices = {0, 1, 2, 3, 4};
      cells.push_back(pyramid);
    }
    {
      CellData<3> wedge;
      wedge.vertices = {5, 6, 7, 0, 2, 4};
      cells.push_back(wedge);
    }
    tria.create_triangulation(vertices, cells, SubCellData());

    for (unsigned int i = 1; i < 3; ++i)
      {
        test_mixed_mesh<3>(tria, 1, i, false);
        test_mixed_mesh<3>(tria, i, 1, true);
      }
  }

  return 0;
}
