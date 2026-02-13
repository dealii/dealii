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



// check that computation of hp-constraints works for FE_SimplexP, FE_PyramidP
// and FE_WedgeP elements by interpolating a function on a pure mesh
// see also the tests in tests/hp/hp_constraints_common.h


#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_pyramid_p.h>
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
  Assert(triangulation.get_reference_cells().size() == 1, ExcInternalError());

  deallog << "Testing interpolation " << fe[0].get_name() << " vs. "
          << fe[1].get_name() << std::endl;


  DoFHandler<dim>    dof_handler(triangulation);
  const unsigned int n_cells      = triangulation.n_active_cells();
  unsigned int       cell_counter = 0;
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell_counter < n_cells / 2)
        {
          cell->set_active_fe_index(0);
        }
      else
        {
          cell->set_active_fe_index(1);
        }
      ++cell_counter;
    }

  dof_handler.distribute_dofs(fe);

  Vector<double> interpolant(dof_handler.n_dofs());
  Vector<float>  error(triangulation.n_active_cells());

  const unsigned int min_degree = std::min(fe[0].degree, fe[1].degree);
  // interpolate monomials of the for x1^a*x2^b*x3^c
  // tets can interpolate a + b + c <= degree exactly
  // wedges can interpolate a + b <= degree and c <= degree exactly
  // pyramids can interpolate a + b + c <= degree exactly
  // so now figure out what to use
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
          if (fe[i].reference_cell().is_simplex() ||
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
test_mesh(const Triangulation<dim> &tria,
          const unsigned int        degree1,
          const unsigned int        degree2)
{
  hp::FECollection<dim> fe;

  const auto reference_cells = tria.get_reference_cells();

  if (reference_cells[0].is_simplex())
    {
      fe.push_back(FE_SimplexP<dim>(degree1));
      fe.push_back(FE_SimplexP<dim>(degree2));
    }
  else if (reference_cells[0] == ReferenceCells::Pyramid)
    {
      fe.push_back(FE_PyramidP<dim>(degree1));
      fe.push_back(FE_PyramidP<dim>(degree2));
    }
  else if (reference_cells[0] == ReferenceCells::Wedge)
    {
      fe.push_back(FE_WedgeP<dim>(degree1));
      fe.push_back(FE_WedgeP<dim>(degree2));
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
    // tet mesh
    tria.clear();
    GridGenerator::subdivided_hyper_cube_with_simplices(tria, 5);

    for (unsigned int i = 1; i < 4; ++i)
      for (unsigned int j = 1; j < 4; ++j)
        test_mesh<3>(tria, i, j);
  }

  {
    // pyramid mesh
    tria.clear();
    GridGenerator::subdivided_hyper_cube_with_pyramids(tria, 5);

    for (unsigned int i = 1; i < 2; ++i)
      for (unsigned int j = 1; j < 2; ++j)
        test_mesh<3>(tria, i, j);
  }

  {
    // wedge mesh
    tria.clear();
    vertices.clear();
    cells.clear();

    vertices.push_back(Point<3>(0., 0., 0.));
    vertices.push_back(Point<3>(1., 0., 0.));
    vertices.push_back(Point<3>(0., 1., 0.));
    vertices.push_back(Point<3>(0., 0., 1.));
    vertices.push_back(Point<3>(1., 0., 1.));
    vertices.push_back(Point<3>(0., 1., 1.));
    vertices.push_back(Point<3>(1., 1., 0.));
    vertices.push_back(Point<3>(1., 1., 1.));
    vertices.push_back(Point<3>(0., 0., 2.));
    vertices.push_back(Point<3>(1., 0., 2.));
    vertices.push_back(Point<3>(0., 1., 2.));
    vertices.push_back(Point<3>(1., 1., 2.));
    {
      CellData<3> wedge;
      wedge.vertices = {0, 1, 2, 3, 4, 5};
      cells.push_back(wedge);
    }
    {
      CellData<3> wedge;
      wedge.vertices = {1, 6, 2, 4, 7, 5};
      cells.push_back(wedge);
    }
    {
      CellData<3> wedge;
      wedge.vertices = {3, 4, 5, 8, 9, 10};
      cells.push_back(wedge);
    }
    {
      CellData<3> wedge;
      wedge.vertices = {4, 7, 5, 9, 11, 10};
      cells.push_back(wedge);
    }
    tria.create_triangulation(vertices, cells, SubCellData());
    tria.refine_global(2);

    for (unsigned int j = 1; j < 3; ++j)
      for (unsigned int i = 1; i < 3; ++i)
        test_mesh<3>(tria, i, j);
  }

  return 0;
}
