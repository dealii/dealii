// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// evaluate polynomials_bernardi_raugel on the reference cell
// watch for the first dim*vertices_per_cell shape functions to
// yield delta functions when evaluated on the nodes
// watch for the remaining faces_per_cell shape functions to
// have magnitude 1 on the face centers

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/polynomials_bernardi_raugel.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor.h>

#include <vector>

#include "../tests.h"

template <int dim>
void
check_poly_q(const PolynomialsBernardiRaugel<dim> &poly)
{
  std::vector<Point<dim>>     points;
  std::vector<Tensor<1, dim>> values;
  std::vector<Tensor<2, dim>> grads;
  std::vector<Tensor<3, dim>> grads2;
  std::vector<Tensor<4, dim>> thirds;
  std::vector<Tensor<5, dim>> fourths;

  // set up reference cell vertices
  if (dim == 2)
    {
      for (unsigned int i = 0; i < 2; ++i)
        {
          for (unsigned int j = 0; j < 2; ++j)
            {
              Point<dim> p;
              p[0] = j;
              p[1] = i;
              points.push_back(p);
            }
        }
    }
  else if (dim == 3)
    {
      for (unsigned int i = 0; i < 2; ++i)
        {
          for (unsigned int j = 0; j < 2; ++j)
            {
              for (unsigned int k = 0; k < 2; ++k)
                {
                  Point<dim> p;
                  p[0] = k;
                  p[1] = j;
                  p[2] = i;
                  points.push_back(p);
                }
            }
        }
    }

  // loop over vertices
  for (unsigned int i = 0; i < points.size(); ++i)
    {
      values.clear();
      grads.clear();
      grads2.clear();
      thirds.clear();
      fourths.clear();
      values.resize(dim * GeometryInfo<dim>::vertices_per_cell +
                    GeometryInfo<dim>::faces_per_cell);

      deallog << "BR1<" << dim << "> point " << i << " (" << points[i][0];
      for (unsigned int d = 1; d < dim; ++d)
        deallog << ", " << points[i][d];
      deallog << ')' << std::endl;

      poly.evaluate(points[i], values, grads, grads2, thirds, fourths);

      // loop through the Q_1^d shape functions
      for (unsigned int j = 0; j < dim * GeometryInfo<dim>::vertices_per_cell;
           ++j)
        {
          deallog << "BR1<" << dim << "> shape fxn " << j << ": ";
          for (unsigned int d = 0; d < dim; ++d)
            deallog << '\t' << values[j][d];
          deallog << std::endl;
        }
    }
}



template <int dim>
void
check_poly_bubble(const PolynomialsBernardiRaugel<dim> &poly)
{
  std::vector<Point<dim>>     points;
  std::vector<Tensor<1, dim>> values;
  std::vector<Tensor<2, dim>> grads;
  std::vector<Tensor<3, dim>> grads2;
  std::vector<Tensor<4, dim>> thirds;
  std::vector<Tensor<5, dim>> fourths;

  // set up reference cell face midpoints
  for (unsigned int i = 0; i < dim; ++i)
    {
      for (unsigned int j = 0; j < 2; ++j)
        {
          Point<dim> p;
          p[0] = 0.5;
          p[1] = 0.5;
          if (dim == 3)
            p[2] = 0.5;
          p[i] = j;
          points.push_back(p);
        }
    }

  // loop over vertices
  for (unsigned int i = 0; i < points.size(); ++i)
    {
      values.clear();
      grads.clear();
      grads2.clear();
      thirds.clear();
      fourths.clear();
      values.resize(dim * GeometryInfo<dim>::vertices_per_cell +
                    GeometryInfo<dim>::faces_per_cell);

      deallog << "BR1<" << dim << "> point "
              << (i + GeometryInfo<dim>::vertices_per_cell) << " ("
              << points[i][0];
      for (unsigned int d = 1; d < dim; ++d)
        deallog << ", " << points[i][d];
      deallog << ')' << std::endl;

      poly.evaluate(points[i], values, grads, grads2, thirds, fourths);

      // loop through the bubble shape functions
      for (unsigned int j = dim * GeometryInfo<dim>::vertices_per_cell;
           j < dim * GeometryInfo<dim>::vertices_per_cell +
                 GeometryInfo<dim>::faces_per_cell;
           ++j)
        {
          deallog << "BR1<" << dim << "> shape fxn " << j << ": ";
          for (unsigned int d = 0; d < dim; ++d)
            deallog << '\t' << values[j][d];
          deallog << std::endl;
        }
    }
}

int
main()
{
  initlog();
  deallog << std::setprecision(3);

  PolynomialsBernardiRaugel<2> p_2d(1);
  check_poly_q(p_2d);
  check_poly_bubble(p_2d);

  PolynomialsBernardiRaugel<3> p_3d(1);
  check_poly_q(p_3d);
  check_poly_bubble(p_3d);
}
