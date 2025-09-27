// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1999 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check accuracy of various quadrature formulas by using them to
// integrate polynomials of increasing degree, and finding the degree
// until which they integrate exactly


#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/grid/reference_cell.h>

#include "../tests.h"

template <int dim>
void
fill_vector(std::vector<Quadrature<dim> *> &quadratures)
{
  quadratures.push_back(new QMidpoint<dim>());
  quadratures.push_back(new QTrapezoid<dim>());
  quadratures.push_back(new QSimpson<dim>());
  quadratures.push_back(new QMilne<dim>());
  quadratures.push_back(new QWeddle<dim>());
  for (unsigned int i = 0; i < 9; ++i)
    {
      quadratures.push_back(new QGauss<dim>(i));
    }
  QMilne<1> q1d;
  quadratures.push_back(new Quadrature<dim>(q1d));
  for (unsigned int i = 2; i < 8; ++i)
    {
      quadratures.push_back(new QGaussLobatto<dim>(i));
    }
  for (unsigned int i = 1; i < 9; ++i)
    {
      quadratures.push_back(new QGaussRadau<dim>(i));
      quadratures.push_back(
        new QGaussRadau<dim>(i, QGaussRadau<dim>::EndPoint::right));
    }
}

template <int dim>
void
check_cells(std::vector<Quadrature<dim> *> &quadratures)
{
  Quadrature<dim> quadrature;
  for (unsigned int n = 0; n < quadratures.size(); ++n)
    {
      quadrature                             = *quadratures[n];
      const std::vector<Point<dim>> &points  = quadrature.get_points();
      const std::vector<double>     &weights = quadrature.get_weights();

      deallog << "Quadrature no." << n;

      unsigned int i              = 0;
      double       quadrature_int = 0;
      double       exact_int      = 0;
      double       err            = 0;

      do
        {
          ++i;

          quadrature_int = 0;
          // Check the polynomial x^i*y^i

          for (unsigned int x = 0; x < quadrature.size(); ++x)
            {
              double f = 1.;
              switch (dim)
                {
                  case 3:
                    f *= std::pow(static_cast<double>(points[x][2]), i * 1.0);
                  case 2:
                    f *= std::pow(static_cast<double>(points[x][1]), i * 1.0);
                  case 1:
                    f *= std::pow(static_cast<double>(points[x][0]), i * 1.0);
                }
              quadrature_int += f * weights[x];
            }

          // the exact integral is 1/(i+1)
          exact_int = 1. / std::pow(static_cast<double>(i + 1), dim);
          err       = std::fabs(quadrature_int - exact_int);
        }
      while (err < 1e-14);
      // Uncomment here for testing
      //      deallog << " (Int " << quadrature_int << ',' << exact_int << ')';
      deallog << " is exact for polynomials of degree " << i - 1 << std::endl;

      if (dim == 1)
        {
          // check the ordering of
          // the quadrature points
          bool in_order = true;
          for (unsigned int x = 1; x < quadrature.size(); ++x)
            {
              if (points[x][0] <= points[x - 1][0])
                in_order = false;
            }
          if (!in_order)
            for (unsigned int x = 0; x < quadrature.size(); ++x)
              deallog << points[x] << std::endl;
        }
    }
}


template <int dim>
void
check_faces(const std::vector<Quadrature<dim - 1> *> &quadratures,
            const bool                                sub)
{
  if (sub)
    deallog.push("subfaces");
  else
    deallog.push("faces");

  for (unsigned int n = 0; n < quadratures.size(); ++n)
    {
      Quadrature<dim> quadrature(
        sub == false ?
          QProjector<dim>::project_to_all_faces(
            ReferenceCells::get_hypercube<dim>(), *quadratures[n]) :
          QProjector<dim>::project_to_all_subfaces(
            ReferenceCells::get_hypercube<dim>(), *quadratures[n]));
      const std::vector<Point<dim>> &points  = quadrature.get_points();
      const std::vector<double>     &weights = quadrature.get_weights();

      deallog << "Quadrature no." << n;

      unsigned int i              = 0;
      long double  quadrature_int = 0;
      double       exact_int      = 0;
      double       err            = 0;

      do
        {
          ++i;

          quadrature_int = 0;
          // Check the polynomial
          // x^i*y^i*z^i

          for (unsigned int x = 0; x < quadrature.size(); ++x)
            {
              long double f = 1.;
              switch (dim)
                {
                  case 3:
                    f *= std::pow((long double)points[x][2], i * 1.0L);
                  case 2:
                    f *= std::pow((long double)points[x][1], i * 1.0L);
                  case 1:
                    f *= std::pow((long double)points[x][0], i * 1.0L);
                }
              quadrature_int += f * weights[x];
            }

          // the exact integral is
          // 1/(i+1)^(dim-1)
          switch (dim)
            {
              case 2:
                exact_int = 2 * (sub ? 2 : 1) * 2 / (double)(i + 1);
                break;
              case 3:
                exact_int =
                  3 * (sub ? (4 + 2 + 2) : 1) * 8 / (double)(i + 1) / (i + 1);
                break;
            }

          err = std::fabs(quadrature_int - exact_int);
        }
      // for comparison: use factor 8 in case
      // of dim==3, as we integrate 8 times
      // over the whole surface (all
      // combinations of face_orientation,
      // face_flip and face_rotation)
      // but 2 in 2d (only two possible
      // orientations) when we are not using
      // subfaces (same as the previous comment)
      while (err < (dim == 3 ? 8 : (sub ? 1 : 2)) * 2e-14);
      // Uncomment here for testing
      //      deallog << " (Int " << quadrature_int << '-' << exact_int << '='
      //      << err << ')';
      deallog << " is exact for polynomials of degree " << i - 1 << std::endl;
    }
  deallog.pop();
}

int
main()
{
  initlog();

  std::vector<Quadrature<1> *> q1;
  std::vector<Quadrature<2> *> q2;
  std::vector<Quadrature<3> *> q3;
  fill_vector(q1);
  fill_vector(q2);
  fill_vector(q3);

  deallog.push("1d");
  check_cells(q1);
  deallog.pop();

  deallog.push("2d");
  check_cells(q2);
  check_faces<2>(q1, false);
  check_faces<2>(q1, true);
  deallog.pop();

  deallog.push("3d");
  check_cells(q3);
  check_faces<3>(q2, false);
  check_faces<3>(q2, true);
  deallog.pop();

  // delete objects again to avoid
  // messages about memory leaks when
  // using purify or other memory
  // checkers
  for (unsigned int i = 0; i < q1.size(); ++i)
    delete q1[i];
  for (unsigned int i = 0; i < q2.size(); ++i)
    delete q2[i];
  for (unsigned int i = 0; i < q3.size(); ++i)
    delete q3[i];
}
