// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



// tests matrix-free face evaluation, matrix-vector products as compared to
// the same implementation with MeshWorker. This example uses a mesh described
// through a complicated manifold

#include <deal.II/base/function.h>

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/grid/grid_tools.h>

#include "../tests.h"

#include "create_mesh.h"
#include "matrix_vector_faces_common.h"


double
f_x(double x_m)
{
  double x   = x_m * 1000.0;
  double y_m = 0.0;

  if (x <= 9.0)
    y_m = 0.001 * (-28.0 + std::min(28.0,
                                    2.8e1 + 6.775070969851e-3 * x * x -
                                      2.124527775800e-3 * x * x * x));

  else if (x > 9.0 && x <= 14.0)
    y_m = 0.001 * (-28.0 + 2.507355893131e1 + 9.754803562315e-1 * x -
                   1.016116352781e-1 * x * x + 1.889794677828e-3 * x * x * x);

  else if (x > 14.0 && x <= 20.0)
    y_m = 0.001 * (-28.0 + 2.579601052357e1 + 8.206693007457e-1 * x -
                   9.055370274339e-2 * x * x + 1.626510569859e-3 * x * x * x);

  else if (x > 20.0 && x <= 30.0)
    y_m = 0.001 * (-28.0 + 4.046435022819e1 - 1.379581654948 * x +
                   1.945884504128e-2 * x * x - 2.070318932190e-4 * x * x * x);

  else if (x > 30.0 && x <= 40.0)
    y_m = 0.001 * (-28.0 + 1.792461334664e1 + 8.743920332081e-1 * x -
                   5.567361123058e-2 * x * x + 6.277731764683e-4 * x * x * x);

  else if (x > 40.0 && x <= 54.0)
    y_m = 0.001 * (-28.0 + std::max(0.0,
                                    5.639011190988e1 - 2.010520359035 * x +
                                      1.644919857549e-2 * x * x +
                                      2.674976141766e-5 * x * x * x));

  else if (x > 54.0)
    y_m = 0.001 * (-28.0);

  return y_m;
}

double
gamma_x(double x_m)
{
  double xmod  = x_m / 0.028;
  double gamma = 1.;
  if (xmod > 6.3 && xmod <= 8.3)
    {
      xmod  = 8.3 - xmod;
      gamma = (-0.02 * std::cos(xmod * numbers::PI) + 1.02);
    }
  else if (xmod > 8.3) // move mesh closer to the wall to get a lower y+ value
                       // at the peak
    {
      xmod  = 9. - xmod;
      gamma = (-0.05 * std::cos(xmod * numbers::PI * 2. / 0.7) + 1.05);
    }
  return 1.0 * gamma;
}



template <int dim>
class PeriodicHillManifold : public ChartManifold<dim>
{
public:
  PeriodicHillManifold()
    : ChartManifold<dim>()
    , h(0.028)
    , x_max(9.0 * h)
    , y_max(2.036 * h)
    , y_FoR(h)
  {}

  virtual std::unique_ptr<Manifold<dim>>
  clone() const override
  {
    return std_cxx14::make_unique<PeriodicHillManifold<dim>>();
  }

  virtual Point<dim>
  push_forward(const Point<dim> &xi) const override
  {
    Point<dim> x = xi;

    const double gamma      = gamma_x(xi[0]);
    const double tanh_gamma = std::tanh(gamma);
    double       s_y =
      std::tanh(gamma * (2.0 * ((xi[1] - y_FoR) / y_max) - 1.0)) / tanh_gamma;
    double t_y =
      std::tanh(gamma * (2.0 * (1.0 - (xi[1] - y_FoR) / y_max) - 1.0)) /
      tanh_gamma;
    if (xi[0] <= x_max / 2.0)
      x[1] =
        y_max / 2.0 * s_y + 4.036 * h / 2.0 + (0.5 * t_y + 0.5) * f_x(xi[0]);
    // y_max/2.0*t_y+4.036*h/2.0
    else if (xi[0] > x_max / 2.0)
      x[1] = y_max / 2.0 * s_y + 4.036 * h / 2.0 +
             (0.5 * t_y + 0.5) * f_x(x_max - xi[0]);

    return x;
  }

  virtual Point<dim>
  pull_back(const Point<dim> &x) const override
  {
    Point<dim> xi = x;

    // y component
    unsigned int iter    = 0;
    unsigned int maxiter = 100;
    double       tol     = 1.0e-14;
    double       eps     = 1.;

    // get a good estimate for Y (pullBack without stretching the cells)
    double Y = 0.0;
    if (x[0] <= x_max / 2.0)
      Y = (x[1] - f_x(x[0]) * (1 + y_FoR / y_max)) / (1.0 - f_x(x[0]) / y_max);
    else if (x[0] > x_max / 2.0)
      Y = (x[1] - f_x(x_max - x[0]) * (1 + y_FoR / y_max)) /
          (1.0 - f_x(x_max - x[0]) / y_max);

    if (x[0] <= x_max / 2.0)
      {
        const double gamma      = gamma_x(x[0]);
        const double tanh_gamma = std::tanh(gamma);
        while (eps > tol && iter < maxiter)
          {
            const double arg =
              gamma * (2.0 * (1.0 - (Y - y_FoR) / y_max) - 1.0);
            const double arg2 = gamma * (2.0 * ((Y - y_FoR) / y_max) - 1.0);
            const double t_y  = std::tanh(arg) / tanh_gamma;
            const double s_y  = std::tanh(arg2) / tanh_gamma;
            const double ts_y = 1.0 / (std::cosh(arg) * std::cosh(arg)) *
                                gamma * (-2.0 / y_max) / tanh_gamma;
            const double ss_y = 1.0 / (std::cosh(arg2) * std::cosh(arg2)) *
                                gamma * (2.0 / y_max) / tanh_gamma;
            const double Yn =
              Y - (y_max / 2.0 * s_y + 4.036 * h / 2.0 +
                   (0.5 * t_y + 0.5) * f_x(x[0]) - x[1]) /
                    (y_max / 2.0 * ss_y + 0.5 * ts_y * f_x(x[0]));

            eps = std::abs(Yn - Y);
            Y   = Yn;
            iter++;
          }
        AssertThrow(iter < maxiter,
                    ExcMessage(
                      "Newton within PullBack did not find the solution. "));
        xi[1] = Y;
      }
    else if (x[0] > x_max / 2.0)
      {
        const double gamma      = gamma_x(x[0]);
        const double tanh_gamma = std::tanh(gamma);
        while (eps > tol && iter < maxiter)
          {
            const double arg =
              gamma * (2.0 * (1.0 - (Y - y_FoR) / y_max) - 1.0);
            const double arg2 = gamma * (2.0 * ((Y - y_FoR) / y_max) - 1.0);
            const double t_y  = std::tanh(arg) / tanh_gamma;
            const double s_y  = std::tanh(arg2) / tanh_gamma;
            const double ts_y = 1.0 / (std::cosh(arg) * std::cosh(arg)) *
                                gamma * (-2.0 / y_max) / tanh_gamma;
            const double ss_y = 1.0 / (std::cosh(arg2) * std::cosh(arg2)) *
                                gamma * (2.0 / y_max) / tanh_gamma;
            const double Yn =
              Y - (y_max / 2.0 * s_y + 4.036 * h / 2.0 +
                   (0.5 * t_y + 0.5) * f_x(x_max - x[0]) - x[1]) /
                    (y_max / 2.0 * ss_y + 0.5 * ts_y * f_x(x_max - x[0]));

            eps = std::abs(Yn - Y);
            Y   = Yn;
            iter++;
          }
        AssertThrow(iter < maxiter,
                    ExcMessage(
                      "Newton within PullBack did not find the solution. "));
        xi[1] = Y;
      }
    return xi;
  }

private:
  const double h = 0.028;

  // data from initial block
  const double x_max = 9.0 * h; // 9.0*h;
  const double y_max = 2.036 * h;

  const double y_FoR = h;
};


template <int dim, int fe_degree>
void
test()
{
  Triangulation<dim> triangulation;
  /* --------------- Generate grid ------------------- */
  const double h = 0.028;
  Point<dim>   coordinates;
  coordinates[0] = 9.0 * h;   // 9.0*h;
  coordinates[1] = 3.036 * h; // 2.036*h;
  if (dim == 3)
    coordinates[2] = 2.25 * h; // 4.5*h;
  std::vector<unsigned int> refinements(dim, 1);
  refinements[0] = 2;

  // start with a cube
  Point<dim> p;
  p[0] = 0.;
  p[1] = h;
  if (dim == 3)
    p[2] = -2.25 * h;

  GridGenerator::subdivided_hyper_rectangle(triangulation,
                                            refinements,
                                            p,
                                            coordinates);

  triangulation.last()->vertex(0)[1] = 0.;
  if (dim == 3)
    triangulation.last()->vertex(4)[1] = 0.;
  // boundary ids for refinements[0] = 2:
  // periodicity in x-direction
  // add 10 to avoid conflicts with dirichlet boundary, which is 0
  triangulation.begin()->face(0)->set_all_boundary_ids(0 + 10);
  triangulation.last()->face(1)->set_all_boundary_ids(1 + 10);

  // periodicity in z-direction, if dim==3
  if (dim == 3)
    {
      triangulation.begin()->face(4)->set_all_boundary_ids(2 + 10);
      triangulation.begin()->face(5)->set_all_boundary_ids(3 + 10);
      triangulation.last()->face(4)->set_all_boundary_ids(2 + 10);
      triangulation.last()->face(5)->set_all_boundary_ids(3 + 10);
    }

  std::vector<
    GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator>>
    periodic_faces;
  GridTools::collect_periodic_faces(
    triangulation, 0 + 10, 1 + 10, 0, periodic_faces);
  if (dim == 3)
    {
      GridTools::collect_periodic_faces(
        triangulation, 2 + 10, 3 + 10, 2, periodic_faces);
    }

  triangulation.add_periodicity(periodic_faces);
  triangulation.set_all_manifold_ids(111);

  static PeriodicHillManifold<dim> manifold;
  triangulation.set_manifold(111, manifold);
  triangulation.refine_global(2);

  FE_DGQ<dim>     fe(fe_degree);
  DoFHandler<dim> dof(triangulation);
  dof.distribute_dofs(fe);
  AffineConstraints<double> constraints;
  constraints.close();

  // cannot test threads right now
  do_test<dim, fe_degree, fe_degree + 1, double>(dof, constraints, false);
}
