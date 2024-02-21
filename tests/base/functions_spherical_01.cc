// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test value, gradient and hessian of a Functions::Spherical as
// compared to radially symmetric function.

#include <deal.II/base/function_spherical.h>
#include <deal.II/base/geometric_utilities.h>

#include "../tests.h"


// reference function f(x) = exp (-Z*r)
template <int dim>
class ExpFunc : public Function<dim>
{
public:
  ExpFunc(const Point<dim> &origin, const double &Z)
    : Function<dim>(1)
    , origin(origin)
    , Z(Z)
  {}

  virtual double
  value(const Point<dim> &point, const unsigned int component = 0) const
  {
    Tensor<1, dim> dist = point - origin;
    const double   r    = dist.norm();
    return std::exp(-Z * r);
  }

  virtual Tensor<1, dim>
  gradient(const Point<dim> &p, const unsigned int component = 0) const
  {
    Tensor<1, dim> dist = p - origin;
    const double   r    = dist.norm();
    Assert(r > 0.0, ExcMessage("r is not positive"));
    dist /= r;
    return -Z * std::exp(-Z * r) * dist;
  }

  virtual SymmetricTensor<2, dim>
  hessian(const Point<dim> &p, const unsigned int component = 0) const
  {
    Tensor<1, dim> dir = p - origin;
    const double   r   = dir.norm();
    Assert(r > 0.0, ExcMessage("r is not positive"));
    dir /= r;
    SymmetricTensor<2, dim> dir_x_dir;
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = i; j < dim; ++j)
        dir_x_dir[i][j] = dir[i] * dir[j];

    return Z * std::exp(-Z * r) *
           ((Z + 1.0 / r) * dir_x_dir - unit_symmetric_tensor<dim>() / r);
  }

private:
  const Point<dim> origin;
  const double     Z;
};


// same as above but using Functions::Spherical
template <int dim>
class ExpFunc2 : public Functions::Spherical<dim>
{
public:
  ExpFunc2(const Point<dim> &origin, const double &Z)
    : Functions::Spherical<dim>(origin)
    , Z(Z)
  {}

private:
  virtual double
  svalue(const std::array<double, dim> &sp, const unsigned int) const
  {
    return std::exp(-Z * sp[0]);
  }

  virtual std::array<double, dim>
  sgradient(const std::array<double, dim> &sp, const unsigned int) const
  {
    std::array<double, dim> res;
    res[0] = -Z * std::exp(-Z * sp[0]);
    for (unsigned int i = 1; i < dim; ++i)
      res[i] = 0.;
    return res;
  }

  virtual std::array<double, 6>
  shessian(const std::array<double, dim> &sp, const unsigned int) const
  {
    std::array<double, 6> res;
    res[0] = Z * Z * std::exp(-Z * sp[0]);
    for (unsigned int i = 1; i < 6; ++i)
      res[i] = 0.;
    return res;
  }

  const double Z;
};

template <int dim>
void
check()
{
  Point<dim>   center;
  const double Z = 2.5;
  center[1]      = 2.0;
  if (dim > 2)
    center[2] = -1.5;

  ExpFunc<dim>  func(center, Z);
  ExpFunc2<dim> func2(center, Z);

  for (double r = 0.1; r < 10; r += 0.35)
    for (double theta = 0; theta < 2 * numbers::PI; theta += numbers::PI / 3.)
      for (double phi = 0.01; phi <= numbers::PI; phi += numbers::PI / 4.)
        {
          std::array<double, dim> sp;
          sp[0]        = r;
          sp[1]        = theta;
          sp[2]        = phi;
          Point<dim> p = GeometricUtilities::Coordinates::from_spherical(sp);
          for (unsigned int i = 0; i < dim; ++i)
            p[i] += center[i];

          // check values:
          const double v1 = func.value(p);
          const double v2 = func2.value(p);
          AssertThrow(std::fabs(v1 - v2) <= std::abs(v1) * 1e-10,
                      ExcInternalError());

          // check gradients:
          const Tensor<1, dim> g1 = func.gradient(p);
          const Tensor<1, dim> g2 = func2.gradient(p);
          const Tensor<1, dim> gd = g1 - g2;
          AssertThrow(gd.norm() <= g1.norm() * 1e-10, ExcInternalError());

          // check hessian:
          const SymmetricTensor<2, dim> h1 = func.hessian(p);
          const SymmetricTensor<2, dim> h2 = func2.hessian(p);
          const SymmetricTensor<2, dim> dh = h1 - h2;
          AssertThrow(dh.norm() <= h1.norm() * 1e-10, ExcInternalError());
        }
  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();

  check<3>();
}
