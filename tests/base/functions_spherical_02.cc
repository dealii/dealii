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


// test value, gradient and Hessian of a Functions::Spherical
// for the case of Px orbital
// (sqrt(3/4pi) sin(phi) cos(theta) == sqrt(3/pi)*x/r)
// multiplied with r*r.

/* MWE in Maxima:

spherical:
f: r*r*cos(theta)*sin(phi);
[diff(f,r), diff(f,theta), diff(f,phi)];
[diff(f,r,2), diff(f,theta,2), diff(f,phi,2), diff(f,r,1,theta,1),
diff(f,r,1,phi,1), diff (f,theta,1,phi,1)];

Cartesian:

assume(r>0);
srule: [r= sqrt(x*x+y*y+z*z)];
srule2: [sqrt(x*x+y*y+z*z) = r, x*x+y*y+z*z = r^2, expand((x*x+y*y+z*z)^2) =
r^4]; f: subst(srule,r*r*x/r); df: fullratsimp([diff(f,x), diff(f,y),
diff(f,z)]); df2 :(subst(srule2,expand(df))); ddf : fullratsimp([diff(f,x,2),
diff(f,y,2), diff(f,z,2), diff(f,x,1,y,1), diff(f,x,1,z,1), diff(f,y,1,z,1)]);
ddf2: (subst(srule2,expand(ddf)));

 */

#include <deal.II/base/function_spherical.h>
#include <deal.II/base/geometric_utilities.h>

#include "../tests.h"



template <int dim>
class RefFunc : public Function<dim>
{
public:
  RefFunc(const Point<dim> origin = Point<dim>())
    : Function<dim>(1)
    , origin(origin)
  {}

  virtual double
  value(const Point<dim> &point, const unsigned int component = 0) const
  {
    Tensor<1, dim> dist = point - origin;
    const double   r    = dist.norm();
    return r * dist[0];
  }

  virtual Tensor<1, dim>
  gradient(const Point<dim> &p, const unsigned int component = 0) const
  {
    Tensor<1, dim> dist = p - origin;
    const double   r    = dist.norm();
    Assert(r > 0.0, ExcMessage("r is not positive"));
    Tensor<1, dim> res;
    const double   x  = dist[0];
    const double   y  = dist[1];
    const double   z  = dist[2];
    const double   x2 = x * x;
    const double   y2 = y * y;
    const double   z2 = z * z;

    res[0] = z2 / r + y2 / r + (2. * x2) / r;
    res[1] = (x * y) / r;
    res[2] = (x * z) / r;

    return res;
  }

  virtual SymmetricTensor<2, dim>
  hessian(const Point<dim> &p, const unsigned int component = 0) const
  {
    const Tensor<1, dim> dist = p - origin;
    const double         r    = dist.norm();
    Assert(r > 0.0, ExcMessage("r is not positive"));
    const double x  = dist[0];
    const double y  = dist[1];
    const double z  = dist[2];
    const double z2 = z * z;
    const double y2 = y * y;
    const double x2 = x * x;
    const double x3 = x2 * x;
    const double z3 = z2 * z;
    const double y3 = y2 * y;
    const double r3 = r * r * r;


    SymmetricTensor<2, dim> res;
    res[0][0] = (3. * x * z2) / r3 + (3 * x * y2) / r3 + (2. * x3) / r3;
    res[1][1] = (x * z2) / r3 + x3 / r3;
    res[2][2] = (x * y2) / r3 + x3 / r3;
    res[0][1] = (y * z2) / r3 + y3 / r3;
    res[0][2] = z3 / r3 + (y2 * z) / r3;
    res[1][2] = -(x * y * z) / r3;

    return res;
  }

private:
  const Point<dim> origin;
};


// same as above but using Functions::Spherical
template <int dim>
class SphFunc : public Functions::Spherical<dim>
{
public:
  SphFunc(const Point<dim> origin = Point<dim>())
    : Functions::Spherical<dim>(origin)
  {}

private:
  virtual double
  svalue(const std::array<double, dim> &sp, const unsigned int) const
  {
    return sp[0] * sp[0] * std::cos(sp[1]) * std::sin(sp[2]);
  }

  virtual std::array<double, dim>
  sgradient(const std::array<double, dim> &sp, const unsigned int) const
  {
    std::array<double, dim> res;
    const double            r     = sp[0];
    const double            theta = sp[1];
    const double            phi   = sp[2];
    res[0]                        = 2. * sin(phi) * r * cos(theta);
    res[1]                        = -sin(phi) * r * r * sin(theta);
    res[2]                        = cos(phi) * r * r * cos(theta);
    return res;
  }

  virtual std::array<double, 6>
  shessian(const std::array<double, dim> &sp, const unsigned int) const
  {
    std::array<double, 6> res;
    const double          r     = sp[0];
    const double          theta = sp[1];
    const double          phi   = sp[2];
    const double          r2    = r * r;
    res[0]                      = 2. * sin(phi) * cos(theta);
    res[1]                      = -sin(phi) * r2 * cos(theta);
    res[2]                      = -sin(phi) * r2 * cos(theta);
    res[3]                      = -2. * sin(phi) * r * sin(theta);
    res[4]                      = 2. * cos(phi) * r * cos(theta);
    res[5]                      = -cos(phi) * r2 * sin(theta);
    return res;
  }
};

template <int dim>
void
check()
{
  Point<dim> center;
  center[1] = 2.0;
  if (dim > 2)
    center[2] = -1.5;

  RefFunc<dim> func(center);
  SphFunc<dim> func2(center);

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
