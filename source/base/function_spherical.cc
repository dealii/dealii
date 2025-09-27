// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/function_spherical.h>
#include <deal.II/base/geometric_utilities.h>
#include <deal.II/base/point.h>

#include <algorithm>
#include <cmath>

DEAL_II_NAMESPACE_OPEN
namespace Functions
{
  // other implementations/notes:
  // https://github.com/apache/commons-math/blob/master/src/main/java/org/apache/commons/math4/geometry/euclidean/threed/SphericalCoordinates.java
  // http://mathworld.wolfram.com/SphericalCoordinates.html

  /*derivation of Hessian in Maxima as function of tensor products of unit
  vectors:

  depends(ur,[theta,phi]);
  depends(utheta,theta);
  depends(uphi,[theta,phi]);
  depends(f,[r,theta,phi]);
  declare([f,r,theta,phi], scalar)$
  dotscrules: true;
  grads(a):=ur.diff(a,r)+(1/r)*uphi.diff(a,phi)+(1/(r*sin(phi)))*utheta.diff(a,theta);


  H : factor(grads(grads(f)));
  H2: subst([diff(ur,theta)=sin(phi)*utheta,
       diff(utheta,theta)=-cos(phi)*uphi-sin(phi)*ur,
       diff(uphi,theta)=cos(phi)*utheta,
       diff(ur,phi)=uphi,
       diff(uphi,phi)=-ur],
      H);
  H3: trigsimp(fullratsimp(H2));


  srules : [diff(f,r)=sg0,
          diff(f,theta)=sg1,
          diff(f,phi)=sg2,
          diff(f,r,2)=sh0,
          diff(f,theta,2)=sh1,
          diff(f,phi,2)=sh2,
          diff(f,r,1,theta,1)=sh3,
          diff(f,r,1,phi,1)=sh4,
          diff(f,theta,1,phi,1)=sh5,
          cos(phi)=cos_phi,
          cos(theta)=cos_theta,
          sin(phi)=sin_phi,
          sin(theta)=sin_theta
        ]$

  c_utheta2 : distrib(subst(srules, ratcoeff(expand(H3), utheta.utheta)));
  c_utheta_ur : (subst(srules, ratcoeff(expand(H3), utheta.ur)));
  (subst(srules, ratcoeff(expand(H3), ur.utheta))) - c_utheta_ur;
  c_utheta_uphi : (subst(srules, ratcoeff(expand(H3), utheta.uphi)));
  (subst(srules, ratcoeff(expand(H3), uphi.utheta))) - c_utheta_uphi;
  c_ur2 : (subst(srules, ratcoeff(expand(H3), ur.ur)));
  c_ur_uphi : (subst(srules, ratcoeff(expand(H3), ur.uphi)));
  (subst(srules, ratcoeff(expand(H3), uphi.ur))) - c_ur_uphi;
  c_uphi2 : (subst(srules, ratcoeff(expand(H3), uphi.uphi)));


  where (used later to do tensor products):

  ur     : [cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi)];
  utheta : [-sin(theta), cos(theta), 0];
  uphi   : [cos(theta)*cos(phi), sin(theta)*cos(phi), -sin(phi)];

  with the following proof of substitution rules above:

  -diff(ur,theta)+sin(phi)*utheta;
  trigsimp(-diff(utheta,theta)-cos(phi)*uphi-sin(phi)*ur);
  -diff(uphi,theta)+cos(phi)*utheta;
  -diff(ur,phi)+uphi;
  -diff(uphi,phi)-ur;

   */

  namespace
  {
    /**
     * Evaluate unit vectors in spherical coordinates
     */
    template <int dim>
    void
    set_unit_vectors(const double    cos_theta,
                     const double    sin_theta,
                     const double    cos_phi,
                     const double    sin_phi,
                     Tensor<1, dim> &unit_r,
                     Tensor<1, dim> &unit_theta,
                     Tensor<1, dim> &unit_phi)
    {
      unit_r[0] = cos_theta * sin_phi;
      unit_r[1] = sin_theta * sin_phi;
      unit_r[2] = cos_phi;

      unit_theta[0] = -sin_theta;
      unit_theta[1] = cos_theta;
      unit_theta[2] = 0.;

      unit_phi[0] = cos_theta * cos_phi;
      unit_phi[1] = sin_theta * cos_phi;
      unit_phi[2] = -sin_phi;
    }


    /**
     * calculates out[i][j] += v*(in1[i]*in2[j]+in1[j]*in2[i])
     */
    template <int dim>
    void
    add_outer_product(SymmetricTensor<2, dim> &out,
                      const double             val,
                      const Tensor<1, dim>    &in1,
                      const Tensor<1, dim>    &in2)
    {
      if (val != 0.)
        for (unsigned int i = 0; i < dim; ++i)
          for (unsigned int j = i; j < dim; ++j)
            out[i][j] += (in1[i] * in2[j] + in1[j] * in2[i]) * val;
    }

    /**
     * calculates out[i][j] += v*in[i]in[j]
     */
    template <int dim>
    void
    add_outer_product(SymmetricTensor<2, dim> &out,
                      const double             val,
                      const Tensor<1, dim>    &in)
    {
      if (val != 0.)
        for (unsigned int i = 0; i < dim; ++i)
          for (unsigned int j = i; j < dim; ++j)
            out[i][j] += val * in[i] * in[j];
    }
  } // namespace



  template <int dim>
  Spherical<dim>::Spherical(const Point<dim>  &p,
                            const unsigned int n_components)
    : Function<dim>(n_components)
    , coordinate_system_offset(p)
  {
    AssertThrow(dim == 3, ExcNotImplemented());
  }



  template <int dim>
  double
  Spherical<dim>::value(const Point<dim>  &p_,
                        const unsigned int component) const
  {
    const Point<dim>              p = p_ - coordinate_system_offset;
    const std::array<double, dim> sp =
      GeometricUtilities::Coordinates::to_spherical(p);
    return svalue(sp, component);
  }



  template <int dim>
  Tensor<1, dim>
  Spherical<dim>::gradient(const Point<dim> & /*p_*/,
                           const unsigned int /*component*/) const

  {
    DEAL_II_NOT_IMPLEMENTED();
    return {};
  }



  template <>
  Tensor<1, 3>
  Spherical<3>::gradient(const Point<3> &p_, const unsigned int component) const
  {
    constexpr int                 dim = 3;
    const Point<dim>              p   = p_ - coordinate_system_offset;
    const std::array<double, dim> sp =
      GeometricUtilities::Coordinates::to_spherical(p);
    const std::array<double, dim> sg = sgradient(sp, component);

    // somewhat backwards, but we need cos/sin's for unit vectors
    const double cos_theta = std::cos(sp[1]);
    const double sin_theta = std::sin(sp[1]);
    const double cos_phi   = std::cos(sp[2]);
    const double sin_phi   = std::sin(sp[2]);

    Tensor<1, dim> unit_r, unit_theta, unit_phi;
    set_unit_vectors(
      cos_theta, sin_theta, cos_phi, sin_phi, unit_r, unit_theta, unit_phi);

    Tensor<1, dim> res;

    if (sg[0] != 0.)
      {
        res += unit_r * sg[0];
      }

    if (sg[1] * sin_phi != 0.)
      {
        Assert(sp[0] != 0., ExcDivideByZero());
        res += unit_theta * sg[1] / (sp[0] * sin_phi);
      }

    if (sg[2] != 0.)
      {
        Assert(sp[0] != 0., ExcDivideByZero());
        res += unit_phi * sg[2] / sp[0];
      }

    return res;
  }



  template <int dim>
  SymmetricTensor<2, dim>
  Spherical<dim>::hessian(const Point<dim> & /*p*/,
                          const unsigned int /*component*/) const
  {
    DEAL_II_NOT_IMPLEMENTED();
    return {};
  }



  template <>
  SymmetricTensor<2, 3>
  Spherical<3>::hessian(const Point<3> &p_, const unsigned int component) const

  {
    constexpr int                 dim = 3;
    const Point<dim>              p   = p_ - coordinate_system_offset;
    const std::array<double, dim> sp =
      GeometricUtilities::Coordinates::to_spherical(p);
    const std::array<double, dim> sg = sgradient(sp, component);
    const std::array<double, 6>   sh = shessian(sp, component);

    // somewhat backwards, but we need cos/sin's for unit vectors
    const double cos_theta = std::cos(sp[1]);
    const double sin_theta = std::sin(sp[1]);
    const double cos_phi   = std::cos(sp[2]);
    const double sin_phi   = std::sin(sp[2]);
    const double r         = sp[0];

    Tensor<1, dim> unit_r, unit_theta, unit_phi;
    set_unit_vectors(
      cos_theta, sin_theta, cos_phi, sin_phi, unit_r, unit_theta, unit_phi);

    const double sin_phi2 = sin_phi * sin_phi;
    const double r2       = r * r;
    Assert(r != 0., ExcDivideByZero());

    const double c_utheta2 =
      sg[0] / r + ((sin_phi != 0.) ? (cos_phi * sg[2]) / (r2 * sin_phi) +
                                       sh[1] / (r2 * sin_phi2) :
                                     0.);
    const double c_utheta_ur =
      ((sin_phi != 0.) ? (r * sh[3] - sg[1]) / (r2 * sin_phi) : 0.);
    const double c_utheta_uphi =
      ((sin_phi != 0.) ? (sh[5] * sin_phi - cos_phi * sg[1]) / (r2 * sin_phi2) :
                         0.);
    const double c_ur2     = sh[0];
    const double c_ur_uphi = (r * sh[4] - sg[2]) / r2;
    const double c_uphi2   = (sh[2] + r * sg[0]) / r2;

    // go through each tensor product
    SymmetricTensor<2, dim> res;

    add_outer_product(res, c_utheta2, unit_theta);

    add_outer_product(res, c_utheta_ur, unit_theta, unit_r);

    add_outer_product(res, c_utheta_uphi, unit_theta, unit_phi);

    add_outer_product(res, c_ur2, unit_r);

    add_outer_product(res, c_ur_uphi, unit_r, unit_phi);

    add_outer_product(res, c_uphi2, unit_phi);

    return res;
  }



  template <int dim>
  std::size_t
  Spherical<dim>::memory_consumption() const
  {
    return sizeof(Spherical<dim>);
  }



  template <int dim>
  double
  Spherical<dim>::svalue(const std::array<double, dim> & /* sp */,
                         const unsigned int /*component*/) const
  {
    AssertThrow(false, ExcNotImplemented());
    return 0.;
  }



  template <int dim>
  std::array<double, dim>
  Spherical<dim>::sgradient(const std::array<double, dim> & /* sp */,
                            const unsigned int /*component*/) const
  {
    AssertThrow(false, ExcNotImplemented());
    return std::array<double, dim>();
  }



  template <int dim>
  std::array<double, 6>
  Spherical<dim>::shessian(const std::array<double, dim> & /* sp */,
                           const unsigned int /*component*/) const
  {
    AssertThrow(false, ExcNotImplemented());
    return std::array<double, 6>();
  }



  // explicit instantiations
  template class Spherical<1>;
  template class Spherical<2>;
  template class Spherical<3>;

} // namespace Functions

DEAL_II_NAMESPACE_CLOSE
