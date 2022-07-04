// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2022 by the deal.II authors
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

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function_signed_distance.h>

#include <algorithm>

DEAL_II_NAMESPACE_OPEN

namespace Functions
{
  namespace SignedDistance
  {
    template <int dim>
    Sphere<dim>::Sphere(const Point<dim> &center, const double radius)
      : center(center)
      , radius(radius)
    {
      Assert(radius > 0, ExcMessage("Radius must be positive."))
    }



    template <int dim>
    double
    Sphere<dim>::value(const Point<dim> & point,
                       const unsigned int component) const
    {
      AssertIndexRange(component, this->n_components);
      (void)component;

      return point.distance(center) - radius;
    }



    template <int dim>
    Tensor<1, dim>
    Sphere<dim>::gradient(const Point<dim> & point,
                          const unsigned int component) const
    {
      AssertIndexRange(component, this->n_components);
      (void)component;

      const Tensor<1, dim> center_to_point = point - center;
      const Tensor<1, dim> grad = center_to_point / center_to_point.norm();
      return grad;
    }



    template <int dim>
    SymmetricTensor<2, dim>
    Sphere<dim>::hessian(const Point<dim> & point,
                         const unsigned int component) const
    {
      AssertIndexRange(component, this->n_components);
      (void)component;

      const Tensor<1, dim> center_to_point = point - center;
      const double         distance        = center_to_point.norm();

      const SymmetricTensor<2, dim> hess =
        unit_symmetric_tensor<dim>() / distance -
        symmetrize(outer_product(center_to_point, center_to_point)) /
          std::pow(distance, 3);

      return hess;
    }



    template <int dim>
    Plane<dim>::Plane(const Point<dim> &point, const Tensor<1, dim> &normal)
      : point_in_plane(point)
      , normal(normal)
    {
      Assert(normal.norm() > 0, ExcMessage("Plane normal must not be 0."))
    }



    template <int dim>
    double
    Plane<dim>::value(const Point<dim> & point,
                      const unsigned int component) const
    {
      AssertIndexRange(component, this->n_components);
      (void)component;

      return normal * (point - point_in_plane);
    }



    template <int dim>
    Tensor<1, dim>
    Plane<dim>::gradient(const Point<dim> &, const unsigned int component) const
    {
      AssertIndexRange(component, this->n_components);
      (void)component;

      return normal;
    }



    template <int dim>
    SymmetricTensor<2, dim>
    Plane<dim>::hessian(const Point<dim> &, const unsigned int component) const
    {
      AssertIndexRange(component, this->n_components);
      (void)component;

      return SymmetricTensor<2, dim>();
    }



    template <int dim>
    Ellipsoid<dim>::Ellipsoid(const Point<dim> &             center,
                              const std::array<double, dim> &radii,
                              const double                   tolerance,
                              const unsigned int             max_iter)
      : center(center)
      , radii(radii)
      , tolerance(tolerance)
      , max_iter(max_iter)
    {
      for (unsigned int d = 0; d < dim; ++d)
        Assert(radii[d] > 0, ExcMessage("All radii must be positive."))
    }



    template <int dim>
    double
    Ellipsoid<dim>::value(const Point<dim> & point,
                          const unsigned int component) const
    {
      AssertIndexRange(component, this->n_components);
      (void)component;

      if (dim == 1)
        return point.distance(center) - radii[0];
      else if (dim == 2)
        return compute_signed_distance_ellipse(point);
      else
        Assert(false, ExcNotImplemented());

      return 0.0;
    }



    template <int dim>
    Tensor<1, dim>
    Ellipsoid<dim>::gradient(const Point<dim> & point,
                             const unsigned int component) const
    {
      AssertIndexRange(component, this->n_components);
      (void)component;

      Tensor<1, dim> grad;
      if (dim == 1)
        grad = point - center;
      else if (dim == 2)
        {
          const Point<dim> point_in_centered_coordinate_system =
            Point<dim>(compute_closest_point_ellipse(point) - center);
          grad = compute_analyical_normal_vector_on_ellipse(
            point_in_centered_coordinate_system);
        }
      else
        AssertThrow(false, ExcNotImplemented());

      if (grad.norm() > 1e-12)
        return grad / grad.norm();
      else
        return grad;
    }



    template <int dim>
    double
    Ellipsoid<dim>::evaluate_ellipsoid(const Point<dim> &point) const
    {
      double val = 0.0;
      for (unsigned int d = 0; d < dim; ++d)
        val += std::pow((point[d] - center[d]) / radii[d], 2);
      return val - 1.0;
    }



    template <int dim>
    Point<dim>
    Ellipsoid<dim>::compute_closest_point_ellipse(const Point<dim> &point) const
    {
      AssertDimension(dim, 2);

      /*
       * Function to compute the closest point on an ellipse (adopted from
       * https://wet-robots.ghost.io/simple-method-for-distance-to-ellipse/ and
       * https://github.com/0xfaded/ellipse_demo):
       *
       * Since the ellipse is symmetric to the two major axes through its
       * center, the point is moved so the center coincides with the origin and
       * into the first quadrant.
       * 1. Choose a point on the ellipse (x), here x = a*cos(pi/4) and y =
       * b*sin(pi/4).
       * 2. Find second point on the ellipse, that has the same distance.
       * 3. Find midpoint on the ellipse (must be closer).
       * 4. Repeat 2.-4. until convergence.
       */
      // get equivalent point in first quadrant of centered ellipse
      const double px      = std::abs(point[0] - center[0]);
      const double py      = std::abs(point[1] - center[1]);
      const double sign_px = std::copysign(1.0, point[0] - center[0]);
      const double sign_py = std::copysign(1.0, point[1] - center[1]);
      // get semi axes radii
      const double &a = radii[0];
      const double &b = radii[1];
      // initial guess (t = angle from x-axis)
      double t = numbers::PI_4;
      double x = a * std::cos(t);
      double y = b * std::sin(t);

      unsigned int iter = 0;
      double       delta_t;
      do
        {
          // compute the ellipse evolute (center of curvature) for the current t
          const double ex = (a * a - b * b) * std::pow(std::cos(t), 3) / a;
          const double ey = (b * b - a * a) * std::pow(std::sin(t), 3) / b;
          // compute distances from current point on ellipse to its evolute
          const double rx = x - ex;
          const double ry = y - ey;
          // compute distances from point to the current evolute
          const double qx = px - ex;
          const double qy = py - ey;
          // compute the curvature radius at the current point on the ellipse
          const double r = std::hypot(rx, ry);
          // compute the distance from evolute to the point
          const double q = std::hypot(qx, qy);
          // compute step size on ellipse
          const double delta_c = r * std::asin((rx * qy - ry * qx) / (r * q));
          // compute approximate angle step
          delta_t = delta_c / std::sqrt(a * a + b * b - x * x - y * y);
          t += delta_t;
          // make sure the angle stays in first quadrant
          t = std::min(numbers::PI_2, std::max(0.0, t));
          x = a * std::cos(t);
          y = b * std::sin(t);
          iter++;
        }
      while (std::abs(delta_t) > tolerance && iter < max_iter);
      AssertIndexRange(iter, max_iter);

      AssertIsFinite(x);
      AssertIsFinite(y);

      return center + Point<dim>(sign_px * x, sign_py * y);
    }



    template <int dim>
    Tensor<1, dim, double>
    Ellipsoid<dim>::compute_analyical_normal_vector_on_ellipse(
      const Point<dim> &) const
    {
      AssertThrow(false, ExcNotImplemented());
      return Tensor<1, dim, double>();
    }



    template <>
    Tensor<1, 2, double>
    Ellipsoid<2>::compute_analyical_normal_vector_on_ellipse(
      const Point<2> &point) const
    {
      const auto &a = radii[0];
      const auto &b = radii[1];
      const auto &x = point[0];
      const auto &y = point[1];
      return Tensor<1, 2, double>({b * x / a, a * y / b});
    }



    template <int dim>
    double
    Ellipsoid<dim>::compute_signed_distance_ellipse(const Point<dim> &) const
    {
      AssertThrow(false, ExcNotImplemented());
      return 0;
    }



    template <>
    double
    Ellipsoid<2>::compute_signed_distance_ellipse(const Point<2> &point) const
    {
      // point corresponds to center
      if (point.distance(center) < tolerance)
        return *std::min_element(radii.begin(), radii.end()) * -1.;

      const Point<2> &closest_point = compute_closest_point_ellipse(point);

      const double distance =
        std::hypot(closest_point[0] - point[0], closest_point[1] - point[1]);

      return evaluate_ellipsoid(point) < 0.0 ? -distance : distance;
    }
  } // namespace SignedDistance
} // namespace Functions

#include "function_signed_distance.inst"

DEAL_II_NAMESPACE_CLOSE
