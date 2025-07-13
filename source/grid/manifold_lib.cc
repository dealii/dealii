// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/table.h>
#include <deal.II/base/tensor.h>

#include <deal.II/fe/mapping.h>
#include <deal.II/fe/mapping_q_internal.h>

#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <deal.II/physics/vector_relations.h>

#include <boost/container/small_vector.hpp>

#include <cmath>
#include <limits>
#include <memory>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  // The pull_back function fails regularly in the compute_chart_points
  // method, and, instead of throwing an exception, returns a point outside
  // the unit cell. The individual coordinates of that point are given by the
  // value below.
  static constexpr double invalid_pull_back_coordinate = 20.0;

  // Rotate a given unit vector u around the axis dir
  // where the angle is given by the length of dir.
  // This is the exponential map for a sphere.
  Tensor<1, 3>
  apply_exponential_map(const Tensor<1, 3> &u, const Tensor<1, 3> &dir)
  {
    const double theta = dir.norm();
    if (theta < 1.e-10)
      {
        return u;
      }
    else
      {
        const Tensor<1, 3> dirUnit = dir / theta;
        const Tensor<1, 3> tmp =
          std::cos(theta) * u + std::sin(theta) * dirUnit;
        return tmp / tmp.norm();
      }
  }

  // Returns the direction to go from v to u
  // projected to the plane perpendicular to the unit vector v.
  // This one is more stable when u and v are nearly equal.
  Tensor<1, 3>
  projected_direction(const Tensor<1, 3> &u, const Tensor<1, 3> &v)
  {
    Tensor<1, 3> ans = u - v;
    ans -= (ans * v) * v;
    return ans; // ans = (u-v) - ((u-v)*v)*v
  }

  // helper function to compute a vector orthogonal to a given one.
  // does nothing unless spacedim == 3.
  template <int spacedim>
  Point<spacedim>
  compute_normal(const Tensor<1, spacedim> & /*vector*/,
                 bool /*normalize*/ = false)
  {
    return {};
  }

  Point<3>
  compute_normal(const Tensor<1, 3> &vector, bool normalize = false)
  {
    Assert(vector.norm_square() != 0.,
           ExcMessage("The direction parameter must not be zero!"));
    Point<3> normal;
    if (std::abs(vector[0]) >= std::abs(vector[1]) &&
        std::abs(vector[0]) >= std::abs(vector[2]))
      {
        normal[1] = -1.;
        normal[2] = -1.;
        normal[0] = (vector[1] + vector[2]) / vector[0];
      }
    else if (std::abs(vector[1]) >= std::abs(vector[0]) &&
             std::abs(vector[1]) >= std::abs(vector[2]))
      {
        normal[0] = -1.;
        normal[2] = -1.;
        normal[1] = (vector[0] + vector[2]) / vector[1];
      }
    else
      {
        normal[0] = -1.;
        normal[1] = -1.;
        normal[2] = (vector[0] + vector[1]) / vector[2];
      }
    if (normalize)
      normal /= normal.norm();
    return normal;
  }
} // namespace internal



// ============================================================
// PolarManifold
// ============================================================

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
template <int dim, int spacedim>
PolarManifold<dim, spacedim>::PolarManifold(const Point<spacedim> center)
  : ChartManifold<dim, spacedim, spacedim>(
      PolarManifold<dim, spacedim>::get_periodicity())
  , center(center)
  , p_center(center)
{}
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS



template <int dim, int spacedim>
std::unique_ptr<Manifold<dim, spacedim>>
PolarManifold<dim, spacedim>::clone() const
{
  return std::make_unique<PolarManifold<dim, spacedim>>(*this);
}



template <int dim, int spacedim>
Tensor<1, spacedim>
PolarManifold<dim, spacedim>::get_periodicity()
{
  Tensor<1, spacedim> periodicity;
  // In two dimensions, theta is periodic.
  // In three dimensions things are a little more complicated, since the only
  // variable that is truly periodic is phi, while theta should be bounded
  // between 0 and pi. There is currently no way to enforce this, so here we
  // only fix periodicity for the last variable, corresponding to theta in 2d
  // and phi in 3d.
  periodicity[spacedim - 1] = 2 * numbers::PI;
  return periodicity;
}



template <int dim, int spacedim>
Point<spacedim>
PolarManifold<dim, spacedim>::push_forward(
  const Point<spacedim> &spherical_point) const
{
  Assert(spherical_point[0] >= 0.0,
         ExcMessage("Negative radius for given point."));
  const double rho   = spherical_point[0];
  const double theta = spherical_point[1];

  Point<spacedim> p;
  if (rho > 1e-10)
    switch (spacedim)
      {
        case 2:
          p[0] = rho * std::cos(theta);
          p[1] = rho * std::sin(theta);
          break;
        case 3:
          {
            const double phi = spherical_point[2];
            p[0]             = rho * std::sin(theta) * std::cos(phi);
            p[1]             = rho * std::sin(theta) * std::sin(phi);
            p[2]             = rho * std::cos(theta);
            break;
          }
        default:
          DEAL_II_NOT_IMPLEMENTED();
      }
  return p + p_center;
}



template <int dim, int spacedim>
Point<spacedim>
PolarManifold<dim, spacedim>::pull_back(
  const Point<spacedim> &space_point) const
{
  const Tensor<1, spacedim> R   = space_point - p_center;
  const double              rho = R.norm();

  Point<spacedim> p;
  p[0] = rho;

  switch (spacedim)
    {
      case 2:
        {
          p[1] = std::atan2(R[1], R[0]);
          if (p[1] < 0)
            p[1] += 2 * numbers::PI;
          break;
        }

      case 3:
        {
          const double z = R[2];
          p[2]           = std::atan2(R[1], R[0]); // phi
          if (p[2] < 0)
            p[2] += 2 * numbers::PI; // phi is periodic
          p[1] = std::atan2(std::sqrt(R[0] * R[0] + R[1] * R[1]), z); // theta
          break;
        }

      default:
        DEAL_II_NOT_IMPLEMENTED();
    }
  return p;
}



template <int dim, int spacedim>
DerivativeForm<1, spacedim, spacedim>
PolarManifold<dim, spacedim>::push_forward_gradient(
  const Point<spacedim> &spherical_point) const
{
  Assert(spherical_point[0] >= 0.0,
         ExcMessage("Negative radius for given point."));
  const double rho   = spherical_point[0];
  const double theta = spherical_point[1];

  DerivativeForm<1, spacedim, spacedim> DX;
  if (rho > 1e-10)
    switch (spacedim)
      {
        case 2:
          {
            DX[0][0] = std::cos(theta);
            DX[0][1] = -rho * std::sin(theta);
            DX[1][0] = std::sin(theta);
            DX[1][1] = rho * std::cos(theta);
            break;
          }

        case 3:
          {
            const double phi = spherical_point[2];
            DX[0][0]         = std::sin(theta) * std::cos(phi);
            DX[0][1]         = rho * std::cos(theta) * std::cos(phi);
            DX[0][2]         = -rho * std::sin(theta) * std::sin(phi);

            DX[1][0] = std::sin(theta) * std::sin(phi);
            DX[1][1] = rho * std::cos(theta) * std::sin(phi);
            DX[1][2] = rho * std::sin(theta) * std::cos(phi);

            DX[2][0] = std::cos(theta);
            DX[2][1] = -rho * std::sin(theta);
            DX[2][2] = 0;
            break;
          }

        default:
          DEAL_II_NOT_IMPLEMENTED();
      }
  return DX;
}



namespace
{
  template <int dim, int spacedim>
  bool
  spherical_face_is_horizontal(
    const typename Triangulation<dim, spacedim>::face_iterator &face,
    const Point<spacedim>                                      &manifold_center)
  {
    // We test whether a face is horizontal by checking that the vertices
    // all have roughly the same distance from the center: If the
    // maximum deviation for the distances from the vertices to the
    // center is less than 1.e-5 of the distance between vertices (as
    // measured by the minimum distance from any of the other vertices
    // to the first vertex), then we call this a horizontal face.
    constexpr unsigned int n_vertices =
      GeometryInfo<spacedim>::vertices_per_face;
    std::array<double, n_vertices>     sqr_distances_to_center;
    std::array<double, n_vertices - 1> sqr_distances_to_first_vertex;
    sqr_distances_to_center[0] =
      (face->vertex(0) - manifold_center).norm_square();
    for (unsigned int i = 1; i < n_vertices; ++i)
      {
        sqr_distances_to_center[i] =
          (face->vertex(i) - manifold_center).norm_square();
        sqr_distances_to_first_vertex[i - 1] =
          (face->vertex(i) - face->vertex(0)).norm_square();
      }
    const auto minmax_sqr_distance =
      std::minmax_element(sqr_distances_to_center.begin(),
                          sqr_distances_to_center.end());
    const auto min_sqr_distance_to_first_vertex =
      std::min_element(sqr_distances_to_first_vertex.begin(),
                       sqr_distances_to_first_vertex.end());

    return (*minmax_sqr_distance.second - *minmax_sqr_distance.first <
            1.e-10 * *min_sqr_distance_to_first_vertex);
  }
} // namespace



template <int dim, int spacedim>
Tensor<1, spacedim>
PolarManifold<dim, spacedim>::normal_vector(
  const typename Triangulation<dim, spacedim>::face_iterator &face,
  const Point<spacedim>                                      &p) const
{
  // Let us first test whether we are on a "horizontal" face
  // (tangential to the sphere).  In this case, the normal vector is
  // easy to compute since it is proportional to the vector from the
  // center to the point 'p'.
  if (spherical_face_is_horizontal<dim, spacedim>(face, p_center))
    {
      // So, if this is a "horizontal" face, then just compute the normal
      // vector as the one from the center to the point 'p', adequately
      // scaled.
      const Tensor<1, spacedim> unnormalized_spherical_normal = p - p_center;
      const Tensor<1, spacedim> normalized_spherical_normal =
        unnormalized_spherical_normal / unnormalized_spherical_normal.norm();
      return normalized_spherical_normal;
    }
  else
    // If it is not a horizontal face, just use the machinery of the
    // base class.
    return Manifold<dim, spacedim>::normal_vector(face, p);

  return Tensor<1, spacedim>();
}



// ============================================================
// SphericalManifold
// ============================================================

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
template <int dim, int spacedim>
SphericalManifold<dim, spacedim>::SphericalManifold(
  const Point<spacedim> center)
  : center(center)
  , p_center(center)
  , polar_manifold(center)
{}
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS



template <int dim, int spacedim>
std::unique_ptr<Manifold<dim, spacedim>>
SphericalManifold<dim, spacedim>::clone() const
{
  return std::make_unique<SphericalManifold<dim, spacedim>>(*this);
}



template <int dim, int spacedim>
Point<spacedim>
SphericalManifold<dim, spacedim>::get_intermediate_point(
  const Point<spacedim> &p1,
  const Point<spacedim> &p2,
  const double           w) const
{
  const double tol = 1e-10;

  if ((p1 - p2).norm_square() < tol * tol || std::abs(w) < tol)
    return p1;
  else if (std::abs(w - 1.0) < tol)
    return p2;

  // If the points are one dimensional then there is no need for anything but
  // a linear combination.
  if (spacedim == 1)
    return Point<spacedim>(w * p2 + (1 - w) * p1);

  const Tensor<1, spacedim> v1 = p1 - p_center;
  const Tensor<1, spacedim> v2 = p2 - p_center;
  const double              r1 = v1.norm();
  const double              r2 = v2.norm();

  Assert(r1 > tol && r2 > tol,
         ExcMessage("p1 and p2 cannot coincide with the center."));

  const Tensor<1, spacedim> e1 = v1 / r1;
  const Tensor<1, spacedim> e2 = v2 / r2;

  // Find the cosine of the angle gamma described by v1 and v2.
  const double cosgamma = e1 * e2;

  // Points are collinear with the center (allow for 8*eps as a tolerance)
  if (cosgamma < -1 + 8. * std::numeric_limits<double>::epsilon())
    return p_center;

  // Points are along a line, in which case e1 and e2 are essentially the same.
  if (cosgamma > 1 - 8. * std::numeric_limits<double>::epsilon())
    return Point<spacedim>(p_center + w * v2 + (1 - w) * v1);

  // Find the angle sigma that corresponds to arclength equal to w. acos
  // should never be undefined because we have ruled out the two special cases
  // above.
  const double sigma = w * std::acos(cosgamma);

  // Normal to v1 in the plane described by v1,v2,and the origin.
  // Since p1 and p2 do not coincide n is not zero and well defined.
  Tensor<1, spacedim> n      = v2 - (v2 * e1) * e1;
  const double        n_norm = n.norm();
  Assert(n_norm > 0,
         ExcInternalError("n should be different from the null vector. "
                          "Probably, this means v1==v2 or v2==0."));

  n /= n_norm;

  // Find the point Q along O,v1 such that
  // P1,V,P2 has measure sigma.
  const Tensor<1, spacedim> P = std::cos(sigma) * e1 + std::sin(sigma) * n;

  // Project this point on the manifold.
  return Point<spacedim>(p_center + (w * r2 + (1.0 - w) * r1) * P);
}



template <int dim, int spacedim>
Tensor<1, spacedim>
SphericalManifold<dim, spacedim>::get_tangent_vector(
  const Point<spacedim> &p1,
  const Point<spacedim> &p2) const
{
  [[maybe_unused]] const double tol = 1e-10;

  Assert(p1 != p2, ExcMessage("p1 and p2 should not concide."));

  const Tensor<1, spacedim> v1 = p1 - p_center;
  const Tensor<1, spacedim> v2 = p2 - p_center;
  const double              r1 = v1.norm();
  const double              r2 = v2.norm();

  Assert(r1 > tol, ExcMessage("p1 cannot coincide with the center."));

  Assert(r2 > tol, ExcMessage("p2 cannot coincide with the center."));

  const Tensor<1, spacedim> e1 = v1 / r1;
  const Tensor<1, spacedim> e2 = v2 / r2;

  // Find the cosine of the angle gamma described by v1 and v2.
  const double cosgamma = e1 * e2;

  Assert(cosgamma > -1 + 8. * std::numeric_limits<double>::epsilon(),
         ExcMessage("p1 and p2 cannot lie on the same diameter and be opposite "
                    "respect to the center."));

  if (cosgamma > 1 - 8. * std::numeric_limits<double>::epsilon())
    return v2 - v1;

  // Normal to v1 in the plane described by v1,v2,and the origin.
  // Since p1 and p2 do not coincide n is not zero and well defined.
  Tensor<1, spacedim> n      = v2 - (v2 * e1) * e1;
  const double        n_norm = n.norm();
  Assert(n_norm > 0,
         ExcInternalError("n should be different from the null vector. "
                          "Probably, this means v1==v2 or v2==0."));

  n /= n_norm;

  // this is the derivative of the geodesic in get_intermediate_point
  // derived with respect to w and inserting w=0.
  const double gamma = std::acos(cosgamma);
  return (r2 - r1) * e1 + r1 * gamma * n;
}



template <int dim, int spacedim>
Tensor<1, spacedim>
SphericalManifold<dim, spacedim>::normal_vector(
  const typename Triangulation<dim, spacedim>::face_iterator &face,
  const Point<spacedim>                                      &p) const
{
  // Let us first test whether we are on a "horizontal" face
  // (tangential to the sphere).  In this case, the normal vector is
  // easy to compute since it is proportional to the vector from the
  // center to the point 'p'.
  if (spherical_face_is_horizontal<dim, spacedim>(face, p_center))
    {
      // So, if this is a "horizontal" face, then just compute the normal
      // vector as the one from the center to the point 'p', adequately
      // scaled.
      const Tensor<1, spacedim> unnormalized_spherical_normal = p - p_center;
      const Tensor<1, spacedim> normalized_spherical_normal =
        unnormalized_spherical_normal / unnormalized_spherical_normal.norm();
      return normalized_spherical_normal;
    }
  else
    // If it is not a horizontal face, just use the machinery of the
    // base class.
    return Manifold<dim, spacedim>::normal_vector(face, p);

  return Tensor<1, spacedim>();
}



template <>
void
SphericalManifold<1, 1>::get_normals_at_vertices(
  const Triangulation<1>::face_iterator &,
  Manifold<1, 1>::FaceVertexNormals &) const
{
  Assert(false, ExcImpossibleInDim(1));
}



template <>
void
SphericalManifold<1, 2>::get_normals_at_vertices(
  const Triangulation<1, 2>::face_iterator &,
  Manifold<1, 2>::FaceVertexNormals &) const
{
  Assert(false, ExcImpossibleInDim(1));
}



template <int dim, int spacedim>
void
SphericalManifold<dim, spacedim>::get_normals_at_vertices(
  const typename Triangulation<dim, spacedim>::face_iterator &face,
  typename Manifold<dim, spacedim>::FaceVertexNormals &face_vertex_normals)
  const
{
  // Let us first test whether we are on a "horizontal" face
  // (tangential to the sphere).  In this case, the normal vector is
  // easy to compute since it is proportional to the vector from the
  // center to the point 'p'.
  if (spherical_face_is_horizontal<dim, spacedim>(face, p_center))
    {
      // So, if this is a "horizontal" face, then just compute the normal
      // vector as the one from the center to the point 'p', adequately
      // scaled.
      for (unsigned int vertex = 0;
           vertex < GeometryInfo<spacedim>::vertices_per_face;
           ++vertex)
        face_vertex_normals[vertex] = face->vertex(vertex) - p_center;
    }
  else
    Manifold<dim, spacedim>::get_normals_at_vertices(face, face_vertex_normals);
}



template <int dim, int spacedim>
void
SphericalManifold<dim, spacedim>::get_new_points(
  const ArrayView<const Point<spacedim>> &surrounding_points,
  const Table<2, double>                 &weights,
  ArrayView<Point<spacedim>>              new_points) const
{
  AssertDimension(new_points.size(), weights.size(0));
  AssertDimension(surrounding_points.size(), weights.size(1));

  do_get_new_points(surrounding_points, make_array_view(weights), new_points);

  return;
}



template <int dim, int spacedim>
Point<spacedim>
SphericalManifold<dim, spacedim>::get_new_point(
  const ArrayView<const Point<spacedim>> &vertices,
  const ArrayView<const double>          &weights) const
{
  // To avoid duplicating all of the logic in get_new_points, simply call it
  // for one position.
  Point<spacedim> new_point;
  do_get_new_points(vertices,
                    weights,
                    make_array_view(&new_point, &new_point + 1));

  return new_point;
}



namespace internal
{
  namespace SphericalManifold
  {
    namespace
    {
      template <int spacedim>
      Point<spacedim>
      do_get_new_point(
        const ArrayView<const Tensor<1, spacedim>> & /*directions*/,
        const ArrayView<const double> & /*distances*/,
        const ArrayView<const double> & /*weights*/,
        const Point<spacedim> & /*candidate_point*/)
      {
        DEAL_II_NOT_IMPLEMENTED();
        return {};
      }

      template <>
      Point<3>
      do_get_new_point(
        const ArrayView<const Tensor<1, 3>>            &directions,
        [[maybe_unused]] const ArrayView<const double> &distances,
        const ArrayView<const double>                  &weights,
        const Point<3>                                 &candidate_point)
      {
        AssertDimension(directions.size(), distances.size());
        AssertDimension(directions.size(), weights.size());

        Point<3>           candidate       = candidate_point;
        const unsigned int n_merged_points = directions.size();
        const double       tolerance       = 1e-10;
        const int          max_iterations  = 10;

        {
          // If the candidate happens to coincide with a normalized
          // direction, we return it. Otherwise, the Hessian would be singular.
          for (unsigned int i = 0; i < n_merged_points; ++i)
            {
              const double squared_distance =
                (candidate - directions[i]).norm_square();
              if (squared_distance < tolerance * tolerance)
                return candidate;
            }

          // check if we only have two points now, in which case we can use the
          // get_intermediate_point function
          if (n_merged_points == 2)
            {
              static const dealii::SphericalManifold<3, 3> unit_manifold;
              Assert(std::abs(weights[0] + weights[1] - 1.0) < 1e-13,
                     ExcMessage("Weights do not sum up to 1"));
              const Point<3> intermediate =
                unit_manifold.get_intermediate_point(Point<3>(directions[0]),
                                                     Point<3>(directions[1]),
                                                     weights[1]);
              return intermediate;
            }

          Tensor<1, 3> vPerp;
          Tensor<2, 2> Hessian;
          Tensor<1, 2> gradient;
          Tensor<1, 2> gradlocal;

          // On success we exit the loop early.
          // Otherwise, we just take the result after max_iterations steps.
          for (unsigned int i = 0; i < max_iterations; ++i)
            {
              // Step 2a: Find new descent direction

              // Get local basis for the estimate candidate
              const Tensor<1, 3> Clocalx = internal::compute_normal(candidate);
              const Tensor<1, 3> Clocaly = cross_product_3d(candidate, Clocalx);

              // For each vertices vector, compute the tangent vector from
              // candidate towards the vertices vector -- its length is the
              // spherical length from candidate to the vertices vector. Then
              // compute its contribution to the Hessian.
              gradient = 0.;
              Hessian  = 0.;
              for (unsigned int i = 0; i < n_merged_points; ++i)
                if (std::abs(weights[i]) > 1.e-15)
                  {
                    vPerp =
                      internal::projected_direction(directions[i], candidate);
                    const double sinthetaSq = vPerp.norm_square();
                    const double sintheta   = std::sqrt(sinthetaSq);
                    if (sintheta < tolerance)
                      {
                        Hessian[0][0] += weights[i];
                        Hessian[1][1] += weights[i];
                      }
                    else
                      {
                        const double costheta = (directions[i]) * candidate;
                        const double theta    = std::atan2(sintheta, costheta);
                        const double sincthetaInv = theta / sintheta;

                        const double cosphi = vPerp * Clocalx;
                        const double sinphi = vPerp * Clocaly;

                        gradlocal[0] = cosphi;
                        gradlocal[1] = sinphi;
                        gradient += (weights[i] * sincthetaInv) * gradlocal;

                        const double wt       = weights[i] / sinthetaSq;
                        const double sinphiSq = sinphi * sinphi;
                        const double cosphiSq = cosphi * cosphi;
                        const double tt       = sincthetaInv * costheta;
                        const double offdiag =
                          cosphi * sinphi * wt * (1.0 - tt);
                        Hessian[0][0] += wt * (cosphiSq + tt * sinphiSq);
                        Hessian[0][1] += offdiag;
                        Hessian[1][0] += offdiag;
                        Hessian[1][1] += wt * (sinphiSq + tt * cosphiSq);
                      }
                  }

              Assert(determinant(Hessian) > tolerance, ExcInternalError());

              const Tensor<2, 2> inverse_Hessian = invert(Hessian);

              const Tensor<1, 2> xDisplocal = inverse_Hessian * gradient;
              const Tensor<1, 3> xDisp =
                xDisplocal[0] * Clocalx + xDisplocal[1] * Clocaly;

              // Step 2b: rotate candidate in direction xDisp for a new
              // candidate.
              const Point<3> candidateOld = candidate;
              candidate =
                Point<3>(internal::apply_exponential_map(candidate, xDisp));

              // Step 2c: return the new candidate if we didn't move
              if ((candidate - candidateOld).norm_square() <
                  tolerance * tolerance)
                break;
            }
        }
        return candidate;
      }
    } // namespace
  }   // namespace SphericalManifold
} // namespace internal



template <int dim, int spacedim>
void
SphericalManifold<dim, spacedim>::do_get_new_points(
  const ArrayView<const Point<spacedim>> &surrounding_points,
  const ArrayView<const double>          &weights,
  ArrayView<Point<spacedim>>              new_points) const
{
  AssertDimension(weights.size(),
                  new_points.size() * surrounding_points.size());
  const unsigned int weight_rows    = new_points.size();
  const unsigned int weight_columns = surrounding_points.size();

  if (surrounding_points.size() == 2)
    {
      for (unsigned int row = 0; row < weight_rows; ++row)
        new_points[row] =
          SphericalManifold<dim, spacedim>::get_intermediate_point(
            surrounding_points[0],
            surrounding_points[1],
            weights[row * weight_columns + 1]);
      return;
    }

  boost::container::small_vector<std::pair<double, Tensor<1, spacedim>>, 100>
    new_candidates(new_points.size());
  boost::container::small_vector<Tensor<1, spacedim>, 100> directions(
    surrounding_points.size(), Point<spacedim>());
  boost::container::small_vector<double, 100> distances(
    surrounding_points.size(), 0.0);
  double max_distance = 0.;
  for (unsigned int i = 0; i < surrounding_points.size(); ++i)
    {
      directions[i] = surrounding_points[i] - p_center;
      distances[i]  = directions[i].norm();

      if (distances[i] != 0.)
        directions[i] /= distances[i];
      else
        Assert(false,
               ExcMessage("One of the vertices coincides with the center. "
                          "This is not allowed!"));

      // Check if an estimate is good enough,
      // this is often the case for sufficiently refined meshes.
      for (unsigned int k = 0; k < i; ++k)
        {
          const double squared_distance =
            (directions[i] - directions[k]).norm_square();
          max_distance = std::max(max_distance, squared_distance);
        }
    }

  // Step 1: Check for some special cases, create simple linear guesses
  // otherwise.
  const double                              tolerance = 1e-10;
  boost::container::small_vector<bool, 100> accurate_point_was_found(
    new_points.size(), false);
  const ArrayView<const Tensor<1, spacedim>> array_directions =
    make_array_view(directions.begin(), directions.end());
  const ArrayView<const double> array_distances =
    make_array_view(distances.begin(), distances.end());
  for (unsigned int row = 0; row < weight_rows; ++row)
    {
      new_candidates[row] =
        guess_new_point(array_directions,
                        array_distances,
                        ArrayView<const double>(&weights[row * weight_columns],
                                                weight_columns));

      // If the candidate is the center, mark it as found to avoid entering
      // the Newton iteration in step 2, which would crash.
      if (new_candidates[row].first == 0.0)
        {
          new_points[row]               = p_center;
          accurate_point_was_found[row] = true;
          continue;
        }

      // If not in 3d, just use the implementation from PolarManifold
      // after we verified that the candidate is not the center.
      if (spacedim < 3)
        new_points[row] = polar_manifold.get_new_point(
          surrounding_points,
          ArrayView<const double>(&weights[row * weight_columns],
                                  weight_columns));
    }

  // In this case, we treated the case that the candidate is the center and
  // obtained the new locations from the PolarManifold object otherwise.
  if constexpr (spacedim < 3)
    return;
  else
    {
      // If all the points are close to each other, we expect the estimate to
      // be good enough. This tolerance was chosen such that the first iteration
      // for a at least three time refined HyperShell mesh with radii .5 and 1.
      // doesn't already succeed.
      if (max_distance < 2e-2)
        {
          for (unsigned int row = 0; row < weight_rows; ++row)
            new_points[row] =
              p_center + new_candidates[row].first * new_candidates[row].second;

          return;
        }

      // Step 2:
      // Do more expensive Newton-style iterations to improve the estimate.

      // Search for duplicate directions and merge them to minimize the cost of
      // the get_new_point function call below.
      boost::container::small_vector<double, 1000> merged_weights(
        weights.size());
      boost::container::small_vector<Tensor<1, spacedim>, 100>
        merged_directions(surrounding_points.size(), Point<spacedim>());
      boost::container::small_vector<double, 100> merged_distances(
        surrounding_points.size(), 0.0);

      unsigned int n_unique_directions = 0;
      for (unsigned int i = 0; i < surrounding_points.size(); ++i)
        {
          bool found_duplicate = false;

          // This inner loop is of $O(N^2)$ complexity, but
          // surrounding_points.size() is usually at most 8 points large.
          for (unsigned int j = 0; j < n_unique_directions; ++j)
            {
              const double squared_distance =
                (directions[i] - directions[j]).norm_square();
              if (!found_duplicate && squared_distance < 1e-28)
                {
                  found_duplicate = true;
                  for (unsigned int row = 0; row < weight_rows; ++row)
                    merged_weights[row * weight_columns + j] +=
                      weights[row * weight_columns + i];
                }
            }

          if (found_duplicate == false)
            {
              merged_directions[n_unique_directions] = directions[i];
              merged_distances[n_unique_directions]  = distances[i];
              for (unsigned int row = 0; row < weight_rows; ++row)
                merged_weights[row * weight_columns + n_unique_directions] =
                  weights[row * weight_columns + i];

              ++n_unique_directions;
            }
        }

      // Search for duplicate weight rows and merge them to minimize the cost of
      // the get_new_point function call below.
      boost::container::small_vector<unsigned int, 100> merged_weights_index(
        new_points.size(), numbers::invalid_unsigned_int);
      for (unsigned int row = 0; row < weight_rows; ++row)
        {
          for (unsigned int existing_row = 0; existing_row < row;
               ++existing_row)
            {
              bool identical_weights = true;

              for (unsigned int weight_index = 0;
                   weight_index < n_unique_directions;
                   ++weight_index)
                if (std::abs(
                      merged_weights[row * weight_columns + weight_index] -
                      merged_weights[existing_row * weight_columns +
                                     weight_index]) > tolerance)
                  {
                    identical_weights = false;
                    break;
                  }

              if (identical_weights)
                {
                  merged_weights_index[row] = existing_row;
                  break;
                }
            }
        }

      // Note that we only use the n_unique_directions first entries in the
      // ArrayView
      const ArrayView<const Tensor<1, spacedim>> array_merged_directions =
        make_array_view(merged_directions.begin(),
                        merged_directions.begin() + n_unique_directions);
      const ArrayView<const double> array_merged_distances =
        make_array_view(merged_distances.begin(),
                        merged_distances.begin() + n_unique_directions);

      for (unsigned int row = 0; row < weight_rows; ++row)
        if (!accurate_point_was_found[row])
          {
            if (merged_weights_index[row] == numbers::invalid_unsigned_int)
              {
                const ArrayView<const double> array_merged_weights(
                  &merged_weights[row * weight_columns], n_unique_directions);
                new_candidates[row].second =
                  internal::SphericalManifold::do_get_new_point(
                    array_merged_directions,
                    array_merged_distances,
                    array_merged_weights,
                    Point<spacedim>(new_candidates[row].second));
              }
            else
              new_candidates[row].second =
                new_candidates[merged_weights_index[row]].second;

            new_points[row] =
              p_center + new_candidates[row].first * new_candidates[row].second;
          }
    }
}



template <int dim, int spacedim>
std::pair<double, Tensor<1, spacedim>>
SphericalManifold<dim, spacedim>::guess_new_point(
  const ArrayView<const Tensor<1, spacedim>> &directions,
  const ArrayView<const double>              &distances,
  const ArrayView<const double>              &weights) const
{
  const double        tolerance = 1e-10;
  double              rho       = 0.;
  Tensor<1, spacedim> candidate;

  // Perform a simple average ...
  double total_weights = 0.;
  for (unsigned int i = 0; i < directions.size(); ++i)
    {
      // if one weight is one, return its direction
      if (std::abs(1 - weights[i]) < tolerance)
        return std::make_pair(distances[i], directions[i]);

      rho += distances[i] * weights[i];
      candidate += directions[i] * weights[i];
      total_weights += weights[i];
    }

  // ... and normalize if the candidate is different from the origin.
  const double norm = candidate.norm();
  if (norm == 0.)
    return std::make_pair(0.0, Point<spacedim>());
  candidate /= norm;
  rho /= total_weights;

  return std::make_pair(rho, candidate);
}



// ============================================================
// CylindricalManifold
// ============================================================
template <int dim, int spacedim>
CylindricalManifold<dim, spacedim>::CylindricalManifold(const unsigned int axis,
                                                        const double tolerance)
  : CylindricalManifold<dim, spacedim>(Point<spacedim>::unit_vector(axis),
                                       Point<spacedim>(),
                                       tolerance)
{
  // do not use static_assert to make dimension-independent programming
  // easier.
  Assert(spacedim == 3,
         ExcMessage("CylindricalManifold can only be used for spacedim==3!"));
}



template <int dim, int spacedim>
CylindricalManifold<dim, spacedim>::CylindricalManifold(
  const Tensor<1, spacedim> &direction,
  const Point<spacedim>     &point_on_axis,
  const double               tolerance)
  : ChartManifold<dim, spacedim, 3>(Tensor<1, 3>({0, 2. * numbers::PI, 0}))
  , normal_direction(internal::compute_normal(direction, true))
  , direction(direction / direction.norm())
  , point_on_axis(point_on_axis)
  , tolerance(tolerance)
  , dxn(cross_product_3d(this->direction, normal_direction))
{
  // do not use static_assert to make dimension-independent programming
  // easier.
  Assert(spacedim == 3,
         ExcMessage("CylindricalManifold can only be used for spacedim==3!"));
}



template <int dim, int spacedim>
std::unique_ptr<Manifold<dim, spacedim>>
CylindricalManifold<dim, spacedim>::clone() const
{
  return std::make_unique<CylindricalManifold<dim, spacedim>>(*this);
}



template <int dim, int spacedim>
Point<spacedim>
CylindricalManifold<dim, spacedim>::get_new_point(
  const ArrayView<const Point<spacedim>> &surrounding_points,
  const ArrayView<const double>          &weights) const
{
  Assert(spacedim == 3,
         ExcMessage("CylindricalManifold can only be used for spacedim==3!"));

  // First check if the average in space lies on the axis.
  Point<spacedim> middle;
  double          average_length = 0.;
  for (unsigned int i = 0; i < surrounding_points.size(); ++i)
    {
      middle += surrounding_points[i] * weights[i];
      average_length += surrounding_points[i].square() * weights[i];
    }
  middle -= point_on_axis;
  const double lambda = middle * direction;

  if ((middle - direction * lambda).square() < tolerance * average_length)
    return point_on_axis + direction * lambda;
  else // If not, using the ChartManifold should yield valid results.
    return ChartManifold<dim, spacedim, 3>::get_new_point(surrounding_points,
                                                          weights);
}



template <int dim, int spacedim>
Point<3>
CylindricalManifold<dim, spacedim>::pull_back(
  const Point<spacedim> &space_point) const
{
  Assert(spacedim == 3,
         ExcMessage("CylindricalManifold can only be used for spacedim==3!"));

  // First find the projection of the given point to the axis.
  const Tensor<1, spacedim> normalized_point = space_point - point_on_axis;
  const double              lambda           = normalized_point * direction;
  const Point<spacedim>     projection = point_on_axis + direction * lambda;
  const Tensor<1, spacedim> p_diff     = space_point - projection;

  // Then compute the angle between the projection direction and
  // another vector orthogonal to the direction vector.
  const double phi = Physics::VectorRelations::signed_angle(normal_direction,
                                                            p_diff,
                                                            /*axis=*/direction);

  // Return distance from the axis, angle and signed distance on the axis.
  return Point<3>(p_diff.norm(), phi, lambda);
}



template <int dim, int spacedim>
Point<spacedim>
CylindricalManifold<dim, spacedim>::push_forward(
  const Point<3> &chart_point) const
{
  Assert(spacedim == 3,
         ExcMessage("CylindricalManifold can only be used for spacedim==3!"));

  // Rotate the orthogonal direction by the given angle
  const double sine_r   = std::sin(chart_point[1]) * chart_point[0];
  const double cosine_r = std::cos(chart_point[1]) * chart_point[0];
  const Tensor<1, spacedim> intermediate =
    normal_direction * cosine_r + dxn * sine_r;

  // Finally, put everything together.
  return point_on_axis + direction * chart_point[2] + intermediate;
}



template <int dim, int spacedim>
DerivativeForm<1, 3, spacedim>
CylindricalManifold<dim, spacedim>::push_forward_gradient(
  const Point<3> &chart_point) const
{
  Assert(spacedim == 3,
         ExcMessage("CylindricalManifold can only be used for spacedim==3!"));

  Tensor<2, 3> derivatives;

  // Rotate the orthogonal direction by the given angle
  const double              sine   = std::sin(chart_point[1]);
  const double              cosine = std::cos(chart_point[1]);
  const Tensor<1, spacedim> intermediate =
    normal_direction * cosine + dxn * sine;

  // avoid compiler warnings
  constexpr int s0 = 0 % spacedim;
  constexpr int s1 = 1 % spacedim;
  constexpr int s2 = 2 % spacedim;

  // derivative w.r.t the radius
  derivatives[s0][s0] = intermediate[s0];
  derivatives[s1][s0] = intermediate[s1];
  derivatives[s2][s0] = intermediate[s2];

  // derivatives w.r.t the angle
  derivatives[s0][s1] = -normal_direction[s0] * sine + dxn[s0] * cosine;
  derivatives[s1][s1] = -normal_direction[s1] * sine + dxn[s1] * cosine;
  derivatives[s2][s1] = -normal_direction[s2] * sine + dxn[s2] * cosine;

  // derivatives w.r.t the direction of the axis
  derivatives[s0][s2] = direction[s0];
  derivatives[s1][s2] = direction[s1];
  derivatives[s2][s2] = direction[s2];

  return derivatives;
}



namespace
{
  template <int dim>
  Tensor<1, dim>
  check_and_normalize(const Tensor<1, dim> &t)
  {
    const double norm = t.norm();
    Assert(norm > 0.0, ExcMessage("The major axis must have a positive norm."));
    return t / norm;
  }
} // namespace



// ============================================================
// EllipticalManifold
// ============================================================
template <int dim, int spacedim>
EllipticalManifold<dim, spacedim>::EllipticalManifold(
  const Point<spacedim>     &center,
  const Tensor<1, spacedim> &major_axis_direction,
  const double               eccentricity)
  : ChartManifold<dim, spacedim, spacedim>(
      EllipticalManifold<dim, spacedim>::get_periodicity())
  , direction(check_and_normalize(major_axis_direction))
  , center(center)
  , eccentricity(eccentricity)
  , cosh_u(1.0 / eccentricity)
  , sinh_u(std::sqrt(cosh_u * cosh_u - 1.0))
{
  // Throws an exception if dim!=2 || spacedim!=2.
  Assert(dim == 2 && spacedim == 2, ExcNotImplemented());
  // Throws an exception if eccentricity is not in range.
  Assert(std::signbit(cosh_u * cosh_u - 1.0) == false,
         ExcMessage(
           "Invalid eccentricity: It must satisfy 0 < eccentricity < 1."));
}



template <int dim, int spacedim>
std::unique_ptr<Manifold<dim, spacedim>>
EllipticalManifold<dim, spacedim>::clone() const
{
  return std::make_unique<EllipticalManifold<dim, spacedim>>(*this);
}



template <int dim, int spacedim>
Tensor<1, spacedim>
EllipticalManifold<dim, spacedim>::get_periodicity()
{
  Tensor<1, spacedim> periodicity;
  // The second elliptical coordinate is periodic, while the first is not.
  // Enforce periodicity on the last variable.
  periodicity[spacedim - 1] = 2.0 * numbers::PI;
  return periodicity;
}



template <int dim, int spacedim>
Point<spacedim>
EllipticalManifold<dim, spacedim>::push_forward(const Point<spacedim> &) const
{
  DEAL_II_NOT_IMPLEMENTED();
  return {};
}



template <>
Point<2>
EllipticalManifold<2, 2>::push_forward(const Point<2> &chart_point) const
{
  const double cs = std::cos(chart_point[1]);
  const double sn = std::sin(chart_point[1]);
  // Coordinates in the reference frame (i.e. major axis direction is
  // x-axis)
  const double x = chart_point[0] * cosh_u * cs;
  const double y = chart_point[0] * sinh_u * sn;
  // Rotate them according to the major axis direction
  const Point<2> p(direction[0] * x - direction[1] * y,
                   direction[1] * x + direction[0] * y);
  return p + center;
}



template <int dim, int spacedim>
Point<spacedim>
EllipticalManifold<dim, spacedim>::pull_back(const Point<spacedim> &) const
{
  DEAL_II_NOT_IMPLEMENTED();
  return {};
}



template <>
Point<2>
EllipticalManifold<2, 2>::pull_back(const Point<2> &space_point) const
{
  // Moving space_point in the reference coordinate system.
  const double x0 = space_point[0] - center[0];
  const double y0 = space_point[1] - center[1];
  const double x  = direction[0] * x0 + direction[1] * y0;
  const double y  = -direction[1] * x0 + direction[0] * y0;
  const double pt0 =
    std::sqrt((x * x) / (cosh_u * cosh_u) + (y * y) / (sinh_u * sinh_u));
  // If the radius is exactly zero, the point coincides with the origin.
  if (pt0 == 0.0)
    {
      return center;
    }
  double cos_eta = x / (pt0 * cosh_u);
  if (cos_eta < -1.0)
    {
      cos_eta = -1.0;
    }
  if (cos_eta > 1.0)
    {
      cos_eta = 1.0;
    }
  const double eta = std::acos(cos_eta);
  const double pt1 = (std::signbit(y) ? 2.0 * numbers::PI - eta : eta);
  return {pt0, pt1};
}



template <int dim, int spacedim>
DerivativeForm<1, spacedim, spacedim>
EllipticalManifold<dim, spacedim>::push_forward_gradient(
  const Point<spacedim> &) const
{
  DEAL_II_NOT_IMPLEMENTED();
  return {};
}



template <>
DerivativeForm<1, 2, 2>
EllipticalManifold<2, 2>::push_forward_gradient(
  const Point<2> &chart_point) const
{
  const double cs = std::cos(chart_point[1]);
  const double sn = std::sin(chart_point[1]);
  Tensor<2, 2> dX;
  dX[0][0] = cosh_u * cs;
  dX[0][1] = -chart_point[0] * cosh_u * sn;
  dX[1][0] = sinh_u * sn;
  dX[1][1] = chart_point[0] * sinh_u * cs;

  // rotate according to the major axis direction
  Tensor<2, 2, double> rot{
    {{+direction[0], -direction[1]}, {direction[1], direction[0]}}};

  return rot * dX;
}



// ============================================================
// FunctionManifold
// ============================================================
template <int dim, int spacedim, int chartdim>
FunctionManifold<dim, spacedim, chartdim>::FunctionManifold(
  const Function<chartdim>  &push_forward_function,
  const Function<spacedim>  &pull_back_function,
  const Tensor<1, chartdim> &periodicity,
  const double               tolerance)
  : ChartManifold<dim, spacedim, chartdim>(periodicity)
  , const_map()
  , push_forward_function(&push_forward_function)
  , pull_back_function(&pull_back_function)
  , tolerance(tolerance)
  , owns_pointers(false)
  , finite_difference_step(0)
{
  AssertDimension(push_forward_function.n_components, spacedim);
  AssertDimension(pull_back_function.n_components, chartdim);
}



template <int dim, int spacedim, int chartdim>
FunctionManifold<dim, spacedim, chartdim>::FunctionManifold(
  std::unique_ptr<Function<chartdim>> push_forward,
  std::unique_ptr<Function<spacedim>> pull_back,
  const Tensor<1, chartdim>          &periodicity,
  const double                        tolerance)
  : ChartManifold<dim, spacedim, chartdim>(periodicity)
  , const_map()
  , push_forward_function(push_forward.release())
  , pull_back_function(pull_back.release())
  , tolerance(tolerance)
  , owns_pointers(true)
  , finite_difference_step(0)
{
  AssertDimension(push_forward_function->n_components, spacedim);
  AssertDimension(pull_back_function->n_components, chartdim);
}



template <int dim, int spacedim, int chartdim>
FunctionManifold<dim, spacedim, chartdim>::FunctionManifold(
  const std::string                                 push_forward_expression,
  const std::string                                 pull_back_expression,
  const Tensor<1, chartdim>                        &periodicity,
  const typename FunctionParser<spacedim>::ConstMap const_map,
  const std::string                                 chart_vars,
  const std::string                                 space_vars,
  const double                                      tolerance,
  const double                                      h)
  : ChartManifold<dim, spacedim, chartdim>(periodicity)
  , const_map(const_map)
  , tolerance(tolerance)
  , owns_pointers(true)
  , push_forward_expression(push_forward_expression)
  , pull_back_expression(pull_back_expression)
  , chart_vars(chart_vars)
  , space_vars(space_vars)
  , finite_difference_step(h)
{
  FunctionParser<chartdim> *pf = new FunctionParser<chartdim>(spacedim, 0.0, h);
  FunctionParser<spacedim> *pb = new FunctionParser<spacedim>(chartdim, 0.0, h);
  pf->initialize(chart_vars, push_forward_expression, const_map);
  pb->initialize(space_vars, pull_back_expression, const_map);
  push_forward_function = pf;
  pull_back_function    = pb;
}



template <int dim, int spacedim, int chartdim>
FunctionManifold<dim, spacedim, chartdim>::~FunctionManifold()
{
  if (owns_pointers == true)
    {
      const Function<chartdim> *pf = push_forward_function;
      push_forward_function        = nullptr;
      delete pf;

      const Function<spacedim> *pb = pull_back_function;
      pull_back_function           = nullptr;
      delete pb;
    }
}



template <int dim, int spacedim, int chartdim>
std::unique_ptr<Manifold<dim, spacedim>>
FunctionManifold<dim, spacedim, chartdim>::clone() const
{
  // This manifold can be constructed either by providing an expression for the
  // push forward and the pull back charts, or by providing two Function
  // objects. In the first case, the push_forward and pull_back functions are
  // created internally in FunctionManifold, and destroyed when this object is
  // deleted. In the second case, the function objects are destroyed if they
  // are passed as pointers upon construction.
  // We need to make sure that our cloned object is constructed in the
  // same way this class was constructed, and that its internal Function
  // pointers point either to the same Function objects used to construct this
  // function or that the newly generated manifold creates internally the
  // push_forward and pull_back functions using the same expressions that were
  // used to construct this class.
  if (!(push_forward_expression.empty() && pull_back_expression.empty()))
    {
      return std::make_unique<FunctionManifold<dim, spacedim, chartdim>>(
        push_forward_expression,
        pull_back_expression,
        this->get_periodicity(),
        const_map,
        chart_vars,
        space_vars,
        tolerance,
        finite_difference_step);
    }
  else
    {
      return std::make_unique<FunctionManifold<dim, spacedim, chartdim>>(
        *push_forward_function,
        *pull_back_function,
        this->get_periodicity(),
        tolerance);
    }
}



template <int dim, int spacedim, int chartdim>
Point<spacedim>
FunctionManifold<dim, spacedim, chartdim>::push_forward(
  const Point<chartdim> &chart_point) const
{
  Vector<double>  pf(spacedim);
  Point<spacedim> result;
  push_forward_function->vector_value(chart_point, pf);
  for (unsigned int i = 0; i < spacedim; ++i)
    result[i] = pf[i];

  if constexpr (running_in_debug_mode())
    {
      Vector<double> pb(chartdim);
      pull_back_function->vector_value(result, pb);
      for (unsigned int i = 0; i < chartdim; ++i)
        Assert(
          (chart_point.norm() > tolerance &&
           (std::abs(pb[i] - chart_point[i]) <
            tolerance * chart_point.norm())) ||
            (std::abs(pb[i] - chart_point[i]) < tolerance),
          ExcMessage(
            "The push forward is not the inverse of the pull back! Bailing out."));
    }

  return result;
}



template <int dim, int spacedim, int chartdim>
DerivativeForm<1, chartdim, spacedim>
FunctionManifold<dim, spacedim, chartdim>::push_forward_gradient(
  const Point<chartdim> &chart_point) const
{
  DerivativeForm<1, chartdim, spacedim> DF;
  for (unsigned int i = 0; i < spacedim; ++i)
    {
      const auto gradient = push_forward_function->gradient(chart_point, i);
      for (unsigned int j = 0; j < chartdim; ++j)
        DF[i][j] = gradient[j];
    }
  return DF;
}



template <int dim, int spacedim, int chartdim>
Point<chartdim>
FunctionManifold<dim, spacedim, chartdim>::pull_back(
  const Point<spacedim> &space_point) const
{
  Vector<double>  pb(chartdim);
  Point<chartdim> result;
  pull_back_function->vector_value(space_point, pb);
  for (unsigned int i = 0; i < chartdim; ++i)
    result[i] = pb[i];
  return result;
}



// ============================================================
// TorusManifold
// ============================================================
template <int dim>
Point<3>
TorusManifold<dim>::pull_back(const Point<3> &p) const
{
  double x     = p[0];
  double z     = p[1];
  double y     = p[2];
  double phi   = std::atan2(y, x);
  double theta = std::atan2(z, std::sqrt(x * x + y * y) - centerline_radius);
  double w =
    std::sqrt(Utilities::fixed_power<2>(y - std::sin(phi) * centerline_radius) +
              Utilities::fixed_power<2>(x - std::cos(phi) * centerline_radius) +
              z * z) /
    inner_radius;
  return {phi, theta, w};
}



template <int dim>
Point<3>
TorusManifold<dim>::push_forward(const Point<3> &chart_point) const
{
  double phi   = chart_point[0];
  double theta = chart_point[1];
  double w     = chart_point[2];

  return {std::cos(phi) * centerline_radius +
            inner_radius * w * std::cos(theta) * std::cos(phi),
          inner_radius * w * std::sin(theta),
          std::sin(phi) * centerline_radius +
            inner_radius * w * std::cos(theta) * std::sin(phi)};
}



template <int dim>
TorusManifold<dim>::TorusManifold(const double centerline_radius,
                                  const double inner_radius)
  : ChartManifold<dim, 3, 3>(Point<3>(2 * numbers::PI, 2 * numbers::PI, 0.0))
  , centerline_radius(centerline_radius)
  , inner_radius(inner_radius)
{
  Assert(centerline_radius > inner_radius,
         ExcMessage("The centerline radius must be greater than the "
                    "inner radius."));
  Assert(inner_radius > 0.0, ExcMessage("The inner radius must be positive."));
}



template <int dim>
std::unique_ptr<Manifold<dim, 3>>
TorusManifold<dim>::clone() const
{
  return std::make_unique<TorusManifold<dim>>(centerline_radius, inner_radius);
}



template <int dim>
DerivativeForm<1, 3, 3>
TorusManifold<dim>::push_forward_gradient(const Point<3> &chart_point) const
{
  DerivativeForm<1, spacedim, spacedim> DX;

  double phi   = chart_point[0];
  double theta = chart_point[1];
  double w     = chart_point[2];

  DX[0][0] = -std::sin(phi) * centerline_radius -
             inner_radius * w * std::cos(theta) * std::sin(phi);
  DX[0][1] = -inner_radius * w * std::sin(theta) * std::cos(phi);
  DX[0][2] = inner_radius * std::cos(theta) * std::cos(phi);

  DX[1][0] = 0;
  DX[1][1] = inner_radius * w * std::cos(theta);
  DX[1][2] = inner_radius * std::sin(theta);

  DX[2][0] = std::cos(phi) * centerline_radius +
             inner_radius * w * std::cos(theta) * std::cos(phi);
  DX[2][1] = -inner_radius * w * std::sin(theta) * std::sin(phi);
  DX[2][2] = inner_radius * std::cos(theta) * std::sin(phi);

  return DX;
}



// ============================================================
// TransfiniteInterpolationManifold
// ============================================================
template <int dim, int spacedim>
TransfiniteInterpolationManifold<dim,
                                 spacedim>::TransfiniteInterpolationManifold()
  : triangulation(nullptr)
  , level_coarse(-1)
{
  AssertThrow(dim > 1, ExcNotImplemented());
}



template <int dim, int spacedim>
TransfiniteInterpolationManifold<dim,
                                 spacedim>::~TransfiniteInterpolationManifold()
{
  if (clear_signal.connected())
    clear_signal.disconnect();
}



template <int dim, int spacedim>
std::unique_ptr<Manifold<dim, spacedim>>
TransfiniteInterpolationManifold<dim, spacedim>::clone() const
{
  auto ptr = new TransfiniteInterpolationManifold<dim, spacedim>();
  if (triangulation)
    ptr->initialize(*triangulation);
  return std::unique_ptr<Manifold<dim, spacedim>>(ptr);
}



template <int dim, int spacedim>
void
TransfiniteInterpolationManifold<dim, spacedim>::initialize(
  const Triangulation<dim, spacedim> &triangulation)
{
  this->triangulation = &triangulation;
  // In case the triangulation is cleared, remove the pointers by a signal:
  clear_signal.disconnect();
  clear_signal = triangulation.signals.clear.connect([&]() -> void {
    this->triangulation = nullptr;
    this->level_coarse  = -1;
  });
  level_coarse = triangulation.last()->level();
  coarse_cell_is_flat.resize(triangulation.n_cells(level_coarse), false);
  quadratic_approximation.clear();

  // In case of dim == spacedim we perform a quadratic approximation in
  // InverseQuadraticApproximation(), thus initialize the unit_points
  // vector with one subdivision to get 3^dim unit_points.
  //
  // In the co-dimension one case (meaning  dim < spacedim) we have to fall
  // back to a simple GridTools::affine_cell_approximation<dim>() which
  // requires 2^dim points, instead. Thus, initialize the QIterated
  // quadrature with no subdivisions.
  const std::vector<Point<dim>> unit_points =
    QIterated<dim>(QTrapezoid<1>(), (dim == spacedim ? 2 : 1)).get_points();
  std::vector<Point<spacedim>> real_points(unit_points.size());

  for (const auto &cell : triangulation.active_cell_iterators())
    {
      bool cell_is_flat = true;
      for (const auto l : cell->line_indices())
        if (cell->line(l)->manifold_id() != cell->manifold_id() &&
            cell->line(l)->manifold_id() != numbers::flat_manifold_id)
          cell_is_flat = false;
      if constexpr (dim > 2)
        for (const auto q : cell->face_indices())
          if (cell->quad(q)->manifold_id() != cell->manifold_id() &&
              cell->quad(q)->manifold_id() != numbers::flat_manifold_id)
            cell_is_flat = false;
      AssertIndexRange(static_cast<unsigned int>(cell->index()),
                       coarse_cell_is_flat.size());
      coarse_cell_is_flat[cell->index()] = cell_is_flat;

      // build quadratic interpolation
      for (unsigned int i = 0; i < unit_points.size(); ++i)
        real_points[i] = push_forward(cell, unit_points[i]);
      quadratic_approximation.emplace_back(real_points, unit_points);
    }
}



namespace
{
  // version for 1d
  template <typename AccessorType>
  Point<AccessorType::space_dimension>
  compute_transfinite_interpolation(const AccessorType &cell,
                                    const Point<1>     &chart_point,
                                    const bool /*cell_is_flat*/)
  {
    return cell.vertex(0) * (1. - chart_point[0]) +
           cell.vertex(1) * chart_point[0];
  }

  // version for 2d
  template <typename AccessorType>
  Point<AccessorType::space_dimension>
  compute_transfinite_interpolation(const AccessorType &cell,
                                    const Point<2>     &chart_point,
                                    const bool          cell_is_flat)
  {
    const unsigned int       dim             = AccessorType::dimension;
    const unsigned int       spacedim        = AccessorType::space_dimension;
    const types::manifold_id my_manifold_id  = cell.manifold_id();
    const Triangulation<dim, spacedim> &tria = cell.get_triangulation();

    // formula see wikipedia
    // https://en.wikipedia.org/wiki/Transfinite_interpolation
    // S(u,v) = (1-v)c_1(u)+v c_3(u) + (1-u)c_2(v) + u c_4(v) -
    //   [(1-u)(1-v)P_0 + u(1-v) P_1 + (1-u)v P_2 + uv P_3]
    const std::array<Point<spacedim>, 4> vertices{
      {cell.vertex(0), cell.vertex(1), cell.vertex(2), cell.vertex(3)}};

    // this evaluates all bilinear shape functions because we need them
    // repeatedly. we will update this values in the complicated case with
    // curved lines below
    std::array<double, 4> weights_vertices{
      {(1. - chart_point[0]) * (1. - chart_point[1]),
       chart_point[0] * (1. - chart_point[1]),
       (1. - chart_point[0]) * chart_point[1],
       chart_point[0] * chart_point[1]}};

    Point<spacedim> new_point;
    if (cell_is_flat)
      for (const unsigned int v : GeometryInfo<2>::vertex_indices())
        new_point += weights_vertices[v] * vertices[v];
    else
      {
        // The second line in the formula tells us to subtract the
        // contribution of the vertices.  If a line employs the same manifold
        // as the cell, we can merge the weights of the line with the weights
        // of the vertex with a negative sign while going through the faces
        // (this is a bit artificial in 2d but it becomes clear in 3d where we
        // avoid looking at the faces' orientation and other complications).

        // add the contribution from the lines around the cell (first line in
        // formula)
        std::array<double, GeometryInfo<2>::vertices_per_face>          weights;
        std::array<Point<spacedim>, GeometryInfo<2>::vertices_per_face> points;
        // note that the views are immutable, but the arrays are not
        const auto weights_view =
          make_array_view(weights.begin(), weights.end());
        const auto points_view = make_array_view(points.begin(), points.end());

        for (unsigned int line = 0; line < GeometryInfo<2>::lines_per_cell;
             ++line)
          {
            const double my_weight =
              (line % 2) ? chart_point[line / 2] : 1 - chart_point[line / 2];
            const double line_point = chart_point[1 - line / 2];

            // Same manifold or invalid id which will go back to the same
            // class -> contribution should be added for the final point,
            // which means that we subtract the current weight from the
            // negative weight applied to the vertex
            const types::manifold_id line_manifold_id =
              cell.line(line)->manifold_id();
            if (line_manifold_id == my_manifold_id ||
                line_manifold_id == numbers::flat_manifold_id)
              {
                weights_vertices[GeometryInfo<2>::line_to_cell_vertices(line,
                                                                        0)] -=
                  my_weight * (1. - line_point);
                weights_vertices[GeometryInfo<2>::line_to_cell_vertices(line,
                                                                        1)] -=
                  my_weight * line_point;
              }
            else
              {
                points[0] =
                  vertices[GeometryInfo<2>::line_to_cell_vertices(line, 0)];
                points[1] =
                  vertices[GeometryInfo<2>::line_to_cell_vertices(line, 1)];
                weights[0] = 1. - line_point;
                weights[1] = line_point;
                new_point +=
                  my_weight * tria.get_manifold(line_manifold_id)
                                .get_new_point(points_view, weights_view);
              }
          }

        // subtract contribution from the vertices (second line in formula)
        for (const unsigned int v : GeometryInfo<2>::vertex_indices())
          new_point -= weights_vertices[v] * vertices[v];
      }

    return new_point;
  }

  // this is replicated from GeometryInfo::face_to_cell_vertices since we need
  // it very often in compute_transfinite_interpolation and the function is
  // performance critical
  static constexpr unsigned int face_to_cell_vertices_3d[6][4] = {{0, 2, 4, 6},
                                                                  {1, 3, 5, 7},
                                                                  {0, 4, 1, 5},
                                                                  {2, 6, 3, 7},
                                                                  {0, 1, 2, 3},
                                                                  {4, 5, 6, 7}};

  // this is replicated from GeometryInfo::face_to_cell_lines since we need it
  // very often in compute_transfinite_interpolation and the function is
  // performance critical
  static constexpr unsigned int face_to_cell_lines_3d[6][4] = {{8, 10, 0, 4},
                                                               {9, 11, 1, 5},
                                                               {2, 6, 8, 9},
                                                               {3, 7, 10, 11},
                                                               {0, 1, 2, 3},
                                                               {4, 5, 6, 7}};

  // version for 3d
  template <typename AccessorType>
  Point<AccessorType::space_dimension>
  compute_transfinite_interpolation(const AccessorType &cell,
                                    const Point<3>     &chart_point,
                                    const bool          cell_is_flat)
  {
    const unsigned int       dim             = AccessorType::dimension;
    const unsigned int       spacedim        = AccessorType::space_dimension;
    const types::manifold_id my_manifold_id  = cell.manifold_id();
    const Triangulation<dim, spacedim> &tria = cell.get_triangulation();

    // Same approach as in 2d, but adding the faces, subtracting the edges, and
    // adding the vertices
    const std::array<Point<spacedim>, 8> vertices{{cell.vertex(0),
                                                   cell.vertex(1),
                                                   cell.vertex(2),
                                                   cell.vertex(3),
                                                   cell.vertex(4),
                                                   cell.vertex(5),
                                                   cell.vertex(6),
                                                   cell.vertex(7)}};

    // store the components of the linear shape functions because we need them
    // repeatedly. we allow for 10 such shape functions to wrap around the
    // first four once again for easier face access.
    double linear_shapes[10];
    for (unsigned int d = 0; d < 3; ++d)
      {
        linear_shapes[2 * d]     = 1. - chart_point[d];
        linear_shapes[2 * d + 1] = chart_point[d];
      }

    // wrap linear shape functions around for access in face loop
    for (unsigned int d = 6; d < 10; ++d)
      linear_shapes[d] = linear_shapes[d - 6];

    std::array<double, 8> weights_vertices;
    for (unsigned int i2 = 0, v = 0; i2 < 2; ++i2)
      for (unsigned int i1 = 0; i1 < 2; ++i1)
        for (unsigned int i0 = 0; i0 < 2; ++i0, ++v)
          weights_vertices[v] =
            (linear_shapes[4 + i2] * linear_shapes[2 + i1]) * linear_shapes[i0];

    Point<spacedim> new_point;
    if (cell_is_flat)
      for (unsigned int v = 0; v < 8; ++v)
        new_point += weights_vertices[v] * vertices[v];
    else
      {
        // identify the weights for the lines to be accumulated (vertex
        // weights are set outside and coincide with the flat manifold case)

        std::array<double, GeometryInfo<3>::lines_per_cell> weights_lines;
        std::fill(weights_lines.begin(), weights_lines.end(), 0.0);

        // start with the contributions of the faces
        std::array<double, GeometryInfo<2>::vertices_per_cell>          weights;
        std::array<Point<spacedim>, GeometryInfo<2>::vertices_per_cell> points;
        // note that the views are immutable, but the arrays are not
        const auto weights_view =
          make_array_view(weights.begin(), weights.end());
        const auto points_view = make_array_view(points.begin(), points.end());

        for (const unsigned int face : GeometryInfo<3>::face_indices())
          {
            const double       my_weight = linear_shapes[face];
            const unsigned int face_even = face - face % 2;

            if (std::abs(my_weight) < 1e-13)
              continue;

            // same manifold or invalid id which will go back to the same class
            // -> face will interpolate from the surrounding lines and vertices
            const types::manifold_id face_manifold_id =
              cell.face(face)->manifold_id();
            if (face_manifold_id == my_manifold_id ||
                face_manifold_id == numbers::flat_manifold_id)
              {
                for (unsigned int line = 0;
                     line < GeometryInfo<2>::lines_per_cell;
                     ++line)
                  {
                    const double line_weight =
                      linear_shapes[face_even + 2 + line];
                    weights_lines[face_to_cell_lines_3d[face][line]] +=
                      my_weight * line_weight;
                  }
                // as to the indices inside linear_shapes: we use the index
                // wrapped around at 2*d, ensuring the correct orientation of
                // the face's coordinate system with respect to the
                // lexicographic indices
                weights_vertices[face_to_cell_vertices_3d[face][0]] -=
                  linear_shapes[face_even + 2] *
                  (linear_shapes[face_even + 4] * my_weight);
                weights_vertices[face_to_cell_vertices_3d[face][1]] -=
                  linear_shapes[face_even + 3] *
                  (linear_shapes[face_even + 4] * my_weight);
                weights_vertices[face_to_cell_vertices_3d[face][2]] -=
                  linear_shapes[face_even + 2] *
                  (linear_shapes[face_even + 5] * my_weight);
                weights_vertices[face_to_cell_vertices_3d[face][3]] -=
                  linear_shapes[face_even + 3] *
                  (linear_shapes[face_even + 5] * my_weight);
              }
            else
              {
                for (const unsigned int v : GeometryInfo<2>::vertex_indices())
                  points[v] = vertices[face_to_cell_vertices_3d[face][v]];
                weights[0] =
                  linear_shapes[face_even + 2] * linear_shapes[face_even + 4];
                weights[1] =
                  linear_shapes[face_even + 3] * linear_shapes[face_even + 4];
                weights[2] =
                  linear_shapes[face_even + 2] * linear_shapes[face_even + 5];
                weights[3] =
                  linear_shapes[face_even + 3] * linear_shapes[face_even + 5];
                new_point +=
                  my_weight * tria.get_manifold(face_manifold_id)
                                .get_new_point(points_view, weights_view);
              }
          }

        // next subtract the contributions of the lines
        const auto weights_view_line =
          make_array_view(weights.begin(), weights.begin() + 2);
        const auto points_view_line =
          make_array_view(points.begin(), points.begin() + 2);
        for (unsigned int line = 0; line < GeometryInfo<3>::lines_per_cell;
             ++line)
          {
            const double line_point =
              (line < 8 ? chart_point[1 - (line % 4) / 2] : chart_point[2]);
            double my_weight = 0.;
            if (line < 8)
              my_weight = linear_shapes[line % 4] * linear_shapes[4 + line / 4];
            else
              {
                const unsigned int subline = line - 8;
                my_weight =
                  linear_shapes[subline % 2] * linear_shapes[2 + subline / 2];
              }
            my_weight -= weights_lines[line];

            if (std::abs(my_weight) < 1e-13)
              continue;

            const types::manifold_id line_manifold_id =
              cell.line(line)->manifold_id();
            if (line_manifold_id == my_manifold_id ||
                line_manifold_id == numbers::flat_manifold_id)
              {
                weights_vertices[GeometryInfo<3>::line_to_cell_vertices(line,
                                                                        0)] -=
                  my_weight * (1. - line_point);
                weights_vertices[GeometryInfo<3>::line_to_cell_vertices(line,
                                                                        1)] -=
                  my_weight * (line_point);
              }
            else
              {
                points[0] =
                  vertices[GeometryInfo<3>::line_to_cell_vertices(line, 0)];
                points[1] =
                  vertices[GeometryInfo<3>::line_to_cell_vertices(line, 1)];
                weights[0] = 1. - line_point;
                weights[1] = line_point;
                new_point -= my_weight * tria.get_manifold(line_manifold_id)
                                           .get_new_point(points_view_line,
                                                          weights_view_line);
              }
          }

        // finally add the contribution of the
        for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
          new_point += weights_vertices[v] * vertices[v];
      }
    return new_point;
  }
} // namespace



template <int dim, int spacedim>
Point<spacedim>
TransfiniteInterpolationManifold<dim, spacedim>::push_forward(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const Point<dim>                                           &chart_point) const
{
  AssertDimension(cell->level(), level_coarse);

  // check that the point is in the unit cell which is the current chart
  // Tolerance 5e-4 chosen that the method also works with manifolds
  // that have some discretization error like SphericalManifold
  Assert(GeometryInfo<dim>::is_inside_unit_cell(chart_point, 5e-4),
         ExcMessage("chart_point is not in unit interval"));

  return compute_transfinite_interpolation(*cell,
                                           chart_point,
                                           coarse_cell_is_flat[cell->index()]);
}



template <int dim, int spacedim>
DerivativeForm<1, dim, spacedim>
TransfiniteInterpolationManifold<dim, spacedim>::push_forward_gradient(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const Point<dim>                                           &chart_point,
  const Point<spacedim> &pushed_forward_chart_point) const
{
  // compute the derivative with the help of finite differences
  DerivativeForm<1, dim, spacedim> grad;
  for (unsigned int d = 0; d < dim; ++d)
    {
      Point<dim>   modified = chart_point;
      const double step     = chart_point[d] > 0.5 ? -1e-8 : 1e-8;

      // avoid checking outside of the unit interval
      modified[d] += step;
      Tensor<1, spacedim> difference =
        compute_transfinite_interpolation(*cell,
                                          modified,
                                          coarse_cell_is_flat[cell->index()]) -
        pushed_forward_chart_point;
      for (unsigned int e = 0; e < spacedim; ++e)
        grad[e][d] = difference[e] / step;
    }
  return grad;
}



template <int dim, int spacedim>
Point<dim>
TransfiniteInterpolationManifold<dim, spacedim>::pull_back(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const Point<spacedim>                                      &point,
  const Point<dim> &initial_guess) const
{
  Point<dim> outside;
  for (unsigned int d = 0; d < dim; ++d)
    outside[d] = internal::invalid_pull_back_coordinate;

  // project the user-given input to unit cell
  Point<dim> chart_point = cell->reference_cell().closest_point(initial_guess);

  // run quasi-Newton iteration with a combination of finite differences for
  // the exact Jacobian and "Broyden's good method". As opposed to the various
  // mapping implementations, this class does not throw exception upon failure
  // as those are relatively expensive and failure occurs quite regularly in
  // the implementation of the compute_chart_points method.
  Tensor<1, spacedim> residual =
    point -
    compute_transfinite_interpolation(*cell,
                                      chart_point,
                                      coarse_cell_is_flat[cell->index()]);
  const double tolerance = 1e-21 * Utilities::fixed_power<2>(cell->diameter());
  double       residual_norm_square = residual.norm_square();
  DerivativeForm<1, dim, spacedim> inv_grad;
  bool                             must_recompute_jacobian = true;
  for (unsigned int i = 0; i < 100; ++i)
    {
      if (residual_norm_square < tolerance)
        {
          // do a final update of the point with the last available Jacobian
          // information. The residual is close to zero due to the check
          // above, but me might improve some of the last digits by a final
          // Newton-like step with step length 1
          Tensor<1, dim> update;
          for (unsigned int d = 0; d < spacedim; ++d)
            for (unsigned int e = 0; e < dim; ++e)
              update[e] += inv_grad[d][e] * residual[d];
          return chart_point + update;
        }

      // every 9 iterations, including the first time around, we create an
      // approximation of the Jacobian with finite differences. Broyden's
      // method usually does not need more than 5-8 iterations, but sometimes
      // we might have had a bad initial guess and then we can accelerate
      // convergence considerably with getting the actual Jacobian rather than
      // using secant-like methods (one gradient calculation in 3d costs as
      // much as 3 more iterations). this usually happens close to convergence
      // and one more step with the finite-differenced Jacobian leads to
      // convergence
      if (must_recompute_jacobian || i % 9 == 0)
        {
          // if the determinant is zero or negative, the mapping is either not
          // invertible or already has inverted and we are outside the valid
          // chart region. Note that the Jacobian here represents the
          // derivative of the forward map and should have a positive
          // determinant since we use properly oriented meshes.
          DerivativeForm<1, dim, spacedim> grad =
            push_forward_gradient(cell,
                                  chart_point,
                                  Point<spacedim>(point - residual));
          if (grad.determinant() <= 0.0)
            return outside;
          inv_grad                = grad.covariant_form();
          must_recompute_jacobian = false;
        }
      Tensor<1, dim> update;
      for (unsigned int d = 0; d < spacedim; ++d)
        for (unsigned int e = 0; e < dim; ++e)
          update[e] += inv_grad[d][e] * residual[d];

      // Line search, accept step if the residual has decreased
      double alpha = 1.;

      // check if point is inside 1.2 times the unit cell to avoid
      // hitting points very far away from valid ones in the manifolds
      while (
        !GeometryInfo<dim>::is_inside_unit_cell(chart_point + alpha * update,
                                                0.2) &&
        alpha > 1e-7)
        alpha *= 0.5;

      const Tensor<1, spacedim> old_residual = residual;
      while (alpha > 1e-4)
        {
          Point<dim>                guess = chart_point + alpha * update;
          const Tensor<1, spacedim> residual_guess =
            point - compute_transfinite_interpolation(
                      *cell, guess, coarse_cell_is_flat[cell->index()]);
          const double residual_norm_new = residual_guess.norm_square();
          if (residual_norm_new < residual_norm_square)
            {
              residual             = residual_guess;
              residual_norm_square = residual_norm_new;
              chart_point += alpha * update;
              break;
            }
          else
            alpha *= 0.5;
        }
      // If alpha got very small, it is likely due to a bad Jacobian
      // approximation with Broyden's method (relatively far away from the
      // zero), which can be corrected by the outer loop when a Newton update
      // is recomputed. The second case is when the Jacobian is actually bad
      // and we should fail as early as possible. Since we cannot really
      // distinguish the two, we must continue here in any case.
      if (alpha <= 1e-4)
        must_recompute_jacobian = true;

      // update the inverse Jacobian with "Broyden's good method" and
      // Sherman-Morrison formula for the update of the inverse, see
      // https://en.wikipedia.org/wiki/Broyden%27s_method
      // J^{-1}_n = J^{-1}_{n-1} + (delta x_n - J^{-1}_{n-1} delta f_n) /
      // (delta x_n^T J_{-1}_{n-1} delta f_n) delta x_n^T J^{-1}_{n-1}

      // switch sign in residual as compared to the formula above because we
      // use a negative definition of the residual with respect to the
      // Jacobian
      const Tensor<1, spacedim> delta_f = old_residual - residual;

      Tensor<1, dim> Jinv_deltaf;
      for (unsigned int d = 0; d < spacedim; ++d)
        for (unsigned int e = 0; e < dim; ++e)
          Jinv_deltaf[e] += inv_grad[d][e] * delta_f[d];

      const Tensor<1, dim> delta_x = alpha * update;

      // prevent division by zero. This number should be scale-invariant
      // because Jinv_deltaf carries no units and x is in reference
      // coordinates.
      if (std::abs(delta_x * Jinv_deltaf) > 1e-12 && !must_recompute_jacobian)
        {
          const Tensor<1, dim> factor =
            (delta_x - Jinv_deltaf) / (delta_x * Jinv_deltaf);
          Tensor<1, spacedim> jac_update;
          for (unsigned int d = 0; d < spacedim; ++d)
            for (unsigned int e = 0; e < dim; ++e)
              jac_update[d] += delta_x[e] * inv_grad[d][e];
          for (unsigned int d = 0; d < spacedim; ++d)
            for (unsigned int e = 0; e < dim; ++e)
              inv_grad[d][e] += factor[e] * jac_update[d];
        }
    }
  return outside;
}



template <int dim, int spacedim>
std::array<unsigned int, 20>
TransfiniteInterpolationManifold<dim, spacedim>::
  get_possible_cells_around_points(
    const ArrayView<const Point<spacedim>> &points) const
{
  // The methods to identify cells around points in GridTools are all written
  // for the active cells, but we are here looking at some cells at the coarse
  // level.
  Assert(triangulation != nullptr, ExcNotInitialized());
  Assert(triangulation->begin_active()->level() >= level_coarse,
         ExcMessage("The manifold was initialized with level " +
                    std::to_string(level_coarse) + " but there are now" +
                    "active cells on a lower level. Coarsening the mesh is " +
                    "currently not supported"));

  // This computes the distance of the surrounding points transformed to the
  // unit cell from the unit cell.
  typename Triangulation<dim, spacedim>::cell_iterator cell =
                                                         triangulation->begin(
                                                           level_coarse),
                                                       endc =
                                                         triangulation->end(
                                                           level_coarse);
  boost::container::small_vector<std::pair<double, unsigned int>, 200>
    distances_and_cells;
  for (; cell != endc; ++cell)
    {
      // only consider cells where the current manifold is attached
      if (&cell->get_manifold() != this)
        continue;

      std::array<Point<spacedim>, GeometryInfo<dim>::vertices_per_cell>
        vertices;
      for (const unsigned int vertex_n : GeometryInfo<dim>::vertex_indices())
        {
          vertices[vertex_n] = cell->vertex(vertex_n);
        }

      // cheap check: if any of the points is not inside a circle around the
      // center of the loop, we can skip the expensive part below (this assumes
      // that the manifold does not deform the grid too much)
      Point<spacedim> center;
      for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
        center += vertices[v];
      center *= 1. / GeometryInfo<dim>::vertices_per_cell;
      double radius_square = 0.;
      for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
        radius_square =
          std::max(radius_square, (center - vertices[v]).norm_square());
      bool inside_circle = true;
      for (unsigned int i = 0; i < points.size(); ++i)
        if ((center - points[i]).norm_square() > radius_square * 1.5)
          {
            inside_circle = false;
            break;
          }
      if (inside_circle == false)
        continue;

      // slightly more expensive search
      double current_distance = 0;
      for (unsigned int i = 0; i < points.size(); ++i)
        {
          Point<dim> point =
            quadratic_approximation[cell->index()].compute(points[i]);
          current_distance += GeometryInfo<dim>::distance_to_unit_cell(point);
        }
      distances_and_cells.push_back(
        std::make_pair(current_distance, cell->index()));
    }
  // no coarse cell could be found -> transformation failed
  AssertThrow(distances_and_cells.size() > 0,
              (typename Mapping<dim, spacedim>::ExcTransformationFailed()));
  std::sort(distances_and_cells.begin(), distances_and_cells.end());
  std::array<unsigned int, 20> cells;
  cells.fill(numbers::invalid_unsigned_int);
  for (unsigned int i = 0; i < distances_and_cells.size() && i < cells.size();
       ++i)
    cells[i] = distances_and_cells[i].second;

  return cells;
}



template <int dim, int spacedim>
typename Triangulation<dim, spacedim>::cell_iterator
TransfiniteInterpolationManifold<dim, spacedim>::compute_chart_points(
  const ArrayView<const Point<spacedim>> &surrounding_points,
  ArrayView<Point<dim>>                   chart_points) const
{
  Assert(surrounding_points.size() == chart_points.size(),
         ExcMessage("The chart points array view must be as large as the "
                    "surrounding points array view."));

  std::array<unsigned int, 20> nearby_cells =
    get_possible_cells_around_points(surrounding_points);

  // This function is nearly always called to place new points on a cell or
  // cell face. In this case, the general structure of the surrounding points
  // is known (i.e., if there are eight surrounding points, then they will
  // almost surely be either eight points around a quadrilateral or the eight
  // vertices of a cube). Hence, making this assumption, we use two
  // optimizations (one for structdim == 2 and one for structdim == 3) that
  // guess the locations of some of the chart points more efficiently than the
  // affine map approximation. The affine map approximation is used whenever
  // we don't have a cheaper guess available.

  // Function that can guess the location of a chart point by assuming that
  // the eight surrounding points are points on a two-dimensional object
  // (either a cell in 2d or the face of a hexahedron in 3d), arranged like
  //
  //     2 - 7 - 3
  //     |       |
  //     4       5
  //     |       |
  //     0 - 6 - 1
  //
  // This function assumes that the first three chart points have been
  // computed since there is no effective way to guess them.
  auto guess_chart_point_structdim_2 = [&](const unsigned int i) -> Point<dim> {
    Assert(surrounding_points.size() == 8 && 2 < i && i < 8,
           ExcMessage("This function assumes that there are eight surrounding "
                      "points around a two-dimensional object. It also assumes "
                      "that the first three chart points have already been "
                      "computed."));
    switch (i)
      {
        case 0:
        case 1:
        case 2:
          DEAL_II_ASSERT_UNREACHABLE();
          break;
        case 3:
          return chart_points[1] + (chart_points[2] - chart_points[0]);
        case 4:
          return 0.5 * (chart_points[0] + chart_points[2]);
        case 5:
          return 0.5 * (chart_points[1] + chart_points[3]);
        case 6:
          return 0.5 * (chart_points[0] + chart_points[1]);
        case 7:
          return 0.5 * (chart_points[2] + chart_points[3]);
        default:
          DEAL_II_ASSERT_UNREACHABLE();
      }

    return {};
  };

  // Function that can guess the location of a chart point by assuming that
  // the eight surrounding points form the vertices of a hexahedron, arranged
  // like
  //
  //         6-------7
  //        /|      /|
  //       /       / |
  //      /  |    /  |
  //     4-------5   |
  //     |   2- -|- -3
  //     |  /    |  /
  //     |       | /
  //     |/      |/
  //     0-------1
  //
  // (where vertex 2 is the back left vertex) we can estimate where chart
  // points 5 - 7 are by computing the height (in chart coordinates) as c4 -
  // c0 and then adding that onto the appropriate bottom vertex.
  //
  // This function assumes that the first five chart points have been computed
  // since there is no effective way to guess them.
  auto guess_chart_point_structdim_3 = [&](const unsigned int i) -> Point<dim> {
    Assert(surrounding_points.size() == 8 && 4 < i && i < 8,
           ExcMessage("This function assumes that there are eight surrounding "
                      "points around a three-dimensional object. It also "
                      "assumes that the first five chart points have already "
                      "been computed."));
    return chart_points[i - 4] + (chart_points[4] - chart_points[0]);
  };

  // Check if we can use the two chart point shortcuts above before we start:
  bool use_structdim_2_guesses = false;
  bool use_structdim_3_guesses = false;
  // note that in the structdim 2 case: 0 - 6 and 2 - 7 should be roughly
  // parallel, while in the structdim 3 case, 0 - 6 and 2 - 7 should be roughly
  // orthogonal. Use the angle between these two vectors to figure out if we
  // should turn on either structdim optimization.
  if (surrounding_points.size() == 8)
    {
      const Tensor<1, spacedim> v06 =
        surrounding_points[6] - surrounding_points[0];
      const Tensor<1, spacedim> v27 =
        surrounding_points[7] - surrounding_points[2];

      // note that we can save a call to sqrt() by rearranging
      const double cosine = scalar_product(v06, v27) /
                            std::sqrt(v06.norm_square() * v27.norm_square());
      if (0.707 < cosine)
        // the angle is less than pi/4, so these vectors are roughly parallel:
        // enable the structdim 2 optimization
        use_structdim_2_guesses = true;
      else if (spacedim == 3)
        // otherwise these vectors are roughly orthogonal: enable the
        // structdim 3 optimization if we are in 3d
        use_structdim_3_guesses = true;
    }
  // we should enable at most one of the optimizations
  Assert((!use_structdim_2_guesses && !use_structdim_3_guesses) ||
           (use_structdim_2_guesses ^ use_structdim_3_guesses),
         ExcInternalError());



  auto compute_chart_point =
    [&](const typename Triangulation<dim, spacedim>::cell_iterator &cell,
        const unsigned int point_index) {
      Point<dim> guess;
      // an optimization: keep track of whether or not we used the quadratic
      // approximation so that we don't call pull_back with the same
      // initial guess twice (i.e., if pull_back fails the first time,
      // don't try again with the same function arguments).
      bool used_quadratic_approximation = false;
      // if we have already computed three points, we can guess the fourth
      // to be the missing corner point of a rectangle
      if (point_index == 3 && surrounding_points.size() >= 8)
        guess = chart_points[1] + (chart_points[2] - chart_points[0]);
      else if (use_structdim_2_guesses && 3 < point_index)
        guess = guess_chart_point_structdim_2(point_index);
      else if (use_structdim_3_guesses && 4 < point_index)
        guess = guess_chart_point_structdim_3(point_index);
      else if (dim == 3 && point_index > 7 && surrounding_points.size() == 26)
        {
          if (point_index < 20)
            guess =
              0.5 * (chart_points[GeometryInfo<dim>::line_to_cell_vertices(
                       point_index - 8, 0)] +
                     chart_points[GeometryInfo<dim>::line_to_cell_vertices(
                       point_index - 8, 1)]);
          else
            guess =
              0.25 * (chart_points[GeometryInfo<dim>::face_to_cell_vertices(
                        point_index - 20, 0)] +
                      chart_points[GeometryInfo<dim>::face_to_cell_vertices(
                        point_index - 20, 1)] +
                      chart_points[GeometryInfo<dim>::face_to_cell_vertices(
                        point_index - 20, 2)] +
                      chart_points[GeometryInfo<dim>::face_to_cell_vertices(
                        point_index - 20, 3)]);
        }
      else
        {
          guess = quadratic_approximation[cell->index()].compute(
            surrounding_points[point_index]);
          used_quadratic_approximation = true;
        }
      chart_points[point_index] =
        pull_back(cell, surrounding_points[point_index], guess);

      // the initial guess may not have been good enough: if applicable,
      // try again with the affine approximation (which is more accurate
      // than the cheap methods used above)
      if (chart_points[point_index][0] ==
            internal::invalid_pull_back_coordinate &&
          !used_quadratic_approximation)
        {
          guess = quadratic_approximation[cell->index()].compute(
            surrounding_points[point_index]);
          chart_points[point_index] =
            pull_back(cell, surrounding_points[point_index], guess);
        }

      if (chart_points[point_index][0] ==
          internal::invalid_pull_back_coordinate)
        {
          for (unsigned int d = 0; d < dim; ++d)
            guess[d] = 0.5;
          chart_points[point_index] =
            pull_back(cell, surrounding_points[point_index], guess);
        }
    };

  // check whether all points are inside the unit cell of the current chart
  for (unsigned int c = 0; c < nearby_cells.size(); ++c)
    {
      typename Triangulation<dim, spacedim>::cell_iterator cell(
        triangulation, level_coarse, nearby_cells[c]);
      bool inside_unit_cell = true;
      for (unsigned int i = 0; i < surrounding_points.size(); ++i)
        {
          compute_chart_point(cell, i);

          // Tolerance 5e-4 chosen that the method also works with manifolds
          // that have some discretization error like SphericalManifold
          if (GeometryInfo<dim>::is_inside_unit_cell(chart_points[i], 5e-4) ==
              false)
            {
              inside_unit_cell = false;
              break;
            }
        }
      if (inside_unit_cell == true)
        {
          return cell;
        }

      // if we did not find a point and this was the last valid cell (the next
      // iterate being the end of the array or an invalid tag), we must stop
      if (c == nearby_cells.size() - 1 ||
          nearby_cells[c + 1] == numbers::invalid_unsigned_int)
        {
          // generate additional information to help debugging why we did not
          // get a point
          std::ostringstream message;
          for (unsigned int b = 0; b <= c; ++b)
            {
              typename Triangulation<dim, spacedim>::cell_iterator cell(
                triangulation, level_coarse, nearby_cells[b]);
              message << "Looking at cell " << cell->id()
                      << " with vertices: " << std::endl;
              for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
                message << cell->vertex(v) << "    ";
              message << std::endl;
              message << "Transformation to chart coordinates: " << std::endl;
              for (unsigned int i = 0; i < surrounding_points.size(); ++i)
                {
                  compute_chart_point(cell, i);
                  message << surrounding_points[i] << " -> " << chart_points[i]
                          << std::endl;
                }
            }

          AssertThrow(false,
                      (typename Mapping<dim, spacedim>::ExcTransformationFailed(
                        message.str())));
        }
    }

  // a valid inversion should have returned a point above. an invalid
  // inversion should have triggered the assertion, so we should never end up
  // here
  DEAL_II_ASSERT_UNREACHABLE();
  return typename Triangulation<dim, spacedim>::cell_iterator();
}



template <int dim, int spacedim>
Point<spacedim>
TransfiniteInterpolationManifold<dim, spacedim>::get_new_point(
  const ArrayView<const Point<spacedim>> &surrounding_points,
  const ArrayView<const double>          &weights) const
{
  boost::container::small_vector<Point<dim>, 100> chart_points(
    surrounding_points.size());
  ArrayView<Point<dim>> chart_points_view =
    make_array_view(chart_points.begin(), chart_points.end());
  const auto cell = compute_chart_points(surrounding_points, chart_points_view);

  const Point<dim> p_chart =
    chart_manifold.get_new_point(chart_points_view, weights);

  return push_forward(cell, p_chart);
}



template <int dim, int spacedim>
void
TransfiniteInterpolationManifold<dim, spacedim>::get_new_points(
  const ArrayView<const Point<spacedim>> &surrounding_points,
  const Table<2, double>                 &weights,
  ArrayView<Point<spacedim>>              new_points) const
{
  Assert(weights.size(0) > 0, ExcEmptyObject());
  AssertDimension(surrounding_points.size(), weights.size(1));

  boost::container::small_vector<Point<dim>, 100> chart_points(
    surrounding_points.size());
  ArrayView<Point<dim>> chart_points_view =
    make_array_view(chart_points.begin(), chart_points.end());
  const auto cell = compute_chart_points(surrounding_points, chart_points_view);

  boost::container::small_vector<Point<dim>, 100> new_points_on_chart(
    weights.size(0));
  chart_manifold.get_new_points(chart_points_view,
                                weights,
                                make_array_view(new_points_on_chart.begin(),
                                                new_points_on_chart.end()));

  for (unsigned int row = 0; row < weights.size(0); ++row)
    new_points[row] = push_forward(cell, new_points_on_chart[row]);
}



// explicit instantiations
#include "grid/manifold_lib.inst"

DEAL_II_NAMESPACE_CLOSE
