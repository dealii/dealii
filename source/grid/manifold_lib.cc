// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include <deal.II/base/std_cxx14/memory.h>
#include <deal.II/base/table.h>
#include <deal.II/base/tensor.h>

#include <deal.II/fe/mapping.h>

#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <cmath>

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
        const Tensor<1, 3> tmp     = cos(theta) * u + sin(theta) * dirUnit;
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

template <int dim, int spacedim>
PolarManifold<dim, spacedim>::PolarManifold(const Point<spacedim> center) :
  ChartManifold<dim, spacedim, spacedim>(
    PolarManifold<dim, spacedim>::get_periodicity()),
  center(center)
{}



template <int dim, int spacedim>
std::unique_ptr<Manifold<dim, spacedim>>
PolarManifold<dim, spacedim>::clone() const
{
  return std_cxx14::make_unique<PolarManifold<dim, spacedim>>(center);
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
          p[0] = rho * cos(theta);
          p[1] = rho * sin(theta);
          break;
        case 3:
          {
            const double phi = spherical_point[2];
            p[0]             = rho * sin(theta) * cos(phi);
            p[1]             = rho * sin(theta) * sin(phi);
            p[2]             = rho * cos(theta);
            break;
          }
        default:
          Assert(false, ExcNotImplemented());
      }
  return p + center;
}



template <int dim, int spacedim>
Point<spacedim>
PolarManifold<dim, spacedim>::pull_back(
  const Point<spacedim> &space_point) const
{
  const Tensor<1, spacedim> R   = space_point - center;
  const double              rho = R.norm();

  Point<spacedim> p;
  p[0] = rho;

  switch (spacedim)
    {
      case 2:
        {
          p[1] = atan2(R[1], R[0]);
          if (p[1] < 0)
            p[1] += 2 * numbers::PI;
          break;
        }

      case 3:
        {
          const double z = R[2];
          p[2]           = atan2(R[1], R[0]); // phi
          if (p[2] < 0)
            p[2] += 2 * numbers::PI;                        // phi is periodic
          p[1] = atan2(sqrt(R[0] * R[0] + R[1] * R[1]), z); // theta
          break;
        }

      default:
        Assert(false, ExcNotImplemented());
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
            DX[0][0] = cos(theta);
            DX[0][1] = -rho * sin(theta);
            DX[1][0] = sin(theta);
            DX[1][1] = rho * cos(theta);
            break;
          }

        case 3:
          {
            const double phi = spherical_point[2];
            DX[0][0]         = sin(theta) * cos(phi);
            DX[0][1]         = rho * cos(theta) * cos(phi);
            DX[0][2]         = -rho * sin(theta) * sin(phi);

            DX[1][0] = sin(theta) * sin(phi);
            DX[1][1] = rho * cos(theta) * sin(phi);
            DX[1][2] = rho * sin(theta) * cos(phi);

            DX[2][0] = cos(theta);
            DX[2][1] = -rho * sin(theta);
            DX[2][2] = 0;
            break;
          }

        default:
          Assert(false, ExcNotImplemented());
      }
  return DX;
}



// ============================================================
// SphericalManifold
// ============================================================

template <int dim, int spacedim>
SphericalManifold<dim, spacedim>::SphericalManifold(
  const Point<spacedim> center) :
  center(center),
  polar_manifold(center)
{}



template <int dim, int spacedim>
std::unique_ptr<Manifold<dim, spacedim>>
SphericalManifold<dim, spacedim>::clone() const
{
  return std_cxx14::make_unique<SphericalManifold<dim, spacedim>>(center);
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

  const Tensor<1, spacedim> v1 = p1 - center;
  const Tensor<1, spacedim> v2 = p2 - center;
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
    return center;

  // Points are along a line, in which case e1 and e2 are essentially the same.
  if (cosgamma > 1 - 8. * std::numeric_limits<double>::epsilon())
    return Point<spacedim>(center + w * v2 + (1 - w) * v1);

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
  return Point<spacedim>(center + (w * r2 + (1.0 - w) * r1) * P);
}



template <int dim, int spacedim>
Tensor<1, spacedim>
SphericalManifold<dim, spacedim>::get_tangent_vector(
  const Point<spacedim> &p1,
  const Point<spacedim> &p2) const
{
  const double tol = 1e-10;
  (void)tol;

  Assert(p1 != p2, ExcMessage("p1 and p2 should not concide."));

  const Tensor<1, spacedim> v1 = p1 - center;
  const Tensor<1, spacedim> v2 = p2 - center;
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
  const Point<spacedim> &                                     p) const
{
  // if the maximum deviation for the distance from the vertices to the center
  // is less than 1.e-5 of the minimum distance to the first vertex, assume we
  // can simply return p-center. otherwise, we compute the normal using
  // get_normal_vector
  constexpr unsigned int n_vertices = GeometryInfo<spacedim>::vertices_per_face;
  std::array<double, n_vertices>     distances_to_center;
  std::array<double, n_vertices - 1> distances_to_first_vertex;
  distances_to_center[0] = (face->vertex(0) - center).norm_square();
  for (unsigned int i = 1; i < n_vertices; ++i)
    {
      distances_to_center[i] = (face->vertex(i) - center).norm_square();
      distances_to_first_vertex[i - 1] =
        (face->vertex(i) - face->vertex(0)).norm_square();
    }
  const auto minmax_distance =
    std::minmax_element(distances_to_center.begin(), distances_to_center.end());
  const auto min_distance_to_first_vertex = std::min_element(
    distances_to_first_vertex.begin(), distances_to_first_vertex.end());

  if (*minmax_distance.second - *minmax_distance.first <
      1.e-10 * *min_distance_to_first_vertex)
    {
      const Tensor<1, spacedim> unnormalized_spherical_normal = p - center;
      const Tensor<1, spacedim> normalized_spherical_normal =
        unnormalized_spherical_normal / unnormalized_spherical_normal.norm();
      return normalized_spherical_normal;
    }
  return Manifold<dim, spacedim>::normal_vector(face, p);
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
  // if the maximum deviation for the distance from the vertices to the center
  // is less than 1.e-5 of the minimum distance to the first vertex, assume we
  // can simply return vertex-center. otherwise, we compute the normal using
  // get_normal_vector
  constexpr unsigned int n_vertices = GeometryInfo<spacedim>::vertices_per_face;
  std::array<double, n_vertices>     distances_to_center;
  std::array<double, n_vertices - 1> distances_to_first_vertex;
  distances_to_center[0] = (face->vertex(0) - center).norm_square();
  for (unsigned int i = 1; i < n_vertices; ++i)
    {
      distances_to_center[i] = (face->vertex(i) - center).norm_square();
      distances_to_first_vertex[i - 1] =
        (face->vertex(i) - face->vertex(0)).norm_square();
    }
  const auto minmax_distance =
    std::minmax_element(distances_to_center.begin(), distances_to_center.end());
  const auto min_distance_to_first_vertex = std::min_element(
    distances_to_first_vertex.begin(), distances_to_first_vertex.end());

  if (*minmax_distance.second - *minmax_distance.first <
      1.e-10 * *min_distance_to_first_vertex)
    {
      for (unsigned int vertex = 0; vertex < n_vertices; ++vertex)
        face_vertex_normals[vertex] = face->vertex(vertex) - center;
    }
  else
    Manifold<dim, spacedim>::get_normals_at_vertices(face, face_vertex_normals);
}



template <int dim, int spacedim>
void
SphericalManifold<dim, spacedim>::get_new_points(
  const ArrayView<const Point<spacedim>> &surrounding_points,
  const Table<2, double> &                weights,
  ArrayView<Point<spacedim>>              new_points) const
{
  AssertDimension(new_points.size(), weights.size(0));
  AssertDimension(surrounding_points.size(), weights.size(1));

  get_new_points(surrounding_points, make_array_view(weights), new_points);

  return;
}



template <int dim, int spacedim>
Point<spacedim>
SphericalManifold<dim, spacedim>::get_new_point(
  const ArrayView<const Point<spacedim>> &vertices,
  const ArrayView<const double> &         weights) const
{
  // To avoid duplicating all of the logic in get_new_points, simply call it
  // for one position.
  Point<spacedim> new_point;
  get_new_points(
    vertices, weights, make_array_view(&new_point, &new_point + 1));

  return new_point;
}



template <int dim, int spacedim>
void
SphericalManifold<dim, spacedim>::get_new_points(
  const ArrayView<const Point<spacedim>> &surrounding_points,
  const ArrayView<const double> &         weights,
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
          get_intermediate_point(surrounding_points[0],
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
      directions[i] = surrounding_points[i] - center;
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
          new_points[row]               = center;
          accurate_point_was_found[row] = true;
          continue;
        }

      // If not in 3D, just use the implementation from PolarManifold
      // after we verified that the candidate is not the center.
      if (spacedim < 3)
        new_points[row] = polar_manifold.get_new_point(
          surrounding_points,
          ArrayView<const double>(&weights[row * weight_columns],
                                  weight_columns));
    }

  // In this case, we treated the case that the candidate is the center and
  // obtained the new locations from the PolarManifold object otherwise.
  if (spacedim < 3)
    return;

  // If all the points are close to each other, we expect the estimate to
  // be good enough. This tolerance was chosen such that the first iteration
  // for a at least three time refined HyperShell mesh with radii .5 and 1.
  // doesn't already succeed.
  if (max_distance < 2e-2)
    {
      for (unsigned int row = 0; row < weight_rows; ++row)
        new_points[row] =
          center + new_candidates[row].first * new_candidates[row].second;

      return;
    }

  // Step 2:
  // Do more expensive Newton-style iterations to improve the estimate.

  // Search for duplicate directions and merge them to minimize the cost of
  // the get_new_point function call below.
  boost::container::small_vector<double, 1000> merged_weights(weights.size());
  boost::container::small_vector<Tensor<1, spacedim>, 100> merged_directions(
    surrounding_points.size(), Point<spacedim>());
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
      for (unsigned int existing_row = 0; existing_row < row; ++existing_row)
        {
          bool identical_weights = true;

          for (unsigned int weight_index = 0;
               weight_index < n_unique_directions;
               ++weight_index)
            if (std::abs(merged_weights[row * weight_columns + weight_index] -
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
  const ArrayView<const double> array_merged_distances = make_array_view(
    merged_distances.begin(), merged_distances.begin() + n_unique_directions);

  for (unsigned int row = 0; row < weight_rows; ++row)
    if (!accurate_point_was_found[row])
      {
        if (merged_weights_index[row] == numbers::invalid_unsigned_int)
          {
            const ArrayView<const double> array_merged_weights(
              &merged_weights[row * weight_columns], n_unique_directions);
            new_candidates[row].second =
              get_new_point(array_merged_directions,
                            array_merged_distances,
                            array_merged_weights,
                            Point<spacedim>(new_candidates[row].second));
          }
        else
          new_candidates[row].second =
            new_candidates[merged_weights_index[row]].second;

        new_points[row] =
          center + new_candidates[row].first * new_candidates[row].second;
      }
}



template <int dim, int spacedim>
std::pair<double, Tensor<1, spacedim>>
SphericalManifold<dim, spacedim>::guess_new_point(
  const ArrayView<const Tensor<1, spacedim>> &directions,
  const ArrayView<const double> &             distances,
  const ArrayView<const double> &             weights) const
{
  const double        tolerance = 1e-10;
  double              rho       = 0.;
  Tensor<1, spacedim> candidate;

  // Perform a simple average ...
  double total_weights = 0.;
  for (unsigned int i = 0; i < directions.size(); i++)
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


namespace
{
  template <int spacedim>
  Point<spacedim>
  do_get_new_point(const ArrayView<const Tensor<1, spacedim>> & /*directions*/,
                   const ArrayView<const double> & /*distances*/,
                   const ArrayView<const double> & /*weights*/,
                   const Point<spacedim> & /*candidate_point*/)
  {
    Assert(false, ExcNotImplemented());
    return Point<spacedim>();
  }

  template <>
  Point<3>
  do_get_new_point(const ArrayView<const Tensor<1, 3>> &directions,
                   const ArrayView<const double> &      distances,
                   const ArrayView<const double> &      weights,
                   const Point<3> &                     candidate_point)
  {
    (void)distances;

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
          SphericalManifold<3, 3> unit_manifold;
          Assert(std::abs(weights[0] + weights[1] - 1.0) < 1e-13,
                 ExcMessage("Weights do not sum up to 1"));
          Point<3> intermediate = unit_manifold.get_intermediate_point(
            Point<3>(directions[0]), Point<3>(directions[1]), weights[1]);
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

          // For each vertices vector, compute the tangent vector from candidate
          // towards the vertices vector -- its length is the spherical length
          // from candidate to the vertices vector.
          // Then compute its contribution to the Hessian.
          gradient = 0.;
          Hessian  = 0.;
          for (unsigned int i = 0; i < n_merged_points; ++i)
            if (std::abs(weights[i]) > 1.e-15)
              {
                vPerp = internal::projected_direction(directions[i], candidate);
                const double sinthetaSq = vPerp.norm_square();
                const double sintheta   = std::sqrt(sinthetaSq);
                if (sintheta < tolerance)
                  {
                    Hessian[0][0] += weights[i];
                    Hessian[1][1] += weights[i];
                  }
                else
                  {
                    const double costheta     = (directions[i]) * candidate;
                    const double theta        = atan2(sintheta, costheta);
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
                    const double offdiag  = cosphi * sinphi * wt * (1.0 - tt);
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

          // Step 2b: rotate candidate in direction xDisp for a new candidate.
          const Point<3> candidateOld = candidate;
          candidate =
            Point<3>(internal::apply_exponential_map(candidate, xDisp));

          // Step 2c: return the new candidate if we didn't move
          if ((candidate - candidateOld).norm_square() < tolerance * tolerance)
            break;
        }
    }
    return candidate;
  }
} // namespace



template <int dim, int spacedim>
Point<spacedim>
SphericalManifold<dim, spacedim>::get_new_point(
  const ArrayView<const Tensor<1, spacedim>> &,
  const ArrayView<const double> &,
  const ArrayView<const double> &,
  const Point<spacedim> &) const
{
  Assert(false, ExcNotImplemented());
  return Point<spacedim>();
}



template <>
Point<3>
SphericalManifold<1, 3>::get_new_point(
  const ArrayView<const Tensor<1, 3>> &directions,
  const ArrayView<const double> &      distances,
  const ArrayView<const double> &      weights,
  const Point<3> &                     candidate_point) const
{
  return do_get_new_point(directions, distances, weights, candidate_point);
}



template <>
Point<3>
SphericalManifold<2, 3>::get_new_point(
  const ArrayView<const Tensor<1, 3>> &directions,
  const ArrayView<const double> &      distances,
  const ArrayView<const double> &      weights,
  const Point<3> &                     candidate_point) const
{
  return do_get_new_point(directions, distances, weights, candidate_point);
}



template <>
Point<3>
SphericalManifold<3, 3>::get_new_point(
  const ArrayView<const Tensor<1, 3>> &directions,
  const ArrayView<const double> &      distances,
  const ArrayView<const double> &      weights,
  const Point<3> &                     candidate_point) const
{
  return do_get_new_point(directions, distances, weights, candidate_point);
}



// ============================================================
// CylindricalManifold
// ============================================================
template <int dim, int spacedim>
CylindricalManifold<dim, spacedim>::CylindricalManifold(
  const unsigned int axis,
  const double       tolerance) :
  CylindricalManifold<dim, spacedim>(Point<spacedim>::unit_vector(axis),
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
  const Point<spacedim> &    point_on_axis,
  const double               tolerance) :
  ChartManifold<dim, spacedim, 3>(Tensor<1, 3>({0, 2. * numbers::PI, 0})),
  normal_direction(internal::compute_normal(direction, true)),
  direction(direction / direction.norm()),
  point_on_axis(point_on_axis),
  tolerance(tolerance)
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
  return std_cxx14::make_unique<CylindricalManifold<dim, spacedim>>(
    direction, point_on_axis, tolerance);
}



template <int dim, int spacedim>
Point<spacedim>
CylindricalManifold<dim, spacedim>::get_new_point(
  const ArrayView<const Point<spacedim>> &surrounding_points,
  const ArrayView<const double> &         weights) const
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
    return Point<spacedim>() + direction * lambda;
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
  const double dot = normal_direction * p_diff;
  const double det = direction * cross_product_3d(normal_direction, p_diff);
  const double phi = std::atan2(det, dot);

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

  // Rotate the orthogonal direction by the given angle.
  // Formula from Section 5.2 in
  // http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/
  // simplified assuming normal_direction and direction are orthogonal
  // and unit vectors.
  const double sine_r           = std::sin(chart_point(1)) * chart_point(0);
  const double cosine_r         = std::cos(chart_point(1)) * chart_point(0);
  const Tensor<1, spacedim> dxn = cross_product_3d(direction, normal_direction);
  const Tensor<1, spacedim> intermediate =
    normal_direction * cosine_r + dxn * sine_r;

  // Finally, put everything together.
  return point_on_axis + direction * chart_point(2) + intermediate;
}



template <int dim, int spacedim>
DerivativeForm<1, 3, spacedim>
CylindricalManifold<dim, spacedim>::push_forward_gradient(
  const Point<3> &chart_point) const
{
  Assert(spacedim == 3,
         ExcMessage("CylindricalManifold can only be used for spacedim==3!"));

  Tensor<2, 3> derivatives;

  // Rotate the orthogonal direction by the given angle.
  // Formula from Section 5.2 in
  // http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/
  // simplified assuming normal_direction and direction are orthogonal
  // and unit vectors.
  const double              sine   = std::sin(chart_point(1));
  const double              cosine = std::cos(chart_point(1));
  const Tensor<1, spacedim> dxn = cross_product_3d(direction, normal_direction);
  const Tensor<1, spacedim> intermediate =
    normal_direction * cosine + dxn * sine;

  // derivative w.r.t the radius
  derivatives[0][0] = intermediate[0];
  derivatives[1][0] = intermediate[1];
  derivatives[2][0] = intermediate[2];

  // derivatives w.r.t the angle
  derivatives[0][1] = -normal_direction[0] * sine + dxn[0] * cosine;
  derivatives[1][1] = -normal_direction[1] * sine + dxn[1] * cosine;
  derivatives[2][1] = -normal_direction[2] * sine + dxn[2] * cosine;

  // derivatives w.r.t the direction of the axis
  derivatives[0][2] = direction[0];
  derivatives[1][2] = direction[1];
  derivatives[2][2] = direction[2];

  return derivatives;
}



// ============================================================
// FunctionManifold
// ============================================================
template <int dim, int spacedim, int chartdim>
FunctionManifold<dim, spacedim, chartdim>::FunctionManifold(
  const Function<chartdim> & push_forward_function,
  const Function<spacedim> & pull_back_function,
  const Tensor<1, chartdim> &periodicity,
  const double               tolerance) :
  ChartManifold<dim, spacedim, chartdim>(periodicity),
  push_forward_function(&push_forward_function),
  pull_back_function(&pull_back_function),
  tolerance(tolerance),
  owns_pointers(false),
  finite_difference_step(0)
{
  AssertDimension(push_forward_function.n_components, spacedim);
  AssertDimension(pull_back_function.n_components, chartdim);
}



template <int dim, int spacedim, int chartdim>
FunctionManifold<dim, spacedim, chartdim>::FunctionManifold(
  const std::string                                 push_forward_expression,
  const std::string                                 pull_back_expression,
  const Tensor<1, chartdim> &                       periodicity,
  const typename FunctionParser<spacedim>::ConstMap const_map,
  const std::string                                 chart_vars,
  const std::string                                 space_vars,
  const double                                      tolerance,
  const double                                      h) :
  ChartManifold<dim, spacedim, chartdim>(periodicity),
  const_map(const_map),
  tolerance(tolerance),
  owns_pointers(true),
  push_forward_expression(push_forward_expression),
  pull_back_expression(pull_back_expression),
  chart_vars(chart_vars),
  space_vars(space_vars),
  finite_difference_step(h)
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
  // deleted. We need to make sure that our cloned object is constructed in the
  // same way this class was constructed, and that its internal Function
  // pointers point either to the same Function objects used to construct this
  // function (owns_pointers == false) or that the newly generated manifold
  // creates internally the push_forward and pull_back functions using the same
  // expressions that were used to construct this class (own_pointers == true).
  if (owns_pointers == true)
    {
      return std_cxx14::make_unique<FunctionManifold<dim, spacedim, chartdim>>(
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
    return std_cxx14::make_unique<FunctionManifold<dim, spacedim, chartdim>>(
      *push_forward_function,
      *pull_back_function,
      this->get_periodicity(),
      tolerance);
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

#ifdef DEBUG
  Vector<double> pb(chartdim);
  pull_back_function->vector_value(result, pb);
  for (unsigned int i = 0; i < chartdim; ++i)
    Assert(
      (chart_point.norm() > tolerance &&
       (std::abs(pb[i] - chart_point[i]) < tolerance * chart_point.norm())) ||
        (std::abs(pb[i] - chart_point[i]) < tolerance),
      ExcMessage(
        "The push forward is not the inverse of the pull back! Bailing out."));
#endif

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
  double x     = p(0);
  double z     = p(1);
  double y     = p(2);
  double phi   = atan2(y, x);
  double theta = atan2(z, std::sqrt(x * x + y * y) - R);
  double w =
    std::sqrt(pow(y - sin(phi) * R, 2.0) + pow(x - cos(phi) * R, 2.0) + z * z) /
    r;
  return Point<3>(phi, theta, w);
}



template <int dim>
Point<3>
TorusManifold<dim>::push_forward(const Point<3> &chart_point) const
{
  double phi   = chart_point(0);
  double theta = chart_point(1);
  double w     = chart_point(2);

  return Point<3>(cos(phi) * R + r * w * cos(theta) * cos(phi),
                  r * w * sin(theta),
                  sin(phi) * R + r * w * cos(theta) * sin(phi));
}



template <int dim>
TorusManifold<dim>::TorusManifold(const double R, const double r) :
  ChartManifold<dim, 3, 3>(Point<3>(2 * numbers::PI, 2 * numbers::PI, 0.0)),
  r(r),
  R(R)
{
  Assert(R > r,
         ExcMessage("Outer radius R must be greater than the inner "
                    "radius r."));
  Assert(r > 0.0, ExcMessage("inner radius must be positive."));
}



template <int dim>
std::unique_ptr<Manifold<dim, 3>>
TorusManifold<dim>::clone() const
{
  return std_cxx14::make_unique<TorusManifold<dim>>(R, r);
}



template <int dim>
DerivativeForm<1, 3, 3>
TorusManifold<dim>::push_forward_gradient(const Point<3> &chart_point) const
{
  DerivativeForm<1, spacedim, spacedim> DX;

  double phi   = chart_point(0);
  double theta = chart_point(1);
  double w     = chart_point(2);

  DX[0][0] = -sin(phi) * R - r * w * cos(theta) * sin(phi);
  DX[0][1] = -r * w * sin(theta) * cos(phi);
  DX[0][2] = r * cos(theta) * cos(phi);

  DX[1][0] = 0;
  DX[1][1] = r * w * cos(theta);
  DX[1][2] = r * sin(theta);

  DX[2][0] = cos(phi) * R + r * w * cos(theta) * cos(phi);
  DX[2][1] = -r * w * sin(theta) * sin(phi);
  DX[2][2] = r * cos(theta) * sin(phi);

  return DX;
}



// ============================================================
// TransfiniteInterpolationManifold
// ============================================================
template <int dim, int spacedim>
TransfiniteInterpolationManifold<dim,
                                 spacedim>::TransfiniteInterpolationManifold() :
  triangulation(nullptr),
  level_coarse(-1)
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
  // in case the triangulatoin is cleared, remove the pointers by a signal
  clear_signal = triangulation.signals.clear.connect([&]() -> void {
    this->triangulation = nullptr;
    this->level_coarse  = -1;
  });
  level_coarse = triangulation.last()->level();
  coarse_cell_is_flat.resize(triangulation.n_cells(level_coarse), false);
  typename Triangulation<dim, spacedim>::active_cell_iterator
    cell = triangulation.begin(level_coarse),
    endc = triangulation.end(level_coarse);
  for (; cell != endc; ++cell)
    {
      bool cell_is_flat = true;
      for (unsigned int l = 0; l < GeometryInfo<dim>::lines_per_cell; ++l)
        if (cell->line(l)->manifold_id() != cell->manifold_id() &&
            cell->line(l)->manifold_id() != numbers::invalid_manifold_id)
          cell_is_flat = false;
      if (dim > 2)
        for (unsigned int q = 0; q < GeometryInfo<dim>::quads_per_cell; ++q)
          if (cell->quad(q)->manifold_id() != cell->manifold_id() &&
              cell->quad(q)->manifold_id() != numbers::invalid_manifold_id)
            cell_is_flat = false;
      AssertIndexRange(static_cast<unsigned int>(cell->index()),
                       coarse_cell_is_flat.size());
      coarse_cell_is_flat[cell->index()] = cell_is_flat;
    }
}



namespace
{
  // version for 1D
  template <typename AccessorType>
  Point<AccessorType::space_dimension>
  compute_transfinite_interpolation(const AccessorType &cell,
                                    const Point<1> &    chart_point,
                                    const bool /*cell_is_flat*/)
  {
    return cell.vertex(0) * (1. - chart_point[0]) +
           cell.vertex(1) * chart_point[0];
  }

  // version for 2D
  template <typename AccessorType>
  Point<AccessorType::space_dimension>
  compute_transfinite_interpolation(const AccessorType &cell,
                                    const Point<2> &    chart_point,
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
      for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_cell; ++v)
        new_point += weights_vertices[v] * vertices[v];
    else
      {
        // The second line in the formula tells us to subtract the
        // contribution of the vertices.  If a line employs the same manifold
        // as the cell, we can merge the weights of the line with the weights
        // of the vertex with a negative sign while going through the faces
        // (this is a bit artificial in 2D but it becomes clear in 3D where we
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
                line_manifold_id == numbers::invalid_manifold_id)
              {
                weights_vertices[GeometryInfo<2>::line_to_cell_vertices(
                  line, 0)] -= my_weight * (1. - line_point);
                weights_vertices[GeometryInfo<2>::line_to_cell_vertices(
                  line, 1)] -= my_weight * line_point;
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
        for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_cell; ++v)
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

  // version for 3D
  template <typename AccessorType>
  Point<AccessorType::space_dimension>
  compute_transfinite_interpolation(const AccessorType &cell,
                                    const Point<3> &    chart_point,
                                    const bool          cell_is_flat)
  {
    const unsigned int       dim             = AccessorType::dimension;
    const unsigned int       spacedim        = AccessorType::space_dimension;
    const types::manifold_id my_manifold_id  = cell.manifold_id();
    const Triangulation<dim, spacedim> &tria = cell.get_triangulation();

    // Same approach as in 2D, but adding the faces, subtracting the edges, and
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

        double weights_lines[GeometryInfo<3>::lines_per_cell];
        for (unsigned int line = 0; line < GeometryInfo<3>::lines_per_cell;
             ++line)
          weights_lines[line] = 0;

        // start with the contributions of the faces
        std::array<double, GeometryInfo<2>::vertices_per_cell>          weights;
        std::array<Point<spacedim>, GeometryInfo<2>::vertices_per_cell> points;
        // note that the views are immutable, but the arrays are not
        const auto weights_view =
          make_array_view(weights.begin(), weights.end());
        const auto points_view = make_array_view(points.begin(), points.end());

        for (unsigned int face = 0; face < GeometryInfo<3>::faces_per_cell;
             ++face)
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
                face_manifold_id == numbers::invalid_manifold_id)
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
                for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_cell;
                     ++v)
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
                line_manifold_id == numbers::invalid_manifold_id)
              {
                weights_vertices[GeometryInfo<3>::line_to_cell_vertices(
                  line, 0)] -= my_weight * (1. - line_point);
                weights_vertices[GeometryInfo<3>::line_to_cell_vertices(
                  line, 1)] -= my_weight * (line_point);
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
        for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
          new_point += weights_vertices[v] * vertices[v];
      }
    return new_point;
  }
} // namespace



template <int dim, int spacedim>
Point<spacedim>
TransfiniteInterpolationManifold<dim, spacedim>::push_forward(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const Point<dim> &                                          chart_point) const
{
  AssertDimension(cell->level(), level_coarse);

  // check that the point is in the unit cell which is the current chart
  // Tolerance 1e-6 chosen that the method also works with
  // SphericalManifold
  Assert(GeometryInfo<dim>::is_inside_unit_cell(chart_point, 1e-6),
         ExcMessage("chart_point is not in unit interval"));

  return compute_transfinite_interpolation(
    *cell, chart_point, coarse_cell_is_flat[cell->index()]);
}



template <int dim, int spacedim>
DerivativeForm<1, dim, spacedim>
TransfiniteInterpolationManifold<dim, spacedim>::push_forward_gradient(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const Point<dim> &                                          chart_point,
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
        compute_transfinite_interpolation(
          *cell, modified, coarse_cell_is_flat[cell->index()]) -
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
  const Point<spacedim> &                                     point,
  const Point<dim> &initial_guess) const
{
  Point<dim> outside;
  for (unsigned int d = 0; d < dim; ++d)
    outside[d] = internal::invalid_pull_back_coordinate;

  // project the user-given input to unit cell
  Point<dim> chart_point =
    GeometryInfo<dim>::project_to_unit_cell(initial_guess);

  // run quasi-Newton iteration with a combination of finite differences for
  // the exact Jacobian and "Broyden's good method". As opposed to the various
  // mapping implementations, this class does not throw exception upon failure
  // as those are relatively expensive and failure occurs quite regularly in
  // the implementation of the compute_chart_points method.
  Tensor<1, spacedim> residual =
    point - compute_transfinite_interpolation(
              *cell, chart_point, coarse_cell_is_flat[cell->index()]);
  const double tolerance = 1e-21 * Utilities::fixed_power<2>(cell->diameter());
  double       residual_norm_square = residual.norm_square();
  DerivativeForm<1, dim, spacedim> inv_grad;
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

      // every 8 iterations, including the first time around, we create an
      // approximation of the Jacobian with finite differences. Broyden's
      // method usually does not need more than 5-8 iterations, but sometimes
      // we might have had a bad initial guess and then we can accelerate
      // convergence considerably with getting the actual Jacobian rather than
      // using secant-like methods (one gradient calculation in 3D costs as
      // much as 3 more iterations). this usually happens close to convergence
      // and one more step with the finite-differenced Jacobian leads to
      // convergence
      if (i % 8 == 0)
        {
          // if the determinant is zero or negative, the mapping is either not
          // invertible or already has inverted and we are outside the valid
          // chart region. Note that the Jacobian here represents the
          // derivative of the forward map and should have a positive
          // determinant since we use properly oriented meshes.
          DerivativeForm<1, dim, spacedim> grad = push_forward_gradient(
            cell, chart_point, Point<spacedim>(point - residual));
          if (grad.determinant() <= 0.0)
            return outside;
          inv_grad = grad.covariant_form();
        }
      Tensor<1, dim> update;
      for (unsigned int d = 0; d < spacedim; ++d)
        for (unsigned int e = 0; e < dim; ++e)
          update[e] += inv_grad[d][e] * residual[d];

      // Line search, accept step if the residual has decreased
      double alpha = 1.;

      // check if point is inside 1.2 times the unit cell to avoid
      // hitting points very far away from valid ones in the manifolds
      while (!GeometryInfo<dim>::is_inside_unit_cell(
               chart_point + alpha * update, 0.2) &&
             alpha > 1e-7)
        alpha *= 0.5;

      const Tensor<1, spacedim> old_residual = residual;
      while (alpha > 1e-7)
        {
          Point<dim> guess = chart_point + alpha * update;
          residual =
            point - compute_transfinite_interpolation(
                      *cell, guess, coarse_cell_is_flat[cell->index()]);
          const double residual_norm_new = residual.norm_square();
          if (residual_norm_new < residual_norm_square)
            {
              residual_norm_square = residual_norm_new;
              chart_point += alpha * update;
              break;
            }
          else
            alpha *= 0.5;
        }
      if (alpha < 1e-7)
        return outside;

      // update the inverse Jacobian with "Broyden's good method" and
      // Sherman-Morrison formula for the update of the inverse, see
      // https://en.wikipedia.org/wiki/Broyden%27s_method
      const Tensor<1, dim> delta_x = alpha * update;

      // switch sign in residual as compared to the wikipedia article because
      // we use a negative definition of the residual with respect to the
      // Jacobian
      const Tensor<1, spacedim> delta_f = old_residual - residual;

      Tensor<1, dim> Jinv_deltaf;
      for (unsigned int d = 0; d < spacedim; ++d)
        for (unsigned int e = 0; e < dim; ++e)
          Jinv_deltaf[e] += inv_grad[d][e] * delta_f[d];
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
                    Utilities::to_string(level_coarse) + " but there are now" +
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
      for (unsigned int vertex_n = 0;
           vertex_n < GeometryInfo<dim>::vertices_per_cell;
           ++vertex_n)
        {
          vertices[vertex_n] = cell->vertex(vertex_n);
        }

      // cheap check: if any of the points is not inside a circle around the
      // center of the loop, we can skip the expensive part below (this assumes
      // that the manifold does not deform the grid too much)
      Point<spacedim> center;
      for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
        center += vertices[v];
      center *= 1. / GeometryInfo<dim>::vertices_per_cell;
      double radius_square = 0.;
      for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
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
            cell->real_to_unit_cell_affine_approximation(points[i]);
          current_distance += GeometryInfo<dim>::distance_to_unit_cell(point);
        }
      distances_and_cells.emplace_back(current_distance, cell->index());
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
  // (either a cell in 2D or the face of a hexahedron in 3D), arranged like
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
          Assert(false, ExcInternalError());
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
          Assert(false, ExcInternalError());
      }

    return Point<dim>();
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
        // structdim 3 optimization if we are in 3D
        use_structdim_3_guesses = true;
    }
  // we should enable at most one of the optimizations
  Assert((!use_structdim_2_guesses && !use_structdim_3_guesses) ||
           (use_structdim_2_guesses ^ use_structdim_3_guesses),
         ExcInternalError());

  // check whether all points are inside the unit cell of the current chart
  for (unsigned int c = 0; c < nearby_cells.size(); ++c)
    {
      typename Triangulation<dim, spacedim>::cell_iterator cell(
        triangulation, level_coarse, nearby_cells[c]);
      bool inside_unit_cell = true;
      for (unsigned int i = 0; i < surrounding_points.size(); ++i)
        {
          Point<dim> guess;
          // an optimization: keep track of whether or not we used the affine
          // approximation so that we don't call pull_back with the same
          // initial guess twice (i.e., if pull_back fails the first time,
          // don't try again with the same function arguments).
          bool used_affine_approximation = false;
          // if we have already computed three points, we can guess the fourth
          // to be the missing corner point of a rectangle
          if (i == 3 && surrounding_points.size() == 8)
            guess = chart_points[1] + (chart_points[2] - chart_points[0]);
          else if (use_structdim_2_guesses && 3 < i)
            guess = guess_chart_point_structdim_2(i);
          else if (use_structdim_3_guesses && 4 < i)
            guess = guess_chart_point_structdim_3(i);
          else
            {
              guess = cell->real_to_unit_cell_affine_approximation(
                surrounding_points[i]);
              used_affine_approximation = true;
            }
          chart_points[i] = pull_back(cell, surrounding_points[i], guess);

          // the initial guess may not have been good enough: if applicable,
          // try again with the affine approximation (which is more accurate
          // than the cheap methods used above)
          if (chart_points[i][0] == internal::invalid_pull_back_coordinate &&
              !used_affine_approximation)
            {
              guess = cell->real_to_unit_cell_affine_approximation(
                surrounding_points[i]);
              chart_points[i] = pull_back(cell, surrounding_points[i], guess);
            }

          // Tolerance 1e-6 chosen that the method also works with
          // SphericalManifold
          if (GeometryInfo<dim>::is_inside_unit_cell(chart_points[i], 1e-6) ==
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
      // iterate being the end of the array or an unvalid tag), we must stop
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
              for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell;
                   ++v)
                message << cell->vertex(v) << "    ";
              message << std::endl;
              message << "Transformation to chart coordinates: " << std::endl;
              for (unsigned int i = 0; i < surrounding_points.size(); ++i)
                {
                  message << surrounding_points[i] << " -> "
                          << pull_back(
                               cell,
                               surrounding_points[i],
                               cell->real_to_unit_cell_affine_approximation(
                                 surrounding_points[i]))
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
  Assert(false, ExcInternalError());
  return typename Triangulation<dim, spacedim>::cell_iterator();
}



template <int dim, int spacedim>
Point<spacedim>
TransfiniteInterpolationManifold<dim, spacedim>::get_new_point(
  const ArrayView<const Point<spacedim>> &surrounding_points,
  const ArrayView<const double> &         weights) const
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
  const Table<2, double> &                weights,
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
  chart_manifold.get_new_points(
    chart_points_view,
    weights,
    make_array_view(new_points_on_chart.begin(), new_points_on_chart.end()));

  for (unsigned int row = 0; row < weights.size(0); ++row)
    new_points[row] = push_forward(cell, new_points_on_chart[row]);
}



// explicit instantiations
#include "manifold_lib.inst"

DEAL_II_NAMESPACE_CLOSE
