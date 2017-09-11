// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2016 by the deal.II authors
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

#include <deal.II/base/tensor.h>
#include <deal.II/base/table.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/fe/mapping.h>
#include <cmath>

DEAL_II_NAMESPACE_OPEN

// ============================================================
// PolarManifold
// ============================================================

template <int dim, int spacedim>
PolarManifold<dim,spacedim>::PolarManifold(const Point<spacedim> center):
  ChartManifold<dim,spacedim,spacedim>(PolarManifold<dim,spacedim>::get_periodicity()),
  center(center)
{}

template <int dim, int spacedim>
Tensor<1,spacedim>
PolarManifold<dim,spacedim>::get_periodicity()
{
  Tensor<1,spacedim> periodicity;
  // In two dimensions, theta is periodic.
  // In three dimensions things are a little more complicated, since the only variable
  // that is truly periodic is phi, while theta should be bounded between
  // 0 and pi. There is currently no way to enforce this, so here we only fix
  // periodicity for the last variable, corresponding to theta in 2d and phi in 3d.
  periodicity[spacedim-1] = 2*numbers::PI;
  return periodicity;
}

template <int dim, int spacedim>
Point<spacedim>
PolarManifold<dim,spacedim>::push_forward(const Point<spacedim> &spherical_point) const
{
  Assert(spherical_point[0] >=0.0,
         ExcMessage("Negative radius for given point."));
  const double rho = spherical_point[0];
  const double theta = spherical_point[1];

  Point<spacedim> p;
  if (rho > 1e-10)
    switch (spacedim)
      {
      case 2:
        p[0] = rho*cos(theta);
        p[1] = rho*sin(theta);
        break;
      case 3:
      {
        const double phi= spherical_point[2];
        p[0] = rho*sin(theta)*cos(phi);
        p[1] = rho*sin(theta)*sin(phi);
        p[2] = rho*cos(theta);
        break;
      }
      default:
        Assert(false, ExcNotImplemented());
      }
  return p+center;
}

template <int dim, int spacedim>
Point<spacedim>
PolarManifold<dim,spacedim>::pull_back(const Point<spacedim> &space_point) const
{
  const Tensor<1,spacedim> R = space_point-center;
  const double rho = R.norm();

  Point<spacedim> p;
  p[0] = rho;

  switch (spacedim)
    {
    case 2:
    {
      p[1] = atan2(R[1],R[0]);
      if (p[1] < 0)
        p[1] += 2*numbers::PI;
      break;
    }

    case 3:
    {
      const double z = R[2];
      p[2] = atan2(R[1],R[0]); // phi
      if (p[2] < 0)
        p[2] += 2*numbers::PI; // phi is periodic
      p[1] = atan2(sqrt(R[0]*R[0]+R[1]*R[1]),z);  // theta
      break;
    }

    default:
      Assert(false, ExcNotImplemented());
    }
  return p;
}

template <int dim, int spacedim>
DerivativeForm<1,spacedim,spacedim>
PolarManifold<dim,spacedim>::push_forward_gradient(const Point<spacedim> &spherical_point) const
{
  Assert(spherical_point[0] >= 0.0,
         ExcMessage("Negative radius for given point."));
  const double rho = spherical_point[0];
  const double theta = spherical_point[1];

  DerivativeForm<1,spacedim,spacedim> DX;
  if (rho > 1e-10)
    switch (spacedim)
      {
      case 2:
      {
        DX[0][0] = cos(theta);
        DX[0][1] = -rho*sin(theta);
        DX[1][0] = sin(theta);
        DX[1][1] = rho*cos(theta);
        break;
      }

      case 3:
      {
        const double phi= spherical_point[2];
        DX[0][0] =      sin(theta)*cos(phi);
        DX[0][1] =  rho*cos(theta)*cos(phi);
        DX[0][2] = -rho*sin(theta)*sin(phi);

        DX[1][0] =     sin(theta)*sin(phi);
        DX[1][1] = rho*cos(theta)*sin(phi);
        DX[1][2] = rho*sin(theta)*cos(phi);

        DX[2][0] =      cos(theta);
        DX[2][1] = -rho*sin(theta);
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
SphericalManifold<dim,spacedim>::SphericalManifold(const Point<spacedim> center):
  center(center)
{}

template <int dim, int spacedim>
Point<spacedim>
SphericalManifold<dim,spacedim>::
get_intermediate_point (const Point<spacedim> &p1,
                        const Point<spacedim> &p2,
                        const double w) const
{
  Assert(w >=0.0 && w <= 1.0,
         ExcMessage("w should be in the range [0.0,1.0]."));

  const double tol = 1e-10;

  if ( p1.distance(p2) < tol || w < tol)
    return p1;
  else if (w > 1.0 - tol)
    return p2;

  // If the points are one dimensional then there is no need for anything but
  // a linear combination.
  if (spacedim == 1)
    return Point<spacedim>(w*p2 + (1-w)*p1);

  const Tensor<1,spacedim> v1 = p1 - center;
  const Tensor<1,spacedim> v2 = p2 - center;
  const double r1 = v1.norm();
  const double r2 = v2.norm();

  Assert(r1 > tol && r2 > tol,
         ExcMessage("p1 and p2 cannot coincide with the center."));

  const Tensor<1,spacedim> e1 = v1/r1;
  const Tensor<1,spacedim> e2 = v2/r2;

  if ((e1 - e2).norm_square() < tol*tol)
    return Point<spacedim>(center + w*v2 + (1-w)*v1);

  // Find the angle gamma described by v1 and v2:
  const double gamma = std::acos(e1*e2);

  // Find the angle sigma that corresponds to arclength equal to w
  const double sigma = w * gamma;

  // Normal to v1 in the plane described by v1,v2,and the origin.
  // Since p1 and p2 do not coincide n is not zero and well defined.
  Tensor<1,spacedim> n = v2 - (v2*e1)*e1;
  Assert( n.norm() > 0,
          ExcInternalError("n should be different from the null vector."
                           "Probably this means v1==v2 or v2==0."));

  n /= n.norm();

  // Find the point Q along O,v1 such that
  // P1,V,P2 has measure sigma.
  const Tensor<1,spacedim> P = std::cos(sigma) * e1 + std::sin(sigma) * n;

  // Project this point on the manifold.
  return Point<spacedim>(center + (w*r2+(1.0-w)*r1)*P);
}

template <int dim, int spacedim>
Tensor<1,spacedim>
SphericalManifold<dim,spacedim>::
get_tangent_vector (const Point<spacedim> &p1,
                    const Point<spacedim> &p2) const
{
  Assert(p1 != p2,
         ExcMessage("p1 and p2 should not concide."));

  const double r1 = (p1 - center).norm();
  const double r2 = (p2 - center).norm();

  const double tolerance = 1e-10;

  Assert(r1 > tolerance,
         ExcMessage("p1 cannot coincide with the center."));

  Assert(r2 > tolerance,
         ExcMessage("p2 cannot coincide with the center."));

  const Tensor<1,spacedim> e1 = (p1 - center)/r1;
  const Tensor<1,spacedim> e2 = (p2 - center)/r2;

  Assert(e1*e2 + 1.0 > tolerance,
         ExcMessage("p1 and p2 cannot lie on the same diameter and be opposite "
                    "respect to the center."));

  // Tangent vector to the unit sphere along the geodesic given by e1 and e2.
  Tensor<1,spacedim> tg = (e2-(e2*e1)*e1);

  // There is a special case if e2*e1==1.0, in which case tg==0
  const double tg_norm = tg.norm();
  if (tg_norm < tolerance)
    return p2-p1;
  else
    tg /= tg_norm;

  const double gamma = std::acos(e1*e2);

  return (r1-r2)*e1 + r1*gamma*tg;
}



template <int dim, int spacedim>
Point<spacedim>
SphericalManifold<dim,spacedim>::
get_new_point (const ArrayView<const Point<spacedim>> &vertices,
               const ArrayView<const double>          &weights) const
{
  const unsigned int n_points = vertices.size();

  double rho = 0.0;
  Tensor<1,spacedim> candidate;
  for (unsigned int i = 0; i<n_points; i++)
    {
      const Tensor<1,spacedim> direction (vertices[i]-center);
      rho += direction.norm()*weights[i];
      candidate += direction*weights[i];
    }

  // Unit norm direction.
  candidate /= candidate.norm();

  return center+rho*candidate;
}

// ============================================================
// CylindricalManifold
// ============================================================

namespace internal
{
  namespace CylindricalManifold
  {
    namespace
    {
      // helper function to compute a vector orthogonal to a given one.
      template <int spacedim>
      Point<spacedim>
      compute_normal(const Tensor<1,spacedim> &vector)
      {
        AssertThrow(vector.norm() != 0.,
                    ExcMessage("The direction parameter must not be zero!"));
        Point<3> normal;
        if (std::abs(vector[0]) >= std::abs(vector[1])
            && std::abs(vector[0]) >= std::abs(vector[2]))
          {
            normal[1]=-1.;
            normal[2]=-1.;
            normal[0]=(vector[1]+vector[2])/vector[0];
          }
        else if (std::abs(vector[1]) >= std::abs(vector[0])
                 && std::abs(vector[1]) >= std::abs(vector[2]))
          {
            normal[0]=-1.;
            normal[2]=-1.;
            normal[1]=(vector[0]+vector[2])/vector[1];
          }
        else
          {
            normal[0]=-1.;
            normal[1]=-1.;
            normal[2]=(vector[0]+vector[1])/vector[2];
          }
        return normal;
      }
    }
  }
}



template <int dim, int spacedim>
CylindricalManifold<dim, spacedim>::CylindricalManifold(const unsigned int axis,
                                                        const double tolerance) :
  CylindricalManifold<dim, spacedim>(Point<3>::unit_vector(axis),
                                     Point<3>(),
                                     tolerance)
{}



template <int dim, int spacedim>
CylindricalManifold<dim, spacedim>::CylindricalManifold(const Point<spacedim> &direction_,
                                                        const Point<spacedim> &point_on_axis_,
                                                        const double tolerance) :
  ChartManifold<dim,3,3>(Tensor<1,3>({0,2.*numbers::PI,0})),
              normal_direction(internal::CylindricalManifold::compute_normal(direction_)),
              direction (direction_/direction_.norm()),
              point_on_axis (point_on_axis_),
              tolerance(tolerance)
{}



template <int dim, int spacedim>
Point<spacedim>
CylindricalManifold<dim,spacedim>::
get_new_point (const ArrayView<const Point<spacedim>> &surrounding_points,
               const ArrayView<const double>          &weights) const
{
  // First check if the average in space lies on the axis.
  Point<spacedim> middle;
  double average_length = 0.;
  for (unsigned int i=0; i<surrounding_points.size(); ++i)
    {
      middle += surrounding_points[i]*weights[i];
      average_length += surrounding_points[i].square()*weights[i];
    }
  middle -= point_on_axis;
  const double lambda = middle*direction;

  if ((middle-direction*lambda).square() < tolerance*average_length)
    return Point<spacedim>()+direction*lambda;
  else // If not, using the ChartManifold should yield valid results.
    return ChartManifold<dim, spacedim, 3>::get_new_point(surrounding_points,
                                                          weights);
}



template <int dim, int spacedim>
Point<3>
CylindricalManifold<dim, spacedim>::pull_back(const Point<spacedim> &space_point) const
{
  // First find the projection of the given point to the axis.
  const Tensor<1,3> normalized_point = space_point-point_on_axis;
  const double lambda = normalized_point*direction;
  const Point<3> projection = point_on_axis + direction*lambda;
  const Tensor<1,3> p_diff = space_point - projection;

  // Then compute the angle between the projection direction and
  // another vector orthogonal to the direction vector.
  const double dot = normal_direction*p_diff;
  const double det = direction*cross_product_3d(normal_direction, p_diff);
  const double phi = std::atan2(det, dot);

  // Return distance from the axis, angle and signed distance on the axis.
  return Point<3> (p_diff.norm(), phi, lambda);
}



template <int dim, int spacedim>
Point<spacedim>
CylindricalManifold<dim, spacedim>::push_forward(const Point<3> &chart_point) const
{
  // Rotate the orthogonal direction by the given angle.
  // Formula from Section 5.2 in
  // http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/
  // simplified assuming normal_direction and direction are orthogonal.
  const double sine = std::sin(chart_point(1));
  const double cosine = std::cos(chart_point(1));
  const double x = normal_direction[0]*cosine
                   -sine*(direction[2]*normal_direction[1]
                          -direction[1]*normal_direction[2]);
  const double y = normal_direction[1]*cosine
                   -sine*(direction[0]*normal_direction[2]
                          -direction[2]*normal_direction[0]);
  const double z = normal_direction[2]*cosine
                   -sine*(direction[1]*normal_direction[0]
                          -direction[0]*normal_direction[1]);

  // Rescale according to the given distance from the axis.
  Point<3> intermediate (x,y,z);
  intermediate *= chart_point(0)/std::sqrt(intermediate.square());

  // Finally, put everything together.
  return intermediate+point_on_axis+direction*chart_point(2);
}



// ============================================================
// FunctionManifold
// ============================================================
template <int dim, int spacedim, int chartdim>
FunctionManifold<dim,spacedim,chartdim>::FunctionManifold
(const Function<chartdim> &push_forward_function,
 const Function<spacedim> &pull_back_function,
 const Tensor<1,chartdim> &periodicity,
 const double tolerance):
  ChartManifold<dim,spacedim,chartdim>(periodicity),
  push_forward_function(&push_forward_function),
  pull_back_function(&pull_back_function),
  tolerance(tolerance),
  owns_pointers(false)
{
  AssertDimension(push_forward_function.n_components, spacedim);
  AssertDimension(pull_back_function.n_components, chartdim);
}



template <int dim, int spacedim, int chartdim>
FunctionManifold<dim,spacedim,chartdim>::FunctionManifold
(const std::string push_forward_expression,
 const std::string pull_back_expression,
 const Tensor<1,chartdim> &periodicity,
 const typename FunctionParser<spacedim>::ConstMap const_map,
 const std::string chart_vars,
 const std::string space_vars,
 const double tolerance,
 const double h) :
  ChartManifold<dim,spacedim,chartdim>(periodicity),
  const_map(const_map),
  tolerance(tolerance),
  owns_pointers(true)
{
  FunctionParser<chartdim> *pf = new FunctionParser<chartdim>(spacedim, 0.0, h);
  FunctionParser<spacedim> *pb = new FunctionParser<spacedim>(chartdim, 0.0, h);
  pf->initialize(chart_vars, push_forward_expression, const_map);
  pb->initialize(space_vars, pull_back_expression, const_map);
  push_forward_function = pf;
  pull_back_function = pb;
}



template <int dim, int spacedim, int chartdim>
FunctionManifold<dim,spacedim,chartdim>::~FunctionManifold()
{
  if (owns_pointers == true)
    {
      const Function<chartdim> *pf = push_forward_function;
      push_forward_function = nullptr;
      delete pf;

      const Function<spacedim> *pb = pull_back_function;
      pull_back_function = nullptr;
      delete pb;
    }
}



template <int dim, int spacedim, int chartdim>
Point<spacedim>
FunctionManifold<dim,spacedim,chartdim>::push_forward(const Point<chartdim> &chart_point) const
{
  Vector<double> pf(spacedim);
  Point<spacedim> result;
  push_forward_function->vector_value(chart_point, pf);
  for (unsigned int i=0; i<spacedim; ++i)
    result[i] = pf[i];

#ifdef DEBUG
  Vector<double> pb(chartdim);
  pull_back_function->vector_value(result, pb);
  for (unsigned int i=0; i<chartdim; ++i)
    Assert((chart_point.norm() > tolerance &&
            (std::abs(pb[i]-chart_point[i]) < tolerance*chart_point.norm())) ||
           (std::abs(pb[i]-chart_point[i]) < tolerance),
           ExcMessage("The push forward is not the inverse of the pull back! Bailing out."));
#endif

  return result;
}



template <int dim, int spacedim, int chartdim>
DerivativeForm<1,chartdim, spacedim>
FunctionManifold<dim,spacedim,chartdim>::push_forward_gradient(const Point<chartdim> &chart_point) const
{
  DerivativeForm<1, chartdim, spacedim> DF;
  for (unsigned int i=0; i<spacedim; ++i)
    {
      const auto gradient = push_forward_function->gradient(chart_point, i);
      for (unsigned int j=0; j<chartdim; ++j)
        DF[i][j] = gradient[j];
    }
  return DF;
}



template <int dim, int spacedim, int chartdim>
Point<chartdim>
FunctionManifold<dim,spacedim,chartdim>::pull_back(const Point<spacedim> &space_point) const
{
  Vector<double> pb(chartdim);
  Point<chartdim> result;
  pull_back_function->vector_value(space_point, pb);
  for (unsigned int i=0; i<chartdim; ++i)
    result[i] = pb[i];
  return result;
}



template <int dim>
Point<3>
TorusManifold<dim>::pull_back(const Point<3> &p) const
{
  double x = p(0);
  double z = p(1);
  double y = p(2);
  double phi = atan2(y, x);
  double theta = atan2(z, std::sqrt(x*x+y*y)-R);
  double w = std::sqrt(pow(y-sin(phi)*R, 2.0)+pow(x-cos(phi)*R, 2.0)+z*z)/r;
  return Point<3>(phi, theta, w);
}



template <int dim>
Point<3>
TorusManifold<dim>::push_forward(const Point<3> &chart_point) const
{
  double phi = chart_point(0);
  double theta = chart_point(1);
  double w = chart_point(2);

  return Point<3>(cos(phi)*R + r*w*cos(theta)*cos(phi),
                  r*w*sin(theta),
                  sin(phi)*R + r*w*cos(theta)*sin(phi)
                 );
}



template <int dim>
TorusManifold<dim>::TorusManifold (const double R, const double r)
  : ChartManifold<dim,3,3> (Point<3>(2*numbers::PI, 2*numbers::PI, 0.0)),
    r(r),
    R(R)
{
  Assert (R>r, ExcMessage("Outer radius R must be greater than the inner "
                          "radius r."));
  Assert (r>0.0, ExcMessage("inner radius must be positive."));
}



template <int dim>
DerivativeForm<1,3,3>
TorusManifold<dim>::push_forward_gradient(const Point<3> &chart_point) const
{
  DerivativeForm<1,spacedim,spacedim> DX;

  double phi = chart_point(0);
  double theta = chart_point(1);
  double w = chart_point(2);

  DX[0][0] = -sin(phi)*R - r*w*cos(theta)*sin(phi);
  DX[0][1] = -r*w*sin(theta)*cos(phi);
  DX[0][2] = r*cos(theta)*cos(phi);

  DX[1][0] = 0;
  DX[1][1] = r*w*cos(theta);
  DX[1][2] = r*sin(theta);

  DX[2][0] = cos(phi)*R + r*w*cos(theta)*cos(phi);
  DX[2][1] = -r*w*sin(theta)*sin(phi);
  DX[2][2] = r*cos(theta)*sin(phi);

  return DX;
}



template <int dim, int spacedim>
TransfiniteInterpolationManifold<dim,spacedim>::TransfiniteInterpolationManifold()
  :
  triangulation(nullptr),
  level_coarse (-1)
{
  AssertThrow(dim > 1, ExcNotImplemented());
}



template <int dim, int spacedim>
void
TransfiniteInterpolationManifold<dim,spacedim>
::initialize(const Triangulation<dim,spacedim> &triangulation)
{
  this->triangulation = &triangulation;
  // in case the triangulatoin is cleared, remove the pointers by a signal
  triangulation.signals.clear.connect
  ([&]() -> void {this->triangulation = nullptr; this->level_coarse = -1;});
  level_coarse = triangulation.last()->level();
  coarse_cell_is_flat.resize(triangulation.n_cells(level_coarse), false);
  typename Triangulation<dim,spacedim>::active_cell_iterator
  cell = triangulation.begin(level_coarse),
  endc = triangulation.end(level_coarse);
  for ( ; cell != endc; ++cell)
    {
      bool cell_is_flat = true;
      for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_cell; ++l)
        if (cell->line(l)->manifold_id() != cell->manifold_id() &&
            cell->line(l)->manifold_id() != numbers::invalid_manifold_id)
          cell_is_flat = false;
      if (dim > 2)
        for (unsigned int q=0; q<GeometryInfo<dim>::quads_per_cell; ++q)
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
                                    const Point<1> &chart_point,
                                    const bool      /*cell_is_flat*/)
  {
    return cell.vertex(0) * (1.-chart_point[0]) + cell.vertex(1) * chart_point[0];
  }

  // version for 2D
  template <typename AccessorType>
  Point<AccessorType::space_dimension>
  compute_transfinite_interpolation(const AccessorType &cell,
                                    const Point<2>     &chart_point,
                                    const bool          cell_is_flat)
  {
    const unsigned int spacedim = AccessorType::space_dimension;
    const types::manifold_id my_manifold_id = cell.manifold_id();

    // formula see wikipedia
    // https://en.wikipedia.org/wiki/Transfinite_interpolation
    // S(u,v) = (1-v)c_1(u)+v c_3(u) + (1-u)c_2(v) + u c_4(v) -
    //   [(1-u)(1-v)P_0 + u(1-v) P_1 + (1-u)v P_2 + uv P_3]
    const std::array<Point<spacedim>, GeometryInfo<2>::vertices_per_cell> vertices
    {{cell.vertex(0), cell.vertex(1), cell.vertex(2), cell.vertex(3)}};

    Point<spacedim> new_point;
    if (cell_is_flat)
      for (unsigned int v=0; v<GeometryInfo<2>::vertices_per_cell; ++v)
        new_point += GeometryInfo<2>::d_linear_shape_function(chart_point, v) *
                     vertices[v];
    else
      {
        // We subtract the contribution of the vertices (second line in formula).
        // If a line applies the same manifold as the cell, we also subtract a
        // weighted part of the vertex, so accumulate the final weight of the
        // the vertices while going through the faces (this is a bit artificial
        // in 2D but it becomes clear in 3D where we avoid looking at the faces'
        // orientation and other complications).
        double weights_vertices[GeometryInfo<2>::vertices_per_cell];
        for (unsigned int v=0; v<GeometryInfo<2>::vertices_per_cell; ++v)
          weights_vertices[v] = -GeometryInfo<2>::d_linear_shape_function(chart_point, v);

        // add the contribution from the lines around the cell (first line in
        // formula)
        std::array<double, GeometryInfo<2>::vertices_per_face> weights;
        std::array<Point<spacedim>, GeometryInfo<2>::vertices_per_face> points;
        // note that the views are immutable, but the arrays are not
        const auto weights_view = make_array_view(weights.begin(), weights.end());
        const auto points_view = make_array_view(points.begin(), points.end());

        for (unsigned int line=0; line<GeometryInfo<2>::lines_per_cell; ++line)
          {
            const double my_weight = line%2 ? chart_point[line/2] : 1-chart_point[line/2];
            const double line_point = chart_point[1-line/2];

            // Same manifold or invalid id which will go back to the same class ->
            // adds to the vertices
            if (cell.line(line)->manifold_id() == my_manifold_id ||
                cell.line(line)->manifold_id() == numbers::invalid_manifold_id)
              {
                weights_vertices[GeometryInfo<2>::line_to_cell_vertices(line,0)]
                += my_weight * (1.-line_point);
                weights_vertices[GeometryInfo<2>::line_to_cell_vertices(line,1)]
                += my_weight * line_point;
              }
            else
              {
                points[0] = vertices[GeometryInfo<2>::line_to_cell_vertices(line,0)];
                points[1] = vertices[GeometryInfo<2>::line_to_cell_vertices(line,1)];
                weights[0] = 1. - line_point;
                weights[1] = line_point;
                new_point += my_weight *
                             cell.line(line)->get_manifold().get_new_point(points_view,
                                                                           weights_view);
              }
          }

        // subtract contribution from the vertices (second line in formula)
        for (unsigned int v=0; v<GeometryInfo<2>::vertices_per_cell; ++v)
          new_point += weights_vertices[v] * vertices[v];
      }

    return new_point;
  }

  // version for 3D
  template <typename AccessorType>
  Point<AccessorType::space_dimension>
  compute_transfinite_interpolation(const AccessorType &cell,
                                    const Point<3> &chart_point,
                                    const bool      cell_is_flat)
  {
    const unsigned int dim = AccessorType::dimension;
    const unsigned int spacedim = AccessorType::space_dimension;
    const types::manifold_id my_manifold_id = cell.manifold_id();

    // Same approach as in 2D, but adding the faces, subtracting the edges, and
    // adding the vertices
    Point<spacedim> new_point;

    if (cell_is_flat)
      for (unsigned int v=0; v<GeometryInfo<3>::vertices_per_cell; ++v)
        new_point += GeometryInfo<3>::d_linear_shape_function(chart_point, v) *
                     cell.vertex(v);
    else
      {
        // identify the weights for the vertices and lines to be accumulated
        double weights_vertices[GeometryInfo<3>::vertices_per_cell];
        for (unsigned int v=0; v<GeometryInfo<3>::vertices_per_cell; ++v)
          weights_vertices[v] = GeometryInfo<3>::d_linear_shape_function(chart_point, v);
        double weights_lines[GeometryInfo<3>::lines_per_cell];
        for (unsigned int line=0; line<GeometryInfo<3>::lines_per_cell; ++line)
          weights_lines[line] = 0;

        // start with the contributions of the faces
        std::array<double, GeometryInfo<2>::vertices_per_cell> weights;
        std::array<Point<spacedim>, GeometryInfo<2>::vertices_per_cell> points;
        // note that the views are immutable, but the arrays are not
        const auto weights_view = make_array_view(weights.begin(), weights.end());
        const auto points_view = make_array_view(points.begin(), points.end());

        for (unsigned int face=0; face<GeometryInfo<3>::faces_per_cell; ++face)
          {
            Point<2> quad_point(chart_point[(face/2+1)%3], chart_point[(face/2+2)%3]);
            const double my_weight = face%2 ? chart_point[face/2] : 1-chart_point[face/2];

            if (std::abs(my_weight) < 1e-13)
              continue;

            // same manifold or invalid id which will go back to the same class
            // -> face will interpolate from the surrounding lines and vertices
            if (cell.face(face)->manifold_id() == my_manifold_id ||
                cell.face(face)->manifold_id() == numbers::invalid_manifold_id)
              {
                for (unsigned int line=0; line<GeometryInfo<2>::lines_per_cell; ++line)
                  {
                    const double line_weight = line%2 ? quad_point[line/2] : 1-quad_point[line/2];
                    weights_lines[GeometryInfo<3>::face_to_cell_lines(face, line)] +=
                      my_weight * line_weight;
                  }
                for (unsigned int v=0; v<GeometryInfo<2>::vertices_per_cell; ++v)
                  weights_vertices[GeometryInfo<3>::face_to_cell_vertices(face, v)]
                  -= GeometryInfo<2>::d_linear_shape_function(quad_point, v) * my_weight;
              }
            else
              {
                for (unsigned int v=0; v<GeometryInfo<2>::vertices_per_cell; ++v)
                  {
                    points[v] = cell.vertex(GeometryInfo<3>::face_to_cell_vertices(face,v));
                    weights[v] = GeometryInfo<2>::d_linear_shape_function(quad_point, v);
                  }
                new_point += my_weight *
                             cell.face(face)->get_manifold().get_new_point(points_view,
                                                                           weights_view);
              }
          }

        // next subtract the contributions of the lines
        for (unsigned int line=0; line<GeometryInfo<3>::lines_per_cell; ++line)
          {
            const double line_point = (line < 8 ? chart_point[1-(line%4)/2] : chart_point[2]);
            double my_weight = 0.;
            if (line < 8)
              {
                const unsigned int subline = line%4;
                my_weight = subline % 2 ? chart_point[subline/2] : 1-chart_point[subline/2];
                my_weight *= line/4 ? chart_point[2] : (1-chart_point[2]);
              }
            else
              {
                Point<2> xy(chart_point[0], chart_point[1]);
                my_weight = GeometryInfo<2>::d_linear_shape_function(xy, line-8);
              }
            my_weight -= weights_lines[line];

            if (std::abs(my_weight) < 1e-13)
              continue;

            if (cell.line(line)->manifold_id() == my_manifold_id ||
                cell.line(line)->manifold_id() == numbers::invalid_manifold_id)
              {
                weights_vertices[GeometryInfo<3>::line_to_cell_vertices(line,0)]
                -= my_weight * (1.-line_point);
                weights_vertices[GeometryInfo<3>::line_to_cell_vertices(line,1)]
                -= my_weight * (line_point);
              }
            else
              {
                points[0] = cell.vertex(GeometryInfo<3>::line_to_cell_vertices(line,0));
                points[1] = cell.vertex(GeometryInfo<3>::line_to_cell_vertices(line,1));
                weights[0] = 1. - line_point;
                weights[1] = line_point;
                new_point -= my_weight *
                             cell.line(line)->get_manifold().get_new_point(points_view,
                                                                           weights_view);
              }
          }

        // finally add the contribution of the
        for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
          new_point += weights_vertices[v] * cell.vertex(v);
      }
    return new_point;
  }
}



template <int dim, int spacedim>
Point<spacedim>
TransfiniteInterpolationManifold<dim,spacedim>
::push_forward(const typename Triangulation<dim,spacedim>::cell_iterator &cell,
               const Point<dim> &chart_point) const
{
  AssertDimension(cell->level(), level_coarse);

  // check that the point is in the unit cell which is the current chart
  // Tolerance 1e-6 chosen that the method also works with
  // SphericalManifold
  Assert(GeometryInfo<dim>::is_inside_unit_cell(chart_point, 1e-6),
         ExcMessage("chart_point is not in unit interval"));

  return compute_transfinite_interpolation(*cell, chart_point,
                                           coarse_cell_is_flat[cell->index()]);
}



template <int dim, int spacedim>
DerivativeForm<1,dim,spacedim>
TransfiniteInterpolationManifold<dim,spacedim>
::push_forward_gradient(const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                        const Point<dim> &chart_point) const
{
  // compute the derivative with the help of finite differences
  Point<spacedim> point = compute_transfinite_interpolation(*cell, chart_point,
                                                            coarse_cell_is_flat[cell->index()]);
  DerivativeForm<1,dim,spacedim> grad;
  for (unsigned int d=0; d<dim; ++d)
    {
      Point<dim> modified = chart_point;
      const double step = chart_point[d] > 0.5 ? -1e-8 : 1e-8;

      // avoid checking outside of the unit interval
      modified[d] += step;
      Tensor<1,spacedim> difference =
        compute_transfinite_interpolation(*cell, modified,
                                          coarse_cell_is_flat[cell->index()]) - point;
      for (unsigned int e=0; e<spacedim; ++e)
        grad[e][d] = difference[e] / step;
    }
  return grad;
}



template <int dim, int spacedim>
Point<dim>
TransfiniteInterpolationManifold<dim,spacedim>
::pull_back(const typename Triangulation<dim,spacedim>::cell_iterator &cell,
            const Point<spacedim> &point) const
{
  Point<dim> outside;
  for (unsigned int d=0; d<dim; ++d)
    outside[d] = 20;

  // initial guess from affine approximation and projection to unit cell
  Point<dim> chart_point =
    GeometryInfo<dim>::project_to_unit_cell(cell->real_to_unit_cell_affine_approximation(point));

  // run Newton iteration. As opposed to the various mapping implementations,
  // this class does not throw exception upon failure as those are relatively
  // expensive and failure occurs quite regularly in the implementation of the
  // compute_chart_points method.
  Tensor<1,spacedim> residual = point - compute_transfinite_interpolation(*cell, chart_point,
                                coarse_cell_is_flat[cell->index()]);
  const double tolerance = 1e-21 * cell->diameter() * cell->diameter();
  double residual_norm_square = residual.norm_square();
  for (unsigned int i=0; i<100; ++i)
    {
      if (residual_norm_square < tolerance)
        return chart_point;

      // if the determinant is zero, the mapping is not invertible and we are
      // outside the valid chart region
      DerivativeForm<1,dim,spacedim> grad = push_forward_gradient(cell, chart_point);
      if (grad.determinant() <= 0.0)
        return outside;
      DerivativeForm<1,dim,spacedim> inv_grad = grad.covariant_form();
      Tensor<1,dim> update;
      for (unsigned int d=0; d<spacedim; ++d)
        for (unsigned int e=0; e<dim; ++e)
          update[e] += inv_grad[d][e] * residual[d];

      // Line search, accept step if the residual has decreased
      double alpha = 1.;
      while (alpha > 1e-7)
        {
          Point<dim> guess = chart_point + alpha*update;
          residual = point - compute_transfinite_interpolation(*cell, guess,
                                                               coarse_cell_is_flat[cell->index()]);
          const double residual_norm_new = residual.norm_square();
          if (residual_norm_new < residual_norm_square)
            {
              residual_norm_square = residual_norm_new;
              chart_point += alpha*update;
              break;
            }
          else
            alpha *= 0.5;
        }
      if (alpha < 1e-7)
        return outside;
    }
  return outside;
}



template <int dim, int spacedim>
std::array<unsigned int, 10>
TransfiniteInterpolationManifold<dim,spacedim>
::get_possible_cells_around_points(const ArrayView<const Point<spacedim>> &points) const
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

  // This computes the distance of the surrouding points transformed to the unit
  // cell from the unit cell.
  typename Triangulation<dim,spacedim>::cell_iterator
  cell = triangulation->begin(level_coarse),
  endc = triangulation->end(level_coarse);
  boost::container::small_vector<std::pair<double, unsigned int>, 200> distances_and_cells;
  for ( ; cell != endc; ++cell)
    {
      // only consider cells where the current manifold is attached
      if (&cell->get_manifold() != this)
        continue;

      std::array<Point<spacedim>, GeometryInfo<dim>::vertices_per_cell> vertices;
      for (unsigned int vertex_n = 0; vertex_n < GeometryInfo<dim>::vertices_per_cell; ++vertex_n)
        {
          vertices[vertex_n] = cell->vertex(vertex_n);
        }

      // cheap check: if any of the points is not inside a circle around the
      // center of the loop, we can skip the expensive part below (this assumes
      // that the manifold does not deform the grid too much)
      Point<spacedim> center;
      for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
        center += vertices[v];
      center *= 1./GeometryInfo<dim>::vertices_per_cell;
      double radius_square = 0.;
      for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
        radius_square = std::max(radius_square, (center-vertices[v]).norm_square());
      bool inside_circle = true;
      for (unsigned int i=0; i<points.size(); ++i)
        if ((center-points[i]).norm_square() > radius_square * 1.5)
          {
            inside_circle = false;
            break;
          }
      if (inside_circle == false)
        continue;

      // slightly more expensive search
      double current_distance = 0;
      for (unsigned int i=0; i<points.size(); ++i)
        {
          Point<dim> point = cell->real_to_unit_cell_affine_approximation(points[i]);
          current_distance += GeometryInfo<dim>::distance_to_unit_cell(point);
        }
      distances_and_cells.emplace_back(current_distance, cell->index());
    }
  // no coarse cell could be found -> transformation failed
  AssertThrow(distances_and_cells.size() > 0,
              (typename Mapping<dim,spacedim>::ExcTransformationFailed()));
  std::sort(distances_and_cells.begin(), distances_and_cells.end());
  std::array<unsigned int,10> cells;
  cells.fill(numbers::invalid_unsigned_int);
  for (unsigned int i=0; i<distances_and_cells.size() && i<10; ++i)
    cells[i] = distances_and_cells[i].second;

  return cells;
}



template <int dim, int spacedim>
typename Triangulation<dim,spacedim>::cell_iterator
TransfiniteInterpolationManifold<dim, spacedim>
::compute_chart_points (const ArrayView<const Point<spacedim>> &surrounding_points,
                        ArrayView<Point<dim>>                   chart_points) const
{
  Assert(surrounding_points.size() == chart_points.size(),
         ExcMessage("The chart points array view must be as large as the "
                    "surrounding points array view."));

  std::array<unsigned int,10> nearby_cells =
    get_possible_cells_around_points(surrounding_points);

  // check whether all points are inside the unit cell of the current chart
  for (unsigned int c=0; c<nearby_cells.size(); ++c)
    {
      AssertThrow(nearby_cells[c] != numbers::invalid_unsigned_int,
                  (typename Mapping<dim,spacedim>::ExcTransformationFailed()));
      typename Triangulation<dim,spacedim>::cell_iterator cell(triangulation,
                                                               level_coarse,
                                                               nearby_cells[c]);
      bool inside_unit_cell = true;
      for (unsigned int i=0; i<surrounding_points.size(); ++i)
        {
          chart_points[i] = pull_back(cell, surrounding_points[i]);

          // Tolerance 1e-6 chosen that the method also works with
          // SphericalManifold
          if (GeometryInfo<dim>::is_inside_unit_cell(chart_points[i],
                                                     1e-6) == false)
            {
              inside_unit_cell = false;
              break;
            }
        }
      if (inside_unit_cell == true)
        {
          return cell;
        }
    }

  // a valid inversion should have returned a point above.
  AssertThrow(false,
              (typename Mapping<dim,spacedim>::ExcTransformationFailed()));
  return typename Triangulation<dim,spacedim>::cell_iterator();
}



template <int dim, int spacedim>
Point<spacedim>
TransfiniteInterpolationManifold<dim, spacedim>
::get_new_point (const ArrayView<const Point<spacedim>> &surrounding_points,
                 const ArrayView<const double>          &weights) const
{
  boost::container::small_vector<Point<dim>, 100> chart_points(surrounding_points.size());
  ArrayView<Point<dim>> chart_points_view = make_array_view(chart_points.begin(),
                                                            chart_points.end());
  const auto cell = compute_chart_points(surrounding_points, chart_points_view);

  const Point<dim> p_chart = chart_manifold.get_new_point (chart_points_view, weights);

  return push_forward(cell, p_chart);
}



template <int dim, int spacedim>
void
TransfiniteInterpolationManifold<dim,spacedim>::
get_new_points (const ArrayView<const Point<spacedim>> &surrounding_points,
                const Table<2,double>                  &weights,
                ArrayView<Point<spacedim> >             new_points) const
{
  Assert(weights.size(0) > 0, ExcEmptyObject());
  AssertDimension(surrounding_points.size(), weights.size(1));

  boost::container::small_vector<Point<dim>, 100> chart_points(surrounding_points.size());
  ArrayView<Point<dim>> chart_points_view = make_array_view(chart_points.begin(),
                                                            chart_points.end());
  const auto cell = compute_chart_points(surrounding_points, chart_points_view);

  boost::container::small_vector<Point<dim>, 100> new_points_on_chart(weights.size(0));
  chart_manifold.get_new_points(chart_points_view,
                                weights,
                                make_array_view(new_points_on_chart.begin(),
                                                new_points_on_chart.end()));

  for (unsigned int row=0; row<weights.size(0); ++row)
    new_points[row] = push_forward(cell, new_points_on_chart[row]);
}



// explicit instantiations
#include "manifold_lib.inst"

DEAL_II_NAMESPACE_CLOSE
