// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2013 by the deal.II authors
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

#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/base/tensor.h>
#include <cmath>



DEAL_II_NAMESPACE_OPEN



template <int dim, int spacedim>
CylinderBoundary<dim,spacedim>::CylinderBoundary (const double radius,
                                                  const unsigned int axis)
  :
  radius(radius),
  direction (get_axis_vector (axis)),
  point_on_axis (Point<spacedim>())
{}


template <int dim, int spacedim>
CylinderBoundary<dim,spacedim>::CylinderBoundary (const double           radius,
                                                  const Point<spacedim> &direction,
                                                  const Point<spacedim> &point_on_axis)
  :
  radius(radius),
  direction (direction / direction.norm()),
  point_on_axis (point_on_axis)
{}


template <int dim, int spacedim>
Point<spacedim>
CylinderBoundary<dim,spacedim>::get_axis_vector (const unsigned int axis)
{
  Assert (axis < dim, ExcIndexRange (axis, 0, dim));

  Point<spacedim> axis_vector;
  axis_vector[axis] = 1;
  return axis_vector;
}



template <int dim, int spacedim>
Point<spacedim>
CylinderBoundary<dim,spacedim>::
get_new_point_on_line (const typename Triangulation<dim,spacedim>::line_iterator &line) const
{
  // compute a proposed new point
  const Point<spacedim> middle = StraightBoundary<dim,spacedim>::get_new_point_on_line (line);

  // we then have to project this
  // point out to the given radius
  // from the axis. to this end, we
  // have to take into account the
  // offset point_on_axis and the
  // direction of the axis
  const Point<spacedim> vector_from_axis = (middle-point_on_axis) -
                                           ((middle-point_on_axis) * direction) * direction;
  // scale it to the desired length
  // and put everything back
  // together, unless we have a point
  // on the axis
  if (vector_from_axis.norm() <= 1e-10 * middle.norm())
    return middle;
  else
    return (vector_from_axis / vector_from_axis.norm() * radius +
            ((middle-point_on_axis) * direction) * direction +
            point_on_axis);
}



template<>
Point<3>
CylinderBoundary<3>::
get_new_point_on_quad (const Triangulation<3>::quad_iterator &quad) const
{
  const Point<3> middle = StraightBoundary<3,3>::get_new_point_on_quad (quad);

  // same algorithm as above
  const unsigned int spacedim = 3;

  const Point<spacedim> vector_from_axis = (middle-point_on_axis) -
                                           ((middle-point_on_axis) * direction) * direction;
  if (vector_from_axis.norm() <= 1e-10 * middle.norm())
    return middle;
  else
    return (vector_from_axis / vector_from_axis.norm() * radius +
            ((middle-point_on_axis) * direction) * direction +
            point_on_axis);
}

template<>
Point<3>
CylinderBoundary<2,3>::
get_new_point_on_quad (const Triangulation<2,3>::quad_iterator &quad) const
{
  const Point<3> middle = StraightBoundary<2,3>::get_new_point_on_quad (quad);

  // same algorithm as above
  const unsigned int spacedim = 3;
  const Point<spacedim> vector_from_axis = (middle-point_on_axis) -
                                           ((middle-point_on_axis) * direction) * direction;
  if (vector_from_axis.norm() <= 1e-10 * middle.norm())
    return middle;
  else
    return (vector_from_axis / vector_from_axis.norm() * radius +
            ((middle-point_on_axis) * direction) * direction +
            point_on_axis);
}


template <int dim, int spacedim>
Point<spacedim>
CylinderBoundary<dim,spacedim>::
get_new_point_on_quad (const typename Triangulation<dim,spacedim>::quad_iterator &) const
{
  Assert (false, ExcImpossibleInDim(dim));
  return Point<spacedim>();
}



template <int dim, int spacedim>
void
CylinderBoundary<dim,spacedim>::get_intermediate_points_on_line (
  const typename Triangulation<dim,spacedim>::line_iterator &line,
  std::vector<Point<spacedim> > &points) const
{
  if (points.size()==1)
    points[0]=get_new_point_on_line(line);
  else
    get_intermediate_points_between_points(line->vertex(0), line->vertex(1), points);
}


template <int dim, int spacedim>
void
CylinderBoundary<dim,spacedim>::get_intermediate_points_between_points (
  const Point<spacedim> &v0,
  const Point<spacedim> &v1,
  std::vector<Point<spacedim> > &points) const
{
  const unsigned int n=points.size();
  Assert(n>0, ExcInternalError());

  // Do a simple linear interpolation followed by projection, using the same
  // algorithm as above
  const std::vector<Point<1> > &line_points = this->get_line_support_points(n);

  for (unsigned int i=0; i<n; ++i)
    {
      const double x = line_points[i+1][0];
      const Point<spacedim> middle = (1-x)*v0 + x*v1;

      const Point<spacedim> vector_from_axis = (middle-point_on_axis) -
                                               ((middle-point_on_axis) * direction) * direction;
      if (vector_from_axis.norm() <= 1e-10 * middle.norm())
        points[i] = middle;
      else
        points[i] = (vector_from_axis / vector_from_axis.norm() * radius +
                     ((middle-point_on_axis) * direction) * direction +
                     point_on_axis);
    }
}



template <>
void
CylinderBoundary<3>::get_intermediate_points_on_quad (
  const Triangulation<3>::quad_iterator &quad,
  std::vector<Point<3> > &points) const
{
  if (points.size()==1)
    points[0]=get_new_point_on_quad(quad);
  else
    {
      unsigned int m=static_cast<unsigned int> (std::sqrt(static_cast<double>(points.size())));
      Assert(points.size()==m*m, ExcInternalError());

      std::vector<Point<3> > lp0(m);
      std::vector<Point<3> > lp1(m);

      get_intermediate_points_on_line(quad->line(0), lp0);
      get_intermediate_points_on_line(quad->line(1), lp1);

      std::vector<Point<3> > lps(m);
      for (unsigned int i=0; i<m; ++i)
        {
          get_intermediate_points_between_points(lp0[i], lp1[i], lps);

          for (unsigned int j=0; j<m; ++j)
            points[i*m+j]=lps[j];
        }
    }
}



template <int dim, int spacedim>
void
CylinderBoundary<dim,spacedim>::get_intermediate_points_on_quad (
  const typename Triangulation<dim,spacedim>::quad_iterator &,
  std::vector<Point<spacedim> > &) const
{
  Assert (false, ExcImpossibleInDim(dim));
}




template <>
void
CylinderBoundary<1>::
get_normals_at_vertices (const Triangulation<1>::face_iterator &,
                         Boundary<1,1>::FaceVertexNormals &) const
{
  Assert (false, ExcImpossibleInDim(1));
}




template <int dim, int spacedim>
void
CylinderBoundary<dim,spacedim>::
get_normals_at_vertices (const typename Triangulation<dim,spacedim>::face_iterator &face,
                         typename Boundary<dim,spacedim>::FaceVertexNormals &face_vertex_normals) const
{
  for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_face; ++v)
    {
      const Point<spacedim> vertex = face->vertex(v);

      const Point<spacedim> vector_from_axis = (vertex-point_on_axis) -
                                               ((vertex-point_on_axis) * direction) * direction;

      face_vertex_normals[v] = (vector_from_axis / vector_from_axis.norm());
    }
}



template <int dim, int spacedim>
double
CylinderBoundary<dim,spacedim>::get_radius () const
{
  return radius;
}


//======================================================================//

template<int dim>
ConeBoundary<dim>::ConeBoundary (const double radius_0,
                                 const double radius_1,
                                 const Point<dim> x_0,
                                 const Point<dim> x_1)
  :
  radius_0 (radius_0),
  radius_1 (radius_1),
  x_0 (x_0),
  x_1 (x_1)
{}



template<int dim>
double ConeBoundary<dim>::get_radius (Point<dim> x) const
{
  for (unsigned int i = 0; i < dim; ++i)
    if ((x_1 (i) - x_0 (i)) != 0)
      return (radius_1 - radius_0) * (x (i) - x_0 (i)) / (x_1 (i) - x_0 (i)) + radius_0;

  return 0;
}



template<int dim>
void
ConeBoundary<dim>::
get_intermediate_points_between_points (const Point<dim> &p0,
                                        const Point<dim> &p1,
                                        std::vector<Point<dim> > &points) const
{
  const unsigned int n = points.size ();
  const Point<dim> axis = x_1 - x_0;

  Assert (n > 0, ExcInternalError ());

  const std::vector<Point<1> > &line_points = this->get_line_support_points(n);

  for (unsigned int i=0; i<n; ++i)
    {
      const double x = line_points[i+1][0];

      // Compute the current point.
      const Point<dim> x_i = (1-x)*p0 + x*p1;
      // To project this point on the boundary of the cone we first compute
      // the orthogonal projection of this point onto the axis of the cone.
      const double c = (x_i - x_0) * axis / axis.square ();
      const Point<dim> x_ip = x_0 + c * axis;
      // Compute the projection of the middle point on the boundary of the
      // cone.
      points[i] = x_ip + get_radius (x_ip) *  (x_i - x_ip) / (x_i - x_ip).norm ();
    }
}

template<int dim>
Point<dim>
ConeBoundary<dim>::
get_new_point_on_line (const typename Triangulation<dim>::line_iterator &line) const
{
  const Point<dim> axis = x_1 - x_0;
  // Compute the middle point of the line.
  const Point<dim> middle = StraightBoundary<dim>::get_new_point_on_line (line);
  // To project it on the boundary of the cone we first compute the orthogonal
  // projection of the middle point onto the axis of the cone.
  const double c = (middle - x_0) * axis / axis.square ();
  const Point<dim> middle_p = x_0 + c * axis;
  // Compute the projection of the middle point on the boundary of the cone.
  return middle_p + get_radius (middle_p) * (middle - middle_p) / (middle - middle_p).norm ();
}



template <>
Point<3>
ConeBoundary<3>::
get_new_point_on_quad (const Triangulation<3>::quad_iterator &quad) const
{
  const int dim = 3;

  const Point<dim> axis = x_1 - x_0;
  // Compute the middle point of the quad.
  const Point<dim> middle = StraightBoundary<3,3>::get_new_point_on_quad (quad);
  // Same algorithm as above: To project it on the boundary of the cone we
  // first compute the orthogonal projection of the middle point onto the axis
  // of the cone.
  const double c = (middle - x_0) * axis / axis.square ();
  const Point<dim> middle_p = x_0 + c * axis;
  // Compute the projection of the middle point on the boundary of the cone.
  return middle_p + get_radius (middle_p) * (middle - middle_p) / (middle - middle_p).norm ();
}



template<int dim>
Point<dim>
ConeBoundary<dim>::
get_new_point_on_quad (const typename Triangulation<dim>::quad_iterator &) const
{
  Assert (false, ExcImpossibleInDim (dim));

  return Point<dim>();
}



template<int dim>
void
ConeBoundary<dim>::
get_intermediate_points_on_line (const typename Triangulation<dim>::line_iterator &line,
                                 std::vector<Point<dim> > &points) const
{
  if (points.size () == 1)
    points[0] = get_new_point_on_line (line);
  else
    get_intermediate_points_between_points (line->vertex (0), line->vertex (1), points);
}




template<>
void
ConeBoundary<3>::
get_intermediate_points_on_quad (const Triangulation<3>::quad_iterator &quad,
                                 std::vector<Point<3> > &points) const
{
  if (points.size () == 1)
    points[0] = get_new_point_on_quad (quad);
  else
    {
      unsigned int n = static_cast<unsigned int> (std::sqrt (static_cast<double> (points.size ())));

      Assert (points.size () == n * n, ExcInternalError ());

      std::vector<Point<3> > points_line_0 (n);
      std::vector<Point<3> > points_line_1 (n);

      get_intermediate_points_on_line (quad->line (0), points_line_0);
      get_intermediate_points_on_line (quad->line (1), points_line_1);

      std::vector<Point<3> > points_line_segment (n);

      for (unsigned int i = 0; i < n; ++i)
        {
          get_intermediate_points_between_points (points_line_0[i],
                                                  points_line_1[i],
                                                  points_line_segment);

          for (unsigned int j = 0; j < n; ++j)
            points[i * n + j] = points_line_segment[j];
        }
    }
}



template <int dim>
void
ConeBoundary<dim>::
get_intermediate_points_on_quad (const typename Triangulation<dim>::quad_iterator &,
                                 std::vector<Point<dim> > &) const
{
  Assert (false, ExcImpossibleInDim (dim));
}




template<>
void
ConeBoundary<1>::
get_normals_at_vertices (const Triangulation<1>::face_iterator &,
                         Boundary<1,1>::FaceVertexNormals &) const
{
  Assert (false, ExcImpossibleInDim (1));
}



template<int dim>
void
ConeBoundary<dim>::
get_normals_at_vertices (const typename Triangulation<dim>::face_iterator &face,
                         typename Boundary<dim>::FaceVertexNormals &face_vertex_normals) const
{
  const Point<dim> axis = x_1 - x_0;

  for (unsigned int vertex = 0; vertex < GeometryInfo<dim>::vertices_per_cell; ++vertex)
    {
      // Compute the orthogonal projection of the vertex onto the axis of the
      // cone.
      const double c = (face->vertex (vertex) - x_0) * axis / axis.square ();
      const Point<dim> vertex_p = x_0 + c * axis;
      // Then compute the vector pointing from the point <tt>vertex_p</tt> on
      // the axis to the vertex.
      const Point<dim> axis_to_vertex = face->vertex (vertex) - vertex_p;

      face_vertex_normals[vertex] = axis_to_vertex / axis_to_vertex.norm ();
    }
}


//======================================================================//

template <int dim, int spacedim>
HyperBallBoundary<dim,spacedim>::HyperBallBoundary (const Point<spacedim> p,
                                                    const double     radius)
  :
  center(p),
  radius(radius),
  compute_radius_automatically(false)
{}



template <int dim, int spacedim>
Point<spacedim>
HyperBallBoundary<dim,spacedim>::get_new_point_on_line (const typename Triangulation<dim,spacedim>::line_iterator &line) const
{
  Point<spacedim> middle = StraightBoundary<dim,spacedim>::get_new_point_on_line (line);

  middle -= center;

  double r=0;
  if (compute_radius_automatically)
    {
      const Point<spacedim> vertex_relative = line->vertex(0) - center;
      r = std::sqrt(vertex_relative.square());
    }
  else
    r=radius;
  // project to boundary
  middle *= r / std::sqrt(middle.square());
  middle += center;
  return middle;
}



template <>
Point<1>
HyperBallBoundary<1,1>::
get_new_point_on_quad (const Triangulation<1,1>::quad_iterator &) const
{
  Assert (false, ExcInternalError());
  return Point<1>();
}


template <>
Point<2>
HyperBallBoundary<1,2>::
get_new_point_on_quad (const Triangulation<1,2>::quad_iterator &) const
{
  Assert (false, ExcInternalError());
  return Point<2>();
}



template <int dim, int spacedim>
Point<spacedim>
HyperBallBoundary<dim,spacedim>::
get_new_point_on_quad (const typename Triangulation<dim,spacedim>::quad_iterator &quad) const
{
  Point<spacedim> middle = StraightBoundary<dim,spacedim>::get_new_point_on_quad (quad);

  middle -= center;

  double r=0;
  if (compute_radius_automatically)
    {
      const Point<spacedim> vertex_relative = quad->vertex(0) - center;
      r = std::sqrt(vertex_relative.square());
    }
  else
    r=radius;
  // project to boundary
  middle *= r / std::sqrt(middle.square());

  middle += center;
  return middle;
}



template <>
void
HyperBallBoundary<1>::get_intermediate_points_on_line (
  const Triangulation<1>::line_iterator &,
  std::vector<Point<1> > &) const
{
  Assert (false, ExcImpossibleInDim(1));
}



template <int dim, int spacedim>
void
HyperBallBoundary<dim,spacedim>::get_intermediate_points_on_line (
  const typename Triangulation<dim,spacedim>::line_iterator &line,
  std::vector<Point<spacedim> > &points) const
{
  if (points.size()==1)
    points[0]=get_new_point_on_line(line);
  else
    get_intermediate_points_between_points(line->vertex(0), line->vertex(1), points);
}



template <int dim, int spacedim>
void
HyperBallBoundary<dim,spacedim>::get_intermediate_points_between_points (
  const Point<spacedim> &p0, const Point<spacedim> &p1,
  std::vector<Point<spacedim> > &points) const
{
  const unsigned int n=points.size();
  Assert(n>0, ExcInternalError());

  const Point<spacedim> v0=p0-center,
                        v1=p1-center;
  const double length=std::sqrt((v1-v0).square());

  double eps=1e-12;
  double r=0;
  if (compute_radius_automatically)
    {
      const Point<spacedim> vertex_relative = p0 - center;
      r = std::sqrt(vertex_relative.square());
    }
  else
    r=radius;


  const double r2=r*r;
  Assert(std::fabs(v0.square()-r2)<eps*r2, ExcInternalError());
  Assert(std::fabs(v1.square()-r2)<eps*r2, ExcInternalError());

  const double alpha=std::acos((v0*v1)/std::sqrt(v0.square()*v1.square()));
  const Point<spacedim> pm=0.5*(v0+v1);

  const double h=std::sqrt(pm.square());

  // n even:  m=n/2,
  // n odd:   m=(n-1)/2
  const std::vector<Point<1> > &line_points = this->get_line_support_points(n);
  const unsigned int m=n/2;
  for (unsigned int i=0; i<m ; ++i)
    {
      const double beta = alpha * (line_points[i+1][0]-0.5);
      const double d=h*std::tan(beta);
      points[i]=pm+d/length*(v1-v0);
      points[n-1-i]=pm-d/length*(v1-v0);
    }

  if ((n+1)%2==0)
    // if the number of parts is even insert the midpoint
    points[(n-1)/2]=pm;


  // project the points from the straight line to the HyperBallBoundary
  for (unsigned int i=0; i<n; ++i)
    {
      points[i] *= r / std::sqrt(points[i].square());
      points[i] += center;
    }
}



template <>
void
HyperBallBoundary<3>::get_intermediate_points_on_quad (
  const Triangulation<3>::quad_iterator &quad,
  std::vector<Point<3> > &points) const
{
  if (points.size()==1)
    points[0]=get_new_point_on_quad(quad);
  else
    {
      unsigned int m=static_cast<unsigned int> (std::sqrt(static_cast<double>(points.size())));
      Assert(points.size()==m*m, ExcInternalError());

      std::vector<Point<3> > lp0(m);
      std::vector<Point<3> > lp1(m);

      get_intermediate_points_on_line(quad->line(0), lp0);
      get_intermediate_points_on_line(quad->line(1), lp1);

      std::vector<Point<3> > lps(m);
      for (unsigned int i=0; i<m; ++i)
        {
          get_intermediate_points_between_points(lp0[i], lp1[i], lps);

          for (unsigned int j=0; j<m; ++j)
            points[i*m+j]=lps[j];
        }
    }
}



template <>
void
HyperBallBoundary<2,3>::get_intermediate_points_on_quad (
  const Triangulation<2,3>::quad_iterator &quad,
  std::vector<Point<3> > &points) const
{
  if (points.size()==1)
    points[0]=get_new_point_on_quad(quad);
  else
    {
      unsigned int m=static_cast<unsigned int> (std::sqrt(static_cast<double>(points.size())));
      Assert(points.size()==m*m, ExcInternalError());

      std::vector<Point<3> > lp0(m);
      std::vector<Point<3> > lp1(m);

      get_intermediate_points_on_line(quad->line(0), lp0);
      get_intermediate_points_on_line(quad->line(1), lp1);

      std::vector<Point<3> > lps(m);
      for (unsigned int i=0; i<m; ++i)
        {
          get_intermediate_points_between_points(lp0[i], lp1[i], lps);

          for (unsigned int j=0; j<m; ++j)
            points[i*m+j]=lps[j];
        }
    }
}



template <int dim, int spacedim>
void
HyperBallBoundary<dim,spacedim>::get_intermediate_points_on_quad (
  const typename Triangulation<dim,spacedim>::quad_iterator &,
  std::vector<Point<spacedim> > &) const
{
  Assert(false, ExcImpossibleInDim(dim));
}



template <int dim, int spacedim>
Tensor<1,spacedim>
HyperBallBoundary<dim,spacedim>::
normal_vector (const typename Triangulation<dim,spacedim>::face_iterator &,
               const Point<spacedim> &p) const
{
  const Tensor<1,spacedim> unnormalized_normal = p-center;
  return unnormalized_normal/unnormalized_normal.norm();
}



template <>
void
HyperBallBoundary<1>::
get_normals_at_vertices (const Triangulation<1>::face_iterator &,
                         Boundary<1,1>::FaceVertexNormals &) const
{
  Assert (false, ExcImpossibleInDim(1));
}

template <>
void
HyperBallBoundary<1,2>::
get_normals_at_vertices (const Triangulation<1,2>::face_iterator &,
                         Boundary<1,2>::FaceVertexNormals &) const
{
  Assert (false, ExcImpossibleInDim(1));
}



template <int dim, int spacedim>
void
HyperBallBoundary<dim,spacedim>::
get_normals_at_vertices (const typename Triangulation<dim,spacedim>::face_iterator &face,
                         typename Boundary<dim,spacedim>::FaceVertexNormals &face_vertex_normals) const
{
  for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_face; ++vertex)
    face_vertex_normals[vertex] = face->vertex(vertex)-center;
}



template <int dim, int spacedim>
Point<spacedim>
HyperBallBoundary<dim,spacedim>::get_center () const
{
  return center;
}



template <int dim, int spacedim>
double
HyperBallBoundary<dim,spacedim>::get_radius () const
{
  Assert(!compute_radius_automatically, ExcRadiusNotSet());
  return radius;
}


/* ---------------------------------------------------------------------- */


template <int dim>
HalfHyperBallBoundary<dim>::HalfHyperBallBoundary (const Point<dim> center,
                                                   const double     radius) :
  HyperBallBoundary<dim> (center, radius)
{}



template <int dim>
Point<dim>
HalfHyperBallBoundary<dim>::
get_new_point_on_line (const typename Triangulation<dim>::line_iterator &line) const
{
  // check whether center of object is at x==x_center, since then it belongs
  // to the plane part of the boundary. however, this is not the case if it is
  // at the outer perimeter
  const Point<dim> line_center = line->center();
  const Point<dim> vertices[2] = { line->vertex(0), line->vertex(1) };

  if ((line_center(0) == this->center(0))
      &&
      ((std::fabs(vertices[0].distance(this->center)-this->radius) >
        1e-5*this->radius)
       ||
       (std::fabs(vertices[1].distance(this->center)-this->radius) >
        1e-5*this->radius)))
    return line_center;
  else
    return HyperBallBoundary<dim>::get_new_point_on_line (line);
}



template <>
Point<1>
HalfHyperBallBoundary<1>::
get_new_point_on_quad (const Triangulation<1>::quad_iterator &) const
{
  Assert (false, ExcInternalError());
  return Point<1>();
}



template <int dim>
Point<dim>
HalfHyperBallBoundary<dim>::
get_new_point_on_quad (const typename Triangulation<dim>::quad_iterator &quad) const
{
  const Point<dim> quad_center = quad->center();
  if (quad_center(0) == this->center(0))
    return quad_center;
  else
    return HyperBallBoundary<dim>::get_new_point_on_quad (quad);
}



template <int dim>
void
HalfHyperBallBoundary<dim>::
get_intermediate_points_on_line (const typename Triangulation<dim>::line_iterator &line,
                                 std::vector<Point<dim> > &points) const
{
  // check whether center of object is at x==0, since then it belongs to the
  // plane part of the boundary
  const Point<dim> line_center = line->center();
  if (line_center(0) == this->center(0))
    return StraightBoundary<dim>::get_intermediate_points_on_line (line, points);
  else
    return HyperBallBoundary<dim>::get_intermediate_points_on_line (line, points);
}



template <int dim>
void
HalfHyperBallBoundary<dim>::
get_intermediate_points_on_quad (const typename Triangulation<dim>::quad_iterator &quad,
                                 std::vector<Point<dim> > &points) const
{
  if (points.size()==1)
    points[0]=get_new_point_on_quad(quad);
  else
    {
      // check whether center of object is at x==0, since then it belongs to
      // the plane part of the boundary
      const Point<dim> quad_center = quad->center();
      if (quad_center(0) == this->center(0))
        StraightBoundary<dim>::get_intermediate_points_on_quad (quad, points);
      else
        HyperBallBoundary<dim>::get_intermediate_points_on_quad (quad, points);
    }
}



template <>
void
HalfHyperBallBoundary<1>::
get_intermediate_points_on_quad (const Triangulation<1>::quad_iterator &,
                                 std::vector<Point<1> > &) const
{
  Assert (false, ExcInternalError());
}



template <>
void
HalfHyperBallBoundary<1>::
get_normals_at_vertices (const Triangulation<1>::face_iterator &,
                         Boundary<1,1>::FaceVertexNormals &) const
{
  Assert (false, ExcImpossibleInDim(1));
}



template <int dim>
void
HalfHyperBallBoundary<dim>::
get_normals_at_vertices (const typename Triangulation<dim>::face_iterator &face,
                         typename Boundary<dim>::FaceVertexNormals &face_vertex_normals) const
{
  // check whether center of object is at x==0, since then it belongs to the
  // plane part of the boundary
  const Point<dim> quad_center = face->center();
  if (quad_center(0) == this->center(0))
    StraightBoundary<dim>::get_normals_at_vertices (face, face_vertex_normals);
  else
    HyperBallBoundary<dim>::get_normals_at_vertices (face, face_vertex_normals);
}


/* ---------------------------------------------------------------------- */



template <int dim>
HyperShellBoundary<dim>::HyperShellBoundary (const Point<dim> &center)
  :
  HyperBallBoundary<dim>(center, 0.)
{
  this->compute_radius_automatically=true;
}


/* ---------------------------------------------------------------------- */




template <int dim>
HalfHyperShellBoundary<dim>::HalfHyperShellBoundary (const Point<dim> &center,
                                                     const double inner_radius,
                                                     const double outer_radius)
  :
  HyperShellBoundary<dim> (center),
  inner_radius (inner_radius),
  outer_radius (outer_radius)
{
  if (dim > 2)
    Assert ((inner_radius >= 0) &&
            (outer_radius > 0) &&
            (outer_radius > inner_radius),
            ExcMessage ("Inner and outer radii must be specified explicitly in 3d."));
}



template <int dim>
Point<dim>
HalfHyperShellBoundary<dim>::
get_new_point_on_line (const typename Triangulation<dim>::line_iterator &line) const
{
  switch (dim)
    {
    // in 2d, first check whether the two end points of the line are on the
    // axis of symmetry. if so, then return the mid point
    case 2:
    {
      if ((line->vertex(0)(0) == this->center(0))
          &&
          (line->vertex(1)(0) == this->center(0)))
        return (line->vertex(0) + line->vertex(1))/2;
      else
        // otherwise we are on the outer or inner part of the shell. proceed
        // as in the base class
        return HyperShellBoundary<dim>::get_new_point_on_line (line);
    }

    // in 3d, a line is a straight line if it is on the symmetry plane and if
    // not both of its end points are on either the inner or outer sphere
    case 3:
    {

      if (((line->vertex(0)(0) == this->center(0))
           &&
           (line->vertex(1)(0) == this->center(0)))
          &&
          !(((std::fabs (line->vertex(0).distance (this->center)
                         - inner_radius) < 1e-12 * outer_radius)
             &&
             (std::fabs (line->vertex(1).distance (this->center)
                         - inner_radius) < 1e-12 * outer_radius))
            ||
            ((std::fabs (line->vertex(0).distance (this->center)
                         - outer_radius) < 1e-12 * outer_radius)
             &&
             (std::fabs (line->vertex(1).distance (this->center)
                         - outer_radius) < 1e-12 * outer_radius))))
        return (line->vertex(0) + line->vertex(1))/2;
      else
        // otherwise we are on the outer or inner part of the shell. proceed
        // as in the base class
        return HyperShellBoundary<dim>::get_new_point_on_line (line);
    }

    default:
      Assert (false, ExcNotImplemented());
    }

  return Point<dim>();
}



template <>
Point<1>
HalfHyperShellBoundary<1>::
get_new_point_on_quad (const Triangulation<1>::quad_iterator &) const
{
  Assert (false, ExcInternalError());
  return Point<1>();
}




template <int dim>
Point<dim>
HalfHyperShellBoundary<dim>::
get_new_point_on_quad (const typename Triangulation<dim>::quad_iterator &quad) const
{
  // if this quad is on the symmetry plane, take the center point and project
  // it outward to the same radius as the centers of the two radial lines
  if ((quad->vertex(0)(0) == this->center(0)) &&
      (quad->vertex(1)(0) == this->center(0)) &&
      (quad->vertex(2)(0) == this->center(0)) &&
      (quad->vertex(3)(0) == this->center(0)))
    {
      const Point<dim> quad_center = (quad->vertex(0) + quad->vertex(1) +
                                      quad->vertex(2) + quad->vertex(3)   )/4;
      const Point<dim> quad_center_offset = quad_center - this->center;


      if (std::fabs (quad->line(0)->center().distance(this->center) -
                     quad->line(1)->center().distance(this->center))
          < 1e-12 * outer_radius)
        {
          // lines 0 and 1 are radial
          const double needed_radius
            = quad->line(0)->center().distance(this->center);

          return (this->center +
                  quad_center_offset/quad_center_offset.norm() * needed_radius);
        }
      else if (std::fabs (quad->line(2)->center().distance(this->center) -
                          quad->line(3)->center().distance(this->center))
               < 1e-12 * outer_radius)
        {
          // lines 2 and 3 are radial
          const double needed_radius
            = quad->line(2)->center().distance(this->center);

          return (this->center +
                  quad_center_offset/quad_center_offset.norm() * needed_radius);
        }
      else
        Assert (false, ExcInternalError());
    }

  // otherwise we are on the outer or inner part of the shell. proceed as in
  // the base class
  return HyperShellBoundary<dim>::get_new_point_on_quad (quad);
}



template <int dim>
void
HalfHyperShellBoundary<dim>::
get_intermediate_points_on_line (const typename Triangulation<dim>::line_iterator &line,
                                 std::vector<Point<dim> > &points) const
{
  switch (dim)
    {
    // in 2d, first check whether the two end points of the line are on the
    // axis of symmetry. if so, then return the mid point
    case 2:
    {
      if ((line->vertex(0)(0) == this->center(0))
          &&
          (line->vertex(1)(0) == this->center(0)))
        StraightBoundary<dim>::get_intermediate_points_on_line (line, points);
      else
        // otherwise we are on the outer or inner part of the shell. proceed
        // as in the base class
        HyperShellBoundary<dim>::get_intermediate_points_on_line (line, points);
      break;
    }

    // in 3d, a line is a straight line if it is on the symmetry plane and if
    // not both of its end points are on either the inner or outer sphere
    case 3:
    {
      if (((line->vertex(0)(0) == this->center(0))
           &&
           (line->vertex(1)(0) == this->center(0)))
          &&
          !(((std::fabs (line->vertex(0).distance (this->center)
                         - inner_radius) < 1e-12 * outer_radius)
             &&
             (std::fabs (line->vertex(1).distance (this->center)
                         - inner_radius) < 1e-12 * outer_radius))
            ||
            ((std::fabs (line->vertex(0).distance (this->center)
                         - outer_radius) < 1e-12 * outer_radius)
             &&
             (std::fabs (line->vertex(1).distance (this->center)
                         - outer_radius) < 1e-12 * outer_radius))))
        StraightBoundary<dim>::get_intermediate_points_on_line (line, points);
      else
        // otherwise we are on the outer or inner part of the shell. proceed
        // as in the base class
        HyperShellBoundary<dim>::get_intermediate_points_on_line (line, points);

      break;
    }

    default:
      Assert (false, ExcNotImplemented());
    }
}



template <int dim>
void
HalfHyperShellBoundary<dim>::
get_intermediate_points_on_quad (const typename Triangulation<dim>::quad_iterator &quad,
                                 std::vector<Point<dim> > &points) const
{
  Assert (dim < 3, ExcNotImplemented());

  // check whether center of object is at x==0, since then it belongs to the
  // plane part of the boundary
  const Point<dim> quad_center = quad->center();
  if (quad_center(0) == this->center(0))
    StraightBoundary<dim>::get_intermediate_points_on_quad (quad, points);
  else
    HyperShellBoundary<dim>::get_intermediate_points_on_quad (quad, points);
}



template <>
void
HalfHyperShellBoundary<1>::
get_intermediate_points_on_quad (const Triangulation<1>::quad_iterator &,
                                 std::vector<Point<1> > &) const
{
  Assert (false, ExcInternalError());
}



template <>
void
HalfHyperShellBoundary<1>::
get_normals_at_vertices (const Triangulation<1>::face_iterator &,
                         Boundary<1,1>::FaceVertexNormals &) const
{
  Assert (false, ExcImpossibleInDim(1));
}





template <int dim>
void
HalfHyperShellBoundary<dim>::
get_normals_at_vertices (const typename Triangulation<dim>::face_iterator &face,
                         typename Boundary<dim>::FaceVertexNormals &face_vertex_normals) const
{
  if (face->center()(0) == this->center(0))
    StraightBoundary<dim>::get_normals_at_vertices (face, face_vertex_normals);
  else
    HyperShellBoundary<dim>::get_normals_at_vertices (face, face_vertex_normals);
}




template <int dim, int spacedim>
TorusBoundary<dim,spacedim>::TorusBoundary (const double R__,
                                            const double r__)
  :
  R(R__),
  r(r__)
{
  Assert (false, ExcNotImplemented());
}



template <>
TorusBoundary<2,3>::TorusBoundary (const double R__,
                                   const double r__)
  :
  R(R__),
  r(r__)
{
  Assert (R>r, ExcMessage("Outer radius must be greater than inner radius."));
}



template <int dim, int spacedim>
double
TorusBoundary<dim,spacedim>::get_correct_angle(const double angle,
                                               const double x,
                                               const double y) const
{
  if (y>=0)
    {
      if (x >=0)
        return angle;

      return numbers::PI-angle;
    }

  if (x <=0)
    return numbers::PI+angle;

  return 2.0*numbers::PI-angle;
}



template <>
Point<3>
TorusBoundary<2,3>::get_real_coord (const Point<2> &surfP) const
{
  const double theta=surfP(0);
  const double phi=surfP(1);

  return Point<3> ((R+r*std::cos(phi))*std::cos(theta),
                   r*std::sin(phi),
                   (R+r*std::cos(phi))*std::sin(theta));
}



template <>
Point<2>
TorusBoundary<2,3>::get_surf_coord(const Point<3> &p) const
{
  const double phi=std::asin(std::abs(p(1))/r);
  const double Rr_2=p(0)*p(0)+p(2)*p(2);

  Point<2> surfP;
  surfP(1)=get_correct_angle(phi,Rr_2-R*R,p(1));//phi

  if (std::abs(p(0))<1.E-5)
    {
      if (p(2)>=0)
        surfP(0) =  numbers::PI*0.5;
      else
        surfP(0) = -numbers::PI*0.5;
    }
  else
    {
      const double theta = std::atan(std::abs(p(2)/p(0)));
      surfP(0)=  get_correct_angle(theta,p(0),p(2));
    }

  return surfP;
}



template <>
Point<3>
TorusBoundary<2,3>::get_new_point_on_line (const Triangulation<2,3>::line_iterator &line) const
{
  //Just get the average
  Point<2>  p0=get_surf_coord(line->vertex(0));
  Point<2>  p1=get_surf_coord(line->vertex(1));

  Point<2>  middle(0,0);

  //Take care for periodic conditions, For instance phi0= 0, phi1= 3/2*Pi
  //middle has to be 7/4*Pi not 3/4*Pi. This also works for -Pi/2 + Pi, middle
  //is 5/4*Pi
  for (unsigned int i=0; i<2; i++)
    if (std::abs(p0(i)-p1(i))> numbers::PI)
      middle(i)=2*numbers::PI;

  middle+=  p0 + p1;
  middle*=0.5;

  Point<3> midReal=get_real_coord(middle);
  return midReal;
}



template <>
Point<3>
TorusBoundary<2,3>::get_new_point_on_quad (const Triangulation<2,3>::quad_iterator &quad) const
{
  //Just get the average
  Point<2> p[4];

  for (unsigned int i=0; i<4; i++)
    p[i]=get_surf_coord(quad->vertex(i));

  Point<2>  middle(0,0);

  //Take care for periodic conditions, see get_new_point_on_line() above
  //For instance phi0= 0, phi1= 3/2*Pi  middle has to be 7/4*Pi not 3/4*Pi
  //This also works for -Pi/2 + Pi + Pi- Pi/2, middle is 5/4*Pi
  for (unsigned int i=0; i<2; i++)
    for (unsigned int j=1; j<4; j++)
      {
        if (std::abs(p[0](i)-p[j](i))> numbers::PI)
          {
            middle(i)+=2*numbers::PI;
          }
      }

  for (unsigned int i=0; i<4; i++)
    middle+=p[i];

  middle*= 0.25;

  return get_real_coord(middle);
}



//Normal field without unit length
template <>
Point<3>
TorusBoundary<2,3>:: get_surf_norm_from_sp(const Point<2> &surfP) const
{

  Point<3> n;
  double theta=surfP[0];
  double phi=surfP[1];

  double f=R+r*std::cos(phi);

  n[0]=r*std::cos(phi)*std::cos(theta)*f;
  n[1]=r*std::sin(phi)*f;
  n[2]=r*std::sin(theta)*std::cos(phi)*f;

  return n;
}



//Normal field without unit length
template <>
Point<3>
TorusBoundary<2,3>::get_surf_norm(const Point<3> &p) const
{

  Point<2> surfP=get_surf_coord(p);
  return get_surf_norm_from_sp(surfP);

}



template<>
void
TorusBoundary<2,3>::
get_intermediate_points_on_line (const Triangulation<2, 3>::line_iterator   &line,
                                 std::vector< Point< 3 > > &points) const
{
  //Almost the same implementation as StraightBoundary<2,3>
  unsigned int npoints=points.size();
  if (npoints==0) return;

  Point<2> p[2];

  for (unsigned int i=0; i<2; i++)
    p[i]=get_surf_coord(line->vertex(i));

  unsigned int offset[2];
  offset[0]=0;
  offset[1]=0;

  //Take care for periodic conditions & negative angles, see
  //get_new_point_on_line() above. Because we dont have a symmetric
  //interpolation (just the middle) we need to add 2*Pi to each almost zero
  //and negative angles.
  for (unsigned int i=0; i<2; i++)
    for (unsigned int j=1; j<2; j++)
      {
        if (std::abs(p[0](i)-p[j](i))> numbers::PI)
          {
            offset[i]++;
            break;
          }
      }

  for (unsigned int i=0; i<2; i++)
    for (unsigned int j=0; j<2; j++)
      if (p[j](i)<1.E-12 ) //Take care for periodic conditions & negative angles
        p[j](i)+=2*numbers::PI*offset[i];


  Point<2>  target;
  const std::vector<Point<1> > &line_points = this->get_line_support_points(npoints);
  for (unsigned int i=0; i<npoints; i++)
    {
      const double x = line_points[i+1][0];
      target=  (1-x)*p[0] + x*p[1];
      points[i]=get_real_coord(target);
    }
}



template<>
void
TorusBoundary<2,3>::
get_intermediate_points_on_quad (const Triangulation< 2, 3 >::quad_iterator &quad,
                                 std::vector< Point< 3 > > &points )const
{
  //Almost the same implementation as  StraightBoundary<2,3>
  const unsigned int n=points.size(),
                     m=static_cast<unsigned int>(std::sqrt(static_cast<double>(n)));
  // is n a square number
  Assert(m*m==n, ExcInternalError());

  Point<2>  p[4];

  for (unsigned int i=0; i<4; i++)
    p[i]=get_surf_coord(quad->vertex(i));

  Point<2>  target;
  unsigned int offset[2];
  offset[0]=0;
  offset[1]=0;

  //Take care for periodic conditions & negative angles, see
  //get_new_point_on_line() above.  Because we dont have a symmetric
  //interpolation (just the middle) we need to add 2*Pi to each almost zero
  //and negative angles.
  for (unsigned int i=0; i<2; i++)
    for (unsigned int j=1; j<4; j++)
      {
        if (std::abs(p[0](i)-p[j](i))> numbers::PI)
          {
            offset[i]++;
            break;
          }
      }

  for (unsigned int i=0; i<2; i++)
    for (unsigned int j=0; j<4; j++)
      if (p[j](i)<1.E-12 ) //Take care for periodic conditions & negative angles
        p[j](i)+=2*numbers::PI*offset[i];

  const std::vector<Point<1> > &line_points = this->get_line_support_points(m);
  for (unsigned int i=0; i<m; ++i)
    {
      const double y=line_points[i+1][0];
      for (unsigned int j=0; j<m; ++j)
        {
          const double x=line_points[j+1][0];
          target=((1-x) * p[0] +
                  x     * p[1]) * (1-y) +
                 ((1-x) * p[2] +
                  x     * p[3]) * y;

          points[i*m+j]=get_real_coord(target);
        }
    }
}



template<>
void
TorusBoundary<2,3>::
get_normals_at_vertices (const Triangulation<2,3 >::face_iterator &face,
                         Boundary<2,3>::FaceVertexNormals &face_vertex_normals) const
{
  for (unsigned int i=0; i<GeometryInfo<2>::vertices_per_face; i++)
    face_vertex_normals[i]=get_surf_norm(face->vertex(i));
}



// explicit instantiations
#include "tria_boundary_lib.inst"

DEAL_II_NAMESPACE_CLOSE
