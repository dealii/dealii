//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <grid/tria_boundary_lib.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <base/tensor.h>
#include <cmath>

DEAL_II_NAMESPACE_OPEN





template <int dim>
CylinderBoundary<dim>::CylinderBoundary (const double radius,
					 const unsigned int axis)
		:
		radius(radius),
		direction (get_axis_vector (axis)),
		point_on_axis (Point<dim>())
{}


template <int dim>
CylinderBoundary<dim>::CylinderBoundary (const double     radius,
					 const Point<dim> direction,
					 const Point<dim> point_on_axis)
		:
		radius(radius),
		direction (direction / direction.norm()),
		point_on_axis (point_on_axis)
{}


template <int dim>
Point<dim>
CylinderBoundary<dim>::get_axis_vector (const unsigned int axis)
{
  Assert (axis < dim, ExcIndexRange (axis, 0, dim));

  Point<dim> axis_vector;
  axis_vector[axis] = 1;
  return axis_vector;
}



template <int dim>
Point<dim>
CylinderBoundary<dim>::
get_new_point_on_line (const typename Triangulation<dim>::line_iterator &line) const
{
				   // compute a proposed new point
  const Point<dim> middle = StraightBoundary<dim>::get_new_point_on_line (line);

				   // we then have to project this
				   // point out to the given radius
				   // from the axis. to this end, we
				   // have to take into account the
				   // offset point_on_axis and the
				   // direction of the axis
  const Point<dim> vector_from_axis = (middle-point_on_axis) -
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


#if deal_II_dimension >=3

template<>
Point<3>
CylinderBoundary<3>::
get_new_point_on_quad (const Triangulation<3>::quad_iterator &quad) const
{
  const Point<3> middle = StraightBoundary<3>::get_new_point_on_quad (quad);

				   // same algorithm as above
  const unsigned int dim = 3;
  const Point<dim> vector_from_axis = (middle-point_on_axis) -
				      ((middle-point_on_axis) * direction) * direction;
  if (vector_from_axis.norm() <= 1e-10 * middle.norm())
    return middle;
  else
    return (vector_from_axis / vector_from_axis.norm() * radius +
	    ((middle-point_on_axis) * direction) * direction +
	    point_on_axis);
}

#endif


template <int dim>
Point<dim>
CylinderBoundary<dim>::
get_new_point_on_quad (const typename Triangulation<dim>::quad_iterator &) const
{
  Assert (false, ExcImpossibleInDim(dim));
  return Point<dim>();
}



template <int dim>
void
CylinderBoundary<dim>::get_intermediate_points_on_line (
  const typename Triangulation<dim>::line_iterator &line,
  std::vector<Point<dim> > &points) const
{
  if (points.size()==1)
    points[0]=get_new_point_on_line(line);
  else
    get_intermediate_points_between_points(line->vertex(0), line->vertex(1), points);
}


template <int dim>
void
CylinderBoundary<dim>::get_intermediate_points_between_points (
  const Point<dim> &v0,
  const Point<dim> &v1,
  std::vector<Point<dim> > &points) const
{
  const unsigned int n=points.size();
  Assert(n>0, ExcInternalError());

				   // Do a simple linear interpolation
				   // followed by projection, using
				   // the same algorithm as above
  const Point<dim> ds = (v1-v0) / (n+1);

  for (unsigned int i=0; i<n; ++i)
    {
      const Point<dim> middle = v0 + (i+1)*ds;

      const Point<dim> vector_from_axis = (middle-point_on_axis) -
					  ((middle-point_on_axis) * direction) * direction;
      if (vector_from_axis.norm() <= 1e-10 * middle.norm())
	points[i] = middle;
      else
	points[i] = (vector_from_axis / vector_from_axis.norm() * radius +
		     ((middle-point_on_axis) * direction) * direction +
		     point_on_axis);
    }
}


#if deal_II_dimension == 3

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

#endif


template <int dim>
void
CylinderBoundary<dim>::get_intermediate_points_on_quad (
  const typename Triangulation<dim>::quad_iterator &,
  std::vector<Point<dim> > &) const
{
  Assert (false, ExcImpossibleInDim(dim));
}



#if deal_II_dimension == 1

template <>
void
CylinderBoundary<1>::
get_normals_at_vertices (const Triangulation<1>::face_iterator &,
			 Boundary<1>::FaceVertexNormals &) const
{
  Assert (false, ExcImpossibleInDim(1));
}

#endif


template <int dim>
void
CylinderBoundary<dim>::
get_normals_at_vertices (const typename Triangulation<dim>::face_iterator &face,
			 typename Boundary<dim>::FaceVertexNormals &face_vertex_normals) const
{
  for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_face; ++v)
    {
      const Point<dim> vertex = face->vertex(v);

      const Point<dim> vector_from_axis = (vertex-point_on_axis) -
					  ((vertex-point_on_axis) * direction) * direction;

      face_vertex_normals[v] = (vector_from_axis / vector_from_axis.norm());
    }
}



template <int dim>
double
CylinderBoundary<dim>::get_radius () const
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
  const Point<dim> dx = (p1 - p0) / (n + 1);

  Assert (n > 0, ExcInternalError ());

  for (unsigned int i = 0; i < n; ++i)
    {
				       // Compute the current point.
      const Point<dim> x_i = p0 + (i + 1) * dx;
				       // To project this point on the
				       // boundary of the cone we first
				       // compute the orthogonal
				       // projection of this point onto
				       // the axis of the cone.
      const double c = (x_i - x_0) * axis / axis.square ();
      const Point<dim> x_ip = x_0 + c * axis;
				       // Compute the projection of
				       // the middle point on the
				       // boundary of the cone.
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
				   // To project it on the boundary of
				   // the cone we first compute the
				   // orthogonal projection of the
				   // middle point onto the axis of
				   // the cone.
  const double c = (middle - x_0) * axis / axis.square ();
  const Point<dim> middle_p = x_0 + c * axis;
				   // Compute the projection of the
				   // middle point on the boundary
				   // of the cone.
  return middle_p + get_radius (middle_p) * (middle - middle_p) / (middle - middle_p).norm ();
}


#if deal_II_dimension == 3

template<int dim>
Point<dim>
ConeBoundary<dim>::
get_new_point_on_quad (const typename Triangulation<dim>::quad_iterator &quad) const
{
  const Point<dim> axis = x_1 - x_0;
				   // Compute the middle point of the
				   // quad.
  const Point<dim> middle = StraightBoundary<3>::get_new_point_on_quad (quad);
				   // Same algorithm as above: To
				   // project it on the boundary of
				   // the cone we first compute the
				   // orthogonal projection of the
				   // middle point onto the axis of
				   // the cone.
  const double c = (middle - x_0) * axis / axis.square ();
  const Point<dim> middle_p = x_0 + c * axis;
				   // Compute the projection of the
				   // middle point on the boundary
				   // of the cone.
  return middle_p + get_radius (middle_p) * (middle - middle_p) / (middle - middle_p).norm ();
}

#else

template<int dim>
Point<dim>
ConeBoundary<dim>::
get_new_point_on_quad (const typename Triangulation<dim>::quad_iterator &) const
{
  Assert (false, ExcImpossibleInDim (dim));

  return Point<dim>();
}

#endif

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


#if deal_II_dimension == 3

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

#endif

template <int dim>
void
ConeBoundary<dim>::
get_intermediate_points_on_quad (const typename Triangulation<dim>::quad_iterator &,
				 std::vector<Point<dim> > &) const
{
  Assert (false, ExcImpossibleInDim (dim));
}


#if deal_II_dimension == 1

template<>
void
ConeBoundary<1>::
get_normals_at_vertices (const Triangulation<1>::face_iterator &,
			 Boundary<1>::FaceVertexNormals &) const
{
  Assert (false, ExcImpossibleInDim (1));
}

#endif

template<int dim>
void
ConeBoundary<dim>::
get_normals_at_vertices (const typename Triangulation<dim>::face_iterator &face,
			 typename Boundary<dim>::FaceVertexNormals &face_vertex_normals) const
{
  const Point<dim> axis = x_1 - x_0;

  for (unsigned int vertex = 0; vertex < GeometryInfo<dim>::vertices_per_cell; ++vertex)
    {
				       // Compute the orthogonal
				       // projection of the vertex onto
				       // the axis of the cone.
      const double c = (face->vertex (vertex) - x_0) * axis / axis.square ();
      const Point<dim> vertex_p = x_0 + c * axis;
				       // Then compute the vector
				       // pointing from the point
				       // <tt>vertex_p</tt> on the axis
				       // to the vertex.
      const Point<dim> axis_to_vertex = face->vertex (vertex) - vertex_p;

      face_vertex_normals[vertex] = axis_to_vertex / axis_to_vertex.norm ();
    }
}


//======================================================================//

template <int dim, int spacedim>
HyperBallBoundary<dim,spacedim>::HyperBallBoundary (const Point<spacedim> p,
					   const double     radius) :
		center(p), radius(radius), compute_radius_automatically(false)
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



#if deal_II_dimension == 1

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

#endif



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


#if deal_II_dimension == 1

template <>
void
HyperBallBoundary<1>::get_intermediate_points_on_line (
  const Triangulation<1>::line_iterator &,
  std::vector<Point<1> > &) const
{
  Assert (false, ExcImpossibleInDim(1));
}

#endif


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
  const double d_alpha=alpha/(n+1);
  const Point<spacedim> pm=0.5*(v0+v1);

  const double h=std::sqrt(pm.square());

  double beta=0;
  unsigned int left_index=0, right_index=0;
  if ((n+1)%2==0)
    {
				       // if the number of
				       // parts is even insert
				       // the midpoint
      left_index=(n-1)/2;
      right_index=left_index;
      points[left_index]=pm;
      --left_index;
      ++right_index;
      beta=d_alpha;
    }
  else
    {
      beta=(n/2+1)*d_alpha-alpha/2.;
      right_index=n/2;
      left_index=n/2-1;
    }

				   // n even:  m=n/2,
				   // n odd:   m=(n-1)/2
  const unsigned int m=n/2;
  for (unsigned int i=0; i<m ; ++i, ++right_index, --left_index, beta+=d_alpha)
    {
      Assert(right_index<n, ExcInternalError());
      Assert(left_index+1>=1, ExcInternalError());

      double d=h*std::tan(beta);
      points[right_index]=pm+d/length*(v1-v0);
      points[left_index]=pm-d/length*(v1-v0);
    }


				   // project the points from the
				   // straight line to the
				   // HyperBallBoundary
  for (unsigned int i=0; i<n; ++i)
    {
      points[i] *= r / std::sqrt(points[i].square());
      points[i] += center;
    }
}

#if deal_II_dimension == 3

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

#endif


template <int dim, int spacedim>
void
HyperBallBoundary<dim,spacedim>::get_intermediate_points_on_quad (
  const typename Triangulation<dim,spacedim>::quad_iterator &,
  std::vector<Point<spacedim> > &) const
{
  Assert(false, ExcImpossibleInDim(dim));
}



#if deal_II_dimension == 1

template <>
void
HyperBallBoundary<1>::
get_normals_at_vertices (const Triangulation<1>::face_iterator &,
			 Boundary<1>::FaceVertexNormals &) const
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

#endif


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
				   // check whether center of object is
				   // at x==0, since then it belongs
				   // to the plane part of the
				   // boundary
  const Point<dim> line_center = line->center();
  if (line_center(0) == this->center(0))
    return line_center;
  else
    return HyperBallBoundary<dim>::get_new_point_on_line (line);
}



#if deal_II_dimension == 1

template <>
Point<1>
HalfHyperBallBoundary<1>::
get_new_point_on_quad (const Triangulation<1>::quad_iterator &) const
{
  Assert (false, ExcInternalError());
  return Point<1>();
}

#endif



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
				   // check whether center of object is
				   // at x==0, since then it belongs
				   // to the plane part of the
				   // boundary
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
				       // check whether center of
				       // object is at x==0, since
				       // then it belongs to the plane
				       // part of the boundary
      const Point<dim> quad_center = quad->center();
      if (quad_center(0) == this->center(0))
	StraightBoundary<dim>::get_intermediate_points_on_quad (quad, points);
      else
	HyperBallBoundary<dim>::get_intermediate_points_on_quad (quad, points);
    }
}


#if deal_II_dimension == 1

template <>
void
HalfHyperBallBoundary<1>::
get_intermediate_points_on_quad (const Triangulation<1>::quad_iterator &,
				 std::vector<Point<1> > &) const
{
  Assert (false, ExcInternalError());
}

#endif


#if deal_II_dimension == 1

template <>
void
HalfHyperBallBoundary<1>::
get_normals_at_vertices (const Triangulation<1>::face_iterator &,
			 Boundary<1>::FaceVertexNormals &) const
{
  Assert (false, ExcImpossibleInDim(1));
}

#endif


template <int dim>
void
HalfHyperBallBoundary<dim>::
get_normals_at_vertices (const typename Triangulation<dim>::face_iterator &face,
			 typename Boundary<dim>::FaceVertexNormals &face_vertex_normals) const
{
				   // check whether center of object is
				   // at x==0, since then it belongs
				   // to the plane part of the
				   // boundary
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
				       // in 2d, first check whether the two
				       // end points of the line are on the
				       // axis of symmetry. if so, then return
				       // the mid point
      case 2:
      {
	if ((line->vertex(0)(0) == this->center(0))
	    &&
	    (line->vertex(1)(0) == this->center(0)))
	  return (line->vertex(0) + line->vertex(1))/2;
	else
					   // otherwise we are on the outer or
					   // inner part of the shell. proceed
					   // as in the base class
	  return HyperShellBoundary<dim>::get_new_point_on_line (line);
      }

					     // in 3d, a line is a straight
					     // line if it is on the symmetry
					     // plane and if not both of its
					     // end points are on either the
					     // inner or outer sphere
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
					   // otherwise we are on the outer or
					   // inner part of the shell. proceed
					   // as in the base class
	  return HyperShellBoundary<dim>::get_new_point_on_line (line);
      }

      default:
	    Assert (false, ExcNotImplemented());
    }

  return Point<dim>();
}



#if deal_II_dimension == 1

template <>
Point<1>
HalfHyperShellBoundary<1>::
get_new_point_on_quad (const Triangulation<1>::quad_iterator &) const
{
  Assert (false, ExcInternalError());
  return Point<1>();
}

#endif



template <int dim>
Point<dim>
HalfHyperShellBoundary<dim>::
get_new_point_on_quad (const typename Triangulation<dim>::quad_iterator &quad) const
{
				   // if this quad is on the symmetry plane,
				   // take the center point and project it
				   // outward to the same radius as the
				   // centers of the two radial lines
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

				   // otherwise we are on the outer or
				   // inner part of the shell. proceed
				   // as in the base class
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
				       // in 2d, first check whether the two
				       // end points of the line are on the
				       // axis of symmetry. if so, then return
				       // the mid point
      case 2:
      {
	if ((line->vertex(0)(0) == this->center(0))
	    &&
	    (line->vertex(1)(0) == this->center(0)))
	  StraightBoundary<dim>::get_intermediate_points_on_line (line, points);
	else
					   // otherwise we are on the outer or
					   // inner part of the shell. proceed
					   // as in the base class
	  HyperShellBoundary<dim>::get_intermediate_points_on_line (line, points);
      }

					     // in 3d, a line is a straight
					     // line if it is on the symmetry
					     // plane and if not both of its
					     // end points are on either the
					     // inner or outer sphere
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
					   // otherwise we are on the outer or
					   // inner part of the shell. proceed
					   // as in the base class
	  HyperShellBoundary<dim>::get_intermediate_points_on_line (line, points);
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

				   // check whether center of object is
				   // at x==0, since then it belongs
				   // to the plane part of the
				   // boundary
  const Point<dim> quad_center = quad->center();
  if (quad_center(0) == this->center(0))
    StraightBoundary<dim>::get_intermediate_points_on_quad (quad, points);
  else
    HyperShellBoundary<dim>::get_intermediate_points_on_quad (quad, points);
}



#if deal_II_dimension == 1

template <>
void
HalfHyperShellBoundary<1>::
get_intermediate_points_on_quad (const Triangulation<1>::quad_iterator &,
				 std::vector<Point<1> > &) const
{
  Assert (false, ExcInternalError());
}

#endif



#if deal_II_dimension == 1

template <>
void
HalfHyperShellBoundary<1>::
get_normals_at_vertices (const Triangulation<1>::face_iterator &,
			 Boundary<1>::FaceVertexNormals &) const
{
  Assert (false, ExcImpossibleInDim(1));
}

#endif



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



// explicit instantiations
template class CylinderBoundary<deal_II_dimension>;
template class ConeBoundary<deal_II_dimension>;
template class HyperBallBoundary<deal_II_dimension>;
template class HalfHyperBallBoundary<deal_II_dimension>;
template class HyperShellBoundary<deal_II_dimension>;
template class HalfHyperShellBoundary<deal_II_dimension>;

#if deal_II_dimension != 3
template class HyperBallBoundary<deal_II_dimension,deal_II_dimension+1>;
#endif

DEAL_II_NAMESPACE_CLOSE

