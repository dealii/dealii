//----------------------------  tria_boundary_lib.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  tria_boundary_lib.cc  ---------------------------


#include <grid/tria_boundary_lib.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <base/tensor.h>
#include <cmath>


template <int dim>
CylinderBoundary<dim>::CylinderBoundary (const double radius) :
		radius(radius)
{};



template <int dim>
Point<dim>
CylinderBoundary<dim>::get_new_point_on_line (const typename Triangulation<dim>::line_iterator &line) const
{
  Point<dim> middle = StraightBoundary<dim>::get_new_point_on_line (line);
				   // project to boundary
  if (dim>=3
      && line->vertex(0).square()-line->vertex(0)(0)*line->vertex(0)(0) >= radius*radius-1.e-12
      && line->vertex(1).square()-line->vertex(1)(0)*line->vertex(1)(0) >= radius*radius-1.e-12)
    {
      const double f = radius / std::sqrt(middle.square()-middle(0)*middle(0));
      for (unsigned int i=1;i<dim;++i)
	middle(i) *= f;
    }
  return middle;
};


template <int dim>
Point<dim>
CylinderBoundary<dim>::
get_new_point_on_quad (const typename Triangulation<dim>::quad_iterator &quad) const
{
  Point<dim> middle = StraightBoundary<dim>::get_new_point_on_quad (quad);
  
				   // project to boundary
  if (dim>=3
      && quad->vertex(0).square()-quad->vertex(0)(0)*quad->vertex(0)(0) >= radius*radius-1.e-12
      && quad->vertex(1).square()-quad->vertex(1)(0)*quad->vertex(1)(0) >= radius*radius-1.e-12
      && quad->vertex(2).square()-quad->vertex(2)(0)*quad->vertex(2)(0) >= radius*radius-1.e-12
      && quad->vertex(3).square()-quad->vertex(3)(0)*quad->vertex(3)(0) >= radius*radius-1.e-12)
      
    {
      const double f = radius / std::sqrt(middle.square()-middle(0)*middle(0));
      for (unsigned int i=1;i<dim;++i)
	middle(i) *= f;
    }
  return middle;
};


template <int dim>
void
CylinderBoundary<dim>::get_intermediate_points_on_line (
  const typename Triangulation<dim>::line_iterator &line,
  typename std::vector<Point<dim> > &points) const
{
  if (points.size()==1)
    points[0]=get_new_point_on_line(line);
  else
    get_intermediate_points_between_points(line->vertex(0), line->vertex(1), points);
}

  
template <int dim>
void
CylinderBoundary<dim>::get_intermediate_points_between_points (
  const Point<dim> &v0, const Point<dim> &v1,
  typename std::vector<Point<dim> > &points) const
{
  const unsigned int n=points.size();
  Assert(n>0, ExcInternalError());
				   // Do a simple linear interpolation
				   // followed by projection. If this
				   // is not sufficient, try to
				   // understand the sophisticated
				   // code in HyperBall later.
  Point<dim> ds = v1-v0;
  ds /= n+1;

  bool scale = (dim>=3
		&& v0.square()-v0(0)*v0(0) >= radius*radius-1.e-12
		&& v1.square()-v1(0)*v1(0) >= radius*radius-1.e-12);
  
  for (unsigned int i=0; i<n; ++i)
    {
      if (i==0)
	points[i] = v0+ds;
      else
	points[i] = points[i-1]+ds;

      if (scale)
	{
	  const double f = radius / std::sqrt(points[i].square()-points[i](0)*points[i](0));
	  for (unsigned int d=1;d<dim;++d)
	    points[i](d) *= f;
	}
    }
}


#if deal_II_dimension == 3

template <>
void
CylinderBoundary<3>::get_intermediate_points_on_quad (
  const Triangulation<3>::quad_iterator &quad,
  std::vector<Point<3> > &points) const
{
  unsigned int m=static_cast<unsigned int> (std::sqrt(static_cast<double>(points.size())));
  Assert(points.size()==m*m, ExcInternalError());

  std::vector<Point<3> > lp3(m);
  std::vector<Point<3> > lp1(m);
  
  get_intermediate_points_on_line(quad->line(3), lp3);
  get_intermediate_points_on_line(quad->line(1), lp1);
  
  std::vector<Point<3> > lps(m);
  for (unsigned int i=0; i<m; ++i)
    {
      get_intermediate_points_between_points(lp3[i], lp1[i], lps);
      
      for (unsigned int j=0; j<m; ++j)
	points[i*m+j]=lps[j];
    }
}

#endif


template <int dim>
void
CylinderBoundary<dim>::get_intermediate_points_on_quad (
  const typename Triangulation<dim>::quad_iterator &,
  typename std::vector<Point<dim> > &) const
{
  Assert(false, Boundary<dim>::ExcFunctionNotUseful(dim));
}



#if deal_II_dimension == 1

template <>
void
CylinderBoundary<1>::
get_normals_at_vertices (const Triangulation<1>::face_iterator &,
			 Boundary<1>::FaceVertexNormals &) const
{
  Assert (false, Boundary<1>::ExcFunctionNotUseful(1));
};

#endif


template <int dim>
void
CylinderBoundary<dim>::
get_normals_at_vertices (const typename Triangulation<dim>::face_iterator &face,
			 typename Boundary<dim>::FaceVertexNormals &face_vertex_normals) const
{
  for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_face; ++vertex)
    {
      face_vertex_normals[vertex] = face->vertex(vertex);
      face_vertex_normals[vertex][0] = 0.;
    }
};



template <int dim>
double
CylinderBoundary<dim>::get_radius () const 
{
  return radius;
};


//======================================================================//

template <int dim>
HyperBallBoundary<dim>::HyperBallBoundary (const Point<dim> p,
					   const double     radius) :
		center(p), radius(radius), compute_radius_automatically(false)
{};



template <int dim>
Point<dim>
HyperBallBoundary<dim>::get_new_point_on_line (const typename Triangulation<dim>::line_iterator &line) const
{
  Point<dim> middle = StraightBoundary<dim>::get_new_point_on_line (line);
  
  middle -= center;

  double r=0;
  if (compute_radius_automatically)
    {
      const Point<dim> vertex_relative = line->vertex(0) - center;
      r = std::sqrt(vertex_relative.square());
    }
  else
    r=radius;
				   // project to boundary
  middle *= r / std::sqrt(middle.square());
  
  middle += center;
  return middle;
};



#if deal_II_dimension == 1

template <>
Point<1>
HyperBallBoundary<1>::
get_new_point_on_quad (const Triangulation<1>::quad_iterator &) const
{
  Assert (false, ExcInternalError());
  return Point<1>();
};

#endif



template <int dim>
Point<dim>
HyperBallBoundary<dim>::
get_new_point_on_quad (const typename Triangulation<dim>::quad_iterator &quad) const
{
  Point<dim> middle = StraightBoundary<dim>::get_new_point_on_quad (quad);
  
  middle -= center;

  double r=0;
  if (compute_radius_automatically)
    {
      const Point<dim> vertex_relative = quad->vertex(0) - center;
      r = std::sqrt(vertex_relative.square());
    }
  else
    r=radius;
				   // project to boundary
  middle *= r / std::sqrt(middle.square());
  
  middle += center;
  return middle;
};


#if deal_II_dimension == 1

template <>
void
HyperBallBoundary<1>::get_intermediate_points_on_line (
  const Triangulation<1>::line_iterator &,
  std::vector<Point<1> > &) const
{
  Assert (false, Boundary<1>::ExcFunctionNotUseful(1));
}

#else


template <int dim>
void
HyperBallBoundary<dim>::get_intermediate_points_on_line (
  const typename Triangulation<dim>::line_iterator &line,
  typename std::vector<Point<dim> > &points) const
{
  if (points.size()==1)
    points[0]=get_new_point_on_line(line);
  else
    get_intermediate_points_between_points(line->vertex(0), line->vertex(1), points);
}

  
template <int dim>
void
HyperBallBoundary<dim>::get_intermediate_points_between_points (
  const Point<dim> &p0, const Point<dim> &p1,
  typename std::vector<Point<dim> > &points) const
{
  const unsigned int n=points.size();
  Assert(n>0, ExcInternalError());

  const Point<dim> v0=p0-center,
		   v1=p1-center;
  const double length=std::sqrt((v1-v0).square());
  
  double eps=1e-12;
  double r=0;
  if (compute_radius_automatically)
    {
      const Point<dim> vertex_relative = p0 - center;
      r = std::sqrt(vertex_relative.square());
    }
  else
    r=radius;

  
  const double r2=r*r;
  Assert(std::fabs(v0.square()-r2)<eps*r2, ExcInternalError());
  Assert(std::fabs(v1.square()-r2)<eps*r2, ExcInternalError());
  
  const double alpha=std::acos((v0*v1)/std::sqrt(v0.square()*v1.square()));
  const double d_alpha=alpha/(n+1);
  const Point<dim> pm=0.5*(v0+v1);
  
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


#endif

#if deal_II_dimension == 3

template <>
void
HyperBallBoundary<3>::get_intermediate_points_on_quad (
  const Triangulation<3>::quad_iterator &quad,
  std::vector<Point<3> > &points) const
{
  unsigned int m=static_cast<unsigned int> (std::sqrt(static_cast<double>(points.size())));
  Assert(points.size()==m*m, ExcInternalError());

  std::vector<Point<3> > lp3(m);
  std::vector<Point<3> > lp1(m);
  
  get_intermediate_points_on_line(quad->line(3), lp3);
  get_intermediate_points_on_line(quad->line(1), lp1);
  
  std::vector<Point<3> > lps(m);
  for (unsigned int i=0; i<m; ++i)
    {
      get_intermediate_points_between_points(lp3[i], lp1[i], lps);
      
      for (unsigned int j=0; j<m; ++j)
	points[i*m+j]=lps[j];
    }
}

#endif


template <int dim>
void
HyperBallBoundary<dim>::get_intermediate_points_on_quad (
  const typename Triangulation<dim>::quad_iterator &,
  typename std::vector<Point<dim> > &) const
{
  Assert(false, Boundary<dim>::ExcFunctionNotUseful(dim));
}



#if deal_II_dimension == 1

template <>
void
HyperBallBoundary<1>::
get_normals_at_vertices (const Triangulation<1>::face_iterator &,
			 Boundary<1>::FaceVertexNormals &) const
{
  Assert (false, Boundary<1>::ExcFunctionNotUseful(1));
};

#endif


template <int dim>
void
HyperBallBoundary<dim>::
get_normals_at_vertices (const typename Triangulation<dim>::face_iterator &face,
			 typename Boundary<dim>::FaceVertexNormals &face_vertex_normals) const
{
  for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_face; ++vertex)
    face_vertex_normals[vertex] = face->vertex(vertex)-center;
};



template <int dim>
Point<dim>
HyperBallBoundary<dim>::get_center () const 
{
  return center;
};



template <int dim>
double
HyperBallBoundary<dim>::get_radius () const 
{
  Assert(!compute_radius_automatically, ExcRadiusNotSet());
  return radius;
};


/* ---------------------------------------------------------------------- */


template <int dim>
HalfHyperBallBoundary<dim>::HalfHyperBallBoundary (const Point<dim> center,
						   const double     radius) :
		HyperBallBoundary<dim> (center, radius)
{};



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
  if (line_center(0) == center(0))
    return line_center;
  else
    return HyperBallBoundary<dim>::get_new_point_on_line (line);
};



#if deal_II_dimension == 1

template <>
Point<1>
HalfHyperBallBoundary<1>::
get_new_point_on_quad (const Triangulation<1>::quad_iterator &) const 
{
  Assert (false, ExcInternalError());
  return Point<1>();
};

#endif



template <int dim>
Point<dim>
HalfHyperBallBoundary<dim>::
get_new_point_on_quad (const typename Triangulation<dim>::quad_iterator &quad) const 
{
  const Point<dim> quad_center = quad->center();
  if (quad_center(0) == center(0))
    return quad_center;
  else
    return HyperBallBoundary<dim>::get_new_point_on_quad (quad);
};



template <int dim>
void
HalfHyperBallBoundary<dim>::
get_intermediate_points_on_line (const typename Triangulation<dim>::line_iterator &line,
				 typename std::vector<Point<dim> > &points) const
{
				   // check whether center of object is
				   // at x==0, since then it belongs
				   // to the plane part of the
				   // boundary
  const Point<dim> line_center = line->center();
  if (line_center(0) == center(0))
    return StraightBoundary<dim>::get_intermediate_points_on_line (line, points);
  else
    return HyperBallBoundary<dim>::get_intermediate_points_on_line (line, points);
};



template <int dim>
void
HalfHyperBallBoundary<dim>::
get_intermediate_points_on_quad (const typename Triangulation<dim>::quad_iterator &quad,
				 typename std::vector<Point<dim> > &points) const
{
				   // check whether center of object is
				   // at x==0, since then it belongs
				   // to the plane part of the
				   // boundary
  const Point<dim> quad_center = quad->center();
  if (quad_center(0) == center(0))
    StraightBoundary<dim>::get_intermediate_points_on_quad (quad, points);
  else
    HyperBallBoundary<dim>::get_intermediate_points_on_quad (quad, points);
};



#if deal_II_dimension == 1

template <>
void
HalfHyperBallBoundary<1>::
get_intermediate_points_on_quad (const Triangulation<1>::quad_iterator &,
				 std::vector<Point<1> > &) const
{
  Assert (false, ExcInternalError());
};

#endif


#if deal_II_dimension == 1

template <>
void
HalfHyperBallBoundary<1>::
get_normals_at_vertices (const Triangulation<1>::face_iterator &,
			 Boundary<1>::FaceVertexNormals &) const
{
  Assert (false, Boundary<1>::ExcFunctionNotUseful(1));
};

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
  if (quad_center(0) == center(0))
    StraightBoundary<dim>::get_normals_at_vertices (face, face_vertex_normals);
  else
    HyperBallBoundary<dim>::get_normals_at_vertices (face, face_vertex_normals);
};


/* ---------------------------------------------------------------------- */



template <int dim>
HyperShellBoundary<dim>::HyperShellBoundary (const Point<dim> &center) :
		HyperBallBoundary<dim>(center, 0.)
{
  compute_radius_automatically=true;
};


/* ---------------------------------------------------------------------- */




template <int dim>
HalfHyperShellBoundary<dim>::HalfHyperShellBoundary (const Point<dim> &center) :
		HyperShellBoundary<dim> (center) 
{};



template <int dim>
Point<dim>
HalfHyperShellBoundary<dim>::
get_new_point_on_line (const typename Triangulation<dim>::line_iterator &line) const 
{
				   // first check whether the two end
				   // points of the line are on the
				   // axis of symmetry. if so, then
				   // return the mid point
  if ((line->vertex(0)(0) == center(0)) &&
      (line->vertex(1)(0) == center(0)))
    return (line->vertex(0) + line->vertex(1))/2;
  

				   // otherwise we are on the outer or
				   // inner part of the shell. proceed
				   // as in the base class
  return HyperShellBoundary<dim>::get_new_point_on_line (line);
};



#if deal_II_dimension == 1

template <>
Point<1>
HalfHyperShellBoundary<1>::
get_new_point_on_quad (const Triangulation<1>::quad_iterator &) const
{
  Assert (false, ExcInternalError());
  return Point<1>();
};

#endif



template <int dim>
Point<dim>
HalfHyperShellBoundary<dim>::
get_new_point_on_quad (const typename Triangulation<dim>::quad_iterator &quad) const
{
				   // same thing as for the new point
				   // on the line
  if ((quad->vertex(0)(0) == center(0)) &&
      (quad->vertex(1)(0) == center(0)) &&
      (quad->vertex(2)(0) == center(0)) &&
      (quad->vertex(3)(0) == center(0)))
    return (quad->vertex(0) + quad->vertex(1) +
	    quad->vertex(2) + quad->vertex(3)   )/4;
  

				   // otherwise we are on the outer or
				   // inner part of the shell. proceed
				   // as in the base class
  return HyperShellBoundary<dim>::get_new_point_on_quad (quad);
};



template <int dim>
void
HalfHyperShellBoundary<dim>::
get_intermediate_points_on_line (const typename Triangulation<dim>::line_iterator &line,
				 typename std::vector<Point<dim> > &points) const
{
				   // check whether center of object is
				   // at x==0, since then it belongs
				   // to the plane part of the
				   // boundary
  const Point<dim> line_center = line->center();
  if (line_center(0) == center(0))
    return StraightBoundary<dim>::get_intermediate_points_on_line (line, points);
  else
    return HyperShellBoundary<dim>::get_intermediate_points_on_line (line, points);
};



template <int dim>
void
HalfHyperShellBoundary<dim>::
get_intermediate_points_on_quad (const typename Triangulation<dim>::quad_iterator &quad,
				 typename std::vector<Point<dim> > &points) const
{
				   // check whether center of object is
				   // at x==0, since then it belongs
				   // to the plane part of the
				   // boundary
  const Point<dim> quad_center = quad->center();
  if (quad_center(0) == center(0))
    StraightBoundary<dim>::get_intermediate_points_on_quad (quad, points);
  else
    HyperShellBoundary<dim>::get_intermediate_points_on_quad (quad, points);
};



#if deal_II_dimension == 1

template <>
void
HalfHyperShellBoundary<1>::
get_intermediate_points_on_quad (const Triangulation<1>::quad_iterator &,
				 std::vector<Point<1> > &) const
{
  Assert (false, ExcInternalError());
};

#endif



#if deal_II_dimension == 1

template <>
void
HalfHyperShellBoundary<1>::
get_normals_at_vertices (const Triangulation<1>::face_iterator &,
			 Boundary<1>::FaceVertexNormals &) const
{
  Assert (false, Boundary<1>::ExcFunctionNotUseful(1));
};

#endif



template <int dim>
void
HalfHyperShellBoundary<dim>::
get_normals_at_vertices (const typename Triangulation<dim>::face_iterator &face,
			 typename Boundary<dim>::FaceVertexNormals &face_vertex_normals) const
{
  if (face->center()(0) == center(0))
    StraightBoundary<dim>::get_normals_at_vertices (face, face_vertex_normals);
  else
    HyperShellBoundary<dim>::get_normals_at_vertices (face, face_vertex_normals);
};



// explicit instantiations
template class HyperBallBoundary<deal_II_dimension>;
template class CylinderBoundary<deal_II_dimension>;
template class HalfHyperBallBoundary<deal_II_dimension>;
template class HyperShellBoundary<deal_II_dimension>;
template class HalfHyperShellBoundary<deal_II_dimension>;
