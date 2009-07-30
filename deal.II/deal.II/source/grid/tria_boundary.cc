//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <base/tensor.h>
#include <grid/tria_boundary.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <cmath>

DEAL_II_NAMESPACE_OPEN



/* -------------------------- Boundary --------------------- */


template <int dim, int spacedim>
Boundary<dim, spacedim>::~Boundary ()
{}



template <int dim, int spacedim>
Point<spacedim>
Boundary<dim, spacedim>::
get_new_point_on_quad (const typename Triangulation<dim, spacedim>::quad_iterator &) const 
{
  Assert (false, ExcPureFunctionCalled());
  return Point<spacedim>();
}



template <int dim, int spacedim>
void
Boundary<dim, spacedim>::
get_intermediate_points_on_line (const typename Triangulation<dim, spacedim>::line_iterator &,
				 std::vector<Point<spacedim> > &) const
{
  Assert (false, ExcPureFunctionCalled());
}



template <int dim, int spacedim>
void
Boundary<dim, spacedim>::
get_intermediate_points_on_quad (const typename Triangulation<dim, spacedim>::quad_iterator &,
				 std::vector<Point<spacedim> > &) const
{
  Assert (false, ExcPureFunctionCalled());
}


#if deal_II_dimension > 1

template <int dim, int spacedim>
Point<spacedim>
Boundary<dim,spacedim>::
get_new_point_on_face (const typename Triangulation<dim,spacedim>::face_iterator &face) const
{
  Assert (dim>1, ExcImpossibleInDim(dim));

  switch (dim)
    {
      case 2:
	    return get_new_point_on_line (face);
      case 3:
	    return get_new_point_on_quad (face);
    }
  
  return Point<spacedim>();
}


template <int dim, int spacedim>
void
Boundary<dim,spacedim>::
get_intermediate_points_on_face (const typename Triangulation<dim,spacedim>::face_iterator &face,
				 std::vector<Point<spacedim> > &points) const
{
  Assert (dim>1, ExcImpossibleInDim(dim));

  switch (dim)
    {
      case 2:
	    get_intermediate_points_on_line (face, points);
	    break;
      case 3:
	    get_intermediate_points_on_quad (face, points);
	    break;
      default:
	    Assert (false, ExcNotImplemented());
    }
}

#else

template <>
Point<1>
Boundary<1,1>::
get_new_point_on_face (const Triangulation<1,1>::face_iterator &) const
{
  Assert (false, ExcImpossibleInDim(1));
  return Point<1>();
}


template <>
void
Boundary<1,1>::
get_intermediate_points_on_face (const Triangulation<1,1>::face_iterator &,
				 std::vector<Point<1> > &) const
{
  Assert (false, ExcImpossibleInDim(1));
}



template <>
Point<2>
Boundary<1,2>::
get_new_point_on_face (const Triangulation<1,2>::face_iterator &) const
{
  Assert (false, ExcImpossibleInDim(1));
  return Point<2>();
}


template <>
void
Boundary<1,2>::
get_intermediate_points_on_face (const Triangulation<1,2>::face_iterator &,
				 std::vector<Point<2> > &) const
{
  Assert (false, ExcImpossibleInDim(1));
}



#endif

template <int dim, int spacedim>
void
Boundary<dim, spacedim>::
get_normals_at_vertices (const typename Triangulation<dim, spacedim>::face_iterator &,
			 FaceVertexNormals                                &) const
{
  Assert (false, ExcPureFunctionCalled());
}



template <int dim, int spacedim>
Point<spacedim>
Boundary<dim, spacedim>::
project_to_surface (const typename Triangulation<dim, spacedim>::line_iterator &,
		    const Point<spacedim>                                &trial_point) const
{
  if (spacedim <= 1)
    return trial_point;
  else
    {
      Assert (false, ExcPureFunctionCalled());
      return Point<spacedim>();
    }
}



template <int dim, int spacedim>
Point<spacedim>
Boundary<dim, spacedim>::
project_to_surface (const typename Triangulation<dim, spacedim>::quad_iterator &,
		    const Point<spacedim>                                &trial_point) const
{
  if (spacedim <= 2)
    return trial_point;
  else
    {
      Assert (false, ExcPureFunctionCalled());
      return Point<spacedim>();
    }
}



template <int dim, int spacedim>
Point<spacedim>
Boundary<dim, spacedim>::
project_to_surface (const typename Triangulation<dim, spacedim>::hex_iterator &,
		    const Point<spacedim>                                &trial_point) const
{
  if (spacedim <= 3)
    return trial_point;
  else
    {
      Assert (false, ExcPureFunctionCalled());
      return Point<spacedim>();
    }
}


/* -------------------------- StraightBoundary --------------------- */


template <int dim, int spacedim>
StraightBoundary<dim, spacedim>::StraightBoundary ()
{}


template <int dim, int spacedim>
Point<spacedim>
StraightBoundary<dim, spacedim>::
get_new_point_on_line (const typename Triangulation<dim, spacedim>::line_iterator &line) const 
{
  return (line->vertex(0) + line->vertex(1)) / 2;
}


namespace
{
				   // compute the new midpoint of a quad --
				   // either of a 2d cell on a manifold in 3d
				   // or of a face of a 3d triangulation in 3d
  template <int dim>
  Point<3>
  compute_new_point_on_quad (const typename Triangulation<dim, 3>::quad_iterator &quad)
  {
				     // generate a new point in the middle of
				     // the face based on the points on the
				     // edges and the vertices.
				     //
				     // there is a pathological situation when
				     // this face is on a straight boundary, but
				     // one of its edges and the face behind it
				     // are not; if that face is refined first,
				     // the new point in the middle of that edge
				     // may not be at the same position as
				     // quad->line(.)->center() would have been,
				     // but would have been moved to the
				     // non-straight boundary. We cater to that
				     // situation by using existing edge
				     // midpoints if available, or center() if
				     // not
				     //
				     // note that this situation can not happen
				     // during mesh refinement, as there the
				     // edges are refined first and only then
				     // the face. thus, the check whether a line
				     // has children does not lead to the
				     // situation where the new face midpoints
				     // have different positions depending on
				     // which of the two cells is refined first.
				     //
				     // the situation where the edges aren't
				     // refined happens when the a higher order
				     // MappingQ requests the midpoint of a
				     // face, though, and it is for these cases
				     // that we need to have the check available
    return (quad->vertex(0) + quad->vertex(1) +
	    quad->vertex(2) + quad->vertex(3) +
	    (quad->line(0)->has_children() ?
	     quad->line(0)->child(0)->vertex(1) :
	     quad->line(0)->center()) +
	    (quad->line(1)->has_children() ?
	     quad->line(1)->child(0)->vertex(1) :
	     quad->line(1)->center()) +
	    (quad->line(2)->has_children() ?
	     quad->line(2)->child(0)->vertex(1) :
	     quad->line(2)->center()) +
	    (quad->line(3)->has_children() ?
	     quad->line(3)->child(0)->vertex(1) :
	     quad->line(3)->center())               ) / 8;
  }
}


#if deal_II_dimension < 3

template <int dim, int spacedim>
Point<spacedim>
StraightBoundary<dim, spacedim>::
get_new_point_on_quad (const typename Triangulation<dim, spacedim>::quad_iterator &) const 
{
  Assert (false, ExcImpossibleInDim(dim));
  return Point<spacedim>();
}

#if deal_II_dimension == 2
template <>
Point<3>
StraightBoundary<2,3>::
get_new_point_on_quad (const Triangulation<2,3>::quad_iterator &quad) const 
{
  return compute_new_point_on_quad<2> (quad);
}

#endif
#else

template <>
Point<3>
StraightBoundary<3>::
get_new_point_on_quad (const Triangulation<3>::quad_iterator &quad) const 
{
  return compute_new_point_on_quad<3> (quad);
}

#endif


#if deal_II_dimension == 1

template <>
void
StraightBoundary<1>::
get_intermediate_points_on_line (const Triangulation<1>::line_iterator &,
				 std::vector<Point<1> > &) const
{
  Assert(false, ExcImpossibleInDim(1));
}


template <>
void
StraightBoundary<1,2>::
get_intermediate_points_on_line (const Triangulation<1,2>::line_iterator &,
				 std::vector<Point<2> > &) const
{
  Assert(false, ExcImpossibleInDim(1));
}
#else


template <int dim, int spacedim>
void
StraightBoundary<dim, spacedim>::
get_intermediate_points_on_line (const typename Triangulation<dim, spacedim>::line_iterator &line,
				 std::vector<Point<spacedim> > &points) const
{
  const unsigned int n=points.size();
  Assert(n>0, ExcInternalError());
  
  const double dx=1./(n+1);
  double x=dx;

  const Point<spacedim> vertices[2] = { line->vertex(0),
				   line->vertex(1) };
  
  for (unsigned int i=0; i<n; ++i, x+=dx)
    points[i] = (1-x)*vertices[0] + x*vertices[1];
}

#endif



#if deal_II_dimension < 3

template <int dim, int spacedim>
void
StraightBoundary<dim, spacedim>::
get_intermediate_points_on_quad (const typename Triangulation<dim, spacedim>::quad_iterator &,
				 std::vector<Point<spacedim> > &) const
{
  Assert(false, ExcImpossibleInDim(dim));
}

#else

template <>
void
StraightBoundary<3>::
get_intermediate_points_on_quad (const Triangulation<3>::quad_iterator &quad,
				 std::vector<Point<3> > &points) const
{
  const unsigned int spacedim = 3;
  
  const unsigned int n=points.size(),
		     m=static_cast<unsigned int>(std::sqrt(static_cast<double>(n)));
				   // is n a square number
  Assert(m*m==n, ExcInternalError());

  const double ds=1./(m+1);
  double y=ds;

  const Point<spacedim> vertices[4] = { quad->vertex(0),
					quad->vertex(1),
					quad->vertex(2),
					quad->vertex(3) };

  for (unsigned int i=0; i<m; ++i, y+=ds)
    {
      double x=ds;
      for (unsigned int j=0; j<m; ++j, x+=ds)
	points[i*m+j]=((1-x) * vertices[0] +
		       x     * vertices[1]) * (1-y) +
		      ((1-x) * vertices[2] +
		       x     * vertices[3]) * y;
    }
}

#endif



#if deal_II_dimension == 1

template <>
void
StraightBoundary<1>::
get_normals_at_vertices (const Triangulation<1>::face_iterator &,
			 Boundary<1>::FaceVertexNormals &) const
{
  Assert (false, ExcImpossibleInDim(1));
}

template <>
void
StraightBoundary<1,2>::
get_normals_at_vertices (const Triangulation<1,2>::face_iterator &,
			 Boundary<1,2>::FaceVertexNormals &) const
{
  Assert (false, ExcImpossibleInDim(1));
}

#endif


#if deal_II_dimension == 2

template <>
void
StraightBoundary<2>::
get_normals_at_vertices (const Triangulation<2>::face_iterator &face,
			 Boundary<2>::FaceVertexNormals &face_vertex_normals) const
{
  const Tensor<1,2> tangent = face->vertex(1) - face->vertex(0);
  for (unsigned int vertex=0; vertex<GeometryInfo<2>::vertices_per_face; ++vertex)
				     // compute normals from tangent
    face_vertex_normals[vertex] = Point<2>(tangent[1],
                                           -tangent[0]);
}

template <>
void
StraightBoundary<2,3>::
get_normals_at_vertices (const Triangulation<2,3>::face_iterator &face,
			 Boundary<2,3>::FaceVertexNormals &face_vertex_normals) const
{
  const Tensor<1,3> tangent = face->vertex(1) - face->vertex(0);
  for (unsigned int vertex=0; vertex<GeometryInfo<2>::vertices_per_face; ++vertex)
				     // compute normals from tangent
    face_vertex_normals[vertex] = Point<3>(tangent[1],
                                           -tangent[0],0);
  Assert(false, ExcNotImplemented());
}

#endif



#if deal_II_dimension == 3

template <>
void
StraightBoundary<3>::
get_normals_at_vertices (const Triangulation<3>::face_iterator &face,
			 Boundary<3>::FaceVertexNormals &face_vertex_normals) const
{
  const unsigned int vertices_per_face = GeometryInfo<3>::vertices_per_face;

  static const unsigned int neighboring_vertices[4][2]=
  { {1,2},{3,0},{0,3},{2,1}};
  for (unsigned int vertex=0; vertex<vertices_per_face; ++vertex)
    {
				       // first define the two tangent
				       // vectors at the vertex by
				       // using the two lines
				       // radiating away from this
				       // vertex
      const Tensor<1,3> tangents[2]
	= { face->vertex(neighboring_vertices[vertex][0])
	      - face->vertex(vertex),
	    face->vertex(neighboring_vertices[vertex][1])
	      - face->vertex(vertex)      };

				       // then compute the normal by
				       // taking the cross
				       // product. since the normal is
				       // not required to be
				       // normalized, no problem here
      cross_product (face_vertex_normals[vertex],
		     tangents[0], tangents[1]);
    };
}

#endif


template <int dim, int spacedim>
Point<spacedim>
StraightBoundary<dim, spacedim>::
project_to_surface (const typename Triangulation<dim, spacedim>::line_iterator &line,
		    const Point<spacedim>                                &trial_point) const
{
  if (spacedim <= 1)
    return trial_point;
  else
    {
				       // find the point that lies on
				       // the line p1--p2. the
				       // formulas pan out to
				       // something rather simple
				       // because the mapping to the
				       // line is linear
      const Point<spacedim> p1 = line->vertex(0),
			    p2 = line->vertex(1);
      const double s = (trial_point-p1)*(p2-p1) / ((p2-p1)*(p2-p1));
      return p1 + s*(p2-p1);
    }
}



template <int dim, int spacedim>
Point<spacedim>
StraightBoundary<dim, spacedim>::
project_to_surface (const typename Triangulation<dim, spacedim>::quad_iterator &,
		    const Point<spacedim>                                &trial_point) const
{
  if (spacedim <= 2)
    return trial_point;
  else
    {
      Assert (false, ExcNotImplemented());
      return Point<spacedim>();
    }
}



template <int dim, int spacedim>
Point<spacedim>
StraightBoundary<dim, spacedim>::
project_to_surface (const typename Triangulation<dim, spacedim>::hex_iterator &,
		    const Point<spacedim>                                &trial_point) const
{
  if (spacedim <= 3)
    return trial_point;
  else
    {
      Assert (false, ExcNotImplemented());
      return Point<spacedim>();
    }
}



// explicit instantiations
template class Boundary<deal_II_dimension>;
template class StraightBoundary<deal_II_dimension>;

#if deal_II_dimension == 1 || deal_II_dimension == 2
template class Boundary<deal_II_dimension, deal_II_dimension+1>;
template class StraightBoundary<deal_II_dimension, deal_II_dimension+1>;
#endif

DEAL_II_NAMESPACE_CLOSE

