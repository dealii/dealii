//----------------------------  mapping_c1.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mapping_c1.cc  ---------------------------


#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary.h>
#include <fe/mapping_c1.h>
#include <cmath>




template <int dim>
MappingC1<dim>::MappingC1 ()
		:
		MappingQ<dim> (3)
{
  Assert (dim > 1, ExcImpossibleInDim(dim));
};


#if deal_II_dimension == 1

template <>
void
MappingC1<1>::add_line_support_points (const Triangulation<1>::cell_iterator &,
				       std::vector<Point<1> > &) const
{
  const unsigned int dim = 1;
  Assert (dim > 1, ExcImpossibleInDim(dim));
};

#endif


#if deal_II_dimension == 2

template <>
void
MappingC1<2>::add_line_support_points (const Triangulation<2>::cell_iterator &cell,
				       std::vector<Point<2> > &a) const
{
  const unsigned int dim = 2;
  std::vector<Point<dim> > line_points (2);
      
				       // loop over each of the lines,
				       // and if it is at the
				       // boundary, then first get the
				       // boundary description and
				       // second compute the points on
				       // it. if not at the boundary,
				       // get the respective points
				       // from another function
  for (unsigned int line_no=0; line_no<GeometryInfo<dim>::lines_per_cell; ++line_no)
    {
      const Triangulation<dim>::line_iterator line = cell->line(line_no);

      if (line->at_boundary())
	{
					   // first get the normal
					   // vectors at the two
					   // vertices of this line
					   // from the boundary
					   // description
	  const Boundary<dim> &boundary
	    = line->get_triangulation().get_boundary(line->boundary_indicator());

	  Boundary<dim>::FaceVertexNormals face_vertex_normals;
	  boundary.get_normals_at_vertices (line, face_vertex_normals);

					   // then transform them into
					   // interpolation points for
					   // a cubic polynomial
					   //
					   // for this, note that if
					   // we describe the boundary
					   // curve as a polynomial in
					   // tangential coordinate
					   // @p{t=0..1} (along the
					   // line) and @p{s} in
					   // normal direction, then
					   // the cubic mapping is
					   // such that @p{s = a*t**3
					   // + b*t**2 + c*t + d}, and
					   // we want to determine the
					   // interpolation points at
					   // @p{t=1/3} and
					   // @p{t=2/3}. Since at
					   // @p{t=0,1} we want a
					   // vertex which is actually
					   // at the boundary, we know
					   // that @p{d=0} and
					   // @p{a=-b-c}. As
					   // side-conditions, we want
					   // that the derivatives at
					   // @p{t=0} and @p{t=1},
					   // i.e. at the vertices
					   // match those returned by
					   // the boundary. We then
					   // have that
					   // @p{s(1/3)=1/27(2b+8c)}
					   // and
					   // @p{s(2/3)=4/27b+10/27c}.
					   //
					   // The task is then first
					   // to determine the
					   // coefficients from the
					   // tangentials. for that,
					   // first rotate the
					   // tangents of @p{s(t)}
					   // into the global
					   // coordinate system. they
					   // are @p{A (1,c)} and @p{A
					   // (1,-b-2c)} with @p{A} the
					   // rotation matrix, since
					   // the tangentials in the
					   // coordinate system
					   // relative to the line are
					   // @p{(1,c)} and @p{(1,-b-2c)}
					   // at the two vertices,
					   // respectively. We then
					   // have to make sure by
					   // matching @p{b,c} that
					   // these tangentials are
					   // orthogonal to the normals
					   // returned by the boundary
					   // object
	  const Tensor<1,2> coordinate_vector = line->vertex(1) - line->vertex(0);
	  const double      h                 = std::sqrt(coordinate_vector * coordinate_vector);
	  Tensor<1,2> coordinate_axis = coordinate_vector;
	  coordinate_axis /= h;

	  const double alpha = std::atan2(coordinate_axis[1], coordinate_axis[0]);
	  const double c = -((face_vertex_normals[0][1] * std::sin(alpha)
			      +face_vertex_normals[0][0] * std::cos(alpha)) /
			     (face_vertex_normals[0][1] * std::cos(alpha)
			      -face_vertex_normals[0][0] * std::sin(alpha)));
	  const double b = ((face_vertex_normals[1][1] * std::sin(alpha)
			      +face_vertex_normals[1][0] * std::cos(alpha)) /
			     (face_vertex_normals[1][1] * std::cos(alpha)
			      -face_vertex_normals[1][0] * std::sin(alpha)))
			   -2*c;
		       

					   // next evaluate the so
					   // determined cubic
					   // polynomial at the points
					   // 1/3 and 2/3, first in
					   // unit coordinates
	  const Point<2> new_unit_points[2] = { Point<2>(1./3., 1./27.*(2*b+8*c)),
						  Point<2>(2./3., 4./27.*b+10./27.*c) };
					   // then transform these
					   // points to real
					   // coordinates by rotating,
					   // scaling and shifting
	  for (unsigned int i=0; i<2; ++i)
	    {
	      Point<2> real_point (std::cos(alpha) * new_unit_points[i][0]
				   - std::sin(alpha) * new_unit_points[i][1],
				   std::sin(alpha) * new_unit_points[i][0]
				   + std::cos(alpha) * new_unit_points[i][1]);
	      real_point *= h;
	      real_point += line->vertex(0);
	      a.push_back (real_point);
	    };
	}
      else
					 // not at boundary
	{
	  static const StraightBoundary<dim> straight_boundary;
	  straight_boundary.get_intermediate_points_on_line (line, line_points);
	  a.insert (a.end(), line_points.begin(), line_points.end());
	};
    };
};

#endif



template <int dim>
void
MappingC1<dim>::add_line_support_points (const typename Triangulation<dim>::cell_iterator &,
					 typename std::vector<Point<dim> > &) const
{
  Assert (false, ExcNotImplemented());
};




#if deal_II_dimension == 1

template <>
void
MappingC1<1>::add_quad_support_points (const Triangulation<1>::cell_iterator &,
				       std::vector<Point<1> > &) const
{
  const unsigned int dim = 1;
  Assert (dim > 2, ExcImpossibleInDim(dim));
};

#endif



#if deal_II_dimension == 2

template <>
void
MappingC1<2>::add_quad_support_points (const Triangulation<2>::cell_iterator &,
				       std::vector<Point<2> > &) const
{
  const unsigned int dim = 2;
  Assert (dim > 2, ExcImpossibleInDim(dim));
};

#endif




template <int dim>
void
MappingC1<dim>::add_quad_support_points (const typename Triangulation<dim>::cell_iterator &,
					 typename std::vector<Point<dim> > &) const
{
  Assert (false, ExcNotImplemented());
};




// explicit instantiations
template class MappingC1<deal_II_dimension>;
