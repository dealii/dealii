//----------------------------  grid_tools.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  grid_tools.cc  ---------------------------



#include <grid/grid_tools.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_dof_accessor.h>

#include <cmath>

#if deal_II_dimension != 1

template <int dim>
double
GridTools::diameter (const Triangulation<dim> &tria)
{
				   // the algorithm used simply
				   // traverses all cells and picks
				   // out the boundary vertices. it
				   // may or may not be faster to
				   // simply get all vectors, don't
				   // mark boundary vertices, and
				   // compute the distances thereof,
				   // but at least as the mesh is
				   // refined, it seems better to
				   // first mark boundary nodes, as
				   // marking is O(N) in the number of
				   // cells/vertices, while computing
				   // the maximal distance is O(N*N)
  const std::vector<Point<dim> > &vertices = tria.get_vertices ();
  std::vector<bool> boundary_vertices (vertices.size(), false);

  typename Triangulation<dim>::active_cell_iterator
    cell = tria.begin_active();
  const typename Triangulation<dim>::active_cell_iterator
    endc = tria.end();
  for (; cell!=endc; ++cell)
    for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
      if (cell->face(face)->at_boundary ())
	for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_face; ++i)
	  boundary_vertices[cell->face(face)->vertex_index(i)] = true;

				   // now traverse the list of
				   // boundary vertices and check
				   // distances. since distances are
				   // symmetric, we only have to check
				   // one half
  double max_distance_sqr = 0;
  std::vector<bool>::const_iterator pi = boundary_vertices.begin();
  const unsigned int N = boundary_vertices.size();
  for (unsigned int i=0; i<N; ++i, ++pi)
    {
      std::vector<bool>::const_iterator pj = pi+1;
      for (unsigned int j=i+1; j<N; ++j, ++pj)
	if ((*pi==true) && (*pj==true) &&
	    ((vertices[i]-vertices[j]).square() > max_distance_sqr))
	  max_distance_sqr = (vertices[i]-vertices[j]).square();
    };

  return std::sqrt(max_distance_sqr);
}


#else

double
GridTools::diameter (const Triangulation<1> &tria)
{
				   // for 1d, simply check the
				   // vertices of the left- and
				   // rightmost coarse grid cell
  Triangulation<1>::cell_iterator leftmost  = tria.begin(0);
  Triangulation<1>::cell_iterator rightmost = tria.begin(0);

  while (!leftmost->at_boundary(0))  leftmost  = leftmost->neighbor(0);
  while (!rightmost->at_boundary(1)) rightmost = rightmost->neighbor(1);

  return std::sqrt((leftmost->vertex(0) - rightmost->vertex(1)).square());
}

#endif



// define some transformations in an anonymous namespace
namespace 
{
  template <int dim>
  class ShiftPoint
  {
    public:
      ShiftPoint (const Point<dim> &shift)
		      :
		      shift(shift)
	{};
      Point<dim> operator() (const Point<dim> p) const
	{
	  return p+shift;
	};
    private:
      const Point<dim> shift;
  };


                                   // the following class is only
                                   // needed in 2d, so avoid trouble
                                   // with compilers warning otherwise
#if deal_II_dimension == 2
  class Rotate2d
  {
    public:
      Rotate2d (const double angle)
		      :
		      angle(angle)
	{};
      Point<2> operator() (const Point<2> p) const
	{
	  return Point<2> (std::cos(angle)*p(0) - std::sin(angle) * p(1),
			   std::sin(angle)*p(0) + std::cos(angle) * p(1));
	};
    private:
      const double angle;
  };
#endif


  template <int dim>
  class ScalePoint
  {
    public:
      ScalePoint (const double factor)
		      :
		      factor(factor)
	{};
      Point<dim> operator() (const Point<dim> p) const
	{
	  return p*factor;
	};
    private:
      const double factor;
  };
}


template <int dim>
void
GridTools::shift (const Point<dim>   &shift_vector,
		  Triangulation<dim> &triangulation)
{
  transform (ShiftPoint<dim>(shift_vector), triangulation);
}


#if deal_II_dimension == 2

void
GridTools::rotate (const double      angle,
		   Triangulation<2> &triangulation)
{
  transform (Rotate2d(angle), triangulation);
}

#endif


template <int dim>
void
GridTools::scale (const double        scaling_factor,
		  Triangulation<dim> &triangulation)
{
  Assert (scaling_factor>0, ExcScalingFactorNotPositive (scaling_factor));
  transform (ScalePoint<dim>(scaling_factor), triangulation);
}



template <int dim, typename Container>
typename Container::active_cell_iterator
GridTools::find_active_cell_around_point (const Container  &container,
                                          const Point<dim> &p)
{
                                   // first find the coarse grid cell
                                   // that contains the point. we can
                                   // only do this by a linear search
  typename Container::cell_iterator cell = container.begin(0);
  for (; cell!=container.end(0); ++cell)
    if (cell->point_inside (p))
      break;

                                   // make sure that we found a cell
                                   // in the coarse grid that contains
                                   // this point. for cases where this
                                   // might happen unexpectedly, see
                                   // the documentation of this
                                   // function
  AssertThrow (cell != container.end(0),
               ExcPointNotFoundInCoarseGrid<dim> (p));

                                   // now do the logarithmic part of
                                   // the algorithm: go from child to
                                   // grandchild
  while (cell->has_children())
    {
      unsigned int c=0;
      for (; c<GeometryInfo<dim>::children_per_cell; ++c)
        if (cell->child(c)->point_inside (p))
          break;

                                       // make sure we found a child
                                       // cell
      AssertThrow (c != GeometryInfo<dim>::children_per_cell,
                   ExcPointNotFound<dim> (p));

                                       // then reset cell to the child
      cell = cell->child(c);
    }

                                   // now that we have a terminal
                                   // cell, return it
  return cell;
}




#if deal_II_dimension != 1
template
double
GridTools::diameter<deal_II_dimension> (const Triangulation<deal_II_dimension> &);
#endif

template
void GridTools::shift<deal_II_dimension> (const Point<deal_II_dimension> &,
					  Triangulation<deal_II_dimension> &);

template
void GridTools::scale<deal_II_dimension> (const double,
					  Triangulation<deal_II_dimension> &);

template
Triangulation<deal_II_dimension>::active_cell_iterator
GridTools::find_active_cell_around_point (const Triangulation<deal_II_dimension> &,
                                          const Point<deal_II_dimension> &p);

template
DoFHandler<deal_II_dimension>::active_cell_iterator
GridTools::find_active_cell_around_point (const DoFHandler<deal_II_dimension> &,
                                          const Point<deal_II_dimension> &p);

template
MGDoFHandler<deal_II_dimension>::active_cell_iterator
GridTools::find_active_cell_around_point (const MGDoFHandler<deal_II_dimension> &,
                                          const Point<deal_II_dimension> &p);
