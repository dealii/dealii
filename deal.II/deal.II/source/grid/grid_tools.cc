//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002, 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------



#include <grid/grid_tools.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <lac/sparsity_pattern.h>
#include <lac/compressed_sparsity_pattern.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_dgq.h>
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



#if deal_II_dimension == 3

template <>
double
GridTools::cell_measure(const std::vector<Point<3> > &all_vertices,
			const int vertex_indices[GeometryInfo<3>::vertices_per_cell])
{
  const double x[8] = { all_vertices[vertex_indices[0]](0),
			all_vertices[vertex_indices[1]](0),
			all_vertices[vertex_indices[2]](0),
			all_vertices[vertex_indices[3]](0),
			all_vertices[vertex_indices[4]](0),
			all_vertices[vertex_indices[5]](0),
			all_vertices[vertex_indices[6]](0),
			all_vertices[vertex_indices[7]](0)   };
  const double y[8] = { all_vertices[vertex_indices[0]](1),
			all_vertices[vertex_indices[1]](1),
			all_vertices[vertex_indices[2]](1),
			all_vertices[vertex_indices[3]](1),
			all_vertices[vertex_indices[4]](1),
			all_vertices[vertex_indices[5]](1),
			all_vertices[vertex_indices[6]](1),
			all_vertices[vertex_indices[7]](1)  };
  const double z[8] = { all_vertices[vertex_indices[0]](2),
			all_vertices[vertex_indices[1]](2),
			all_vertices[vertex_indices[2]](2),
			all_vertices[vertex_indices[3]](2),
			all_vertices[vertex_indices[4]](2),
			all_vertices[vertex_indices[5]](2),
			all_vertices[vertex_indices[6]](2),
			all_vertices[vertex_indices[7]](2)  };

/*
  Get the computation of the measure by the same little Maple script
  as in TriaObjectAccessor<3, 3>::barycenter(), see tria_accessor.cc.
*/
  double s1, s2, s3;
  
  s3 = -x[2]*y[7]*z[3]/12.0+x[1]*y[4]*z[0]/12.0-x[0]*y[3]*z[7]/12.0-z[1]*x
       [2]*y[0]/12.0+x[5]*y[4]*z[0]/12.0+x[1]*y[5]*z[0]/12.0-x[1]*y[2]*z[0]/12.0-y[2]*
       x[3]*z[7]/12.0-x[1]*y[0]*z[4]/12.0-y[1]*x[0]*z[2]/12.0-y[0]*x[7]*z[3]/12.0+y[5]
       *x[4]*z[7]/12.0-x[0]*y[3]*z[4]/12.0-y[1]*x[0]*z[3]/12.0-z[0]*x[7]*y[4]/12.0+x
       [1]*y[0]*z[3]/12.0-x[2]*y[3]*z[0]/12.0-y[1]*x[2]*z[5]/12.0;
  s2 = s3+y[0]*x[3]*z[4]/12.0+y[1]*x[2]*z[0]/12.0+y[2]*x[7]*z[3]/12.0+z[1]*
       x[4]*y[0]/12.0-y[5]*x[4]*z[0]/12.0-y[1]*x[5]*z[0]/12.0+x[5]*y[7]*z[4]/12.0+y[1]
       *x[3]*z[0]/12.0+x[0]*y[4]*z[7]/12.0+z[0]*x[7]*y[3]/12.0-z[2]*x[3]*y[0]/12.0+y
       [1]*x[0]*z[5]/12.0-x[1]*y[5]*z[2]/12.0+y[1]*x[0]*z[4]/12.0+z[1]*x[2]*y[5]/12.0-
       x[5]*y[0]*z[4]/12.0+y[2]*x[3]*z[1]/12.0+x[1]*y[0]*z[2]/12.0;
  s3 = -z[1]*x[0]*y[4]/12.0+z[1]*x[0]*y[2]/12.0-x[0]*y[7]*z[4]/12.0+z[0]*x
       [4]*y[7]/12.0+x[2]*y[3]*z[7]/12.0+y[2]*x[6]*z[3]/12.0+y[1]*x[6]*z[2]/12.0-y[2]*
       x[1]*z[3]/12.0+x[1]*y[2]*z[5]/12.0+y[0]*x[3]*z[7]/12.0+y[2]*x[3]*z[0]/12.0-y[0]
       *x[4]*z[3]/12.0-y[2]*x[0]*z[3]/12.0-y[1]*x[4]*z[0]/12.0+z[1]*x[5]*y[0]/12.0+x
       [1]*y[2]*z[6]/12.0-x[1]*y[0]*z[5]/12.0-x[2]*y[3]*z[1]/12.0;
  s1 = s3-z[0]*x[3]*y[4]/12.0+z[1]*x[0]*y[3]/12.0-z[1]*x[5]*y[2]/12.0+z[1]*
       x[2]*y[6]/12.0-z[2]*x[3]*y[1]/12.0+z[2]*x[3]*y[6]/12.0+x[2]*y[1]*z[3]/12.0-z[1]
       *x[3]*y[0]/12.0+x[2]*y[3]*z[6]/12.0-z[0]*x[3]*y[7]/12.0+y[5]*x[0]*z[4]/12.0-z
       [1]*x[6]*y[2]/12.0-z[2]*x[6]*y[3]/12.0+z[2]*x[0]*y[3]/12.0+z[2]*x[3]*y[7]/12.0+
       z[2]*x[1]*y[3]/12.0+y[1]*x[5]*z[2]/12.0-z[2]*x[7]*y[3]/12.0+s2;
  s3 = x[0]*y[7]*z[3]/12.0-z[1]*x[0]*y[5]/12.0-x[1]*y[3]*z[0]/12.0-x[1]*y
       [6]*z[2]/12.0-x[2]*y[6]*z[3]/12.0-x[5]*y[4]*z[7]/12.0+z[0]*x[4]*y[3]/12.0+x[0]*
       y[4]*z[3]/12.0-x[6]*y[5]*z[7]/12.0-x[6]*y[7]*z[2]/12.0+z[1]*x[6]*y[5]/12.0+y[6]
       *x[4]*z[7]/12.0+y[1]*x[5]*z[6]/12.0+x[1]*y[6]*z[5]/12.0-y[5]*x[7]*z[4]/12.0+y
       [0]*x[7]*z[4]/12.0-y[0]*x[4]*z[7]/12.0-z[5]*x[0]*y[4]/12.0;
  s2 = s3+z[6]*x[3]*y[7]/12.0+z[6]*x[7]*y[4]/12.0-z[6]*x[4]*y[7]/12.0-z[6]*
       x[7]*y[3]/12.0+y[6]*x[5]*z[7]/12.0-y[6]*x[7]*z[5]/12.0-x[2]*y[5]*z[6]/12.0+z[6]
       *x[2]*y[7]/12.0+z[6]*x[7]*y[5]/12.0-z[6]*x[5]*y[7]/12.0-z[6]*x[7]*y[2]/12.0+x
       [6]*y[7]*z[5]/12.0+z[5]*x[7]*y[4]/12.0-z[5]*x[4]*y[7]/12.0+y[5]*x[1]*z[4]/12.0-
       y[5]*x[6]*z[4]/12.0-y[5]*x[4]*z[1]/12.0+y[5]*x[4]*z[6]/12.0;
  s3 = y[6]*x[7]*z[3]/12.0+x[6]*y[5]*z[2]/12.0-x[5]*y[6]*z[2]/12.0-x[6]*y
       [2]*z[5]/12.0+x[5]*y[2]*z[6]/12.0+x[6]*y[2]*z[7]/12.0+z[5]*x[4]*y[0]/12.0-y[6]*
       x[2]*z[7]/12.0+x[2]*y[6]*z[5]/12.0+y[6]*x[7]*z[2]/12.0-y[6]*x[3]*z[7]/12.0-x[7]
       *y[4]*z[3]/12.0-x[6]*y[4]*z[7]/12.0-x[5]*y[1]*z[4]/12.0+x[6]*y[7]*z[4]/12.0+x
       [7]*y[3]*z[4]/12.0-z[1]*x[5]*y[6]/12.0+x[3]*y[4]*z[7]/12.0+x[5]*y[6]*z[4]/12.0;

				   // the vertices entered into the
				   // above field and the vertices
				   // used by the maple script use
				   // different orderings (i.e. one is
				   // the inside-out orientation of
				   // the other one). the measure is
				   // thus negative. Thus take the
				   // negative value of the result
  return -(s3+x[5]*y[4]*z[1]/12.0-x[5]*y[4]*z[6]/12.0-x[3]*y[7]*z[4]/12.0+x[6]*
    y[3]*z[7]/12.0-y[1]*x[6]*z[5]/12.0-z[5]*x[1]*y[4]/12.0+z[5]*x[6]*y[4]/12.0+z[5]
    *x[4]*y[1]/12.0-z[5]*x[4]*y[6]/12.0-x[6]*y[7]*z[3]/12.0-x[4]*y[3]*z[7]/12.0+x
    [4]*y[7]*z[3]/12.0-x[1]*y[5]*z[6]/12.0-y[6]*x[7]*z[4]/12.0-y[1]*x[2]*z[6]/12.0-
    y[2]*x[3]*z[6]/12.0+s2+s1+x[2]*y[0]*z[3]/12.0);
}

#else

template <int dim>
double
GridTools::cell_measure(const std::vector<Point<dim> > &all_vertices,
			const int [GeometryInfo<dim>::vertices_per_cell])
{
  Assert(false, ExcNotImplemented());
  return 0.;
}

#endif



// define some transformations in an anonymous namespace
#ifdef DEAL_II_ANON_NAMESPACE_BOGUS_WARNING
  namespace TRANS
#else
  namespace
#endif
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
#ifdef DEAL_II_ANON_NAMESPACE_BOGUS_WARNING
  transform (TRANS::ShiftPoint<dim>(shift_vector), triangulation);
#else
  transform (ShiftPoint<dim>(shift_vector), triangulation);
#endif  
}


#if deal_II_dimension == 2

void
GridTools::rotate (const double      angle,
		   Triangulation<2> &triangulation)
{
#ifdef DEAL_II_ANON_NAMESPACE_BOGUS_WARNING
  transform (TRANS::Rotate2d(angle), triangulation);
#else
  transform (Rotate2d(angle), triangulation);
#endif  
}

#endif


template <int dim>
void
GridTools::scale (const double        scaling_factor,
		  Triangulation<dim> &triangulation)
{
  Assert (scaling_factor>0, ExcScalingFactorNotPositive (scaling_factor));
#ifdef DEAL_II_ANON_NAMESPACE_BOGUS_WARNING
  transform (TRANS::ScalePoint<dim>(scaling_factor), triangulation);
#else
  transform (ScalePoint<dim>(scaling_factor), triangulation);
#endif  
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



template <int dim>
void
GridTools::
partition_triangulation (const unsigned int  n_partitions,
                         Triangulation<dim> &triangulation)
{
  Assert (n_partitions > 0, ExcInvalidNumberOfPartitions(n_partitions));

                                   // check for an easy return
  if (n_partitions == 1)
    {
      for (typename Triangulation<dim>::active_cell_iterator
             cell = triangulation.begin_active();
           cell != triangulation.end(); ++cell)
        cell->set_subdomain_id (0);
      return;
    }

                                   // we decompose the domain by first
                                   // generating the connection graph of all
                                   // cells with their neighbors, and then
                                   // passing this graph off to METIS. To make
                                   // things a little simpler and more
                                   // general, we let the function
                                   // DoFTools:make_flux_sparsity_pattern
                                   // function generate the connection graph
                                   // for us and reuse the SparsityPattern
                                   // data structure for the connection
                                   // graph. The connection structure of the
                                   // mesh is obtained by using a fake
                                   // piecewise constant finite element
                                   //
                                   // Since in 3d the generation of a
                                   // sparsity pattern can be expensive, we
                                   // take the detour of the compressed
                                   // sparsity pattern, which is a little
                                   // slower but more efficient in terms of
                                   // memory
  FE_DGQ<dim> fake_q0(0);
  DoFHandler<dim> dof_handler (triangulation);
  dof_handler.distribute_dofs (fake_q0);
  Assert (dof_handler.n_dofs() == triangulation.n_active_cells(),
          ExcInternalError());
  
  CompressedSparsityPattern csp (dof_handler.n_dofs(),
                                 dof_handler.n_dofs());
  DoFTools::make_flux_sparsity_pattern (dof_handler, csp);
  
  SparsityPattern sparsity_pattern;
  sparsity_pattern.copy_from (csp);

                                   // partition this connection graph and get
                                   // back a vector of indices, one per degree
                                   // of freedom (which is associated with a
                                   // cell)
  std::vector<unsigned int> partition_indices (triangulation.n_active_cells());
  sparsity_pattern.partition (n_partitions,  partition_indices);

                                   // finally loop over all cells and set the
                                   // subdomain ids. for this, get the DoF
                                   // index of each cell and extract the
                                   // subdomain id from the vector obtained
                                   // above
  std::vector<unsigned int> dof_indices(1);
  for (typename DoFHandler<dim>::active_cell_iterator
         cell = dof_handler.begin_active();
       cell != dof_handler.end(); ++cell)
    {
      cell->get_dof_indices(dof_indices);
      Assert (dof_indices[0] < triangulation.n_active_cells(),
              ExcInternalError());
      Assert (partition_indices[dof_indices[0]] < n_partitions,
              ExcInternalError());
      
      cell->set_subdomain_id (partition_indices[dof_indices[0]]);
    }
}



template <int dim>
void
GridTools::
get_subdomain_association (const Triangulation<dim>  &triangulation,
                           std::vector<unsigned int> &subdomain)
{
  Assert (subdomain.size() == triangulation.n_active_cells(),
          ExcDimensionMismatch (subdomain.size(),
                                triangulation.n_active_cells()));
  unsigned int index = 0;
  for (typename Triangulation<dim>::active_cell_iterator
         cell = triangulation.begin_active();
       cell!=triangulation.end(); ++cell, ++index)
    subdomain[index] = cell->subdomain_id();

  Assert (index == subdomain.size(), ExcInternalError());
}



template <int dim>
unsigned int
GridTools::
count_cells_with_subdomain_association (const Triangulation<dim> &triangulation,
                                        const unsigned int        subdomain)
{
  unsigned int count = 0;
  for (typename Triangulation<dim>::active_cell_iterator
         cell = triangulation.begin_active();
       cell!=triangulation.end(); ++cell)
    if (cell->subdomain_id() == subdomain)
      ++count;

  Assert (count != 0, ExcNonExistentSubdomain(subdomain));

  return count;
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

template
void
GridTools::partition_triangulation (const unsigned int,
                                    Triangulation<deal_II_dimension> &);

template
void
GridTools::
get_subdomain_association (const Triangulation<deal_II_dimension>  &,
                           std::vector<unsigned int> &);

template
unsigned int
GridTools::
count_cells_with_subdomain_association (const Triangulation<deal_II_dimension> &,
                                        const unsigned int        );

