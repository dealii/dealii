//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <base/memory_consumption.h>
#include <grid/tria.h>
#include <grid/tria_levels.h>
#include <grid/tria_faces.h>
#include <grid/tria_boundary.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <base/geometry_info.h>
#include <grid/grid_tools.h>
#include <grid/magic_numbers.h>
#include <lac/vector.h>

#include <algorithm>
#include <numeric>
#include <map>
#include <cmath>
#include <iostream>
#include <boost/array.hpp>

DEAL_II_NAMESPACE_OPEN

// initialize the @p{straight_boundary} pointer of the triangulation
// class. for the reasons why it is done like this, see the
// documentation of that member variable
namespace 
{
  const StraightBoundary<deal_II_dimension> dummy_straight_boundary;
}


template <int dim>
const StraightBoundary<dim> *
Triangulation<dim>::straight_boundary = &dummy_straight_boundary;

template <int dim>
const unsigned int
Triangulation<dim>::dimension;


template <int dim>
Triangulation<dim>::Triangulation (const MeshSmoothing smooth_grid) :
		Subscriptor (),
		faces(NULL),
		smooth_grid(smooth_grid)
{
				   // set default boundary for all
				   // possible components
  for (unsigned int i=0;i<255;++i)
    boundary[i] = straight_boundary;
}


template <int dim>
Triangulation<dim>::Triangulation (const Triangulation<dim> &) 
		:
		Subscriptor ()
				   // do not set any subscriptors;
				   // anyway, calling this constructor
				   // is an error!
{
  Assert (false, ExcInternalError());
}


template <int dim>
Triangulation<dim>::~Triangulation ()
{
  for (unsigned int i=0; i<levels.size(); ++i)
    delete levels[i];
  levels.clear ();
  delete faces;
}


template <int dim>
void Triangulation<dim>::clear () 
{
  clear_despite_subscriptions();
}


template <int dim>
void
Triangulation<dim>::set_boundary (const unsigned int number,
				  const Boundary<dim>& boundary_object)
{
  Assert(number<255, ExcIndexRange(number,0,255));
  
  boundary[number] = &boundary_object;
}



template <int dim>
void
Triangulation<dim>::set_boundary (const unsigned int number)
{
  set_boundary (number, *straight_boundary);
}



template <int dim>
const Boundary<dim> &
Triangulation<dim>::get_boundary (const unsigned int number) const 
{
  Assert(number<255, ExcIndexRange(number,0,255));
  
  return *(boundary[number]);
}


#if deal_II_dimension == 1

template <>
std::vector<unsigned char>
Triangulation<1>::get_boundary_indicators () const 
{
  return std::vector<unsigned char> (0);
}

#else

template <int dim>
std::vector<unsigned char>
Triangulation<dim>::get_boundary_indicators () const 
{
  std::vector<bool> bi_exists(255, false);
  active_cell_iterator cell=begin_active();
  for (; cell!=end(); ++cell)
    for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
      if (cell->at_boundary(face))
	bi_exists[cell->face(face)->boundary_indicator()]=true;

  const unsigned int n_bi=
    std::count(bi_exists.begin(), bi_exists.end(), true);

  std::vector<unsigned char> boundary_indicators(n_bi);
  unsigned int bi_counter=0;
  for (unsigned int i=0; i<bi_exists.size(); ++i)
    if (bi_exists[i]==true)
      boundary_indicators[bi_counter++]=i;

  return boundary_indicators;
}

#endif

/*--------- Put the next functions a bit out-or-order to avoid use before
  --------- explicit specialization, which is not allowed.                */

#if deal_II_dimension == 1

template <>
Triangulation<1>::cell_iterator
Triangulation<1>::begin (const unsigned int level) const
{
  return begin_line (level);
}


template <>
Triangulation<1>::raw_cell_iterator
Triangulation<1>::end () const
{
  return end_line ();
}

#endif


#if deal_II_dimension == 2

template <>
Triangulation<2>::cell_iterator
Triangulation<2>::begin (const unsigned int level) const
{
  return begin_quad (level);
}


template <>
Triangulation<2>::raw_cell_iterator
Triangulation<2>::end () const
{
  return end_quad ();
}

#endif


#if deal_II_dimension == 3

template <>
Triangulation<3>::cell_iterator
Triangulation<3>::begin (const unsigned int level) const
{
  return begin_hex (level);
}


template <>
Triangulation<3>::raw_cell_iterator
Triangulation<3>::end () const
{
  return end_hex ();
}

#endif


/*-----------------------------------------------------------------*/


template <int dim>
void Triangulation<dim>::copy_triangulation (const Triangulation<dim> &old_tria)
{
  Assert (vertices.size() == 0, ExcTriangulationNotEmpty());
  Assert (levels.size () == 0, ExcTriangulationNotEmpty());
  Assert (faces == NULL, ExcTriangulationNotEmpty());

  Assert (old_tria.levels.size() != 0, ExcInternalError());
  Assert (old_tria.vertices.size() != 0, ExcInternalError());
  Assert (dim == 1 || old_tria.faces != NULL, ExcInternalError());
  
  
				   // copy normal elements
  vertices      = old_tria.vertices;
  vertices_used = old_tria.vertices_used;
  smooth_grid   = old_tria.smooth_grid;
				   // do a shallow copy of the list of
				   // RefinementListeners
  refinement_listeners = old_tria.refinement_listeners;
  
  faces         = new internal::Triangulation::TriaFaces<dim>(*old_tria.faces);

  for (unsigned i=0;i<255;++i)
    boundary[i] = old_tria.boundary[i];

  levels.reserve (old_tria.levels.size());
  for (unsigned int level=0; level<old_tria.levels.size(); ++level)
    levels.push_back (new internal::Triangulation::TriaLevel<dim>(*old_tria.levels[level]));

  number_cache = old_tria.number_cache;
  
				   // note that we need not copy the
				   // subscriptor!
}


#if deal_II_dimension == 1

template <>
void Triangulation<1>::create_triangulation_compatibility (
  const std::vector<Point<1> >    &v,
  const std::vector<CellData<1> > &cells,
  const SubCellData               &subcelldata)
{
				   // no reordering required
  create_triangulation(v, cells, subcelldata);
}


template <>
void Triangulation<1>::create_triangulation (const std::vector<Point<1> >    &v,
					     const std::vector<CellData<1> > &cells,
					     const SubCellData &subcelldata)
{
				   // note: since no boundary
				   // information can be given in one
				   // dimension, the @p{subcelldata}
				   // field is ignored. (only used for
				   // error checking, which is a good
				   // idea in any case)
  
  const unsigned int dim=1;
  
  Assert (vertices.size() == 0, ExcTriangulationNotEmpty());
  Assert (levels.size() == 0, ExcTriangulationNotEmpty());
  Assert (faces == NULL, ExcTriangulationNotEmpty());
  
				   // check that no forbidden arrays
				   // are used
  Assert (subcelldata.check_consistency(dim), ExcInternalError());

				   // copy vertices
  vertices = v;
  vertices_used = std::vector<bool> (v.size(), true);
    
				   // store the indices of the lines
				   // which are adjacent to a given
				   // vertex
  std::vector<std::vector<int> > lines_at_vertex (v.size());

				   // reserve enough space
  levels.push_back (new internal::Triangulation::TriaLevel<dim>);
  levels[0]->internal::Triangulation::TriaLevel<0>
    ::reserve_space (cells.size(), dim);
  levels[0]->lines.reserve_space (cells.size());
  
				   // make up cells
  raw_line_iterator next_free_line = begin_raw_line ();
  for (unsigned int cell=0; cell<cells.size(); ++cell) 
    {
      while (next_free_line->used())
	++next_free_line;
      
      next_free_line->set (internal::Triangulation
                           ::Line (cells[cell].vertices[0],
                                   cells[cell].vertices[1]));
      next_free_line->set_used_flag ();
      next_free_line->set_material_id (cells[cell].material_id);
      next_free_line->clear_user_pointer ();
      next_free_line->set_subdomain_id (0);
      
				       // note that this cell is
				       // adjacent to these vertices
      lines_at_vertex[cells[cell].vertices[0]].push_back (cell);
      lines_at_vertex[cells[cell].vertices[1]].push_back (cell);
    }


				   // some security tests
  unsigned int boundary_nodes = 0;
  for (unsigned int i=0; i<lines_at_vertex.size(); ++i)
    switch (lines_at_vertex[i].size()) 
      {
	case 1:
					       // this vertex has only
					       // one adjacent line
	      ++boundary_nodes;
	      break;
	case 2:
	      break;
	default:
					       // a node must have one
					       // or two adjacent
					       // lines
	      clear_despite_subscriptions ();

	      AssertThrow (false, ExcInternalError());
      }

				   // assert there are no more than
				   // two boundary nodes
  if (boundary_nodes != 2)
    {
      clear_despite_subscriptions();
      
      AssertThrow (false, ExcInternalError());
    }
  


				   // update neighborship info
  active_line_iterator line = begin_active_line ();
				   // for all lines
  for (; line!=end(); ++line)
				     // for each of the two vertices
    for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
				       // if first cell adjacent to
				       // this vertex is the present
				       // one, then the neighbor is
				       // the second adjacent cell and
				       // vice versa
      if (lines_at_vertex[line->vertex_index(vertex)][0] == line->index())
	if (lines_at_vertex[line->vertex_index(vertex)].size() == 2) 
	  {
	    const cell_iterator neighbor (const_cast<Triangulation<1>*>(this),
					  0,              // level
					  lines_at_vertex[line->vertex_index(vertex)][1]);
	    line->set_neighbor (vertex, neighbor);
	  }
	else
					   // no second adjacent cell
					   // entered -> cell at
					   // boundary
	  line->set_neighbor (vertex, end());
      else
					 // present line is not first
					 // adjacent one -> first
					 // adjacent one is neighbor
	{
	  const cell_iterator neighbor (const_cast<Triangulation<1>*>(this),
					0,              // level
					lines_at_vertex[line->vertex_index(vertex)][0]);
	  line->set_neighbor (vertex, neighbor);
	}


				   // re-compute numbers of lines, etc
  update_number_cache ();  
}

#endif


#if deal_II_dimension == 2

template <>
void Triangulation<2>::create_triangulation_compatibility (
  const std::vector<Point<2> >    &v,
  const std::vector<CellData<2> > &c,
  const SubCellData               &subcelldata)
{
  const bool output=false;

				   // copy cell data as it will be
				   // modified in the following
  std::vector<CellData<2> > cells(c);
  
  for (unsigned int cell=0; cell<cells.size(); ++cell)
    {
      if (output)
	{
	  std::cout << "cell" << cell << ": ";
	  for (unsigned int vertex=0; vertex<4; ++vertex)
	    std::cout << cells[cell].vertices[vertex] << " ";
	}

      std::swap(cells[cell].vertices[2],cells[cell].vertices[3]);
      

      if (output)
	{
	  std::cout << "  after swap: ";
	  for (unsigned int vertex=0; vertex<4; ++vertex)
	    std::cout << cells[cell].vertices[vertex] << " ";
	  std::cout << std::endl;
	}
    }

				   // subcelldata.boundary_lines refer
				   // to the vertices in v. they don't
				   // need to be changed
  
  create_triangulation(v, cells, subcelldata);
}

template <>
void Triangulation<2>::create_triangulation (const std::vector<Point<2> >    &v,
					     const std::vector<CellData<2> > &cells,
					     const SubCellData               &subcelldata)
{
  const unsigned int dim=2;

  Assert (vertices.size() == 0, ExcTriangulationNotEmpty());
  Assert (levels.size() == 0, ExcTriangulationNotEmpty());
  Assert (faces == NULL, ExcTriangulationNotEmpty());
				   // check that no forbidden arrays
				   // are used
  Assert (subcelldata.check_consistency(dim), ExcInternalError());

				   // copy vertices
  vertices = v;
  vertices_used = std::vector<bool> (v.size(), true);

				   // make up a list of the needed
				   // lines each line is a pair of
				   // vertices. The list is kept
				   // sorted and it is guaranteed that
				   // each line is inserted only once.
				   // While the key of such an entry
				   // is the pair of vertices, the
				   // thing it points to is an
				   // iterator pointing to the line
				   // object itself. In the first run,
				   // these iterators are all invalid
				   // ones, but they are filled
				   // afterwards
  std::map<std::pair<int,int>,line_iterator> needed_lines;
  for (unsigned int cell=0; cell<cells.size(); ++cell)
    {
      for (unsigned int vertex=0; vertex<4; ++vertex)
	if ( ! (cells[cell].vertices[vertex]<vertices.size()))
	  {
                                             // store the number of
                                             // vertices, as this
                                             // information will be
                                             // deleted by clear()
                                             // though we'd like to
                                             // include it in the
                                             // error message
            const unsigned int n_vertices = vertices.size();

	    clear_despite_subscriptions();
	    
	    Assert (false,
		    ExcInvalidVertexIndex (cell, cells[cell].vertices[vertex],
                                           n_vertices));
	  }


      for (unsigned int line=0; line<GeometryInfo<dim>::faces_per_cell; ++line)
	{
				       // given a line vertex number
				       // (0,1) on a specific line we
				       // get the cell vertex number
				       // (0-4) through the
				       // line_to_cell_vertices
				       // function
	  std::pair<int,int> line_vertices(
	    cells[cell].vertices[GeometryInfo<dim>::line_to_cell_vertices(line, 0)],
	    cells[cell].vertices[GeometryInfo<dim>::line_to_cell_vertices(line, 1)]);

					   // assert that the line was
					   // not already inserted in
					   // reverse order. This
					   // happens in spite of the
					   // vertex rotation above,
					   // if the sense of the cell
					   // was incorrect.
					   //
					   // Here is what usually
					   // happened when this
					   // exception is thrown:
					   // consider these two cells
					   // and the vertices
					   //  3---4---5
					   //  |   |   |
					   //  0---1---2
					   // If in the input vector
					   // the two cells are given
					   // with vertices <0 1 4 3>
					   // and <4 1 2 5>, in the
					   // first cell the middle
					   // line would have
					   // direction 1->4, while in
					   // the second it would be
					   // 4->1.  This will cause
					   // the exception.
	  if (! (needed_lines.find(std::make_pair(line_vertices.second,
						  line_vertices.first))
		 ==
		 needed_lines.end()))
	    {
	      clear_despite_subscriptions();
	      
	      AssertThrow (false,
			   ExcGridHasInvalidCell(cell));
	    }
		  
					   // insert line, with
					   // invalid iterator if line
					   // already exists, then
					   // nothing bad happens here
	  needed_lines[line_vertices] = end_line();
	}
    }


				   // check that every vertex has at
				   // least two adjacent lines
  if (true) 
    {
      std::vector<unsigned short int> vertex_touch_count (v.size(), 0);
      std::map<std::pair<int,int>,line_iterator>::iterator i;
      for (i=needed_lines.begin(); i!=needed_lines.end(); i++) 
	{
					   // touch the vertices of
					   // this line
	  ++vertex_touch_count[i->first.first];
	  ++vertex_touch_count[i->first.second];
	}

				       // assert minimum touch count
				       // is at least two. if not so,
				       // then clean triangulation and
				       // exit with an exception
      if ( ! (* (std::min_element(vertex_touch_count.begin(),
				  vertex_touch_count.end())) >= 2))
	{
	  clear_despite_subscriptions();
	  
	  AssertThrow (false, ExcGridHasInvalidVertices());
	}
    }
	
  				   // reserve enough space
  levels.push_back (new internal::Triangulation::TriaLevel<dim>);
  faces = new internal::Triangulation::TriaFaces<dim>;
  levels[0]->internal::Triangulation::TriaLevel<0>
    ::reserve_space (cells.size(), dim);
  faces->lines.reserve_space (needed_lines.size());
  levels[0]->quads.reserve_space (cells.size());

				   // make up lines
  if (true) 
    {
      raw_line_iterator line = begin_raw_line();
      std::map<std::pair<int,int>,line_iterator>::iterator i;
      for (i = needed_lines.begin(); line!=end_line(); ++line, ++i) 
	{
	  line->set (internal::Triangulation::Line(i->first.first,
                                                   i->first.second));
	  line->set_used_flag ();
	  line->clear_user_flag ();
	  line->clear_user_pointer ();
	  i->second = line;
	}
    }


				   // store for each line index
				   // the adjacent cells
  std::map<int,std::vector<cell_iterator> > adjacent_cells;

				   // finally make up cells
  if (true) 
    {
      raw_cell_iterator cell = begin_raw_quad();
      for (unsigned int c=0; c<cells.size(); ++c, ++cell)
	{
	  line_iterator lines[GeometryInfo<dim>::lines_per_cell];
	  for (unsigned int line=0; line<GeometryInfo<dim>::lines_per_cell; ++line)
	    lines[line]=needed_lines[std::make_pair(
	      cells[c].vertices[GeometryInfo<dim>::line_to_cell_vertices(line, 0)],
	      cells[c].vertices[GeometryInfo<dim>::line_to_cell_vertices(line, 1)])];
	  
	  cell->set (internal::Triangulation::Quad(lines[0]->index(),
                                                   lines[1]->index(),
                                                   lines[2]->index(),
                                                   lines[3]->index()));
	  
	  cell->set_used_flag ();
	  cell->set_material_id (cells[c].material_id);
	  cell->clear_user_pointer ();
	  cell->set_subdomain_id (0);
	  
					   // note that this cell is
					   // adjacent to the four
					   // lines
	  for (unsigned int line=0; line<GeometryInfo<dim>::lines_per_cell; ++line)
	    adjacent_cells[lines[line]->index()].push_back (cell);
	}
    }


  for (line_iterator line=begin_line(); line!=end_line(); ++line) 
    {
      const unsigned int n_adj_cells = adjacent_cells[line->index()].size();
				       // assert that every line has
				       // one or two adjacent cells
      if (! ((n_adj_cells >= 1) &&
	     (n_adj_cells <= 2)))
	{
	  clear_despite_subscriptions();
	  
	  AssertThrow (false, ExcInternalError());
	}

				       // if only one cell: line is at
				       // boundary -> give it the
				       // boundary indicator zero by
				       // default
      if (n_adj_cells == 1)
	line->set_boundary_indicator (0);
      else
					 // interior line -> 255
      	line->set_boundary_indicator (255);
    }

				   // set boundary indicators where
				   // given
  std::vector<CellData<1> >::const_iterator boundary_line
    = subcelldata.boundary_lines.begin();
  std::vector<CellData<1> >::const_iterator end_boundary_line
    = subcelldata.boundary_lines.end();
  for (; boundary_line!=end_boundary_line; ++boundary_line) 
    {
      line_iterator line;
      std::pair<int,int> line_vertices(std::make_pair(boundary_line->vertices[0],
						      boundary_line->vertices[1]));
      if (needed_lines.find(line_vertices) != needed_lines.end())
					 // line found in this
					 // direction
	line = needed_lines[line_vertices];
      else 
	{
					   // look whether it exists
					   // in reverse direction
	  std::swap (line_vertices.first, line_vertices.second);
	  if (needed_lines.find(line_vertices) != needed_lines.end())
	    line = needed_lines[line_vertices];
	  else 
	    {
					       // line does not exist
	      clear_despite_subscriptions();

	      AssertThrow (false, ExcLineInexistant(line_vertices.first,
						    line_vertices.second));
	    }
	}

				       // assert that we only set
				       // boundary info once
      if (line->boundary_indicator() != 0 &&
	  line->boundary_indicator() != 255)
	{
	  clear_despite_subscriptions();
	  
	  AssertThrow (false, ExcMultiplySetLineInfoOfLine(
	    line_vertices.first, line_vertices.second));
	}
      
				       // Assert that only exterior lines
				       // are given a boundary indicator
      if (line->boundary_indicator() == 255)
	{
	  clear_despite_subscriptions();
	  
	  AssertThrow (false, ExcInteriorLineCantBeBoundary());
	}

      line->set_boundary_indicator (boundary_line->material_id);
    }


				   // finally update neighborship info
  for (cell_iterator cell=begin(); cell!=end(); ++cell)
    for (unsigned int side=0; side<4; ++side)
      if (adjacent_cells[cell->line(side)->index()][0] == cell)
					 // first adjacent cell is
					 // this one
	{
	  if (adjacent_cells[cell->line(side)->index()].size() == 2)
					     // there is another
					     // adjacent cell
	    cell->set_neighbor (side,
				adjacent_cells[cell->line(side)->index()][1]);
	}
				   // first adjacent cell is not this
				   // one, -> it must be the neighbor
				   // we are looking for
      else
	cell->set_neighbor (side,
			    adjacent_cells[cell->line(side)->index()][0]);


				   // re-compute numbers of lines, etc
  update_number_cache ();  
}

#endif


#if deal_II_dimension == 3

/**
 * Invent an object which compares two internal::Triangulation::Quad against
 * each other. This comparison is needed in order to establish a map of Quads
 * to iterators in the Triangulation<3>::create_triangulation function.
 *
 * Since this comparison is not canonical, we do not include it into
 * the general internal::Triangulation::Quad class.
 */
struct QuadComparator
{
    inline bool operator () (const internal::Triangulation::Quad &q1,
                             const internal::Triangulation::Quad &q2) const
      {
					 // here is room to optimize
					 // the repeated equality test
					 // of the previous lines, but
					 // I don't care at present
	if ((q1.line(0) < q2.line(0))          ||
	    ((q1.line(0) == q2.line(0)) &&
	     (q1.line(1) <  q2.line(1)))       ||
	    ((q1.line(0) == q2.line(0)) &&
	     (q1.line(1) == q2.line(1)) &&
	     (q1.line(2) <  q2.line(2)))       ||
	    ((q1.line(0) == q2.line(0)) &&
	     (q1.line(1) == q2.line(1)) &&
	     (q1.line(2) == q2.line(2)) &&
	     (q1.line(3) <  q2.line(3))))
	  return true;
	else
	  return false;
      }
};


template <>
void
Triangulation<3>::create_triangulation_compatibility (
  const std::vector<Point<3> >    &v,
  const std::vector<CellData<3> > &c,
  const SubCellData               &subc)
{
  const bool output1=false;
  const bool output2=false;

				   // copy cell and subcell data as
				   // they will be modified in the
				   // following
  std::vector<CellData<3> > cells(c);
  SubCellData subcelldata(subc);
  
  unsigned int tmp[GeometryInfo<3>::vertices_per_cell];
  for (unsigned int cell=0; cell<cells.size(); ++cell)
    {
      if (output1)
	{
	  std::cout << "cell" << cell << ": ";
	  for (unsigned int vertex=0;
	       vertex<GeometryInfo<3>::vertices_per_cell; ++vertex)
	    std::cout << cells[cell].vertices[vertex] << " ";
	}

      for (unsigned int i=0; i<GeometryInfo<3>::vertices_per_cell; ++i)
	tmp[i]=cells[cell].vertices[i];
      for (unsigned int i=0; i<GeometryInfo<3>::vertices_per_cell; ++i)
	cells[cell].vertices[GeometryInfo<3>::ucd_to_deal[i]]=tmp[i];

      if (output1)
	{
	  std::cout << "  after swap: ";
	  for (unsigned int vertex=0;
	       vertex<GeometryInfo<3>::vertices_per_cell; ++vertex)
	    std::cout << cells[cell].vertices[vertex] << " ";
	  std::cout << std::endl;
	  std::cout << "vertices:" << std::endl;
	  for (unsigned int vertex=0;
	       vertex<GeometryInfo<3>::vertices_per_cell; ++vertex)
	    std::cout << v[cells[cell].vertices[vertex]] << std::endl;
	}
    }

  if (output2)
    {
      std::cout << "vertices:" << std::endl;
      for (unsigned int i=0; i<vertices.size(); ++i)
	std::cout << vertices[i] << std::endl;
    }
				   // now points in boundary quads
  std::vector<CellData<2> >::iterator boundary_quad
    = subcelldata.boundary_quads.begin();
  std::vector<CellData<2> >::iterator end_quad
    = subcelldata.boundary_quads.end();
  for (unsigned int quad_no=0; boundary_quad!=end_quad; ++boundary_quad, ++quad_no)
    {
      if (output2)
	std::cout << quad_no << ": " << boundary_quad->vertices[0]
		  << ",  " << boundary_quad->vertices[1]
		  << ",  " << boundary_quad->vertices[2]
		  << ",  " << boundary_quad->vertices[3];
      
      std::swap(boundary_quad->vertices[2], boundary_quad->vertices[3]);

      if (output2)
	std::cout << ", after swap: " << boundary_quad->vertices[0]
		  << ",  " << boundary_quad->vertices[1]
		  << ",  " << boundary_quad->vertices[2]
		  << ",  " << boundary_quad->vertices[3] << std::endl;
    }  

  create_triangulation(v, cells, subcelldata);
}


template <>
void
Triangulation<3>::create_triangulation (const std::vector<Point<3> >    &v,
                                        const std::vector<CellData<3> > &cells,
                                        const SubCellData               &subcelldata)
{
  const unsigned int dim=3;

  Assert (vertices.size() == 0, ExcTriangulationNotEmpty());
  Assert (levels.size() == 0, ExcTriangulationNotEmpty());
  Assert (faces == NULL, ExcTriangulationNotEmpty());
				   // check that no forbidden arrays
				   // are used
  Assert (subcelldata.check_consistency(dim), ExcInternalError());

				   // copy vertices
  vertices = v;
  vertices_used = std::vector<bool> (v.size(), true);

				   // check that all cells have
				   // positive volume. if not call the
				   // invert_all_cells_of_negative_grid
				   // and reorder_cells function of
				   // GridReordering before creating
				   // the triangulation
  for (unsigned int cell_no=0; cell_no<cells.size(); ++cell_no)
    if (GridTools::cell_measure(vertices, cells[cell_no].vertices) < 0)
      {
	clear_despite_subscriptions();
	
	AssertThrow (false, ExcGridHasInvalidCell(cell_no));
      }

				   ///////////////////////////////////////
				   // first set up some collections of data
				   //
				   // make up a list of the needed
				   // lines
				   //
				   // each line is a pair of
				   // vertices. The list is kept
				   // sorted and it is guaranteed that
				   // each line is inserted only once.
				   // While the key of such an entry
				   // is the pair of vertices, the
				   // thing it points to is an
				   // iterator pointing to the line
				   // object itself. In the first run,
				   // these iterators are all invalid
				   // ones, but they are filled
				   // afterwards same applies for the
				   // quads
  std::map<std::pair<int,int>,line_iterator> needed_lines;
  for (unsigned int cell=0; cell<cells.size(); ++cell)
    {

				       // check whether vertex indices
				       // are valid ones
      for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
	if (! (cells[cell].vertices[vertex] < vertices.size()))
	  {
                                             // store the number of
                                             // vertices, as this
                                             // information will be
                                             // deleted by clear()
                                             // though we'd like to
                                             // include it in the
                                             // error message
            const unsigned int n_vertices = vertices.size();

	    clear_despite_subscriptions();            
	    
	    AssertThrow (false,
			 ExcInvalidVertexIndex (cell, cells[cell].vertices[vertex],
						n_vertices));
	  }      

      
      for (unsigned int line=0; line<GeometryInfo<dim>::lines_per_cell; ++line)
	{
				       // given a line vertex number
				       // (0,1) on a specific line we
				       // get the cell vertex number
				       // (0-7) through the
				       // line_to_cell_vertices
				       // function
	  std::pair<int,int> line_vertices(
	    cells[cell].vertices[GeometryInfo<dim>::line_to_cell_vertices(line, 0)],
	    cells[cell].vertices[GeometryInfo<dim>::line_to_cell_vertices(line, 1)]);
	  
					   // if that line was already inserted
					   // in reverse order do nothing, else
					   // insert the line
	  if ( (needed_lines.find(std::make_pair(line_vertices.second,
						  line_vertices.first))
		 ==
		 needed_lines.end()))
	    {
					       // insert line, with
					       // invalid iterator. if line
					       // already exists, then
					       // nothing bad happens here
	      needed_lines[line_vertices] = end_line();
	    }
	}
    }


				   /////////////////////////////////
				   // now for some sanity-checks:
				   //
				   // check that every vertex has at
				   // least tree adjacent lines
  if (true) 
    {
      std::vector<unsigned short int> vertex_touch_count (v.size(), 0);
      std::map<std::pair<int,int>,line_iterator>::iterator i;
      for (i=needed_lines.begin(); i!=needed_lines.end(); i++) 
	{
					   // touch the vertices of
					   // this line
	  ++vertex_touch_count[i->first.first];
	  ++vertex_touch_count[i->first.second];
	}

 				       // assert minimum touch count
 				       // is at least three. if not so,
				       // then clean triangulation and
				       // exit with an exception
      if (! (* (std::min_element(vertex_touch_count.begin(),
				 vertex_touch_count.end())) >= 3))
	{
	  clear_despite_subscriptions();

	  AssertThrow (false, ExcGridHasInvalidVertices());
	}
    }


				   ///////////////////////////////////
				   // actually set up data structures
				   // for the lines
  				   // reserve enough space
  levels.push_back (new internal::Triangulation::TriaLevel<dim>);
  faces = new internal::Triangulation::TriaFaces<dim>;
  levels[0]->internal::Triangulation::TriaLevel<0>
    ::reserve_space (cells.size(), dim);
  faces->lines.reserve_space (needed_lines.size());

				   // make up lines
  if (true) 
    {
      raw_line_iterator line = begin_raw_line();
      std::map<std::pair<int,int>,line_iterator>::iterator i;
      for (i = needed_lines.begin(); line!=end_line(); ++line, ++i) 
	{
	  line->set (internal::Triangulation::Line(i->first.first,
                                                   i->first.second));
	  line->set_used_flag ();
	  line->clear_user_flag ();
	  line->clear_user_pointer ();

					   // now set the iterator for
					   // this line
	  i->second = line;
	}
    }


				   ///////////////////////////////////////////
				   // make up the quads of this triangulation
				   //
				   // same thing: the iterators are
				   // set to the invalid value at
				   // first, we only collect the data
				   // now

				   // the bool array stores, wether the lines
				   // are in the standard orientation or not

				   // note that QuadComparator is a
				   // class declared and defined in
				   // this file
  std::map<internal::Triangulation::Quad,std::pair<quad_iterator,boost::array<bool,GeometryInfo<dim>::lines_per_face> >,QuadComparator>
    needed_quads;
  for (unsigned int cell=0; cell<cells.size(); ++cell) 
    {
				       // the faces are quads which
				       // consist of four numbers
				       // denoting the index of the
				       // four lines bounding the
				       // quad. we can get this index
				       // by asking @p{needed_lines}
				       // for an iterator to this
				       // line, dereferencing it and
				       // thus return an iterator into
				       // the @p{lines} array of the
				       // triangulation, which is
				       // already set up. we can then
				       // ask this iterator for its
				       // index within the present
				       // level (the level is zero, of
				       // course)
				       //
				       // to make things easier, we
				       // don't create the lines
				       // (pairs of their vertex
				       // indices) in place, but
				       // before they are really
				       // needed.
      std::pair<int,int> line_list[GeometryInfo<dim>::lines_per_cell],
	inverse_line_list[GeometryInfo<dim>::lines_per_cell];
      unsigned int face_line_list[GeometryInfo<dim>::lines_per_face];
      boost::array<bool,GeometryInfo<dim>::lines_per_face> orientation;
      
      for (unsigned int line=0; line<GeometryInfo<dim>::lines_per_cell; ++line)
	{
	  line_list[line]=std::pair<int,int> (
	    cells[cell].vertices[GeometryInfo<dim>::line_to_cell_vertices(line, 0)],
	    cells[cell].vertices[GeometryInfo<dim>::line_to_cell_vertices(line, 1)]);
	  inverse_line_list[line]=std::pair<int,int> (
	    cells[cell].vertices[GeometryInfo<dim>::line_to_cell_vertices(line, 1)],
	    cells[cell].vertices[GeometryInfo<dim>::line_to_cell_vertices(line, 0)]);
	}
      
      for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
        {
					   // set up a list of the lines to be
					   // used for this face. check the
					   // direction for each line
					   //
					   // given a face line number (0-3) on
					   // a specific face we get the cell
					   // line number (0-11) through the
					   // face_to_cell_lines function
	  for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_face; ++l)
	    if (needed_lines.find (inverse_line_list[GeometryInfo<dim>::
						     face_to_cell_lines(face,l)]) == needed_lines.end())
	      {
		face_line_list[l]=needed_lines[line_list[GeometryInfo<dim>::
							 face_to_cell_lines(face,l)]]->index();
		orientation[l]=true;
	      }
	    else
	      {
		face_line_list[l]=needed_lines[inverse_line_list[GeometryInfo<dim>::
								 face_to_cell_lines(face,l)]]->index();
		orientation[l]=false;
	      }
	    
	  
	  internal::Triangulation::Quad
            quad(face_line_list[0],
		 face_line_list[1],
		 face_line_list[2],
		 face_line_list[3]);
	  
                                           // insert quad, with
                                           // invalid iterator
                                           //
                                           // if quad already exists,
                                           // then nothing bad happens
                                           // here, as this will then
                                           // simply become an
                                           // interior face of the
                                           // triangulation. however,
                                           // we will run into major
                                           // trouble if the face was
                                           // already inserted in the
                                           // opposite
                                           // direction. there are
                                           // really only two
                                           // orientations for a face
                                           // to be in, since the edge
                                           // directions are already
                                           // set. thus, vertex 0 is
                                           // the one from which two
                                           // edges originate, and
                                           // vertex 3 is the one to
                                           // which they converge. we
                                           // are then left with
                                           // orientations 0-1-2-3 and
                                           // 2-3-0-1 for the order of
                                           // lines. the
                                           // corresponding quad can
                                           // be easily constructed by
                                           // exchanging lines. we do
                                           // so here, just to check
                                           // that that flipped quad
                                           // isn't already in the
                                           // triangulation. if it is,
                                           // then don't insert the
                                           // new one and instead
                                           // later set the
                                           // face_orientation flag
          const internal::Triangulation::Quad
	    test_quad_1(quad.line(2), quad.line(3),
			quad.line(0), quad.line(1)),//face_orientation=false, face_flip=false, face_rotation=false
	    test_quad_2(quad.line(0), quad.line(1),
			quad.line(3), quad.line(2)),//face_orientation=false, face_flip=false, face_rotation=true
	    test_quad_3(quad.line(3), quad.line(2),
			quad.line(1), quad.line(0)),//face_orientation=false, face_flip=true,  face_rotation=false
	    test_quad_4(quad.line(1), quad.line(0),
			quad.line(2), quad.line(3)),//face_orientation=false, face_flip=true,  face_rotation=true
	    test_quad_5(quad.line(2), quad.line(3),
			quad.line(1), quad.line(0)),//face_orientation=true,  face_flip=false, face_rotation=true
	    test_quad_6(quad.line(1), quad.line(0),
			quad.line(3), quad.line(2)),//face_orientation=true,  face_flip=true,  face_rotation=false
	    test_quad_7(quad.line(3), quad.line(2),
			quad.line(0), quad.line(1));//face_orientation=true,  face_flip=true,  face_rotation=true
          if (needed_quads.find (test_quad_1) == needed_quads.end() &&
	      needed_quads.find (test_quad_2) == needed_quads.end() &&
	      needed_quads.find (test_quad_3) == needed_quads.end() &&
	      needed_quads.find (test_quad_4) == needed_quads.end() &&
	      needed_quads.find (test_quad_5) == needed_quads.end() &&
	      needed_quads.find (test_quad_6) == needed_quads.end() &&
	      needed_quads.find (test_quad_7) == needed_quads.end())
            needed_quads[quad] = std::make_pair(end_quad(),orientation);
        }
    }


				   /////////////////////////////////
				   // enter the resulting quads into
				   // the arrays of the Triangulation
				   //
				   // first reserve enough space
  faces->quads.reserve_space (needed_quads.size());
  
  if (true) 
    {
      raw_quad_iterator quad = begin_raw_quad();
      std::map<internal::Triangulation::Quad,std::pair<quad_iterator,boost::array<bool,GeometryInfo<dim>::lines_per_face> >,QuadComparator>
        ::iterator q;
      for (q = needed_quads.begin(); quad!=end_quad(); ++quad, ++q)
	{
	  quad->set (q->first);
	  quad->set_used_flag ();
	  quad->clear_user_flag ();
	  quad->clear_user_pointer ();
					   // set the line orientation
	  quad->set_line_orientation(0,q->second.second[0]);
	  quad->set_line_orientation(1,q->second.second[1]);
	  quad->set_line_orientation(2,q->second.second[2]);
	  quad->set_line_orientation(3,q->second.second[3]);
	  

					   // now set the iterator for
					   // this quad
	  q->second.first = quad;
	}
    }

				   /////////////////////////////////
				   // finally create the cells
  levels[0]->hexes.reserve_space (cells.size());  

				   // store for each quad index the
				   // adjacent cells
  std::map<int,std::vector<cell_iterator> > adjacent_cells;

				   // finally make up cells
  if (true) 
    {
      raw_cell_iterator cell = begin_raw_hex();
      for (unsigned int c=0; c<cells.size(); ++c, ++cell)
	{
					   // first find for each of
					   // the cells the quad
					   // iterator of the
					   // respective faces.
					   //
					   // to this end, set up the
					   // lines of this cell and
					   // find the quads that are
					   // bounded by these lines;
					   // these are then the faces
					   // of the present cell
	  std::pair<int,int> line_list[GeometryInfo<dim>::lines_per_cell],
	    inverse_line_list[GeometryInfo<dim>::lines_per_cell];
	  unsigned int face_line_list[4];
	  for (unsigned int line=0; line<GeometryInfo<dim>::lines_per_cell; ++line)
	    {
	      line_list[line]=std::make_pair(
		cells[c].vertices[GeometryInfo<dim>::line_to_cell_vertices(line, 0)],
		cells[c].vertices[GeometryInfo<dim>::line_to_cell_vertices(line, 1)]);
	      inverse_line_list[line]=std::pair<int,int> (
		cells[c].vertices[GeometryInfo<dim>::line_to_cell_vertices(line, 1)],
		cells[c].vertices[GeometryInfo<dim>::line_to_cell_vertices(line, 0)]);
	    }
	  
					   // get the iterators
					   // corresponding to the
					   // faces. also store
					   // whether they are
					   // reversed or not
	  quad_iterator face_iterator[GeometryInfo<dim>::faces_per_cell];
          bool face_orientation[GeometryInfo<dim>::faces_per_cell];
          bool face_flip[GeometryInfo<dim>::faces_per_cell];
          bool face_rotation[GeometryInfo<dim>::faces_per_cell];
          for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	    {
	      for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_face; ++l)
		if (needed_lines.find (inverse_line_list[GeometryInfo<dim>::
							 face_to_cell_lines(face,l)]) == needed_lines.end())
		  face_line_list[l]=needed_lines[line_list[GeometryInfo<dim>::
							   face_to_cell_lines(face,l)]]->index();
		else
		  face_line_list[l]=needed_lines[inverse_line_list[GeometryInfo<dim>::
								   face_to_cell_lines(face,l)]]->index();
	      
	      internal::Triangulation::Quad
		quad(face_line_list[0],
		     face_line_list[1],
		     face_line_list[2],
		     face_line_list[3]);
	      
	      if (needed_quads.find (quad) != needed_quads.end())
		{
						   // face is in standard
						   // orientation (and not
						   // flipped or rotated). this
						   // must be true for at least
						   // one of the two cells
						   // containing this face
						   // (i.e. for the cell which
						   // originally inserted the
						   // face)
		  face_iterator[face] = needed_quads[quad].first;
		  face_orientation[face] = true;
		  face_flip[face]=false;
		  face_rotation[face]=false;
		}
	      else
		{
						   // face must be available in
						   // reverse order
						   // then. construct all
						   // possibilities and check
						   // them one after the other
		  const internal::Triangulation::Quad
		    test_quad_1(quad.line(2), quad.line(3),
				quad.line(0), quad.line(1)),//face_orientation=false, face_flip=false, face_rotation=false
		    test_quad_2(quad.line(0), quad.line(1),
				quad.line(3), quad.line(2)),//face_orientation=false, face_flip=false, face_rotation=true
		    test_quad_3(quad.line(3), quad.line(2),
				quad.line(1), quad.line(0)),//face_orientation=false, face_flip=true,  face_rotation=false
		    test_quad_4(quad.line(1), quad.line(0),
				quad.line(2), quad.line(3)),//face_orientation=false, face_flip=true,  face_rotation=true
		    test_quad_5(quad.line(2), quad.line(3),
				quad.line(1), quad.line(0)),//face_orientation=true,  face_flip=false, face_rotation=true
		    test_quad_6(quad.line(1), quad.line(0),
				quad.line(3), quad.line(2)),//face_orientation=true,  face_flip=true,  face_rotation=false
		    test_quad_7(quad.line(3), quad.line(2),
				quad.line(0), quad.line(1));//face_orientation=true,  face_flip=true,  face_rotation=true
		  if (needed_quads.find (test_quad_1) != needed_quads.end())
		    {
		      face_iterator[face] = needed_quads[test_quad_1].first;
		      face_orientation[face] = false;
		      face_flip[face]=false;
		      face_rotation[face]=false;
		    }
		  else if (needed_quads.find (test_quad_2) != needed_quads.end())
		    {
		      face_iterator[face] = needed_quads[test_quad_2].first;
		      face_orientation[face] = false;
		      face_flip[face]=false;
		      face_rotation[face]=true;
		    }
		  else if (needed_quads.find (test_quad_3) != needed_quads.end())
		    {
		      face_iterator[face] = needed_quads[test_quad_3].first;
		      face_orientation[face] = false;
		      face_flip[face]=true;
		      face_rotation[face]=false;
		    }
		  else if (needed_quads.find (test_quad_4) != needed_quads.end())
		    {
		      face_iterator[face] = needed_quads[test_quad_4].first;
		      face_orientation[face] = false;
		      face_flip[face]=true;
		      face_rotation[face]=true;
		    }
		  else if (needed_quads.find (test_quad_5) != needed_quads.end())
		    {
		      face_iterator[face] = needed_quads[test_quad_5].first;
		      face_orientation[face] = true;
		      face_flip[face]=false;
		      face_rotation[face]=true;
		    }
		  else if (needed_quads.find (test_quad_6) != needed_quads.end())
		    {
		      face_iterator[face] = needed_quads[test_quad_6].first;
		      face_orientation[face] = true;
		      face_flip[face]=true;
		      face_rotation[face]=false;
		    }
		  else if (needed_quads.find (test_quad_7) != needed_quads.end())
		    {
		      face_iterator[face] = needed_quads[test_quad_7].first;
		      face_orientation[face] = true;
		      face_flip[face]=true;
		      face_rotation[face]=true;
		    }

		  else
						       // we didn't find the
						       // face in any direction,
						       // so something went
						       // wrong above
		    Assert(false,ExcInternalError());
		  
		}
	    }// for all faces
	  
					   // make the cell out of
					   // these iterators
	  cell->set (internal::Triangulation
                     ::Hexahedron(face_iterator[0]->index(),
                                  face_iterator[1]->index(),
                                  face_iterator[2]->index(),
                                  face_iterator[3]->index(),
                                  face_iterator[4]->index(),
                                  face_iterator[5]->index()));
	  
	  cell->set_used_flag ();
	  cell->set_material_id (cells[c].material_id);
	  cell->clear_user_flag ();
	  cell->clear_user_pointer ();
	  cell->set_subdomain_id (0);

                                           // set orientation flag for
                                           // each of the faces
	  for (unsigned int quad=0; quad<GeometryInfo<dim>::faces_per_cell; ++quad)
	    {
	      cell->set_face_orientation (quad, face_orientation[quad]);
	      cell->set_face_flip (quad, face_flip[quad]);
	      cell->set_face_rotation (quad, face_rotation[quad]);
	    }
	  
          
					   // note that this cell is
					   // adjacent to the six
					   // quads
	  for (unsigned int quad=0; quad<GeometryInfo<dim>::faces_per_cell; ++quad)
	    adjacent_cells[face_iterator[quad]->index()].push_back (cell);

#ifdef DEBUG          
					   // make some checks on the
					   // lines and their
					   // ordering

					   // first map all cell lines
					   // to the two face lines
					   // which should
					   // coincide. all face lines
					   // are included with a cell
					   // line number (0-11)
					   // key. At the end all keys
					   // will be included twice
					   // (for each of the two
					   // coinciding lines once)
	  std::multimap<unsigned int, std::pair<unsigned int, unsigned int> >
	    cell_to_face_lines;
	  for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	    for (unsigned int line=0; line<GeometryInfo<dim>::lines_per_face; ++line)
	      cell_to_face_lines.insert(
		std::pair<unsigned int, std::pair<unsigned int, unsigned int> > (
		  GeometryInfo<dim>::face_to_cell_lines(face,line),
		  std::pair<unsigned int, unsigned int> (face,line)));
	  std::multimap<unsigned int, std::pair<unsigned int, unsigned int> >::const_iterator
	    map_iter=cell_to_face_lines.begin();
	  
	  for (; map_iter!=cell_to_face_lines.end(); ++map_iter)
	    {
	      const unsigned int cell_line=map_iter->first;
	      const unsigned int face1=map_iter->second.first;
	      const unsigned int line1=map_iter->second.second;
	      ++map_iter;
	      Assert(map_iter!=cell_to_face_lines.end(), ExcInternalErrorOnCell(c));
	      Assert(map_iter->first==cell_line, ExcInternalErrorOnCell(c));
	      const unsigned int face2=map_iter->second.first;
	      const unsigned int line2=map_iter->second.second;

					       // check that the pair
					       // of lines really
					       // coincide. Take care
					       // about the face
					       // orientation;
	      Assert (face_iterator[face1]->line(GeometryInfo<dim>::standard_to_real_face_line(
						   line1,
						   face_orientation[face1],
						   face_flip[face1],
						   face_rotation[face1])) ==
		      face_iterator[face2]->line(GeometryInfo<dim>::standard_to_real_face_line(
						   line2,
						   face_orientation[face2],
						   face_flip[face2],
						   face_rotation[face2])),
		      ExcInternalErrorOnCell(c));
	    }
#endif
	}
    }


				   /////////////////////////////////////////
				   // find those quads which are at the
				   // boundary and mark them appropriately
  for (quad_iterator quad=begin_quad(); quad!=end_quad(); ++quad)
    {
      const unsigned int n_adj_cells = adjacent_cells[quad->index()].size();
				       // assert that every quad has
				       // one or two adjacent cells
      if (! ((n_adj_cells >= 1) &&
	     (n_adj_cells <= 2)))
	{
	  clear_despite_subscriptions();
	  
	  AssertThrow (false, ExcInternalError());
	}

				       // if only one cell: quad is at
				       // boundary -> give it the
				       // boundary indicator zero by
				       // default
      if (n_adj_cells == 1)
	quad->set_boundary_indicator (0);
      else
					 // interior quad -> 255
      	quad->set_boundary_indicator (255);
    }

				   /////////////////////////////////////////
				   // next find those lines which are at
				   // the boundary and mark all others as
				   // interior ones
				   //
				   // for this: first mark all lines
				   // as interior
  for (line_iterator line=begin_line(); line!=end_line(); ++line)
    line->set_boundary_indicator (255);
				   // next reset all lines bounding
				   // boundary quads as on the
				   // boundary also. note that since
				   // we are in 3d, there are cases
				   // where one or more lines of a
				   // quad that is not on the
				   // boundary, are actually boundary
				   // lines. they will not be marked
				   // when visiting this
				   // face. however, since we do not
				   // support dim-2 dimensional
				   // boundaries (i.e. internal lines
				   // constituting boundaries), every
				   // such line is also part of a face
				   // that is actually on the
				   // boundary, so sooner or later we
				   // get to mark that line for being
				   // on the boundary
  for (quad_iterator quad=begin_quad(); quad!=end_quad(); ++quad)
    if (quad->at_boundary())
      for (unsigned int l=0; l<4; ++l)
	quad->line(l)->set_boundary_indicator (0);
  
				   // Check that all cells have
				   // positive volume. As the input
				   // CellData already passed this
				   // test we don't really expect
				   // negative measures here
  for (active_cell_iterator cell=begin_active(); cell!=end(); ++cell)
    Assert(cell->measure()>0, ExcInternalError());
      
				   ///////////////////////////////////////
				   // now set boundary indicators
				   // where given
				   //
				   // first do so for lines
  std::vector<CellData<1> >::const_iterator boundary_line
    = subcelldata.boundary_lines.begin();
  std::vector<CellData<1> >::const_iterator end_boundary_line
    = subcelldata.boundary_lines.end();
  for (; boundary_line!=end_boundary_line; ++boundary_line)
    {
      line_iterator line;
      std::pair <int, int> line_vertices(std::make_pair(boundary_line->vertices[0],
							boundary_line->vertices[1]));
      if (needed_lines.find(line_vertices) != needed_lines.end())
					 // line found in this
					 // direction
	line = needed_lines[line_vertices];
 	  	  
      else
	{
					   // look wether it exists in
					   // reverse direction
	  std::swap (line_vertices.first, line_vertices.second);
	  if (needed_lines.find(line_vertices) != needed_lines.end())
	    line = needed_lines[line_vertices];
	  else
	    {
	      clear_despite_subscriptions ();
					       // line does not exist
	      AssertThrow (false, ExcLineInexistant(line_vertices.first,
						    line_vertices.second));
	      line = end_line();
	    }
	}
				       // Assert that only exterior
				       // lines are given a boundary
				       // indicator
      if (!line->at_boundary())
	{
	  clear_despite_subscriptions ();

	  AssertThrow (false, ExcInteriorLineCantBeBoundary());
	}

                                       // and make sure that we don't
                                       // attempt to reset the
                                       // boundary indicator to a
                                       // different than the
                                       // previously set value
      if (line->boundary_indicator() != 0)
	{
	  if (line->boundary_indicator() != boundary_line->material_id)
	    {
	      clear_despite_subscriptions ();

	      AssertThrow (false, ExcMessage (
		"Duplicate boundary lines are only allowed "
		"if they carry the same boundary indicator."));
	    }
	}
      line->set_boundary_indicator (boundary_line->material_id);
    }

  
				   // now go on with boundary faces
  std::vector<CellData<2> >::const_iterator boundary_quad
    = subcelldata.boundary_quads.begin();
  std::vector<CellData<2> >::const_iterator end_boundary_quad
    = subcelldata.boundary_quads.end();
  for (; boundary_quad!=end_boundary_quad; ++boundary_quad)
    {
      quad_iterator quad;
      line_iterator line[4];
 
				       // first find the lines that
				       // are made up of the given
				       // vertices, then build up a
				       // quad from these lines
				       // finally use the find
				       // function of the map template
				       // to find the quad
      for (unsigned int i=0; i<4; ++i)
	{
	  std::pair<int, int> line_vertices(
	    boundary_quad->vertices[GeometryInfo<dim-1>::line_to_cell_vertices(i,0)],
	    boundary_quad->vertices[GeometryInfo<dim-1>::line_to_cell_vertices(i,1)]);
	  
					   // check whether line
					   // already exists
	  if (needed_lines.find(line_vertices) != needed_lines.end())
	    line[i] = needed_lines[line_vertices];
	  else
					     // look wether it exists
					     // in reverse direction
	    {
	      std::swap (line_vertices.first, line_vertices.second);
	      if (needed_lines.find(line_vertices) != needed_lines.end())
		line[i] = needed_lines[line_vertices];
	      else
		{
		  clear_despite_subscriptions ();
						   // line does not
						   // exist
		  AssertThrow (false, ExcLineInexistant(line_vertices.first,
							line_vertices.second));
		  line[i] = end_line();
		}
	    }
	}
      
 
				       // Set up 2 quads that are
				       // built up from the lines for
				       // reasons of comparison to
				       // needed_quads.  The second
				       // quad is the reversed version
				       // of the first quad in order
				       // find the quad regardless of
				       // its orientation.  This is
				       // introduced for convenience
				       // and because boundary quad
				       // orientation does not carry
				       // any information.
      internal::Triangulation::Quad
        quad_compare_1(line[0]->index(), line[1]->index(),
                       line[2]->index(), line[3]->index());
      internal::Triangulation::Quad
        quad_compare_2(line[2]->index(), line[3]->index(),
                       line[0]->index(), line[1]->index());
      
				       // try to find the quad with
				       // lines situated as
				       // constructed above.  if it
				       // could not be found, rotate
				       // the boundary lines 3 times
				       // until it is found or it does
				       // not exist.

				       // mapping from counterclock to
				       // lexicographic ordering of
				       // quad lines
      static const unsigned int lex2cclock[4]={3,1,0,2};
				       // copy lines from
				       // lexicographic to
				       // counterclock ordering, as
				       // rotation is much simpler in
				       // counterclock ordering
      line_iterator line_counterclock[4];
      for (unsigned int i=0; i<4; ++i)
	line_counterclock[lex2cclock[i]]=line[i];
      unsigned int n_rotations=0;
      bool not_found_quad_1;
      while ( (not_found_quad_1=(needed_quads.find(quad_compare_1) == needed_quads.end())) && 
	      (                  needed_quads.find(quad_compare_2) == needed_quads.end()) &&
	      (n_rotations<4))
 	{
					   // use the rotate defined
					   // in <algorithms>
 	  rotate(line_counterclock, line_counterclock+1, line_counterclock+4);
					   // update the quads with
					   // rotated lines (i runs in
					   // lexicographic ordering)
	  for (unsigned int i=0; i<4; ++i)
	    {
	      quad_compare_1.set_line(i,       line_counterclock[lex2cclock[i]]->index());
	      quad_compare_2.set_line((i+2)%4, line_counterclock[lex2cclock[i]]->index());
	    }
 	      
 	  ++n_rotations;
 	}
 
      if (n_rotations==4)
	{
	  clear_despite_subscriptions();

 	                        // quad does not exist
 	  AssertThrow (false, ExcQuadInexistant(line[0]->index(), line[1]->index(),
 	                                        line[2]->index(), line[3]->index()));
 	  quad = end_quad();
 	}
 
				       // exception not thrown,
				       // therefore assign the quad
				       // appropriately
      if (not_found_quad_1)
	quad = needed_quads[quad_compare_2].first;
      else
	quad = needed_quads[quad_compare_1].first;

				       // check whether this face is
				       // really an exterior one
      if (!quad->at_boundary())
	{
	  clear_despite_subscriptions ();

	  AssertThrow(false, ExcInteriorQuadCantBeBoundary());
	}

                                       // and make sure that we don't
                                       // attempt to reset the
                                       // boundary indicator to a
                                       // different than the
                                       // previously set value
      if (quad->boundary_indicator() != 0)
	{
	  if (quad->boundary_indicator() != boundary_quad->material_id)
	    {
	      clear_despite_subscriptions ();
	  
	      AssertThrow (false, ExcMessage (
		"Duplicate boundary quads are only allowed "
		"if they carry the same boundary indicator."));
	    }
	}
      quad->set_boundary_indicator (boundary_quad->material_id); 
    }


				   /////////////////////////////////////////
				   // finally update neighborship info
  for (cell_iterator cell=begin(); cell!=end(); ++cell)
    for (unsigned int face=0; face<6; ++face)
      if (adjacent_cells[cell->quad(face)->index()][0] == cell)
					 // first adjacent cell is
					 // this one
	{
	  if (adjacent_cells[cell->quad(face)->index()].size() == 2)
					     // there is another
					     // adjacent cell
	    cell->set_neighbor (face,
				adjacent_cells[cell->quad(face)->index()][1]);
	}
				   // first adjacent cell is not this
				   // one, -> it must be the neighbor
				   // we are looking for
      else
	cell->set_neighbor (face,
			    adjacent_cells[cell->quad(face)->index()][0]);


				   // re-compute numbers of lines, etc
  update_number_cache ();  
}

#endif


#if deal_II_dimension == 1

template <>
void Triangulation<1>::distort_random (const double factor,
				       const bool   keep_boundary)
{
				   // this function is mostly
				   // equivalent to that for the
				   // general dimensional case the
				   // only difference being the
				   // correction for split faces which
				   // is not necessary in 1D
				   //
				   // if you change something here,
				   // don't forget to do so there as
				   // well

  const unsigned int dim = 1;
  
				   // find the smallest length of the
				   // lines adjacent to the
				   // vertex. take the initial value
				   // to be larger than anything that
				   // might be found: the diameter of
				   // the triangulation, here computed
				   // by adding up the diameters of
				   // the coarse grid cells.
  double almost_infinite_length = 0;
  for (cell_iterator cell=begin(0); cell!=end(0); ++cell)
    almost_infinite_length += cell->diameter();
  
  std::vector<double> minimal_length (vertices.size(),
				      almost_infinite_length);
				   // also note if a vertex is at
				   // the boundary
  std::vector<bool>   at_boundary (vertices.size(), false);
  
  for (active_line_iterator line=begin_active_line();
       line != end_line(); ++line)
    {
      if (keep_boundary && line->at_boundary())
	{
	  at_boundary[line->vertex_index(0)] = true;
	  at_boundary[line->vertex_index(1)] = true;
	}
      
      minimal_length[line->vertex_index(0)]
	= std::min(line->diameter(),
		   minimal_length[line->vertex_index(0)]);
      minimal_length[line->vertex_index(1)]
	= std::min(line->diameter(),
		   minimal_length[line->vertex_index(1)]);
    }


  const unsigned int n_vertices = vertices.size();
  Point<dim> shift_vector;
  
  for (unsigned int vertex=0; vertex<n_vertices; ++vertex) 
    {
				       // ignore this vertex if we
				       // whall keep the boundary and
				       // this vertex *is* at the
				       // boundary
      if (keep_boundary && at_boundary[vertex])
	continue;
      
				       // first compute a random shift
				       // vector
      for (unsigned int d=0; d<dim; ++d)
	shift_vector(d) = std::rand()*1.0/RAND_MAX;

      shift_vector *= factor * minimal_length[vertex] /
		      std::sqrt(shift_vector.square());

				       // finally move the vertex
      vertices[vertex] += shift_vector;
    }
}

#endif


template <int dim>
void Triangulation<dim>::distort_random (const double factor,
					 const bool   keep_boundary)
{
//TODO:[?]Implement the random distortion in Triangulation for hanging nodes as well
//    Hanging nodes need to be reset to the correct mean value
//    at the end, which is simple for 2D but difficult for 3D. Maybe take
//    a look at how we get to the original location of the point in the
//    execute_refinement function and copy the relevant lines.
  
				   // this function is mostly
				   // equivalent to that for the
				   // general dimensional case the
				   // only difference being the
				   // correction for split faces which
				   // is not necessary in 1D
				   //
				   // if you change something here,
				   // don't forget to do so there as
				   // well
  
				   // find the smallest length of the
				   // lines adjacent to the
				   // vertex. take the initial value
				   // to be larger than anything that
				   // might be found: the diameter of
				   // the triangulation, here
				   // estimated by adding up the
				   // diameters of the coarse grid
				   // cells.
  double almost_infinite_length = 0;
  for (cell_iterator cell=begin(0); cell!=end(0); ++cell)
    almost_infinite_length += cell->diameter();
  
  std::vector<double> minimal_length (vertices.size(),
				      almost_infinite_length);

				   // also note if a vertex is at the
				   // boundary
  std::vector<bool>   at_boundary (vertices.size(), false);
  
  for (active_line_iterator line=begin_active_line();
       line != end_line(); ++line)
    {
      if (keep_boundary && line->at_boundary())
	{
	  at_boundary[line->vertex_index(0)] = true;
	  at_boundary[line->vertex_index(1)] = true;
	}
      
      minimal_length[line->vertex_index(0)]
	= std::min(line->diameter(),
		   minimal_length[line->vertex_index(0)]);
      minimal_length[line->vertex_index(1)]
	= std::min(line->diameter(),
		   minimal_length[line->vertex_index(1)]);
    }


  const unsigned int n_vertices = vertices.size();
  Point<dim> shift_vector;
  
  for (unsigned int vertex=0; vertex<n_vertices; ++vertex) 
    {
				       // ignore this vertex if we
				       // whall keep the boundary and
				       // this vertex *is* at the
				       // boundary
      if (keep_boundary && at_boundary[vertex])
	continue;
      
				       // first compute a random shift
				       // vector
      for (unsigned int d=0; d<dim; ++d)
	shift_vector(d) = std::rand()*1.0/RAND_MAX;

      shift_vector *= factor * minimal_length[vertex] /
		      std::sqrt(shift_vector.square());

				       // finally move the vertex
      vertices[vertex] += shift_vector;
    }


				   // finally correct hanging nodes
				   // again. The following is not
				   // necessary for 1D
  active_cell_iterator cell = begin_active(),
		       endc = end();
  for (; cell!=endc; ++cell) 
    for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
      if (cell->face(face)->has_children() &&
	  !cell->face(face)->at_boundary())
					 // this lines has children,
					 // thus there are restricted
					 // nodes
	{
					   // not implemented at
					   // present for dim=3 or
					   // higher
	  Assert (dim<=2, ExcNotImplemented());

					   // compute where the common
					   // point of the two child
					   // lines will lie and reset
					   // it to the correct value
	  vertices[cell->face(face)->child(0)->vertex_index(1)]
	    = (cell->face(face)->vertex(0) +
	       cell->face(face)->vertex(1)) / 2;
	}
}



template <int dim>
void Triangulation<dim>::set_all_refine_flags ()
{
  active_cell_iterator cell = begin_active(),
		       endc = end();

  for (; cell != endc; ++cell)
    {
      cell->clear_coarsen_flag();
      cell->set_refine_flag ();
    }
}



template <int dim>
void Triangulation<dim>::refine_global (const unsigned int times)
{
  for (unsigned int i=0; i<times; ++i)
    {
      set_all_refine_flags();
      execute_coarsening_and_refinement ();
    }
}



template <int dim>
void Triangulation<dim>::save_refine_flags (std::vector<bool> &v) const
{
  v.resize (n_active_cells(), false);
  std::vector<bool>::iterator  i = v.begin();
  active_cell_iterator cell = begin_active(),
		       endc = end();
  for (; cell!=endc; ++cell, ++i)
    *i = cell->refine_flag_set();
}



template <int dim>
void Triangulation<dim>::save_refine_flags (std::ostream &out) const
{
  std::vector<bool> v;
  save_refine_flags (v);
  write_bool_vector (mn_tria_refine_flags_begin, v, mn_tria_refine_flags_end,
		     out);
}



template <int dim>
void Triangulation<dim>::load_refine_flags (std::istream &in)
{
  std::vector<bool> v;
  read_bool_vector (mn_tria_refine_flags_begin, v, mn_tria_refine_flags_end,
		    in);
  load_refine_flags (v);
}



template <int dim>
void Triangulation<dim>::load_refine_flags (const std::vector<bool> &v)
{
  AssertThrow (v.size() == n_active_cells(), ExcGridReadError());
  
  active_cell_iterator cell = begin_active(),
		       endc = end();
  std::vector<bool>::const_iterator i = v.begin();
  for (; cell!=endc; ++cell, ++i)
    if (*i == true)
      cell->set_refine_flag();
    else
      cell->clear_refine_flag();
}



template <int dim>
void Triangulation<dim>::save_coarsen_flags (std::vector<bool> &v) const
{
  v.resize (n_active_cells(), false);
  std::vector<bool>::iterator  i = v.begin();
  active_cell_iterator cell = begin_active(),
		       endc = end();
  for (; cell!=endc; ++cell, ++i)
    *i = cell->coarsen_flag_set();
}



template <int dim>
void Triangulation<dim>::save_coarsen_flags (std::ostream &out) const
{
  std::vector<bool> v;
  save_coarsen_flags (v);
  write_bool_vector (mn_tria_coarsen_flags_begin, v, mn_tria_coarsen_flags_end,
		     out);
}



template <int dim>
void Triangulation<dim>::load_coarsen_flags (std::istream &in)
{
  std::vector<bool> v;
  read_bool_vector (mn_tria_coarsen_flags_begin, v, mn_tria_coarsen_flags_end,
		    in);
  load_coarsen_flags (v);
}



template <int dim>
void Triangulation<dim>::load_coarsen_flags (const std::vector<bool> &v)
{
  Assert (v.size() == n_active_cells(), ExcGridReadError());
  
  active_cell_iterator cell = begin_active(),
		       endc = end();
  std::vector<bool>::const_iterator i = v.begin();
  for (; cell!=endc; ++cell, ++i)
    if (*i == true)
      cell->set_coarsen_flag();
    else
      cell->clear_coarsen_flag();
}


#if deal_II_dimension == 1

template <>
void Triangulation<1>::clear_user_pointers ()
{
  cell_iterator cell = begin(),
		endc = end();
  for (; cell!=endc; ++cell)
    cell->clear_user_pointer ();
}



template <>
void Triangulation<1>::clear_user_flags ()
{
  clear_user_flags_line();
}



template <>
void Triangulation<1>::clear_user_flags_quad ()
{}


template <>
void Triangulation<1>::clear_user_flags_hex ()
{}


#endif


#if deal_II_dimension == 2

template <>
void Triangulation<2>::clear_user_pointers ()
{
  line_iterator line = begin_line(),
		endl = end_line();
  for (; line!=endl; ++line)
    line->clear_user_pointer ();

  cell_iterator cell = begin(),
		endc = end();
  for (; cell!=endc; ++cell)
    cell->clear_user_pointer ();
}



template <>
void Triangulation<2>::clear_user_flags ()
{
  clear_user_flags_line ();
  clear_user_flags_quad ();  
}



template <>
void Triangulation<2>::clear_user_flags_hex ()
{}


#endif


#if deal_II_dimension == 3

template <>
void Triangulation<3>::clear_user_pointers ()
{
  line_iterator line = begin_line(),
		endl = end_line();
  for (; line!=endl; ++line)
    line->clear_user_pointer ();

  quad_iterator quad = begin_quad(),
		endq = end_quad();
  for (; quad!=endq; ++quad)
    quad->clear_user_pointer ();

  cell_iterator cell = begin(),
		endc = end();
  for (; cell!=endc; ++cell)
    cell->clear_user_pointer ();
}



template <>
void Triangulation<3>::clear_user_flags ()
{
  clear_user_flags_line ();
  clear_user_flags_quad ();
  clear_user_flags_hex ();
}



#endif


template <int dim>
void Triangulation<dim>::clear_user_flags_line ()
{
  line_iterator line = begin_line(),
		endl = end_line();
  for (; line!=endl; ++line)
    line->clear_user_flag ();
}



template <int dim>
void Triangulation<dim>::clear_user_flags_quad ()
{
  quad_iterator quad = begin_quad(),
		endq = end_quad();
  for (; quad!=endq; ++quad)
    quad->clear_user_flag ();
}



template <int dim>
void Triangulation<dim>::clear_user_flags_hex ()
{
  hex_iterator hex  = begin_hex(),
               endh = end_hex();
  for (; hex!=endh; ++hex)
    hex->clear_user_flag ();
}



template <int dim>
void Triangulation<dim>::save_user_flags (std::ostream &out) const
{
  save_user_flags_line (out);
  
  if (dim>=2)
    save_user_flags_quad (out);
  
  if (dim>=3)
    save_user_flags_hex (out);

  if (dim >= 4)
    Assert (false, ExcNotImplemented());
}



template <int dim>
void Triangulation<dim>::save_user_flags (std::vector<bool> &v) const
{
				   // clear vector and append
				   // all the stuff later on
  v.clear ();

  std::vector<bool> tmp;

  save_user_flags_line (tmp);
  v.insert (v.end(), tmp.begin(), tmp.end());

  if (dim >= 2)
    {
      save_user_flags_quad (tmp);
      v.insert (v.end(), tmp.begin(), tmp.end());
    }
  
  if (dim >= 3)
    {
      save_user_flags_hex (tmp);
      v.insert (v.end(), tmp.begin(), tmp.end());
    }      

  if (dim >= 4)
    Assert (false, ExcNotImplemented());
}



template <int dim>
void Triangulation<dim>::load_user_flags (std::istream &in)
{
  load_user_flags_line (in);
  
  if (dim>=2)
    load_user_flags_quad (in);
  
  if (dim>=3)
    load_user_flags_hex (in);

  if (dim >= 4)
    Assert (false, ExcNotImplemented());
}



template <int dim>
void Triangulation<dim>::load_user_flags (const std::vector<bool> &v)
{
  Assert (v.size() == n_lines()+n_quads()+n_hexs(), ExcInternalError());
  std::vector<bool> tmp;

				   // first extract the flags
				   // belonging to lines
  tmp.insert (tmp.end(),
	      v.begin(), v.begin()+n_lines());
				   // and set the lines
  load_user_flags_line (tmp);

  if (dim >= 2)
    {
      tmp.clear ();
      tmp.insert (tmp.end(),
		  v.begin()+n_lines(), v.begin()+n_lines()+n_quads());
      load_user_flags_quad (tmp);
    }
  
  if (dim >= 3)
    {
      tmp.clear();
      tmp.insert (tmp.end(),
		  v.begin()+n_lines()+n_quads(), v.begin()+n_lines()+n_quads()+n_hexs());
      load_user_flags_hex (tmp);
    }      

  if (dim >= 4)
    Assert (false, ExcNotImplemented());
}



template <int dim>
void Triangulation<dim>::save_user_flags_line (std::vector<bool> &v) const
{
  v.resize (n_lines(), false);
  std::vector<bool>::iterator  i = v.begin();
  line_iterator line = begin_line(),
		endl = end_line();
  for (; line!=endl; ++line, ++i)
    *i = line->user_flag_set();
}



template <int dim>
void Triangulation<dim>::save_user_flags_line (std::ostream &out) const
{
  std::vector<bool> v;
  save_user_flags_line (v);
  write_bool_vector (mn_tria_line_user_flags_begin, v, mn_tria_line_user_flags_end,
		     out);
}



template <int dim>
void Triangulation<dim>::load_user_flags_line (std::istream &in)
{
  std::vector<bool> v;
  read_bool_vector (mn_tria_line_user_flags_begin, v, mn_tria_line_user_flags_end,
		    in);
  load_user_flags_line (v);
}



template <int dim>
void Triangulation<dim>::load_user_flags_line (const std::vector<bool> &v)
{
  Assert (v.size() == n_lines(), ExcGridReadError());
  
  line_iterator line = begin_line(),
		endl = end_line();
  std::vector<bool>::const_iterator i = v.begin();
  for (; line!=endl; ++line, ++i)
    if (*i == true)
      line->set_user_flag();
    else
      line->clear_user_flag();
}


#if deal_II_dimension == 1

template <>
void Triangulation<1>::save_user_flags_quad (std::ostream &) const
{
  Assert (false, ExcImpossibleInDim(1));
}



template <>
void Triangulation<1>::save_user_flags_quad (std::vector<bool> &) const
{
  Assert (false, ExcImpossibleInDim(1));
}



template <>
void Triangulation<1>::load_user_flags_quad (std::istream &)
{
  Assert (false, ExcImpossibleInDim(1));
}



template <>
void Triangulation<1>::load_user_flags_quad (const std::vector<bool> &)
{
  Assert (false, ExcImpossibleInDim(1));
}



template <>
void Triangulation<1>::save_user_flags_hex (std::ostream &) const
{
  Assert (false, ExcImpossibleInDim(1));
}



template <>
void Triangulation<1>::save_user_flags_hex (std::vector<bool> &) const
{
  Assert (false, ExcImpossibleInDim(1));
}



template <>
void Triangulation<1>::load_user_flags_hex (std::istream &)
{
  Assert (false, ExcImpossibleInDim(1));
}



template <>
void Triangulation<1>::load_user_flags_hex (const std::vector<bool> &)
{
  Assert (false, ExcImpossibleInDim(1));
}

#endif


template <int dim>
void Triangulation<dim>::save_user_flags_quad (std::vector<bool> &v) const
{
  v.resize (n_quads(), false);
  std::vector<bool>::iterator  i = v.begin();
  quad_iterator quad = begin_quad(),
		endq = end_quad();
  for (; quad!=endq; ++quad, ++i)
    *i = quad->user_flag_set();
}



template <int dim>
void Triangulation<dim>::save_user_flags_quad (std::ostream &out) const
{
  std::vector<bool> v;
  save_user_flags_quad (v);
  write_bool_vector (mn_tria_quad_user_flags_begin, v, mn_tria_quad_user_flags_end,
		     out);
}



template <int dim>
void Triangulation<dim>::load_user_flags_quad (std::istream &in)
{
  std::vector<bool> v;
  read_bool_vector (mn_tria_quad_user_flags_begin, v, mn_tria_quad_user_flags_end,
		    in);
  load_user_flags_quad (v);
}



template <int dim>
void Triangulation<dim>::load_user_flags_quad (const std::vector<bool> &v)
{
  Assert (v.size() == n_quads(), ExcGridReadError());
  
  quad_iterator quad = begin_quad(),
		endq = end_quad();
  std::vector<bool>::const_iterator i = v.begin();
  for (; quad!=endq; ++quad, ++i)
    if (*i == true)
      quad->set_user_flag();
    else
      quad->clear_user_flag();
}


#if deal_II_dimension == 2

template <>
void Triangulation<2>::save_user_flags_hex (std::ostream &) const
{
  Assert (false, ExcImpossibleInDim(2));
}



template <>
void Triangulation<2>::save_user_flags_hex (std::vector<bool> &) const
{
  Assert (false, ExcImpossibleInDim(2));
}



template <>
void Triangulation<2>::load_user_flags_hex (std::istream &)
{
  Assert (false, ExcImpossibleInDim(2));
}



template <>
void Triangulation<2>::load_user_flags_hex (const std::vector<bool> &)
{
  Assert (false, ExcImpossibleInDim(2));
}


#endif


template <int dim>
void Triangulation<dim>::save_user_flags_hex (std::vector<bool> &v) const
{
  v.resize (n_hexs(), false);
  std::vector<bool>::iterator  i = v.begin();
  hex_iterator hex = begin_hex(),
	      endh = end_hex();
  for (; hex!=endh; ++hex, ++i)
    *i = hex->user_flag_set();
}



template <int dim>
void Triangulation<dim>::save_user_flags_hex (std::ostream &out) const
{
  std::vector<bool> v;
  save_user_flags_hex (v);
  write_bool_vector (mn_tria_hex_user_flags_begin, v, mn_tria_hex_user_flags_end,
		     out);
}



template <int dim>
void Triangulation<dim>::load_user_flags_hex (std::istream &in)
{
  std::vector<bool> v;
  read_bool_vector (mn_tria_hex_user_flags_begin, v, mn_tria_hex_user_flags_end,
		    in);
  load_user_flags_hex (v);
}



template <int dim>
void Triangulation<dim>::load_user_flags_hex (const std::vector<bool> &v)
{
  Assert (v.size() == n_hexs(), ExcGridReadError());
  
  hex_iterator hex = begin_hex(),
	      endh = end_hex();
  std::vector<bool>::const_iterator i = v.begin();
  for (; hex!=endh; ++hex, ++i)
    if (*i == true)
      hex->set_user_flag();
    else
      hex->clear_user_flag();
}



template <int dim>
void Triangulation<dim>::save_user_pointers (std::vector<void *> &v) const
{
				   // clear vector and append all the
				   // stuff later on
  v.clear ();

  std::vector<void *> tmp;

  save_user_pointers_line (tmp);
  v.insert (v.end(), tmp.begin(), tmp.end());

  if (dim >= 2)
    {
      save_user_pointers_quad (tmp);
      v.insert (v.end(), tmp.begin(), tmp.end());
    }
  
  if (dim >= 3)
    {
      save_user_pointers_hex (tmp);
      v.insert (v.end(), tmp.begin(), tmp.end());
    }      

  if (dim >= 4)
    Assert (false, ExcNotImplemented());
}



template <int dim>
void Triangulation<dim>::load_user_pointers (const std::vector<void *> &v)
{
  Assert (v.size() == n_lines()+n_quads()+n_hexs(), ExcInternalError());
  std::vector<void *> tmp;

				   // first extract the pointers
				   // belonging to lines
  tmp.insert (tmp.end(),
	      v.begin(), v.begin()+n_lines());
				   // and set the lines
  load_user_pointers_line (tmp);

  if (dim >= 2)
    {
      tmp.clear ();
      tmp.insert (tmp.end(),
		  v.begin()+n_lines(), v.begin()+n_lines()+n_quads());
      load_user_pointers_quad (tmp);
    }
  
  if (dim >= 3)
    {
      tmp.clear ();
      tmp.insert (tmp.end(),
		  v.begin()+n_lines()+n_quads(), v.begin()+n_lines()+n_quads()+n_hexs());
      load_user_pointers_hex (tmp);
    }      

  if (dim >= 4)
    Assert (false, ExcNotImplemented());
}



template <int dim>
void Triangulation<dim>::save_user_pointers_line (std::vector<void *> &v) const
{
  v.resize (n_lines(), 0);
  std::vector<void *>::iterator  i = v.begin();
  line_iterator line = begin_line(),
		endl = end_line();
  for (; line!=endl; ++line, ++i)
    *i = line->user_pointer();
}



template <int dim>
void Triangulation<dim>::load_user_pointers_line (const std::vector<void *> &v)
{
  Assert (v.size() == n_lines(), ExcGridReadError());
  
  line_iterator line = begin_line(),
		endl = end_line();
  std::vector<void *>::const_iterator i = v.begin();
  for (; line!=endl; ++line, ++i)
    line->set_user_pointer(*i);
}


#if deal_II_dimension == 1


template <>
void Triangulation<1>::save_user_pointers_quad (std::vector<void *> &) const
{
  Assert (false, ExcImpossibleInDim(1));
}



template <>
void Triangulation<1>::load_user_pointers_quad (const std::vector<void *> &)
{
  Assert (false, ExcImpossibleInDim(1));
}



template <>
void Triangulation<1>::save_user_pointers_hex (std::vector<void *> &) const
{
  Assert (false, ExcImpossibleInDim(1));
}



template <>
void Triangulation<1>::load_user_pointers_hex (const std::vector<void *> &)
{
  Assert (false, ExcImpossibleInDim(1));
}

#endif


template <int dim>
void Triangulation<dim>::save_user_pointers_quad (std::vector<void *> &v) const
{
  v.resize (n_quads(), 0);
  std::vector<void *>::iterator  i = v.begin();
  quad_iterator quad = begin_quad(),
		endq = end_quad();
  for (; quad!=endq; ++quad, ++i)
    *i = quad->user_pointer();
}



template <int dim>
void Triangulation<dim>::load_user_pointers_quad (const std::vector<void *> &v)
{
  Assert (v.size() == n_quads(), ExcGridReadError());
  
  quad_iterator quad = begin_quad(),
		endq = end_quad();
  std::vector<void *>::const_iterator i = v.begin();
  for (; quad!=endq; ++quad, ++i)
    quad->set_user_pointer(*i);
}


#if deal_II_dimension == 2



template <>
void Triangulation<2>::save_user_pointers_hex (std::vector<void *> &) const
{
  Assert (false, ExcImpossibleInDim(2));
}



template <>
void Triangulation<2>::load_user_pointers_hex (const std::vector<void *> &)
{
  Assert (false, ExcImpossibleInDim(2));
}


#endif


template <int dim>
void Triangulation<dim>::save_user_pointers_hex (std::vector<void *> &v) const
{
  v.resize (n_hexs(), 0);
  std::vector<void *>::iterator  i = v.begin();
  hex_iterator hex = begin_hex(),
	      endh = end_hex();
  for (; hex!=endh; ++hex, ++i)
    *i = hex->user_pointer();
}



template <int dim>
void Triangulation<dim>::load_user_pointers_hex (const std::vector<void *> &v)
{
  Assert (v.size() == n_hexs(), ExcGridReadError());
  
  hex_iterator hex = begin_hex(),
	      endh = end_hex();
  std::vector<void *>::const_iterator i = v.begin();
  for (; hex!=endh; ++hex, ++i)
    hex->set_user_pointer(*i);
}



/*------------------------ Iterator functions ------------------------*/


#if deal_II_dimension == 1

template <>
Triangulation<1>::raw_cell_iterator
Triangulation<1>::begin_raw (const unsigned int level) const
{
  return begin_raw_line (level);
}



template <>
Triangulation<1>::active_cell_iterator
Triangulation<1>::begin_active (const unsigned int level) const
{
  return begin_active_line (level);
}



template <>
Triangulation<1>::raw_cell_iterator
Triangulation<1>::last_raw () const
{
  return last_raw_line ();
}



template <>
Triangulation<1>::raw_cell_iterator
Triangulation<1>::last_raw (const unsigned int level) const
{
  return last_raw_line (level);
}



template <>
Triangulation<1>::cell_iterator
Triangulation<1>::last () const
{
  return last_line ();
}



template <>
Triangulation<1>::cell_iterator
Triangulation<1>::last (const unsigned int level) const
{
  return last_line (level);
}



template <>
Triangulation<1>::active_cell_iterator
Triangulation<1>::last_active () const
{
  return last_active_line ();
}



template <>
Triangulation<1>::active_cell_iterator
Triangulation<1>::last_active (const unsigned int level) const
{
  return last_active_line (level);
}



template <>
Triangulation<1>::raw_face_iterator
Triangulation<1>::begin_raw_face () const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}



template <>
Triangulation<1>::face_iterator
Triangulation<1>::begin_face () const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}



template <>
Triangulation<1>::active_face_iterator
Triangulation<1>::begin_active_face () const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}



template <>
Triangulation<1>::raw_face_iterator
Triangulation<1>::end_face () const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}



template <>
Triangulation<1>::raw_face_iterator
Triangulation<1>::last_raw_face () const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}



template <>
Triangulation<1>::face_iterator
Triangulation<1>::last_face () const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}



template <>
Triangulation<1>::active_face_iterator
Triangulation<1>::last_active_face () const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}



template <>
Triangulation<1>::raw_quad_iterator
Triangulation<1>::begin_raw_quad (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}



template <>
Triangulation<1>::quad_iterator
Triangulation<1>::begin_quad (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}



template <>
Triangulation<1>::active_quad_iterator
Triangulation<1>::begin_active_quad (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}



template <>
Triangulation<1>::raw_quad_iterator
Triangulation<1>::end_quad () const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}



template <>
Triangulation<1>::raw_quad_iterator
Triangulation<1>::last_raw_quad (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}



template <>
Triangulation<1>::raw_quad_iterator
Triangulation<1>::last_raw_quad () const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}



template <>
Triangulation<1>::quad_iterator
Triangulation<1>::last_quad (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}



template <>
Triangulation<1>::quad_iterator
Triangulation<1>::last_quad () const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}



template <>
Triangulation<1>::active_quad_iterator
Triangulation<1>::last_active_quad (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}



template <>
Triangulation<1>::active_quad_iterator
Triangulation<1>::last_active_quad () const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}



template <>
Triangulation<1>::raw_hex_iterator
Triangulation<1>::begin_raw_hex (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}



template <>
Triangulation<1>::hex_iterator
Triangulation<1>::begin_hex (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}



template <>
Triangulation<1>::active_hex_iterator
Triangulation<1>::begin_active_hex (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}



template <>
Triangulation<1>::raw_hex_iterator
Triangulation<1>::end_hex () const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}



template <>
Triangulation<1>::raw_hex_iterator
Triangulation<1>::last_raw_hex (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}



template <>
Triangulation<1>::raw_hex_iterator
Triangulation<1>::last_raw_hex () const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}



template <>
Triangulation<1>::hex_iterator
Triangulation<1>::last_hex (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}



template <>
Triangulation<1>::hex_iterator
Triangulation<1>::last_hex () const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}



template <>
Triangulation<1>::active_hex_iterator
Triangulation<1>::last_active_hex (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}



template <>
Triangulation<1>::active_hex_iterator
Triangulation<1>::last_active_hex () const
{
  Assert (false, ExcImpossibleInDim(1));
  return 0;
}

#endif


#if deal_II_dimension == 2

template <>
Triangulation<2>::raw_cell_iterator
Triangulation<2>::begin_raw (const unsigned int level) const
{
  return begin_raw_quad (level);
}



template <>
Triangulation<2>::active_cell_iterator
Triangulation<2>::begin_active (const unsigned int level) const
{
  return begin_active_quad (level);
}



template <>
Triangulation<2>::raw_cell_iterator
Triangulation<2>::last_raw () const
{
  return last_raw_quad ();
}



template <>
Triangulation<2>::raw_cell_iterator
Triangulation<2>::last_raw (const unsigned int level) const
{
  return last_raw_quad (level);
}



template <>
Triangulation<2>::cell_iterator
Triangulation<2>::last () const
{
  return last_quad ();
}



template <>
Triangulation<2>::cell_iterator
Triangulation<2>::last (const unsigned int level) const
{
  return last_quad (level);
}



template <>
Triangulation<2>::active_cell_iterator
Triangulation<2>::last_active () const
{
  return last_active_quad ();
}



template <>
Triangulation<2>::active_cell_iterator
Triangulation<2>::last_active (const unsigned int level) const
{
  return last_active_quad (level);
}



template <>
Triangulation<2>::raw_face_iterator
Triangulation<2>::begin_raw_face () const
{
  return begin_raw_line ();
}



template <>
Triangulation<2>::face_iterator
Triangulation<2>::begin_face () const
{
  return begin_line ();
}



template <>
Triangulation<2>::active_face_iterator
Triangulation<2>::begin_active_face () const
{
  return begin_active_line ();
}



template <>
Triangulation<2>::raw_face_iterator
Triangulation<2>::end_face () const
{
  return end_line ();
}



template <>
Triangulation<2>::raw_face_iterator
Triangulation<2>::last_raw_face () const
{
  return last_raw_line ();
}



template <>
Triangulation<2>::face_iterator
Triangulation<2>::last_face () const
{
  return last_line ();
}



template <>
Triangulation<2>::active_face_iterator
Triangulation<2>::last_active_face () const
{
  return last_active_line ();
}



template <>
Triangulation<2>::raw_hex_iterator
Triangulation<2>::begin_raw_hex (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(2));
  return 0;
}



template <>
Triangulation<2>::hex_iterator
Triangulation<2>::begin_hex (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(2));
  return 0;
}



template <>
Triangulation<2>::active_hex_iterator
Triangulation<2>::begin_active_hex (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(2));
  return 0;
}



template <>
Triangulation<2>::raw_hex_iterator
Triangulation<2>::end_hex () const
{
  Assert (false, ExcImpossibleInDim(2));
  return 0;
}



template <>
Triangulation<2>::raw_hex_iterator
Triangulation<2>::last_raw_hex (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(2));
  return 0;
}



template <>
Triangulation<2>::raw_hex_iterator
Triangulation<2>::last_raw_hex () const
{
  Assert (false, ExcImpossibleInDim(2));
  return 0;
}



template <>
Triangulation<2>::hex_iterator
Triangulation<2>::last_hex (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(2));
  return 0;
}



template <>
Triangulation<2>::hex_iterator
Triangulation<2>::last_hex () const
{
  Assert (false, ExcImpossibleInDim(2));
  return 0;
}



template <>
Triangulation<2>::active_hex_iterator
Triangulation<2>::last_active_hex (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(2));
  return 0;
}



template <>
Triangulation<2>::active_hex_iterator
Triangulation<2>::last_active_hex () const
{
  Assert (false, ExcImpossibleInDim(2));
  return 0;
}

#endif


#if deal_II_dimension == 3

template <>
Triangulation<3>::raw_cell_iterator
Triangulation<3>::begin_raw (const unsigned int level) const
{
  return begin_raw_hex (level);
}



template <>
Triangulation<3>::active_cell_iterator
Triangulation<3>::begin_active (const unsigned int level) const
{
  return begin_active_hex (level);
}



template <>
Triangulation<3>::raw_cell_iterator
Triangulation<3>::last_raw () const
{
  return last_raw_hex ();
}



template <>
Triangulation<3>::raw_cell_iterator
Triangulation<3>::last_raw (const unsigned int level) const
{
  return last_raw_hex (level);
}



template <>
Triangulation<3>::cell_iterator
Triangulation<3>::last () const
{
  return last_hex ();
}



template <>
Triangulation<3>::cell_iterator
Triangulation<3>::last (const unsigned int level) const
{
  return last_hex (level);
}



template <>
Triangulation<3>::active_cell_iterator
Triangulation<3>::last_active () const
{
  return last_active_hex ();
}



template <>
Triangulation<3>::active_cell_iterator
Triangulation<3>::last_active (const unsigned int level) const
{
  return last_active_hex (level);
}



template <>
Triangulation<3>::raw_face_iterator
Triangulation<3>::begin_raw_face () const
{
  return begin_raw_quad ();
}



template <>
Triangulation<3>::face_iterator
Triangulation<3>::begin_face () const
{
  return begin_quad ();
}



template <>
Triangulation<3>::active_face_iterator
Triangulation<3>::begin_active_face () const
{
  return begin_active_quad ();
}



template <>
Triangulation<3>::raw_face_iterator
Triangulation<3>::end_face () const
{
  return end_quad ();
}



template <>
Triangulation<3>::raw_face_iterator
Triangulation<3>::last_raw_face () const
{
  return last_raw_quad ();
}



template <>
Triangulation<3>::face_iterator
Triangulation<3>::last_face () const
{
  return last_quad ();
}



template <>
Triangulation<3>::active_face_iterator
Triangulation<3>::last_active_face () const
{
  return last_active_quad ();
}


#endif


template <int dim>
typename Triangulation<dim>::raw_line_iterator
Triangulation<dim>::begin_raw_line (const unsigned int level) const
{
  Assert (level == 0, ExcFacesHaveNoLevel());
  return raw_line_iterator (const_cast<Triangulation<dim>*>(this),
			    0,
			    0);
}

#if deal_II_dimension == 1

template <>
Triangulation<1>::raw_line_iterator
Triangulation<1>::begin_raw_line (const unsigned int level) const
{
  Assert (level<levels.size(), ExcInvalidLevel(level));

  if (levels[level]->lines.cells.size() == 0)
    return end_line ();

  return raw_line_iterator (const_cast<Triangulation<1>*>(this),
			    level,
			    0);
}

#endif

#if deal_II_dimension == 3

template <>
Triangulation<3>::raw_quad_iterator
Triangulation<3>::begin_raw_quad (const unsigned int level) const
{
  Assert (level == 0, ExcFacesHaveNoLevel())

  return raw_quad_iterator (const_cast<Triangulation<3>*>(this),
			    0,
			    0);
}

#endif

#if deal_II_dimension == 2

template <>
Triangulation<2>::raw_quad_iterator
Triangulation<2>::begin_raw_quad (const unsigned int level) const
{
  Assert (level<levels.size(), ExcInvalidLevel(level));

  if (levels[level]->quads.cells.size() == 0)
    return end_quad();
  
  return raw_quad_iterator (const_cast<Triangulation<2>*>(this),
			    level,
			    0);
}

#endif


template <int dim>
typename Triangulation<dim>::raw_hex_iterator
Triangulation<dim>::begin_raw_hex (const unsigned int level) const
{
  Assert (level<levels.size(), ExcInvalidLevel(level));

  if (levels[level]->hexes.cells.size() == 0)
    return end_hex();
  
  return raw_hex_iterator (const_cast<Triangulation<dim>*>(this),
			   level,
			   0);
}



template <int dim>
typename Triangulation<dim>::line_iterator
Triangulation<dim>::begin_line (const unsigned int level) const
{
  				   // level is checked in begin_raw
  raw_line_iterator ri = begin_raw_line (level);
  if (ri.state() != IteratorState::valid)
    return ri;
  while (ri->used() == false)
    if ((++ri).state() != IteratorState::valid)
      return ri;
  return ri;
}



template <int dim>
typename Triangulation<dim>::quad_iterator
Triangulation<dim>::begin_quad (const unsigned int level) const
{
  				   // level is checked in begin_raw
  raw_quad_iterator ri = begin_raw_quad (level);
  if (ri.state() != IteratorState::valid)
    return ri;
  while (ri->used() == false)
    if ((++ri).state() != IteratorState::valid)
      return ri;
  return ri;
}



template <int dim>
typename Triangulation<dim>::hex_iterator
Triangulation<dim>::begin_hex (const unsigned int level) const
{
  				   // level is checked in begin_raw
  raw_hex_iterator ri = begin_raw_hex (level);
  if (ri.state() != IteratorState::valid)
    return ri;
  while (ri->used() == false)
    if ((++ri).state() != IteratorState::valid)
      return ri;
  return ri;
}



template <int dim>
typename Triangulation<dim>::active_line_iterator
Triangulation<dim>::begin_active_line (const unsigned int level) const
{
  				   // level is checked in begin_raw
  line_iterator i = begin_line (level);
  if (i.state() != IteratorState::valid)
    return i;
  while (i->has_children())
    if ((++i).state() != IteratorState::valid)
      return i;
  return i;
}



template <int dim>
typename Triangulation<dim>::active_quad_iterator
Triangulation<dim>::begin_active_quad (const unsigned int level) const
{
  				   // level is checked in begin_raw
  quad_iterator i = begin_quad (level);
  if (i.state() != IteratorState::valid)
    return i;
  while (i->has_children())
    if ((++i).state() != IteratorState::valid)
      return i;
  return i;
}



template <int dim>
typename Triangulation<dim>::active_hex_iterator
Triangulation<dim>::begin_active_hex (const unsigned int level) const
{
  				   // level is checked in begin_raw
  hex_iterator i = begin_hex (level);
  if (i.state() != IteratorState::valid)
    return i;
  while (i->has_children())
    if ((++i).state() != IteratorState::valid)
      return i;
  return i;
}



template <int dim>
typename Triangulation<dim>::raw_line_iterator
Triangulation<dim>::end_line () const
{
  return raw_line_iterator (const_cast<Triangulation<dim>*>(this),
			    -1,
			    -1);
}



template <int dim>
typename Triangulation<dim>::raw_quad_iterator
Triangulation<dim>::end_quad () const
{
  return raw_quad_iterator (const_cast<Triangulation<dim>*>(this),
			    -1,
			    -1);
}



template <int dim>
typename Triangulation<dim>::raw_hex_iterator
Triangulation<dim>::end_hex () const
{
  return raw_hex_iterator (const_cast<Triangulation<dim>*>(this),
			   -1,
			   -1);
}



template <int dim>
typename Triangulation<dim>::raw_line_iterator
Triangulation<dim>::last_raw_line (const unsigned int level) const
{
  Assert (level == 0, ExcFacesHaveNoLevel());

  return raw_line_iterator (const_cast<Triangulation<dim>*>(this),
			    0,
			    faces->lines.cells.size()-1);
  
}


#if deal_II_dimension == 1

template <>
Triangulation<1>::raw_line_iterator
Triangulation<1>::last_raw_line (const unsigned int level) const
{
  Assert (level<levels.size(), ExcInvalidLevel(level));
  Assert (levels[level]->lines.cells.size() != 0,
	  ExcEmptyLevel (level));

  return raw_line_iterator (const_cast<Triangulation<1>*>(this),
			    level,
			    levels[level]->lines.cells.size()-1);
}

#endif


template <int dim>
typename Triangulation<dim>::raw_quad_iterator
Triangulation<dim>::last_raw_quad (const unsigned int level) const
{
  Assert (level == 0, ExcFacesHaveNoLevel());
  return raw_quad_iterator (const_cast<Triangulation<dim>*>(this),
			    0,
			    faces->quads.cells.size()-1);
}


#if deal_II_dimension == 2

template <>
Triangulation<2>::raw_quad_iterator
Triangulation<2>::last_raw_quad (const unsigned int level) const
{
  Assert (level<levels.size(),
	  ExcInvalidLevel(level));
  Assert (levels[level]->quads.cells.size() != 0,
	  ExcEmptyLevel (level));
  return raw_quad_iterator (const_cast<Triangulation<2>*>(this),
			    level,
			    levels[level]->quads.cells.size()-1);
}

#endif


template <int dim>
typename Triangulation<dim>::raw_hex_iterator
Triangulation<dim>::last_raw_hex (const unsigned int level) const
{
  Assert (level<levels.size(),
	  ExcInvalidLevel(level));
  Assert (levels[level]->hexes.cells.size() != 0,
	  ExcEmptyLevel (level));

  return raw_hex_iterator (const_cast<Triangulation<dim>*>(this),
			   level,
			   levels[level]->hexes.cells.size()-1);
}


template <int dim>
typename Triangulation<dim>::raw_line_iterator
Triangulation<dim>::last_raw_line () const
{
  if (dim == 1)
    return last_raw_line (levels.size()-1);
  else
    return last_raw_line (0);
}


template <int dim>
typename Triangulation<dim>::raw_quad_iterator
Triangulation<dim>::last_raw_quad () const
{
  if (dim == 2)
    return last_raw_quad (levels.size()-1);
  else
    return last_raw_quad (0);
}


template <int dim>
typename Triangulation<dim>::raw_hex_iterator
Triangulation<dim>::last_raw_hex () const
{
  return last_raw_hex (levels.size()-1);
}


template <int dim>
typename Triangulation<dim>::line_iterator
Triangulation<dim>::last_line (const unsigned int level) const
{
  				   // level is checked in last_raw
  raw_line_iterator ri = last_raw_line(level);
  if (ri->used()==true)
    return ri;
  while ((--ri).state() == IteratorState::valid)
    if (ri->used()==true)
      return ri;
  return ri;
}


template <int dim>
typename Triangulation<dim>::quad_iterator
Triangulation<dim>::last_quad (const unsigned int level) const
{
  				   // level is checked in last_raw
  raw_quad_iterator ri = last_raw_quad(level);
  if (ri->used()==true)
    return ri;
  while ((--ri).state() == IteratorState::valid)
    if (ri->used()==true)
      return ri;
  return ri;
}


template <int dim>
typename Triangulation<dim>::hex_iterator
Triangulation<dim>::last_hex (const unsigned int level) const
{
  				   // level is checked in last_raw
  raw_hex_iterator ri = last_raw_hex(level);
  if (ri->used()==true)
    return ri;
  while ((--ri).state() == IteratorState::valid)
    if (ri->used()==true)
      return ri;
  return ri;
}


template <int dim>
typename Triangulation<dim>::line_iterator
Triangulation<dim>::last_line () const
{
  if (dim == 1)
    return last_line (levels.size()-1);
  else
    return last_line (0);
}


template <int dim>
typename Triangulation<dim>::quad_iterator
Triangulation<dim>::last_quad () const
{
  if (dim == 2)
    return last_quad (levels.size()-1);
  else
    return last_quad (0);
}


template <int dim>
typename Triangulation<dim>::hex_iterator
Triangulation<dim>::last_hex () const
{
  return last_hex (levels.size()-1);
}


template <int dim>
typename Triangulation<dim>::active_line_iterator
Triangulation<dim>::last_active_line (const unsigned int level) const
{
				   // level is checked in last_raw
  line_iterator i = last_line(level);
  if (i->has_children()==false)
    return i;
  while ((--i).state() == IteratorState::valid)
    if (i->has_children()==false)
      return i;
  return i;
}


template <int dim>
typename Triangulation<dim>::active_quad_iterator
Triangulation<dim>::last_active_quad (const unsigned int level) const
{
				   // level is checked in last_raw
  quad_iterator i = last_quad(level);
  if (i->has_children()==false)
    return i;
  while ((--i).state() == IteratorState::valid)
    if (i->has_children()==false)
      return i;
  return i;
}


template <int dim>
typename Triangulation<dim>::active_hex_iterator
Triangulation<dim>::last_active_hex (const unsigned int level) const
{
				   // level is checked in last_raw
  hex_iterator i = last_hex(level);
  if (i->has_children()==false)
    return i;
  while ((--i).state() == IteratorState::valid)
    if (i->has_children()==false)
      return i;
  return i;
}


template <int dim>
typename Triangulation<dim>::active_line_iterator
Triangulation<dim>::last_active_line () const
{
  if (dim == 1)
    return last_active_line (levels.size()-1);
  else
    return last_active_line (0);
}


template <int dim>
typename Triangulation<dim>::active_quad_iterator
Triangulation<dim>::last_active_quad () const
{
  if (dim == 2)
    return last_active_quad (levels.size()-1);
  else
    return last_active_quad (0);
}


template <int dim>
typename Triangulation<dim>::active_hex_iterator
Triangulation<dim>::last_active_hex () const
{
  return last_active_hex (levels.size()-1);
}


template <int dim>
typename Triangulation<dim>::raw_cell_iterator
Triangulation<dim>::end_raw (const unsigned int level) const
{
  return (level == levels.size()-1 ?
	  end() :
	  begin_raw (level+1));
}


template <int dim>
typename Triangulation<dim>::cell_iterator
Triangulation<dim>::end (const unsigned int level) const
{
  return (level == levels.size()-1 ?
	  cell_iterator(end()) :
	  begin (level+1));
}


template <int dim>
typename Triangulation<dim>::active_cell_iterator
Triangulation<dim>::end_active (const unsigned int level) const
{
  return (level == levels.size()-1 ?
	  active_cell_iterator(end()) :
	  begin_active (level+1));
}


template <int dim>
typename Triangulation<dim>::raw_face_iterator
Triangulation<dim>::end_raw_face () const
{
  return end_face();
}


template <int dim>
typename Triangulation<dim>::active_face_iterator
Triangulation<dim>::end_active_face () const
{
  return active_face_iterator(end_face());
}


template <int dim>
typename Triangulation<dim>::raw_line_iterator
Triangulation<dim>::end_raw_line (const unsigned int level) const
{
  Assert (dim == 1 || level == 0, ExcFacesHaveNoLevel());
  if (dim == 1)
    return (level == levels.size()-1 ?
	    end_line() :
	    begin_raw_line (level+1));
  else
    return end_line();
}


template <int dim>
typename Triangulation<dim>::line_iterator
Triangulation<dim>::end_line (const unsigned int level) const
{
  Assert (dim == 1 || level == 0, ExcFacesHaveNoLevel());
  if (dim == 1)
    return (level == levels.size()-1 ?
	    line_iterator(end_line()) :
	    begin_line (level+1));
  else
    return line_iterator(end_line());
}


template <int dim>
typename Triangulation<dim>::active_line_iterator
Triangulation<dim>::end_active_line (const unsigned int level) const
{
  Assert (dim == 1 || level == 0, ExcFacesHaveNoLevel());
  if (dim == 1)
    return (level == levels.size()-1 ?
	    active_line_iterator(end_line()) :
	    begin_active_line (level+1));
  else
    return active_line_iterator(end_line());
}


template <int dim>
typename Triangulation<dim>::raw_quad_iterator
Triangulation<dim>::end_raw_quad (const unsigned int level) const
{
  Assert (dim == 2 || level == 0, ExcFacesHaveNoLevel());
  if (dim == 2)
    return (level == levels.size()-1 ?
	    end_quad() :
	    begin_raw_quad (level+1));
  else
    return end_quad();
}


template <int dim>
typename Triangulation<dim>::quad_iterator
Triangulation<dim>::end_quad (const unsigned int level) const
{
  Assert (dim == 2 || level == 0, ExcFacesHaveNoLevel());
  if (dim == 2)
    return (level == levels.size()-1 ?
	    quad_iterator(end_quad()) :
	    begin_quad (level+1));
  else
    return quad_iterator(end_quad());
  
}


template <int dim>
typename Triangulation<dim>::active_quad_iterator
Triangulation<dim>::end_active_quad (const unsigned int level) const
{
  Assert(dim == 2 || level == 0, ExcFacesHaveNoLevel());
  if (dim == 2)
    return (level == levels.size()-1 ?
	    active_quad_iterator(end_quad()) :
	    begin_active_quad (level+1));
  else
    return active_quad_iterator(end_quad());
}


template <int dim>
typename Triangulation<dim>::raw_hex_iterator
Triangulation<dim>::end_raw_hex (const unsigned int level) const
{
  return (level == levels.size()-1 ?
	  end_hex() :
	  begin_raw_hex (level+1));
}


template <int dim>
typename Triangulation<dim>::hex_iterator
Triangulation<dim>::end_hex (const unsigned int level) const
{
  return (level == levels.size()-1 ?
	  hex_iterator(end_hex()) :
	  begin_hex (level+1));
}


template <int dim>
typename Triangulation<dim>::active_hex_iterator
Triangulation<dim>::end_active_hex (const unsigned int level) const
{
  return (level == levels.size()-1 ?
	  active_hex_iterator(end_hex()) :
	  begin_active_hex (level+1));
}


template <int dim>
unsigned int Triangulation<dim>::n_cells () const
{
  unsigned int n=0;
  for (unsigned int l=0; l<levels.size(); ++l)
    n += n_cells (l);
  return n;
}


template <int dim>
unsigned int Triangulation<dim>::n_active_cells () const
{
  unsigned int n=0;
  for (unsigned int l=0; l<levels.size(); ++l)
    n += n_active_cells (l);
  return n;
}


template <int dim>
unsigned int Triangulation<dim>::n_faces () const
{
  Assert (dim<=3, ExcNotImplemented());
  if (dim==2) return n_lines();
  if (dim==3) return n_quads();
  return 0;
}


template <int dim>
unsigned int Triangulation<dim>::n_active_faces () const
{
  Assert (dim<=3, ExcNotImplemented());
  if (dim==2) return n_active_lines();
  if (dim==3) return n_active_quads();
  return 0;
}


//TODO:[GK] Implement this like above and remove specialization declaration in tria.h
//TODO:[GK] Add a remark to Triangulation doc telling 1d-Triangulations have no face

#if deal_II_dimension == 1

template <>
unsigned int Triangulation<1>::n_raw_cells (const unsigned int level) const
{
  return n_raw_lines (level);
}

template <>
unsigned int Triangulation<1>::n_cells (const unsigned int level) const
{
  return n_lines (level);
}


template <>
unsigned int Triangulation<1>::n_active_cells (const unsigned int level) const
{
  return n_active_lines (level);
}


#endif


#if deal_II_dimension == 2

template <>
unsigned int Triangulation<2>::n_raw_cells (const unsigned int level) const
{
  return n_raw_quads (level);
}


template <>
unsigned int Triangulation<2>::n_cells (const unsigned int level) const
{
  return n_quads (level);
}


template <>
unsigned int Triangulation<2>::n_active_cells (const unsigned int level) const
{
  return n_active_quads (level);
}

#endif


#if deal_II_dimension == 3

template <>
unsigned int Triangulation<3>::n_raw_cells (const unsigned int level) const
{
  return n_raw_hexs (level);
}


template <>
unsigned int Triangulation<3>::n_cells (const unsigned int level) const
{
  return n_hexs (level);
}


template <>
unsigned int Triangulation<3>::n_active_cells (const unsigned int level) const
{
  return n_active_hexs (level);
}

#endif


template <int dim>
unsigned int Triangulation<dim>::n_lines () const
{
  return number_cache.n_lines;
}


#if deal_II_dimension == 1

template <>
unsigned int Triangulation<1>::n_raw_lines (const unsigned int level) const
{
  Assert(level < n_levels(), ExcIndexRange(level,0,n_levels()));
  return levels[level]->lines.cells.size();
}


template <>
unsigned int Triangulation<1>::n_raw_lines () const
{
  Assert(false, ExcNotImplemented());
  return 0;
}

#endif


template <int dim>
unsigned int Triangulation<dim>::n_raw_lines (const unsigned int) const
{
  Assert(false, ExcFacesHaveNoLevel());
  return 0;
}


template <int dim>
unsigned int Triangulation<dim>::n_raw_lines () const
{
  return faces->lines.cells.size();
}


template <int dim>
unsigned int Triangulation<dim>::n_lines (const unsigned int level) const
{
  Assert (level < number_cache.n_lines_level.size(),
	  ExcIndexRange (level, 0, number_cache.n_lines_level.size()));
  Assert (dim == 1, ExcFacesHaveNoLevel());
  return number_cache.n_lines_level[level];
}


template <int dim>
unsigned int Triangulation<dim>::n_active_lines () const
{
  return number_cache.n_active_lines;
}


template <int dim>
unsigned int Triangulation<dim>::n_active_lines (const unsigned int level) const
{
  Assert (level < number_cache.n_lines_level.size(),
	  ExcIndexRange (level, 0, number_cache.n_lines_level.size()));
  Assert (dim == 1, ExcFacesHaveNoLevel());
  
  return number_cache.n_active_lines_level[level];
}


#if deal_II_dimension == 1

template <>
unsigned int Triangulation<1>::n_quads () const
{
  return 0;
}


template <>
unsigned int Triangulation<1>::n_quads (const unsigned int) const
{
  return 0;
}


template <>
unsigned int Triangulation<1>::n_raw_quads (const unsigned int) const
{
  return 0;
}


template <>
unsigned int Triangulation<1>::n_raw_hexs (const unsigned int) const
{
  return 0;
}


template <>
unsigned int Triangulation<1>::n_active_quads (const unsigned int) const
{
  return 0;
}


template <>
unsigned int Triangulation<1>::n_active_quads () const
{
  return 0;
}


#endif


template <int dim>
unsigned int Triangulation<dim>::n_quads () const
{
  return number_cache.n_quads;
}


template <int dim>
unsigned int Triangulation<dim>::n_quads (const unsigned int level) const
{
  Assert (dim == 2, ExcFacesHaveNoLevel());
  Assert (level < number_cache.n_quads_level.size(),
	  ExcIndexRange (level, 0, number_cache.n_quads_level.size()));  
  return number_cache.n_quads_level[level];
}

#if deal_II_dimension == 2


template <>
unsigned int Triangulation<2>::n_raw_quads (const unsigned int level) const
{
  Assert(level < n_levels(), ExcIndexRange(level,0,n_levels()));
  return levels[level]->quads.cells.size();
}

#endif

#if deal_II_dimension == 3

template <>
unsigned int Triangulation<3>::n_raw_quads (const unsigned int) const
{
  Assert(false, ExcFacesHaveNoLevel());
  return 0;
}

#endif




template <int dim>
unsigned int Triangulation<dim>::n_raw_quads () const
{
  Assert (false, ExcNotImplemented());
  return 0;
}


#if deal_II_dimension == 3

template <>
unsigned int Triangulation<3>::n_raw_quads () const
{
  return faces->quads.cells.size();
}

#endif

template <int dim>
unsigned int Triangulation<dim>::n_active_quads () const
{
  return number_cache.n_active_quads;
}


template <int dim>
unsigned int Triangulation<dim>::n_active_quads (const unsigned int level) const
{
  Assert (level < number_cache.n_quads_level.size(),
	  ExcIndexRange (level, 0, number_cache.n_quads_level.size()));
  Assert (dim == 2, ExcFacesHaveNoLevel());
  
  return number_cache.n_active_quads_level[level];
}


#if deal_II_dimension < 3

template <int dim>
unsigned int Triangulation<dim>::n_hexs () const
{
  return 0;
}



template <int dim>
unsigned int Triangulation<dim>::n_hexs (const unsigned int) const
{
  return 0;
}



template <int dim>
unsigned int Triangulation<dim>::n_raw_hexs (const unsigned int) const
{
  return 0;
}


template <int dim>
unsigned int Triangulation<dim>::n_active_hexs () const
{
  return 0;
}



template <int dim>
unsigned int Triangulation<dim>::n_active_hexs (const unsigned int) const
{
  return 0;
}


#else


template <int dim>
unsigned int Triangulation<dim>::n_hexs () const
{
  return number_cache.n_hexes;
}



template <int dim>
unsigned int Triangulation<dim>::n_hexs (const unsigned int level) const
{
  Assert (level < number_cache.n_hexes_level.size(),
	  ExcIndexRange (level, 0, number_cache.n_hexes_level.size()));
  
  return number_cache.n_hexes_level[level];
}



template <int dim>
unsigned int Triangulation<dim>::n_raw_hexs (const unsigned int level) const
{
  Assert(level < n_levels(), ExcIndexRange(level,0,n_levels()));
  return levels[level]->hexes.cells.size();
}


template <int dim>
unsigned int Triangulation<dim>::n_active_hexs () const
{
  return number_cache.n_active_hexes;
}



template <int dim>
unsigned int Triangulation<dim>::n_active_hexs (const unsigned int level) const
{
  Assert (level < number_cache.n_hexes_level.size(),
	  ExcIndexRange (level, 0, number_cache.n_hexes_level.size()));
  
  return number_cache.n_active_hexes_level[level];
}


#endif


template <int dim>
unsigned int Triangulation<dim>::n_levels () const
{
  if (levels.size() == 0)
    return 0;
				   // check whether there are cells on
				   // the highest levels (there need
				   // not be, since they might all
				   // have been coarsened away)
  raw_cell_iterator cell = last_raw (levels.size()-1),
		    endc = end();
  for (; cell!=endc; --cell)
    if (cell->used())
				       // return level of most refined
				       // existing cell (+1 because of
				       // counting conventions)
      return cell->level()+1;

				   // no cells at all?
  Assert (false, ExcInternalError());
				   // just to make the compiler happy:
  return 0;
}



template <int dim>
unsigned int
Triangulation<dim>::n_vertices () const 
{
  return vertices.size();
}



template <int dim>
const std::vector<Point<dim> > &
Triangulation<dim>::get_vertices () const 
{
  return vertices;
}



template <int dim>
unsigned int
Triangulation<dim>::n_used_vertices () const 
{
  return std::count_if (vertices_used.begin(), vertices_used.end(),
			std::bind2nd (std::equal_to<bool>(), true));
}



template <int dim>
const std::vector<bool> &
Triangulation<dim>::get_used_vertices () const 
{
  return vertices_used;
}


#if deal_II_dimension == 1

template <>
unsigned int Triangulation<1>::max_adjacent_cells () const
{
  return 2;
}

#endif


template <int dim>
unsigned int Triangulation<dim>::max_adjacent_cells () const
{
  cell_iterator cell = begin(0),
		endc = (levels.size() > 1 ? begin(1) : cell_iterator(end()));
				   // store the largest index of the
				   // vertices used on level 0
  unsigned int max_vertex_index = 0;
  for (; cell!=endc; ++cell)
    for (unsigned vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
      if (cell->vertex_index(vertex) > (signed int)max_vertex_index)
	max_vertex_index = cell->vertex_index(vertex);

				   // store the number of times a cell
				   // touches a vertex. An unsigned
				   // int should suffice, even for
				   // larger dimensions
  std::vector<unsigned short int> usage_count (max_vertex_index+1, 0);
				   // touch a vertex's usage count
				   // everytime we find an adjacent
				   // element
  for (cell=begin(); cell!=endc; ++cell)
    for (unsigned vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
      ++usage_count[cell->vertex_index(vertex)];

  return std::max (GeometryInfo<dim>::vertices_per_cell,
		   static_cast<unsigned int>(*std::max_element (usage_count.begin(),
								usage_count.end())));
}



template <int dim>
void
Triangulation<dim>::execute_coarsening_and_refinement ()
{
  prepare_coarsening_and_refinement ();

				   // Inform RefinementListeners
                                   // about beginning of refinement.
  typename std::list<RefinementListener *>::iterator ref_listener =
      refinement_listeners.begin (),
      end_listener = refinement_listeners.end ();
  for (; ref_listener != end_listener; ++ref_listener)
      (*ref_listener)->pre_refinement_notification (*this);

  execute_coarsening();
  execute_refinement();

				   // Inform RefinementListeners
                                   // about end of refinement.
  for (ref_listener = refinement_listeners.begin ();
       ref_listener != end_listener; ++ref_listener)
      (*ref_listener)->post_refinement_notification (*this);
}


template<int dim>
void
Triangulation<dim>::clear_despite_subscriptions()
{
				   // This is the former function
				   // clear without the assertion in
				   // the beginning.
  for (unsigned int i=0; i<levels.size(); ++i)
    delete levels[i];
  levels.clear ();
  
  delete faces;
  faces = NULL;

  vertices.clear ();
  vertices_used.clear ();
  
  for (unsigned int i=0; i<255; ++i)
    boundary[i] = straight_boundary;
  
  number_cache = internal::Triangulation::NumberCache<dim>();
}


#if deal_II_dimension == 1

template <>
void
Triangulation<1>::execute_refinement ()
{
  const unsigned int dim = 1;
  
				   // check whether a new level is
				   // needed we have to check for this
				   // on the highest level only (on
				   // this, all used cells are also
				   // active, so we only have to check
				   // for this)
  if (true)
    {
      raw_cell_iterator cell = begin_active (levels.size()-1),
			endc = end();
      for (; cell != endc; ++cell)
	if (cell->used())
	  if (cell->refine_flag_set()) 
	    {
	      levels.push_back (new internal::Triangulation::TriaLevel<dim>);
	      break;
	    }
    }


				   // check how much space is needed
				   // on every level we need not check
				   // the highest level since either -
				   // on the highest level no cells
				   // are flagged for refinement -
				   // there are, but
				   // prepare_refinement added another
				   // empty level
  unsigned int needed_vertices = 0;
  for (int level=levels.size()-2; level>=0; --level)
    {
				       // count number of flagged
				       // cells on this level
      unsigned int flagged_cells = 0;
      active_cell_iterator acell = begin_active(level),
			   aendc = begin_active(level+1);
      for (; acell!=aendc; ++acell) 
	if (acell->refine_flag_set())
	  ++flagged_cells;

				       // count number of used cells
				       // on the next higher level
      const unsigned int used_cells
	=  std::count_if (levels[level+1]->lines.used.begin(),
			  levels[level+1]->lines.used.end(),
			  std::bind2nd (std::equal_to<bool>(), true));

				       // reserve space for the
				       // used_cells cells already
				       // existing on the next higher
				       // level as well as for the
				       // 2*flagged_cells that will be
				       // created on that level
      levels[level+1]->
	internal::Triangulation::TriaLevel<0>
        ::reserve_space (used_cells+
                         GeometryInfo<1>::children_per_cell *
                         flagged_cells, 1);
				       // reserve space for
				       // 2*flagged_cells new lines on
				       // the next higher level
      levels[level+1]->lines.
	reserve_space (GeometryInfo<1>::children_per_cell*flagged_cells);
      
      needed_vertices += flagged_cells;
    }

				   // add to needed vertices how many
				   // vertices are already in use
  needed_vertices += std::count_if (vertices_used.begin(), vertices_used.end(),
				    std::bind2nd (std::equal_to<bool>(), true));
				   // if we need more vertices: create
				   // them, if not: leave the array as
				   // is, since shrinking is not
				   // really possible because some of
				   // the vertices at the end may be
				   // in use
  if (needed_vertices > vertices.size())
    {
      vertices.resize (needed_vertices, Point<1>());
      vertices_used.resize (needed_vertices, false);
    }


				   // Do REFINEMENT
				   // on every level; exclude highest
				   // level as above

				   // index of next unused vertex
  unsigned int next_unused_vertex = 0;
  
  for (int level=levels.size()-2; level>=0; --level) 
    {
      
      active_cell_iterator cell = begin_active(level),
			   endc = begin_active(level+1);
      
      raw_cell_iterator next_unused_cell = begin_raw (level+1);

      for (; (cell!=endc) && (cell->level()==level); ++cell) 
	if (cell->refine_flag_set()) 
	  {
					     // clear refinement flag
	    cell->clear_refine_flag ();

					     // search for next unused
					     // vertex
	    while (vertices_used[next_unused_vertex] == true)
	      ++next_unused_vertex;
	    Assert (next_unused_vertex < vertices.size(),
		    ExcTooFewVerticesAllocated());
	    
					     // first insert new vertex
	    Point<1> new_point = cell->vertex(0);
	    new_point += cell->vertex(1);
	    new_point /= 2.0;
	    
	    vertices[next_unused_vertex] = new_point;
	    vertices_used[next_unused_vertex] = true;

					     // search for next two
					     // unused cell (++ takes
					     // care of the end of the
					     // vector)
	    raw_cell_iterator first_child, second_child;
	    while (next_unused_cell->used() == true)
	      ++next_unused_cell;
	    first_child = next_unused_cell;
	    first_child->set_used_flag ();
	    first_child->clear_user_pointer ();
	    ++next_unused_cell;
	    Assert (next_unused_cell->used() == false,
		    ExcCellShouldBeUnused());
	    second_child = next_unused_cell;
	    second_child->set_used_flag ();
	    second_child->clear_user_pointer ();

					     // insert first child
	    cell->set_children (first_child->index());
	    first_child->clear_children ();
	    first_child->set (internal::Triangulation
                              ::Line (cell->vertex_index(0),
                                      next_unused_vertex));
	    first_child->set_material_id (cell->material_id());
	    first_child->set_subdomain_id (cell->subdomain_id());	    
	    
					     // reset neighborship info (refer
					     // to
					     // internal::Triangulation::TriaLevel<0>
					     // for details)
	    first_child->set_neighbor (1, second_child);
	    if (cell->neighbor(0).state() != IteratorState::valid)
	      first_child->set_neighbor (0, cell->neighbor(0));
	    else
	      if (cell->neighbor(0)->active())
		{
						   // since the
						   // neighbors level
						   // is always
						   // <=level, if the
						   // cell is active,
						   // then there are
						   // no cells to the
						   // left which may
						   // want to know
						   // about this new
						   // child cell.
		  Assert (cell->neighbor(0)->level() <= cell->level(),
			  ExcInternalError());
		  first_child->set_neighbor (0, cell->neighbor(0));
		}
	      else
						 // left neighbor is
						 // refined
		{
						   // set neighbor to
						   // cell on same
						   // level
		  first_child->set_neighbor (0, cell->neighbor(0)->child(1));

						   // reset neighbor
						   // info of all
						   // right descendant
						   // of the left
						   // neighbor of cell
		  cell_iterator left_neighbor = cell->neighbor(0);
		  while (left_neighbor->has_children())
		    {
		      left_neighbor = left_neighbor->child(1);
		      left_neighbor->set_neighbor (1, first_child);
		    }
		}
	    
					     // insert second child
	    second_child->clear_children ();
	    second_child->set (internal::Triangulation
                               ::Line (next_unused_vertex,
                                       cell->vertex_index(1)));
	    second_child->set_neighbor (0, first_child);
	    second_child->set_material_id (cell->material_id());
	    second_child->set_subdomain_id (cell->subdomain_id());	    
	    if (cell->neighbor(1).state() != IteratorState::valid)
	      second_child->set_neighbor (1, cell->neighbor(1));
	    else
	      if (cell->neighbor(1)->active())
		{
		  Assert (cell->neighbor(1)->level() <= cell->level(),
			  ExcInternalError());
		  second_child->set_neighbor (1, cell->neighbor(1));
		}
	      else
						 // right neighbor is
						 // refined same as
						 // above
		{
		  second_child->set_neighbor (1, cell->neighbor(1)->child(0));
		  
		  cell_iterator right_neighbor = cell->neighbor(1);
		  while (right_neighbor->has_children())
		    {
		      right_neighbor = right_neighbor->child(0);
		      right_neighbor->set_neighbor (0, second_child);
		    }
		}
	  }      
    }


				   // re-compute number of lines
  update_number_cache ();
  
#ifdef DEBUG
  for (unsigned int level=0; level<levels.size(); ++level) 
    levels[level]->lines.monitor_memory (1);

				   // check whether really all
				   // refinement flags are reset (also
				   // of previously non-active cells
				   // which we may not have
				   // touched. If the refinement flag
				   // of a non-active cell is set,
				   // something went wrong since the
				   // cell-accessors should have
				   // caught this)
  cell_iterator cell = begin(),
		endc = end();
  while (cell != endc)
    Assert (!(cell++)->refine_flag_set(), ExcInternalError ());  
#endif
}

#endif


#if deal_II_dimension == 2

template <>
void
Triangulation<2>::execute_refinement ()
{
  const unsigned int dim = 2;
  
				   // check whether a new level is
				   // needed we have to check for this
				   // on the highest level only (on
				   // this, all used cells are also
				   // active, so we only have to check
				   // for this)
  if (true)
    {
      raw_cell_iterator cell = begin_active (levels.size()-1),
			endc = end();
      for (; cell != endc; ++cell)
	if (cell->used())
	  if (cell->refine_flag_set()) 
	    {
	      levels.push_back (new internal::Triangulation::TriaLevel<dim>);
	      break;
	    }
    }


				   // first clear user flags for
				   // lines; we're going to use them
				   // to flag which lines need
				   // refinement
  for (line_iterator line=begin_line(); line!=end_line(); ++line)
    line->clear_user_flag();

				   // check how much space is needed
				   // on every level we need not check
				   // the highest level since either
				   // - on the highest level no cells
				   //   are flagged for refinement
				   // - there are, but prepare_refinement
				   //   added another empty level
  unsigned int needed_vertices = 0;
  unsigned int needed_lines  = 0;
  for (int level=levels.size()-2; level>=0; --level)
    {
      				       // count number of flagged
      				       // cells on this level and
      				       // compute how many new
      				       // vertices and new lines will
      				       // be needed
      unsigned int flagged_cells = 0;
      active_cell_iterator cell = begin_active(level),
			   endc = begin_active(level+1);
      for (; cell!=endc; ++cell) 
	if (cell->refine_flag_set()) 
	  {
	    ++flagged_cells;

					     // new vertex at center
					     // of cell is needed in
					     // any case
	    ++needed_vertices;
					     //	also the four inner
					     //	lines
	    needed_lines += 4;
	    
					     // mark all faces (lines)
					     // for refinement;
					     // checking locally
					     // whether the neighbor
					     // would also like to
					     // refine them is rather
					     // difficult for lines so
					     // we only flag them and
					     // after visiting all
					     // cells, we decide which
					     // lines need refinement;
	    for (unsigned int line_no=0; line_no<GeometryInfo<dim>::faces_per_cell;
		 ++line_no)
	      {
		line_iterator line = cell->line(line_no);
		
		if (line->has_children() == false) 
		  line->set_user_flag ();
	      }
	  }

      
      
      				       // count number of used cells
      				       // on the next higher level
      const unsigned int used_cells
	= std::count_if (levels[level+1]->quads.used.begin(),
			 levels[level+1]->quads.used.end(),
			 std::bind2nd (std::equal_to<bool>(), true));
      

				       // reserve space for the
				       // used_cells cells already
				       // existing on the next higher
				       // level as well as for the
				       // 4*flagged_cells that will be
				       // created on that level
      levels[level+1]->
	internal::Triangulation::TriaLevel<0>
        ::reserve_space (used_cells+4*flagged_cells, 2);

      				       // reserve space for
      				       // 4*flagged_cells
				       // new quads on the next higher
				       // level
      levels[level+1]->quads.
	reserve_space (4*flagged_cells);
    }

				   // now count the lines which
				   // were flagged for refinement
  for (line_iterator line=begin_line(); line!=end_line(); ++line)
    if (line->user_flag_set())
      {
	Assert (line->has_children() == false, ExcInternalError());
	needed_lines    += 2;
	needed_vertices += 1;
      }

				   // reserve space for
				   // needed_lines new lines
  faces->lines.
    reserve_space (needed_lines);
  
				   // add to needed vertices how many
				   // vertices are already in use
  needed_vertices += std::count_if (vertices_used.begin(), vertices_used.end(),
				    std::bind2nd (std::equal_to<bool>(), true));
				   // if we need more vertices: create
				   // them, if not: leave the array as
				   // is, since shrinking is not
				   // really possible because some of
				   // the vertices at the end may be
				   // in use
  if (needed_vertices > vertices.size())
    {
      vertices.resize (needed_vertices, Point<dim>());
      vertices_used.resize (needed_vertices, false);
    }


				   // Do REFINEMENT  
				   // on every level; exclude highest
				   // level as above

				   //  index of next unused vertex
  unsigned int next_unused_vertex = 0;

				   // first for lines
				   //
				   // only active objects can be
				   // refined further
  if (true)
    {
      active_line_iterator line = begin_active_line(),
			   endl = end_line();
      raw_line_iterator next_unused_line = begin_raw_line ();
      
      for (; line!=endl; ++line)
	if (line->user_flag_set())
	  {
					     // this line needs to be
					     // refined

					     // find the next unused
					     // vertex and set it
					     // appropriately
	    while (vertices_used[next_unused_vertex] == true)
	      ++next_unused_vertex;
	    Assert (next_unused_vertex < vertices.size(),
		    ExcTooFewVerticesAllocated());
	    vertices_used[next_unused_vertex] = true;

	    if (line->at_boundary())
	      vertices[next_unused_vertex]
		= boundary[line->boundary_indicator()]->get_new_point_on_line (line);
	    else
	      vertices[next_unused_vertex]
		= (line->vertex(0) + line->vertex(1)) / 2;
	  
					     // now that we created
					     // the right point, make
					     // up the two child lines
					     // (++ takes care of the
					     // end of the vector)
	    while (next_unused_line->used() == true)
	      ++next_unused_line;
					     // there should always be
					     // two consecutive unused
					     // lines, such that the
					     // children of a line
					     // will be consecutive.
					     // then set the child
					     // pointer of the present
					     // line
	    line->set_children (next_unused_line->index());
	  
					     // set the two new lines
	    raw_line_iterator children[2] = { next_unused_line,
					      ++next_unused_line };
					     // some tests; if any of
					     // the iterators should
					     // be invalid, then
					     // already dereferencing
					     // will fail
	    Assert (children[0]->used() == false, ExcCellShouldBeUnused());
	    Assert (children[1]->used() == false, ExcCellShouldBeUnused());
	  
	    children[0]->set (internal::Triangulation
                              ::Line(line->vertex_index(0),
                                     next_unused_vertex));
	    children[1]->set (internal::Triangulation
                              ::Line(next_unused_vertex,
                                     line->vertex_index(1)));
    
	    children[0]->set_used_flag();
	    children[1]->set_used_flag();
	    children[0]->clear_children();
	    children[1]->clear_children();
	    children[0]->clear_user_pointer();
	    children[1]->clear_user_pointer();
	    children[0]->clear_user_flag();
	    children[1]->clear_user_flag();

	    children[0]->set_boundary_indicator (line->boundary_indicator());
	    children[1]->set_boundary_indicator (line->boundary_indicator());
	    
					     // finally clear flag
					     // indicating the need
					     // for refinement
	    line->clear_user_flag ();
	  }
    }

				   // Now set up the new cells
  for (int level=0; level<static_cast<int>(levels.size())-1; ++level) 
    {
      
      active_cell_iterator cell = begin_active(level),
			   endc = begin_active(level+1);
      
      raw_line_iterator next_unused_line = begin_raw_line ();
      raw_cell_iterator next_unused_cell = begin_raw (level+1);

      for (; cell!=endc; ++cell) 
	if (cell->refine_flag_set()) 
	  {
					     // clear refinement flag
	    cell->clear_refine_flag ();

					     // do some additional
					     // checks.
#ifdef DEBUG
	    for (unsigned int neighbor=0;
		 neighbor<GeometryInfo<dim>::faces_per_cell; ++neighbor)
	      if (cell->neighbor(neighbor).state() == IteratorState::valid)
		Assert (((cell->neighbor(neighbor)->level() == cell->level()) &&
			 (cell->neighbor(neighbor)->coarsen_flag_set() == false))  ||
			((cell->neighbor(neighbor)->level() == cell->level()-1) &&
			 (cell->neighbor(neighbor)->refine_flag_set() == true)),
			ExcInternalError());
#endif
	    
/* For the refinement process: since we go the levels up from the lowest, there
   are (unlike above) only two possibilities: a neighbor cell is on the same
   level or one level up (in both cases, it may or may not be refined later on,
   but we don't care here).
   
   First:
   Set up an array of the 3x3 vertices, which are distributed on the cell
   (the array consists of indices into the @p{vertices} std::vector
   
   2--7--3
   |  |  |
   4--8--5
   |  |  |
   0--6--1
	
   Second:  
   Set up an array of the new lines (the array consists of iterator pointers
   into the lines arrays)
   
   .-6-.-7-.         The directions are:  .->-.->-.
   1   9   3                              ^   ^   ^
   .-10.11-.                             .->-.->-.
   0   8   2                              ^   ^   ^
   .-4-.-5-.                              .->-.->-.

     
   Third:
   Set up an array of neighbors:
   
      6  7
    .--.--.
   1|  |  |3
    .--.--.
   0|  |  |2
    .--.--.
     4   5

   We need this array for two reasons: first to get the lines which will
   bound the four subcells (if the neighboring cell is refined, these
   lines already exist), and second to update neighborship information.
   Since if a neighbor is not refined, its neighborship record only
   points to the present, unrefined, cell rather than the children we
   are presently creating, we only need the neighborship information
   if the neighbor cells are refined. In all other cases, we store
   the unrefined neighbor address

   We also need for every neighbor (if refined) which number among its
   neighbors the present (unrefined) cell has, since that number is to
   be replaced and because that also is the number of the subline which
   will be the interface between that neighbor and the to be created cell.
   We will store this number (between 0 and 3) in the field
   @p{neighbors_neighbor}.

   It would be sufficient to use the children of the common line to the
   neighbor, if we only wanted to get the new sublines and the new vertex,
   but because we need to update the neighborship information of the
   two refined subcells of the neighbor, we need to search these anyway.

   Convention:
   The created children are numbered like this:

   .--.--.
   |2 . 3|
   .--.--.
   |0 | 1|
   .--.--.
*/
	    
					     // find the next unused
					     // vertex and set it
					     // appropriately
	    while (vertices_used[next_unused_vertex] == true)
	      ++next_unused_vertex;
	    Assert (next_unused_vertex < vertices.size(),
		    ExcTooFewVerticesAllocated());
	    vertices_used[next_unused_vertex] = true;

					     // collect the indices
					     // all vertices
	    int new_vertices[9] = {cell->vertex_index(0),
				   cell->vertex_index(1),
				   cell->vertex_index(2),
				   cell->vertex_index(3),
				   cell->line(0)->child(0)->vertex_index(1),
				   cell->line(1)->child(0)->vertex_index(1),
				   cell->line(2)->child(0)->vertex_index(1),
				   cell->line(3)->child(0)->vertex_index(1),
				   next_unused_vertex};

					     // new vertex is placed
					     // at the arithmetic mean
					     // of all 8 neighboring
					     // points.
	    Point<2> new_point(0,0);
	    for (unsigned int i=0; i<8; ++i)
	      new_point +=  vertices[new_vertices[i]];
	    new_point /= 8.0;
	    
	    vertices[new_vertices[8]] = new_point;

					     // Now the lines:

					     // lines 0-7 already
					     // exist, create only the
					     // four interior lines
					     // 8-11
	    raw_line_iterator new_lines[12];
	    unsigned int l=0;
	    for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell;
		 ++face_no)
	      for (unsigned int child_no=0;
		   child_no<GeometryInfo<dim>::subfaces_per_face; ++child_no, ++l)
		new_lines[l]=cell->line(face_no)->child(child_no);
	    Assert(l==8, ExcInternalError());
		
	    for (; l<12; ++l)
	      {
		while (next_unused_line->used() == true)
		  ++next_unused_line;
		new_lines[l] = next_unused_line;
		++next_unused_line;

		Assert (new_lines[l]->used() == false,
                        ExcCellShouldBeUnused());
	      }

	    new_lines[8] ->set (internal::Triangulation
                                ::Line(new_vertices[6], new_vertices[8]));
	    new_lines[9] ->set (internal::Triangulation
                                ::Line(new_vertices[8], new_vertices[7]));
	    new_lines[10]->set (internal::Triangulation
                                ::Line(new_vertices[4], new_vertices[8]));
	    new_lines[11]->set (internal::Triangulation
                                ::Line(new_vertices[8], new_vertices[5]));
	
	    for (l=8; l<12; ++l)
	      {
		new_lines[l]->set_used_flag();
		new_lines[l]->clear_user_flag();
		new_lines[l]->clear_user_pointer();
		new_lines[l]->clear_children();
						 // interior line
		new_lines[l]->set_boundary_indicator(255);
	      }

					     // Now add the four new
					     // cells!
	    
					     // search for next unused
					     // cell the four children
					     // have to be put into
					     // the array
					     // consecutively
	    while (next_unused_cell->used() == true)
	      ++next_unused_cell;

	    const unsigned int n_children=GeometryInfo<dim>::children_per_cell;
	    raw_cell_iterator subcells[n_children];
	    for (unsigned int i=0; i<n_children; ++i) 
	      {
		Assert (next_unused_cell->used() == false,
			ExcCellShouldBeUnused());
		subcells[i] = next_unused_cell;
		++next_unused_cell;
	      }


	    cell->set_children (subcells[0]->index());

	    Assert(n_children==4, ExcNotImplemented());
	    subcells[0]->set (internal::Triangulation
                              ::Quad(new_lines[0]->index(),
                                     new_lines[8]->index(),
                                     new_lines[4]->index(),
                                     new_lines[10]->index()));
	    subcells[1]->set (internal::Triangulation
                              ::Quad(new_lines[8]->index(),
                                     new_lines[2]->index(),
                                     new_lines[5]->index(),
                                     new_lines[11]->index()));
	    subcells[2]->set (internal::Triangulation
                              ::Quad(new_lines[1]->index(),
                                     new_lines[9]->index(),
                                     new_lines[10]->index(),
                                     new_lines[6]->index()));
	    subcells[3]->set (internal::Triangulation
                              ::Quad(new_lines[9]->index(),
                                     new_lines[3]->index(),
                                     new_lines[11]->index(),
                                     new_lines[7]->index()));

	    for (unsigned int i=0; i<n_children; ++i)
	      {
		subcells[i]->set_used_flag();
		subcells[i]->clear_user_flag();
		subcells[i]->clear_user_pointer();
		subcells[i]->clear_children();
						 // inherit material
						 // properties		
		subcells[i]->set_material_id (cell->material_id());
		subcells[i]->set_subdomain_id (cell->subdomain_id());
	      }
	    
					     // now the only thing
					     // still to be done is
					     // setting neighborship
					     // information.
					     //
					     // to do so, first
					     // collect the iterators
					     // pointing to the 4x2
					     // neighbors of this
					     // cell.
					     //
					     // note that in case the
					     // neighboring cell is
					     // not refined, the
					     // neighbor iterators
					     // point to the common
					     // mother cell. the same
					     // applies if there is no
					     // neighbor: the
					     // iterators are past the
					     // end
	    cell_iterator neighbor_cells[n_children][2];
	    for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell;
		 ++face)
	      {
		const cell_iterator neighbor = cell->neighbor(face);

						 // if no neighbor
		if (neighbor.state() != IteratorState::valid)
		  for (unsigned int child_face=0;
		       child_face<GeometryInfo<dim>::subfaces_per_face;
		       ++child_face)
		    neighbor_cells[face][child_face] = neighbor;
		
		else
						   // neighbor exists
		  {
						     // neighbor's
						     // level must not
						     // be higher
						     // (else
						     // something went
						     // wrong when
						     // constructing
						     // either of the
						     // two cells) and
						     // not lower
						     // since then
						     // this cell
						     // should not
						     // have been
						     // refined.
		    Assert (neighbor->level() == cell->level(),
			    ExcInternalError());

						     // now there are
						     // two
						     // possibilities:
						     // either the
						     // neighbor has
						     // no children or
						     // it has
						     // children. these
						     // must be
						     // terminal then.
		    if (!neighbor->has_children())
		      for (unsigned int child_face=0;
			   child_face<GeometryInfo<dim>::subfaces_per_face;
			   ++child_face)
			neighbor_cells[face][child_face] = neighbor;
		    else
						       // neighbor has
						       // children
		      {
							 // first find
							 // the face
							 // of the
							 // neighbor
							 // adjacent
							 // to which
							 // the
							 // present
							 // cell is
			const unsigned int nb_nb = cell->neighbor_of_neighbor(face);
			Assert (nb_nb<GeometryInfo<dim>::faces_per_cell,
				ExcInternalError());

							 // and set
							 // the
							 // neighbors
							 // accordingly
			for (unsigned int c=0;
			     c<GeometryInfo<dim>::subfaces_per_face; ++c)
			  {
			    neighbor_cells[face][c]
			      = neighbor->child(GeometryInfo<dim>::
                                                child_cell_on_face(nb_nb, c));
			    
			    Assert (neighbor_cells[face][c].state() ==
				    IteratorState::valid,
				    ExcInternalError());
			    Assert (!neighbor_cells[face][c]->has_children(),
				    ExcInternalError());
			  }
		      }
		  }
	      }

					     // now we've got all
					     // neighbors, so set them
					     // in the new cells
	    subcells[0]->set_neighbor (0, neighbor_cells[0][0]);
	    subcells[0]->set_neighbor (1, subcells[1]);
	    subcells[0]->set_neighbor (2, neighbor_cells[2][0]);
	    subcells[0]->set_neighbor (3, subcells[2]);

	    subcells[1]->set_neighbor (0, subcells[0]);
	    subcells[1]->set_neighbor (1, neighbor_cells[1][0]);
	    subcells[1]->set_neighbor (2, neighbor_cells[2][1]);
	    subcells[1]->set_neighbor (3, subcells[3]);

	    subcells[2]->set_neighbor (0, neighbor_cells[0][1]);
	    subcells[2]->set_neighbor (1, subcells[3]);
	    subcells[2]->set_neighbor (2, subcells[0]);
	    subcells[2]->set_neighbor (3, neighbor_cells[3][0]);

	    subcells[3]->set_neighbor (0, subcells[2]);
	    subcells[3]->set_neighbor (1, neighbor_cells[1][1]);
	    subcells[3]->set_neighbor (2, subcells[1]);
	    subcells[3]->set_neighbor (3, neighbor_cells[3][1]);

					     // now we need to set the
					     // neighbors' neighborship
					     // information; this is
					     // only necessary if the
					     // neighboring cell is
					     // refined, i.e. is on
					     // the same level as the
					     // new children of the
					     // present cell
	    for (unsigned int nb=0; nb<GeometryInfo<dim>::faces_per_cell; ++nb)
	      for (unsigned int subface=0;
		   subface<GeometryInfo<dim>::subfaces_per_face; ++subface)
		if ((neighbor_cells[nb][subface].state() ==
		     IteratorState::valid) &&
		    (neighbor_cells[nb][subface]->level() ==
		     cell->level()+1))
		  {
						     // ok, the
						     // neighbor is a
						     // refined one
						     // and we need to
						     // set one of the
						     // new children
						     // as its
						     // neighbor
		    const cell_iterator neighbor = neighbor_cells[nb][subface];

						     // find which
						     // neighbor
						     // pointer is to
						     // be reset; this
						     // pointer still
						     // points to the
						     // present cell
		    unsigned int face;
		    for (face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
		      if (neighbor->neighbor(face) == cell)
			break;

		    Assert (face<GeometryInfo<dim>::faces_per_cell,
			    ExcInternalError());

		    neighbor->set_neighbor(face, subcells[
		      GeometryInfo<dim>::child_cell_on_face(nb, subface)]);
		  }

					     // note that the
					     // refinement flag was
					     // already cleared at the
					     // beginning of this loop
	    
	  }
    }

				   // re-compute number of lines and
				   // quads
  update_number_cache ();


#ifdef DEBUG
  for (unsigned int level=0; level<levels.size(); ++level) 
    levels[level]->quads.monitor_memory (2);

				   // check whether really all
				   // refinement flags are reset (also
				   // of previously non-active cells
				   // which we may not have
				   // touched. If the refinement flag
				   // of a non-active cell is set,
				   // something went wrong since the
				   // cell-accessors should have
				   // caught this)
  line_iterator line = begin_line(),
		endl = end_line();
  while (line != endl)
    Assert (!(line++)->user_flag_set(), ExcInternalError ());

  cell_iterator cell = begin(),
		endc = end();
  while (cell != endc)
    Assert (!(cell++)->refine_flag_set(), ExcInternalError ());
#endif
}

#endif


#if deal_II_dimension == 3

template <>
void
Triangulation<3>::execute_refinement ()
{
  const unsigned int dim = 3;

				   // check whether a new level is
				   // needed we have to check for this
				   // on the highest level only (on
				   // this, all used cells are also
				   // active, so we only have to check
				   // for this)
  if (true)
    {
      raw_cell_iterator cell = begin_active (levels.size()-1),
			endc = end();
      for (; cell != endc; ++cell)
	if (cell->used())
	  if (cell->refine_flag_set()) 
	    {
	      levels.push_back (new internal::Triangulation::TriaLevel<dim>);
	      break;
	    }
    }


				   // first clear user flags for quads
				   // and lines; we're going to use them
				   // to flag which lines and quads
				   // need refinement
  for (line_iterator line=begin_line(); line!=end_line(); ++line)
    line->clear_user_flag();
  for (quad_iterator quad=begin_quad(); quad!=end_quad(); ++quad)
    quad->clear_user_flag();


				   // check how much space is needed
				   // on every level
				   // we need not check the highest
				   // level since either
				   // - on the highest level no cells
				   //   are flagged for refinement
				   // - there are, but prepare_refinement
				   //   added another empty level which
				   //   then is the highest level
  unsigned int needed_vertices = 0;
  unsigned int needed_lines    = 0;
  unsigned int needed_quads    = 0;
  int level;
  for (level=levels.size()-2; level>=0; --level)
    {
      				       // count number of flagged
      				       // cells on this level and
      				       // compute how many new
      				       // vertices and new lines will
      				       // be needed
      unsigned int flagged_cells = 0;

      
      active_cell_iterator acell = begin_active(level),
			   aendc = begin_active(level+1);
      for (; acell!=aendc; ++acell) 
	if (acell->refine_flag_set()) 
	  {
	    ++flagged_cells;

					     // new vertex at center
					     // of cell is needed in
					     // any case
	    ++needed_vertices;
					     //	also the six inner
					     //	lines
	    needed_lines += 6;
					     // and the 12 inner quads
	    needed_quads += 12;
	    
					     // mark all faces and
					     // lines for refinement;
					     // checking locally
					     // whether the neighbor
					     // would also like to
					     // refine them is rather
					     // difficult for lines so
					     // we only flag them and
					     // after visiting all
					     // cells, we decide which
					     // lines need refinement;
					     // same for the quads
	    for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell;
		 ++face)
	      {
		face_iterator aface = acell->face(face);
		
		if (aface->has_children() == false) 
		  {
		    aface->set_user_flag ();
		    for (unsigned int line=0; line<4; ++line)
		      if (aface->line(line)->has_children() == false)
			aface->line(line)->set_user_flag ();
		  }
	      }
	  }

				       // count number of used cells on
				       // the next higher level
      const unsigned int used_cells
	= std::count_if (levels[level+1]->hexes.used.begin(),
			 levels[level+1]->hexes.used.end(),
			 std::bind2nd (std::equal_to<bool>(), true));


				       // reserve space for the
				       // used_cells cells already
				       // existing on the next higher
				       // level as well as for the
				       // 8*flagged_cells that will be
				       // created on that level
      levels[level+1]->
	internal::Triangulation::TriaLevel<0>
        ::reserve_space (used_cells+8*flagged_cells, 3);
      				       // reserve space for
      				       // 8*flagged_cells
				       // new hexes on the next higher
				       // level
      levels[level+1]->hexes.
	reserve_space (8*flagged_cells);
    }

				   // now count the quads and
				   // lines which were flagged for
				   // refinement
  for (quad_iterator quad=begin_quad(); quad!=end_quad(); ++quad)
    if (quad->user_flag_set())
      {
	Assert (quad->has_children() == false, ExcInternalError());
	needed_quads    += 4;
	needed_lines    += 4;
	needed_vertices += 1;
      }
  
  for (line_iterator line=begin_line(); line!=end_line(); ++line)
    if (line->user_flag_set())
      {
	Assert (line->has_children() == false, ExcInternalError());
	needed_lines += 2;
	needed_vertices += 1;
      }

				   // reserve space for
				   // needed_lines new lines
  faces->lines.
    reserve_space (needed_lines);
				   // reserve space for
				   // needed_quads new quads
  faces->quads.
    reserve_space (needed_quads);

				   // add to needed vertices how many
				   // vertices are already in use
  needed_vertices += std::count_if (vertices_used.begin(), vertices_used.end(),
				    std::bind2nd (std::equal_to<bool>(), true));
				   // if we need more vertices: create
				   // them, if not: leave the array as
				   // is, since shrinking is not
				   // really possible because some of
				   // the vertices at the end may be
				   // in use
  if (needed_vertices > vertices.size())
    {
      vertices.resize (needed_vertices, Point<dim>());
      vertices_used.resize (needed_vertices, false);
    }


				   ///////////////////////////////////////////
				   // Before we start with the actual
				   // refinement, we do some sanity
				   // checks if in debug
				   // mode. especially, we try to
				   // catch the notorious problem with
				   // lines being twice refined,
				   // i.e. there are cells adjacent at
				   // one line ("around the edge", but
				   // not at a face), with two cells
				   // differing by more than one
				   // refinement level
				   //
				   // this check is very simple to
				   // implement here, since we have
				   // all lines flagged if they shall
				   // be refined
#ifdef DEBUG
  for (active_cell_iterator cell=begin_active(); cell!=end(); ++cell)
    if (!cell->refine_flag_set())
      for (unsigned int line=0; line<GeometryInfo<dim>::lines_per_cell; ++line)
	if (cell->line(line)->has_children())
	  for (unsigned int c=0; c<2; ++c)
	    Assert (cell->line(line)->child(c)->user_flag_set() == false,
		    ExcInternalError());
#endif

				   ///////////////////////////////////////////
				   // Do refinement on every level
				   //
				   // To make life a bit easier, we
				   // first refine those lines and
				   // quads that were flagged for
				   // refinement and then compose the
				   // newly to be created cells.
				   //
				   // index of next unused vertex
  unsigned int next_unused_vertex = 0;

				   // first for lines
				   // only active objects can be
				   // refined further;
  if (true)
    {
      active_line_iterator line = begin_active_line(),
			   endl = end_line();
      raw_line_iterator next_unused_line = begin_raw_line ();
  

      for (; line!=endl; ++line)
	if (line->user_flag_set())
	  {
					     // this line needs to be
					     // refined

					     // find the next unused
					     // vertex and set it
					     // appropriately
	    while (vertices_used[next_unused_vertex] == true)
	      ++next_unused_vertex;
	    Assert (next_unused_vertex < vertices.size(),
		    ExcTooFewVerticesAllocated());
	    vertices_used[next_unused_vertex] = true;

	    if (line->at_boundary())
	      vertices[next_unused_vertex]
		= boundary[line->boundary_indicator()]->get_new_point_on_line (line);
	    else
	      vertices[next_unused_vertex]
		= (line->vertex(0) + line->vertex(1)) / 2;
	  
					     // now that we created
					     // the right point, make
					     // up the two child lines
					     // (++ takes care of the
					     // end of the vector)
	    while (next_unused_line->used() == true)
	      ++next_unused_line;
					     // there should always be
					     // two consecutive unused
					     // lines, such that the
					     // children of a line
					     // will be consecutive.
					     // then set the child
					     // pointer of the present
					     // line
	    line->set_children (next_unused_line->index());
	  
					     // set the two new lines
	    raw_line_iterator children[2] = { next_unused_line,
					      ++next_unused_line };
					     // some tests; if any of
					     // the iterators should
					     // be invalid, then
					     // already dereferencing
					     // will fail
	    Assert (children[0]->used() == false, ExcCellShouldBeUnused());
	    Assert (children[1]->used() == false, ExcCellShouldBeUnused());
	  
	    children[0]->set (internal::Triangulation
                              ::Line(line->vertex_index(0),
                                     next_unused_vertex));
	    children[1]->set (internal::Triangulation
                              ::Line(next_unused_vertex,
                                     line->vertex_index(1)));
    
	    children[0]->set_used_flag();
	    children[1]->set_used_flag();
	    children[0]->clear_children();
	    children[1]->clear_children();
	    children[0]->clear_user_pointer();
	    children[1]->clear_user_pointer();
	    children[0]->clear_user_flag();
	    children[1]->clear_user_flag();

	    children[0]->set_boundary_indicator (line->boundary_indicator());
	    children[1]->set_boundary_indicator (line->boundary_indicator());
	    
					     // finally clear flag
					     // indicating the need
					     // for refinement
	    line->clear_user_flag ();
	  }
    }


				   ///////////////////////////////////////
				   // now refine marked quads
				   ///////////////////////////////////////
				   //
				   // only active objects can be
				   // refined further;
  if (true)
    {
      active_quad_iterator quad = begin_active_quad(),
			   endq = end_quad();
      raw_line_iterator next_unused_line = begin_raw_line ();
      raw_quad_iterator next_unused_quad = begin_raw_quad ();
      
      for (; quad!=endq; ++quad)
	if (quad->user_flag_set())
	  {
					     // this quad needs to be
					     // refined

					     // find the next unused
					     // vertex and set it
					     // appropriately
	    while (vertices_used[next_unused_vertex] == true)
	      ++next_unused_vertex;
	    Assert (next_unused_vertex < vertices.size(),
		    ExcTooFewVerticesAllocated());
	    vertices_used[next_unused_vertex] = true;
	    
	    if (quad->at_boundary()) 
	      vertices[next_unused_vertex]
		= boundary[quad->boundary_indicator()]->get_new_point_on_quad (quad);
	    else
					       // it might be that the
					       // quad itself is not
					       // at the boundary, but
					       // that one of its lines
					       // actually is. in this
					       // case, the newly
					       // created vertices at
					       // the centers of the
					       // lines are not
					       // necessarily the mean
					       // values of the
					       // adjacent vertices,
					       // so do not compute
					       // the new vertex as
					       // the mean value of
					       // the 4 vertices of
					       // the face, but rather
					       // as a weighted mean
					       // value of the 8
					       // vertices which we
					       // already have (the
					       // four old ones, and
					       // the four ones
					       // inserted as middle
					       // points for the four
					       // lines). summing up
					       // some more points is
					       // generally cheaper
					       // than first asking
					       // whether one of the
					       // lines is at the
					       // boundary
					       //
					       // note that the exact
					       // weights are chosen
					       // such as to minimize
					       // the distortion of
					       // the four new quads
					       // from the optimal
					       // shape; their
					       // derivation and
					       // values is copied
					       // over from the
					       // @p{MappingQ::set_laplace_on_vector}
					       // function
	      vertices[next_unused_vertex]
		= (quad->vertex(0) + quad->vertex(1) +
		   quad->vertex(2) + quad->vertex(3) +
		   3*(quad->line(0)->child(0)->vertex(1) +
		      quad->line(1)->child(0)->vertex(1) +
		      quad->line(2)->child(0)->vertex(1) +
		      quad->line(3)->child(0)->vertex(1))   ) / 16;
	  
					     // now that we created
					     // the right point, make
					     // up the four lines
					     // interior to the quad
					     // (++ takes care of the
					     // end of the vector)
	    raw_line_iterator new_lines[4];

	    for (unsigned int i=0; i<4; ++i)
	      {
		while (next_unused_line->used() == true)
		  ++next_unused_line;
		new_lines[i] = next_unused_line;
		++next_unused_line;

		Assert (new_lines[i]->used() == false,
                        ExcCellShouldBeUnused());
	      }

					     // set the data of the
					     // four lines.
					     // first collect the
					     // indices of the five
					     // vertices:
					     // *--3--*
					     // |  |  |
					     // 0--4--1
					     // |  |  |
					     // *--2--*
					     // the lines are numbered
					     // as follows:
					     // *--*--*
					     // |  1  |
					     // *2-*-3*
					     // |  0  |
					     // *--*--*
	    const unsigned int vertex_indices[5]
	      = { quad->line(0)->child(0)->vertex_index(1),
		  quad->line(1)->child(0)->vertex_index(1),
		  quad->line(2)->child(0)->vertex_index(1),
		  quad->line(3)->child(0)->vertex_index(1),
		  next_unused_vertex 
	      };
	    
	    new_lines[0]->set (internal::Triangulation
                               ::Line(vertex_indices[2], vertex_indices[4]));
	    new_lines[1]->set (internal::Triangulation
                               ::Line(vertex_indices[4], vertex_indices[3]));
	    new_lines[2]->set (internal::Triangulation
                               ::Line(vertex_indices[0], vertex_indices[4]));
	    new_lines[3]->set (internal::Triangulation
                               ::Line(vertex_indices[4], vertex_indices[1]));

	    for (unsigned int i=0; i<4; ++i)
	      {
		new_lines[i]->set_used_flag();
		new_lines[i]->clear_user_flag();
		new_lines[i]->clear_user_pointer();
		new_lines[i]->clear_children();
		new_lines[i]->set_boundary_indicator(quad->boundary_indicator());
	      }


					     // now for the
					     // quads. again, first
					     // collect some data
					     // about the indices of
					     // the lines, with the
					     // following numbering:
					     //   .-6-.-7-.
					     //   1   9   3
					     //   .-10.11-.
					     //   0   8   2
					     //   .-4-.-5-.

					     // child 0 and 1 of a line are
					     // switched if the line orientation
					     // is false. set up a miniature
					     // table, indicating which child to
					     // take for line orientations false
					     // and true. first index: child
					     // index in standard orientation,
					     // second index: line orientation
	    unsigned int index[2][2]=
	      {{1,0},   // child 0, line_orientation=false and true
	       {0,1}};  // child 1, line_orientation=false and true
		    
	    const unsigned int line_indices[12]
	      = { quad->line(0)->child(index[0][quad->line_orientation(0)])->index(),
		  quad->line(0)->child(index[1][quad->line_orientation(0)])->index(),
		  quad->line(1)->child(index[0][quad->line_orientation(1)])->index(),
		  quad->line(1)->child(index[1][quad->line_orientation(1)])->index(),
		  quad->line(2)->child(index[0][quad->line_orientation(2)])->index(),
		  quad->line(2)->child(index[1][quad->line_orientation(2)])->index(),
		  quad->line(3)->child(index[0][quad->line_orientation(3)])->index(),
		  quad->line(3)->child(index[1][quad->line_orientation(3)])->index(),
		  new_lines[0]->index(),
		  new_lines[1]->index(),
		  new_lines[2]->index(),
		  new_lines[3]->index() 
	      };
	    
					     // find some space for
					     // the four newly to be
					     // created quads.  note
					     // that there should
					     // always be four
					     // consecutive free slots
					     // for them
	    raw_quad_iterator new_quads[4];

	    while (next_unused_quad->used() == true)
	      ++next_unused_quad;

	    new_quads[0] = next_unused_quad;
	    Assert (new_quads[0]->used() == false, ExcCellShouldBeUnused());

	    ++next_unused_quad;
	    new_quads[1] = next_unused_quad;
	    Assert (new_quads[1]->used() == false, ExcCellShouldBeUnused());

	    ++next_unused_quad;
	    new_quads[2] = next_unused_quad;
	    Assert (new_quads[2]->used() == false, ExcCellShouldBeUnused());

	    ++next_unused_quad;
	    new_quads[3] = next_unused_quad;
	    Assert (new_quads[3]->used() == false, ExcCellShouldBeUnused());

					     // note these quads as
					     // children to the
					     // present one
	    quad->set_children (new_quads[0]->index());

	    new_quads[0]->set (internal::Triangulation
                               ::Quad(line_indices[0],
                                      line_indices[8],
                                      line_indices[4],
                                      line_indices[10]));
	    new_quads[1]->set (internal::Triangulation
                               ::Quad(line_indices[8],
                                      line_indices[2],
                                      line_indices[5],
                                      line_indices[11]));
	    new_quads[2]->set (internal::Triangulation
                               ::Quad(line_indices[1],
                                      line_indices[9],
                                      line_indices[10],
                                      line_indices[6]));
	    new_quads[3]->set (internal::Triangulation
                               ::Quad(line_indices[9],
                                      line_indices[3],
                                      line_indices[11],
                                      line_indices[7]));
	    for (unsigned int i=0; i<4; ++i)
	      {
		new_quads[i]->set_used_flag();
		new_quads[i]->clear_user_flag();
		new_quads[i]->clear_user_pointer();
		new_quads[i]->clear_children();
		new_quads[i]->set_boundary_indicator (quad->boundary_indicator());
						 // set all line orientations to
						 // true, change this after the
						 // loop, as we have to consider
						 // different lines for each
						 // child
		for (unsigned int j=0; j<GeometryInfo<dim>::lines_per_face; ++j)
		  new_quads[i]->set_line_orientation(j,true);
	      }
					     // now set the line orientation of
					     // children of outer lines
					     // correctly, the lines in the
					     // interior of the refined quad are
					     // automatically oriented
					     // conforming to the standard
	    new_quads[0]->set_line_orientation(0,quad->line_orientation(0));
	    new_quads[0]->set_line_orientation(2,quad->line_orientation(2));
	    new_quads[1]->set_line_orientation(1,quad->line_orientation(1));
	    new_quads[1]->set_line_orientation(2,quad->line_orientation(2));
	    new_quads[2]->set_line_orientation(0,quad->line_orientation(0));
	    new_quads[2]->set_line_orientation(3,quad->line_orientation(3));
	    new_quads[3]->set_line_orientation(1,quad->line_orientation(1));
	    new_quads[3]->set_line_orientation(3,quad->line_orientation(3));
				       
	    
					     // finally clear flag
					     // indicating the need
					     // for refinement
	    quad->clear_user_flag ();
	  }
    }

				   ///////////////////////////////////
				   // Now, finally, set up the new
				   // cells
				   ///////////////////////////////////
  for (unsigned int level=0; level!=levels.size()-1; ++level)
    {
				       // only active objects can be
				       // refined further; remember
				       // that we won't operate on the
				       // finest level, so
				       // begin_*(level+1) is allowed
      active_hex_iterator hex  = begin_active_hex(level),
			  endh = begin_active_hex(level+1);
      raw_line_iterator next_unused_line = begin_raw_line ();
      raw_quad_iterator next_unused_quad = begin_raw_quad ();
      raw_hex_iterator  next_unused_hex  = begin_raw_hex (level+1);

      for (; hex!=endh; ++hex)
	if (hex->refine_flag_set())
	  {
					     // do some additional
					     // checks.
#ifdef DEBUG
	    for (unsigned int neighbor=0;
		 neighbor<GeometryInfo<dim>::faces_per_cell; ++neighbor)
	      if (hex->neighbor(neighbor).state() == IteratorState::valid)
		Assert (((hex->neighbor(neighbor)->level() == hex->level()) &&
			 (hex->neighbor(neighbor)->coarsen_flag_set() == false))  ||
			((hex->neighbor(neighbor)->level() == hex->level()-1) &&
			 (hex->neighbor(neighbor)->refine_flag_set() == true)),
			ExcInternalError());
#endif
					     // this hex needs to be
					     // refined
	    
					     // clear flag indicating
					     // the need for
					     // refinement. do it here
					     // already, since we
					     // can't do it anymore
					     // once the cell has
					     // children
	    hex->clear_refine_flag ();

					     // find the next unused
					     // vertex and set it
					     // appropriately
	    while (vertices_used[next_unused_vertex] == true)
	      ++next_unused_vertex;
	    Assert (next_unused_vertex < vertices.size(),
		    ExcTooFewVerticesAllocated());
	    vertices_used[next_unused_vertex] = true;
	    
					     // the new vertex is
					     // definitely in the
					     // interior, so we need
					     // not worry about the
					     // boundary.  let it be
					     // the average of the 26
					     // vertices surrounding
					     // it. weight these
					     // vertices in the same
					     // way as they are
					     // weighted in the
					     // @p{MappingQ::set_laplace_on_hex_vector}
					     // function, and like the
					     // new vertex at the
					     // center of the quad is
					     // weighted (see above)
	    vertices[next_unused_vertex] = Point<dim>();
					     // first add corners of hex
	    for (unsigned int vertex=0;
		 vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
	      vertices[next_unused_vertex] += hex->vertex(vertex) / 128;
					     // now add center of lines
	    for (unsigned int line=0;
		 line<GeometryInfo<dim>::lines_per_cell; ++line)
	      vertices[next_unused_vertex] += hex->line(line)->child(0)->vertex(1) *
					      7./192.;
					     // finally add centers of
					     // faces. note that
					     // vertex 3 of child 0 is an
					     // invariant with respect
					     // to the face
					     // orientation, flip and rotation
	    for (unsigned int face=0;
		 face<GeometryInfo<dim>::faces_per_cell; ++face)
	      vertices[next_unused_vertex] += hex->face(face)->child(0)->vertex(3) *
					      1./12.;

					     // now that we created
					     // the right point, make
					     // up the six lines
					     // interior to the hex
					     // (++ takes care of the
					     // end of the vector)
	    raw_line_iterator new_lines[6];

	    for (unsigned int i=0; i<6; ++i)
	      {
		while (next_unused_line->used() == true)
		  ++next_unused_line;
		new_lines[i] = next_unused_line;
		++next_unused_line;

		Assert (new_lines[i]->used() == false,
                        ExcCellShouldBeUnused());
	      }

					     // set the data of the
					     // six lines.  first
					     // collect the indices of
					     // the seven vertices
					     // (consider the two
					     // planes to be crossed
					     // to form the planes
					     // cutting the hex in two
					     // vertically and
					     // horizontally)
					     //     *--3--*   *--5--*
					     //    /  /  /    |  |  |
					     //   0--6--1     0--6--1
					     //  /  /  /      |  |  |
					     // *--2--*       *--4--*
					     // the lines are numbered
					     // as follows:
					     //     *--*--*   *--*--*
					     //    /  1  /    |  5  |
					     //   *2-*-3*     *2-*-3*
					     //  /  0  /      |  4  |
					     // *--*--*       *--*--*
                                             //
                                             // note that asking for child 0 and
                                             // vertex 3 within that is
                                             // invariant with respect to the
                                             // face orientation, rotation and
                                             // flip, so we do not have to ask
                                             // here
	    const unsigned int vertex_indices[7]
	      = { hex->face(0)->child(0)->vertex_index(3),
		  hex->face(1)->child(0)->vertex_index(3),
		  hex->face(2)->child(0)->vertex_index(3),
		  hex->face(3)->child(0)->vertex_index(3),
		  hex->face(4)->child(0)->vertex_index(3),
		  hex->face(5)->child(0)->vertex_index(3),
		  next_unused_vertex 
	      };
	    
	    new_lines[0]->set (internal::Triangulation
                               ::Line(vertex_indices[2], vertex_indices[6]));
	    new_lines[1]->set (internal::Triangulation
                               ::Line(vertex_indices[6], vertex_indices[3]));
	    new_lines[2]->set (internal::Triangulation
                               ::Line(vertex_indices[0], vertex_indices[6]));
	    new_lines[3]->set (internal::Triangulation
                               ::Line(vertex_indices[6], vertex_indices[1]));
	    new_lines[4]->set (internal::Triangulation
                               ::Line(vertex_indices[4], vertex_indices[6]));
	    new_lines[5]->set (internal::Triangulation
                               ::Line(vertex_indices[6], vertex_indices[5]));

	    for (unsigned int i=0; i<6; ++i)
	      {
		new_lines[i]->set_used_flag();
		new_lines[i]->clear_user_flag();
		new_lines[i]->clear_user_pointer();
		new_lines[i]->clear_children();
						 // interior line
		new_lines[i]->set_boundary_indicator(255);
	      }


					     // now for the
					     // quads. again, first
					     // collect some data
					     // about the indices of
					     // the lines, with the
					     // following numbering:
					     // (note that face 0 and
					     // 1 each are shown twice
					     // for better
					     // readability)
	    
					     // face 0: left plane
					     //       *            *
					     //      /|           /| 
					     //     * |          * |
					     //    /| *         /| *
					     //   * 1/|        * |3|
					     //   | * |        | * |
					     //   |/| *        |2| *
					     //   * 0/         * |/
					     //   | *          | *
					     //   |/           |/ 
					     //   *            *
					     // face 1: right plane
					     //       *            *
					     //      /|           /| 
					     //     * |          * |
					     //    /| *         /| *
					     //   * 5/|        * |7|
					     //   | * |        | * |
					     //   |/| *        |6| *
					     //   * 4/         * |/
					     //   | *          | *
					     //   |/           |/
					     //   *            *
					     // face 2: front plane
					     //   (note: x,y exchanged)
					     //   *---*---*
					     //   |   11  |
					     //   *-8-*-9-*
					     //   |   10  |
					     //   *---*---*
					     // face 3: back plane
					     //   (note: x,y exchanged)
					     //   *---*---*
					     //   |   15  | 
					     //   *12-*-13*
					     //   |   14  | 
					     //   *---*---*
					     // face 4: bottom plane
					     //       *---*---*
					     //      /  17   /
					     //     *18-*-19*
					     //    /   16  /
					     //   *---*---*
					     // face 5: top plane
					     //       *---*---*
					     //      /  21   /
					     //     *22-*-23*
					     //    /   20  /
					     //   *---*---*
					     // middle planes
					     //     *---*---*   *---*---*
					     //    /  25   /    |   29  |
					     //   *26-*-27*     *26-*-27*
					     //  /   24  /      |   28  |
					     // *---*---*       *---*---*
                                             //
                                             // this time we have to
                                             // take into account
                                             // whether the different
                                             // faces are oriented
                                             // correctly or in the
                                             // opposite direction, so
                                             // store that up front

					     // TODO: shorten this
	    
					     // face_orientation
            const bool f_or[6]
              = { hex->face_orientation (0),
                  hex->face_orientation (1),
                  hex->face_orientation (2),
                  hex->face_orientation (3),
                  hex->face_orientation (4),
                  hex->face_orientation (5) };
                    
					     // face_flip
            const bool f_fl[6]
              = { hex->face_flip (0),
                  hex->face_flip (1),
                  hex->face_flip (2),
                  hex->face_flip (3),
                  hex->face_flip (4),
                  hex->face_flip (5) };
                    
					     // face_rotation
            const bool f_ro[6]
              = { hex->face_rotation (0),
                  hex->face_rotation (1),
                  hex->face_rotation (2),
                  hex->face_rotation (3),
                  hex->face_rotation (4),
                  hex->face_rotation (5) };

					     // set up a list of line iterators
					     // first. from this, construct
					     // lists of line_indices and
					     // line orientations later on
	    const raw_line_iterator lines[30]
	      = {
		    hex->face(0)->child(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[0],f_fl[0],f_ro[0]))
                    ->line(GeometryInfo<dim>::standard_to_real_face_line(1,f_or[0],f_fl[0],f_ro[0])),        //0
		    hex->face(0)->child(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[0],f_fl[0],f_ro[0]))
                    ->line(GeometryInfo<dim>::standard_to_real_face_line(0,f_or[0],f_fl[0],f_ro[0])),        //1
		    hex->face(0)->child(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[0],f_fl[0],f_ro[0]))
                    ->line(GeometryInfo<dim>::standard_to_real_face_line(3,f_or[0],f_fl[0],f_ro[0])),        //2
		    hex->face(0)->child(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[0],f_fl[0],f_ro[0]))
                    ->line(GeometryInfo<dim>::standard_to_real_face_line(2,f_or[0],f_fl[0],f_ro[0])),        //3

		    hex->face(1)->child(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[1],f_fl[1],f_ro[1]))
                    ->line(GeometryInfo<dim>::standard_to_real_face_line(1,f_or[1],f_fl[1],f_ro[1])),        //4
		    hex->face(1)->child(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[1],f_fl[1],f_ro[1]))
                    ->line(GeometryInfo<dim>::standard_to_real_face_line(0,f_or[1],f_fl[1],f_ro[1])),        //5
		    hex->face(1)->child(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[1],f_fl[1],f_ro[1]))
                    ->line(GeometryInfo<dim>::standard_to_real_face_line(3,f_or[1],f_fl[1],f_ro[1])),        //6
		    hex->face(1)->child(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[1],f_fl[1],f_ro[1]))
                    ->line(GeometryInfo<dim>::standard_to_real_face_line(2,f_or[1],f_fl[1],f_ro[1])),        //7

		    hex->face(2)->child(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[2],f_fl[2],f_ro[2]))
                    ->line(GeometryInfo<dim>::standard_to_real_face_line(1,f_or[2],f_fl[2],f_ro[2])),        //8
		    hex->face(2)->child(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[2],f_fl[2],f_ro[2]))
                    ->line(GeometryInfo<dim>::standard_to_real_face_line(0,f_or[2],f_fl[2],f_ro[2])),        //9
		    hex->face(2)->child(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[2],f_fl[2],f_ro[2]))
                    ->line(GeometryInfo<dim>::standard_to_real_face_line(3,f_or[2],f_fl[2],f_ro[2])),        //10
		    hex->face(2)->child(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[2],f_fl[2],f_ro[2]))
                    ->line(GeometryInfo<dim>::standard_to_real_face_line(2,f_or[2],f_fl[2],f_ro[2])),        //11

		    hex->face(3)->child(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[3],f_fl[3],f_ro[3]))
                    ->line(GeometryInfo<dim>::standard_to_real_face_line(1,f_or[3],f_fl[3],f_ro[3])),        //12
		    hex->face(3)->child(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[3],f_fl[3],f_ro[3]))
                    ->line(GeometryInfo<dim>::standard_to_real_face_line(0,f_or[3],f_fl[3],f_ro[3])),        //13
		    hex->face(3)->child(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[3],f_fl[3],f_ro[3]))
                    ->line(GeometryInfo<dim>::standard_to_real_face_line(3,f_or[3],f_fl[3],f_ro[3])),        //14
		    hex->face(3)->child(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[3],f_fl[3],f_ro[3]))
                    ->line(GeometryInfo<dim>::standard_to_real_face_line(2,f_or[3],f_fl[3],f_ro[3])),        //15

		    hex->face(4)->child(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[4],f_fl[4],f_ro[4]))
                    ->line(GeometryInfo<dim>::standard_to_real_face_line(1,f_or[4],f_fl[4],f_ro[4])),        //16
		    hex->face(4)->child(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[4],f_fl[4],f_ro[4]))
                    ->line(GeometryInfo<dim>::standard_to_real_face_line(0,f_or[4],f_fl[4],f_ro[4])),        //17
		    hex->face(4)->child(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[4],f_fl[4],f_ro[4]))
                    ->line(GeometryInfo<dim>::standard_to_real_face_line(3,f_or[4],f_fl[4],f_ro[4])),        //18
		    hex->face(4)->child(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[4],f_fl[4],f_ro[4]))
                    ->line(GeometryInfo<dim>::standard_to_real_face_line(2,f_or[4],f_fl[4],f_ro[4])),        //19

		    hex->face(5)->child(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[5],f_fl[5],f_ro[5]))
                    ->line(GeometryInfo<dim>::standard_to_real_face_line(1,f_or[5],f_fl[5],f_ro[5])),        //20
		    hex->face(5)->child(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[5],f_fl[5],f_ro[5]))
                    ->line(GeometryInfo<dim>::standard_to_real_face_line(0,f_or[5],f_fl[5],f_ro[5])),        //21
		    hex->face(5)->child(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[5],f_fl[5],f_ro[5]))
                    ->line(GeometryInfo<dim>::standard_to_real_face_line(3,f_or[5],f_fl[5],f_ro[5])),        //22
		    hex->face(5)->child(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[5],f_fl[5],f_ro[5]))
                    ->line(GeometryInfo<dim>::standard_to_real_face_line(2,f_or[5],f_fl[5],f_ro[5])),        //23

		    new_lines[0],                       //24
		    new_lines[1],                       //25
		    new_lines[2],                       //26
		    new_lines[3],                       //27
		    new_lines[4],                       //28
		    new_lines[5]                        //29
	      };

	    unsigned int line_indices[30];
	    for (unsigned int i=0; i<30; ++i)
	      line_indices[i]=lines[i]->index();
	    
					     // the orientation of lines for the
					     // inner quads is quite tricky. as
					     // these lines are newly created
					     // ones and thus have no parents,
					     // they cannot inherit this
					     // property. set up an array and
					     // fill it with the respective
					     // values
	    bool line_orientation[30];
					     // note: for the first 24 lines
					     // (inner lines of the outer quads)
					     // the following holds: the second
					     // vertex of the even lines in
					     // standard orientation is the
					     // vertex in the middle of the
					     // quad, whereas for odd lines the
					     // first vertex is the same middle
					     // vertex.
	    for (unsigned int i=0; i<24; ++i)
	      if (lines[i]->vertex_index((i+1)%2)==static_cast<int>(vertex_indices[i/4]))
		line_orientation[i]=true;
	      else
		{
						   // it must be the other way
						   // round then
		  Assert(lines[i]->vertex_index(i%2)==static_cast<int>(vertex_indices[i/4]),
			 ExcInternalError());
		  line_orientation[i]=false;
		}
					     // for the last 6 lines the line
					     // orientation is always true,
					     // since they were just constructed
					     // that way
	    for (unsigned int i=24; i<30; ++i)
	      line_orientation[i]=true;

	    
					     // find some space for
					     // the 12 newly to be
					     // created quads.
	    raw_quad_iterator new_quads[12];

	    for (unsigned int i=0; i<12; ++i)
	      {
		while (next_unused_quad->used() == true)
		  ++next_unused_quad;
		new_quads[i] = next_unused_quad;
		++next_unused_quad;

		Assert (new_quads[i]->used() == false,
                        ExcCellShouldBeUnused());
	      }

					     // set up the 12 quads,
					     // numbered as follows
					     // (left quad numbering,
					     // right line numbering
					     // extracted from above)
					     //
					     //      *          *
					     //     /|        21|
					     //    * |        * 15
					     //  y/|3*      20| *  
					     //  * |/|      * |/| 
					     //  |2* |x    11 * 14
					     //  |/|1*      |/| * 
					     //  * |/       * |17 
					     //  |0*       10 *
					     //  |/         |16
					     //  *          *
					     //
					     //  x
					     //  *---*---*      *22-*-23*
					     //  | 5 | 7 |      1  29   5
					     //  *---*---*      *26-*-27*
					     //  | 4 | 6 |      0  28   4
					     //  *---*---*y     *18-*-19*
					     //
					     //       y
					     //      *----*----*      *-12-*-13-*
					     //     / 10 / 11 /      3   25    7
					     //    *----*----*      *-26-*-27-*
					     //   / 8  / 9  /      2   24    6
					     //  *----*----*x     *--8-*--9-*
	    
	    new_quads[0]->set (internal::Triangulation
                               ::Quad(line_indices[10],
                                      line_indices[28],
                                      line_indices[16],
                                      line_indices[24]));
	    new_quads[1]->set (internal::Triangulation
                               ::Quad(line_indices[28],
                                      line_indices[14],
                                      line_indices[17],
                                      line_indices[25]));
	    new_quads[2]->set (internal::Triangulation
                               ::Quad(line_indices[11],
                                      line_indices[29],
                                      line_indices[24],
                                      line_indices[20]));
	    new_quads[3]->set (internal::Triangulation
                               ::Quad(line_indices[29],
                                      line_indices[15],
                                      line_indices[25],
                                      line_indices[21]));
	    new_quads[4]->set (internal::Triangulation
                               ::Quad(line_indices[18],
                                      line_indices[26],
                                      line_indices[0],
                                      line_indices[28]));
	    new_quads[5]->set (internal::Triangulation
                               ::Quad(line_indices[26],
                                      line_indices[22],
                                      line_indices[1],
                                      line_indices[29]));
	    new_quads[6]->set (internal::Triangulation
                               ::Quad(line_indices[19],
                                      line_indices[27],
                                      line_indices[28],
                                      line_indices[4]));
	    new_quads[7]->set (internal::Triangulation
                               ::Quad(line_indices[27],
                                      line_indices[23],
                                      line_indices[29],
                                      line_indices[5]));
	    new_quads[8]->set (internal::Triangulation
                               ::Quad(line_indices[2],
                                      line_indices[24],
                                      line_indices[8],
                                      line_indices[26]));
	    new_quads[9]->set (internal::Triangulation
                               ::Quad(line_indices[24],
                                      line_indices[6],
                                      line_indices[9],
                                      line_indices[27]));
	    new_quads[10]->set (internal::Triangulation
                                ::Quad(line_indices[3],
                                       line_indices[25],
                                       line_indices[26],
                                       line_indices[12]));
	    new_quads[11]->set (internal::Triangulation
                                ::Quad(line_indices[25],
                                       line_indices[7],
                                       line_indices[27],
                                       line_indices[13]));

					     // set flags
	    for (unsigned int i=0; i<12; ++i)
	      {
		new_quads[i]->set_used_flag();
		new_quads[i]->clear_user_flag();
		new_quads[i]->clear_user_pointer();
		new_quads[i]->clear_children();
						 // interior quad
		new_quads[i]->set_boundary_indicator (255);
						 // set all line orientation
						 // flags to true, chnage this
						 // afterwards, if necessary
		for (unsigned int j=0; j<GeometryInfo<dim>::lines_per_face; ++j)
		  new_quads[i]->set_line_orientation(j,true);
	      }  
					     // now reset the line_orientation
					     // flags of outer lines as they
					     // cannot be set in a loop (at
					     // least not easily)
	    new_quads[0]->set_line_orientation(0,line_orientation[10]);
	    new_quads[0]->set_line_orientation(2,line_orientation[16]);

	    new_quads[1]->set_line_orientation(1,line_orientation[14]);
	    new_quads[1]->set_line_orientation(2,line_orientation[17]);

	    new_quads[2]->set_line_orientation(0,line_orientation[11]);
	    new_quads[2]->set_line_orientation(3,line_orientation[20]);

	    new_quads[3]->set_line_orientation(1,line_orientation[15]);
	    new_quads[3]->set_line_orientation(3,line_orientation[21]);

	    new_quads[4]->set_line_orientation(0,line_orientation[18]);
	    new_quads[4]->set_line_orientation(2,line_orientation[0]);

	    new_quads[5]->set_line_orientation(1,line_orientation[22]);
	    new_quads[5]->set_line_orientation(2,line_orientation[1]);

	    new_quads[6]->set_line_orientation(0,line_orientation[19]);
	    new_quads[6]->set_line_orientation(3,line_orientation[4]);

	    new_quads[7]->set_line_orientation(1,line_orientation[23]);
	    new_quads[7]->set_line_orientation(3,line_orientation[5]);

	    new_quads[8]->set_line_orientation(0,line_orientation[2]);
	    new_quads[8]->set_line_orientation(2,line_orientation[8]);

	    new_quads[9]->set_line_orientation(1,line_orientation[6]);
	    new_quads[9]->set_line_orientation(2,line_orientation[9]);

	    new_quads[10]->set_line_orientation(0,line_orientation[3]);
	    new_quads[10]->set_line_orientation(3,line_orientation[12]);

	    new_quads[11]->set_line_orientation(1,line_orientation[7]);
	    new_quads[11]->set_line_orientation(3,line_orientation[13]);


					     /////////////////////////////////
					     // create the eight new hexes
					     //
					     // again first collect
					     // some data.  here, we
					     // need the indices of a
					     // whole lotta
					     // quads. they are
					     // numbered as follows:
					     //
					     // planes in the interior
					     // of the old hex:
					     //      *
					     //     /|
					     //    * |
					     //   /|3*  *---*---*      *----*----*
					     //  * |/|  | 5 | 7 |     / 10 / 11 /
					     //  |2* |  *---*---*    *----*----*
					     //  |/|1*  | 4 | 6 |   / 8  / 9  /
					     //  * |/   *---*---*y *----*----*x
					     //  |0*
					     //  |/
					     //  *
					     //
					     // children of the faces
					     // of the old hex
					     //      *-------*        *-------*
					     //     /|25   27|       /34   35/|
					     //    15|       |      /       /19
					     //   /  |       |     /32   33/  |
					     //  *   |24   26|    *-------*18 |
					     //  1413*-------*    |21   23| 17*
					     //  |  /30   31/     |       |  /
					     //  12/       /      |       |16
					     //  |/28   29/       |20   22|/
					     //  *-------*        *-------*
                                             //
                                             // note that we have to
                                             // take care of the
                                             // orientation of
                                             // faces.
					     //
					     // TODO: simplify this
	    const unsigned int quad_indices[36]
	      = {
		    new_quads[0]->index(),     //0
		    new_quads[1]->index(),
		    new_quads[2]->index(),
		    new_quads[3]->index(),
		    new_quads[4]->index(),
		    new_quads[5]->index(),
		    new_quads[6]->index(),
		    new_quads[7]->index(),
		    new_quads[8]->index(),
		    new_quads[9]->index(),
		    new_quads[10]->index(),
		    new_quads[11]->index(),    //11

		    hex->face(0)->child_index(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[0],f_fl[0],f_ro[0])),  //12
		    hex->face(0)->child_index(GeometryInfo<dim>::standard_to_real_face_vertex(1,f_or[0],f_fl[0],f_ro[0])),
		    hex->face(0)->child_index(GeometryInfo<dim>::standard_to_real_face_vertex(2,f_or[0],f_fl[0],f_ro[0])),
		    hex->face(0)->child_index(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[0],f_fl[0],f_ro[0])),

		    hex->face(1)->child_index(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[1],f_fl[1],f_ro[1])),  //16
		    hex->face(1)->child_index(GeometryInfo<dim>::standard_to_real_face_vertex(1,f_or[1],f_fl[1],f_ro[1])),
		    hex->face(1)->child_index(GeometryInfo<dim>::standard_to_real_face_vertex(2,f_or[1],f_fl[1],f_ro[1])),
		    hex->face(1)->child_index(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[1],f_fl[1],f_ro[1])),

		    hex->face(2)->child_index(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[2],f_fl[2],f_ro[2])),  //20
		    hex->face(2)->child_index(GeometryInfo<dim>::standard_to_real_face_vertex(1,f_or[2],f_fl[2],f_ro[2])),
		    hex->face(2)->child_index(GeometryInfo<dim>::standard_to_real_face_vertex(2,f_or[2],f_fl[2],f_ro[2])),
		    hex->face(2)->child_index(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[2],f_fl[2],f_ro[2])),

		    hex->face(3)->child_index(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[3],f_fl[3],f_ro[3])),  //24
		    hex->face(3)->child_index(GeometryInfo<dim>::standard_to_real_face_vertex(1,f_or[3],f_fl[3],f_ro[3])),
		    hex->face(3)->child_index(GeometryInfo<dim>::standard_to_real_face_vertex(2,f_or[3],f_fl[3],f_ro[3])),
		    hex->face(3)->child_index(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[3],f_fl[3],f_ro[3])),

		    hex->face(4)->child_index(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[4],f_fl[4],f_ro[4])),  //28
		    hex->face(4)->child_index(GeometryInfo<dim>::standard_to_real_face_vertex(1,f_or[4],f_fl[4],f_ro[4])),
		    hex->face(4)->child_index(GeometryInfo<dim>::standard_to_real_face_vertex(2,f_or[4],f_fl[4],f_ro[4])),
		    hex->face(4)->child_index(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[4],f_fl[4],f_ro[4])),

		    hex->face(5)->child_index(GeometryInfo<dim>::standard_to_real_face_vertex(0,f_or[5],f_fl[5],f_ro[5])),  //32
		    hex->face(5)->child_index(GeometryInfo<dim>::standard_to_real_face_vertex(1,f_or[5],f_fl[5],f_ro[5])),
		    hex->face(5)->child_index(GeometryInfo<dim>::standard_to_real_face_vertex(2,f_or[5],f_fl[5],f_ro[5])),
		    hex->face(5)->child_index(GeometryInfo<dim>::standard_to_real_face_vertex(3,f_or[5],f_fl[5],f_ro[5]))
	      };


					     // find some space for
					     // the eight newly to be
					     // created hexes.  note
					     // that there should
					     // always be eight
					     // consecutive free slots
					     // for them
	    raw_hex_iterator new_hexes[8];

	    while (next_unused_hex->used() == true)
	      ++next_unused_hex;

	    for (unsigned int i=0; i<8; ++i)
	      {
		new_hexes[i] = next_unused_hex;
		Assert (new_hexes[i]->used() == false,
                        ExcCellShouldBeUnused());
		++next_unused_hex;
	      }

					     // note these hexes as
					     // children to the
					     // present cell
	    hex->set_children (new_hexes[0]->index());

					     // bottom children
	    new_hexes[0]->set (internal::Triangulation
                               ::Hexahedron(quad_indices[12],
                                            quad_indices[0],
                                            quad_indices[20],
                                            quad_indices[4],
                                            quad_indices[28],
                                            quad_indices[8]));
	    new_hexes[1]->set (internal::Triangulation
                               ::Hexahedron(quad_indices[0],
                                            quad_indices[16],
                                            quad_indices[22],
                                            quad_indices[6],
                                            quad_indices[29],
                                            quad_indices[9]));
	    new_hexes[2]->set (internal::Triangulation
                               ::Hexahedron(quad_indices[13],
                                            quad_indices[1],
                                            quad_indices[4],
                                            quad_indices[24],
                                            quad_indices[30],
                                            quad_indices[10]));
	    new_hexes[3]->set (internal::Triangulation
                               ::Hexahedron(quad_indices[1],
                                            quad_indices[17],
                                            quad_indices[6],
                                            quad_indices[26],
                                            quad_indices[31],
                                            quad_indices[11]));
	    
					     // top children
	    new_hexes[4]->set (internal::Triangulation
                               ::Hexahedron(quad_indices[14],
                                            quad_indices[2],
                                            quad_indices[21],
                                            quad_indices[5],
                                            quad_indices[8],
                                            quad_indices[32]));
	    new_hexes[5]->set (internal::Triangulation
                               ::Hexahedron(quad_indices[2],
                                            quad_indices[18],
                                            quad_indices[23],
                                            quad_indices[7],
                                            quad_indices[9],
                                            quad_indices[33]));
	    new_hexes[6]->set (internal::Triangulation
                               ::Hexahedron(quad_indices[15],
                                            quad_indices[3],
                                            quad_indices[5],
                                            quad_indices[25],
                                            quad_indices[10],
                                            quad_indices[34]));
	    new_hexes[7]->set (internal::Triangulation
                               ::Hexahedron(quad_indices[3],
                                            quad_indices[19],
                                            quad_indices[7],
                                            quad_indices[27],
                                            quad_indices[11],
                                            quad_indices[35]));
            

	    for (unsigned int i=0; i<8; ++i)
	      {
		new_hexes[i]->set_used_flag();
		new_hexes[i]->clear_user_flag();
		new_hexes[i]->clear_user_pointer();
		new_hexes[i]->clear_children();
						 // inherit material
						 // properties
		new_hexes[i]->set_material_id (hex->material_id());
		new_hexes[i]->set_subdomain_id (hex->subdomain_id());
              }
            
                                             // and set face
                                             // orientation
                                             // flags. note that
                                             // new faces in the
                                             // interior of the
                                             // mother cell always
                                             // have a correctly
                                             // oriented face, but
                                             // the ones on the
                                             // outer faces will
                                             // inherit this flag
                                             //
                                             // set the flag to true
                                             // for all faces
                                             // initially, then go the
                                             // other way round and
                                             // reset faces that are
                                             // at the boundary of the
                                             // mother cube
					     //
					     // the same is true for the
					     // face_flip and face_rotation
					     // flags. however, the latter two
					     // are set to false by default as
					     // this is the standard value
	    for (unsigned int i=0; i<8; ++i)
              for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
		{
		  new_hexes[i]->set_face_orientation(f, true);
		  new_hexes[i]->set_face_flip(f, false);
		  new_hexes[i]->set_face_rotation(f, false);
		}
	    
            for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
              for (unsigned int s=0; s<GeometryInfo<dim>::subfaces_per_face; ++s)
                {
		  new_hexes[GeometryInfo<dim>::child_cell_on_face(f,s)]
		    ->set_face_orientation(f, hex->face_orientation(f));
		  new_hexes[GeometryInfo<dim>::child_cell_on_face(f,s)]
		    ->set_face_flip(f, hex->face_flip(f));
		  new_hexes[GeometryInfo<dim>::child_cell_on_face(f,s)]
		    ->set_face_rotation(f, hex->face_rotation(f));
		}
            
#ifdef DEBUG
					     // check consistency
					     // against
					     // GeometryInfo<3>::child_cell_on_face
	    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
	      for (unsigned int s=0; s<GeometryInfo<dim>::subfaces_per_face; ++s) 
		Assert(hex->face(f)->child(s)==hex->child(
		  GeometryInfo<dim>::child_cell_on_face(
		    f, s, hex->face_orientation(f), hex->face_flip(f), hex->face_rotation(f)))->face(f), ExcInternalError());
#endif

					     /////////////////////////////////
					     // now the only thing still
					     // to be done is setting
					     // neighborship
					     // information.
					     //
					     // to do so, first
					     // collect the iterators
					     // pointing to the 6x4
					     // neighbors of this
					     // cell.
					     //
					     // note that in case the
					     // neighboring cell is
					     // not refined, the
					     // neighbor iterators
					     // point to the common
					     // mother cell. the same
					     // applies if there is no
					     // neighbor: the
					     // iterators are past the
					     // end
	    cell_iterator neighbor_cells[6][4];
	    for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell;
		 ++face)
	      {
		const cell_iterator neighbor = hex->neighbor(face);

						 // if no neighbor
		if (neighbor.state() != IteratorState::valid)
		  for (unsigned int child_face=0;
		       child_face<GeometryInfo<dim>::subfaces_per_face;
		       ++child_face)
		    neighbor_cells[face][child_face] = neighbor;
		
		else
						   // neighbor exists
		  {
						     // neighbor's
						     // level must not
						     // be higher
						     // (else
						     // something went
						     // wrong when
						     // constructing
						     // either of the
						     // two cells) and
						     // not lower
						     // since then
						     // this cell
						     // should not
						     // have been
						     // refined.
		    Assert (neighbor->level() == hex->level(),
			    ExcInternalError());

						     // now there are
						     // two
						     // possibilities:
						     // either the
						     // neighbor has
						     // no children or
						     // it has
						     // children. these
						     // must be
						     // terminal then.
		    if (!neighbor->has_children())
		      for (unsigned int child_face=0;
			   child_face<GeometryInfo<dim>::subfaces_per_face;
			   ++child_face)
			neighbor_cells[face][child_face] = neighbor;
		    else
						       // neighbor has
						       // children;
						       // now it's
						       // getting
						       // complicated
		      {
							 // first find
							 // the face
							 // of the
							 // neighbor
							 // adjacent
							 // to which
							 // the
							 // present
							 // cell is
			const unsigned int nb_nb = hex->neighbor_of_neighbor(face);
			Assert (nb_nb<GeometryInfo<dim>::faces_per_cell,
				ExcInternalError());

							 // now the
							 // four child
							 // cells of
							 // neighbor
							 // adjacent
							 // to the
							 // present
							 // cell can
							 // be
							 // obtained
							 // by a
							 // function
							 // of
							 // GeometryInfo. however,
							 // if our face or the
							 // neighbor's
							 // face has
							 // non-standard
							 // orientation, flip or rotation,
							 // then we
							 // run into
							 // trouble
							 // and have
							 // to swap
							 // subfaces
							 // to account
							 // for
							 // that.
                                                         //
							 // it is quite
							 // difficult to find
							 // out, which neighbors
							 // child is adjacent to
							 // a given child of our
							 // current cell. it is
							 // easier to determine
							 // the children
							 // adjacent to a given
							 // subface in both
							 // cases. here, we
							 // consider that for
							 // the neighbor, when
							 // we set the neighbors
							 // of our children we
							 // will consider it for
							 // the children of our
							 // current cell.

			for (unsigned int c=0;
			     c<GeometryInfo<dim>::subfaces_per_face; ++c)
			  {
			    neighbor_cells[face][c]
			      = neighbor->child(GeometryInfo<dim>::
						child_cell_on_face(nb_nb, c,
								   neighbor->face_orientation(nb_nb),
								   neighbor->face_flip(nb_nb),
								   neighbor->face_rotation(nb_nb)));
			    
			    Assert (neighbor_cells[face][c].state() ==
				    IteratorState::valid,
				    ExcInternalError());
			    Assert (!neighbor_cells[face][c]->has_children(),
				    ExcInternalError());
			  }
		      }
		  }
	      }

					     // now we've got all
					     // neighbors, so set them
					     // in the new cells
	    
					     // the subface ordering
					     // is as follows
					     // face 0 and 1:
					     //        * 
					     //       /| 
					     //      * | 
					     //    y/|3* 
					     //    * |/| 
					     //    |2* |x
					     //    |/|1* 
					     //    * |/  
					     //    |0*   
					     //    |/    
					     //    *
					     // face 2 and 3:
					     //  x
					     //  *---*---*
					     //  | 1 | 3 |
					     //  *---*---*
					     //  | 0 | 2 |
					     //  *---*---*y
					     // face 4 and 5:
					     //       y
					     //      *---*---*
					     //     / 2 / 3 /
					     //    *---*---*
					     //   / 0 / 1 /
					     //  *---*---*x

					     // Step 1: set all 'outer neighbors'
	    for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	      for (unsigned int sub=0; sub<GeometryInfo<dim>::subfaces_per_face; ++sub)
		new_hexes[GeometryInfo<dim>::
			  child_cell_on_face(face,sub,
					     hex->face_orientation(face),
					     hex->face_flip(face),
					     hex->face_rotation(face))]
		  ->set_neighbor(face,neighbor_cells[face][sub]);
					     // Step 2: set all 'inner
					     // neighbors', i.e. neighbors that
					     // are itself children of our
					     // current cell

					     // TODO: do we really want to test
					     // (assert) this static information
					     // for each cell?
	    new_hexes[0]->set_neighbor (1, new_hexes[1]);
	    new_hexes[0]->set_neighbor (3, new_hexes[2]);
	    new_hexes[0]->set_neighbor (5, new_hexes[4]);
	    Assert(GeometryInfo<dim>::child_cell_on_face(0,0)==0, ExcInternalError());
	    Assert(GeometryInfo<dim>::child_cell_on_face(2,0)==0, ExcInternalError());
	    Assert(GeometryInfo<dim>::child_cell_on_face(4,0)==0, ExcInternalError());

	    new_hexes[1]->set_neighbor (0, new_hexes[0]);
	    new_hexes[1]->set_neighbor (3, new_hexes[3]);
	    new_hexes[1]->set_neighbor (5, new_hexes[5]);
	    Assert(GeometryInfo<dim>::child_cell_on_face(1,0)==1, ExcInternalError());
	    Assert(GeometryInfo<dim>::child_cell_on_face(2,2)==1, ExcInternalError());
	    Assert(GeometryInfo<dim>::child_cell_on_face(4,1)==1, ExcInternalError());

	    new_hexes[2]->set_neighbor (1, new_hexes[3]);
	    new_hexes[2]->set_neighbor (2, new_hexes[0]);
	    new_hexes[2]->set_neighbor (5, new_hexes[6]);
	    Assert(GeometryInfo<dim>::child_cell_on_face(0,1)==2, ExcInternalError());
	    Assert(GeometryInfo<dim>::child_cell_on_face(3,0)==2, ExcInternalError());
	    Assert(GeometryInfo<dim>::child_cell_on_face(4,2)==2, ExcInternalError());

	    new_hexes[3]->set_neighbor (0, new_hexes[2]);
	    new_hexes[3]->set_neighbor (2, new_hexes[1]);
	    new_hexes[3]->set_neighbor (5, new_hexes[7]);
	    Assert(GeometryInfo<dim>::child_cell_on_face(1,1)==3, ExcInternalError());
	    Assert(GeometryInfo<dim>::child_cell_on_face(3,2)==3, ExcInternalError());
	    Assert(GeometryInfo<dim>::child_cell_on_face(4,3)==3, ExcInternalError());

	    new_hexes[4]->set_neighbor (1, new_hexes[5]);
	    new_hexes[4]->set_neighbor (3, new_hexes[6]);
	    new_hexes[4]->set_neighbor (4, new_hexes[0]);
	    Assert(GeometryInfo<dim>::child_cell_on_face(0,2)==4, ExcInternalError());
	    Assert(GeometryInfo<dim>::child_cell_on_face(2,1)==4, ExcInternalError());
	    Assert(GeometryInfo<dim>::child_cell_on_face(5,0)==4, ExcInternalError());

	    new_hexes[5]->set_neighbor (0, new_hexes[4]);
	    new_hexes[5]->set_neighbor (3, new_hexes[7]);
	    new_hexes[5]->set_neighbor (4, new_hexes[1]);
	    Assert(GeometryInfo<dim>::child_cell_on_face(1,2)==5, ExcInternalError());
	    Assert(GeometryInfo<dim>::child_cell_on_face(2,3)==5, ExcInternalError());
	    Assert(GeometryInfo<dim>::child_cell_on_face(5,1)==5, ExcInternalError());

	    new_hexes[6]->set_neighbor (1, new_hexes[7]);
	    new_hexes[6]->set_neighbor (2, new_hexes[4]);
	    new_hexes[6]->set_neighbor (4, new_hexes[2]);
	    Assert(GeometryInfo<dim>::child_cell_on_face(0,3)==6, ExcInternalError());
	    Assert(GeometryInfo<dim>::child_cell_on_face(3,1)==6, ExcInternalError());
	    Assert(GeometryInfo<dim>::child_cell_on_face(5,2)==6, ExcInternalError());

	    new_hexes[7]->set_neighbor (0, new_hexes[6]);
	    new_hexes[7]->set_neighbor (2, new_hexes[5]);
	    new_hexes[7]->set_neighbor (4, new_hexes[3]);
	    Assert(GeometryInfo<dim>::child_cell_on_face(1,3)==7, ExcInternalError());
	    Assert(GeometryInfo<dim>::child_cell_on_face(3,3)==7, ExcInternalError());
	    Assert(GeometryInfo<dim>::child_cell_on_face(5,3)==7, ExcInternalError());


					     // now we need to set the
					     // neighbors' neighborship
					     // information; this is
					     // only necessary if the
					     // neighboring cell is
					     // refined, i.e. is on
					     // the same level as the
					     // new children of the
					     // present cell
	    for (unsigned int nb=0; nb<GeometryInfo<dim>::faces_per_cell; ++nb)
	      for (unsigned int subface=0;
		   subface<GeometryInfo<dim>::subfaces_per_face; ++subface)
		if ((neighbor_cells[nb][subface].state() ==
		     IteratorState::valid) &&
		    (neighbor_cells[nb][subface]->level() ==
		     hex->level()+1))
		  {
						     // ok, the
						     // neighbor is a
						     // refined one
						     // and we need to
						     // set one of the
						     // new children
						     // as its
						     // neighbor
		    const cell_iterator neighbor = neighbor_cells[nb][subface];

						     // find which
						     // neighbor
						     // pointer is to
						     // be reset; this
						     // pointer still
						     // points to the
						     // present cell
		    unsigned int face;
		    for (face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
		      if (neighbor->neighbor(face) == hex)
			break;

		    Assert (face<GeometryInfo<dim>::faces_per_cell,
			    ExcInternalError());
						     // check, that the face we
						     // found is consistent with
						     // the information obtained
						     // by the
						     // neighbor_of_neighbor()
						     // function
		    Assert (face==hex->neighbor_of_neighbor(nb),
			    ExcInternalError());

                                                     // then figure out which of the new
                                                     // cells points to this neighbor.
						     //
                                                     // We have considered
                                                     // face-orientation, -flip
                                                     // and -rotation issues in
                                                     // the construction of the
                                                     // neighbor_cells array,
                                                     // now we have to consider
                                                     // it for this cell as
                                                     // well.
		    int c = GeometryInfo<dim>::child_cell_on_face(nb,subface,
								  hex->face_orientation(nb),
								  hex->face_flip(nb),
								  hex->face_rotation(nb));
		    neighbor->set_neighbor(face,new_hexes[c]);
		  }


					     // note that the
					     // refinement flag was
					     // already cleared at the
					     // beginning of this loop
	  }
    }

				   // re-compute number of lines and
				   // quads
  update_number_cache ();


#ifdef DEBUG
  for (unsigned int level=0; level<levels.size(); ++level) 
    levels[level]->hexes.monitor_memory (3);

				   // check whether really all
				   // refinement flags are reset (also
				   // of previously non-active cells
				   // which we may not have
				   // touched. If the refinement flag
				   // of a non-active cell is set,
				   // something went wrong since the
				   // cell-accessors should have
				   // caught this)
  line_iterator line = begin_line(),
		endl = end_line();
  while (line != endl)
    Assert (!(line++)->user_flag_set(), ExcInternalError ());

  quad_iterator quad = begin_quad(),
		endq = end_quad();
  while (quad != endq)
    Assert (!(quad++)->user_flag_set(), ExcInternalError ());

  cell_iterator cell = begin(),
		endc = end();
  while (cell != endc)
    Assert (!(cell++)->refine_flag_set(), ExcInternalError ());
#endif
}


#endif


template <int dim>
void Triangulation<dim>::execute_coarsening ()
{
  				   // loop over all cells. Flag all
  				   // cells of which all children are
  				   // flagged for
				   // coarsening and delete the childrens'
				   // flags. In effect, only those
				   // cells are flagged of which originally
				   // all children were flagged and for which
				   // all children are on the same refinement
				   // level. For flagging, the user flags are
				   // used, to avoid confusion and because
				   // non-active cells can't be flagged for
				   // coarsening. Note that because of the
				   // effects of @p{prepare_coarsening}, of a
				   // cell either all or no children must
				   // be flagged for coarsening, so it is
				   // ok to only check the first child
  clear_user_flags ();

  cell_iterator cell = begin(),
		endc = end();
  for (; cell!=endc; ++cell) 
    if (!cell->active())
      if (cell->child(0)->coarsen_flag_set())
	{
	  cell->set_user_flag();
	  for (unsigned int child=0;
	       child<GeometryInfo<dim>::children_per_cell; ++child)
	    {
	      Assert (cell->child(child)->coarsen_flag_set(),
		      ExcInternalError());
	      cell->child(child)->clear_coarsen_flag();
	    }
	}


				   // now do the actual coarsening
				   // step. Since the loop goes over
				   // used cells we only need not
				   // worry about deleting some cells
				   // since the ++operator will then
				   // just hop over them if we should
				   // hit one. Do the loop in the
				   // reverse way since we may only
				   // delete some cells if their
				   // neighbors have already been
				   // deleted (if the latter are on a
				   // higher level for example)
				   //
				   // if there is only one level,
				   // there can not be anything to do
  if (levels.size() >= 2)
    for (cell = last(levels.size()-2); cell!=endc; --cell)
      if (cell->user_flag_set())
					 // use a separate function,
					 // since this is dimension
					 // specific
	delete_children (cell);

  				   // re-compute number of lines and
  				   // quads
  update_number_cache ();

  				   // in principle no user flags
  				   // should be
				   // set any more at this point
#if DEBUG
  for (cell=begin(); cell!=endc; ++cell)
    Assert (cell->user_flag_set() == false, ExcInternalError());
#endif
}


template <int dim>
void Triangulation<dim>::prepare_refinement_dim_dependent () 
{}


#if deal_II_dimension == 3

template <>
void Triangulation<3>::prepare_refinement_dim_dependent () 
{
  const unsigned int dim = 3;
  
				   // first clear flags on lines,
				   // since we need them to determine
				   // which lines will be refined
  clear_user_flags_line();

				   // variable to store whether the
				   // mesh was changed in the present
				   // loop and in the whole process
  bool mesh_changed      = false;

  do
    {
      mesh_changed = false;

				       // for this following, we need to know
				       // which cells are going to be
				       // coarsened, if we had to make a
				       // decision. the following function
				       // sets these flags:
      fix_coarsen_flags ();

      
				       // flag those lines that will
				       // be refined
      for (active_cell_iterator cell=begin_active(); cell!=end(); ++cell)
	if (cell->refine_flag_set())
	  for (unsigned int line=0; line<GeometryInfo<dim>::lines_per_cell; ++line)
					     // if the line is not yet
					     // refined, it will be in
					     // the process
	    if (!cell->line(line)->has_children())
	      cell->line(line)->set_user_flag();


				       // now check whether there are
				       // cells with lines that are
				       // more than once refined or
				       // that will be more than once
				       // refined. The first thing
				       // should never be the case, in
				       // the second case we flag the
				       // cell for refinement
      for (active_cell_iterator cell=last_active(); cell!=end(); --cell)
	for (unsigned int line=0; line<GeometryInfo<dim>::lines_per_cell; ++line)
	  {
	    if (cell->line(line)->has_children())
	      {
						 // if this line is
						 // refined, its
						 // children should
						 // not have further
						 // children
						 //
						 // however, if any of
						 // the children is
						 // flagged for
						 // further
						 // refinement, we
						 // need to refine
						 // this cell also (at
						 // least, if the cell
						 // is not already
						 // flagged)
		bool offending_line_found = false;
		
		for (unsigned int c=0; c<2; ++c)
		  {
		    Assert (cell->line(line)->child(c)->has_children() == false,
			    ExcInternalError());
		    
		    if (cell->line(line)->child(c)->user_flag_set () &&
			!cell->refine_flag_set())
		      {
							 // tag this
							 // cell for
							 // refinement
			cell->clear_coarsen_flag ();
			cell->set_refine_flag();
			
							 // note that
							 // we have
							 // changed
							 // the grid
			offending_line_found = true;
			
							 // it may save us several
							 // loop iterations if we
							 // flag all lines of
							 // this cell now (and not
							 // at the outset of the
							 // next iteration) for
							 // refinement
			for (unsigned int line=0;
			     line<GeometryInfo<dim>::lines_per_cell; ++line)
			  if (!cell->line(line)->has_children())
			    cell->line(line)->set_user_flag();
			
			break;
		      }
		  }

		if (offending_line_found)
		  {
		    mesh_changed = true;
		    break;
		  }
	      }
	  }
	

				       // there is another thing here:
				       // if any of the lines will be
				       // refined, then we may not
				       // coarsen the present cell
      for (active_cell_iterator cell=last_active(); cell!=end(); --cell)
	{
	  if (cell->coarsen_flag_set())
	    for (unsigned int line=0; line<GeometryInfo<dim>::lines_per_cell; ++line)
	      {
		if (cell->line(line)->user_flag_set())
		  {
		    cell->clear_coarsen_flag ();
		    mesh_changed = true;
		
		    goto next_cell;
		  }

						 // similarly, if any of the lines
						 // *is* already refined, we may
						 // not coarsen the current
						 // cell. however, there's a
						 // catch: if the line is refined,
						 // but the cell behind it is
						 // going to be coarsened, then
						 // the situation changes. if we
						 // forget this second condition,
						 // the refine_and_coarsen_3d test
						 // will start to fail. note that
						 // to know which cells are going
						 // to be coarsened, the call for
						 // fix_coarsen_flags above is
						 // necessary
						 //
						 // the problem is that finding
						 // all cells that are behind an
						 // edge in 3d is somewhat of a
						 // pain and worst of all
						 // introduces a quadratic
						 // behavior in this algorithm. on
						 // the other hand, not many cells
						 // have their coarsen flag set
						 // usually, and fixing
						 // refine_and_coarsen_3d is a
						 // somewhat important case
		if (cell->line(line)->has_children())
		  {
		    bool cant_be_coarsened = false;

						     // loop over all cells of this
						     // level to find neighbors of
						     // this cell and edge
		    for (cell_iterator edge_neighbor=begin(cell->level());
			 ((edge_neighbor != end(cell->level()))
			  &&
			  (cant_be_coarsened == false));
			 ++edge_neighbor)
		      if (edge_neighbor != cell)
			for (unsigned int e=0; e<GeometryInfo<dim>::lines_per_cell; ++e)
			  if (edge_neighbor->line(e) == cell->line(line))
			    {
							       // this is a cell
							       // that is adjacent
							       // to the present
							       // cell across this
							       // edge. so treat
							       // it, but only if
							       // it is actually
							       // refined or will
							       // be refined
			      if (! (edge_neighbor->has_children()
				     ||
				     (!edge_neighbor->has_children() &&
				      edge_neighbor->refine_flag_set())))
				break;

							       // figure out if
							       // the neighbor is
							       // going to be
							       // coarsened. as a
							       // post-condition
							       // of the call to
							       // fix_coarsen_flags(),
							       // either all
							       // children of a
							       // cell must be
							       // flagged for
							       // coarsening, or
							       // none may. above
							       // we delete some
							       // coarsen flags,
							       // and in the next
							       // call to
							       // fix_coarsen_flags()
							       // the flags to all
							       // siblings will be
							       // removed. we will
							       // check here if
							       // still all
							       // children have
							       // that flag set
			      unsigned int n_children_flagged = 0;
			      for (unsigned int c=0; c<GeometryInfo<dim>::children_per_cell; ++c)
				if ((edge_neighbor->child(c)->has_children() == false)
				    &&
				    edge_neighbor->child(c)->coarsen_flag_set())
				  ++n_children_flagged;
			  
							       // now, if not all
							       // children are
							       // flagged, then
							       // the neighboring
							       // cell isn't going
							       // to be
							       // coarsened. that
							       // means that the
							       // common edge
							       // isn't going to
							       // be coarsened and
							       // that we can't
							       // coarsen the
							       // present cell
			      if (n_children_flagged !=
				  GeometryInfo<dim>::children_per_cell)
				cant_be_coarsened = true;


							       // neighbor was
							       // found. no reason
							       // to keep looping
							       // over edges of
							       // the possible
							       // edge_neighbor
			      break;
			    }

		    if (cant_be_coarsened == true)
		      {
			cell->clear_coarsen_flag ();
			mesh_changed = true;

			goto next_cell;
		      }
		  }
	      }
	  next_cell:
	  ;
	}
    }
  while (mesh_changed == true);
}

#endif


template <int dim>
void Triangulation<dim>::fix_coarsen_flags ()
{
				   // loop over all cells. Flag all
				   // cells of which all children are
				   // flagged for coarsening and
				   // delete the childrens'
				   // flags. Also delete all flags of
				   // cells for which not all children
				   // of a cell are flagged. In
				   // effect, only those cells are
				   // flagged of which originally all
				   // children were flagged and for
				   // which all children are on the
				   // same refinement level. For
				   // flagging, the user flags are
				   // used, to avoid confusion and
				   // because non-active cells can't
				   // be flagged for coarsening
				   //
				   // In effect, all coarsen flags are
				   // turned into user flags of the
				   // mother cell if coarsening is
				   // possible or deleted
				   // otherwise. Coarsen flags of
				   // cells with no mother cell,
				   // i.e. on the coarsest level are
				   // deleted explicitly.
  clear_user_flags ();
      
  cell_iterator cell = begin(),
		endc = end();
  for (; cell!=endc; ++cell) 
    {
				       // nothing to do if we are
				       // already on the finest level;
				       // if we are on the coarsest
				       // level, delete coarsen flag
				       // since no coarsening possible
      if (cell->active()) 
	{
	  if (cell->level() == 0)
	    cell->clear_coarsen_flag();
	  continue;
	}
	  
      unsigned int flagged_children = 0;
      for (unsigned int child=0; child<GeometryInfo<dim>::children_per_cell; ++child)
	if (cell->child(child)->active() &&
	    cell->child(child)->coarsen_flag_set()) 
	  {
	    ++flagged_children;
					     // clear flag since we
					     // don't need it anymore
	    cell->child(child)->clear_coarsen_flag();
	  }
	  
				       // flag this cell for
				       // coarsening if all children
				       // were flagged
      if (flagged_children == GeometryInfo<dim>::children_per_cell)
	cell->set_user_flag();
    }
      
				   // in principle no coarsen flags
				   // should be set any more at this
				   // point
#if DEBUG
  for (cell=begin(); cell!=endc; ++cell)
    Assert (cell->coarsen_flag_set() == false, ExcInternalError());
#endif

				   // now loop over all cells which have the
				   // user flag set. their children were
				   // flagged for coarsening. set the coarsen
				   // flag again if we are sure that none of
				   // the neighbors of these children are
				   // refined, or will be refined, since then
				   // we would get a two-level jump in
				   // refinement. on the other hand, if one of
				   // the children's neighbors has their user
				   // flag set, then we know that its children
				   // will go away by coarsening, and we will
				   // be ok.
				   //
				   // note on the other hand that we do allow
				   // level-2 jumps in refinement between
				   // neighbors in 1d, so this whole procedure
				   // is only necessary if we are not in 1d
				   //
				   // since we remove some coarsening/user
				   // flags in the process, we have to work
				   // from the finest level to the coarsest
				   // one, since we occasionally inspect user
				   // flags of cells on finer levels and need
				   // to be sure that these flags are final
  for (cell=last(); cell!=endc; --cell)
    if (cell->user_flag_set())
      {
	bool coarsening_allowed = true;

	if (dim > 1)
	  {
	    for (unsigned int c=0;
		 (c<GeometryInfo<dim>::children_per_cell) && (coarsening_allowed==true);
		 ++c)
	      for (unsigned int n=0; n<GeometryInfo<dim>::faces_per_cell; ++n)
		{
		  const cell_iterator child_neighbor = cell->child(c)->neighbor(n);
		  if ((child_neighbor.state() == IteratorState::valid)
		      &&
		      (child_neighbor->level()==cell->level()+1)
		      &&
		      ((child_neighbor->has_children()
			&&
			!child_neighbor->user_flag_set())
		       ||
		       (child_neighbor->has_children()
			&&
			child_neighbor->refine_flag_set())))
		    {
		      coarsening_allowed = false;
		      break;
		    }
		}
	  }
	
					 // if allowed: tag the
					 // children for coarsening
	if (coarsening_allowed)
	  for (unsigned int c=0; c<GeometryInfo<dim>::children_per_cell; ++c)
	    {
	      Assert (cell->child(c)->refine_flag_set()==false,
		      ExcInternalError());
	      
	      cell->child(c)->set_coarsen_flag();
	    }
      }

				   // clear all user flags again, now that we
				   // don't need them any more
  clear_user_flags ();
}



template <int dim>
bool Triangulation<dim>::prepare_coarsening_and_refinement ()
{
      
				   // save the flags to determine
				   // whether something was changed in
				   // the course of this function
  std::vector<bool> flags_before[2];
  save_coarsen_flags (flags_before[0]);
  save_refine_flags (flags_before[1]);

				   // do nothing in 1d, except setting
				   // the coarsening flags correctly
  if (dim == 1)
    {
      fix_coarsen_flags ();

      std::vector<bool> flags_after[2];
      save_coarsen_flags (flags_after[0]);
      save_refine_flags (flags_after[1]);

      return ((flags_before[0] != flags_after[0]) ||
	      (flags_before[1] != flags_after[1]));
    }

				   // for all other dimensions

				   // save the flags at the outset of
				   // each loop. we do so in order to
				   // find out whether something was
				   // changed in the present loop, in
				   // which case we would have to
				   // re-run the loop. the other
				   // possibility to find this out
				   // would be to set a flag
				   // @p{something_changed} to true
				   // each time we change something.
				   // however, sometimes one change in
				   // one of the parts of the loop is
				   // undone by another one, so we
				   // might end up in an endless
				   // loop. we could be tempted to
				   // break this loop at an arbitrary
				   // number of runs, but that would
				   // not be a clean solution, since
				   // we would either have to  
				   // 1/ break the loop too early, in which
				   //    case the promise that a second
				   //    call to this function immediately
				   //    after the first one does not
				   //    change anything, would be broken,
				   // or
				   // 2/ we do as many loops as there are
				   //    levels. we know that information
				   //    is transported over one level
				   //    in each run of the loop, so this
				   //    is enough. Unfortunately, each
				   //    loop is rather expensive, so
				   //    we chose the way presented here
  std::vector<bool> flags_before_loop[2] = {flags_before[0],
                                            flags_before[1]};

				   // now for what is done in each
				   // loop: we have to fulfill several
				   // tasks at the same time, namely
				   // several mesh smoothing
				   // algorithms and mesh
				   // regularisation, by which we mean
				   // that the next mesh fulfills
				   // several requirements such as no
				   // double refinement at each face
				   // or line, etc.
				   //
				   // since doing these things at once
				   // seems almost impossible (in the
				   // first year of this library, they
				   // were done in two functions, one
				   // for refinement and one for
				   // coarsening, and most things
				   // within these were done at once,
				   // so the code was rather
				   // impossible to join into this,
				   // only, function), we do them one
				   // after each other. the order in
				   // which we do them is such that
				   // the important tasks, namely
				   // regularisation, are done last
				   // and the least important things
				   // are done the first. the
				   // following order is chosen:
				   //
				   // 0/ Only if coarsest_level_1 or
				   //    patch_level_1 is set:
				   //    clear all coarsen flags on level 1
				   //    to avoid level 0 cells being
				   //    created by coarsening.
				   // 1/ do not coarsen a cell if
				   //    'most of the neighbors' will be
				   //    refined after the step. This is
				   //    to prevent occurence of
				   //    unrefined islands.
				   // 2/ eliminate refined islands in the
				   //    interior and at the boundary. since
				   //    they don't do much harm besides
				   //    increasing the number of degrees
				   //    of freedom, doing this has a
				   //    rather low priority.
  				   // 3/ limit the level difference of
				   //    neighboring cells at each vertex.
  				   // 4/ eliminate unrefined islands. this
				   //    has higher priority since this
				   //    diminishes the approximation
				   //    properties not only of the unrefined
				   //    island, but also of the surrounding
				   //    patch.
				   // 5/ ensure patch level 1. Then the
				   //    triangulation consists of patches,
				   //    i.e. of cells that are
				   //    refined once. It follows that if
				   //    least one of the children of a cell
				   //    is or will be refined than all children
				   //    need to be refined. This step
				   //    only sets refinement flags and does
				   //    not set coarsening flags.
				   //    If the path_level_1 flag is set, then
				   //    eliminate_unrefined_islands,
				   //    eliminate_refined_inner_islands and
				   //    eliminate_refined_boundary_islands will
				   //    be fulfilled automatically and do not
				   //    need to be enforced separately.
  				   // 6/ take care of the requirement that no
				   //    double refinement is done at each face
  				   // 7/ take care that no double refinement
				   //    is done at each line in 3d or higher
				   //    dimensions.
				   // 8/ make sure that all children of each
				   //    cell are either flagged for coarsening
				   //    or none of the children is
				   //
				   // For some of these steps, it is
				   // known that they
				   // interact. Namely, it is not
				   // possible to guarantee that after
				   // step 6 another step 5 would have
				   // no effect; the same holds for
				   // the opposite order and also when
				   // taking into account step
				   // 7. however, it is important to
				   // guarantee that step five or six
				   // do not undo something that step
				   // 5 did, and step 7 not something
				   // of step 6, otherwise the
				   // requirements will not be
				   // satisfied even if the loop
				   // terminates. this is accomplished
				   // by the fact that steps 5 and 6
				   // only *add* refinement flags and
				   // delete coarsening flags
				   // (therefore, step 6 can't undo
				   // something that step 4 already
				   // did), and step 7 only deletes
				   // coarsening flags, never adds
				   // some. step 7 needs also take
				   // care that it won't tag cells for
				   // refinement for which some
				   // neighbors are more refined or
				   // will be refined.
  bool mesh_changed_in_this_loop = false;
  do
    {

      				       //////////////////////////////////////
				       // STEP 0:
				       //    Only if coarsest_level_1 or
				       //    patch_level_1 is set:
				       //    clear all coarsen flags on level 1
				       //    to avoid level 0 cells being
				       //    created by coarsening.
      if (((smooth_grid & coarsest_level_1) ||
	   (smooth_grid & patch_level_1)) && n_levels()>=2)
	{
	  typename Triangulation<dim>::active_cell_iterator
	    cell=begin_active(1),
	    endc=end_active(1);
	  
	  for (; cell!=endc; ++cell)
	    cell->clear_coarsen_flag();
	}


				       //////////////////////////////////////
				       // STEP 1:
				       //    do not coarsen a cell if 'most of
				       //    the neighbors' will be refined after
				       //    the step. This is to prevent the
				       //    occurence of unrefined islands.
				       //    If patch_level_1 is set, this will
				       //    be automatically fulfilled.
      if (smooth_grid & do_not_produce_unrefined_islands &&
	  !(smooth_grid & patch_level_1))
	{
	  cell_iterator       cell;
	  const cell_iterator endc = end();
	  
	  for (cell=begin(); cell!=endc; ++cell)
	    {
	      if (!cell->active())
		{
						   // count the children whose
						   // coarsen_flags are set
		  unsigned int n_childrens_coarsen_flags_set=0;
		  for (unsigned int c=0;
		       c<GeometryInfo<dim>::children_per_cell; ++c)
		    if (cell->child(c)->active() &&
			cell->child(c)->coarsen_flag_set())
		      ++n_childrens_coarsen_flags_set;

						   // only do
						   // something if all
						   // children are
						   // flagged for
						   // coarsening since
						   // only then are
						   // they coarsened
						   // anyway.
		  if (n_childrens_coarsen_flags_set==
		      GeometryInfo<dim>::children_per_cell)
		    {
		      unsigned int n_neighbors=0;
						       // count all
						       // neighbors
						       // that will be
						       // refined
						       // after the
						       // next step
		      unsigned int count=0;
		      for (unsigned int n=0;
			   n<GeometryInfo<dim>::faces_per_cell; ++n)
			{
			  const cell_iterator neighbor = cell->neighbor(n);
			  if (neighbor.state() == IteratorState::valid)
			    {
			      ++n_neighbors;

			      bool not_active_neighbor_will_be_coarsened=false;
			      unsigned int
				n_neighbors_childrens_coarsen_flags_set=0;
			      if (!neighbor->active())
				{
				  for (unsigned int c=0;
				       c<GeometryInfo<dim>::children_per_cell; ++c)
				    if (neighbor->child(c)->active() &&
					neighbor->child(c)->coarsen_flag_set())
				      ++n_neighbors_childrens_coarsen_flags_set;

				  if (n_neighbors_childrens_coarsen_flags_set
				      ==GeometryInfo<dim>::children_per_cell)
				    not_active_neighbor_will_be_coarsened=true;
				}


			      if ((neighbor->active() &&
				   neighbor->refine_flag_set()) ||
				  !not_active_neighbor_will_be_coarsened)
				++count;
			    }
			}
		      
		      if ((dim==1 && count==n_neighbors) ||
			  (dim>1 && (count==n_neighbors ||
				     (count>=n_neighbors-1 &&
				      n_neighbors==
				      GeometryInfo<dim>::faces_per_cell))))
			for (unsigned int c=0;
			     c<GeometryInfo<dim>::children_per_cell; ++c)
			  cell->child(c)->clear_coarsen_flag();
		    }
		  
		}  // if (!cell->active())
	    }  // for (all cells)
	} // if (smooth_grid & ...)


				       //////////////////////////////////////
				       // STEP 2:
				       //    eliminate refined islands in the
				       //    interior and at the boundary. since
				       //    they don't do much harm besides
				       //    increasing the number of degrees of
				       //    freedom, doing this has a rather low
				       //    priority.
				       //    If patch_level_1 is set, this will
				       //    be automatically fulfilled.
      if (smooth_grid & (eliminate_refined_inner_islands |
			 eliminate_refined_boundary_islands) &&
	  !(smooth_grid & patch_level_1)) 
	{
	  cell_iterator       cell;
	  const cell_iterator endc = end();
	  
	  for (cell=begin(); cell!=endc; ++cell)
	    if (!cell->active() || (cell->active() && cell->refine_flag_set()))
	      {
						 // check whether all
						 // children are
						 // active, i.e. not
						 // refined
						 // themselves. This
						 // is a precondition
						 // that the children
						 // may be coarsened
						 // away. If the cell
						 // is only flagged
						 // for refinement,
						 // then all future
						 // children will be
						 // active
		bool all_children_active = true;
		if (!cell->active())
		  for (unsigned int c=0; c<GeometryInfo<dim>::children_per_cell; ++c)
		    if (!cell->child(c)->active()) 
		      {
			all_children_active = false;
			break;
		      }

		if (all_children_active) 
		  {
						     // count number
						     // of refined and
						     // unrefined
						     // neighbors of
						     // cell.
						     // neighbors on
						     // lower levels
						     // are counted as
						     // unrefined
						     // since they can
						     // only get to
						     // the same level
						     // as this cell
						     // by the next
						     // refinement
						     // cycle
		    unsigned int unrefined_neighbors = 0,
				     total_neighbors = 0;

		    for (unsigned int n=0; n<GeometryInfo<dim>::faces_per_cell; ++n) 
		      {
			const cell_iterator neighbor = cell->neighbor(n);
			if (neighbor.state() == IteratorState::valid)
			  ++total_neighbors;

			if (neighbor.state() == IteratorState::valid)
			  {
			    if ((neighbor->active() &&
				 !neighbor->refine_flag_set()) ||
				(neighbor->level() == cell->level()-1))
			      ++unrefined_neighbors;
			    else
			      if (!neighbor->active())
								 // maybe this cell's
								 // children will be
								 // coarsened
				{
				  unsigned int tagged_children = 0;
				  for (unsigned int c=0;
				       c<GeometryInfo<dim>::children_per_cell;
				       ++c)
				    if (neighbor->child(c)->coarsen_flag_set())
				      ++tagged_children;
				  if (tagged_children ==
				      GeometryInfo<dim>::children_per_cell)
				    ++unrefined_neighbors;
				}
			  }
		      }
		    
						     // if all
						     // neighbors
						     // unrefined:
						     // mark this cell
						     // for coarsening
						     // or don't
						     // refine if
						     // marked for
						     // that
						     //
						     // also do the
						     // distinction
						     // between the
						     // two versions
						     // of the
						     // eliminate_refined_*_islands
						     // flag
						     //
						     // the last check
						     // is whether
						     // there are any
						     // neighbors at
						     // all. if not
						     // so, then we
						     // are (e.g.) on
						     // the coarsest
						     // grid with one
						     // cell, for
						     // which, of
						     // course, we do
						     // not remove the
						     // refine flag.
		    if ((unrefined_neighbors == total_neighbors)
			&&
			(((unrefined_neighbors==GeometryInfo<dim>::faces_per_cell) &&
			  (smooth_grid & eliminate_refined_inner_islands)) ||
			 ((unrefined_neighbors<GeometryInfo<dim>::faces_per_cell) &&
			  (smooth_grid & eliminate_refined_boundary_islands)) )
			&&
			(total_neighbors != 0))
		      {
			if (!cell->active())
			  for (unsigned int c=0;
			       c<GeometryInfo<dim>::children_per_cell; ++c)
			    {
			      cell->child(c)->clear_refine_flag ();
			      cell->child(c)->set_coarsen_flag ();
			    }
			else
			  cell->clear_refine_flag();
		      }
		  }
	      }
	}

				       //////////////////////////////////////
				       // STEP 3:
				       //    limit the level difference of
				       //    neighboring cells at each vertex.
      if (smooth_grid & limit_level_difference_at_vertices) 
	{
					   // store highest level one
					   // of the cells adjacent to
					   // a vertex belongs to
	  std::vector<int> vertex_level (vertices.size(), 0);
	  active_cell_iterator cell = begin_active(),
			       endc = end();
	  for (; cell!=endc; ++cell)
	    for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell;
		 ++vertex)
	      if (cell->refine_flag_set())
		vertex_level[cell->vertex_index(vertex)]
		  = std::max (vertex_level[cell->vertex_index(vertex)],
			      cell->level()+1);
	      else
		vertex_level[cell->vertex_index(vertex)]
		  = std::max (vertex_level[cell->vertex_index(vertex)],
			      cell->level());

					   // loop over all cells in
					   // reverse order. do so
					   // because we can then
					   // update the vertex levels
					   // on the and maybe already
					   // flag additional cells in
					   // this loop
	  for (cell=last_active(); cell != endc; --cell)
	    if (cell->refine_flag_set() == false) 
	      for (unsigned int vertex=0;
		   vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
		if (vertex_level[cell->vertex_index(vertex)] >
		    cell->level()+1)
		  {
						     // refine cell
						     // and update
						     // vertex levels
		    cell->clear_coarsen_flag();
		    cell->set_refine_flag();

		    for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell;
			 ++v)
		      vertex_level[cell->vertex_index(v)]
			= std::max (vertex_level[cell->vertex_index(v)],
				    cell->level()+1);

						     // now that we
						     // fixed this
						     // cell, we can
						     // safely leave
						     // this inner
						     // loop.
		    break;
		  }	  
	}

				       /////////////////////////////////////
				       // STEP 4:      
				       //    eliminate unrefined
				       //    islands. this has higher
				       //    priority since this
				       //    diminishes the
				       //    approximation properties
				       //    not only of the unrefined
				       //    island, but also of the
				       //    surrounding patch.
      if (smooth_grid & eliminate_unrefined_islands)
	{
	  active_cell_iterator cell = begin_active(),
			       endc = end();
	  for (; cell!=endc; ++cell) 
	    {
					       // if cell is already
					       // flagged for
					       // refinement: nothing
					       // to do anymore
	      if (cell->refine_flag_set())
		continue;
		  
	      unsigned int refined_neighbors = 0,
			 unrefined_neighbors = 0;
	      for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
		if (!cell->at_boundary(face))
		  {
						     // neighbor may
						     // only be on the
						     // same level or
						     // one level
						     // below because
						     // of the
						     // regularisation
						     // above
		    Assert ((cell->neighbor_level(face)==cell->level()) ||
			    (cell->neighbor_level(face)==cell->level()-1),
			    ExcInternalError());
		    if ((cell->neighbor_level(face) == cell->level()) &&
			(cell->neighbor(face)->refine_flag_set() ||
			 cell->neighbor(face)->has_children()))
		      ++refined_neighbors;
		    else
		      ++unrefined_neighbors;
		  }

	      if (unrefined_neighbors < refined_neighbors)
		{
		  if (cell->coarsen_flag_set())
		    cell->clear_coarsen_flag();
		  cell->set_refine_flag ();
		}
	    }
	}


				       /////////////////////////////////
				       // STEP 5:
				       //    ensure patch level 1.
				       //
				       //    Introduce some terminology:
				       //    - a cell that is refined
				       //      once is a patch of
				       //      level 1 simply called patch.
				       //    - a cell that is globally
				       //      refined twice is called
				       //      a patch of level 2.
				       //    - patch level n says that
				       //      the triangulation consists
				       //      of patches of level n.
				       //      This makes sense only
				       //      if the grid is already at
				       //      least n times globally
				       //      refined.
				       //
				       //    E.g. from patch level 1
				       //    follows: if at least one
				       //    of the children of a cell
				       //    is or will be refined
				       //    than enforce all
				       //    children to be
				       //    refined.

				       //    This step 4 only
				       //    sets refinement flags and
				       //    does not set coarsening
				       //    flags.
      if (smooth_grid & patch_level_1) 
	{
					   // An important assumption
					   // (A) is that before
					   // calling this function
					   // the grid was already of
					   // patch level 1.

					   // loop over all cells
					   // whose children are all
					   // active.  (By assumption
					   // (A) either all or none
					   // of the children are
					   // active).  If the refine
					   // flag of at least one of
					   // the children is set then
					   // set_refine_flag and
					   // clear_coarsen_flag of
					   // all children.
	  unsigned int n_children=GeometryInfo<dim>::children_per_cell;
	  for (cell_iterator cell = begin(); cell != end(); ++cell)
	    if (!cell->active() && cell->child(0)->active())
	      {
						 // cell is found to
						 // be a patch
		bool any_refine_flag_set=false;
		for (unsigned int i=0; i<n_children; ++i)
		  {
		    cell_iterator child=cell->child(i);
						     // check
						     // consistency:
						     // cell is really
						     // a patch,
						     // i.e. no child
						     // is refined.
		    Assert(child->active(), ExcInternalError());
		    
		    if (child->refine_flag_set())
		      {
			any_refine_flag_set=true;
			break;
		      }
		  }
		if (any_refine_flag_set)
		  {
		    for (unsigned int i=0; i<n_children; ++i)
		      {
			cell_iterator child=cell->child(i);

			child->clear_coarsen_flag();
			child->set_refine_flag();
		      }
		  }
	      }

					   // Loop over all patches of
					   // level 2, i.e. over all
					   // cells whose
					   // grandchildren are all
					   // active.  Coarsen the
					   // children (and remove the
					   // grandchildren) only if
					   // all cell->grandchild(i)
	                                   //       ->coarsen_flag_set()
					   // are set.
	  for (cell_iterator cell = begin(); cell != end(); ++cell)	    
	    if (!cell->active() &&
		!cell->child(0)->active() &&
		cell->child(0)->child(0)->active())
	      {
						 // count all coarsen
						 // flags of the
						 // grandchildren.
		unsigned int n_coarsen_flags=0;
						 // cell is not a
						 // patch (of level 1)
						 // as it has a
						 // grandchild.  Is
						 // cell a patch of
						 // level 2??
						 // Therefore: find
						 // out whether all
						 // cell->child(i) are
						 // patches
		for (unsigned int c=0; c<n_children; ++c)
		  {
		    cell_iterator child=cell->child(c);
						     // check
						     // consistency:
						     // cell is not a
						     // patch of level
						     // 1.
		    Assert(!child->active(), ExcInternalError());
		    
		    if (child->child(0)->active())
		      {
							 // child is
							 // found to
							 // be a patch
			for (unsigned int cc=0; cc<n_children; ++cc)
			  {
			    cell_iterator grand_child=child->child(cc);
							     // check
							     // consistency:
							     // child is
							     // a patch
			    Assert(grand_child->active(), ExcInternalError());

			    if (grand_child->coarsen_flag_set())
			      ++n_coarsen_flags;
			  }
		      }
		  }
		
		if (n_coarsen_flags!=n_children*n_children)
		  {
						     // clear all
						     // grandchildren's
						     // coarsen_flags
		    for (unsigned int c=0; c<n_children; ++c)
		      {
			cell_iterator child=cell->child(c);
			Assert(!child->active(), ExcInternalError());
		    
			if (child->child(0)->active())
			  {
			    for (unsigned int cc=0; cc<n_children; ++cc)
			      {
				cell_iterator grand_child=child->child(cc);
				Assert(grand_child->active(), ExcInternalError());
				
				grand_child->clear_coarsen_flag();
			      }
			  }
		      } 
		  }
	      }
	}
      
      

				       /////////////////////////////////
				       // STEP 6:
				       //    take care of the requirement that no
				       //    double refinement is done at each face
      for (active_cell_iterator cell = last_active(); cell != end(); --cell)
	if (cell->refine_flag_set() == true) 
	  {
					     // loop over neighbors of cell
	    for (unsigned int i=0; i<GeometryInfo<dim>::faces_per_cell; ++i)
	      if (cell->neighbor(i).state() == IteratorState::valid)
		{
						   // regularisation?
		  if ((cell->neighbor_level(i) == cell->level()-1)
		      &&
		      (cell->neighbor(i)->refine_flag_set() == false))
		    {
		      cell->neighbor(i)->clear_coarsen_flag();
		      cell->neighbor(i)->set_refine_flag();
		    }
		  else
		    if ((cell->neighbor_level(i) == cell->level())
			&&
			(cell->neighbor(i)->coarsen_flag_set() == true))
						       // if this cell will
						       // be refined and the
						       // neighbor may or may
						       // not be coarsened
						       // (depending on whether
						       // all children of its
						       // mother cell are tagged
						       // for coarsening), then
						       // disallow coarsening.
						       // to do so, it suffices
						       // to delete the coarsen
						       // flag from one child,
						       // namely our present
						       // neighbor
		      cell->neighbor(i)->clear_coarsen_flag ();
		}
	  }

				       //////////////////////////////////////
				       // STEP 7:
				       //    take care that no double refinement
				       //    is done at each line in 3d or higher
				       //    dimensions.
      prepare_refinement_dim_dependent ();
      
				       //////////////////////////////////////
				       // STEP 8:
				       //    make sure that all children of each
				       //    cell are either flagged for coarsening
				       //    or none of the children is
      fix_coarsen_flags ();

				       // get the refinement and coarsening
				       // flags
      std::vector<bool> flags_after_loop[2];
      save_coarsen_flags (flags_after_loop[0]);
      save_refine_flags (flags_after_loop[1]);

				       // find out whether something was
				       // changed in this loop
      mesh_changed_in_this_loop
	= ((flags_before_loop[0] != flags_after_loop[0]) ||
	   (flags_before_loop[1] != flags_after_loop[1]));

				       // set the flags for the next loop
				       // already
      flags_before_loop[0] = flags_after_loop[0];
      flags_before_loop[1] = flags_after_loop[1];
    }
  while (mesh_changed_in_this_loop);

  
				   // find out whether something was really
				   // changed in this function. Note that
				   // @p{flags_before_loop} represents the
				   // state after the last loop, i.e.
				   // the present state
  return ((flags_before[0] != flags_before_loop[0]) ||
	  (flags_before[1] != flags_before_loop[1]));
}


#if deal_II_dimension == 1

template <>
void Triangulation<1>::delete_children (cell_iterator &cell)
{
  const unsigned int dim=1;
				   // first we need to reset the
				   // neighbor pointers of the
				   // neighbors of this cell's
				   // children to this cell. This is
				   // different for one dimension,
				   // since there neighbors can have a
				   // refinement level differing from
				   // that of this cell's children by
				   // more than one level.

  Assert (!cell->child(0)->has_children() && !cell->child(1)->has_children(),
	  ExcInternalError());
  
				   // first do it for the cells to the
				   // left
  if (cell->neighbor(0).state() == IteratorState::valid)
    if (cell->neighbor(0)->has_children())
      {
	cell_iterator neighbor = cell->neighbor(0);
	Assert (neighbor->level() == cell->level(), ExcInternalError());
	
					 // right child
	neighbor = neighbor->child(1);
	while (1)
	  {
	    Assert (neighbor->neighbor(1) == cell->child(0),
		    ExcInternalError());
	    neighbor->set_neighbor (1, cell);
	    
					     // move on to further
					     // children on the
					     // boundary between this
					     // cell and its neighbor
	    if (neighbor->has_children())
	      neighbor = neighbor->child(1);
	    else
	      break;
	  }
      }

  				   // now do it for the cells to the
				   // left
  if (cell->neighbor(1).state() == IteratorState::valid)
    if (cell->neighbor(1)->has_children())
      {
	cell_iterator neighbor = cell->neighbor(1);
	Assert (neighbor->level() == cell->level(), ExcInternalError());
	
					 // left child
	neighbor = neighbor->child(0);
	while (1)
	  {
	    Assert (neighbor->neighbor(0) == cell->child(1),
		    ExcInternalError());
	    neighbor->set_neighbor (0, cell);
	    
					     // move on to further
					     // children on the
					     // boundary between this
					     // cell and its neighbor
	    if (neighbor->has_children())
	      neighbor = neighbor->child(0);
	    else
	      break;
	  }
      }


				   // delete the vertex which will not
				   // be needed anymore. This vertex
				   // is the second of the second line
				   // of the first child
  vertices_used[cell->child(0)->vertex_index(1)] = false;

				   // invalidate children.  clear user
				   // pointers, to avoid that they may
				   // appear at unwanted places later
				   // on...
  for (unsigned int child=0; child<GeometryInfo<dim>::children_per_cell; ++child)
    {
      cell->child(child)->clear_user_pointer();
      cell->child(child)->clear_user_flag();
      cell->child(child)->clear_used_flag();
    }


				   // delete pointer to children
  cell->set_children (-1);
  cell->clear_user_flag();
}

#endif


#if deal_II_dimension == 2

template <>
void Triangulation<2>::delete_children (cell_iterator &cell)
{
  const unsigned int dim=2;
				   // first we need to reset the
				   // neighbor pointers of the
				   // neighbors of this cell's
				   // children to this cell. This is
				   // different for one dimension,
				   // since there neighbors can have a
				   // refinement level differing from
				   // that of this cell's children by
				   // more than one level.  For two or
				   // more dimensions, the neighbors
				   // of the children may only be on
				   // the same level or on the level
				   // of this cell (the case that the
				   // neighbors are more refined than
				   // the children was eliminated in
				   // @p{prepare_coarsening}
  for (unsigned int child=0; child<GeometryInfo<dim>::children_per_cell; ++child)
    for (unsigned int n=0; n<GeometryInfo<dim>::faces_per_cell;
	 ++n)
      {
	const cell_iterator neighbor = cell->child(child)->neighbor(n);
					 // do nothing if at boundary
	if (neighbor.state() != IteratorState::valid)
	  continue;
	
	Assert ((neighbor->level()==cell->level()) ||
		(neighbor->level()==cell->level()+1),
		ExcInternalError());
	
					 // if the neighbor's level is
					 // the same as that of
					 // @p{cell}, then it's
					 // neighbor pointers points
					 // to this cell rather than
					 // to this cell's child. In
					 // that case we need not do
					 // anything.  If the neighbor
					 // is refined as often as are
					 // the children, we need to
					 // reset those neigbor
					 // pointers that point to the
					 // child of this cell; when
					 // resetting the neighbor
					 // pointers of neighbors of
					 // one of the children, we
					 // will also reset the
					 // neighbor pointers other
					 // children to the present
					 // cell, but this does no
					 // harm since we delete the
					 // children afterwards anyway
	if (neighbor->level() == cell->level()+1)
	  for (unsigned int neighbor_neighbor=0;
	       neighbor_neighbor<GeometryInfo<dim>::faces_per_cell;
	       ++neighbor_neighbor)
	    if (neighbor->neighbor(neighbor_neighbor) == cell->child(child))
	      neighbor->set_neighbor(neighbor_neighbor, cell);
      }

				   // delete the vertex which will not
				   // be needed anymore. This vertex
				   // is the second of the second line
				   // of the first child
  vertices_used[cell->child(0)->line(1)->vertex_index(1)] = false;

				   // clear user pointers, to avoid
				   // that they may appear at unwanted
				   // places later on...
  cell->child(0)->line(1)->clear_user_pointer();
  cell->child(0)->line(3)->clear_user_pointer();
  cell->child(3)->line(0)->clear_user_pointer();
  cell->child(3)->line(2)->clear_user_pointer();
  
				   // same for user flags
  cell->child(0)->line(1)->clear_user_flag();
  cell->child(0)->line(3)->clear_user_flag();
  cell->child(3)->line(0)->clear_user_flag();
  cell->child(3)->line(2)->clear_user_flag();
  
				   // delete the four interior lines
  cell->child(0)->line(1)->clear_used_flag();
  cell->child(0)->line(3)->clear_used_flag();
  cell->child(3)->line(0)->clear_used_flag();
  cell->child(3)->line(2)->clear_used_flag();

				   // for the four faces: if the
				   // neighbor does not itself need
				   // the subfaces, delete them. note
				   // that since dim>1 the level of a
				   // neighbor is either one less or
				   // the same as that of cell
  for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face) 
    if ((cell->neighbor(face).state() != IteratorState::valid) ||
	(cell->neighbor(face)->level() == cell->level()-1) ||
	((cell->neighbor(face)->level() == cell->level()) &&
	 !cell->neighbor(face)->has_children()))
      {
					 // delete middle vertex
	vertices_used[cell->face(face)->child(0)->vertex_index(1)] = false;
					 // delete the two subfaces
	for (unsigned int subface=0;
	     subface<GeometryInfo<dim>::subfaces_per_face; ++subface)
	  {
	    cell->face(face)->child(subface)->clear_user_pointer ();
	    cell->face(face)->child(subface)->clear_user_flag ();
	    cell->face(face)->child(subface)->clear_used_flag ();
	  }
	
	cell->face(face)->clear_children();
      }
  
				   // invalidate children
  for (unsigned int child=0; child<GeometryInfo<dim>::children_per_cell; ++child)
    {
      cell->child(child)->clear_user_pointer();
      cell->child(child)->clear_user_flag();
      cell->child(child)->clear_used_flag();
    }


				   // delete pointer to children
  cell->set_children (-1);
  cell->clear_user_flag();
}

#endif


#if deal_II_dimension == 3


template <>
void Triangulation<3>::delete_children (cell_iterator &cell)
{
  const unsigned int dim=3;
				   // first we need to reset the
				   // neighbor pointers of the
				   // neighbors of this cell's
				   // children to this cell. This is
				   // different for one dimension,
				   // since there neighbors can have a
				   // refinement level differing from
				   // that of this cell's children by
				   // more than one level.  For two or
				   // more dimensions, the neighbors
				   // of the children may only be on
				   // the same level or on the level
				   // of this cell (the case that the
				   // neighbors are more refined than
				   // the children was eliminated in
				   // @p{prepare_coarsening}
  for (unsigned int child=0; child<GeometryInfo<dim>::children_per_cell; ++child)
    for (unsigned int n=0; n<GeometryInfo<dim>::faces_per_cell; ++n)
      {
	const cell_iterator neighbor = cell->child(child)->neighbor(n);
					 // do nothing if at boundary
	if (neighbor.state() != IteratorState::valid)
	  continue;
	
	Assert ((neighbor->level()==cell->level()) ||
		(neighbor->level()==cell->level()+1),
		ExcInternalError());
	
					 // if the neighbor's level is
					 // the same as that of
					 // @p{cell}, then it's
					 // neighbor pointers points
					 // to this cell rather than
					 // to this cell's child. In
					 // that case we need not do
					 // anything.  If the neighbor
					 // is refined as often as are
					 // the children, we need to
					 // reset those neigbor
					 // pointers that point to the
					 // child of this cell; when
					 // resetting the neighbor
					 // pointers of neighbors of
					 // one of the children, we
					 // will also reset the
					 // neighbor pointers other
					 // children to the present
					 // cell, but this does no
					 // harm since we delete the
					 // children afterwards anyway
	if (neighbor->level() == cell->level()+1)
	  for (unsigned int neighbor_neighbor=0;
	       neighbor_neighbor<GeometryInfo<dim>::faces_per_cell;
	       ++neighbor_neighbor)
	    if (neighbor->neighbor(neighbor_neighbor) == cell->child(child))
	      neighbor->set_neighbor(neighbor_neighbor, cell);
      }

				   // delete the vertex which will not
				   // be needed anymore. This vertex
				   // is the vertex at the heart of
				   // this cell, which is the sixth of
				   // the first child
  vertices_used[cell->child(0)->vertex_index(7)] = false;

				   ///////////////////////////////////////
				   // delete interior quads and lines
				   //
				   // first set up a list of these
				   // line's and quad's indices
  const quad_iterator interior_quads[12]
    = {  cell->child(0)->face(1),
	 cell->child(2)->face(1),
	 cell->child(4)->face(1),
	 cell->child(6)->face(1),

	 cell->child(0)->face(3),
	 cell->child(1)->face(3),
	 cell->child(4)->face(3),
	 cell->child(5)->face(3),

	 cell->child(0)->face(5),
	 cell->child(1)->face(5),
	 cell->child(2)->face(5),
	 cell->child(3)->face(5) 
    };

  const line_iterator interior_lines[6]
    = {  cell->child(0)->line(5),
	 cell->child(0)->line(7),
	 cell->child(0)->line(11),

	 cell->child(7)->line(0),
	 cell->child(7)->line(2),
	 cell->child(7)->line(8),
    };

  				   // clear user pointers, to avoid that
				   // they may appear at unwanted places
				   // later on...
				   // same for user flags, then finally
				   // delete thes quads and lines
  for (unsigned int q=0; q<12; ++q)
    {
      interior_quads[q]->clear_user_pointer();
      interior_quads[q]->clear_user_flag();
      interior_quads[q]->clear_used_flag();
    }

  for (unsigned int l=0; l<6; ++l)
    {
      interior_lines[l]->clear_user_pointer();
      interior_lines[l]->clear_user_flag();
      interior_lines[l]->clear_used_flag();
    }

				   // for the six faces: if the
				   // neighbor does not itself need
				   // the subfaces, delete them. note
				   // that since dim>1 the level of a
				   // neighbor is either one less or
				   // the same as that of cell
  for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face) 
    if ((cell->neighbor(face).state() != IteratorState::valid) ||
	(cell->neighbor(face)->level() == cell->level()-1) ||
	((cell->neighbor(face)->level() == cell->level()) &&
	 !cell->neighbor(face)->has_children()))
      {
	quad_iterator quad = cell->face(face);
	
					 // delete middle vertex
	vertices_used[quad->child(0)->vertex_index(3)] = false;

	const line_iterator interior_lines[4]
	  = { quad->child(0)->line(1),
	      quad->child(0)->line(3),
	      quad->child(3)->line(0),
	      quad->child(3)->line(2) 
	  };

					 // delete interior lines
	for (unsigned int l=0; l<4; ++l)
	  {
	    interior_lines[l]->clear_user_pointer ();
	    interior_lines[l]->clear_user_flag ();
	    interior_lines[l]->clear_used_flag ();
	  }
	
					 // delete the four subfaces
	for (unsigned int subface=0;
	     subface<GeometryInfo<dim>::subfaces_per_face; ++subface)
	  {
	    quad->child(subface)->clear_user_pointer ();
	    quad->child(subface)->clear_user_flag ();
	    quad->child(subface)->clear_used_flag ();
	  }
	
	quad->clear_children();
      }

				   // invalidate children
  for (unsigned int child=0; child<GeometryInfo<dim>::children_per_cell; ++child)
    {
      cell->child(child)->clear_user_pointer();
      cell->child(child)->clear_user_flag();

      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
	{
					   // set flags denoting deviations from
					   // standard orientation of faces back
					   // to initialization values
	  cell->child(child)->set_face_orientation (f, true);
	  cell->child(child)->set_face_flip(f,false);
	  cell->child(child)->set_face_rotation(f,false);
	}
      
      cell->child(child)->clear_used_flag();
    }


				   // delete pointer to children
  cell->clear_children ();
  cell->clear_user_flag();

				   // now there still are the 12 lines
				   // of this hex which are refined
				   // and which may need
				   // coarsening. however, it is not
				   // so easy to decide whether they
				   // are still needed, since it does
				   // not suffice to ask the
				   // neighbors. we also need to ask
				   // those cells, which are "around
				   // the corner".
				   //
				   // to do so: first set up a list of
				   // 12 pairs of line_iterators and
				   // flags which denote whether the
				   // line's children are still needed
				   //
				   // we default to: "is not needed"
				   // because in this case, if we make
				   // an error in the code below, some
				   // lines will be deleted that in
				   // fact are needed. this will
				   // eventually be caught somewhen,
				   // because there are many checks
				   // whether an iterator points to
				   // something used. the opposite
				   // case, that we do not delete
				   // lines that are no more used, is
				   // more severe and causes a memory
				   // leak which is probably
				   // impossible to find.
  const std::pair<line_iterator,bool> line_is_needed_pairs[12]
    = {   std::make_pair(cell->line(0), false),
	  std::make_pair(cell->line(1), false),
	  std::make_pair(cell->line(2), false),
	  std::make_pair(cell->line(3), false),
	  std::make_pair(cell->line(4), false),
	  std::make_pair(cell->line(5), false),
	  std::make_pair(cell->line(6), false),
	  std::make_pair(cell->line(7), false),
	  std::make_pair(cell->line(8), false),
	  std::make_pair(cell->line(9), false),
	  std::make_pair(cell->line(10), false),
	  std::make_pair(cell->line(11), false)  };
  
				   // if in debug mode: make sure that
				   // none of the lines of this cell
				   // is twice refined; else, deleting
				   // this cell's children will result
				   // in an invalid state. also check
				   // that each of the lines for which
				   // we consider deleting the
				   // children in fact has children
				   // (the bits/coarsening_3d test
				   // tripped over this initially)
  for (unsigned int line=0; line<12; ++line)
    {
      Assert (cell->line(line)->has_children(),
              ExcInternalError());
      for (unsigned int c=0; c<2; ++c)
	Assert (!cell->line(line)->child(c)->has_children(),
		ExcInternalError());
    }
  

				   // next make a map out of this for
				   // simpler access to the flag
				   // associated with a line
  std::map<line_iterator,bool> line_is_needed (&line_is_needed_pairs[0],
					       &line_is_needed_pairs[12]);
  
				   // then ask each neighbor and their
				   // neighbors
  for (unsigned int nb=0; nb<GeometryInfo<dim>::faces_per_cell; ++nb) 
    {
      const cell_iterator neighbor = cell->neighbor(nb);
				       // do nothing if at boundary
      if (neighbor.state() != IteratorState::valid)
	continue;

      Assert (neighbor->level() == cell->level(),
	      ExcInternalError());

				       // if the neighbor itself has
				       // children, then the four
				       // lines of the common face are
				       // definitely needed
      if (neighbor->has_children())
	{
	  for (unsigned int i=0; i<4; ++i)
	    {
	      Assert (line_is_needed.find(cell->face(nb)->line(i))
		      != line_is_needed.end(),
		      ExcInternalError());
	      
	      line_is_needed[cell->face(nb)->line(i)] = true;
	    }
					   // ok for this neighbor
	  continue;
	}

				       // if the neighbor is not
				       // refined, then it may still
				       // be that one of his neighbors
				       // may need one of our
				       // lines. the present cell is
				       // also one of his neighbors,
				       // but this one is not any more
				       // refined
      for (unsigned int nb_nb=0; nb_nb<GeometryInfo<dim>::faces_per_cell;
	   ++nb_nb)
	{
	  const cell_iterator neighbor_neighbor = neighbor->neighbor(nb_nb);
					   // do nothing if at
					   // boundary
	  if (neighbor_neighbor.state() != IteratorState::valid)
	    continue;

	  Assert ((neighbor_neighbor->level() == cell->level()) ||
		  (neighbor_neighbor->level() == cell->level()-1),
		  ExcInternalError());
	  
					   // also do nothing if the
					   // neighbor is at a lower
					   // level (but then it
					   // should not be refined)
	  if (neighbor_neighbor->level() == cell->level()-1)
	    {
	      Assert (!neighbor_neighbor->has_children(),
		      ExcInternalError());
	      continue;
	    }

					   // neighbor's neighbor is
					   // on same level. if it has
					   // children, then look
					   // closer
	  if (neighbor_neighbor->has_children())
					     // if any of those cell's
					     // lines is one of those
					     // that we are interested
					     // in, then flag it
	    for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_cell; ++l)
	      {
		line_iterator line = neighbor_neighbor->line(l);

		if (line_is_needed.find(line) != line_is_needed.end())
		  line_is_needed[line] = true;
	      }
	}
    }


				   // now, if the lines are not marked
				   // as needed, we may delete their
				   // children and the midpoint
  std::map<line_iterator,bool>::iterator line_and_flag;
  for (line_and_flag=line_is_needed.begin();
       line_and_flag!=line_is_needed.end(); ++line_and_flag)
    if  (line_and_flag->second == false)
      {
	line_iterator line = line_and_flag->first;

	vertices_used[line->child(0)->vertex_index(1)] = false;

	for (unsigned int child=0; child<2; ++child)
	  {
	    line->child(child)->clear_user_pointer();
	    line->child(child)->clear_user_flag();
	    line->child(child)->clear_used_flag();
	  }

	line->clear_children();
      }
}

#endif


template <int dim>
void Triangulation<dim>::write_bool_vector (const unsigned int  magic_number1,
					    const std::vector<bool> &v,
					    const unsigned int  magic_number2,
					    std::ostream            &out)
{
  const unsigned int N = v.size();
  unsigned char *flags = new unsigned char[N/8+1];
  for (unsigned int i=0; i<N/8+1; ++i) flags[i]=0;
  
  for (unsigned int position=0; position<N; ++position)
    flags[position/8] |= (v[position] ? (1<<(position%8)) : 0);

  AssertThrow (out, ExcIO());
  
				   // format:
				   // 0. magic number
				   // 1. number of flags
				   // 2. the flags
				   // 3. magic number
  out << magic_number1 << ' ' << N << std::endl;
  for (unsigned int i=0; i<N/8+1; ++i) 
    out << static_cast<unsigned int>(flags[i]) << ' ';
  
  out << std::endl << magic_number2 << std::endl;
  
  delete[] flags;

  AssertThrow (out, ExcIO());
}


template <int dim>
void Triangulation<dim>::read_bool_vector (const unsigned int  magic_number1,
					   std::vector<bool>       &v,
					   const unsigned int  magic_number2,
					   std::istream            &in)
{
  AssertThrow (in, ExcIO());

  unsigned int magic_number;
  in >> magic_number;
  AssertThrow (magic_number==magic_number1, ExcGridReadError());

  unsigned int N;
  in >> N;
  v.resize (N);

  unsigned char *flags = new unsigned char[N/8+1];
  unsigned short int tmp;
  for (unsigned int i=0; i<N/8+1; ++i) 
    {
      in >> tmp;
      flags[i] = tmp;
    }

  for (unsigned int position=0; position!=N; ++position)
    v[position] = (flags[position/8] & (1<<(position%8)));

  in >> magic_number;
  AssertThrow (magic_number==magic_number2, ExcGridReadError());

  delete[] flags;

  AssertThrow (in, ExcIO());
}



template <int dim>
unsigned int
Triangulation<dim>::memory_consumption () const 
{
  unsigned int mem = 0;
  mem += MemoryConsumption::memory_consumption(levels);
  for (unsigned int i=0; i<levels.size(); ++i)
    mem += MemoryConsumption::memory_consumption (*levels[i]);
  mem += MemoryConsumption::memory_consumption (vertices);
  mem += MemoryConsumption::memory_consumption (vertices_used);
  mem += sizeof(boundary);
  mem += sizeof(smooth_grid);
  mem += MemoryConsumption::memory_consumption (number_cache);
  mem += MemoryConsumption::memory_consumption (faces);
  mem += MemoryConsumption::memory_consumption (*faces);  

  return mem;
}

  
#if deal_II_dimension == 1

template <>
void Triangulation<1>::update_number_cache_lines () 
{
				   ///////////////////////////////////
				   // update the number of lines
				   // on the different levels in
				   // the cache
  number_cache.n_lines_level.resize (levels.size());
  number_cache.n_lines = 0;
  for (unsigned int level=0; level<levels.size(); ++level)
    {
				       // count lines on this level
      number_cache.n_lines_level[level] = 0;
      if (levels[level]->lines.cells.size() != 0) 
	{
	  line_iterator line = begin_line (level),
			endc = (level == levels.size()-1 ?
				line_iterator(end_line()) :
				begin_line (level+1));
	  for (; line!=endc; ++line)
	    ++number_cache.n_lines_level[level];
	}
      
				       // update total number of lines
      number_cache.n_lines += number_cache.n_lines_level[level];
    }

				   // do the update for the number of
				   // active lines as well
  number_cache.n_active_lines_level.resize (levels.size());
  number_cache.n_active_lines = 0;
  for (unsigned int level=0; level<levels.size(); ++level)
    {
				       // count lines on this level
      number_cache.n_active_lines_level[level] = 0;
      if (levels[level]->lines.cells.size() != 0) 
	{
	  active_line_iterator line = begin_active_line (level),
			       endc = end_active_line (level);
	  for (; line!=endc; ++line)
	    ++number_cache.n_active_lines_level[level];
	}
      
				       // update total number of lines
      number_cache.n_active_lines += number_cache.n_active_lines_level[level];
    }
}

#endif


template <int dim>
void Triangulation<dim>::update_number_cache_lines () 
{
				   ///////////////////////////////////
				   // update the number of lines
				   // as @p dim>1, there are no levels
				   // here to take care of.
  number_cache.n_lines = 0;
  if (faces->lines.cells.size() != 0) 
    {
      line_iterator line = begin_line (),
		    endc = end_line();
      for (; line!=endc; ++line)
	++number_cache.n_lines;
    }
				   // do the update for the number of
				   // active lines as well
  number_cache.n_active_lines = 0;
  if (faces->lines.cells.size() != 0) 
    {
      active_line_iterator line = begin_active_line (),
			   endc = end_line ();
      for (; line!=endc; ++line)
	++number_cache.n_active_lines;
    }
}


#if deal_II_dimension == 1

template <>
void Triangulation<1>::update_number_cache_quads () 
{
  Assert (false, ExcInternalError());
}

#endif

#if deal_II_dimension == 2

template <>
void Triangulation<2>::update_number_cache_quads () 
{
  				   ///////////////////////////////////
				   // update the number of quads
				   // on the different levels in
				   // the cache
  number_cache.n_quads_level.resize (levels.size());
  number_cache.n_quads = 0;
  for (unsigned int level=0; level<levels.size(); ++level)
    {
				       // count quads on this level
      number_cache.n_quads_level[level] = 0;
      if (levels[level]->quads.cells.size() != 0) 
	{
	  quad_iterator quad = begin_quad (level),
			endc = (level == levels.size()-1 ?
				quad_iterator(end_quad()) :
				begin_quad (level+1));
	  for (; quad!=endc; ++quad)
	    ++number_cache.n_quads_level[level];
	}
      
				       // update total number of quads
      number_cache.n_quads += number_cache.n_quads_level[level];
    }

				   // do the update for the number
				   // of active quads as well
  number_cache.n_active_quads_level.resize (levels.size());
  number_cache.n_active_quads = 0;
  for (unsigned int level=0; level<levels.size(); ++level)
    {
				       // count quads on this level
      number_cache.n_active_quads_level[level] = 0;
      if (levels[level]->quads.cells.size() != 0) 
	{
	  active_quad_iterator quad = begin_active_quad (level),
			       endc = end_active_quad (level);
	  for (; quad!=endc; ++quad)
	    ++number_cache.n_active_quads_level[level];
	}
      
				       // update total number of quads
      number_cache.n_active_quads += number_cache.n_active_quads_level[level];
    }
}

#endif


#if deal_II_dimension > 2

template <int dim>
void Triangulation<dim>::update_number_cache_quads () 
{
  				   ///////////////////////////////////
				   // update the number of quads
				   // as @p dim>2, there are no levels
				   // here to take care of.
  number_cache.n_quads = 0;
  if (faces->quads.cells.size() != 0) 
    {
      quad_iterator quad = begin_quad (),
		    endc = end_quad();
      for (; quad!=endc; ++quad)
	++number_cache.n_quads;
    }
  
				   // do the update for the number
				   // of active quads as well
  number_cache.n_active_quads = 0;
  
  if (faces->quads.cells.size() != 0) 
    {
      active_quad_iterator quad = begin_active_quad (),
			   endc = end_quad ();
      for (; quad!=endc; ++quad)
	++number_cache.n_active_quads;
    }
}

#endif


#if deal_II_dimension == 1

template <>
void Triangulation<1>::update_number_cache_hexes () 
{
  Assert (false, ExcInternalError());
}

#endif


#if deal_II_dimension == 2

template <>
void Triangulation<2>::update_number_cache_hexes () 
{
  Assert (false, ExcInternalError());
}

#endif

#if deal_II_dimension == 3

template <int dim>
void Triangulation<dim>::update_number_cache_hexes () 
{
    				   ///////////////////////////////////
				   // update the number of hexes
				   // on the different levels in
				   // the cache
  number_cache.n_hexes_level.resize (levels.size());
  number_cache.n_hexes = 0;
  for (unsigned int level=0; level<levels.size(); ++level)
    {
				       // count hexes on this level
      number_cache.n_hexes_level[level] = 0;
      if (levels[level]->hexes.cells.size() != 0) 
	{
	  hex_iterator hex = begin_hex (level),
		      endc = (level == levels.size()-1 ?
			      hex_iterator(end_hex()) :
			      begin_hex (level+1));
	  for (; hex!=endc; ++hex)
	    ++number_cache.n_hexes_level[level];
	}
      
				       // update total number of hexes
      number_cache.n_hexes += number_cache.n_hexes_level[level];
    }

				   // do the update for the number
				   // of active hexes as well
  number_cache.n_active_hexes_level.resize (levels.size());
  number_cache.n_active_hexes = 0;
  for (unsigned int level=0; level<levels.size(); ++level)
    {
				       // count hexes on this level
      number_cache.n_active_hexes_level[level] = 0;
      if (levels[level]->hexes.cells.size() != 0) 
	{
	  active_hex_iterator hex = begin_active_hex (level),
			     endc = end_active_hex (level);
	  for (; hex!=endc; ++hex)
	    ++number_cache.n_active_hexes_level[level];
	}
      
				       // update total number of hexes
      number_cache.n_active_hexes += number_cache.n_active_hexes_level[level];
    }
}

#endif


template <int dim>
void Triangulation<dim>::update_number_cache () 
{
				   // play a little bit with switch
				   // statements without 'break's
  switch (dim)
    {
      case 3:
	    update_number_cache_hexes();
      case 2:
	    update_number_cache_quads();
      case 1:
	    update_number_cache_lines();
	    break;
	    
      default:
	    Assert (false, ExcInternalError());
    }
}



template<int dim>
Triangulation<dim>::RefinementListener::~RefinementListener ()
{}




template<int dim>
void Triangulation<dim>::
RefinementListener::pre_refinement_notification (const Triangulation<dim> &)
{}



template<int dim>
void Triangulation<dim>::
RefinementListener::post_refinement_notification (const Triangulation<dim> &)
{}



template<int dim>
void
Triangulation<dim>::add_refinement_listener (RefinementListener &listener) const
{
    refinement_listeners.push_back (&listener);
}



template<int dim>
void
Triangulation<dim>::remove_refinement_listener (RefinementListener &listener) const
{
    typename std::list<RefinementListener *>::iterator p =
	find (refinement_listeners.begin (), refinement_listeners.end (), &listener);
    Assert (p != refinement_listeners.end (), ExcInternalError ());
    refinement_listeners.erase (p);
}


// explicit instantiations
template class Triangulation<deal_II_dimension>;

DEAL_II_NAMESPACE_CLOSE

