/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998, 1999 */


#include <grid/tria.h>
#include <grid/tria_levels.h>
#include <grid/tria_boundary.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/geometry_info.h>
#include <basic/magic_numbers.h>
#include <lac/vector.h>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <map>
#include <cmath>



template <int dim>
const StraightBoundary<dim>& Triangulation<dim>::straight_boundary
= StraightBoundary<dim>();




template <int dim>
Triangulation<dim>::Triangulation (const MeshSmoothing smooth_grid) :
		Subscriptor (),
		smooth_grid(smooth_grid)
{
				   // set default boundary for all possible components
  for (unsigned int i=0;i<255;++i)
    {
      boundary[i] = &straight_boundary;
      straight_boundary.subscribe();
    }
};



template <int dim>
Triangulation<dim>::Triangulation (const Triangulation<dim> &) 
		:
		Subscriptor ()      // do not set any subscriptors; anyway,
			       // calling this constructor is an error!
{
  Assert (false, ExcInternalError());
};



template <int dim>
Triangulation<dim>::~Triangulation ()
{
  for (unsigned int i=0; i<levels.size(); ++i)
    delete levels[i];
  levels.clear ();

  for (unsigned int i=0;i<255;++i)
    boundary[i]->unsubscribe ();
};



template <int dim>
void Triangulation<dim>::clear () 
{
				   // only allow this operation if there
				   // are no subscribers any more
  Assert (n_subscriptions() == 0, ExcInternalError());
  
  for (unsigned int i=0; i<levels.size(); ++i)
    delete levels[i];
  levels.clear ();

  vertices.clear ();
  vertices_used.clear ();
  
  for (unsigned int i=0; i<255; ++i)
    {
      boundary[i]->unsubscribe ();
      boundary[i] = &straight_boundary;
      straight_boundary.subscribe ();
    };

  number_cache = TriaNumberCache<dim>();
};


template <int dim>
void
Triangulation<dim>::set_boundary (const unsigned int number,
				  const Boundary<dim>& boundary_object)
{
  Assert(number<255, ExcIndexRange(number,0,255));
  
  boundary[number]->unsubscribe ();
  boundary[number] = &boundary_object;
  boundary_object.subscribe();
};



template <int dim>
const Boundary<dim> &
Triangulation<dim>::get_boundary (const unsigned int number) const 
{
  Assert(number<255, ExcIndexRange(number,0,255));
  
  return *(boundary[number]);
};




/*--------- Put the next functions a bit out-or-order to avoid use before
  --------- explicit specialization, which is not allowed.                */

#if deal_II_dimension == 1
template <>
TriaDimensionInfo<1>::cell_iterator
Triangulation<1>::begin (const unsigned int level) const {
  return begin_line (level);
};



template <>
TriaDimensionInfo<1>::raw_cell_iterator
Triangulation<1>::end () const {
  return end_line ();
};

#endif


#if deal_II_dimension == 2

template <>
TriaDimensionInfo<2>::cell_iterator
Triangulation<2>::begin (const unsigned int level) const {
  return begin_quad (level);
};



template <>
TriaDimensionInfo<2>::raw_cell_iterator
Triangulation<2>::end () const {
  return end_quad ();
};

#endif



#if deal_II_dimension == 3

template <>
TriaDimensionInfo<3>::cell_iterator
Triangulation<3>::begin (const unsigned int level) const {
  return begin_hex (level);
};



template <>
TriaDimensionInfo<3>::raw_cell_iterator
Triangulation<3>::end () const {
  return end_hex ();
};

#endif




/*-----------------------------------------------------------------*/



template <int dim>
void Triangulation<dim>::copy_triangulation (const Triangulation<dim> &old_tria)
{
  Assert (vertices.size() == 0, ExcTriangulationNotEmpty());
  Assert (levels.size () == 0, ExcTriangulationNotEmpty());

  Assert (old_tria.levels.size() != 0, ExcInternalError());
  Assert (old_tria.vertices.size() != 0, ExcInternalError());
  
				   // copy normal elements
  vertices      = old_tria.vertices;
  vertices_used = old_tria.vertices_used;
  smooth_grid   = old_tria.smooth_grid;

  for (unsigned i=0;i<255;++i)
    {
      boundary[i]->unsubscribe ();
      boundary[i]      = old_tria.boundary[i];
      boundary[i]->subscribe ();
    }

  levels.reserve (old_tria.levels.size());
  for (unsigned int level=0; level<old_tria.levels.size(); ++level)
    levels.push_back (new TriangulationLevel<dim>(*old_tria.levels[level]));

  number_cache = old_tria.number_cache;
  
				   // note that we need not copy the
				   // subscriptor!
};



template <int dim>
void Triangulation<dim>::block_write (ostream &out) const 
{
				   // write out the level information
  out << levels.size() << endl << '[';
  for (unsigned int level=0; level<levels.size(); ++level)
    levels[level]->block_write (out);
  out << ']';
  
				   // write out the vertices
  out << '[';
  out.write (reinterpret_cast<const char*>(vertices.begin()),
	     reinterpret_cast<const char*>(vertices.end())
	     - reinterpret_cast<const char*>(vertices.begin()));
  out << ']';

  write_bool_vector (0, vertices_used, 0, out);
};



template <int dim>
void Triangulation<dim>::block_read (istream &) 
{
  Assert (false, ExcNotImplemented());
};


  

#if deal_II_dimension == 1

template <>
void Triangulation<1>::create_triangulation (const vector<Point<1> >    &v,
					     const vector<CellData<1> > &cells,
					     const SubCellData &subcelldata) {
				   // note: since no boundary information
				   // can be given in one dimension, the
				   // #subcelldata# field is ignored. (only
				   // used for error checking, which is a
				   // good idea in any case)
  
  const unsigned int dim=1;
  
  Assert (vertices.size() == 0, ExcTriangulationNotEmpty());
  Assert (levels.size() == 0, ExcTriangulationNotEmpty());
				   // check that no forbidden arrays are used
  Assert (subcelldata.check_consistency(dim), ExcInternalError());

				   // copy vertices
  vertices = v;
  vertices_used = vector<bool> (v.size(), true);
    
				   // store the indices of the lines which
				   // are adjacent to a given vertex
  vector<vector<int> > lines_at_vertex (v.size());

				   // reserve enough space
  levels.push_back (new TriangulationLevel<dim>);
  levels[0]->TriangulationLevel<0>::reserve_space (cells.size(), dim);
  levels[0]->TriangulationLevel<1>::reserve_space (cells.size());
  
				   // make up cells
  raw_line_iterator next_free_line = begin_raw_line ();
  for (unsigned int cell=0; cell<cells.size(); ++cell) 
    {
      while (next_free_line->used())
	++next_free_line;
      
      next_free_line->set (Line (cells[cell].vertices[0], cells[cell].vertices[1]));
      next_free_line->set_used_flag ();
      next_free_line->set_material_id (cells[cell].material_id);
      next_free_line->clear_user_pointer ();
      
				       // note that this cell
				       // is adjacent to these vertices
      lines_at_vertex[cells[cell].vertices[0]].push_back (cell);
      lines_at_vertex[cells[cell].vertices[1]].push_back (cell);
    };


				   // some security tests
  unsigned int boundary_nodes = 0;
  for (unsigned int i=0; i<lines_at_vertex.size(); ++i)
    switch (lines_at_vertex[i].size()) 
      {
	case 1:      // this vertex has only one adjacent line
	      ++boundary_nodes;
	      break;
	case 2:
	      break;
	default:     // a node must have one or two adjacent lines
	      AssertThrow (false, ExcInternalError());
      };

				   // assert there are no more than two boundary
				   // nodes
  AssertThrow (boundary_nodes == 2, ExcInternalError());


  				   // update neighborship info
  active_line_iterator line = begin_active_line ();
				   // for all lines
  for (; line!=end(); ++line)
				     // for each of the two vertices
    for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
				       // if first cell adjacent to this
				       // vertex is the present one, then
				       // the neighbor is the second adjacent
				       // cell and vice versa
      if (lines_at_vertex[line->vertex_index(vertex)][0] == line->index())
	if (lines_at_vertex[line->vertex_index(vertex)].size() == 2) 
	  {
	    const cell_iterator neighbor (const_cast<Triangulation<1>*>(this),
					  0,              // level
					  lines_at_vertex[line->vertex_index(vertex)][1]);
	    line->set_neighbor (vertex, neighbor);
	  }
	else
					   // no second adjacent cell entered
					   // -> cell at boundary
	  line->set_neighbor (vertex, end());
      else
					 // present line is not first adjacent
					 // one -> first adjacent one is neighbor
	{
	  const cell_iterator neighbor (const_cast<Triangulation<1>*>(this),
					0,              // level
					lines_at_vertex[line->vertex_index(vertex)][0]);
	  line->set_neighbor (vertex, neighbor);
	};


				   // re-compute numbers of lines, etc
  update_number_cache ();  
};

#endif



#if deal_II_dimension == 2

template <>
void Triangulation<2>::create_triangulation (const vector<Point<2> >    &v,
					     const vector<CellData<2> > &c,
					     const SubCellData          &subcelldata) {
  const unsigned int dim=2;

  Assert (vertices.size() == 0, ExcTriangulationNotEmpty());
  Assert (levels.size() == 0, ExcTriangulationNotEmpty());
				   // check that no forbidden arrays are used
  Assert (subcelldata.check_consistency(dim), ExcInternalError());

				   // copy vertices
  vertices = v;
  vertices_used = vector<bool> (v.size(), true);

				   // copy cells. This is needed since we
				   // may need to change entries
  vector<CellData<2> > cells(c);

  
				   // make up a list of the needed lines
				   // each line is a pair of vertices. The list
				   // is kept sorted and it is guaranteed
				   // that each line is inserted only once.
				   // While the key of such an entry is the
				   // pair of vertices, the thing it points
				   // to is an iterator pointing to the line
				   // object itself. In the first run, these
				   // iterators are all invalid ones, but they
				   // are filled afterwards
  map<pair<int,int>,line_iterator> needed_lines;
  for (unsigned int cell=0; cell<cells.size(); ++cell)
    {
      for (unsigned int vertex=0; vertex<4; ++vertex)
	AssertThrow ((0<=cells[cell].vertices[vertex]) &&
		     (cells[cell].vertices[vertex]<static_cast<signed int>(vertices.size())),
		     ExcInvalidVertexIndex (cell, cells[cell].vertices[vertex], vertices.size()));
      
      pair<int,int> line_vertices[4] = {   // note the order of the vertices
	    make_pair (cells[cell].vertices[0], cells[cell].vertices[1]),
	    make_pair (cells[cell].vertices[1], cells[cell].vertices[2]),
	    make_pair (cells[cell].vertices[3], cells[cell].vertices[2]),
	    make_pair (cells[cell].vertices[0], cells[cell].vertices[3])};

				       // note the following: if the sense
				       // of the vertices of a cell is correct,
				       // but the vertices are given in an
				       // order which makes the sense of one line
				       // ambiguous when viewed from the two
				       // adjacent cells, we can heal this by
				       // shifting the vertex indices of one
				       // cell by two (diagonally exchanging
				       // the two vertices from which the
				       // four lines originate and to which
				       // they converge).
				       // If two lines are wrong, we could heal
				       // this by rotating by one or three
				       // vertices, but deciding this is
				       // difficult and not implemented.
      for (unsigned int line=0; line<4; ++line)
	if (needed_lines.find(make_pair(line_vertices[line].second,
					line_vertices[line].first))
	    !=
	    needed_lines.end())
	  {
					     // rotate vertex numbers
	    swap (cells[cell].vertices[0], cells[cell].vertices[2]);
	    swap (cells[cell].vertices[1], cells[cell].vertices[3]);
					     // remake lines
	    line_vertices[0]
	      = make_pair (cells[cell].vertices[0], cells[cell].vertices[1]);
	    line_vertices[1]
	      = make_pair (cells[cell].vertices[1], cells[cell].vertices[2]);
	    line_vertices[2]
	      = make_pair (cells[cell].vertices[0], cells[cell].vertices[3]);
	    line_vertices[3]
	      = make_pair (cells[cell].vertices[3], cells[cell].vertices[2]);
					     // allow for only one such
					     // rotation
	    break;
	  };
      
      
      for (unsigned int line=0; line<4; ++line)
	{
					   // assert that the line was not
					   // already inserted in reverse
					   // order. This happens in spite of
					   // the vertex rotation above, if the
					   // sense of the cell was incorrect.
					   //
					   // Here is what usually happened when
					   // this exception is thrown:
					   // consider these two cells
					   // and the vertices
					   //  3---4---5
					   //  |   |   |
					   //  0---1---2
					   // If in the input vector the
					   // two cells are given with
					   // vertices <0 1 4 3> and
					   // <4 1 2 5>, in the first cell
					   // the middle line would have
					   // direction 1->4, while in
					   // the second it would be 4->1.
					   // This will cause the exception.
	  AssertThrow (needed_lines.find(make_pair(line_vertices[line].second,
						   line_vertices[line].first))
		       ==
		       needed_lines.end(),
		       ExcGridHasInvalidCell(cell));
		  
					   // insert line, with invalid iterator
					   // if line already exists, then
					   // nothing bad happens here
	  needed_lines[line_vertices[line]] = end_line();
	};
    };
  

				   // check the every vertex has
				   // at least two adjacent lines
  if (true) 
    {
      vector<unsigned short int> vertex_touch_count (v.size(), 0);
      map<pair<int,int>,line_iterator>::iterator i;
      for (i=needed_lines.begin(); i!=needed_lines.end(); i++) 
	{
					   // touch the vertices of this line
	  ++vertex_touch_count[i->first.first];
	  ++vertex_touch_count[i->first.second];
	};

				       // assert minimum touch count is at
				       // least two
      AssertThrow (* (min_element(vertex_touch_count.begin(),
				  vertex_touch_count.end())) >= 2,
		   ExcGridHasInvalidVertices());
    };
	
  				   // reserve enough space
  levels.push_back (new TriangulationLevel<dim>);
  levels[0]->TriangulationLevel<0>::reserve_space (cells.size(), dim);
  levels[0]->TriangulationLevel<1>::reserve_space (needed_lines.size());
  levels[0]->TriangulationLevel<2>::reserve_space (cells.size());

				   // make up lines
  if (true) 
    {
      raw_line_iterator line = begin_raw_line();
      map<pair<int,int>,line_iterator>::iterator i;
      for (i = needed_lines.begin(); line!=end_line(); ++line, ++i) 
	{
	  line->set (Line(i->first.first, i->first.second));
	  line->set_used_flag ();
	  line->clear_user_flag ();
	  line->clear_user_pointer ();
	  i->second = line;
	};
    };


				   // store for each line index
				   // the adjacent cells
  map<int,vector<cell_iterator> > adjacent_cells;

				   // finally make up cells
  if (true) 
    {
      raw_cell_iterator cell = begin_raw_quad();
      for (unsigned int c=0; c<cells.size(); ++c, ++cell)
	{
					   // list of iterators of lines
	  const line_iterator lines[4] = {
		needed_lines[make_pair(cells[c].vertices[0], cells[c].vertices[1])],
		needed_lines[make_pair(cells[c].vertices[1], cells[c].vertices[2])],
		needed_lines[make_pair(cells[c].vertices[3], cells[c].vertices[2])],
		needed_lines[make_pair(cells[c].vertices[0], cells[c].vertices[3])]};
	  
	  cell->set (Quad(lines[0]->index(),
			  lines[1]->index(),
			  lines[2]->index(),
			  lines[3]->index()));
	  
	  cell->set_used_flag ();
	  cell->set_material_id (cells[c].material_id);
	  cell->clear_user_pointer ();
					   // note that this cell is adjacent
					   // to the four lines
	  for (unsigned int line=0; line<4; ++line)
	    adjacent_cells[lines[line]->index()].push_back (cell);
	  
					   // make some checks on the vertices
					   // and their ordering
	  Assert (lines[0]->vertex_index(0) == lines[3]->vertex_index(0),
		  ExcInternalErrorOnCell(c));
	  Assert (lines[0]->vertex_index(1) == lines[1]->vertex_index(0),
		  ExcInternalErrorOnCell(c));
	  Assert (lines[1]->vertex_index(1) == lines[2]->vertex_index(1),
		  ExcInternalErrorOnCell(c));
	  Assert (lines[2]->vertex_index(0) == lines[3]->vertex_index(1),
		  ExcInternalErrorOnCell(c));
	};
    };
  

  for (line_iterator line=begin_line(); line!=end_line(); ++line) 
    {
      const unsigned int n_adj_cells = adjacent_cells[line->index()].size();
				       // assert that every line has one or
				       // two adjacent cells
      AssertThrow ((n_adj_cells >= 1) &&
		   (n_adj_cells <= 2),
		   ExcInternalError());

				       // if only one cell: line is at
				       // boundary -> give it the boundary
				       // indicator zero by default
      if (n_adj_cells == 1)
	line->set_boundary_indicator (0);
      else
					 // interior line -> 255
      	line->set_boundary_indicator (255);
    };

				   // set boundary indicators where given
  vector<CellData<1> >::const_iterator boundary_line
    = subcelldata.boundary_lines.begin();
  vector<CellData<1> >::const_iterator end_boundary_line
    = subcelldata.boundary_lines.end();
  for (; boundary_line!=end_boundary_line; ++boundary_line) 
    {
      line_iterator line;
      pair<int,int> line_vertices(make_pair(boundary_line->vertices[0],
					    boundary_line->vertices[1]));
      if (needed_lines.find(line_vertices) != needed_lines.end())
					 // line found in this direction
	line = needed_lines[line_vertices];
      else 
	{
					   // look whether it exists in
					   // reverse direction
	  swap (line_vertices.first, line_vertices.second);
	  if (needed_lines.find(line_vertices) != needed_lines.end())
	    line = needed_lines[line_vertices];
	  else 
	    {
					       // line does not exist
	      AssertThrow (false, ExcLineInexistant(line_vertices.first,
						    line_vertices.second));
	      line = end_line();
	    };
	};

				       // Assert that only exterior lines
				       // are given a boundary indicator
      AssertThrow (line->boundary_indicator() == 0,
		   ExcInteriorLineCantBeBoundary());

      line->set_boundary_indicator (boundary_line->material_id);
    };



  

				   // finally update neighborship info
  for (cell_iterator cell=begin(); cell!=end(); ++cell)
    for (unsigned int side=0; side<4; ++side)
      if (adjacent_cells[cell->line(side)->index()][0] == cell)
					 // first adjacent cell is this one
	{
	  if (adjacent_cells[cell->line(side)->index()].size() == 2)
					     // there is another adjacent cell
	    cell->set_neighbor (side,
				adjacent_cells[cell->line(side)->index()][1]);
	}
				   // first adjacent cell is not this one,
				   // -> it must be the neighbor we are
				   // looking for
      else
	cell->set_neighbor (side,
			    adjacent_cells[cell->line(side)->index()][0]);


				   // re-compute numbers of lines, etc
  update_number_cache ();  
};

#endif


#if deal_II_dimension == 3

/**
 * Invent an object which compares two Quads against each other. This
 * comparison is needed in order to establish a map of Quads to iterators
 * in the Triangulation<3>::create_triangulation function.
 *
 * Since this comparison is not canonical, we do not include it into the
 * general Quad class.
 */
struct QuadComparator {
    bool operator () (const Quad &q1, const Quad &q2) {
				       // here is room to optimize the
				       // repeated equality test of
				       // the previous lines, but I don't
				       // care at present
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
    };
};



template <>
void Triangulation<3>::create_triangulation (const vector<Point<3> >    &v,
					     const vector<CellData<3> > &c,
					     const SubCellData          &subcelldata) {
  const unsigned int dim=3;

  Assert (vertices.size() == 0, ExcTriangulationNotEmpty());
  Assert (levels.size() == 0, ExcTriangulationNotEmpty());
				   // check that no forbidden arrays are used
  Assert (subcelldata.check_consistency(dim), ExcInternalError());

				   // copy vertices
  vertices = v;
  vertices_used = vector<bool> (v.size(), true);

				   // copy cells. This is needed since we
				   // may need to change entries
  vector<CellData<3> > cells(c);

				   ///////////////////////////////////////
				   // first set up some collections of data
				   //
				   // make up a list of the needed lines
				   //
				   // each line is a pair of vertices. The list
				   // is kept sorted and it is guaranteed
				   // that each line is inserted only once.
				   // While the key of such an entry is the
				   // pair of vertices, the thing it points
				   // to is an iterator pointing to the line
				   // object itself. In the first run, these
				   // iterators are all invalid ones, but they
				   // are filled afterwards
				   // same applies for the quads
  map<pair<int,int>,line_iterator> needed_lines;
  for (unsigned int cell=0; cell<cells.size(); ++cell)
    {

				       // check whether vertex
				       // indices are valid ones
      for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
	AssertThrow ((0<=cells[cell].vertices[vertex]) &&
		     (cells[cell].vertices[vertex]<static_cast<signed int>(vertices.size())),
		     ExcInvalidVertexIndex (cell, cells[cell].vertices[vertex], vertices.size()));
      
      pair<int,int> line_vertices[12] = {   // note the order of the vertices
					     // front face
	    make_pair (cells[cell].vertices[0], cells[cell].vertices[1]),
	    make_pair (cells[cell].vertices[1], cells[cell].vertices[2]),
	    make_pair (cells[cell].vertices[3], cells[cell].vertices[2]),
	    make_pair (cells[cell].vertices[0], cells[cell].vertices[3]),
					     // back face
	    make_pair (cells[cell].vertices[4], cells[cell].vertices[5]),
	    make_pair (cells[cell].vertices[5], cells[cell].vertices[6]),
	    make_pair (cells[cell].vertices[7], cells[cell].vertices[6]),
	    make_pair (cells[cell].vertices[4], cells[cell].vertices[7]),
					     // connects of front and back face
	    make_pair (cells[cell].vertices[0], cells[cell].vertices[4]),
	    make_pair (cells[cell].vertices[1], cells[cell].vertices[5]),
	    make_pair (cells[cell].vertices[2], cells[cell].vertices[6]),
	    make_pair (cells[cell].vertices[3], cells[cell].vertices[7])
      };

				       // in the 2d code, some tests were performed
				       // which may heal a problem with quads that
				       // are rotated such that the lines don't
				       // fit snuggly any more. I don't know how
				       // to do this in 3d also, so I leve it for
				       // future student generations.
				       //
				       // however, check that the line does not
				       // exist in the other direction
      for (unsigned int line=0; line<12; ++line)
	{
					   // assert that the line was not
					   // already inserted in reverse
					   // order.
	  AssertThrow (needed_lines.find(make_pair(line_vertices[line].second,
						   line_vertices[line].first))
		       ==
		       needed_lines.end(),
		       ExcGridHasInvalidCell(cell));
		  
					   // insert line, with invalid iterator
					   // if line already exists, then
					   // nothing bad happens here
	  needed_lines[line_vertices[line]] = end_line();
	};
    };
  

				   /////////////////////////////////
				   // now for some sanity-checks:
				   //
				   // check the every vertex has
				   // at least two adjacent lines
  if (true) 
    {
      vector<unsigned short int> vertex_touch_count (v.size(), 0);
      map<pair<int,int>,line_iterator>::iterator i;
      for (i=needed_lines.begin(); i!=needed_lines.end(); i++) 
	{
					   // touch the vertices of this line
	  ++vertex_touch_count[i->first.first];
	  ++vertex_touch_count[i->first.second];
	};

				       // assert minimum touch count is at
				       // least two
      AssertThrow (* (min_element(vertex_touch_count.begin(),
			     vertex_touch_count.end())) >= 2,
	      ExcGridHasInvalidVertices());
    };


				   ///////////////////////////////////
				   // actually set up data structures
				   // for the lines
  				   // reserve enough space
  levels.push_back (new TriangulationLevel<dim>);
  levels[0]->TriangulationLevel<0>::reserve_space (cells.size(), dim);
  levels[0]->TriangulationLevel<1>::reserve_space (needed_lines.size());

				   // make up lines
  if (true) 
    {
      raw_line_iterator line = begin_raw_line();
      map<pair<int,int>,line_iterator>::iterator i;
      for (i = needed_lines.begin(); line!=end_line(); ++line, ++i) 
	{
	  line->set (Line(i->first.first, i->first.second));
	  line->set_used_flag ();
	  line->clear_user_flag ();
	  line->clear_user_pointer ();

					   // now set the iterator for this
					   // line
	  i->second = line;
	};
    };


				   ///////////////////////////////////////////
				   // make up the quads of this triangulation
				   //
				   // same thing: the iterators are set
				   // to the invalid value at first, we only
				   // collect the data now

				   // note that QuadComparator is a class
				   // declared and defined in this file
  map<Quad,quad_iterator,QuadComparator> needed_quads;
  for (unsigned int cell=0; cell<cells.size(); ++cell) 
    {
				       // the faces are quads which consist
				       // of four numbers denoting the index
				       // of the four lines bounding th
				       // quad. we can get this index by
				       // asking #needed_lines# for an
				       // iterator to this line, dereferencing
				       // it and thus return an iterator into
				       // the #lines# array of the
				       // triangulation, which is already set
				       // up. we can then ask this iterator
				       // for its index within the present
				       // level (the level is zero, of course)
				       //
				       // to make things easier, we don't
				       // creare the lines (pairs of their
				       // vertex indices) in place, but before
				       // they are really needed. This is just
				       // copied from above.
      pair<int,int> line_list[12] = {   // note the order of the vertices
					     // front face
	    make_pair (cells[cell].vertices[0], cells[cell].vertices[1]),
	    make_pair (cells[cell].vertices[1], cells[cell].vertices[2]),
	    make_pair (cells[cell].vertices[3], cells[cell].vertices[2]),
	    make_pair (cells[cell].vertices[0], cells[cell].vertices[3]),
					     // back face
	    make_pair (cells[cell].vertices[4], cells[cell].vertices[5]),
	    make_pair (cells[cell].vertices[5], cells[cell].vertices[6]),
	    make_pair (cells[cell].vertices[7], cells[cell].vertices[6]),
	    make_pair (cells[cell].vertices[4], cells[cell].vertices[7]),
					     // connects of front and back face
	    make_pair (cells[cell].vertices[0], cells[cell].vertices[4]),
	    make_pair (cells[cell].vertices[1], cells[cell].vertices[5]),
	    make_pair (cells[cell].vertices[2], cells[cell].vertices[6]),
	    make_pair (cells[cell].vertices[3], cells[cell].vertices[7])
      };

      Quad faces[6]
	= {
					       // front face
	      Quad (needed_lines[line_list[0]]->index(),
		    needed_lines[line_list[1]]->index(),
		    needed_lines[line_list[2]]->index(),
		    needed_lines[line_list[3]]->index()),
					       // back face
	      Quad (needed_lines[line_list[4]]->index(),
		    needed_lines[line_list[5]]->index(),
		    needed_lines[line_list[6]]->index(),
		    needed_lines[line_list[7]]->index()),
					       // bottom face
	      Quad (needed_lines[line_list[0]]->index(),
		    needed_lines[line_list[9]]->index(),
		    needed_lines[line_list[4]]->index(),
		    needed_lines[line_list[8]]->index()),
					       // right face
	      Quad (needed_lines[line_list[9]]->index(),
		    needed_lines[line_list[5]]->index(),
		    needed_lines[line_list[10]]->index(),
		    needed_lines[line_list[1]]->index()),
					       // top face
	      Quad (needed_lines[line_list[2]]->index(),
		    needed_lines[line_list[10]]->index(),
		    needed_lines[line_list[6]]->index(),
		    needed_lines[line_list[11]]->index()),
					       // left face
	      Quad (needed_lines[line_list[8]]->index(),
		    needed_lines[line_list[7]]->index(),
		    needed_lines[line_list[11]]->index(),
		    needed_lines[line_list[3]]->index())    };

				       // in the 2d code, some tests were performed
				       // which may heal a problem with hexes that
				       // are rotated such that the quads don't
				       // fit snuggly any more. I don't know how
				       // to do this in here, so I leve it for
				       // future student generations.
      for (unsigned int quad=0; quad<6; ++quad)
					 // insert quad, with invalid iterator
					 // if quad already exists, then
					 // nothing bad happens here
	  needed_quads[faces[quad]] = end_quad();
    };
  

				   /////////////////////////////////
				   // enter the resulting quads into
				   // the arrays of the Triangulation
				   //
				   // first reserve enough space
  levels[0]->TriangulationLevel<2>::reserve_space (needed_quads.size());
  if (true) 
    {
      raw_quad_iterator quad = begin_raw_quad();
      map<Quad,quad_iterator>::iterator q;
      for (q = needed_quads.begin(); quad!=end_quad(); ++quad, ++q)
	{
	  quad->set (q->first);
	  quad->set_used_flag ();
	  quad->clear_user_flag ();
	  quad->clear_user_pointer ();

					   // now set the iterator for this
					   // quad
	  q->second = quad;
	};
    };

				   /////////////////////////////////
				   // finally create the cells
  levels[0]->TriangulationLevel<3>::reserve_space (cells.size());  

				   // store for each quad index
				   // the adjacent cells
  map<int,vector<cell_iterator> > adjacent_cells;

				   // finally make up cells
  if (true) 
    {
      raw_cell_iterator cell = begin_raw_hex();
      for (unsigned int c=0; c<cells.size(); ++c, ++cell)
	{
					   // first find for each of the
					   // cells the quad iterator of
					   // the respective faces.
					   //
					   // to this end, set up the
					   // lines of this cell and find
					   // the quads that are bounded by
					   // these lines; these are then
					   // the faces of the present cell
	  pair<int,int> line_list[12] = {   // note the order of the vertices
						 // front face
		make_pair (cells[c].vertices[0], cells[c].vertices[1]),
		make_pair (cells[c].vertices[1], cells[c].vertices[2]),
		make_pair (cells[c].vertices[3], cells[c].vertices[2]),
		make_pair (cells[c].vertices[0], cells[c].vertices[3]),
						 // back face
		make_pair (cells[c].vertices[4], cells[c].vertices[5]),
		make_pair (cells[c].vertices[5], cells[c].vertices[6]),
		make_pair (cells[c].vertices[7], cells[c].vertices[6]),
		make_pair (cells[c].vertices[4], cells[c].vertices[7]),
						 // connects of front and back face
		make_pair (cells[c].vertices[0], cells[c].vertices[4]),
		make_pair (cells[c].vertices[1], cells[c].vertices[5]),
		make_pair (cells[c].vertices[2], cells[c].vertices[6]),
		make_pair (cells[c].vertices[3], cells[c].vertices[7])
	  };

	  Quad faces[6]
	    = {
						   // front face
		  Quad (needed_lines[line_list[0]]->index(),
			needed_lines[line_list[1]]->index(),
			needed_lines[line_list[2]]->index(),
			needed_lines[line_list[3]]->index()),
						   // back face
		  Quad (needed_lines[line_list[4]]->index(),
			needed_lines[line_list[5]]->index(),
			needed_lines[line_list[6]]->index(),
			needed_lines[line_list[7]]->index()),
						   // bottom face
		  Quad (needed_lines[line_list[0]]->index(),
			needed_lines[line_list[9]]->index(),
			needed_lines[line_list[4]]->index(),
			needed_lines[line_list[8]]->index()),
						   // right face
		  Quad (needed_lines[line_list[9]]->index(),
			needed_lines[line_list[5]]->index(),
			needed_lines[line_list[10]]->index(),
			needed_lines[line_list[1]]->index()),
						   // top face
		  Quad (needed_lines[line_list[2]]->index(),
			needed_lines[line_list[10]]->index(),
			needed_lines[line_list[6]]->index(),
			needed_lines[line_list[11]]->index()),
						   // left face
		  Quad (needed_lines[line_list[8]]->index(),
			needed_lines[line_list[7]]->index(),
			needed_lines[line_list[11]]->index(),
			needed_lines[line_list[3]]->index())    };

					   // get the iterators corresponding
					   // to the faces
	  const quad_iterator face_iterator[6] = {
		needed_quads[faces[0]],
		needed_quads[faces[1]],
		needed_quads[faces[2]],
		needed_quads[faces[3]],
		needed_quads[faces[4]],
		needed_quads[faces[5]]};

					   // make the cell out of these
					   // iterators
	  cell->set (Hexahedron(face_iterator[0]->index(),
				face_iterator[1]->index(),
				face_iterator[2]->index(),
				face_iterator[3]->index(),
				face_iterator[4]->index(),
				face_iterator[5]->index()));
	  
	  cell->set_used_flag ();
	  cell->set_material_id (cells[c].material_id);
	  cell->clear_user_flag ();
	  cell->clear_user_pointer ();
					   // note that this cell is adjacent
					   // to the four lines
	  for (unsigned int quad=0; quad<6; ++quad)
	    adjacent_cells[face_iterator[quad]->index()].push_back (cell);
	  
					   // make some checks on the lines
					   // and their ordering; if the
					   // lines are right, so are the
					   // vertices
	  Assert (face_iterator[0]->line(0) == face_iterator[2]->line(0),
		  ExcInternalErrorOnCell(c));
	  Assert (face_iterator[0]->line(1) == face_iterator[3]->line(3),
		  ExcInternalErrorOnCell(c));
	  Assert (face_iterator[0]->line(2) == face_iterator[4]->line(0),
		  ExcInternalErrorOnCell(c));
	  Assert (face_iterator[0]->line(3) == face_iterator[5]->line(3),
		  ExcInternalErrorOnCell(c));

	  Assert (face_iterator[1]->line(0) == face_iterator[2]->line(2),
		  ExcInternalErrorOnCell(c));
	  Assert (face_iterator[1]->line(1) == face_iterator[3]->line(1),
		  ExcInternalErrorOnCell(c));
	  Assert (face_iterator[1]->line(2) == face_iterator[4]->line(2),
		  ExcInternalErrorOnCell(c));
	  Assert (face_iterator[1]->line(3) == face_iterator[5]->line(1),
		  ExcInternalErrorOnCell(c));

	  Assert (face_iterator[2]->line(1) == face_iterator[3]->line(0),
		  ExcInternalErrorOnCell(c));
	  Assert (face_iterator[3]->line(2) == face_iterator[4]->line(1),
		  ExcInternalErrorOnCell(c));
	  Assert (face_iterator[4]->line(3) == face_iterator[5]->line(2),
		  ExcInternalErrorOnCell(c));
	  Assert (face_iterator[5]->line(0) == face_iterator[2]->line(3),
		  ExcInternalErrorOnCell(c));
	};
    };
  

				   /////////////////////////////////////////
				   // find those quads which are at the
				   // boundary and mark them appropriately
  for (quad_iterator quad=begin_quad(); quad!=end_quad(); ++quad)
    {
      const unsigned int n_adj_cells = adjacent_cells[quad->index()].size();
				       // assert that every quad has one or
				       // two adjacent cells
      AssertThrow ((n_adj_cells >= 1) &&
		   (n_adj_cells <= 2),
		   ExcInternalError());

				       // if only one cell: quad is at
				       // boundary -> give it the boundary
				       // indicator zero by default
      if (n_adj_cells == 1)
	quad->set_boundary_indicator (0);
      else
					 // interior quad -> 255
      	quad->set_boundary_indicator (255);
    };

				   /////////////////////////////////////////
				   // next find those lines which are at
				   // the boundary and mark all others as
				   // interior ones
				   //
				   // for this: first mark all lines as
				   // interior
  for (line_iterator line=begin_line(); line!=end_line(); ++line)
    line->set_boundary_indicator (255);
				   // next reset all lines bounding boundary
				   // quads as on the boundary also
  for (quad_iterator quad=begin_quad(); quad!=end_quad(); ++quad)
    if (quad->at_boundary())
      for (unsigned int l=0; l<4; ++l)
	quad->line(l)->set_boundary_indicator (0);
  

				   ///////////////////////////////////////
				   // now set boundary indicators
				   // where given
				   //
				   // not yet implemented
  AssertThrow ((subcelldata.boundary_lines.size() == 0) &&
	       (subcelldata.boundary_quads.size() == 0),
	       ExcInternalError());
  

				   /////////////////////////////////////////
				   // finally update neighborship info
  for (cell_iterator cell=begin(); cell!=end(); ++cell)
    for (unsigned int face=0; face<6; ++face)
      if (adjacent_cells[cell->quad(face)->index()][0] == cell)
					 // first adjacent cell is this one
	{
	  if (adjacent_cells[cell->quad(face)->index()].size() == 2)
					     // there is another adjacent cell
	    cell->set_neighbor (face,
				adjacent_cells[cell->quad(face)->index()][1]);
	}
				   // first adjacent cell is not this one,
				   // -> it must be the neighbor we are
				   // looking for
      else
	cell->set_neighbor (face,
			    adjacent_cells[cell->quad(face)->index()][0]);


				   // re-compute numbers of lines, etc
  update_number_cache ();  
};

#endif




#if deal_II_dimension == 1

template <>
void Triangulation<1>::distort_random (const double factor,
				       const bool   keep_boundary) {
				   // this function is mostly equivalent to
				   // that for the general dimensional case
				   // the only difference being the correction
				   // for split faces which is not necessary
				   // in 1D
				   //
				   // if you change something here, don't
				   // forget to do so there as well

  const unsigned int dim = 1;
  
				   // find the smallest length of the lines
				   // adjacent to the vertex
  vector<double>             minimal_length (vertices.size(), 1e308);
				   // also note if a vertex is at
				   // the boundary
  vector<bool>               at_boundary (vertices.size(), false);
  
  for (active_line_iterator line=begin_active_line();
       line != end_line(); ++line)
    {
      if (keep_boundary && line->at_boundary())
	{
	  at_boundary[line->vertex_index(0)] = true;
	  at_boundary[line->vertex_index(1)] = true;
	};
      
      minimal_length[line->vertex_index(0)]
	= min(line->diameter(), minimal_length[line->vertex_index(0)]);
      minimal_length[line->vertex_index(1)]
	= min(line->diameter(), minimal_length[line->vertex_index(1)]);
    };


  const unsigned int n_vertices = vertices.size();
  Point<dim> shift_vector;
  
  for (unsigned int vertex=0; vertex<n_vertices; ++vertex) 
    {
				       // ignore this vertex if we whall keep
				       // the boundary and this vertex *is* at
				       // the boundary
      if (keep_boundary && at_boundary[vertex])
	continue;
      
				       // first compute a random shift vector
      for (unsigned int d=0; d<dim; ++d)
	shift_vector(d) = rand()*1.0/RAND_MAX;

      shift_vector *= factor * minimal_length[vertex] /
		      sqrt(shift_vector.square());

				       // finally move the vertex
      vertices[vertex] += shift_vector;
    };
};

#endif



template <int dim>
void Triangulation<dim>::distort_random (const double factor,
					 const bool   keep_boundary) {
				   // this function is mostly equivalent to
				   // that for the general dimensional case
				   // the only difference being the correction
				   // for split faces which is not necessary
				   // in 1D
				   //
				   // if you change something here, don't
				   // forget to do so there as well
  
				   // find the smallest length of the lines
				   // adjecent to the vertex
  vector<double>             minimal_length (vertices.size(), 1e308);
				   // also note if a vertex is at
				   // the boundary
  vector<bool>               at_boundary (vertices.size(), false);
  
  for (active_line_iterator line=begin_active_line();
       line != end_line(); ++line)
    {
      if (keep_boundary && line->at_boundary())
	{
	  at_boundary[line->vertex_index(0)] = true;
	  at_boundary[line->vertex_index(1)] = true;
	};
      
      minimal_length[line->vertex_index(0)]
	= min(line->diameter(), minimal_length[line->vertex_index(0)]);
      minimal_length[line->vertex_index(1)]
	= min(line->diameter(), minimal_length[line->vertex_index(1)]);
    };


  const unsigned int n_vertices = vertices.size();
  Point<dim> shift_vector;
  
  for (unsigned int vertex=0; vertex<n_vertices; ++vertex) 
    {
				       // ignore this vertex if we whall keep
				       // the boundary and this vertex *is* at
				       // the boundary
      if (keep_boundary && at_boundary[vertex])
	continue;
      
				       // first compute a random shift vector
      for (unsigned int d=0; d<dim; ++d)
	shift_vector(d) = rand()*1.0/RAND_MAX;

      shift_vector *= factor * minimal_length[vertex] /
		      sqrt(shift_vector.square());

				       // finally move the vertex
      vertices[vertex] += shift_vector;
    };

  
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
					   // not implemented at present
					   // for dim=3 or higher
	  Assert (dim<=2, ExcInternalError());

					   // compute where the common
					   // point of the two child lines
					   // will lie and
					   // reset it to the correct value
	  vertices[cell->face(face)->child(0)->vertex_index(1)]
	    = (cell->face(face)->vertex(0) +
	       cell->face(face)->vertex(1)) / 2;
	}
};




template <int dim>
void Triangulation<dim>::set_all_refine_flags () {
  active_cell_iterator cell = begin_active(),
		       endc = end();

  for (; cell != endc; ++cell)
    {
      cell->clear_coarsen_flag();
      cell->set_refine_flag ();
    };
};



template <int dim>
void Triangulation<dim>::refine_global (const unsigned int times) {
  for (unsigned int i=0; i<times; ++i)
    {
      set_all_refine_flags();
      execute_refinement ();
    };
};



template <int dim>
void Triangulation<dim>::save_refine_flags (vector<bool> &v) const {
  v.resize (n_active_cells(), false);
  vector<bool>::iterator  i = v.begin();
  active_cell_iterator cell = begin_active(),
		       endc = end();
  for (; cell!=endc; ++cell, ++i)
    *i = cell->refine_flag_set();
};



template <int dim>
void Triangulation<dim>::save_refine_flags (ostream &out) const {
  vector<bool> v;
  save_refine_flags (v);
  write_bool_vector (mn_tria_refine_flags_begin, v, mn_tria_refine_flags_end,
		     out);
};



template <int dim>
void Triangulation<dim>::load_refine_flags (istream &in) {
  vector<bool> v;
  read_bool_vector (mn_tria_refine_flags_begin, v, mn_tria_refine_flags_end,
		    in);
  load_refine_flags (v);
};



template <int dim>
void Triangulation<dim>::load_refine_flags (const vector<bool> &v) {
  Assert (v.size() == n_active_cells(), ExcGridReadError());
  
  active_cell_iterator cell = begin_active(),
		       endc = end();
  vector<bool>::const_iterator i = v.begin();
  for (; cell!=endc; ++cell, ++i)
    if (*i == true)
      cell->set_refine_flag();
    else
      cell->clear_refine_flag();
};



template <int dim>
void Triangulation<dim>::save_coarsen_flags (vector<bool> &v) const {
  v.resize (n_active_cells(), false);
  vector<bool>::iterator  i = v.begin();
  active_cell_iterator cell = begin_active(),
		       endc = end();
  for (; cell!=endc; ++cell, ++i)
    *i = cell->coarsen_flag_set();
};



template <int dim>
void Triangulation<dim>::save_coarsen_flags (ostream &out) const {
  vector<bool> v;
  save_coarsen_flags (v);
  write_bool_vector (mn_tria_coarsen_flags_begin, v, mn_tria_coarsen_flags_end,
		     out);
};



template <int dim>
void Triangulation<dim>::load_coarsen_flags (istream &in) {
  vector<bool> v;
  read_bool_vector (mn_tria_coarsen_flags_begin, v, mn_tria_coarsen_flags_end,
		    in);
  load_coarsen_flags (v);
};



template <int dim>
void Triangulation<dim>::load_coarsen_flags (const vector<bool> &v) {
  Assert (v.size() == n_active_cells(), ExcGridReadError());
  
  active_cell_iterator cell = begin_active(),
		       endc = end();
  vector<bool>::const_iterator i = v.begin();
  for (; cell!=endc; ++cell, ++i)
    if (*i == true)
      cell->set_coarsen_flag();
    else
      cell->clear_coarsen_flag();
};



#if deal_II_dimension == 1

template <>
void Triangulation<1>::clear_user_pointers () {
  cell_iterator cell = begin(),
		endc = end();
  for (; cell!=endc; ++cell)
    cell->clear_user_pointer ();
};



template <>
void Triangulation<1>::clear_user_flags () {
  cell_iterator cell = begin(),
		endc = end();
  for (; cell!=endc; ++cell)
    cell->clear_user_flag ();
};

#endif


#if deal_II_dimension == 2

template <>
void Triangulation<2>::clear_user_pointers () {
  line_iterator line = begin_line(),
		endl = end_line();
  for (; line!=endl; ++line)
    line->clear_user_pointer ();

  cell_iterator cell = begin(),
		endc = end();
  for (; cell!=endc; ++cell)
    cell->clear_user_pointer ();
};



template <>
void Triangulation<2>::clear_user_flags () {
  line_iterator line = begin_line(),
		endl = end_line();
  for (; line!=endl; ++line)
    line->clear_user_flag ();

  cell_iterator cell = begin(),
		endc = end();
  for (; cell!=endc; ++cell)
    cell->clear_user_flag ();
};


#endif




#if deal_II_dimension == 3

template <>
void Triangulation<3>::clear_user_pointers () {
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
};



template <>
void Triangulation<3>::clear_user_flags () {
  line_iterator line = begin_line(),
		endl = end_line();
  for (; line!=endl; ++line)
    line->clear_user_flag ();

  quad_iterator quad = begin_quad(),
		endq = end_quad();
  for (; quad!=endq; ++quad)
    quad->clear_user_flag ();

  cell_iterator cell = begin(),
		endc = end();
  for (; cell!=endc; ++cell)
    cell->clear_user_flag ();
};


#endif




template <int dim>
void Triangulation<dim>::save_user_flags (ostream &out) const {
  save_user_flags_line (out);
  
  if (dim>=2)
    save_user_flags_quad (out);
  
  if (dim>=3)
    save_user_flags_hex (out);

  if (dim >= 4)
    Assert (false, ExcNotImplemented());
};



template <int dim>
void Triangulation<dim>::save_user_flags (vector<bool> &v) const {
				   // clear vector and append
				   // all the stuff later on
  v.clear ();

  vector<bool> tmp;

  save_user_flags_line (tmp);
  v.insert (v.end(), tmp.begin(), tmp.end());

  if (dim >= 2)
    {
      save_user_flags_quad (tmp);
      v.insert (v.end(), tmp.begin(), tmp.end());
    };
  
  if (dim >= 3)
    {
      save_user_flags_hex (tmp);
      v.insert (v.end(), tmp.begin(), tmp.end());
    };      

  if (dim >= 4)
    Assert (false, ExcNotImplemented());
};



template <int dim>
void Triangulation<dim>::load_user_flags (istream &in)
{
  load_user_flags_line (in);
  
  if (dim>=2)
    load_user_flags_quad (in);
  
  if (dim>=3)
    load_user_flags_hex (in);

  if (dim >= 4)
    Assert (false, ExcNotImplemented());
};



template <int dim>
void Triangulation<dim>::load_user_flags (const vector<bool> &v)
{
  Assert (v.size() == n_lines()+n_quads()+n_hexs(), ExcInternalError());
  vector<bool> tmp;

				   // first extract the flags belonging
				   // to lines
  Assert (v.size() > n_lines(), ExcInternalError());
  tmp.insert (tmp.end(),
	      v.begin(), v.begin()+n_lines());
				   // and set the lines
  load_user_flags_line (tmp);
  tmp.clear ();

  if (dim >= 2)
    {
      Assert (v.size() > n_lines()+n_quads(), ExcInternalError());
      tmp.insert (tmp.end(),
		  v.begin()+n_lines(), v.begin()+n_lines()+n_quads());
      load_user_flags_quad (tmp);
    };
  
  if (dim >= 3)
    {
      tmp.insert (tmp.end(),
		  v.begin()+n_lines()+n_quads(), v.begin()+n_lines()+n_quads()+n_hexs());
      load_user_flags_hex (tmp);
    };      

  if (dim >= 4)
    Assert (false, ExcNotImplemented());
};



template <int dim>
void Triangulation<dim>::save_user_flags_line (vector<bool> &v) const {
  v.resize (n_lines(), false);
  vector<bool>::iterator  i = v.begin();
  line_iterator line = begin_line(),
		endl = end_line();
  for (; line!=endl; ++line, ++i)
    *i = line->user_flag_set();
};



template <int dim>
void Triangulation<dim>::save_user_flags_line (ostream &out) const {
  vector<bool> v;
  save_user_flags_line (v);
  write_bool_vector (mn_tria_line_user_flags_begin, v, mn_tria_line_user_flags_end,
		     out);
};



template <int dim>
void Triangulation<dim>::load_user_flags_line (istream &in) {
  vector<bool> v;
  read_bool_vector (mn_tria_line_user_flags_begin, v, mn_tria_line_user_flags_end,
		    in);
  load_user_flags_line (v);
};



template <int dim>
void Triangulation<dim>::load_user_flags_line (const vector<bool> &v) {
  Assert (v.size() == n_lines(), ExcGridReadError());
  
  line_iterator line = begin_line(),
		endl = end_line();
  vector<bool>::const_iterator i = v.begin();
  for (; line!=endl; ++line, ++i)
    if (*i == true)
      line->set_user_flag();
    else
      line->clear_user_flag();
};



#if deal_II_dimension == 1

template <>
void Triangulation<1>::save_user_flags_quad (ostream &) const {
  Assert (false, ExcFunctionNotUseful());
};



template <>
void Triangulation<1>::save_user_flags_quad (vector<bool> &) const {
  Assert (false, ExcFunctionNotUseful());
};


template <>
void Triangulation<1>::load_user_flags_quad (istream &) {
  Assert (false, ExcFunctionNotUseful());
};


template <>
void Triangulation<1>::load_user_flags_quad (const vector<bool> &) {
  Assert (false, ExcFunctionNotUseful());
};


template <>
void Triangulation<1>::save_user_flags_hex (ostream &) const {
  Assert (false, ExcFunctionNotUseful());
};



template <>
void Triangulation<1>::save_user_flags_hex (vector<bool> &) const {
  Assert (false, ExcFunctionNotUseful());
};


template <>
void Triangulation<1>::load_user_flags_hex (istream &) {
  Assert (false, ExcFunctionNotUseful());
};


template <>
void Triangulation<1>::load_user_flags_hex (const vector<bool> &) {
  Assert (false, ExcFunctionNotUseful());
};

#endif




template <int dim>
void Triangulation<dim>::save_user_flags_quad (vector<bool> &v) const {
  v.resize (n_quads(), false);
  vector<bool>::iterator  i = v.begin();
  quad_iterator quad = begin_quad(),
		endq = end_quad();
  for (; quad!=endq; ++quad, ++i)
    *i = quad->user_flag_set();
};



template <int dim>
void Triangulation<dim>::save_user_flags_quad (ostream &out) const {
  vector<bool> v;
  save_user_flags_quad (v);
  write_bool_vector (mn_tria_quad_user_flags_begin, v, mn_tria_quad_user_flags_end,
		     out);
};



template <int dim>
void Triangulation<dim>::load_user_flags_quad (istream &in) {
  vector<bool> v;
  read_bool_vector (mn_tria_quad_user_flags_begin, v, mn_tria_quad_user_flags_end,
		    in);
  load_user_flags_quad (v);
};



template <int dim>
void Triangulation<dim>::load_user_flags_quad (const vector<bool> &v) {
  Assert (v.size() == n_quads(), ExcGridReadError());
  
  quad_iterator quad = begin_quad(),
		endq = end_quad();
  vector<bool>::const_iterator i = v.begin();
  for (; quad!=endq; ++quad, ++i)
    if (*i == true)
      quad->set_user_flag();
    else
      quad->clear_user_flag();
};



#if deal_II_dimension == 2

template <>
void Triangulation<2>::save_user_flags_hex (ostream &) const {
  Assert (false, ExcFunctionNotUseful());
};



template <>
void Triangulation<2>::save_user_flags_hex (vector<bool> &) const {
  Assert (false, ExcFunctionNotUseful());
};


template <>
void Triangulation<2>::load_user_flags_hex (istream &) {
  Assert (false, ExcFunctionNotUseful());
};


template <>
void Triangulation<2>::load_user_flags_hex (const vector<bool> &) {
  Assert (false, ExcFunctionNotUseful());
};


#endif




template <int dim>
void Triangulation<dim>::save_user_flags_hex (vector<bool> &v) const {
  v.resize (n_hexs(), false);
  vector<bool>::iterator  i = v.begin();
  hex_iterator hex = begin_hex(),
	      endh = end_hex();
  for (; hex!=endh; ++hex, ++i)
    *i = hex->user_flag_set();
};



template <int dim>
void Triangulation<dim>::save_user_flags_hex (ostream &out) const {
  vector<bool> v;
  save_user_flags_hex (v);
  write_bool_vector (mn_tria_hex_user_flags_begin, v, mn_tria_hex_user_flags_end,
		     out);
};



template <int dim>
void Triangulation<dim>::load_user_flags_hex (istream &in) {
  vector<bool> v;
  read_bool_vector (mn_tria_hex_user_flags_begin, v, mn_tria_hex_user_flags_end,
		    in);
  load_user_flags_hex (v);
};



template <int dim>
void Triangulation<dim>::load_user_flags_hex (const vector<bool> &v) {
  Assert (v.size() == n_hexs(), ExcGridReadError());
  
  hex_iterator hex = begin_hex(),
	      endh = end_hex();
  vector<bool>::const_iterator i = v.begin();
  for (; hex!=endh; ++hex, ++i)
    if (*i == true)
      hex->set_user_flag();
    else
      hex->clear_user_flag();
};



/*------------------------ Iterator functions ------------------------*/


#if deal_II_dimension == 1

template <>
TriaDimensionInfo<1>::raw_cell_iterator
Triangulation<1>::begin_raw (const unsigned int level) const {
  return begin_raw_line (level);
};



template <>
TriaDimensionInfo<1>::active_cell_iterator
Triangulation<1>::begin_active (const unsigned int level) const {
  return begin_active_line (level);
};



template <>
TriaDimensionInfo<1>::raw_cell_iterator
Triangulation<1>::last_raw () const {
  return last_raw_line ();
};



template <>
TriaDimensionInfo<1>::raw_cell_iterator
Triangulation<1>::last_raw (const unsigned int level) const {
  return last_raw_line (level);
};



template <>
TriaDimensionInfo<1>::cell_iterator
Triangulation<1>::last () const {
  return last_line ();
};



template <>
TriaDimensionInfo<1>::cell_iterator
Triangulation<1>::last (const unsigned int level) const {
  return last_line (level);
};



template <>
TriaDimensionInfo<1>::active_cell_iterator
Triangulation<1>::last_active () const {
  return last_active_line ();
};



template <>
TriaDimensionInfo<1>::active_cell_iterator
Triangulation<1>::last_active (const unsigned int level) const {
  return last_active_line (level);
};



template <>
TriaDimensionInfo<1>::raw_face_iterator
Triangulation<1>::begin_raw_face (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<1>::face_iterator
Triangulation<1>::begin_face (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<1>::active_face_iterator
Triangulation<1>::begin_active_face (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<1>::raw_face_iterator
Triangulation<1>::end_face () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<1>::raw_face_iterator
Triangulation<1>::last_raw_face () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<1>::raw_face_iterator
Triangulation<1>::last_raw_face (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<1>::face_iterator
Triangulation<1>::last_face () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<1>::face_iterator
Triangulation<1>::last_face (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<1>::active_face_iterator
Triangulation<1>::last_active_face () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<1>::active_face_iterator
Triangulation<1>::last_active_face (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};


template <>
TriaDimensionInfo<1>::raw_quad_iterator
Triangulation<1>::begin_raw_quad (unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<1>::quad_iterator
Triangulation<1>::begin_quad (unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<1>::active_quad_iterator
Triangulation<1>::begin_active_quad (unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<1>::raw_quad_iterator
Triangulation<1>::end_quad () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<1>::raw_quad_iterator
Triangulation<1>::last_raw_quad (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<1>::raw_quad_iterator
Triangulation<1>::last_raw_quad () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<1>::quad_iterator
Triangulation<1>::last_quad (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<1>::quad_iterator
Triangulation<1>::last_quad () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<1>::active_quad_iterator
Triangulation<1>::last_active_quad (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<1>::active_quad_iterator
Triangulation<1>::last_active_quad () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<1>::raw_hex_iterator
Triangulation<1>::begin_raw_hex (unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<1>::hex_iterator
Triangulation<1>::begin_hex (unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<1>::active_hex_iterator
Triangulation<1>::begin_active_hex (unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<1>::raw_hex_iterator
Triangulation<1>::end_hex () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<1>::raw_hex_iterator
Triangulation<1>::last_raw_hex (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<1>::raw_hex_iterator
Triangulation<1>::last_raw_hex () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<1>::hex_iterator
Triangulation<1>::last_hex (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<1>::hex_iterator
Triangulation<1>::last_hex () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<1>::active_hex_iterator
Triangulation<1>::last_active_hex (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<1>::active_hex_iterator
Triangulation<1>::last_active_hex () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};

#endif




#if deal_II_dimension == 2

template <>
TriaDimensionInfo<2>::raw_cell_iterator
Triangulation<2>::begin_raw (const unsigned int level) const {
  return begin_raw_quad (level);
};



template <>
TriaDimensionInfo<2>::active_cell_iterator
Triangulation<2>::begin_active (const unsigned int level) const {
  return begin_active_quad (level);
};



template <>
TriaDimensionInfo<2>::raw_cell_iterator
Triangulation<2>::last_raw () const {
  return last_raw_quad ();
};



template <>
TriaDimensionInfo<2>::raw_cell_iterator
Triangulation<2>::last_raw (const unsigned int level) const {
  return last_raw_quad (level);
};



template <>
TriaDimensionInfo<2>::cell_iterator
Triangulation<2>::last () const {
  return last_quad ();
};



template <>
TriaDimensionInfo<2>::cell_iterator
Triangulation<2>::last (const unsigned int level) const {
  return last_quad (level);
};



template <>
TriaDimensionInfo<2>::active_cell_iterator
Triangulation<2>::last_active () const {
  return last_active_quad ();
};



template <>
TriaDimensionInfo<2>::active_cell_iterator
Triangulation<2>::last_active (const unsigned int level) const {
  return last_active_quad (level);
};







template <>
TriaDimensionInfo<2>::raw_face_iterator
Triangulation<2>::begin_raw_face (const unsigned int level) const {
  return begin_raw_line (level);
};



template <>
TriaDimensionInfo<2>::face_iterator
Triangulation<2>::begin_face (const unsigned int level) const {
  return begin_line (level);
};



template <>
TriaDimensionInfo<2>::active_face_iterator
Triangulation<2>::begin_active_face (const unsigned int level) const {
  return begin_active_line (level);
};



template <>
TriaDimensionInfo<2>::raw_face_iterator
Triangulation<2>::end_face () const {
  return end_line ();
};



template <>
TriaDimensionInfo<2>::raw_face_iterator
Triangulation<2>::last_raw_face () const {
  return last_raw_line ();
};



template <>
TriaDimensionInfo<2>::raw_face_iterator
Triangulation<2>::last_raw_face (const unsigned int level) const {
  return last_raw_line (level);
};



template <>
TriaDimensionInfo<2>::face_iterator
Triangulation<2>::last_face () const {
  return last_line ();
};



template <>
TriaDimensionInfo<2>::face_iterator
Triangulation<2>::last_face (const unsigned int level) const {
  return last_line (level);
};



template <>
TriaDimensionInfo<2>::active_face_iterator
Triangulation<2>::last_active_face () const {
  return last_active_line ();
};



template <>
TriaDimensionInfo<2>::active_face_iterator
Triangulation<2>::last_active_face (const unsigned int level) const {
  return last_active_line (level);
};


template <>
TriaDimensionInfo<2>::raw_hex_iterator
Triangulation<2>::begin_raw_hex (unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<2>::hex_iterator
Triangulation<2>::begin_hex (unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<2>::active_hex_iterator
Triangulation<2>::begin_active_hex (unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<2>::raw_hex_iterator
Triangulation<2>::end_hex () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<2>::raw_hex_iterator
Triangulation<2>::last_raw_hex (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<2>::raw_hex_iterator
Triangulation<2>::last_raw_hex () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<2>::hex_iterator
Triangulation<2>::last_hex (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<2>::hex_iterator
Triangulation<2>::last_hex () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<2>::active_hex_iterator
Triangulation<2>::last_active_hex (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
TriaDimensionInfo<2>::active_hex_iterator
Triangulation<2>::last_active_hex () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};


#endif



#if deal_II_dimension == 3

template <>
TriaDimensionInfo<3>::raw_cell_iterator
Triangulation<3>::begin_raw (const unsigned int level) const {
  return begin_raw_hex (level);
};



template <>
TriaDimensionInfo<3>::active_cell_iterator
Triangulation<3>::begin_active (const unsigned int level) const {
  return begin_active_hex (level);
};



template <>
TriaDimensionInfo<3>::raw_cell_iterator
Triangulation<3>::last_raw () const {
  return last_raw_hex ();
};



template <>
TriaDimensionInfo<3>::raw_cell_iterator
Triangulation<3>::last_raw (const unsigned int level) const {
  return last_raw_hex (level);
};



template <>
TriaDimensionInfo<3>::cell_iterator
Triangulation<3>::last () const {
  return last_hex ();
};



template <>
TriaDimensionInfo<3>::cell_iterator
Triangulation<3>::last (const unsigned int level) const {
  return last_hex (level);
};



template <>
TriaDimensionInfo<3>::active_cell_iterator
Triangulation<3>::last_active () const {
  return last_active_hex ();
};



template <>
TriaDimensionInfo<3>::active_cell_iterator
Triangulation<3>::last_active (const unsigned int level) const {
  return last_active_hex (level);
};







template <>
TriaDimensionInfo<3>::raw_face_iterator
Triangulation<3>::begin_raw_face (const unsigned int level) const {
  return begin_raw_quad (level);
};



template <>
TriaDimensionInfo<3>::face_iterator
Triangulation<3>::begin_face (const unsigned int level) const {
  return begin_quad (level);
};



template <>
TriaDimensionInfo<3>::active_face_iterator
Triangulation<3>::begin_active_face (const unsigned int level) const {
  return begin_active_quad (level);
};



template <>
TriaDimensionInfo<3>::raw_face_iterator
Triangulation<3>::end_face () const {
  return end_quad ();
};



template <>
TriaDimensionInfo<3>::raw_face_iterator
Triangulation<3>::last_raw_face () const {
  return last_raw_quad ();
};



template <>
TriaDimensionInfo<3>::raw_face_iterator
Triangulation<3>::last_raw_face (const unsigned int level) const {
  return last_raw_quad (level);
};



template <>
TriaDimensionInfo<3>::face_iterator
Triangulation<3>::last_face () const {
  return last_quad ();
};



template <>
TriaDimensionInfo<3>::face_iterator
Triangulation<3>::last_face (const unsigned int level) const {
  return last_quad (level);
};



template <>
TriaDimensionInfo<3>::active_face_iterator
Triangulation<3>::last_active_face () const {
  return last_active_quad ();
};



template <>
TriaDimensionInfo<3>::active_face_iterator
Triangulation<3>::last_active_face (const unsigned int level) const {
  return last_active_quad (level);
};



#endif




template <int dim>
typename TriaDimensionInfo<dim>::raw_line_iterator
Triangulation<dim>::begin_raw_line (unsigned int level) const {
  Assert (level<levels.size(), ExcInvalidLevel(level));
  
  if (levels[level]->lines.lines.size() == 0)
    return end_line ();
  
  return raw_line_iterator (const_cast<Triangulation<dim>*>(this),
			    level,
			    0);
};




template <int dim>
typename TriaDimensionInfo<dim>::raw_quad_iterator
Triangulation<dim>::begin_raw_quad (unsigned int level) const {
  Assert (level<levels.size(), ExcInvalidLevel(level));

  if (levels[level]->quads.quads.size() == 0)
    return end_quad();
  
  return raw_quad_iterator (const_cast<Triangulation<dim>*>(this),
			    level,
			    0);
};



template <int dim>
typename TriaDimensionInfo<dim>::raw_hex_iterator
Triangulation<dim>::begin_raw_hex (unsigned int level) const {
  Assert (level<levels.size(), ExcInvalidLevel(level));

  if (levels[level]->hexes.hexes.size() == 0)
    return end_hex();
  
  return raw_hex_iterator (const_cast<Triangulation<dim>*>(this),
			   level,
			   0);
};




template <int dim>
typename TriaDimensionInfo<dim>::line_iterator
Triangulation<dim>::begin_line (unsigned int level) const {
  				   // level is checked in begin_raw
  raw_line_iterator ri = begin_raw_line (level);
  if (ri.state() != valid)
    return ri;
  while (ri->used() == false)
    if ((++ri).state() != valid)
      return ri;
  return ri;
};




template <int dim>
typename TriaDimensionInfo<dim>::quad_iterator
Triangulation<dim>::begin_quad (unsigned int level) const {
  				   // level is checked in begin_raw
  raw_quad_iterator ri = begin_raw_quad (level);
  if (ri.state() != valid)
    return ri;
  while (ri->used() == false)
    if ((++ri).state() != valid)
      return ri;
  return ri;
};



template <int dim>
typename TriaDimensionInfo<dim>::hex_iterator
Triangulation<dim>::begin_hex (unsigned int level) const {
  				   // level is checked in begin_raw
  raw_hex_iterator ri = begin_raw_hex (level);
  if (ri.state() != valid)
    return ri;
  while (ri->used() == false)
    if ((++ri).state() != valid)
      return ri;
  return ri;
};




template <int dim>
typename TriaDimensionInfo<dim>::active_line_iterator
Triangulation<dim>::begin_active_line (unsigned int level) const {
  				   // level is checked in begin_raw
  line_iterator i = begin_line (level);
  if (i.state() != valid)
    return i;
  while (i->has_children())
    if ((++i).state() != valid)
      return i;
  return i;
};



template <int dim>
typename TriaDimensionInfo<dim>::active_quad_iterator
Triangulation<dim>::begin_active_quad (unsigned int level) const {
  				   // level is checked in begin_raw
  quad_iterator i = begin_quad (level);
  if (i.state() != valid)
    return i;
  while (i->has_children())
    if ((++i).state() != valid)
      return i;
  return i;
};



template <int dim>
typename TriaDimensionInfo<dim>::active_hex_iterator
Triangulation<dim>::begin_active_hex (unsigned int level) const {
  				   // level is checked in begin_raw
  hex_iterator i = begin_hex (level);
  if (i.state() != valid)
    return i;
  while (i->has_children())
    if ((++i).state() != valid)
      return i;
  return i;
};



template <int dim>
typename TriaDimensionInfo<dim>::raw_line_iterator
Triangulation<dim>::end_line () const {
  return raw_line_iterator (const_cast<Triangulation<dim>*>(this),
			    -1,
			    -1);
};




template <int dim>
typename TriaDimensionInfo<dim>::raw_quad_iterator
Triangulation<dim>::end_quad () const {
  return raw_quad_iterator (const_cast<Triangulation<dim>*>(this),
			    -1,
			    -1);
};



template <int dim>
typename TriaDimensionInfo<dim>::raw_hex_iterator
Triangulation<dim>::end_hex () const {
  return raw_hex_iterator (const_cast<Triangulation<dim>*>(this),
			   -1,
			   -1);
};



template <int dim>
typename TriaDimensionInfo<dim>::raw_line_iterator
Triangulation<dim>::last_raw_line (const unsigned int level) const {
  Assert (level<levels.size(), ExcInvalidLevel(level));
  Assert (levels[level]->lines.lines.size() != 0,
	  ExcEmptyLevel (level));
  
  return raw_line_iterator (const_cast<Triangulation<dim>*>(this),
			    level,
			    levels[level]->lines.lines.size()-1);
};



template <int dim>
typename TriaDimensionInfo<dim>::raw_quad_iterator
Triangulation<dim>::last_raw_quad (const unsigned int level) const {
  Assert (level<levels.size(),
	  ExcInvalidLevel(level));
  Assert (levels[level]->quads.quads.size() != 0,
	  ExcEmptyLevel (level));

  return raw_quad_iterator (const_cast<Triangulation<dim>*>(this),
			    level,
			    levels[level]->quads.quads.size()-1);
};



template <int dim>
typename TriaDimensionInfo<dim>::raw_hex_iterator
Triangulation<dim>::last_raw_hex (const unsigned int level) const {
  Assert (level<levels.size(),
	  ExcInvalidLevel(level));
  Assert (levels[level]->hexes.hexes.size() != 0,
	  ExcEmptyLevel (level));

  return raw_hex_iterator (const_cast<Triangulation<dim>*>(this),
			   level,
			   levels[level]->hexes.hexes.size()-1);
};



template <int dim>
typename TriaDimensionInfo<dim>::raw_line_iterator
Triangulation<dim>::last_raw_line () const {
  return last_raw_line (levels.size()-1);
};




template <int dim>
typename TriaDimensionInfo<dim>::raw_quad_iterator
Triangulation<dim>::last_raw_quad () const {
  return last_raw_quad (levels.size()-1);
};



template <int dim>
typename TriaDimensionInfo<dim>::raw_hex_iterator
Triangulation<dim>::last_raw_hex () const {
  return last_raw_hex (levels.size()-1);
};



template <int dim>
typename TriaDimensionInfo<dim>::line_iterator
Triangulation<dim>::last_line (const unsigned int level) const {
  				   // level is checked in begin_raw
  raw_line_iterator ri = last_raw_line(level);
  if (ri->used()==true)
    return ri;
  while ((--ri).state() == valid)
    if (ri->used()==true)
      return ri;
  return ri;
};



template <int dim>
typename TriaDimensionInfo<dim>::quad_iterator
Triangulation<dim>::last_quad (const unsigned int level) const {
  				   // level is checked in begin_raw
  raw_quad_iterator ri = last_raw_quad(level);
  if (ri->used()==true)
    return ri;
  while ((--ri).state() == valid)
    if (ri->used()==true)
      return ri;
  return ri;
};



template <int dim>
typename TriaDimensionInfo<dim>::hex_iterator
Triangulation<dim>::last_hex (const unsigned int level) const {
  				   // level is checked in begin_raw
  raw_hex_iterator ri = last_raw_hex(level);
  if (ri->used()==true)
    return ri;
  while ((--ri).state() == valid)
    if (ri->used()==true)
      return ri;
  return ri;
};



template <int dim>
typename TriaDimensionInfo<dim>::line_iterator
Triangulation<dim>::last_line () const {
  return last_line (levels.size()-1);
};



template <int dim>
typename TriaDimensionInfo<dim>::quad_iterator
Triangulation<dim>::last_quad () const {
  return last_quad (levels.size()-1);
};



template <int dim>
typename TriaDimensionInfo<dim>::hex_iterator
Triangulation<dim>::last_hex () const {
  return last_hex (levels.size()-1);
};




template <int dim>
typename TriaDimensionInfo<dim>::active_line_iterator
Triangulation<dim>::last_active_line (const unsigned int level) const {
				   // level is checked in begin_raw
  line_iterator i = last_line(level);
  if (i->has_children()==false)
    return i;
  while ((--i).state() == valid)
    if (i->has_children()==false)
      return i;
  return i;
};



template <int dim>
typename TriaDimensionInfo<dim>::active_quad_iterator
Triangulation<dim>::last_active_quad (const unsigned int level) const {
				   // level is checked in begin_raw
  quad_iterator i = last_quad(level);
  if (i->has_children()==false)
    return i;
  while ((--i).state() == valid)
    if (i->has_children()==false)
      return i;
  return i;
};



template <int dim>
typename TriaDimensionInfo<dim>::active_hex_iterator
Triangulation<dim>::last_active_hex (const unsigned int level) const {
				   // level is checked in begin_raw
  hex_iterator i = last_hex(level);
  if (i->has_children()==false)
    return i;
  while ((--i).state() == valid)
    if (i->has_children()==false)
      return i;
  return i;
};




template <int dim>
typename TriaDimensionInfo<dim>::active_line_iterator
Triangulation<dim>::last_active_line () const {
  return last_active_line (levels.size()-1);
};




template <int dim>
typename TriaDimensionInfo<dim>::active_quad_iterator
Triangulation<dim>::last_active_quad () const {
  return last_active_quad (levels.size()-1);
};



template <int dim>
typename TriaDimensionInfo<dim>::active_hex_iterator
Triangulation<dim>::last_active_hex () const {
  return last_active_hex (levels.size()-1);
};




template <int dim>
typename TriaDimensionInfo<dim>::raw_cell_iterator
Triangulation<dim>::end_raw (const unsigned int level) const {
  return (level == levels.size()-1 ?
	  end() :
	  begin_raw (level+1));
};



template <int dim>
typename TriaDimensionInfo<dim>::cell_iterator
Triangulation<dim>::end (const unsigned int level) const {
  return (level == levels.size()-1 ?
	  cell_iterator(end()) :
	  begin (level+1));
};



template <int dim>
typename TriaDimensionInfo<dim>::active_cell_iterator
Triangulation<dim>::end_active (const unsigned int level) const {
  return (level == levels.size()-1 ?
	  active_cell_iterator(end()) :
	  begin_active (level+1));
};



template <int dim>
typename TriaDimensionInfo<dim>::raw_face_iterator
Triangulation<dim>::end_raw_face (const unsigned int level) const {
  return (level == levels.size()-1 ?
	  end_face() :
	  begin_raw_face (level+1));
};



template <int dim>
typename TriaDimensionInfo<dim>::face_iterator
Triangulation<dim>::end_face (const unsigned int level) const {
  return (level == levels.size()-1 ?
	  face_iterator(end_face()) :
	  begin_face (level+1));
};



template <int dim>
typename TriaDimensionInfo<dim>::active_face_iterator
Triangulation<dim>::end_active_face (const unsigned int level) const {
  return (level == levels.size()-1 ?
	  active_face_iterator(end_face()) :
	  begin_active_face (level+1));
};



template <int dim>
typename TriaDimensionInfo<dim>::raw_line_iterator
Triangulation<dim>::end_raw_line (const unsigned int level) const {
  return (level == levels.size()-1 ?
	  end_line() :
	  begin_raw_line (level+1));
};



template <int dim>
typename TriaDimensionInfo<dim>::line_iterator
Triangulation<dim>::end_line (const unsigned int level) const {
  return (level == levels.size()-1 ?
	  line_iterator(end_line()) :
	  begin_line (level+1));
};



template <int dim>
typename TriaDimensionInfo<dim>::active_line_iterator
Triangulation<dim>::end_active_line (const unsigned int level) const {
  return (level == levels.size()-1 ?
	  active_line_iterator(end_line()) :
	  begin_active_line (level+1));
};




template <int dim>
typename TriaDimensionInfo<dim>::raw_quad_iterator
Triangulation<dim>::end_raw_quad (const unsigned int level) const {
  return (level == levels.size()-1 ?
	  end_quad() :
	  begin_raw_quad (level+1));
};



template <int dim>
typename TriaDimensionInfo<dim>::quad_iterator
Triangulation<dim>::end_quad (const unsigned int level) const {
  return (level == levels.size()-1 ?
	  quad_iterator(end_quad()) :
	  begin_quad (level+1));
};



template <int dim>
typename TriaDimensionInfo<dim>::active_quad_iterator
Triangulation<dim>::end_active_quad (const unsigned int level) const {
  return (level == levels.size()-1 ?
	  active_quad_iterator(end_quad()) :
	  begin_active_quad (level+1));
};





template <int dim>
typename TriaDimensionInfo<dim>::raw_hex_iterator
Triangulation<dim>::end_raw_hex (const unsigned int level) const {
  return (level == levels.size()-1 ?
	  end_hex() :
	  begin_raw_hex (level+1));
};



template <int dim>
typename TriaDimensionInfo<dim>::hex_iterator
Triangulation<dim>::end_hex (const unsigned int level) const {
  return (level == levels.size()-1 ?
	  hex_iterator(end_hex()) :
	  begin_hex (level+1));
};



template <int dim>
typename TriaDimensionInfo<dim>::active_hex_iterator
Triangulation<dim>::end_active_hex (const unsigned int level) const {
  return (level == levels.size()-1 ?
	  active_hex_iterator(end_hex()) :
	  begin_active_hex (level+1));
};





template <int dim>
unsigned int Triangulation<dim>::n_cells () const {
  unsigned int n=0;
  for (unsigned int l=0; l<levels.size(); ++l)
    n += n_cells (l);
  return n;
};



template <int dim>
unsigned int Triangulation<dim>::n_active_cells () const {
  unsigned int n=0;
  for (unsigned int l=0; l<levels.size(); ++l)
    n += n_active_cells (l);
  return n;
};



#if deal_II_dimension == 1

template <>
unsigned int Triangulation<1>::n_active_cells (const unsigned int level) const {
  return n_active_lines (level);
};


template <>
unsigned int Triangulation<1>::n_cells (const unsigned int level) const {
  return n_lines (level);
};

#endif



#if deal_II_dimension == 2

template <>
unsigned int Triangulation<2>::n_cells (const unsigned int level) const {
  return n_quads (level);
};



template <>
unsigned int Triangulation<2>::n_active_cells (const unsigned int level) const {
  return n_active_quads (level);
};

#endif



#if deal_II_dimension == 3

template <>
unsigned int Triangulation<3>::n_cells (const unsigned int level) const {
  return n_hexs (level);
};



template <>
unsigned int Triangulation<3>::n_active_cells (const unsigned int level) const {
  return n_active_hexs (level);
};

#endif




template <int dim>
unsigned int Triangulation<dim>::n_lines () const {
  return number_cache.n_lines;
};



template <int dim>
unsigned int Triangulation<dim>::n_lines (const unsigned int level) const {
  Assert (level < number_cache.n_lines_level.size(),
	  ExcIndexRange (level, 0, number_cache.n_lines_level.size()));
  
  return number_cache.n_lines_level[level];
};



template <int dim>
unsigned int Triangulation<dim>::n_active_lines () const {
  return number_cache.n_active_lines;
};



template <int dim>
unsigned int Triangulation<dim>::n_active_lines (const unsigned int level) const {
  Assert (level < number_cache.n_lines_level.size(),
	  ExcIndexRange (level, 0, number_cache.n_lines_level.size()));
  
  return number_cache.n_active_lines_level[level];
};



#if deal_II_dimension == 1

template <>
unsigned int Triangulation<1>::n_quads () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
unsigned int Triangulation<1>::n_quads (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
unsigned int Triangulation<1>::n_active_quads (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};


template <>
unsigned int Triangulation<1>::n_active_quads () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};


#endif


template <int dim>
unsigned int Triangulation<dim>::n_quads () const {
  return number_cache.n_quads;
};



template <int dim>
unsigned int Triangulation<dim>::n_quads (const unsigned int level) const {
  Assert (level < number_cache.n_quads_level.size(),
	  ExcIndexRange (level, 0, number_cache.n_quads_level.size()));
  
  return number_cache.n_quads_level[level];
};



template <int dim>
unsigned int Triangulation<dim>::n_active_quads () const {
  return number_cache.n_active_quads;
};



template <int dim>
unsigned int Triangulation<dim>::n_active_quads (const unsigned int level) const {
  Assert (level < number_cache.n_quads_level.size(),
	  ExcIndexRange (level, 0, number_cache.n_quads_level.size()));
  
  return number_cache.n_active_quads_level[level];
};



#if deal_II_dimension < 3

template <>
unsigned int Triangulation<1>::n_hexs () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
unsigned int Triangulation<1>::n_hexs (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
unsigned int Triangulation<1>::n_active_hexs () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
unsigned int Triangulation<1>::n_active_hexs (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
unsigned int Triangulation<2>::n_hexs () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
unsigned int Triangulation<2>::n_hexs (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
unsigned int Triangulation<2>::n_active_hexs () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <>
unsigned int Triangulation<2>::n_active_hexs (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};

#endif


template <int dim>
unsigned int Triangulation<dim>::n_hexs () const {
  return number_cache.n_hexes;
};



template <int dim>
unsigned int Triangulation<dim>::n_hexs (const unsigned int level) const {
  Assert (level < number_cache.n_hexes_level.size(),
	  ExcIndexRange (level, 0, number_cache.n_hexes_level.size()));
  
  return number_cache.n_hexes_level[level];
};



template <int dim>
unsigned int Triangulation<dim>::n_active_hexs () const {
  return number_cache.n_active_hexes;
};



template <int dim>
unsigned int Triangulation<dim>::n_active_hexs (const unsigned int level) const {
  Assert (level < number_cache.n_hexes_level.size(),
	  ExcIndexRange (level, 0, number_cache.n_hexes_level.size()));
  
  return number_cache.n_active_hexes_level[level];
};




template <int dim>
unsigned int Triangulation<dim>::n_levels () const {
  if (levels.size() == 0)
    return 0;
				   // check whether there are
				   // cells on the highest levels
				   // (there need not be, since they
				   // might all have been coarsened
				   // away)
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
};



template <int dim>
unsigned int
Triangulation<dim>::n_vertices () const 
{
  return vertices.size();
};



template <int dim>
unsigned int
Triangulation<dim>::n_used_vertices () const 
{
  return count_if (vertices_used.begin(), vertices_used.end(),
		   bind2nd (equal_to<bool>(), true));
};



#if deal_II_dimension == 1

template <>
unsigned int Triangulation<1>::max_adjacent_cells () const {
  return 2;
};

#endif



template <int dim>
unsigned int Triangulation<dim>::max_adjacent_cells () const {
  cell_iterator cell = begin(0),
		endc = (levels.size() > 1 ? begin(1) : cell_iterator(end()));
				   // store the largest index of the vertices
				   // used on level 0
  unsigned int max_vertex_index = 0;
  for (; cell!=endc; ++cell)
    for (unsigned vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
      if (cell->vertex_index(vertex) > (signed int)max_vertex_index)
	max_vertex_index = cell->vertex_index(vertex);

				   // store the number of times a cell touches
				   // a vertex. An unsigned int should suffice,
				   // even for larger dimensions
  vector<unsigned short int> usage_count (max_vertex_index+1, 0);
				   // touch a vertex's usage count everytime
				   // we find an adjacent element
  for (cell=begin(); cell!=endc; ++cell)
    for (unsigned vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
      ++usage_count[cell->vertex_index(vertex)];

  return max (GeometryInfo<dim>::vertices_per_cell,
	      static_cast<unsigned int>(*(max_element (usage_count.begin(),
						       usage_count.end()))));
};





template <int dim>
template <typename number>
void Triangulation<dim>::refine (const Vector<number> &criteria,
				 const double         threshold)
{
  Assert (criteria.size() == n_active_cells(),
	  ExcInvalidVectorSize(criteria.size(), n_active_cells()));

				   // nothing to do; especially we
				   // do not want to flag with zero
				   // error since then we may get
				   // into conflict with coarsening
				   // in some cases
  if (threshold==0)
    return;
  
  active_cell_iterator cell = begin_active();
  const unsigned int n_cells = criteria.size();
  
  for (unsigned int index=0; index<n_cells; ++cell, ++index)
    if (fabs(criteria(index)) >= threshold)
      cell->set_refine_flag();
};



template <int dim>
template <typename number>
void Triangulation<dim>::coarsen (const Vector<number> &criteria,
				  const double         threshold)
{
  Assert (criteria.size() == n_active_cells(),
	  ExcInvalidVectorSize(criteria.size(), n_active_cells()));

  active_cell_iterator cell = begin_active();
  const unsigned int n_cells = criteria.size();
  
  for (unsigned int index=0; index<n_cells; ++cell, ++index)
    if (fabs(criteria(index)) <= threshold)
      cell->set_coarsen_flag();
};



template <int dim>
template <typename number>
void
Triangulation<dim>::refine_and_coarsen_fixed_number (const Vector<number> &criteria,
						     const double         top_fraction,
						     const double         bottom_fraction)
{
				   // correct number of cells is
				   // checked in #refine#
  Assert ((top_fraction>=0) && (top_fraction<=1), ExcInvalidParameterValue());
  Assert ((bottom_fraction>=0) && (bottom_fraction<=1), ExcInvalidParameterValue());
  Assert (top_fraction+bottom_fraction <= 1, ExcInvalidParameterValue());

  const int refine_cells=static_cast<int>(top_fraction*criteria.size());
  const int coarsen_cells=static_cast<int>(bottom_fraction*criteria.size());

  if (refine_cells || coarsen_cells)
    {
      Vector<number> tmp(criteria);
      if (refine_cells)
	{
	  nth_element (tmp.begin(), tmp.begin()+refine_cells,
		       tmp.end(),
		       greater<double>());
	  refine (criteria, *(tmp.begin() + refine_cells));
	}

      if (coarsen_cells)
	{
	  nth_element (tmp.begin(), tmp.begin()+tmp.size()-coarsen_cells,
		       tmp.end(),
		       greater<double>());
	  coarsen (criteria, *(tmp.begin() + tmp.size() - coarsen_cells));
	}
    }
};



static
inline
double sqr(double a) {
  return a*a;
};



template <int dim>
template <typename number>
void
Triangulation<dim>::refine_and_coarsen_fixed_fraction (const Vector<number> &criteria,
						       const double         top_fraction,
						       const double         bottom_fraction) {
				   // correct number of cells is
				   // checked in #refine#
  Assert ((top_fraction>=0) && (top_fraction<=1), ExcInvalidParameterValue());
  Assert ((bottom_fraction>=0) && (bottom_fraction<=1), ExcInvalidParameterValue());
  Assert (top_fraction+bottom_fraction <= 1, ExcInvalidParameterValue());

				   // let tmp be the cellwise square of the
				   // error, which is what we have to sum
				   // up and compare with
				   // #fraction_of_error*total_error#.
  Vector<number> tmp(criteria);
  const double total_error = tmp.l1_norm();

  Vector<float> partial_sums(criteria.size());
  
				   // sort the largest criteria to the
				   // beginning of the vector
  sort (tmp.begin(), tmp.end(), greater<double>());
  partial_sum (tmp.begin(), tmp.end(), partial_sums.begin());

				   // compute thresholds
  const Vector<float>::const_iterator
    q = lower_bound (partial_sums.begin(), partial_sums.end(),
		     top_fraction*total_error),
    p = upper_bound (partial_sums.begin(), partial_sums.end(),
		     total_error*(1-bottom_fraction));
  
  double bottom_threshold = tmp(p != partial_sums.end() ?
				p-partial_sums.begin() :
				criteria.size()-1),
	 top_threshold    = tmp(q-partial_sums.begin());

				   // in some rare cases it may happen that
				   // both thresholds are the same (e.g. if
				   // there are many cells with the same
				   // error indicator). That would mean that
				   // all cells will be flagged for
				   // refinement or coarsening, but some will
				   // be flagged for both, namely those for
				   // which the indicator equals the
				   // thresholds. This is forbidden, however.
				   //
				   // In some rare cases with very few cells
				   // we also could get integer round off
				   // errors and get problems with
				   // the top and bottom fractions.
				   //
				   // In these case we arbitrarily reduce the
				   // bottom threshold by one permille below
				   // the top threshold
  if (bottom_threshold>=top_threshold)
    bottom_threshold = 0.999*top_threshold;
  
				   // actually flag cells
  if (top_threshold < *max_element(criteria.begin(), criteria.end()))
    refine (criteria, top_threshold);
  if (bottom_threshold > *min_element(criteria.begin(), criteria.end()))
    coarsen (criteria, bottom_threshold);

  prepare_coarsening_and_refinement ();
};





template <int dim>
void Triangulation<dim>::execute_coarsening_and_refinement () {
  prepare_coarsening_and_refinement ();

  execute_coarsening();
  execute_refinement();
};




#if deal_II_dimension == 1

template <>
void Triangulation<1>::execute_refinement () {
  const unsigned int dim = 1;
  
				   // check whether a new level is needed
				   // we have to check for this on the
				   // highest level only (on this, all
				   // used cells are also active, so we
				   // only have to check for this)
  if (true)
    {
      raw_cell_iterator cell = begin_active (levels.size()-1),
			endc = end();
      for (; cell != endc; ++cell)
	if (cell->used())
	  if (cell->refine_flag_set()) 
	    {
	      levels.push_back (new TriangulationLevel<dim>);
	      break;
	    };
    };
  



				   // check how much space is needed
				   // on every level
				   // we need not check the highest
				   // level since either
				   // - on the highest level no cells
				   //   are flagged for refinement
				   // - there are, but prepare_refinement
				   //   added another empty level
  unsigned int needed_vertices = 0;
  for (int level=levels.size()-2; level>=0; --level)
    {
				       // count number of flagged cells on
				       // this level
      unsigned int flagged_cells = 0;
      active_cell_iterator acell = begin_active(level),
			   aendc = begin_active(level+1);
      for (; acell!=aendc; ++acell) 
	if (acell->refine_flag_set())
	  ++flagged_cells;

				       // count number of used cells on
				       // the next higher level
      const unsigned int used_cells
	=  count_if (levels[level+1]->lines.used.begin(),
		     levels[level+1]->lines.used.end(),
		     bind2nd (equal_to<bool>(), true));

				       // reserve space for the used_cells
				       // cells already existing on the next
				       // higher level as well as for the
				       // 2*flagged_cells that will be created
				       // on that level
      levels[level+1]->
	TriangulationLevel<0>::reserve_space (used_cells+
					      GeometryInfo<1>::children_per_cell *
					      flagged_cells, 1);
				       // reserve space for 2*flagged_cells
				       // new lines on the next higher
				       // level
      levels[level+1]->
	TriangulationLevel<1>::reserve_space (GeometryInfo<1>::children_per_cell*flagged_cells);
      
      needed_vertices += flagged_cells;
    };

				   // add to needed vertices how
				   // many vertices are already in use
  needed_vertices += count_if (vertices_used.begin(), vertices_used.end(),
			       bind2nd (equal_to<bool>(), true));
				   // if we need more vertices: create them,
				   // if not: leave the array as is, since
				   // shrinking is not really possible because
				   // some of the vertices at the end may be
				   // in use
  if (needed_vertices > vertices.size())
    {
      vertices.resize (needed_vertices, Point<1>());
      vertices_used.resize (needed_vertices, false);
    };


				   // Do REFINEMENT
				   // on every level; exclude highest level as
				   // above

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

					     // search for next unused vertex      
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

					     // search for next two unused cell
					     // (++ takes care of the end of
					     // the vector)
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
	    first_child->set (Line (cell->vertex_index(0), next_unused_vertex));
	    first_child->set_material_id (cell->material_id());
	    
					     // reset neighborship info
					     // (refer to \Ref{TriangulationLevel<0>}
					     // for details)
	    first_child->set_neighbor (1, second_child);
	    if (cell->neighbor(0).state() != valid)
	      first_child->set_neighbor (0, cell->neighbor(0));
	    else
	      if (cell->neighbor(0)->active())
		{
						   // since the neighbors level is
						   // always <=level, if the
						   // cell is active, then there
						   // are no cells to the left which
						   // may want to know about this
						   // new child cell.
		  Assert (cell->neighbor(0)->level() <= cell->level(),
			  ExcInternalError());
		  first_child->set_neighbor (0, cell->neighbor(0));
		}
	      else
						 // left neighbor is refined
		{
						   // set neighbor to cell on
						   // same level
		  first_child->set_neighbor (0, cell->neighbor(0)->child(1));

						   // reset neighbor info of
						   // all right descendant of the
						   // left neighbor of cell
		  cell_iterator left_neighbor = cell->neighbor(0);
		  while (left_neighbor->has_children())
		    {
		      left_neighbor = left_neighbor->child(1);
		      left_neighbor->set_neighbor (1, first_child);
		    };
		};
	    
					     // insert second child
	    second_child->clear_children ();
	    second_child->set (Line (next_unused_vertex, cell->vertex_index(1)));
	    second_child->set_neighbor (0, first_child);
	    second_child->set_material_id (cell->material_id());
	    if (cell->neighbor(1).state() != valid)
	      second_child->set_neighbor (1, cell->neighbor(1));
	    else
	      if (cell->neighbor(1)->active())
		{
		  Assert (cell->neighbor(1)->level() <= cell->level(),
			  ExcInternalError());
		  second_child->set_neighbor (1, cell->neighbor(1));
		}
	      else
						 // right neighbor is refined
						 // same as above
		{
		  second_child->set_neighbor (1, cell->neighbor(1)->child(0));
		  
		  cell_iterator right_neighbor = cell->neighbor(1);
		  while (right_neighbor->has_children())
		    {
		      right_neighbor = right_neighbor->child(0);
		      right_neighbor->set_neighbor (0, second_child);
		    };
		};
	  };      
    };


				   // re-compute number of lines
  update_number_cache ();
  
#ifdef DEBUG
  for (unsigned int level=0; level<levels.size(); ++level) 
    levels[level]->monitor_memory (1);

				   // check whether really all refinement flags
				   // are reset (also of previously non-active
				   // cells which we may not have touched. If
				   // the refinement flag of a non-active cell
				   // is set, something went wrong since the
				   // cell-accessors should have caught this)
  cell_iterator cell = begin(),
		endc = end();
  while (cell != endc)
    Assert (!(cell++)->refine_flag_set(), ExcInternalError ());  
#endif
};

#endif



#if deal_II_dimension == 2

template <>
void Triangulation<2>::execute_refinement () {
  const unsigned int dim = 2;
  
				   // check whether a new level is needed
				   // we have to check for this on the
				   // highest level only (on this, all
				   // used cells are also active, so we
				   // only have to check for this)
  if (true)
    {
      raw_cell_iterator cell = begin_active (levels.size()-1),
			endc = end();
      for (; cell != endc; ++cell)
	if (cell->used())
	  if (cell->refine_flag_set()) 
	    {
	      levels.push_back (new TriangulationLevel<dim>);
	      break;
	    };
    };



				   // check how much space is needed
				   // on every level
				   // we need not check the highest
				   // level since either
				   // - on the highest level no cells
				   //   are flagged for refinement
				   // - there are, but prepare_refinement
				   //   added another empty level
  unsigned int needed_vertices = 0;
  for (int level=levels.size()-2; level>=0; --level)
    {
      				       // count number of flagged cells on
				       // this level and compute how many
				       // new vertices and new lines will
				       // be needed
      unsigned int flagged_cells = 0;
      unsigned int needed_lines  = 0;
      active_cell_iterator acell = begin_active(level),
			   aendc = begin_active(level+1);
      for (; acell!=aendc; ++acell) 
	if (acell->refine_flag_set()) 
	  {
	    ++flagged_cells;

					     // new vertex at center of cell
					     // is needed in any case
	    ++needed_vertices;
					     //	also the four inner lines
	    needed_lines += 4;
	    
					     // for all neighbors of
					     // this cell
	    for (unsigned int nb=0; nb<GeometryInfo<dim>::faces_per_cell; ++nb) 
	      {
		const cell_iterator neighbor = acell->neighbor(nb);
						 // if cell is at boundary
		if (neighbor.state() != valid) 
		  {
						     // new midpoint vertex
						     // necessary
		    ++needed_vertices;
						     // also two new lines
		    needed_lines += 2;
		    
		    continue;
		  };
						 // there is a neighbor. There
						 // are three cases:
						 // 1 nb is on same level and
						 //   not refined (subcases:
						 //   flagged for refinement
						 //   or not)
						 // 2 nb is on same level and
						 //   refined (->no additional
						 //   vertices and lines needed)
						 // 3 nb is one level down
						 //   (but will be refined)
		if ((neighbor->level() == acell->level()) &&
		    (neighbor->active() == true))
		  {
						     // case 1
		    if (((neighbor->refine_flag_set() == true) &&
			 (acell->index() < neighbor->index()))
							  // case 1a
							  // we need one more vertex
							  // and two more lines, but
							  // we must only count them
							  // once. Convention: count
							  // them for the cell with
							  // the lower index
			||
			(neighbor->refine_flag_set() == false))
						       // case 1b
		      {
			++needed_vertices;
			needed_lines += 2;
		      };

		    continue;
		  };

		if ((neighbor->level() == acell->level()) &&
		    (neighbor->active() == false))
						   // case 2
		  continue;
		
		if (neighbor->level() == acell->level()-1)
						   // case 3
		  {
		    ++needed_vertices;
		    needed_lines += 2;
		    
		    continue;
		  };

		Assert (false, ExcUncaughtState());
	      };
	  };
      
      				       // count number of used cells on
				       // the next higher level
      const unsigned int used_cells
	= count_if (levels[level+1]->quads.used.begin(),
		levels[level+1]->quads.used.end(),
		bind2nd (equal_to<bool>(), true));
      
      
				       // reserve space for the used_cells
				       // cells already existing on the next
				       // higher level as well as for the
				       // 4*flagged_cells that will be created
				       // on that level
      levels[level+1]->
	TriangulationLevel<0>::reserve_space (used_cells+4*flagged_cells, 2);
				       // reserve space for needed_lines
				       // new lines
      levels[level+1]->
	TriangulationLevel<1>::reserve_space (needed_lines);
      				       // reserve space for 4*flagged_cells
				       // new quads on the next higher
				       // level
      levels[level+1]->
	TriangulationLevel<2>::reserve_space (4*flagged_cells);
    };

				   // add to needed vertices how
				   // many vertices are already in use
  needed_vertices += count_if (vertices_used.begin(), vertices_used.end(),
			       bind2nd (equal_to<bool>(), true));
				   // if we need more vertices: create them,
				   // if not: leave the array as is, since
				   // shrinking is not really possible because
				   // some of the vertices at the end may be
				   // in use
  if (needed_vertices > vertices.size())
    {
      vertices.resize (needed_vertices, Point<dim>());
      vertices_used.resize (needed_vertices, false);
    };


				   // Do REFINEMENT
				   // on every level; exclude highest level as
				   // above

				   //  index of next unused vertex
  unsigned int next_unused_vertex = 0;
  
  for (int level=0; level<static_cast<int>(levels.size())-1; ++level) 
    {
      
      active_cell_iterator cell = begin_active(level),
			   endc = begin_active(level+1);
      
      raw_line_iterator next_unused_line = begin_raw_line (level+1);
      raw_cell_iterator next_unused_cell = begin_raw (level+1);

      for (; (cell!=endc) && (cell->level()==level); ++cell) 
	if (cell->refine_flag_set()) 
	  {
					     // clear refinement flag
	    cell->clear_refine_flag ();

					     // do some additional checks.
#ifdef DEBUG
	    for (unsigned int neighbor=0;
		 neighbor<GeometryInfo<dim>::faces_per_cell; ++neighbor)
	      if (cell->neighbor(neighbor).state() == valid)
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
   (the array consists of indices into the #vertices# vector
   
     6--5--4
     |  |  |
     7--8--3
     |  |  |
     0--1--2
	
   Second:  
   Set up an array of the new lines (the array consists of iterator pointers
   into the lines arrays)
   
     .-5-.-4-.         The directions are:  .->-.->-.
     6   9   3                              ^   ^   ^
     .-10.11- .                             .->-.->-.
     7   8   2                              ^   ^   ^
     .-0-.-1-.                              .->-.->-.

   Please note that since the children of line are created in the direction of
   that line, the lines 4,5 and 6,7 are created in the wrong time order. This
   has the consequence that if n be the next free line number before the
   refinement process, the line numbered with 4 above will get index n+5,
   while the line number 5 above will get the index n+4. The same applies
   to the lines 6 and 7.
     
   Third:
   Set up an array of neighbors:
   
       5   4
      .--.--.
     6|  |  |3
      .--.--.
     7|  |  |2
      .--.--.
       0   1

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
   #neighbors_neighbor#.

   It would be sufficient to use the children of the common line to the
   neighbor, if we only wanted to get the new sublines and the new vertex,
   but because we need to update the neighborship information of the
   two refined subcells of the neighbor, we need to search these anyway.

   Convention:
   The created children are numbered like this:

     .--.--.
     |3 . 2|
     .--.--.
     |0 | 1|
     .--.--.
   */

	    int               new_vertices[9] = {cell->vertex_index(0), -1,
						 cell->vertex_index(1), -1,
						 cell->vertex_index(2), -1,
						 cell->vertex_index(3), -1,
						 -1};
	    raw_line_iterator new_lines[12];
	    cell_iterator     neighbors[8] = {cell->neighbor(0),
					      cell->neighbor(0),
					      cell->neighbor(1),
					      cell->neighbor(1),
					      cell->neighbor(2),
					      cell->neighbor(2),
					      cell->neighbor(3),
					      cell->neighbor(3)};
	    int               neighbors_neighbor[8] = {-1,-1,-1,-1,-1,-1,-1,-1};

					     // remember: the #i#th line
					     // is the common line to the
					     // #i#th neighbor
	    for (unsigned int nb=0; nb<4; ++nb)
	      {
		bool neighbor_refined=false;
		if (cell->neighbor(nb).state() == valid)
		  if (cell->neighbor(nb)->active() == false)
						   // (ask in two if-statements,
						   // since otherwise both
						   // conditions would be executed,
						   // but the second will throw an
						   // error if the first fails!)
		    neighbor_refined=true;
		
		if (neighbor_refined)
		  {
						     // neighbor exists and is
						     // refined
						     // ->the common line has
						     // two children which
						     // we can use.
		    cell_iterator neighbor = cell->neighbor(nb);
		    for (unsigned int nb_nb=0; nb_nb<4; ++nb_nb)
		      if (neighbor->neighbor(nb_nb)==cell)
							 // this cell is the nb_nb-th
							 // neighbor or neighbor(nb)
			{
			  neighbors_neighbor[2*nb] = neighbors_neighbor[2*nb+1] = nb_nb;
							   // vertex 1 of child 0
							   // is always the interior
							   // one
			  new_vertices[2*nb+1] = neighbor->line(nb_nb)
						 ->child(0)->vertex_index(1);

			  if (nb < 2) 
			    {
			      new_lines[2*nb]  = neighbor->line(nb_nb)->child(0);
			      new_lines[2*nb+1]= neighbor->line(nb_nb)->child(1);
			    } else {
							       // lines 2 and 3 have
							       // opposite sense
			      new_lines[2*nb]  = neighbor->line(nb_nb)->child(1);
			      new_lines[2*nb+1]= neighbor->line(nb_nb)->child(0);
			    };

							   // finally find out which
							   // are the two neighbor
							   // subcells, adjacent to
							   // the two sublines
			  const int child_mapping[4][2] = {{0,1},{1,2},{3,2},{0,3}};     
			  if (nb < 2) 
			    {
			      neighbors[2*nb]  = neighbor->child(child_mapping[nb_nb][0]);
			      neighbors[2*nb+1]= neighbor->child(child_mapping[nb_nb][1]);
			    } else {
			      neighbors[2*nb]  = neighbor->child(child_mapping[nb_nb][1]);
			      neighbors[2*nb+1]= neighbor->child(child_mapping[nb_nb][0]);
			    };
			};
		  }
		else 
	    
						   // neighboring cell either
						   // does not exist or is
						   // not refined -> we need a
						   // new vertex and two new lines
		  {
						     // search for next unused vertex
		    while (vertices_used[next_unused_vertex] == true)
		      ++next_unused_vertex;
		    Assert (next_unused_vertex < vertices.size(),
			    ExcTooFewVerticesAllocated());

						     // where shall we put the new
						     // vertex?
		    Point<2> new_point;
		    
		    face_iterator face=cell->line(nb);
		    
		    if ( face->boundary_indicator() != 255 )
		      {
						       // boundary vertex
		      new_point = boundary[face->boundary_indicator()]->
				  get_new_point_on_line (face);
		      } else {
							 // vertex between two
							 // normal cells
			new_point = vertices[new_vertices[2*nb]];
			new_point += vertices[new_vertices[(2*nb+2)%8]];
			new_point /= 2.0;
		      };
		    
		    new_vertices[nb*2+1] = next_unused_vertex;
		    vertices[new_vertices[nb*2+1]] = new_point;
		    vertices_used[new_vertices[nb*2+1]] = true;

						     // search for next unused line
						     // (++ takes care of the end of
						     // the vector)
		    while (next_unused_line->used() == true)
		      ++next_unused_line;

		    cell->line(nb)->set_children (next_unused_line->index());
		    
		    if (nb<2) 
		      {
			new_lines[nb*2] = next_unused_line;
			++next_unused_line;
			Assert (next_unused_line->used() == false,
				ExcCellShouldBeUnused());
			new_lines[nb*2+1] = next_unused_line;

			new_lines[nb*2]->set(Line(new_vertices[2*nb],
						  new_vertices[2*nb+1]));
			new_lines[nb*2]->set_used_flag ();
			new_lines[nb*2]->clear_children ();
			new_lines[nb*2]->clear_user_pointer ();
			
			new_lines[nb*2+1]->set(Line(new_vertices[2*nb+1],
						    new_vertices[(2*nb+2)%8]));
			new_lines[nb*2+1]->set_used_flag ();
			new_lines[nb*2+1]->clear_children ();
			new_lines[nb*2+1]->clear_user_pointer ();
		      } else {
			new_lines[nb*2+1] = next_unused_line;
			++next_unused_line;
			Assert (next_unused_line->used() == false,
				ExcCellShouldBeUnused());
			new_lines[nb*2] = next_unused_line;

			new_lines[nb*2]->set(Line(new_vertices[2*nb+1],
						    new_vertices[2*nb]));
			new_lines[nb*2]->set_used_flag ();
			new_lines[nb*2]->clear_children ();
			new_lines[nb*2]->clear_user_pointer ();

			new_lines[nb*2+1]->set(Line(new_vertices[(2*nb+2)%8],
						    new_vertices[2*nb+1]));
			new_lines[nb*2+1]->set_used_flag ();
			new_lines[nb*2+1]->clear_children ();
			new_lines[nb*2+1]->clear_user_pointer ();
		      };
		  };
	      };
	    
					     // add new vertex in the middle
					     // search for next unused
					     // vertex
	    while (vertices_used[next_unused_vertex] == true)
	      ++next_unused_vertex;
	    Assert (next_unused_vertex < vertices.size(),
		    ExcTooFewVerticesAllocated());

					     // new vertex is placed at the
					     // arithmetic mean of all 8
					     // neighboring points.
	    Point<2> new_point(0,0);
	    for (unsigned int i=0; i<8; ++i)
	      new_point +=  vertices[new_vertices[i]];
	    new_point /= 8.0;
	    
	    new_vertices[8] = next_unused_vertex;
	    vertices[new_vertices[8]] = new_point;
	    vertices_used[new_vertices[8]] = true;
	    
					     // add the 4 inner lines
	    
					     // search for next unused line
	    while (next_unused_line->used() == true)
	      ++next_unused_line;
	    new_lines[8] = next_unused_line;
	    new_lines[8]->set(Line(new_vertices[1],
				   new_vertices[8]));
	    new_lines[8]->set_used_flag ();
	    new_lines[8]->clear_children ();
	    new_lines[8]->clear_user_pointer ();
	    
	    while (next_unused_line->used() == true)
	      ++next_unused_line;
	    new_lines[9] = next_unused_line;
	    new_lines[9]->set(Line(new_vertices[8],
				   new_vertices[5]));
	    new_lines[9]->set_used_flag ();
	    new_lines[9]->clear_children ();
	    new_lines[9]->clear_user_pointer ();

	    while (next_unused_line->used() == true)
	      ++next_unused_line;
	    new_lines[10] = next_unused_line;
	    new_lines[10]->set(Line(new_vertices[7],
				    new_vertices[8]));
	    new_lines[10]->set_used_flag ();
	    new_lines[10]->clear_children ();
	    new_lines[10]->clear_user_pointer ();
	    
	    while (next_unused_line->used() == true)
	      ++next_unused_line;
	    new_lines[11] = next_unused_line;
	    new_lines[11]->set(Line(new_vertices[8],
				    new_vertices[3]));
	    new_lines[11]->set_used_flag ();
	    new_lines[11]->clear_children ();
	    new_lines[11]->clear_user_pointer ();
					     // set the boundary indicators of
					     // the outer cells.
	    new_lines[0]->set_boundary_indicator (cell->line(0)->boundary_indicator());
	    new_lines[1]->set_boundary_indicator (cell->line(0)->boundary_indicator());
	    new_lines[2]->set_boundary_indicator (cell->line(1)->boundary_indicator());
	    new_lines[3]->set_boundary_indicator (cell->line(1)->boundary_indicator());
	    new_lines[4]->set_boundary_indicator (cell->line(2)->boundary_indicator());
	    new_lines[5]->set_boundary_indicator (cell->line(2)->boundary_indicator());
	    new_lines[6]->set_boundary_indicator (cell->line(3)->boundary_indicator());
	    new_lines[7]->set_boundary_indicator (cell->line(3)->boundary_indicator());
					     // inner cells have boundary
					     // indicator 255
	    new_lines[8]->set_boundary_indicator (255);
	    new_lines[9]->set_boundary_indicator (255);
	    new_lines[10]->set_boundary_indicator (255);
	    new_lines[11]->set_boundary_indicator (255);

	    
					     // finally add the four new cells!
	    
					     // search for next unused cell
					     // the four children have to be put
					     // into the array consecutively
	    while (next_unused_cell->used() == true)
	      ++next_unused_cell;

	    raw_cell_iterator subcells[4];
	    for (unsigned int i=0; i<4; ++i) 
	      {
		Assert (next_unused_cell->used() == false,
			ExcCellShouldBeUnused());
		subcells[i] = next_unused_cell;
		++next_unused_cell;
	      };
	    
	    
	    
	    cell->set_children (subcells[0]->index());
	    
	    subcells[0]->set (Quad(new_lines[0]->index(),  new_lines[8]->index(),
				   new_lines[10]->index(), new_lines[7]->index()));
	    subcells[0]->set_used_flag();
	    subcells[0]->clear_children();
	    subcells[0]->clear_user_pointer ();

	    subcells[1]->set (Quad(new_lines[1]->index(),  new_lines[2]->index(),
				   new_lines[11]->index(), new_lines[8]->index()));
	    subcells[1]->set_used_flag();
	    subcells[1]->clear_children();
	    subcells[1]->clear_user_pointer ();


	    subcells[2]->set (Quad(new_lines[11]->index(),  new_lines[3]->index(),
				   new_lines[4]->index(), new_lines[9]->index()));
	    subcells[2]->set_used_flag();
	    subcells[2]->clear_children();
	    subcells[2]->clear_user_pointer ();
	    

	    subcells[3]->set (Quad(new_lines[10]->index(),  new_lines[9]->index(),
				   new_lines[5]->index(), new_lines[6]->index()));
	    subcells[3]->set_used_flag();
	    subcells[3]->clear_children();
	    subcells[3]->clear_user_pointer ();
	    
					     // finally set neighborship info of
					     // external cells
					     // (neighbor_mapping is the mapping
					     // between the 8 neighbors and the
					     // adjacent new cells in the interior)
	    const int neighbor_mapping[8] = {0,1, 1,2, 2,3, 3,0};
	    
	    for (unsigned int nb=0; nb<8; ++nb)
	      if (neighbors[nb].state() == valid)
		if (neighbors[nb]->level() == level+1)
						   // neighbor is refined cell
		  neighbors[nb]->set_neighbor(neighbors_neighbor[nb],
					      subcells[neighbor_mapping[nb]]);

					     // and neighbarship of new cells
	    subcells[0]->set_neighbor (0, neighbors[0]);
	    subcells[0]->set_neighbor (1, subcells[1]);
	    subcells[0]->set_neighbor (2, subcells[3]);
	    subcells[0]->set_neighbor (3, neighbors[7]);

	    subcells[1]->set_neighbor (0, neighbors[1]);
	    subcells[1]->set_neighbor (1, neighbors[2]);
	    subcells[1]->set_neighbor (2, subcells[2]);
	    subcells[1]->set_neighbor (3, subcells[0]);

	    subcells[2]->set_neighbor (0, subcells[1]);
	    subcells[2]->set_neighbor (1, neighbors[3]);
	    subcells[2]->set_neighbor (2, neighbors[4]);
	    subcells[2]->set_neighbor (3, subcells[3]);

	    subcells[3]->set_neighbor (0, subcells[0]);
	    subcells[3]->set_neighbor (1, subcells[2]);
	    subcells[3]->set_neighbor (2, neighbors[5]);
	    subcells[3]->set_neighbor (3, neighbors[6]);

	    subcells[0]->set_material_id (cell->material_id());
	    subcells[1]->set_material_id (cell->material_id());
	    subcells[2]->set_material_id (cell->material_id());
	    subcells[3]->set_material_id (cell->material_id());
	  };
    };

				   // re-compute number of lines and quads
  update_number_cache ();


#ifdef DEBUG
  for (unsigned int level=0; level<levels.size(); ++level) 
    levels[level]->monitor_memory (2);

				   // check whether really all refinement flags
				   // are reset (also of previously non-active
				   // cells which we may not have touched. If
				   // the refinement flag of a non-active cell
				   // is set, something went wrong since the
				   // cell-accessors should have caught this)
  cell_iterator cell = begin(),
		endc = end();
  while (cell != endc)
    Assert (!(cell++)->refine_flag_set(), ExcInternalError ());
#endif
};

#endif



#if deal_II_dimension == 3

template <>
void Triangulation<3>::execute_refinement () {
  const unsigned int dim = 3;

				   // check whether a new level is needed
				   // we have to check for this on the
				   // highest level only (on this, all
				   // used cells are also active, so we
				   // only have to check for this)
  if (true)
    {
      raw_cell_iterator cell = begin_active (levels.size()-1),
			endc = end();
      for (; cell != endc; ++cell)
	if (cell->used())
	  if (cell->refine_flag_set()) 
	    {
	      levels.push_back (new TriangulationLevel<dim>);
	      break;
	    };
    };


				   // first clear user flags for
				   // quads and lines; we're gonna
				   // use them to flag which lines
				   // and quads need refinement
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
  for (int level=levels.size()-2; level>=0; --level)
    {
      				       // count number of flagged cells on
				       // this level and compute how many
				       // new vertices and new lines will
				       // be needed
      unsigned int flagged_cells = 0;
      unsigned int needed_lines  = 0;
      unsigned int needed_quads  = 0;
      
      active_cell_iterator acell = begin_active(level),
			   aendc = begin_active(level+1);
      for (; acell!=aendc; ++acell) 
	if (acell->refine_flag_set()) 
	  {
	    ++flagged_cells;

					     // new vertex at center of cell
					     // is needed in any case
	    ++needed_vertices;
					     //	also the six inner lines
	    needed_lines += 6;
					     // and the 12 inner quads
	    needed_quads += 12;
	    
					     // mark all faces and lines for
					     // refinement; checking locally
					     // whether the neighbor would
					     // also like to refine them is
					     // rather difficult for lines
					     // so we only flag them and after
					     // visiting all cells, we decide
					     // which lines need refinement;
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
		  };
	      };
	  };

				       // now count the quads and lines which
				       // were flagged for refinement
      for (quad_iterator quad=begin_quad(); quad!=end_quad(); ++quad)
	if (quad->user_flag_set())
	  {
	    Assert (quad->has_children() == false, ExcInternalError());
	    needed_quads    += 4;
	    needed_lines    += 4;
	    needed_vertices += 1;
	  };

      for (line_iterator line=begin_line(); line!=end_line(); ++line)
	if (line->user_flag_set())
	  {
	    Assert (line->has_children() == false, ExcInternalError());
	    needed_lines += 2;
	    needed_vertices += 1;
	  };
      
      
      				       // count number of used cells on
				       // the next higher level
      const unsigned int used_cells
	= count_if (levels[level+1]->quads.used.begin(),
		levels[level+1]->quads.used.end(),
		bind2nd (equal_to<bool>(), true));
      
      
				       // reserve space for the used_cells
				       // cells already existing on the next
				       // higher level as well as for the
				       // 8*flagged_cells that will be created
				       // on that level
      levels[level+1]->
	TriangulationLevel<0>::reserve_space (used_cells+8*flagged_cells, 3);
				       // reserve space for needed_lines
				       // new lines
      levels[level+1]->
	TriangulationLevel<1>::reserve_space (needed_lines);
				       // reserve space for needed_quads
				       // new quads
      levels[level+1]->
	TriangulationLevel<2>::reserve_space (needed_quads);
      				       // reserve space for 8*flagged_cells
				       // new hexes on the next higher
				       // level
      levels[level+1]->
	TriangulationLevel<3>::reserve_space (8*flagged_cells);
    };

				   // add to needed vertices how
				   // many vertices are already in use
  needed_vertices += count_if (vertices_used.begin(), vertices_used.end(),
			       bind2nd (equal_to<bool>(), true));
				   // if we need more vertices: create them,
				   // if not: leave the array as is, since
				   // shrinking is not really possible because
				   // some of the vertices at the end may be
				   // in use
  if (needed_vertices > vertices.size())
    {
      vertices.resize (needed_vertices, Point<dim>());
      vertices_used.resize (needed_vertices, false);
    };


				   ///////////////////////////////////////////
				   // Before we start with the actual
				   // refinement, we do some sanity checks if
				   // in debug mode. especially, we try to
				   // catch the notorious problem with
				   // lines being twice refined, i.e. there
				   // are cells adjacent at one line ("around
				   // the edge", but not at a face), with
				   // two cells differing by more than one
				   // refinement level
				   //
				   // this check is very simple to implement
				   // here, since we have all lines flagged if
				   // they shall be refined
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
				   // To make life a bit easier, we first
				   // refine those lines and quads that were
				   // flagged for refinement and then compose
				   // the newly to be created cells.
				   //
				   //  index of next unused vertex
  unsigned int next_unused_vertex = 0;

				   // first for lines
  for (unsigned int level=0; level!=levels.size()-1; ++level)
    {
				       // only active objects can be refined
				       // further; remember that we won't
				       // operate on the finest level, so
				       // begin_*(level+1) is allowed
      active_line_iterator line = begin_active_line(level),
			   endl = begin_active_line(level+1);
      raw_line_iterator next_unused_line = begin_raw_line (level+1);

      for (; line!=endl; ++line)
	if (line->user_flag_set())
	  {
					     // this line needs to be refined

					     // find the next unused vertex
					     // and set it appropriately
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
	  
					     // now that we created the right
					     // point, make up the two child
					     // lines
					     // (++ takes care of the end of
					     // the vector)
	    while (next_unused_line->used() == true)
	      ++next_unused_line;
					     // there should always be two
					     // consecutive unused lines, such
					     // that the children of a line
					     // will be consecutive.
					     // then set the child pointer of
					     // the present line
	    line->set_children (next_unused_line->index());
	  
					     // set the two new lines
	    raw_line_iterator children[2] = { next_unused_line,
					      ++next_unused_line };
					     // some tests; if any of the
					     // iterators should be invalid,
					     // then already dereferencing
					     // will fail
	    Assert (children[0]->used() == false, ExcCellShouldBeUnused());
	    Assert (children[1]->used() == false, ExcCellShouldBeUnused());
	  
	    children[0]->set (Line(line->vertex_index(0),
				   next_unused_vertex));
	    children[1]->set (Line(next_unused_vertex,
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
	    
					     // finally clear flag indicating
					     // the need for refinement
	    line->clear_user_flag ();
	  };
    };


				   ///////////////////////////////////////
				   // now refine marked quads
				   ///////////////////////////////////////
  for (unsigned int level=0; level!=levels.size()-1; ++level)
    {
				       // only active objects can be refined
				       // further; remember that we won't
				       // operate on the finest level, so
				       // begin_*(level+1) is allowed
      active_quad_iterator quad = begin_active_quad(level),
			   endq = begin_active_quad(level+1);
      raw_line_iterator next_unused_line = begin_raw_line (level+1);
      raw_quad_iterator next_unused_quad = begin_raw_quad (level+1);

      for (; quad!=endq; ++quad)
	if (quad->user_flag_set())
	  {
					     // this quad needs to be refined

					     // find the next unused vertex
					     // and set it appropriately
	    while (vertices_used[next_unused_vertex] == true)
	      ++next_unused_vertex;
	    Assert (next_unused_vertex < vertices.size(),
		    ExcTooFewVerticesAllocated());
	    vertices_used[next_unused_vertex] = true;
	    
	    if (quad->at_boundary()) 
	      vertices[next_unused_vertex]
		= boundary[quad->boundary_indicator()]->get_new_point_on_quad (quad);
	    else
	      vertices[next_unused_vertex]
		= (quad->vertex(0) + quad->vertex(1) +
		   quad->vertex(2) + quad->vertex(3)) / 4;
	  
					     // now that we created the right
					     // point, make up the four
					     // lines interior to the quad
					     // (++ takes care of the end of
					     // the vector)
	    raw_line_iterator new_lines[4];

	    for (unsigned int i=0; i<4; ++i)
	      {
		while (next_unused_line->used() == true)
		  ++next_unused_line;
		new_lines[i] = next_unused_line;
		++next_unused_line;

		Assert (new_lines[i]->used() == false, ExcCellShouldBeUnused());
	      };

					     // set the data of the four lines.
					     // first collect the indices of
					     // the five vertices:
					     // *--2--*
					     // |  |  |
					     // 3--4--1
					     // |  |  |
					     // *--0--*
					     // the lines are numbered as follows:
					     // *--*--*
					     // |  2  |
					     // *3-*-1*
					     // |  0  |
					     // *--*--*
	    const unsigned int vertex_indices[5]
	      = { quad->line(0)->child(0)->vertex_index(1),
		  quad->line(1)->child(0)->vertex_index(1),
		  quad->line(2)->child(0)->vertex_index(1),
		  quad->line(3)->child(0)->vertex_index(1),
		  next_unused_vertex 
	      };
	    
	    new_lines[0]->set (Line(vertex_indices[0], vertex_indices[4]));
	    new_lines[1]->set (Line(vertex_indices[4], vertex_indices[1]));
	    new_lines[2]->set (Line(vertex_indices[4], vertex_indices[2]));
	    new_lines[3]->set (Line(vertex_indices[3], vertex_indices[4]));

	    for (unsigned int i=0; i<4; ++i)
	      {
		new_lines[i]->set_used_flag();
		new_lines[i]->clear_user_flag();
		new_lines[i]->clear_user_pointer();
		new_lines[i]->clear_children();
		new_lines[i]->set_boundary_indicator(quad->boundary_indicator());
	      };



					     // now for the quads. again, first
					     // collect some data about the
					     // indices of the lines, with
					     // the following numbering:
					     // *5-*-4*
					     // 6  10 3
					     // *11*-9*
					     // 7  8  2
					     // *0-*-1*
	    const unsigned int line_indices[12]
	      = { quad->line(0)->child(0)->index(),
		  quad->line(0)->child(1)->index(),
		  quad->line(1)->child(0)->index(),
		  quad->line(1)->child(1)->index(),
		  quad->line(2)->child(1)->index(),
		  quad->line(2)->child(0)->index(),
		  quad->line(3)->child(1)->index(),
		  quad->line(3)->child(0)->index(),
		  new_lines[0]->index(),
		  new_lines[1]->index(),
		  new_lines[2]->index(),
		  new_lines[3]->index() 
	      };
	    
					     // find some space for the four
					     // newly to be created quads.
					     // note that there should
					     // always be four consecutive
					     // free slots for them
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

					     // note these quads as children
					     // to the present one
	    quad->set_children (new_quads[0]->index());

	    new_quads[0]->set (Quad(line_indices[0],
				    line_indices[8],
				    line_indices[11],
				    line_indices[7]));
	    new_quads[1]->set (Quad(line_indices[1],
				    line_indices[2],
				    line_indices[9],
				    line_indices[8]));
	    new_quads[2]->set (Quad(line_indices[9],
				    line_indices[3],
				    line_indices[4],
				    line_indices[10]));
	    new_quads[3]->set (Quad(line_indices[11],
				    line_indices[10],
				    line_indices[5],
				    line_indices[6]));
	    for (unsigned int i=0; i<4; ++i)
	      {
		new_quads[i]->set_used_flag();
		new_quads[i]->clear_user_flag();
		new_quads[i]->clear_user_pointer();
		new_quads[i]->clear_children();
		new_quads[i]->set_boundary_indicator (quad->boundary_indicator());
	      };  
	    
					     // finally clear flag indicating
					     // the need for refinement
	    quad->clear_user_flag ();
	  };
    };

				   ///////////////////////////////////
				   // Now, finally, set up the new
				   // cells
				   ///////////////////////////////////
  for (unsigned int level=0; level!=levels.size()-1; ++level)
    {
				       // only active objects can be refined
				       // further; remember that we won't
				       // operate on the finest level, so
				       // begin_*(level+1) is allowed
      active_hex_iterator hex  = begin_active_hex(level),
			  endh = begin_active_hex(level+1);
      raw_line_iterator next_unused_line = begin_raw_line (level+1);
      raw_quad_iterator next_unused_quad = begin_raw_quad (level+1);
      raw_hex_iterator  next_unused_hex  = begin_raw_hex (level+1);

      for (; hex!=endh; ++hex)
	if (hex->refine_flag_set())
	  {
					     // do some additional checks.
#ifdef DEBUG
	    for (unsigned int neighbor=0;
		 neighbor<GeometryInfo<dim>::faces_per_cell; ++neighbor)
	      if (hex->neighbor(neighbor).state() == valid)
		Assert (((hex->neighbor(neighbor)->level() == hex->level()) &&
			 (hex->neighbor(neighbor)->coarsen_flag_set() == false))  ||
			((hex->neighbor(neighbor)->level() == hex->level()-1) &&
			 (hex->neighbor(neighbor)->refine_flag_set() == true)),
			ExcInternalError());
#endif
					     // this hex needs to be refined
	    
					     // clear flag indicating
					     // the need for refinement. do it
					     // here already, since we can't
					     // do it anymore once the cell
					     // has children
	    hex->clear_refine_flag ();

					     // find the next unused vertex
					     // and set it appropriately
	    while (vertices_used[next_unused_vertex] == true)
	      ++next_unused_vertex;
	    Assert (next_unused_vertex < vertices.size(),
		    ExcTooFewVerticesAllocated());
	    vertices_used[next_unused_vertex] = true;
	    
					     // the new vertex is definitely in
					     // the interior, so we need not
					     // worry about the boundary.
					     // let it be the average of
					     // the 26 vertices surrounding
					     // it
	    vertices[next_unused_vertex] = Point<dim>();
					     // first add corners of hex
	    for (unsigned int vertex=0;
		 vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
	      vertices[next_unused_vertex] += hex->vertex(vertex);
					     // now add center of lines
	    for (unsigned int line=0; line<12; ++line)
	      vertices[next_unused_vertex] += hex->line(line)->child(0)->vertex(1);
					     // finally add centers of faces
	    for (unsigned int face=0;
		 face<GeometryInfo<dim>::faces_per_cell; ++face)
	      vertices[next_unused_vertex] += hex->face(face)->child(0)->vertex(2);

					     // compute average
	    vertices[next_unused_vertex] /= 26;
	    
	  
					     // now that we created the right
					     // point, make up the six
					     // lines interior to the quad
					     // (++ takes care of the end of
					     // the vector)
	    raw_line_iterator new_lines[6];

	    for (unsigned int i=0; i<6; ++i)
	      {
		while (next_unused_line->used() == true)
		  ++next_unused_line;
		new_lines[i] = next_unused_line;
		++next_unused_line;

		Assert (new_lines[i]->used() == false, ExcCellShouldBeUnused());
	      };

					     // set the data of the six lines.
					     // first collect the indices of
					     // the seven vertices (consider
					     // the two planes to be crossed
					     // to form the planes cutting
					     // the hex in two vertically
					     // and horizontally)
					     //     *--2--*   *--5--*
					     //    /  /  /    |  |  |
					     //   3--6--1     3--6--1
					     //  /  /  /      |  |  |
					     // *--0--*       *--4--*
					     // the lines are numbered as follows:
					     //     *--*--*   *--*--*
					     //    /  2  /    |  5  |
					     //   *3-*-1*     *3-*-1*
					     //  /  0  /      |  4  |
					     // *--*--*       *--*--*
	    const unsigned int vertex_indices[7]
	      = { hex->face(0)->child(0)->vertex_index(2),
		  hex->face(3)->child(0)->vertex_index(2),
		  hex->face(1)->child(0)->vertex_index(2),
		  hex->face(5)->child(0)->vertex_index(2),
		  hex->face(2)->child(0)->vertex_index(2),
		  hex->face(4)->child(0)->vertex_index(2),
		  next_unused_vertex 
	      };
	    
	    new_lines[0]->set (Line(vertex_indices[0], vertex_indices[6]));
	    new_lines[1]->set (Line(vertex_indices[6], vertex_indices[1]));
	    new_lines[2]->set (Line(vertex_indices[6], vertex_indices[2]));
	    new_lines[3]->set (Line(vertex_indices[3], vertex_indices[6]));
	    new_lines[4]->set (Line(vertex_indices[4], vertex_indices[6]));
	    new_lines[5]->set (Line(vertex_indices[6], vertex_indices[5]));

	    for (unsigned int i=0; i<6; ++i)
	      {
		new_lines[i]->set_used_flag();
		new_lines[i]->clear_user_flag();
		new_lines[i]->clear_user_pointer();
		new_lines[i]->clear_children();
						 // interior line
		new_lines[i]->set_boundary_indicator(255);
	      };



					     // now for the quads. again, first
					     // collect some data about the
					     // indices of the lines, with
					     // the following numbering:
					     // front plane    *---*---*
					     //                |   2   |
					     //                *3--*--1*
					     //                |   0   |
					     //                *---*---*
					     //
					     // middle plane   *9--*--8*
					     //                10  14  7 
					     //                *15-*-13*
					     //                11  12  6 
					     //                *4--*--5*
					     //
					     // back plane     *---*---*
					     //                |   18  | 
					     //                *19-*-17*
					     //                |   16  | 
					     //                *---*---*
					     //
					     // left plane (the left-to-right
					     // planes are displayed twice each,
					     // for better readability of the
					     // indices; the left part is
					     // already determined by the
					     // pictures above)
					     //                  *            *
					     //                 /|           /| 
					     //                * |          * |
					     //               /| *         /| *
					     //              *10/|        * |21
					     //              | * |        | * |
					     //              |/| *        |20 *
					     //              *11/         * |/
					     //              | *          | *
					     //              |/           |/ 
					     //              *            *
					     //
					     // middle plane
					     //                  *            *
					     //                 /|          23| 
					     //                * 18         * |
					     //               /| *        22| *
					     //              *14/16       * |25
					     //              | * |        | * |
					     //              2/| *        |24 *
					     //              *12/         * |27
					     //              | *          | *
					     //              0/           |26
					     //              *            *
					     //
					     //
					     // right plane
					     //                  *            *
					     //                 /|           /| 
					     //                * |          * |
					     //               /| *         /| *
					     //              * 6/|        * |29
					     //              | * |        | * |
					     //              |/| *        |28 *
					     //              * 7/         * |/
					     //              | *          | *
					     //              |/           |/
					     //              *            *
	       
	    const unsigned int line_indices[30]
	      = {
		    hex->face(0)->child(0)->line_index(1),   //0
		    hex->face(0)->child(1)->line_index(2),   //1
		    hex->face(0)->child(2)->line_index(3),   //2
		    hex->face(0)->child(3)->line_index(0),   //3

		    hex->face(2)->child(0)->line_index(2),   //4
		    hex->face(2)->child(1)->line_index(2),   //5
		    hex->face(3)->child(0)->line_index(1),   //6
		    hex->face(3)->child(3)->line_index(1),   //7

		    hex->face(4)->child(1)->line_index(2),   //8
		    hex->face(4)->child(0)->line_index(2),   //9
		    hex->face(5)->child(3)->line_index(1),   //10
		    hex->face(5)->child(0)->line_index(1),   //11

		    new_lines[4]->index(),                   //12
		    new_lines[1]->index(),                   //13
		    new_lines[5]->index(),                   //14
		    new_lines[3]->index(),                   //15

		    hex->face(1)->child(0)->line_index(1),   //16
		    hex->face(1)->child(1)->line_index(2),   //17
		    hex->face(1)->child(2)->line_index(3),   //18
		    hex->face(1)->child(3)->line_index(0),   //19

		    hex->face(5)->child(0)->line_index(2),   //20
		    hex->face(5)->child(1)->line_index(2),   //21
		    hex->face(4)->child(0)->line_index(1),   //22
		    hex->face(4)->child(3)->line_index(1),   //23

		    new_lines[0]->index(),                   //24
		    new_lines[2]->index(),                   //25
		    hex->face(2)->child(0)->line_index(1),   //26
		    hex->face(2)->child(3)->line_index(1),   //27

		    hex->face(3)->child(0)->line_index(2),   //28
		    hex->face(3)->child(1)->line_index(2)    //29
	      };
	    
					     // find some space for the 12
					     // newly to be created quads.
	    raw_quad_iterator new_quads[12];

	    for (unsigned int i=0; i<12; ++i)
	      {
		while (next_unused_quad->used() == true)
		  ++next_unused_quad;
		new_quads[i] = next_unused_quad;
		++next_unused_quad;

		Assert (new_quads[i]->used() == false, ExcCellShouldBeUnused());
	      };

					     // set up the 12 quads, numbered
					     // as follows (shown are the three
					     // planes cutting the hex in two):
					     //
					     //  *-----*-----*
					     //  |  3  |  2  |
					     //  |     |     |
					     //  *-----*-----*
					     //  |     |     |
					     //  |  0  |  1  |
					     //  *-----*-----*
					     //
					     //       *----*----*
					     //      / 7  / 6  /
					     //     *----*----*
					     //    / 4  / 5  /
					     //   *----*----*
					     //
					     //
					     //             *
					     //            /|
					     //           / |
					     //          *10|
					     //         /|  *
					     //        / | /|
					     //       *  |/ |
					     //       |11*  |
					     //       | /| 9*
					     //       |/ | /
					     //       *  |/
					     //       |8 *
					     //       | /
					     //       |/
					     //       *
	    new_quads[0]->set (Quad(line_indices[4],
				    line_indices[12],
				    line_indices[15],
				    line_indices[11]));
	    new_quads[1]->set (Quad(line_indices[5],
				    line_indices[6],
				    line_indices[13],
				    line_indices[12]));
	    new_quads[2]->set (Quad(line_indices[13],
				    line_indices[7],
				    line_indices[8],
				    line_indices[14]));
	    new_quads[3]->set (Quad(line_indices[15],
				    line_indices[14],
				    line_indices[9],
				    line_indices[10]));
	    new_quads[4]->set (Quad(line_indices[3],
				    line_indices[24],
				    line_indices[15],
				    line_indices[20]));
	    new_quads[5]->set (Quad(line_indices[1],
				    line_indices[28],
				    line_indices[13],
				    line_indices[24]));
	    new_quads[6]->set (Quad(line_indices[13],
				    line_indices[29],
				    line_indices[17],
				    line_indices[25]));
	    new_quads[7]->set (Quad(line_indices[15],
				    line_indices[25],
				    line_indices[19],
				    line_indices[21]));
	    new_quads[8]->set (Quad(line_indices[26],
				    line_indices[12],
				    line_indices[24],
				    line_indices[0]));
	    new_quads[9]->set (Quad(line_indices[27],
				    line_indices[16],
				    line_indices[25],
				    line_indices[12]));
	    new_quads[10]->set (Quad(line_indices[25],
				    line_indices[18],
				    line_indices[23],
				    line_indices[14]));
	    new_quads[11]->set (Quad(line_indices[24],
				    line_indices[14],
				    line_indices[22],
				    line_indices[2]));
	    for (unsigned int i=0; i<12; ++i)
	      {
		new_quads[i]->set_used_flag();
		new_quads[i]->clear_user_flag();
		new_quads[i]->clear_user_pointer();
		new_quads[i]->clear_children();
						 // interior quad
		new_quads[i]->set_boundary_indicator (255);
	      };  


					     /////////////////////////////////
					     // create the eight new hexes
					     //
					     // again first collect some data.
					     // here, we need the indices of a
					     // whole lotta quads. they are
					     // numbered as follows:
					     //
					     // planes in the interior of
					     // the old hex:
					     //  *-----*-----*
					     //  |  3  |  2  |
					     //  |     |     |
					     //  *-----*-----*
					     //  |     |     |
					     //  |  0  |  1  |
					     //  *-----*-----*
					     //
					     //       *----*----*
					     //      / 7  / 6  /
					     //     *----*----*
					     //    / 4  / 5  /
					     //   *----*----*
					     //
					     //
					     //             *
					     //            /|
					     //           / |
					     //          *10|
					     //         /|  *
					     //        / | /|
					     //       *  |/ |
					     //       |11*  |
					     //       | /| 9*
					     //       |/ | /
					     //       *  |/
					     //       |8 *
					     //       | /
					     //       |/
					     //       *
					     //
					     // children of the faces of the
					     // old hex
					     //      *-------*        *-------*
					     //     /|19   18|       /31   30/|
					     //    34|       |      /       /26
					     //   /  |       |     /28   29/  |
					     //  *   |16   17|    *-------*27 |
					     //  3533*-------*    |15   14| 25*
					     //  |  /23   22/     |       |  /
					     //  32/       /      |       |24
					     //  |/20   21/       |12   13|/
					     //  *-------*        *-------*
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

		    hex->face(0)->child_index(0),  //12
		    hex->face(0)->child_index(1),
		    hex->face(0)->child_index(2),
		    hex->face(0)->child_index(3),

		    hex->face(1)->child_index(0),  //16
		    hex->face(1)->child_index(1),
		    hex->face(1)->child_index(2),
		    hex->face(1)->child_index(3),

		    hex->face(2)->child_index(0),  //20
		    hex->face(2)->child_index(1),
		    hex->face(2)->child_index(2),
		    hex->face(2)->child_index(3),

		    hex->face(3)->child_index(0),  //24
		    hex->face(3)->child_index(1),
		    hex->face(3)->child_index(2),
		    hex->face(3)->child_index(3),

		    hex->face(4)->child_index(0),  //28
		    hex->face(4)->child_index(1),
		    hex->face(4)->child_index(2),
		    hex->face(4)->child_index(3),

		    hex->face(5)->child_index(0),  //32
		    hex->face(5)->child_index(1),
		    hex->face(5)->child_index(2),
		    hex->face(5)->child_index(3)
	      };
	    
	    
					     // find some space for the eight
					     // newly to be created hexes.
					     // note that there should
					     // always be eight consecutive
					     // free slots for them
	    raw_hex_iterator new_hexes[8];

	    while (next_unused_hex->used() == true)
	      ++next_unused_hex;

	    for (unsigned int i=0; i<8; ++i)
	      {
		new_hexes[i] = next_unused_hex;
		Assert (new_hexes[i]->used() == false, ExcCellShouldBeUnused());
		++next_unused_hex;
	      };

					     // note these hexes as children
					     // to the present cell
	    hex->set_children (new_hexes[0]->index());

					     // front children
	    new_hexes[0]->set (Hexahedron(quad_indices[12],
					  quad_indices[0],
					  quad_indices[20],
					  quad_indices[8],
					  quad_indices[4],
					  quad_indices[32]));
	    new_hexes[1]->set (Hexahedron(quad_indices[13],
					  quad_indices[1],
					  quad_indices[21],
					  quad_indices[24],
					  quad_indices[5],
					  quad_indices[8]));
	    new_hexes[2]->set (Hexahedron(quad_indices[14],
					  quad_indices[2],
					  quad_indices[5],
					  quad_indices[27],
					  quad_indices[29],
					  quad_indices[11]));
	    new_hexes[3]->set (Hexahedron(quad_indices[15],
					  quad_indices[3],
					  quad_indices[4],
					  quad_indices[11],
					  quad_indices[28],
					  quad_indices[35]));
	    
					     // back children
	    new_hexes[4]->set (Hexahedron(quad_indices[0],
					  quad_indices[16],
					  quad_indices[23],
					  quad_indices[9],
					  quad_indices[7],
					  quad_indices[33]));
	    new_hexes[5]->set (Hexahedron(quad_indices[1],
					  quad_indices[17],
					  quad_indices[22],
					  quad_indices[25],
					  quad_indices[6],
					  quad_indices[9]));
	    new_hexes[6]->set (Hexahedron(quad_indices[2],
					  quad_indices[18],
					  quad_indices[6],
					  quad_indices[26],
					  quad_indices[30],
					  quad_indices[10]));
	    new_hexes[7]->set (Hexahedron(quad_indices[3],
					  quad_indices[19],
					  quad_indices[7],
					  quad_indices[10],
					  quad_indices[31],
					  quad_indices[34]));
	    
	    
	    for (unsigned int i=0; i<8; ++i)
	      {
		new_hexes[i]->set_used_flag();
		new_hexes[i]->clear_user_flag();
		new_hexes[i]->clear_user_pointer();
		new_hexes[i]->clear_children();
						 // inherit material properties
		new_hexes[i]->set_material_id (hex->material_id());
	      };


					     /////////////////////////////////
					     // now the only thing still to be
					     // done is setting neighborship
					     // information.
					     //
					     // to do so, first collect the
					     // iterators pointing to
					     // the 6x4 neighbors of this
					     // cell.
					     //
					     // note that in case the
					     // neighboring cell is not refined,
					     // the neighbor iterators point
					     // to the common moter cell. the
					     // same applies if there is no
					     // neighbor: the iterators are
					     // past the end
	    cell_iterator neighbor_cells[6][4];
	    for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell;
		 ++face)
	      {
		cell_iterator neighbor = hex->neighbor(face);

						 // if no neighbor
		if (neighbor.state() != valid)
		  for (unsigned int child_face=0;
		       child_face<GeometryInfo<dim>::subfaces_per_face;
		       ++child_face)
		    neighbor_cells[face][child_face] = neighbor;
		
		else
						   // neighbor exists
		  {
						     // neighbor's level must
						     // not be higher (else
						     // something went wrong
						     // when constructing
						     // either of the two
						     // cells) and not
						     // lower since then
						     // this cell should not
						     // have been refined.
		    Assert (neighbor->level() == hex->level(),
			    ExcInternalError());

						     // now there are two
						     // possibilities: either
						     // the nieghbor has no
						     // children or it has
						     // children. these
						     // must be terminal then.
		    if (!neighbor->has_children())
		      for (unsigned int child_face=0;
			   child_face<GeometryInfo<dim>::subfaces_per_face;
			   ++child_face)
			neighbor_cells[face][child_face] = neighbor;
		    else
						       // neighbor has
						       // children; now it's
						       // getting complicated
		      {
							 // first find the face
							 // of the neighbor
							 // adjacent to which
							 // the present cell is
			unsigned int nb_nb;
			for (nb_nb=0; nb_nb<GeometryInfo<dim>::faces_per_cell;
			     ++nb_nb)
			  if (neighbor->neighbor(nb_nb) == hex)
			    break;
			Assert (nb_nb<GeometryInfo<dim>::faces_per_cell,
				ExcInternalError());

							 // now the four child
							 // cells of neighbor
							 // adjacent to the
							 // present cell can
							 // be obtained by a
							 // function of
							 // GeometryInfo
			for (unsigned int c=0;
			     c<GeometryInfo<dim>::subfaces_per_face; ++c)
			  {
			    neighbor_cells[face][c]
			      = neighbor->child(GeometryInfo<dim>::child_cell_on_face(nb_nb, c));
			    
			    Assert (neighbor_cells[face][c].state() == valid,
				    ExcInternalError());
			    Assert (!neighbor_cells[face][c]->has_children(),
				    ExcInternalError());
			  };
		      };
		  };
	      };

					     // now we've got all neighbors, so set
					     // them in the new cells
	    new_hexes[0]->set_neighbor (0, neighbor_cells[0][0]);
	    new_hexes[0]->set_neighbor (1, new_hexes[4]);
	    new_hexes[0]->set_neighbor (2, neighbor_cells[2][0]);
	    new_hexes[0]->set_neighbor (3, new_hexes[1]);
	    new_hexes[0]->set_neighbor (4, new_hexes[3]);
	    new_hexes[0]->set_neighbor (5, neighbor_cells[5][0]);

	    new_hexes[1]->set_neighbor (0, neighbor_cells[0][1]);
	    new_hexes[1]->set_neighbor (1, new_hexes[5]);
	    new_hexes[1]->set_neighbor (2, neighbor_cells[2][1]);
	    new_hexes[1]->set_neighbor (3, neighbor_cells[3][0]);
	    new_hexes[1]->set_neighbor (4, new_hexes[2]);
	    new_hexes[1]->set_neighbor (5, new_hexes[0]);

	    new_hexes[2]->set_neighbor (0, neighbor_cells[0][2]);
	    new_hexes[2]->set_neighbor (1, new_hexes[6]);
	    new_hexes[2]->set_neighbor (2, new_hexes[1]);
	    new_hexes[2]->set_neighbor (3, neighbor_cells[3][3]);
	    new_hexes[2]->set_neighbor (4, neighbor_cells[4][1]);
	    new_hexes[2]->set_neighbor (5, new_hexes[3]);

	    new_hexes[3]->set_neighbor (0, neighbor_cells[0][3]);
	    new_hexes[3]->set_neighbor (1, new_hexes[7]);
	    new_hexes[3]->set_neighbor (2, new_hexes[0]);
	    new_hexes[3]->set_neighbor (3, new_hexes[2]);
	    new_hexes[3]->set_neighbor (4, neighbor_cells[4][0]);
	    new_hexes[3]->set_neighbor (5, neighbor_cells[5][3]);

	    new_hexes[4]->set_neighbor (0, new_hexes[0]);
	    new_hexes[4]->set_neighbor (1, neighbor_cells[1][0]);
	    new_hexes[4]->set_neighbor (2, neighbor_cells[2][3]);
	    new_hexes[4]->set_neighbor (3, new_hexes[5]);
	    new_hexes[4]->set_neighbor (4, new_hexes[7]);
	    new_hexes[4]->set_neighbor (5, neighbor_cells[5][1]);

	    new_hexes[5]->set_neighbor (0, new_hexes[1]);
	    new_hexes[5]->set_neighbor (1, neighbor_cells[1][1]);
	    new_hexes[5]->set_neighbor (2, neighbor_cells[2][2]);
	    new_hexes[5]->set_neighbor (3, neighbor_cells[3][1]);
	    new_hexes[5]->set_neighbor (4, new_hexes[6]);
	    new_hexes[5]->set_neighbor (5, new_hexes[4]);

	    new_hexes[6]->set_neighbor (0, new_hexes[2]);
	    new_hexes[6]->set_neighbor (1, neighbor_cells[1][2]);
	    new_hexes[6]->set_neighbor (2, new_hexes[5]);
	    new_hexes[6]->set_neighbor (3, neighbor_cells[3][2]);
	    new_hexes[6]->set_neighbor (4, neighbor_cells[4][2]);
	    new_hexes[6]->set_neighbor (5, new_hexes[7]);

	    new_hexes[7]->set_neighbor (0, new_hexes[3]);
	    new_hexes[7]->set_neighbor (1, neighbor_cells[1][3]);
	    new_hexes[7]->set_neighbor (2, new_hexes[4]);
	    new_hexes[7]->set_neighbor (3, new_hexes[6]);
	    new_hexes[7]->set_neighbor (4, neighbor_cells[4][3]);
	    new_hexes[7]->set_neighbor (5, neighbor_cells[5][2]);


					     // now we need to set the neighbors
					     // neighborship information; this
					     // is only necessary if the
					     // neighboring cell is refined,
					     // i.e. is on the same level as
					     // the new children of the
					     // present cell
	    for (unsigned int nb=0; nb<GeometryInfo<dim>::faces_per_cell; ++nb)
	      for (unsigned int subface=0;
		   subface<GeometryInfo<dim>::subfaces_per_face; ++subface)
		if ((neighbor_cells[nb][subface].state() == valid) &&
		    (neighbor_cells[nb][subface]->level() ==
		     hex->level()+1))
		  {
						     // ok, the neighbor is
						     // a refined one and we
						     // need to set one of
						     // the new children as
						     // its neighbor
		    cell_iterator neighbor = neighbor_cells[nb][subface];

						     // find which neighbor
						     // pointer is to be reset;
						     // this pointer still
						     // points to the present
						     // cell
		    unsigned int face;
		    for (face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
		      if (neighbor->neighbor(face) == hex)
			break;

		    Assert (face<GeometryInfo<dim>::faces_per_cell,
			    ExcInternalError());

						     // now neighbor->neighbor(face)
						     // needs to be reset. do a
						     // large switch statement to
						     // find out which of the
						     // eight children is next
						     // to the face/subface we
						     // are presently at
		    static const unsigned int child_on_face_and_subface[6][4]
		      = {  
			    {0, 1, 2, 3},
			    {4, 5, 6, 7},
			    {0, 1, 5, 4},
			    {1, 5, 6, 2},
			    {3, 2, 6, 7},
			    {0, 4, 7, 3}
		      };
		    neighbor->set_neighbor(face,
					   new_hexes[child_on_face_and_subface[nb][subface]]);
		  };

		    
	    
					     // note that the refinement
					     // flag was already cleared at
					     // the beginning of this loop
	  };
    };

				   // re-compute number of lines and quads
  update_number_cache ();


#ifdef DEBUG
  for (unsigned int level=0; level<levels.size(); ++level) 
    levels[level]->monitor_memory (3);

				   // check whether really all refinement flags
				   // are reset (also of previously non-active
				   // cells which we may not have touched. If
				   // the refinement flag of a non-active cell
				   // is set, something went wrong since the
				   // cell-accessors should have caught this)
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
};



#endif




template <int dim>
void Triangulation<dim>::execute_coarsening () {
  				   // loop over all cells. Flag all cells of
				   // which all children are flagged for
				   // coarsening and delete the childrens'
				   // flags. In effect, only those
				   // cells are flagged of which originally
				   // all children were flagged and for which
				   // all children are on the same refinement
				   // level. For flagging, the user flags are
				   // used, to avoid confusion and because
				   // non-active cells can't be flagged for
				   // coarsening. Note that because of the
				   // effects of #prepare_coarsening#, of a
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
	    };
	};
  

				   // now do the actual coarsening step. Since
				   // the loop goes over used cells only we need
				   // not worry about deleting some cells since
				   // the ++operator will then just hop over
				   // them if we should hit one. Do the loop
				   // in the reverse way since we may only
				   // delete some cells if their neighbors
				   // have already been deleted (if the latter
				   // are on a higher level for example)
				   //
				   // if there is only one level, there can not
				   // be anything to do
  if (levels.size() >= 2)
    for (cell = last(levels.size()-2); cell!=endc; --cell)
      if (cell->user_flag_set())
					 // use a separate function, since this
					 // is dimension specific
	delete_children (cell);

  				   // re-compute number of lines and quads
  update_number_cache ();

  				   // in principle no user flags should be
				   // set any more at this point
#if DEBUG
  for (cell=begin(); cell!=endc; ++cell)
    Assert (cell->user_flag_set() == false, ExcInternalError());
#endif
};





template <int dim>
void Triangulation<dim>::prepare_refinement_dim_dependent () 
{};



#if deal_II_dimension == 3

template <>
void Triangulation<3>::prepare_refinement_dim_dependent () 
{
  const unsigned int dim = 3;
  
				   // first clear flags on lines, since
				   // we need them to determine which
				   // lines will be refined
  for (line_iterator line=begin_line(); line!=end_line(); ++line)
    line->clear_user_flag();

				   // variables to store whether the mesh was
				   // changed in the present loop and in the
				   // whole process
  bool mesh_changed      = false;

  do
    {
      mesh_changed = false;

				       // flag those lines that will be refined
      for (active_cell_iterator cell=begin_active(); cell!=end(); ++cell)
	if (cell->refine_flag_set())
	  for (unsigned int line=0; line<GeometryInfo<dim>::lines_per_cell; ++line)
					     // if the line is not yet refined,
					     // it will be in the process
	    if (!cell->line(line)->has_children())
	      cell->line(line)->set_user_flag();


				       // now check whether there are cells with
				       // lines that are more than once refined
				       // or that will be more than once
				       // refined. The first thing should never
				       // be the case, in the second case we
				       // flag the cell for refinement
      for (active_cell_iterator cell=last_active(); cell!=end(); --cell)
	for (unsigned int line=0; line<GeometryInfo<dim>::lines_per_cell; ++line)
	  {
	    if (cell->line(line)->has_children())
	      {
						 // if this line is refined, its
						 // children should not have
						 // further children
						 //
						 // however, if any of the
						 // children is flagged for
						 // further refinement, we need
						 // to refine this cell also
						 // (at least, if the cell is
						 // not already flagged)
		bool offending_line_found = false;
		
		for (unsigned int c=0; c<2; ++c)
		  {
		    Assert (cell->line(line)->child(c)->has_children() == false,
			    ExcInternalError());
		    
		    if (cell->line(line)->child(c)->user_flag_set () &&
			!cell->refine_flag_set())
		      {
							 // tag this cell for
							 // refinement
			cell->clear_coarsen_flag ();
			cell->set_refine_flag();
			
							 // note that we have
							 // changed the grid
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
		      };
		  };

		if (offending_line_found)
		  {
		    mesh_changed = true;
		    break;
		  };
	      };

						 // there is another thing here:
						 // if any of the lines if refined,
						 // we may not coarsen this cell.
						 // this also holds true if the
						 // line is not yet refined, but
						 // will be
						 //
						 // this is not totally true,
						 // since the neighbors' children
						 // may also be all coarsened,
						 // but we do not catch these
						 // aspects here; in effect, we
						 // disallow to coarsen sharp
						 // edges where the refinement
						 // level decreases from each cell
						 // to the next
	    if (cell->line(line)->has_children() ||
		cell->line(line)->user_flag_set())
	      if (cell->coarsen_flag_set())
		{
		  cell->clear_coarsen_flag ();
		  mesh_changed = true;
		  
		  break;
		};
	  };
    }
  while (mesh_changed == true);
};

#endif



template <int dim>
void Triangulation<dim>::fix_coarsen_flags () {
				   // loop over all cells. Flag all cells of
				   // which all children are flagged for
				   // coarsening and delete the childrens'
				   // flags. Also delete all flags of cells
				   // for which not all children of a cell
				   // are flagged. In effect, only those
				   // cells are flagged of which originally
				   // all children were flagged and for which
				   // all children are on the same refinement
				   // level. For flagging, the user flags are
				   // used, to avoid confusion and because
				   // non-active cells can't be flagged for
				   // coarsening
				   //
				   // In effect, all coarsen flags are turned
				   // into user flags of the mother cell if
				   // coarsening is possible or deleted
				   // otherwise. Coarsen flags of cells with
				   // no mother cell, i.e. on the
				   // coarsest level are deleted explicitely.
  clear_user_flags ();
				   // number of active children of #cell#.
				   // number of children of #cell# which are
				   // flagged for coarsening
  unsigned int flagged_children;
      
  cell_iterator cell = begin(),
		endc = end();
  for (; cell!=endc; ++cell) 
    {
				       // nothing to do if we are already on
				       // the finest level; if we are on the
				       // coarsest level, delete coarsen flag
				       // since no coarsening possible
      if (cell->active()) 
	{
	  if (cell->level() == 0)
	    cell->clear_coarsen_flag();
	  continue;
	};
	  
      flagged_children = 0;
      for (unsigned int child=0; child<GeometryInfo<dim>::children_per_cell; ++child)
	if (cell->child(child)->active() &&
	    cell->child(child)->coarsen_flag_set()) 
	  {
	    ++flagged_children;
					     // clear flag since we don't need
					     // it anymore
	    cell->child(child)->clear_coarsen_flag();
	  };
	  
				       // flag this cell for coarsening if all
				       // children were flagged
      if (flagged_children == GeometryInfo<dim>::children_per_cell)
	cell->set_user_flag();
    };
      
				   // in principle no coarsen flags should be
				   // set any more at this point
#if DEBUG
  for (cell=begin(); cell!=endc; ++cell)
    Assert (cell->coarsen_flag_set() == false, ExcInternalError());
#endif

				   // revert change of flags: use coarsen
				   // flags again and delete user flags
  for (cell=last(); cell!=endc; --cell)
    if (cell->user_flag_set())
      {
	cell->clear_user_flag();

					 // find out whether the children
					 // of this cell may be flagged
					 // for refinement
	bool coarsening_allowed = true;
	for (unsigned int c=0; c<GeometryInfo<dim>::children_per_cell; ++c)
	  for (unsigned int n=0; n<GeometryInfo<dim>::faces_per_cell; ++n)
	    {
	      const cell_iterator child_neighbor = cell->child(c)->neighbor(n);
	      if (child_neighbor.state() == valid)
		if (child_neighbor->has_children() ||
						     // neighbor has children,
						     // then only allow coarsening
						     // if this neighbor will be
						     // coarsened as well. however,
						     // I don't see the consequences
						     // of this case, so simply
						     // disallow coarsening if
						     // any neighbor is more
						     // refined. maybe someone
						     // else will want to do this
						     // some time
		    child_neighbor->refine_flag_set())
		  coarsening_allowed = false;
	    };
	
					 // if allowed: tag the children
					 // for coarsening
	if (coarsening_allowed)
	  for (unsigned int c=0; c<GeometryInfo<dim>::children_per_cell; ++c)
	    {
	      Assert (cell->child(c)->refine_flag_set()==false,
		    ExcInternalError());
	      
	      cell->child(c)->set_coarsen_flag();
	    };
      };
};
  



template <int dim>
bool Triangulation<dim>::prepare_coarsening_and_refinement () {
      
				   // save the flags to determine whether
				   // something was changed in the course
				   // of this function
  vector<bool> flags_before[2];
  save_coarsen_flags (flags_before[0]);
  save_refine_flags (flags_before[1]);


				   // do nothing in 1d, except setting the
				   // coarsening flags correctly
  if (dim == 1)
    {
      fix_coarsen_flags ();

      vector<bool> flags_after[2];
      save_coarsen_flags (flags_after[0]);
      save_refine_flags (flags_after[1]);

      return ((flags_before[0] != flags_after[0]) ||
	      (flags_before[1] != flags_after[1]));
    };

				   // for all other dimensions

				   // save the flags at the outset of
				   // each loop. we do so in order to
				   // find out whether something was
				   // changed in the present loop, in
				   // which case we would have to re-run
				   // the loop. the other possibility to
				   // find this out would be to set a
				   // flag #something_changed# to true
				   // each time we change something.
				   // however, sometimes one change in
				   // one of the parts of the loop is
				   // undone by another one, so we might
				   // end up in an endless loop. we
				   // could be tempted to break this loop
				   // at an arbitrary number of runs,
				   // but that would not be a clean
				   // solution, since we would either have
				   // to
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
  vector<bool> flags_before_loop[2] = {flags_before[0],
				       flags_before[1]};

				   // now for what is done in each loop: we
				   // have to fulfill several tasks at the
				   // same time, namely several mesh smoothing
				   // algorithms and mesh regularisation, by
				   // which we mean that the next mesh fulfills
				   // several requirements such as no double
				   // refinement at each face or line, etc.
				   //
				   // since doing these things at once seems
				   // almost impossible (in the first year of
				   // this library, they were done in two
				   // functions, one for refinement and
				   // one for coarsening, and most things
				   // within these were done at once, so the
				   // code was rather impossible to join into
				   // this, only, function), we do them one
				   // after each other. the order in which
				   // we do them is such that the important
				   // tasks, namely regularisation, are done
				   // last and the least important things
				   // are done the first. the following
				   // order is chosen:
				   //
				   // 1/ eliminate refined islands in the
				   //    interior and at the boundary. since
				   //    they don't do much harm besides
				   //    increasing the number of degrees of
				   //    freedom, doing this has a rather low
				   //    priority.
  				   // 2/ limit the level difference of
				   //    neighboring cells at each vertex.
  				   // 3/ eliminate unrefined islands. this
				   //    has higher priority since this
				   //    diminishes the approximation
				   //    properties not only of the unrefined
				   //    island, but also of the surrounding
				   //    patch.
  				   // 4/ take care of the requirement that no
				   //    double refinement is done at each face
  				   // 5/ take care that no double refinement
				   //    is done at each line in 3d or higher
				   //    dimensions.
				   // 6/ make sure that all children of each
				   //    cell are either flagged for coarsening
				   //    or none of the children is
				   //
				   // For some of these steps, it is known that
				   // they interact. Namely, it is not possible
				   // to guarantee that after step 5 another
				   // step 4 would have no effect; the same
				   // holds for the opposite order and also
				   // when taking into account step 6. however,
				   // it is important to guarantee that step
				   // five or six do not undo something that
				   // step 4 did, and step 6 not something of
				   // step 5, otherwise the requirements will
				   // not be satisfied even if the loop
				   // terminates. this is accomplished by
				   // the fact that steps 4 and 5 only *add*
				   // refinement flags and delete coarsening
				   // flags (therefore, step 5
				   // can't undo something that step 4 already
				   // did), and step 6 only deletes coarsening
				   // flags, never adds some. step 6 needs also
				   // take care that it won't tag cells for
				   // refinement for which some neighbors
				   // are more refined or will be refined.
  bool mesh_changed_in_this_loop = false;
  do
    {
				       //////////////////////////////////////
				       // STEP 1:
				       //    eliminate refined islands in the
				       //    interior and at the boundary. since
				       //    they don't do much harm besides
				       //    increasing the number of degrees of
				       //    freedom, doing this has a rather low
				       //    priority.
      if (smooth_grid & (eliminate_refined_inner_islands |
			 eliminate_refined_boundary_islands)) 
	{
	  cell_iterator       cell;
	  const cell_iterator endc = end();
	  
	  for (cell=begin(); cell!=endc; ++cell)
	    if (!cell->active() || (cell->active() && cell->refine_flag_set()))
	      {
						 // check whether all children are
						 // active, i.e. not refined
						 // themselves. This is a precondition
						 // that the children may be coarsened
						 // away. If the cell is only flagged
						 // for refinement, then all future
						 // children will be active
		bool all_children_active = true;
		if (!cell->active())
		  for (unsigned int c=0; c<GeometryInfo<dim>::children_per_cell; ++c)
		    if (!cell->child(c)->active()) 
		      {
			all_children_active = false;
			break;
		      };

		if (all_children_active) 
		  {
						     // count number of refined and
						     // unrefined neighbors of cell.
						     // neighbors on lower levels
						     // are counted as unrefined since
						     // they can only get to the same
						     // level as this cell by the
						     // next refinement cycle
		    unsigned int unrefined_neighbors = 0,
				     total_neighbors = 0;

		    for (unsigned int n=0; n<GeometryInfo<dim>::faces_per_cell; ++n) 
		      {
			const cell_iterator neighbor = cell->neighbor(n);
			if (neighbor.state() == valid)
			  ++total_neighbors;

			if (neighbor.state() == valid)
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
			      };
		      };
		
						     // if all neighbors unrefined:
						     // mark this
						     // cell for coarsening or don't
						     // refine if marked for that
						     //
						     // also do the distinction
						     // between the two versions of
						     // the eliminate_refined_*_islands
						     // flag
		    if ((unrefined_neighbors == total_neighbors) &&
			(((unrefined_neighbors==GeometryInfo<dim>::faces_per_cell) &&
			  (smooth_grid & eliminate_refined_inner_islands)) ||
			 ((unrefined_neighbors<GeometryInfo<dim>::faces_per_cell) &&
			  (smooth_grid & eliminate_refined_boundary_islands)) ))
		      if (!cell->active())
			for (unsigned int c=0;
			     c<GeometryInfo<dim>::children_per_cell; ++c)
			  {
			    cell->child(c)->clear_refine_flag ();
			    cell->child(c)->set_coarsen_flag ();
			  }
		      else 
			cell->clear_refine_flag();
		  };
	      };
	};

				       //////////////////////////////////////
				       // STEP 2:
				       //    limit the level difference of
				       //    neighboring cells at each vertex.
      if (smooth_grid & limit_level_difference_at_vertices) 
	{
					   // store highest level one of the cells
					   // adjacent to a vertex belongs to
	  vector<int> vertex_level (vertices.size(), 0);
	  active_cell_iterator cell = begin_active(),
			       endc = end();
	  for (; cell!=endc; ++cell)
	    for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell;
		 ++vertex)
	      if (cell->refine_flag_set())
		vertex_level[cell->vertex_index(vertex)]
		  = max (vertex_level[cell->vertex_index(vertex)],
			 cell->level()+1);
	      else
		vertex_level[cell->vertex_index(vertex)]
		  = max (vertex_level[cell->vertex_index(vertex)],
			 cell->level());

					   // loop over all cells in reverse
					   // order. do so because we can then
					   // update the vertex levels on the
					   // and maybe already flag additional
					   // cells in this loop
	  for (cell=last_active(); cell != endc; --cell)
	    if (cell->refine_flag_set() == false) 
	      for (unsigned int vertex=0;
		   vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
		if (vertex_level[cell->vertex_index(vertex)] >
		    cell->level()+1)
		  {
						     // refine cell and
						     // update vertex levels
		    cell->clear_coarsen_flag();
		    cell->set_refine_flag();

		    for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell;
			 ++v)
		      vertex_level[cell->vertex_index(v)]
			= max (vertex_level[cell->vertex_index(v)],
			       cell->level()+1);

						     // now that we fixed this cell,
						     // we can safely leave this
						     // inner loop.
		    break;
		  };	  
	};

				       /////////////////////////////////////
				       // STEP 3:
				       //    eliminate unrefined islands. this
				       //    has higher priority since this
				       //    diminishes the approximation
				       //    properties not only of the unrefined
				       //    island, but also of the surrounding
				       //    patch.
      if (smooth_grid & eliminate_unrefined_islands)
	{
	  active_cell_iterator cell = begin_active(),
			       endc = end();
	  for (; cell!=endc; ++cell) 
	    {
					       // if cell is already flagged
					       // for refinement: nothing to
					       // do anymore
	      if (cell->refine_flag_set())
		continue;
		  
	      unsigned int refined_neighbors = 0,
			 unrefined_neighbors = 0;
	      for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
		if (!cell->at_boundary(face))
		  {
						     // neighbor may only be on
						     // the same level or one
						     // level below because of
						     // the regularisation above
		    Assert ((cell->neighbor_level(face)==cell->level()) ||
			    (cell->neighbor_level(face)==cell->level()-1),
			    ExcInternalError());
		    if ((cell->neighbor_level(face) == cell->level()) &&
			(cell->neighbor(face)->refine_flag_set() ||
			 cell->neighbor(face)->has_children()))
		      ++refined_neighbors;
		    else
		      ++unrefined_neighbors;
		  };

	      if (unrefined_neighbors < refined_neighbors)
		{
		  if (cell->coarsen_flag_set())
		    cell->clear_coarsen_flag();
		  cell->set_refine_flag ();
		};
	    };
	};

				       /////////////////////////////////
				       // STEP 4:
				       //    take care of the requirement that no
				       //    double refinement is done at each face
      for (active_cell_iterator cell = last_active(); cell != end(); --cell)
	if (cell->refine_flag_set() == true) 
	  {
					     // loop over neighbors of cell
	    for (unsigned int i=0; i<GeometryInfo<dim>::faces_per_cell; ++i)
	      if (cell->neighbor(i).state() == valid)
		{
						   // regularisation?
		  if ((cell->neighbor_level(i) == cell->level()-1)
		      &&
		      (cell->neighbor(i)->refine_flag_set() == false))
		    {
		      if (cell->neighbor(i)->coarsen_flag_set())
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
		};
	  };
      
				       //////////////////////////////////////
				       // STEP 5:
				       //    take care that no double refinement
				       //    is done at each line in 3d or higher
				       //    dimensions.
      prepare_refinement_dim_dependent ();
      
				       //////////////////////////////////////
				       // STEP 6:
				       //    make sure that all children of each
				       //    cell are either flagged for coarsening
				       //    or none of the children is
      fix_coarsen_flags ();

				       // get the refinement and coarsening
				       // flags
      vector<bool> flags_after_loop[2];
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
				   // #flags_before_loop# represents the
				   // state after the last loop, i.e.
				   // the present state
  return ((flags_before[0] != flags_before_loop[0]) ||
	  (flags_before[1] != flags_before_loop[1]));
};



#if deal_II_dimension == 1

template <>
void Triangulation<1>::delete_children (cell_iterator &cell) {
  const unsigned int dim=1;
				   // first we need to reset the
				   // neighbor pointers of the neighbors
				   // of this cell's children to this
				   // cell. This is different for one
				   // dimension, since there neighbors
				   // can have a refinement level
				   // differing from that of this cell's
				   // children by more than one level.

  Assert (!cell->child(0)->has_children() && !cell->child(1)->has_children(),
	  ExcInternalError());
  
				   // first do it for the cells to the
				   // left
  if (cell->neighbor(0).state() == valid)
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
	    
					     // move on to further children
					     // on the boundary between this
					     // cell and its neighbor
	    if (neighbor->has_children())
	      neighbor = neighbor->child(1);
	    else
	      break;
	  };
      };

  				   // now do it for the cells to the
				   // left
  if (cell->neighbor(1).state() == valid)
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
	    
					     // move on to further children
					     // on the boundary between this
					     // cell and its neighbor
	    if (neighbor->has_children())
	      neighbor = neighbor->child(0);
	    else
	      break;
	  };
      };

	
				   // delete the vertex which will not be
				   // needed anymore. This vertex is the
				   // second of the second line of the
				   // first child
  vertices_used[cell->child(0)->vertex_index(1)] = false;

				   // invalidate children.
				   // clear user pointers, to avoid that
				   // they may appear at unwanted places
				   // later on...
  for (unsigned int child=0; child<GeometryInfo<dim>::children_per_cell; ++child)
    {
      cell->child(child)->clear_user_pointer();
      cell->child(child)->clear_user_flag();
      cell->child(child)->clear_used_flag();
    };
  
  
				   // delete pointer to children
  cell->set_children (-1);
  cell->clear_user_flag();
};

#endif


#if deal_II_dimension == 2

template <>
void Triangulation<2>::delete_children (cell_iterator &cell) {
  const unsigned int dim=2;
				   // first we need to reset the
				   // neighbor pointers of the neighbors
				   // of this cell's children to this
				   // cell. This is different for one
				   // dimension, since there neighbors
				   // can have a refinement level
				   // differing from that of this cell's
				   // children by more than one level.
				   // For two or more dimensions, the
				   // neighbors of the children may only
				   // be on the same level or on the level
				   // of this cell (the case that the
				   // neighbors are more refined than
				   // the children was eliminated in
				   // #prepare_coarsening#
  for (unsigned int child=0; child<GeometryInfo<dim>::children_per_cell; ++child)
    for (unsigned int n=0; n<GeometryInfo<dim>::faces_per_cell;
	 ++n)
      {
	const cell_iterator neighbor = cell->child(child)->neighbor(n);
					 // do nothing if at boundary
	if (neighbor.state() != valid)
	  continue;
	
	Assert ((neighbor->level()==cell->level()) ||
		(neighbor->level()==cell->level()+1),
		ExcInternalError());
	
					 // if the neighbor's level is the
					 // same as that of #cell#, then
					 // it's neighbor pointers points
					 // to this cell rather than to
					 // this cell's child. In that
					 // case we need not do anything.
					 // If the neighbor is refined
					 // as often as are the children,
					 // we need to reset those neigbor
					 // pointers that point to the
					 // child of this cell; when
					 // resetting the neighbor
					 // pointers of neighbors of one
					 // of the children, we will also
					 // reset the neighbor pointers
					 // other children to the present
					 // cell, but this does no harm
					 // since we delete the children
					 // afterwards anyway
	if (neighbor->level() == cell->level()+1)
	  for (unsigned int neighbor_neighbor=0;
	       neighbor_neighbor<GeometryInfo<dim>::faces_per_cell;
	       ++neighbor_neighbor)
	    if (neighbor->neighbor(neighbor_neighbor) == cell->child(child))
	      neighbor->set_neighbor(neighbor_neighbor, cell);
      };

				   // delete the vertex which will not be
				   // needed anymore. This vertex is the
				   // second of the second line of the
				   // first child
  vertices_used[cell->child(0)->line(1)->vertex_index(1)] = false;

				   // clear user pointers, to avoid that
				   // they may appear at unwanted places
				   // later on...
  cell->child(0)->line(1)->clear_user_pointer();
  cell->child(0)->line(2)->clear_user_pointer();
  cell->child(2)->line(0)->clear_user_pointer();
  cell->child(2)->line(3)->clear_user_pointer();
  
				   // same for user flags
  cell->child(0)->line(1)->clear_user_flag();
  cell->child(0)->line(2)->clear_user_flag();
  cell->child(2)->line(0)->clear_user_flag();
  cell->child(2)->line(3)->clear_user_flag();
  
				   // delete the four interior lines
  cell->child(0)->line(1)->clear_used_flag();
  cell->child(0)->line(2)->clear_used_flag();
  cell->child(2)->line(0)->clear_used_flag();
  cell->child(2)->line(3)->clear_used_flag();

				   // for the four faces: if the neighbor
				   // does not itself need the subfaces,
				   // delete them. note that since dim>1
				   // the level of a neighbor is either
				   // one less or the same as that of
				   // cell
  for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face) 
    if ((cell->neighbor(face).state() != valid) ||
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
	  };
	
	cell->face(face)->clear_children();
      };
  
				   // invalidate children
  for (unsigned int child=0; child<GeometryInfo<dim>::children_per_cell; ++child)
    {
      cell->child(child)->clear_user_pointer();
      cell->child(child)->clear_user_flag();
      cell->child(child)->clear_used_flag();
    };
  
  
				   // delete pointer to children
  cell->set_children (-1);
  cell->clear_user_flag();
};

#endif



#if deal_II_dimension == 3


template <>
void Triangulation<3>::delete_children (cell_iterator &cell) {
  const unsigned int dim=3;
				   // first we need to reset the
				   // neighbor pointers of the neighbors
				   // of this cell's children to this
				   // cell. This is different for one
				   // dimension, since there neighbors
				   // can have a refinement level
				   // differing from that of this cell's
				   // children by more than one level.
				   // For two or more dimensions, the
				   // neighbors of the children may only
				   // be on the same level or on the level
				   // of this cell (the case that the
				   // neighbors are more refined than
				   // the children was eliminated in
				   // #prepare_coarsening#
  for (unsigned int child=0; child<GeometryInfo<dim>::children_per_cell; ++child)
    for (unsigned int n=0; n<GeometryInfo<dim>::faces_per_cell; ++n)
      {
	const cell_iterator neighbor = cell->child(child)->neighbor(n);
					 // do nothing if at boundary
	if (neighbor.state() != valid)
	  continue;
	
	Assert ((neighbor->level()==cell->level()) ||
		(neighbor->level()==cell->level()+1),
		ExcInternalError());
	
					 // if the neighbor's level is the
					 // same as that of #cell#, then
					 // it's neighbor pointers points
					 // to this cell rather than to
					 // this cell's child. In that
					 // case we need not do anything.
					 // If the neighbor is refined
					 // as often as are the children,
					 // we need to reset those neigbor
					 // pointers that point to the
					 // child of this cell; when
					 // resetting the neighbor
					 // pointers of neighbors of one
					 // of the children, we will also
					 // reset the neighbor pointers
					 // other children to the present
					 // cell, but this does no harm
					 // since we delete the children
					 // afterwards anyway
	if (neighbor->level() == cell->level()+1)
	  for (unsigned int neighbor_neighbor=0;
	       neighbor_neighbor<GeometryInfo<dim>::faces_per_cell;
	       ++neighbor_neighbor)
	    if (neighbor->neighbor(neighbor_neighbor) == cell->child(child))
	      neighbor->set_neighbor(neighbor_neighbor, cell);
      };

				   // delete the vertex which will not be
				   // needed anymore. This vertex is the
				   // vertex at the heart of this cell,
				   // which is the sixth of the 
				   // first child
  vertices_used[cell->child(0)->vertex_index(6)] = false;

				   ///////////////////////////////////////
				   // delete interior quads and lines
				   //
				   // first set up a list of these
				   // line's and quad's indices
  const quad_iterator interior_quads[12]
    = {  cell->child(0)->face(1),
	 cell->child(1)->face(1),
	 cell->child(2)->face(1),
	 cell->child(3)->face(1),

	 cell->child(0)->face(3),
	 cell->child(3)->face(3),
	 cell->child(4)->face(3),
	 cell->child(7)->face(3),

	 cell->child(0)->face(4),
	 cell->child(1)->face(4),
	 cell->child(4)->face(4),
	 cell->child(5)->face(4) 
    };

  const line_iterator interior_lines[6]
    = {  cell->child(0)->line(10),
	 cell->child(4)->line(10),

	 cell->child(0)->line(6),
	 cell->child(1)->line(6),

	 cell->child(0)->line(5),
	 cell->child(3)->line(5) 
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
    };

  for (unsigned int l=0; l<6; ++l)
    {
      interior_lines[l]->clear_user_pointer();
      interior_lines[l]->clear_user_flag();
      interior_lines[l]->clear_used_flag();
    };

				   // for the six faces: if the neighbor
				   // does not itself need the subfaces,
				   // delete them. note that since dim>1
				   // the level of a neighbor is either
				   // one less or the same as that of
				   // cell
  for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face) 
    if ((cell->neighbor(face).state() != valid) ||
	(cell->neighbor(face)->level() == cell->level()-1) ||
	((cell->neighbor(face)->level() == cell->level()) &&
	 !cell->neighbor(face)->has_children()))
      {
	quad_iterator quad = cell->face(face);
	
					 // delete middle vertex
	vertices_used[quad->child(0)->vertex_index(2)] = false;

	const line_iterator interior_lines[4]
	  = { quad->child(0)->line(1),
	      quad->child(0)->line(2),
	      quad->child(2)->line(0),
	      quad->child(2)->line(3) 
	  };

					 // delete interior lines
	for (unsigned int l=0; l<4; ++l)
	  {
	    interior_lines[l]->clear_user_pointer ();
	    interior_lines[l]->clear_user_flag ();
	    interior_lines[l]->clear_used_flag ();
	  };
	
					 // delete the four subfaces
	for (unsigned int subface=0;
	     subface<GeometryInfo<dim>::subfaces_per_face; ++subface)
	  {
	    quad->child(subface)->clear_user_pointer ();
	    quad->child(subface)->clear_user_flag ();
	    quad->child(subface)->clear_used_flag ();
	  };
	
	quad->clear_children();
      };

				   // invalidate children
  for (unsigned int child=0; child<GeometryInfo<dim>::children_per_cell; ++child)
    {
      cell->child(child)->clear_user_pointer();
      cell->child(child)->clear_user_flag();
      cell->child(child)->clear_used_flag();
    };
  
  
				   // delete pointer to children
  cell->clear_children ();
  cell->clear_user_flag();

				   // now there still are the 12 lines of
				   // this hex which are refined and which
				   // may need coarsening. however, it is not
				   // so easy to decide whether they are still
				   // needed, since it does not suffice to
				   // ask the neighbors. we also need to ask
				   // those cells, which are "around the
				   // corner".
				   //
				   // to do so: first set up a list of 12 pairs
				   // of line_iterators and flags which denote
				   // whether the line's children are still
				   // needed
				   //
				   // we default to: "is not needed" because in
				   // this case, if we make an error in the code
				   // below, some lines will be deleted that in
				   // fact are needed. this will eventually be
				   // caught somewhen, because there are many
				   // checks whether an iterator points to
				   // something used. the opposite case, that
				   // we do not delete lines that are no more
				   // used, is more severe and causes a memory
				   // leak which is probably impossible to
				   // find.
  const pair<line_iterator,bool> line_is_needed_pairs[12]
    = { make_pair(cell->line(0), false),
	make_pair(cell->line(1), false),
	make_pair(cell->line(2), false),
	make_pair(cell->line(3), false),
	make_pair(cell->line(4), false),
	make_pair(cell->line(5), false),
	make_pair(cell->line(6), false),
	make_pair(cell->line(7), false),
	make_pair(cell->line(8), false),
	make_pair(cell->line(9), false),
	make_pair(cell->line(10), false),
	make_pair(cell->line(11), false)  };
				   // if in debug mode: make sure that
				   // none of the lines of this cell is
				   // twice refined; else, deleting this
				   // cell's children will result
				   // in an invalid state
  for (unsigned int line=0; line<12; ++line)
    for (unsigned int c=0; c<2; ++c)
      Assert (!cell->line(line)->child(c)->has_children(),
	      ExcInternalError());

				   // next make a map out of this for simpler
				   // access to the flag associated with a
				   // line
  map<line_iterator,bool> line_is_needed (&line_is_needed_pairs[0],
					  &line_is_needed_pairs[12]);
  
				   // then ask each neighbor and their neighbors
  for (unsigned int nb=0; nb<GeometryInfo<dim>::faces_per_cell; ++nb) 
    {
      const cell_iterator neighbor = cell->neighbor(nb);
					 // do nothing if at boundary
      if (neighbor.state() != valid)
	continue;

      Assert (neighbor->level() == cell->level(),
	      ExcInternalError());

				       // if the neighbor itself has
				       // children, then the four lines
				       // of the common face are definitely
				       // needed
      if (neighbor->has_children())
	{
	  for (unsigned int i=0; i<4; ++i)
	    {
	      Assert (line_is_needed.find(cell->face(nb)->line(i))
		      != line_is_needed.end(),
		      ExcInternalError());
	      
	      line_is_needed[cell->face(nb)->line(i)] = true;
	    };
					   // ok for this neighbor
	  continue;
	};

				       // if the neighbor is not refined,
				       // then it may still be that one of
				       // his neighbors may need one of our
				       // lines. the present cell is also
				       // one of his neighbors, but this one
				       // is not any more refined
      for (unsigned int nb_nb=0; nb_nb<GeometryInfo<dim>::faces_per_cell;
	   ++nb_nb)
	{
	  const cell_iterator neighbor_neighbor = neighbor->neighbor(nb_nb);
					   // do nothing if at boundary
	  if (neighbor_neighbor.state() != valid)
	    continue;

	  Assert ((neighbor_neighbor->level() == cell->level()) ||
		  (neighbor_neighbor->level() == cell->level()-1),
		  ExcInternalError());
	  
					   // also do nothing if the neighbor
					   // is at a lower level (but then it
					   // should not be refined)
	  if (neighbor_neighbor->level() == cell->level()-1)
	    {
	      Assert (!neighbor_neighbor->has_children(),
		      ExcInternalError());
	      continue;
	    };

					   // neighbor's neighbor is on
					   // same level. if it has children,
					   // then look closer
	  if (neighbor_neighbor->has_children())
					     // if any of those cell's lines
					     // is one of those that we are
					     // interested in, then flag it
	    for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_cell; ++l)
	      {
		line_iterator line = neighbor_neighbor->line(l);

		if (line_is_needed.find(line) != line_is_needed.end())
		  line_is_needed[line] = true;
	      };
	};
    };

  
				   // now, if the lines are not marked as
				   // needed, we may delete their children
				   // and the midpoint
  map<line_iterator,bool>::iterator line_and_flag;
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
	  };

	line->clear_children();
      };
};

#endif



template <int dim>
void Triangulation<dim>::write_bool_vector (const unsigned int  magic_number1,
					    const vector<bool> &v,
					    const unsigned int  magic_number2,
					    ostream            &out) {
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
  out << magic_number1 << ' ' << N << endl;
  for (unsigned int i=0; i<N/8+1; ++i) 
    out << static_cast<unsigned int>(flags[i]) << ' ';
  
  out << endl << magic_number2 << endl;
  
  delete[] flags;

  AssertThrow (out, ExcIO());
};



template <int dim>
void Triangulation<dim>::read_bool_vector (const unsigned int  magic_number1,
					   vector<bool>       &v,
					   const unsigned int  magic_number2,
					   istream            &in) {
  AssertThrow (in, ExcIO());

  unsigned int magic_number;
  in >> magic_number;
  Assert (magic_number==magic_number1, ExcGridReadError());

  unsigned int N;
  in >> N;
  v.resize (N);

  unsigned char *flags = new unsigned char[N/8+1];
  unsigned short int tmp;
  for (unsigned int i=0; i<N/8+1; ++i) 
    {
      in >> tmp;
      flags[i] = tmp;
    };

  for (unsigned int position=0; position!=N; ++position)
    v[position] = (flags[position/8] & (1<<(position%8)));

  in >> magic_number;
  Assert (magic_number==magic_number2, ExcGridReadError());

  delete[] flags;

  AssertThrow (in, ExcIO());
};




template <int dim>
void Triangulation<dim>::update_number_cache_lines () 
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
      if (levels[level]->lines.lines.size() != 0) 
	{
	  line_iterator line = begin_line (level),
			endc = (level == levels.size()-1 ?
				line_iterator(end_line()) :
				begin_line (level+1));
	  for (; line!=endc; ++line)
	    ++number_cache.n_lines_level[level];
	};
      
				       // update total number of lines
      number_cache.n_lines += number_cache.n_lines_level[level];
    };

				   // do the update for the number
				   // of active lines as well
  number_cache.n_active_lines_level.resize (levels.size());
  number_cache.n_active_lines = 0;
  for (unsigned int level=0; level<levels.size(); ++level)
    {
				       // count lines on this level
      number_cache.n_active_lines_level[level] = 0;
      if (levels[level]->lines.lines.size() != 0) 
	{
	  active_line_iterator line = begin_active_line (level),
			       endc = end_active_line (level);
	  for (; line!=endc; ++line)
	    ++number_cache.n_active_lines_level[level];
	};
      
				       // update total number of lines
      number_cache.n_active_lines += number_cache.n_active_lines_level[level];
    };
};



#if deal_II_dimension == 1

template <>
void Triangulation<1>::update_number_cache_quads () 
{
  Assert (false, ExcInternalError());
};

#endif

#if deal_II_dimension >= 2

template <int dim>
void Triangulation<dim>::update_number_cache_quads () 
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
      if (levels[level]->quads.quads.size() != 0) 
	{
	  quad_iterator quad = begin_quad (level),
			endc = (level == levels.size()-1 ?
				quad_iterator(end_quad()) :
				begin_quad (level+1));
	  for (; quad!=endc; ++quad)
	    ++number_cache.n_quads_level[level];
	};
      
				       // update total number of quads
      number_cache.n_quads += number_cache.n_quads_level[level];
    };

				   // do the update for the number
				   // of active quads as well
  number_cache.n_active_quads_level.resize (levels.size());
  number_cache.n_active_quads = 0;
  for (unsigned int level=0; level<levels.size(); ++level)
    {
				       // count quads on this level
      number_cache.n_active_quads_level[level] = 0;
      if (levels[level]->quads.quads.size() != 0) 
	{
	  active_quad_iterator quad = begin_active_quad (level),
			       endc = end_active_quad (level);
	  for (; quad!=endc; ++quad)
	    ++number_cache.n_active_quads_level[level];
	};
      
				       // update total number of quads
      number_cache.n_active_quads += number_cache.n_active_quads_level[level];
    };
};

#endif


#if deal_II_dimension == 1

template <>
void Triangulation<1>::update_number_cache_hexes () 
{
  Assert (false, ExcInternalError());
};

#endif


#if deal_II_dimension == 2

template <>
void Triangulation<2>::update_number_cache_hexes () 
{
  Assert (false, ExcInternalError());
};

#endif

#if deal_II_dimension >= 3

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
      if (levels[level]->hexes.hexes.size() != 0) 
	{
	  hex_iterator hex = begin_hex (level),
		      endc = (level == levels.size()-1 ?
			      hex_iterator(end_hex()) :
			      begin_hex (level+1));
	  for (; hex!=endc; ++hex)
	    ++number_cache.n_hexes_level[level];
	};
      
				       // update total number of hexes
      number_cache.n_hexes += number_cache.n_hexes_level[level];
    };

				   // do the update for the number
				   // of active hexes as well
  number_cache.n_active_hexes_level.resize (levels.size());
  number_cache.n_active_hexes = 0;
  for (unsigned int level=0; level<levels.size(); ++level)
    {
				       // count hexes on this level
      number_cache.n_active_hexes_level[level] = 0;
      if (levels[level]->hexes.hexes.size() != 0) 
	{
	  active_hex_iterator hex = begin_active_hex (level),
			     endc = end_active_hex (level);
	  for (; hex!=endc; ++hex)
	    ++number_cache.n_active_hexes_level[level];
	};
      
				       // update total number of hexes
      number_cache.n_active_hexes += number_cache.n_active_hexes_level[level];
    };
};

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
    };
};





// explicit instantiations
template class Triangulation<deal_II_dimension>;

template void Triangulation<deal_II_dimension>
::refine (const Vector<float> &, const double);

template void Triangulation<deal_II_dimension>
::refine (const Vector<double> &, const double);

template void Triangulation<deal_II_dimension>
::coarsen (const Vector<float> &, const double);

template void Triangulation<deal_II_dimension>
::coarsen (const Vector<double> &, const double);


template void Triangulation<deal_II_dimension>
::refine_and_coarsen_fixed_number (const Vector<double> &,
				   const double         top_fraction,
				   const double         bottom_fraction);

template void Triangulation<deal_II_dimension>
::refine_and_coarsen_fixed_number (const Vector<float> &criteria,
				   const double         top_fraction,
				   const double         bottom_fraction);

template void Triangulation<deal_II_dimension>
::refine_and_coarsen_fixed_fraction (const Vector<double> &criteria,
				     const double         top_fraction,
				     const double         bottom_fraction);

template void Triangulation<deal_II_dimension>
::refine_and_coarsen_fixed_fraction (const Vector<float> &criteria,
				     const double         top_fraction,
				     const double         bottom_fraction);

