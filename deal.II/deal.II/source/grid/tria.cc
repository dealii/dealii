/* $Id$ */

#include <grid/tria_boundary.h>
#include <grid/tria_iterator.h>
#include <grid/tria.h>
#include <basic/magic_numbers.h>
#include <basic/data_io.h>
#include <lac/dvector.h>
#include <iostream>
#include <algorithm>
#include <map>
#include <cmath>




template <int dim>
Triangulation<dim>::Triangulation () {
  static StraightBoundary<dim> default_boundary;
  boundary = &default_boundary;
  
  levels.push_back (new TriangulationLevel<dim>);
};


template <int dim>
Triangulation<dim>::~Triangulation () {
  for (unsigned int i=0; i<levels.size(); ++i)
    delete levels[i];

  levels.resize (0);
};



template <int dim>
void Triangulation<dim>::set_boundary (const Boundary<dim> *boundary_object) {
  boundary = boundary_object;
};




void Triangulation<1>::create_triangulation (const vector<Point<1> >    &v,
					     const vector<CellData<1> > &cells,
					     const SubCellData          &) {
				   // note: since no boundary information
				   // can be given in one dimension, the
				   // #subcelldata# field is ignored.
  
  const unsigned int dim=1;
  
  Assert (vertices.size() == 0, ExcTriangulationNotEmpty());
  Assert (n_lines() == 0, ExcTriangulationNotEmpty());

				   // copy vertices
  vertices = v;
  vertices_used = vector<bool> (v.size(), true);
    
				   // store the indices of the lines which
				   // are adjacent to a given vertex
  vector<vector<int> > lines_at_vertex (v.size());

				   // reserve enough space
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
      
				       // note that this cell
				       // is adjacent to these vertices
      lines_at_vertex[cells[cell].vertices[0]].push_back (cell);
      lines_at_vertex[cells[cell].vertices[1]].push_back (cell);
    };


#ifdef DEBUG
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
	      Assert (false, ExcInternalError());
      };

				   // assert there are no more than two boundary
				   // nodes
  Assert (boundary_nodes == 2, ExcInternalError());
#endif


  				   // update neighborship info
  active_line_iterator line = begin_active_line ();
				   // for all lines
  for (; line!=end(); ++line)
				     // for each of the two vertices
    for (unsigned int vertex=0; vertex<(1<<dim); ++vertex)
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
};



void Triangulation<2>::create_triangulation (const vector<Point<2> >    &v,
					     const vector<CellData<2> > &c,
					     const SubCellData          &subcelldata) {
  const unsigned int dim=2;

  Assert (vertices.size() == 0, ExcTriangulationNotEmpty());
  Assert (n_lines() == 0, ExcTriangulationNotEmpty());
  Assert (n_quads() == 0, ExcTriangulationNotEmpty());

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
#ifdef DEBUG
				       // in debug mode: check whether vertex
				       // indices are valid ones
      for (unsigned int vertex=0; vertex<4; ++vertex)
	Assert ((0<=cells[cell].vertices[vertex]) &&
		(cells[cell].vertices[vertex]<(signed int)vertices.size()),
		ExcInvalidVertexIndex (cell, cells[cell].vertices[vertex], vertices.size()));
#endif
      
      pair<int,int> line_vertices[4] = {   // note the order of the vertices
	    make_pair (cells[cell].vertices[0], cells[cell].vertices[1]),
	    make_pair (cells[cell].vertices[1], cells[cell].vertices[2]),
	    make_pair (cells[cell].vertices[0], cells[cell].vertices[3]),
	    make_pair (cells[cell].vertices[3], cells[cell].vertices[2])};

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
	    line_vertices[0] = make_pair (cells[cell].vertices[0], cells[cell].vertices[1]);
	    line_vertices[1] = make_pair (cells[cell].vertices[1], cells[cell].vertices[2]);
	    line_vertices[2] = make_pair (cells[cell].vertices[0], cells[cell].vertices[3]);
	    line_vertices[3] = make_pair (cells[cell].vertices[3], cells[cell].vertices[2]);
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
	  Assert (needed_lines.find(make_pair(line_vertices[line].second,
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
  

#ifdef DEBUG
				   // in debug mode: check the every vertex has
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
      Assert (* (min_element(vertex_touch_count.begin(),
			     vertex_touch_count.end())) >= 2,
	      ExcGridHasInvalidVertices());
    };
#endif
	
  				   // reserve enough space
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
      Assert ((n_adj_cells >= 1) &&
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
	      Assert (false, ExcLineInexistant(line_vertices.first,
					       line_vertices.second));
	      line = end_line();
	    };
	};

				       // Assert that only exterior lines
				       // are sgiven a boundary indicator
      Assert (line->boundary_indicator() == 0,
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
};



void Triangulation<1>::create_hypercube (const double left,
					 const double right) {
  Assert (vertices.size() == 0, ExcTriangulationNotEmpty());
  Assert (n_lines() == 0, ExcTriangulationNotEmpty());

  const Point<1> vertices[2] = { Point<1>(left), Point<1>(right) };
  vector<CellData<1> > cells (1, CellData<1>());
  cells[0].vertices[0] = 0;
  cells[0].vertices[1] = 1;
  cells[0].material_id = 0;
  
  create_triangulation (vector<Point<1> >(&vertices[0], &vertices[2]),
			cells,
			SubCellData());       // no boundary information
};



  
void Triangulation<2>::create_hypercube (const double left,
					 const double right) {
  Assert (vertices.size() == 0, ExcTriangulationNotEmpty());
  Assert (n_lines() == 0, ExcTriangulationNotEmpty());
  Assert (n_quads() == 0, ExcTriangulationNotEmpty());

  const Point<2> vertices[4] = { Point<2>(left,left),
				 Point<2>(right,left),
				 Point<2>(right,right),
				 Point<2>(left,right)  };
  const int cell_vertices[1][4] = { { 0,1,2,3 } };
  vector<CellData<2> > cells (1, CellData<2>());
  for (unsigned int j=0; j<4; ++j)
    cells[0].vertices[j] = cell_vertices[0][j];
  cells[0].material_id = 0;
  
  create_triangulation (vector<Point<2> >(&vertices[0], &vertices[4]),
			cells,
			SubCellData());       // no boundary information
};



void Triangulation<1>::create_hyper_L (const double, const double) {
  Assert (false, ExcInternalError());
};



void Triangulation<2>::create_hyper_L (const double a, const double b) {
  const unsigned int dim=2;
  const Point<dim> vertices[8] = { Point<dim> (a,a),
				   Point<dim> ((a+b)/2,a),
				   Point<dim> (b,a),
				   Point<dim> (a,(a+b)/2),
				   Point<dim> ((a+b)/2,(a+b)/2),
				   Point<dim> (b,(a+b)/2),
				   Point<dim> (a,b),
				   Point<dim> ((a+b)/2,b)  };
  const int cell_vertices[3][4] = {{0, 1, 4, 3},
				   {1, 2, 5, 4},
				   {3, 4, 7, 6}};

  vector<CellData<2> > cells (3, CellData<2>());
  
  for (unsigned int i=0; i<3; ++i) 
    {
      for (unsigned int j=0; j<4; ++j)
	cells[i].vertices[j] = cell_vertices[i][j];
      cells[i].material_id = 0;
    };
  
  create_triangulation (vector<Point<dim> >(&vertices[0], &vertices[8]),
			cells,
			SubCellData());       // no boundary information
};



void Triangulation<1>::create_hyper_ball (const Point<1>, const double) {
  Assert (false, ExcInternalError());
};



void Triangulation<2>::create_hyper_ball (const Point<2> p, const double radius) {
  const unsigned int dim=2;
  const double a = 0.453;         // equilibrate cell sizes
  const Point<dim> vertices[8] = { p+Point<dim>(-1,-1)*radius/sqrt(2),
				   p+Point<dim>(+1,-1)*radius/sqrt(2),
				   p+Point<dim>(-1,-1)*radius/sqrt(2)*a,
				   p+Point<dim>(+1,-1)*radius/sqrt(2)*a,
				   p+Point<dim>(-1,+1)*radius/sqrt(2)*a,
				   p+Point<dim>(+1,+1)*radius/sqrt(2)*a,
				   p+Point<dim>(-1,+1)*radius/sqrt(2),
				   p+Point<dim>(+1,+1)*radius/sqrt(2) };
  const int cell_vertices[5][4] = {{0, 1, 3, 2},
				   {0, 2, 4, 6},
				   {2, 3, 5, 4},
				   {1, 7, 5, 3},
				   {6, 4, 5, 7}};

  vector<CellData<2> > cells (5, CellData<2>());
  
  for (unsigned int i=0; i<5; ++i) 
    {
      for (unsigned int j=0; j<4; ++j)
	cells[i].vertices[j] = cell_vertices[i][j];
      cells[i].material_id = 0;
    };
  
  create_triangulation (vector<Point<dim> >(&vertices[0], &vertices[8]),
			cells,
			SubCellData());       // no boundary information
};




template <int dim>
void Triangulation<dim>::set_all_refine_flags () {
  active_cell_iterator cell = begin_active(),
		       endc = end();

  for (; cell != endc; ++cell)
    cell->set_refine_flag ();
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
void Triangulation<dim>::save_refine_flags (ostream &out) const {
  unsigned int N = n_active_cells();
  active_cell_iterator cell = begin_active(),
		       endc = end();

  unsigned char *flags = new unsigned char[N/8+1];
  for (unsigned int i=0; i<N/8+1; ++i) flags[i]=0;
  
  for (unsigned int position=0; cell!=endc; ++cell, ++position)
    flags[position/8] |= (cell->refine_flag_set() ? (1<<(position%8)) : 0);

				   // format:
				   // 0. magic number 0xabcc
				   // 1. number of active cells
				   // 2. the flags
				   // 3. magic number 0xabcd
  out << mn_tria_refine_flags_begin << " " << N << endl;
  for (unsigned int i=0; i<N/8+1; ++i) 
    out << static_cast<unsigned int>(flags[i]) << " ";
  
  out << endl;
  out << mn_tria_refine_flags_end << endl;

  delete[] flags;
};



template <int dim>
void Triangulation<dim>::load_refine_flags (istream &in) {
  unsigned int magic_number;
  in >> magic_number;
  Assert (magic_number==mn_tria_refine_flags_begin, ExcGridReadError());

  unsigned int N;
  in >> N;
  Assert (N==n_active_cells(), ExcGridsDoNotMatch(N, n_active_cells()));

  unsigned char *flags = new unsigned char[N/8+1];
  unsigned short int tmp;
  for (unsigned int i=0; i<N/8+1; ++i) 
    {
      in >> tmp;
      flags[i] = tmp;
    };
  
  
  active_cell_iterator cell = begin_active(),
		       endc = end();
  unsigned int position=0;
  for (; cell!=endc; ++cell, ++position)
    if (flags[position/8] & (1<<(position%8)))
      cell->set_refine_flag();
    else
      cell->clear_refine_flag();

  Assert (position==N, ExcGridReadError());

  in >> magic_number;
  Assert (magic_number==mn_tria_refine_flags_end, ExcGridReadError());

  delete[] flags;
};




void Triangulation<1>::clear_user_flags () {
  cell_iterator cell = begin(),
		endc = end();
  for (; cell!=endc; ++cell)
    cell->clear_user_flag ();
};



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




void Triangulation<1>::save_user_flags (ostream &out) const {
  save_user_flags_line (out);
};



void Triangulation<2>::save_user_flags (ostream &out) const {
  save_user_flags_line (out);
  save_user_flags_quad (out);
};



template <int dim>
void Triangulation<dim>::save_user_flags_line (ostream &out) const {
  unsigned int N = n_lines();
  line_iterator line = begin_line(),
		endc = end_line();

  unsigned char *flags = new unsigned char[N/8+1];
  for (unsigned int i=0; i<N/8+1; ++i) flags[i]=0;
  
  for (unsigned int position=0; line!=endc; ++line, ++position)
    flags[position/8] |= (line->user_flag_set() ? (1<<(position%8)) : 0);

				   // format:
				   // 0. magic number 0xabce
				   // 1. number of lines
				   // 2. the flags
				   // 3. magic number 0xabcf
  out << mn_tria_line_user_flags_begin << " " << N << endl;
  for (unsigned int i=0; i<N/8+1; ++i) 
    out << static_cast<unsigned int>(flags[i]) << " ";
  
  out << endl;
  out << mn_tria_line_user_flags_end << endl;

  delete[] flags;
};



template <int dim>
void Triangulation<dim>::load_user_flags_line (istream &in) {
  unsigned int magic_number;
  in >> magic_number;
  Assert (magic_number==mn_tria_line_user_flags_begin, ExcGridReadError());

  unsigned int N;
  in >> N;
  Assert (N==n_lines(), ExcGridsDoNotMatch(N, n_lines()));

  unsigned char *flags = new unsigned char[N/8+1];
  unsigned short int tmp;
  for (unsigned int i=0; i<N/8+1; ++i) 
    {
      in >> tmp;
      flags[i] = tmp;
    };
  
  
  line_iterator line = begin_line(),
		endc = end_line();
  unsigned int position=0;
  for (; line!=endc; ++line, ++position)
    if (flags[position/8] & (1<<(position%8)))
      line->set_user_flag();
    else
      line->clear_user_flag();

  Assert (position==N, ExcGridReadError());

  in >> magic_number;
  Assert (magic_number==mn_tria_line_user_flags_end, ExcGridReadError());

  delete[] flags;
};



void Triangulation<1>::save_user_flags_quad (ostream &) const {
  Assert (false, ExcFunctionNotUseful());
};



template <int dim>
void Triangulation<dim>::save_user_flags_quad (ostream &out) const {
  unsigned int N = n_quads();
  quad_iterator quad = begin_quad(),
		endc = end_quad();

  unsigned char *flags = new unsigned char[N/8+1];
  for (unsigned int i=0; i<N/8+1; ++i) flags[i]=0;
  
  for (unsigned int position=0; quad!=endc; ++quad, ++position)
    flags[position/8] |= (quad->user_flag_set() ? (1<<(position%8)) : 0);

				   // format:
				   // 0. magic number 0xabce
				   // 1. number of quads
				   // 2. the flags
				   // 3. magic number 0xabcf
  out << mn_tria_quad_user_flags_begin << " " << N << endl;
  for (unsigned int i=0; i<N/8+1; ++i) 
    out << static_cast<unsigned int>(flags[i]) << " ";
  
  out << endl;
  out << mn_tria_quad_user_flags_end << endl;

  delete[] flags;
};



void Triangulation<1>::load_user_flags_quad (istream &) {
  Assert (false, ExcFunctionNotUseful());
};



template <int dim>
void Triangulation<dim>::load_user_flags_quad (istream &in) {
  unsigned int magic_number;
  in >> magic_number;
  Assert (magic_number==mn_tria_quad_user_flags_begin, ExcGridReadError());

  unsigned int N;
  in >> N;
  Assert (N==n_quads(), ExcGridsDoNotMatch(N, n_quads()));

  unsigned char *flags = new unsigned char[N/8+1];
  unsigned short int tmp;
  for (unsigned int i=0; i<N/8+1; ++i) 
    {
      in >> tmp;
      flags[i] = tmp;
    };
  
  
  quad_iterator quad = begin_quad(),
		endc = end_quad();
  unsigned int position=0;
  for (; quad!=endc; ++quad, ++position)
    if (flags[position/8] & (1<<(position%8)))
      quad->set_user_flag();
    else
      quad->clear_user_flag();

  Assert (position==N, ExcGridReadError());

  in >> magic_number;
  Assert (magic_number==mn_tria_quad_user_flags_end, ExcGridReadError());

  delete[] flags;
};



TriaDimensionInfo<1>::raw_cell_iterator
Triangulation<1>::begin_raw (const unsigned int level) const {
  return begin_raw_line (level);
};



TriaDimensionInfo<1>::cell_iterator
Triangulation<1>::begin (const unsigned int level) const {
  return begin_line (level);
};



TriaDimensionInfo<1>::active_cell_iterator
Triangulation<1>::begin_active (const unsigned int level) const {
  return begin_active_line (level);
};



TriaDimensionInfo<1>::raw_cell_iterator
Triangulation<1>::end () const {
  return end_line ();
};



TriaDimensionInfo<1>::raw_cell_iterator
Triangulation<1>::last_raw () const {
  return last_raw_line ();
};



TriaDimensionInfo<1>::raw_cell_iterator
Triangulation<1>::last_raw (const unsigned int level) const {
  return last_raw_line (level);
};



TriaDimensionInfo<1>::cell_iterator
Triangulation<1>::last () const {
  return last_line ();
};



TriaDimensionInfo<1>::cell_iterator
Triangulation<1>::last (const unsigned int level) const {
  return last_line (level);
};



TriaDimensionInfo<1>::active_cell_iterator
Triangulation<1>::last_active () const {
  return last_active_line ();
};



TriaDimensionInfo<1>::active_cell_iterator
Triangulation<1>::last_active (const unsigned int level) const {
  return last_active_line (level);
};





TriaDimensionInfo<2>::raw_cell_iterator
Triangulation<2>::begin_raw (const unsigned int level) const {
  return begin_raw_quad (level);
};



TriaDimensionInfo<2>::cell_iterator
Triangulation<2>::begin (const unsigned int level) const {
  return begin_quad (level);
};



TriaDimensionInfo<2>::active_cell_iterator
Triangulation<2>::begin_active (const unsigned int level) const {
  return begin_active_quad (level);
};



TriaDimensionInfo<2>::raw_cell_iterator
Triangulation<2>::end () const {
  return end_quad ();
};



TriaDimensionInfo<2>::raw_cell_iterator
Triangulation<2>::last_raw () const {
  return last_raw_quad ();
};



TriaDimensionInfo<2>::raw_cell_iterator
Triangulation<2>::last_raw (const unsigned int level) const {
  return last_raw_quad (level);
};



TriaDimensionInfo<2>::cell_iterator
Triangulation<2>::last () const {
  return last_quad ();
};



TriaDimensionInfo<2>::cell_iterator
Triangulation<2>::last (const unsigned int level) const {
  return last_quad (level);
};



TriaDimensionInfo<2>::active_cell_iterator
Triangulation<2>::last_active () const {
  return last_active_quad ();
};



TriaDimensionInfo<2>::active_cell_iterator
Triangulation<2>::last_active (const unsigned int level) const {
  return last_active_quad (level);
};







TriaDimensionInfo<1>::raw_face_iterator
Triangulation<1>::begin_raw_face (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



TriaDimensionInfo<1>::face_iterator
Triangulation<1>::begin_face (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



TriaDimensionInfo<1>::active_face_iterator
Triangulation<1>::begin_active_face (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



TriaDimensionInfo<1>::raw_face_iterator
Triangulation<1>::end_face () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



TriaDimensionInfo<1>::raw_face_iterator
Triangulation<1>::last_raw_face () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



TriaDimensionInfo<1>::raw_face_iterator
Triangulation<1>::last_raw_face (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



TriaDimensionInfo<1>::face_iterator
Triangulation<1>::last_face () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



TriaDimensionInfo<1>::face_iterator
Triangulation<1>::last_face (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



TriaDimensionInfo<1>::active_face_iterator
Triangulation<1>::last_active_face () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



TriaDimensionInfo<1>::active_face_iterator
Triangulation<1>::last_active_face (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};




TriaDimensionInfo<2>::raw_face_iterator
Triangulation<2>::begin_raw_face (const unsigned int level) const {
  return begin_raw_line (level);
};



TriaDimensionInfo<2>::face_iterator
Triangulation<2>::begin_face (const unsigned int level) const {
  return begin_line (level);
};



TriaDimensionInfo<2>::active_face_iterator
Triangulation<2>::begin_active_face (const unsigned int level) const {
  return begin_active_line (level);
};



TriaDimensionInfo<2>::raw_face_iterator
Triangulation<2>::end_face () const {
  return end_line ();
};



TriaDimensionInfo<2>::raw_face_iterator
Triangulation<2>::last_raw_face () const {
  return last_raw_line ();
};



TriaDimensionInfo<2>::raw_face_iterator
Triangulation<2>::last_raw_face (const unsigned int level) const {
  return last_raw_line (level);
};



TriaDimensionInfo<2>::face_iterator
Triangulation<2>::last_face () const {
  return last_line ();
};



TriaDimensionInfo<2>::face_iterator
Triangulation<2>::last_face (const unsigned int level) const {
  return last_line (level);
};



TriaDimensionInfo<2>::active_face_iterator
Triangulation<2>::last_active_face () const {
  return last_active_line ();
};



TriaDimensionInfo<2>::active_face_iterator
Triangulation<2>::last_active_face (const unsigned int level) const {
  return last_active_line (level);
};





template <int dim>
typename TriaDimensionInfo<dim>::raw_line_iterator
Triangulation<dim>::begin_raw_line (unsigned int level) const {
  Assert (level<levels.size(),
	  ExcInvalidLevel(level));
  if (levels[level]->lines.lines.size() == 0)
    return end_line ();
  
  return raw_line_iterator (const_cast<Triangulation<dim>*>(this),
			    level,
			    0);
};




TriaDimensionInfo<1>::raw_quad_iterator
Triangulation<1>::begin_raw_quad (unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <int dim>
typename TriaDimensionInfo<dim>::raw_quad_iterator
Triangulation<dim>::begin_raw_quad (unsigned int level) const {
  Assert (level<levels.size(),
	  ExcInvalidLevel(level));
  
  if (levels[level]->quads.quads.size() == 0)
    return end_quad();
  
  return raw_quad_iterator (const_cast<Triangulation<dim>*>(this),
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




TriaDimensionInfo<1>::quad_iterator
Triangulation<1>::begin_quad (unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
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



TriaDimensionInfo<1>::active_quad_iterator
Triangulation<1>::begin_active_quad (unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
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
typename TriaDimensionInfo<dim>::raw_line_iterator
Triangulation<dim>::end_line () const {
  return raw_line_iterator (const_cast<Triangulation<dim>*>(this),
			    -1,
			    -1);
};




TriaDimensionInfo<1>::raw_quad_iterator
Triangulation<1>::end_quad () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <int dim>
typename TriaDimensionInfo<dim>::raw_quad_iterator
Triangulation<dim>::end_quad () const {
  return raw_quad_iterator (const_cast<Triangulation<dim>*>(this),
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





TriaDimensionInfo<1>::raw_quad_iterator
Triangulation<1>::last_raw_quad (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
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
typename TriaDimensionInfo<dim>::raw_line_iterator
Triangulation<dim>::last_raw_line () const {
  return last_raw_line (levels.size()-1);
};




TriaDimensionInfo<1>::raw_quad_iterator
Triangulation<1>::last_raw_quad () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <int dim>
typename TriaDimensionInfo<dim>::raw_quad_iterator
Triangulation<dim>::last_raw_quad () const {
  return last_raw_quad (levels.size()-1);
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




TriaDimensionInfo<1>::quad_iterator
Triangulation<1>::last_quad (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
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
typename TriaDimensionInfo<dim>::line_iterator
Triangulation<dim>::last_line () const {
  return last_line (levels.size()-1);
};



TriaDimensionInfo<1>::quad_iterator
Triangulation<1>::last_quad () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};




template <int dim>
typename TriaDimensionInfo<dim>::quad_iterator
Triangulation<dim>::last_quad () const {
  return last_quad (levels.size()-1);
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



TriaDimensionInfo<1>::active_quad_iterator
Triangulation<1>::last_active_quad (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
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
typename TriaDimensionInfo<dim>::active_line_iterator
Triangulation<dim>::last_active_line () const {
  return last_active_line (levels.size()-1);
};



TriaDimensionInfo<1>::active_quad_iterator
Triangulation<1>::last_active_quad () const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};




template <int dim>
typename TriaDimensionInfo<dim>::active_quad_iterator
Triangulation<dim>::last_active_quad () const {
  return last_active_quad (levels.size()-1);
};







template <int dim>
unsigned int Triangulation<dim>::n_cells () const {
  unsigned int n=0;
  for (unsigned int l=0; l<levels.size(); ++l)
    n += n_cells (l);
  return n;
};



unsigned int Triangulation<1>::n_cells (const unsigned int level) const {
  return n_lines (level);
};



unsigned int Triangulation<2>::n_cells (const unsigned int level) const {
  return n_quads (level);
};



template <int dim>
unsigned int Triangulation<dim>::n_active_cells () const {
  unsigned int n=0;
  for (unsigned int l=0; l<levels.size(); ++l)
    n += n_active_cells (l);
  return n;
};



unsigned int Triangulation<1>::n_active_cells (const unsigned int level) const {
  return n_active_lines (level);
};



unsigned int Triangulation<2>::n_active_cells (const unsigned int level) const {
  return n_active_quads (level);
};



template <int dim>
unsigned int Triangulation<dim>::n_lines () const {
  unsigned int n=0;
  for (unsigned int l=0; l<levels.size(); ++l)
    n += n_lines (l);
  return n;
};



template <int dim>
unsigned int Triangulation<dim>::n_lines (const unsigned int level) const {
  if (levels[level]->lines.lines.size() == 0)
    return 0;

				   // only evaluate begin_/end_line if there
				   // are lines.
  line_iterator line = begin_line (level),
		endc = (level == levels.size()-1 ?
			line_iterator(end_line()) :
			begin_line (level+1));
  unsigned int n=0;
  for (; line!=endc; ++line)
    ++n;
  return n;
};



template <int dim>
unsigned int Triangulation<dim>::n_active_lines () const {
  unsigned int n=0;
  for (unsigned int l=0; l<levels.size(); ++l)
    n += n_active_lines (l);
  return n;
};



template <int dim>
unsigned int Triangulation<dim>::n_active_lines (const unsigned int level) const {
  if (levels[level]->lines.lines.size() == 0)
    return 0;

				   // only evaluate begin_/end_line if there
				   // are lines.
  active_line_iterator line = begin_active_line (level),
		       endc = (level == levels.size()-1 ?
			       active_line_iterator(end_line()) :
			       begin_active_line (level+1));
  unsigned int n=0;
  for (; line!=endc; ++line)
    ++n;
  return n;
};



template <int dim>
unsigned int Triangulation<dim>::n_quads () const {
  unsigned int n=0;
  for (unsigned int l=0; l<levels.size(); ++l)
    n += n_quads (l);
  return n;
};



unsigned int Triangulation<1>::n_quads (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <int dim>
unsigned int Triangulation<dim>::n_quads (const unsigned int level) const {
  if (levels[level]->quads.quads.size() == 0)
    return 0;

				   // only evaluate begin_/end_quad if there
				   // are quads.
  quad_iterator quad = begin_quad (level),
		endc = (level == levels.size()-1 ?
			quad_iterator(end_quad()) :
			begin (level+1));
  unsigned int n=0;
  for (; quad!=endc; ++quad)
    ++n;
  return n;
};



template <int dim>
unsigned int Triangulation<dim>::n_active_quads () const {
  unsigned int n=0;
  for (unsigned int l=0; l<levels.size(); ++l)
    n += n_active_quads (l);
  return n;
};



unsigned int Triangulation<1>::n_active_quads (const unsigned int) const {
  Assert (false, ExcFunctionNotUseful());
  return 0;
};



template <int dim>
unsigned int Triangulation<dim>::n_active_quads (const unsigned int level) const {
  if (levels[level]->quads.quads.size() == 0)
    return 0;

				   // only evaluate begin_/end_quad if there
				   // are quads.
  active_quad_iterator quad = begin_active_quad (level),
		       endc = (level == levels.size()-1 ?
			       active_quad_iterator(end_quad()) :
			       begin_active_quad (level+1));
  unsigned int n=0;
  for (; quad!=endc; ++quad)
    ++n;
  return n;
};




template <int dim>
unsigned int Triangulation<dim>::n_levels () const {
  return levels.size();
};



unsigned int Triangulation<1>::max_adjacent_cells () const {
  return 2;
};



template <int dim>
unsigned int Triangulation<dim>::max_adjacent_cells () const {
  cell_iterator cell = begin(0),
		endc = (levels.size() > 1 ? begin(1) : cell_iterator(end()));
				   // store the largest index of the vertices
				   // used on level 0
  unsigned int max_vertex_index = 0;
  for (; cell!=endc; ++cell)
    for (unsigned vertex=0; vertex<(1<<dim); ++vertex)
      if (cell->vertex_index(vertex) > (signed int)max_vertex_index)
	max_vertex_index = cell->vertex_index(vertex);

				   // store the number of times a cell touches
				   // a vertex. An unsigned int should suffice,
				   // even for larger dimensions
  vector<unsigned short int> usage_count (max_vertex_index+1, 0);
				   // touch a vertex's usage count everytime
				   // we find an adjacent element
  for (cell=begin(); cell!=endc; ++cell)
    for (unsigned vertex=0; vertex<(1<<dim); ++vertex)
      ++usage_count[cell->vertex_index(vertex)];

  return max (static_cast<unsigned int>(1<<dim),
	      static_cast<unsigned int>(*(max_element (usage_count.begin(),
						       usage_count.end()))));
};





template <int dim>
void Triangulation<dim>::print_gnuplot (ostream &out) const {
  for (unsigned int i=0; i<levels.size(); ++i)
    print_gnuplot (out, i);
};



template <int dim>
void Triangulation<dim>::print_gnuplot (ostream &out, const unsigned int level) const {
  active_cell_iterator cell = begin_active (level),
		       endc = (level == levels.size()-1 ?
			       active_cell_iterator(end()) :
			       begin_active (level+1));

  out << "#Active cells on level " << level
      << ": " << n_active_cells (level) << endl;
  for (; cell != endc; ++cell)
    print_gnuplot (out, cell);
};



void Triangulation<1>::print_gnuplot (ostream &out,
				      const TriaDimensionInfo<1>::active_cell_iterator & cell) const {
  out << cell->vertex(0) << " " << cell->level() << endl
      << cell->vertex(1) << " " << cell->level() << endl
      << endl;
};



void Triangulation<2>::print_gnuplot (ostream &out,
				      const TriaDimensionInfo<2>::active_cell_iterator & cell) const {
  out << cell->vertex(0) << " " << cell->level() << endl
      << cell->vertex(1) << " " << cell->level() << endl
      << cell->vertex(2) << " " << cell->level() << endl
      << cell->vertex(3) << " " << cell->level() << endl
      << cell->vertex(0) << " " << cell->level() << endl
      << endl  // double new line for gnuplot 3d plots
      << endl;
};




template <int dim>
void Triangulation<dim>::refine (const dVector &criteria,
				 const double   threshold) {
  Assert (criteria.n() == n_active_cells(),
	  ExcInvalidVectorSize(criteria.n(), n_active_cells()));

  active_cell_iterator cell = begin_active();
  const unsigned int n_cells = criteria.n();
  
  for (unsigned int index=0; index<n_cells; ++cell, ++index)
    if (criteria(index) > threshold)
      cell->set_refine_flag();
};



template <int dim>
void Triangulation<dim>::refine_fixed_number (const dVector &criteria,
					      const double   fraction) {
  dVector tmp(criteria);
  nth_element (tmp.begin(),
	       tmp.begin()+static_cast<int>(fraction*tmp.n()),
	       tmp.end(),
	       greater<double>());

  refine (criteria, *(tmp.begin() +
		      static_cast<int>(fraction*tmp.n())));
};




void Triangulation<1>::execute_refinement () {
  prepare_refinement ();

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
      unsigned int used_cells = 0;
      count_if (levels[level+1]->lines.used.begin(),
		levels[level+1]->lines.used.end(),
		bind2nd (equal_to<bool>(), true),
		used_cells);

				       // reserve space for the used_cells
				       // cells already existing on the next
				       // higher level as well as for the
				       // 2*flagged_cells that will be created
				       // on that level
      levels[level+1]->
	TriangulationLevel<0>::reserve_space (used_cells+2*flagged_cells, 1);
				       // reserve space for 2*flagged_cells
				       // new lines on the next higher
				       // level
      levels[level+1]->
	TriangulationLevel<1>::reserve_space (2*flagged_cells);
      
      needed_vertices += flagged_cells;
    };

				   // add to needed vertices how
				   // many vertices are already in use
  count_if (vertices_used.begin(), vertices_used.end(),
	    bind2nd (equal_to<bool>(), true),
	    needed_vertices);

  vertices.resize (needed_vertices, Point<1>());
  vertices_used.resize (needed_vertices, false);


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
	    ++next_unused_cell;
	    Assert (next_unused_cell->used() == false,
		    ExcCellShouldBeUnused());
	    second_child = next_unused_cell;
	    second_child->set_used_flag ();

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
						 // since the neighbors level is
						 // always <=level, if the
						 // cell is active, then there
						 // are no cells to the left which
						 // may want to know about this
						 // new child cell.
		first_child->set_neighbor (0, cell->neighbor(0));
	      else
						 // left neighbor is refined
		{
						   // set neighbor to cell on
						   // same level
		  first_child->set_neighbor (0, cell->neighbor(0)->child(1));

						   // reset neighbor info of
						   // leftmost descendant of the
						   // left neighbor of cell
		  cell_iterator left_neighbor = cell->neighbor(0);
		  while (left_neighbor->active() == false)
		    left_neighbor = left_neighbor->child(1);
		  left_neighbor->set_neighbor(1, first_child);
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
		second_child->set_neighbor (1, cell->neighbor(1));
	      else
						 // right neighbor is refined
		{
		  second_child->set_neighbor (1, cell->neighbor(1)->child(0));
		  
		  cell_iterator right_neighbor = cell->neighbor(1);
		  while (right_neighbor->active() == false)
		    right_neighbor = right_neighbor->child(0);
		  right_neighbor->set_neighbor(0, second_child);
		};
	  };      
    };

#ifdef DEBUG
  for (unsigned int level=0; level<levels.size(); ++level) 
    levels[level]->monitor_memory (1);
#endif
};






void Triangulation<2>::execute_refinement () {
  prepare_refinement ();

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
	    for (unsigned int nb=0; nb<4; ++nb) 
	      {
		cell_iterator neighbor = acell->neighbor(nb);
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
      
				       // count number of used lines on
				       // the next higher level
      unsigned int used_lines = 0;
      count_if (levels[level+1]->lines.used.begin(),
		levels[level+1]->lines.used.end(),
		bind2nd (equal_to<bool>(), true),
		used_lines);

      				       // count number of used cells on
				       // the next higher level
      unsigned int used_cells = 0;
      count_if (levels[level+1]->quads.used.begin(),
		levels[level+1]->quads.used.end(),
		bind2nd (equal_to<bool>(), true),
		used_cells);
      
      
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
				       // new lines on the next higher
				       // level
      levels[level+1]->
	TriangulationLevel<2>::reserve_space (4*flagged_cells);
    };

				   // add to needed vertices how
				   // many vertices are already in use
  count_if (vertices_used.begin(), vertices_used.end(),
	    bind2nd (equal_to<bool>(), true),
	    needed_vertices);

				   // reserve enough space for all vertices
  vertices.resize (needed_vertices, Point<2>());
  vertices_used.resize (needed_vertices, false);



				   // Do REFINEMENT
				   // on every level; exclude highest level as
				   // above

				   //  index of next unused vertex
  unsigned int next_unused_vertex = 0;
  
  for (int level=0; level<(int)levels.size()-1; ++level) 
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
		    if (cell->at_boundary(nb)) 
		      {
							 // boundary vertex
			const Point<2> *neighbors[2] =
			{&vertices[new_vertices[2*nb]],
			 &vertices[new_vertices[(2*nb+2)%8]]};
			new_point = boundary->in_between (neighbors);
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
			
			new_lines[nb*2+1]->set(Line(new_vertices[2*nb+1],
						    new_vertices[(2*nb+2)%8]));
			new_lines[nb*2+1]->set_used_flag ();
			new_lines[nb*2+1]->clear_children ();
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

			new_lines[nb*2+1]->set(Line(new_vertices[(2*nb+2)%8],
						    new_vertices[2*nb+1]));
			new_lines[nb*2+1]->set_used_flag ();
			new_lines[nb*2+1]->clear_children ();			
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
	    
	    while (next_unused_line->used() == true)
	      ++next_unused_line;
	    new_lines[9] = next_unused_line;
	    new_lines[9]->set(Line(new_vertices[8],
				   new_vertices[5]));
	    new_lines[9]->set_used_flag ();
	    new_lines[9]->clear_children ();

	    while (next_unused_line->used() == true)
	      ++next_unused_line;
	    new_lines[10] = next_unused_line;
	    new_lines[10]->set(Line(new_vertices[7],
				    new_vertices[8]));
	    new_lines[10]->set_used_flag ();
	    new_lines[10]->clear_children ();

	    while (next_unused_line->used() == true)
	      ++next_unused_line;
	    new_lines[11] = next_unused_line;
	    new_lines[11]->set(Line(new_vertices[8],
				    new_vertices[3]));
	    new_lines[11]->set_used_flag ();
	    new_lines[11]->clear_children ();

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


	    subcells[1]->set (Quad(new_lines[1]->index(),  new_lines[2]->index(),
				   new_lines[11]->index(), new_lines[8]->index()));
	    subcells[1]->set_used_flag();
	    subcells[1]->clear_children();


	    subcells[2]->set (Quad(new_lines[11]->index(),  new_lines[3]->index(),
				   new_lines[4]->index(), new_lines[9]->index()));
	    subcells[2]->set_used_flag();
	    subcells[2]->clear_children();


	    subcells[3]->set (Quad(new_lines[10]->index(),  new_lines[9]->index(),
				   new_lines[5]->index(), new_lines[6]->index()));
	    subcells[3]->set_used_flag();
	    subcells[3]->clear_children();

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
#ifdef DEBUG
  for (unsigned int level=0; level<levels.size(); ++level) 
    levels[level]->monitor_memory (2);
#endif
};

  


template <int dim>
void Triangulation<dim>::prepare_refinement () {
				   // make sure no two adjacent active cells
				   // have refinement levels differing
				   // with more than one.
				   // Precondition: on the old grid,
				   // there are no such cells.
  if (dim>=2) 
    {
      active_cell_iterator cell = last_active(),
			   endc = end();
				       // loop over active cells
      for (; cell != endc; --cell)
	if (cell->refine_flag_set() == true)
					   // loop over neighbors of cell
	  for (unsigned int i=0; i<(2*dim); ++i)
	    if (cell->neighbor(i).state() == valid)
	      if ((cell->neighbor_level(i) == cell->level()-1)
		  &&
		  (cell->neighbor(i)->refine_flag_set() == false))
		cell->neighbor(i)->set_refine_flag();
    };


				   // check whether a new level is needed
  raw_cell_iterator cell = begin_active (levels.size()-1),
		    endc = end();
  for (; cell != endc; ++cell)
    if (cell->refine_flag_set()==true) 
      {
	levels.push_back (new TriangulationLevel<dim>);
	return;
      };
};





void TriangulationLevel<0>::reserve_space (const unsigned int total_cells,
					   const unsigned int dimension) {
  refine_flags.reserve (total_cells);
  refine_flags.insert (refine_flags.end(),
		       total_cells - refine_flags.size(),
		       false);

  neighbors.reserve (total_cells*(2*dimension));
  neighbors.insert (neighbors.end(),
		    total_cells*(2*dimension) - neighbors.size(),
		    make_pair(-1,-1));
};




void TriangulationLevel<0>::monitor_memory (const unsigned int true_dimension) const {
//  Assert (refine_flags.size() == refine_flags.capacity() ||
//	  refine_flags.size()<256,
//	  ExcMemoryWasted ("refine_flags",
//			   refine_flags.size(), refine_flags.capacity()));
//  Assert (neighbors.size() ==  neighbors.capacity() ||
//	  neighbors.size()<256,
//	  ExcMemoryWasted ("neighbors",
//			   neighbors.size(), neighbors.capacity()));
  Assert (2*true_dimension*refine_flags.size() == neighbors.size(),
	  ExcMemoryInexact (refine_flags.size(), neighbors.size()));
};



void TriangulationLevel<1>::reserve_space (const unsigned int new_lines) {
  vector<bool>::iterator u = lines.used.begin(),
			 e = lines.used.end();

  unsigned int used_lines = 0;
  for (; u!=e; ++u)
    ++used_lines;

  unsigned int new_size = used_lines + new_lines;

//  cout << "  lines: pre: siz=" << lines.lines.size() << ", cap=" << lines.lines.capacity();
  lines.lines.reserve (new_size);
//  cout << " inter: siz=" << lines.lines.size() << ", cap=" << lines.lines.capacity()
//       << " (newsize=" << new_size << ")";
  lines.lines.insert (lines.lines.end(), new_size-lines.lines.size(), Line());
//  cout << " post: siz=" << lines.lines.size() << ", cap=" << lines.lines.capacity() << endl;
  
//  cout << "  used : pre: siz=" << lines.used.size() << ", cap=" << lines.used.capacity();
  lines.used.reserve (new_size);
//  cout << " inter: siz=" << lines.used.size() << ", cap=" << lines.used.capacity()
//       << " (newsize=" << new_size << ")";
  lines.used.insert (lines.used.end(), new_size-lines.used.size(), false);
//  cout << " post: siz=" << lines.used.size() << ", cap=" << lines.used.capacity() << endl;
  
  lines.user_flags.reserve (new_size);
  lines.user_flags.insert (lines.user_flags.end(),
			   new_size-lines.user_flags.size(), false);
  
  lines.children.reserve (new_size);
  lines.children.insert (lines.children.end(), new_size-lines.children.size(),
			 -1);

  lines.material_id.reserve (new_size);
  lines.material_id.insert (lines.material_id.end(),
			    new_size-lines.material_id.size(),
			    255);

};




void TriangulationLevel<1>::monitor_memory (const unsigned int true_dimension) const {
//  Assert (lines.lines.size() == lines.lines.capacity() ||
//	  lines.lines.size()<256,
//	  ExcMemoryWasted ("lines",
//			   lines.lines.size(), lines.lines.capacity()));
//  Assert (lines.children.size() == lines.children.capacity() ||
//	  lines.children.size()<256,
//	  ExcMemoryWasted ("children",
//			   lines.children.size(), lines.children.capacity()));
//  Assert (lines.used.size() == lines.used.capacity() ||
//	  lines.used.size()<256,
//	  ExcMemoryWasted ("used",
//			   lines.used.size(), lines.used.capacity()));
//  Assert (lines.user_flags.size() == lines.user_flags.capacity() ||
//	  lines.user_flags.size()<256,
//	  ExcMemoryWasted ("user_flags",
//			   lines.user_flags.size(), lines.user_flags.capacity()));
  Assert (lines.lines.size() == lines.used.size(),
	  ExcMemoryInexact (lines.lines.size(), lines.used.size()));
  Assert (lines.lines.size() == lines.user_flags.size(),
	  ExcMemoryInexact (lines.lines.size(), lines.user_flags.size()));
  Assert (lines.lines.size() == lines.children.size(),
	  ExcMemoryInexact (lines.lines.size(), lines.children.size()));
  Assert (lines.lines.size() == lines.material_id.size(),
	  ExcMemoryInexact (lines.lines.size(), lines.material_id.size()));
  Assert (lines.used[lines.used.size()-1]==true ,
	  ExcUnusedMemoryAtEnd());

  TriangulationLevel<0>::monitor_memory (true_dimension);
};



void TriangulationLevel<2>::reserve_space (const unsigned int new_quads) {
  vector<bool>::iterator u = quads.used.begin(),
			 e = quads.used.end();

  unsigned int used_quads = 0;
  for (; u!=e; ++u)
    ++used_quads;

  unsigned int new_size = used_quads + new_quads;
  
  quads.quads.reserve (new_size);
  quads.quads.insert (quads.quads.end(), new_size-quads.quads.size(), Quad());
  
  quads.used.reserve (new_size);
  quads.used.insert (quads.used.end(), new_size-quads.used.size(), false);
  
  quads.user_flags.reserve (new_size);
  quads.user_flags.insert (quads.user_flags.end(),
			   new_size-quads.user_flags.size(), false);
  
  quads.children.reserve (new_size);
  quads.children.insert (quads.children.end(), new_size-quads.children.size(),
			 -1);

  quads.material_id.reserve (new_size);
  quads.material_id.insert (quads.material_id.end(),
			    new_size-quads.material_id.size(),
			    255);
};



void TriangulationLevel<2>::monitor_memory (const unsigned int true_dimension) const {
//  Assert (quads.quads.size() == quads.quads.capacity() ||
//	  quads.quads.size()<256,
//	  ExcMemoryWasted ("quads",
//			   quads.quads.size(), quads.quads.capacity()));
//  Assert (quads.children.size() == quads.children.capacity() ||
//	  quads.children.size()<256,
//	  ExcMemoryWasted ("children",
//			   quads.children.size(), quads.children.capacity()));
//  Assert (quads.used.size() == quads.used.capacity() ||
//	  quads.used.size()<256,
//	  ExcMemoryWasted ("used",
//			   quads.used.size(), quads.used.capacity()));
//  Assert (quads.user_flags.size() == quads.user_flags.capacity() ||
//	  quads.user_flags.size()<256,
//	  ExcMemoryWasted ("user_flags",
//			   quads.user_flags.size(), quads.user_flags.capacity()));
  Assert (quads.quads.size() == quads.used.size(),
	  ExcMemoryInexact (quads.quads.size(), quads.used.size()));
  Assert (quads.quads.size() == quads.user_flags.size(),
	  ExcMemoryInexact (quads.quads.size(), quads.user_flags.size()));
  Assert (quads.quads.size() == quads.children.size(),
	  ExcMemoryInexact (quads.quads.size(), quads.children.size()));
  Assert (quads.quads.size() == quads.material_id.size(),
	  ExcMemoryInexact (quads.quads.size(), quads.material_id.size()));
  Assert (quads.used[quads.used.size()-1]==true ,
	  ExcUnusedMemoryAtEnd());

  TriangulationLevel<1>::monitor_memory (true_dimension);
};






// explicit instantiations
template class Triangulation<1>;
template class Triangulation<2>;

