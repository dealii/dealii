//----------------------------  grid_in.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  grid_in.cc  ---------------------------


#include <grid/grid_in.h>
#include <grid/tria.h>
#include <grid/grid_reordering.h>

#include <map>
#include <algorithm>





#include <fstream>


template <int dim>
GridIn<dim>::GridIn () :
		tria(0) {};


template <int dim>
void GridIn<dim>::attach_triangulation (Triangulation<dim> &t)
{
  tria = &t;
};


template <int dim>
void GridIn<dim>::read_ucd (istream &in)
{
  Assert (tria != 0, ExcNoTriangulationSelected());
  Assert ((1<=dim) && (dim<=2), ExcNotImplemented());

  AssertThrow (in, ExcIO());
  
				   // skip comments at start of file
  skip_comment_lines (in, '#');


  unsigned int n_vertices;
  unsigned int n_cells;
  int dummy;

  in >> n_vertices
     >> n_cells
     >> dummy         // number of data vectors
     >> dummy         // cell data
     >> dummy;        // model data

				   // set up array of vertices
  vector<Point<dim> >     vertices (n_vertices);
				   // set up mapping between numbering
				   // in ucd-file (key) and in the
				   // vertices vector
  map<int,int> vertex_indices;
  
  for (unsigned int vertex=0; vertex<n_vertices; ++vertex) 
    {
      int vertex_number;
      double x[3];

				       // read vertex
      in >> vertex_number
	 >> x[0] >> x[1] >> x[2];

				       // store vertex
      for (unsigned int d=0; d<dim; ++d)
	vertices[vertex](d) = x[d];
				       // store mapping; note that
				       // vertices_indices[i] is automatically
				       // created upon first usage
      vertex_indices[vertex_number] = vertex;
    };

				   // set up array of cells
  vector<CellData<dim> > cells;
  SubCellData            subcelldata;

  for (unsigned int cell=0; cell<n_cells; ++cell) 
    {
      string cell_type;
      int material_id;
      
      in >> dummy          // cell number
	 >> material_id;
      in >> cell_type;

      if (((cell_type == "line") && (dim == 1)) ||
	  ((cell_type == "quad") && (dim == 2)))
					 // found a cell
	{
					   // allocate and read indices
	  cells.push_back (CellData<dim>());
	  for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
	    in >> cells.back().vertices[i];
	  cells.back().material_id = material_id;

					   // transform from ucd to
					   // consecutive numbering
	  for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
	    if (vertex_indices.find (cells.back().vertices[i]) != vertex_indices.end())
					       // vertex with this index exists
	      cells.back().vertices[i] = vertex_indices[cells.back().vertices[i]];
	    else 
	      {
						 // no such vertex index
		AssertThrow (false, ExcInvalidVertexIndex(cell, cells.back().vertices[i]));
		cells.back().vertices[i] = -1;
	      };
	}
      else
	if ((cell_type == "line") && (dim == 2))
					   // boundary info
	  {
	    subcelldata.boundary_lines.push_back (CellData<1>());
	    in >> subcelldata.boundary_lines.back().vertices[0]
	       >> subcelldata.boundary_lines.back().vertices[1];
	    subcelldata.boundary_lines.back().material_id = material_id;

					     // transform from ucd to
					     // consecutive numbering
	    for (unsigned int i=0; i<2; ++i)
	      if (vertex_indices.find (subcelldata.boundary_lines.back().vertices[i]) !=
		  vertex_indices.end())
						 // vertex with this index exists
		subcelldata.boundary_lines.back().vertices[i]
		  = vertex_indices[subcelldata.boundary_lines.back().vertices[i]];
	    else 
	      {
						 // no such vertex index
		AssertThrow (false,
			     ExcInvalidVertexIndex(cell,
						   subcelldata.boundary_lines.back().vertices[i]));
		subcelldata.boundary_lines.back().vertices[i] = -1;
	      };
	  }
	else
					   // cannot read this
	  AssertThrow (false, ExcUnknownIdentifier(cell_type));
    };

				   // check that no forbidden arrays are used
  Assert (subcelldata.check_consistency(dim), ExcInternalError());

  AssertThrow (in, ExcIO());

				   // do some clean-up on vertices
  delete_unused_vertices (vertices, cells, subcelldata);
  tria->create_triangulation (vertices, cells, subcelldata);
};



template <int dim>
void GridIn<dim>::read_dbmesh (istream &in)
{
  Assert (tria != 0, ExcNoTriangulationSelected());
  Assert (dim==2, ExcNotImplemented());

  AssertThrow (in, ExcIO());
  
				   // skip comments at start of file
  skip_comment_lines (in, '#');


				   // first read in identifier string
  string line;
  getline (in, line);

  AssertThrow (line=="MeshVersionFormatted 0",
	       ExcInvalidDBMESHInput(line));

  skip_empty_lines (in);

				   // next read dimension
  getline (in, line);
  AssertThrow (line=="Dimension", ExcInvalidDBMESHInput(line));
  unsigned int dimension;
  in >> dimension;
  AssertThrow (dimension == dim, ExcDBMESHWrongDimension(dimension));
  skip_empty_lines (in);

				   // now there are a lot of fields of
				   // which we don't know the exact
				   // meaning and which are far from
				   // being properly documented in the
				   // manual. we skip everything until
				   // we find a comment line with the
				   // string "# END". at some point in
				   // the future, someone may have the
				   // knowledge to parse and interpret
				   // the other fields in between as
				   // well...
  while (getline(in,line), line.find("# END")==string::npos);
  skip_empty_lines (in);


				   // now read vertices
  getline (in, line);
  AssertThrow (line=="Vertices", ExcInvalidDBMESHInput(line));
  
  unsigned int n_vertices;
  double dummy;
  
  in >> n_vertices;
  vector<Point<dim> >     vertices (n_vertices);
  for (unsigned int vertex=0; vertex<n_vertices; ++vertex)
    {
				       // read vertex coordinates
      for (unsigned int d=0; d<dim; ++d)
	in >> vertices[vertex][d];
				       // read Ref phi_i, whatever that may be
      in >> dummy;
    };
  AssertThrow (in, ExcInvalidDBMeshFormat());

  skip_empty_lines(in);

				   // read edges. we ignore them at
				   // present, so just read them and
				   // discard the input
  getline (in, line);
  AssertThrow (line=="Edges", ExcInvalidDBMESHInput(line));
  
  unsigned int n_edges;
  in >> n_edges;
  for (unsigned int edge=0; edge<n_edges; ++edge)
    {
				       // read vertex indices
      in >> dummy >> dummy;
				       // read Ref phi_i, whatever that may be
      in >> dummy;
    };
  AssertThrow (in, ExcInvalidDBMeshFormat());

  skip_empty_lines(in);



				   // read cracked edges (whatever
				   // that may be). we ignore them at
				   // present, so just read them and
				   // discard the input
  getline (in, line);
  AssertThrow (line=="CrackedEdges", ExcInvalidDBMESHInput(line));
  
  in >> n_edges;
  for (unsigned int edge=0; edge<n_edges; ++edge)
    {
				       // read vertex indices
      in >> dummy >> dummy;
				       // read Ref phi_i, whatever that may be
      in >> dummy;
    };
  AssertThrow (in, ExcInvalidDBMeshFormat());

  skip_empty_lines(in);


				   // now read cells.
				   // set up array of cells
  getline (in, line);
  AssertThrow (line=="Quadrilaterals", ExcInvalidDBMESHInput(line));

  vector<CellData<dim> > cells;
  SubCellData            subcelldata;
  unsigned int n_cells;
  in >> n_cells;
  for (unsigned int cell=0; cell<n_cells; ++cell) 
    {
				       // read in vertex numbers. they
				       // are 1-based, so subtract one
      cells.push_back (CellData<dim>());
      for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
	{
	  in >> cells.back().vertices[i];
	  
	  AssertThrow ((cells.back().vertices[i] >= 1)
		       &&
		       (static_cast<unsigned int>(cells.back().vertices[i]) <= vertices.size()),
		       ExcInvalidVertexIndex(cell, cells.back().vertices[i]));
		  
	  --cells.back().vertices[i];
	};

				       // read and discard Ref phi_i
      in >> dummy;
    };
  AssertThrow (in, ExcInvalidDBMeshFormat());

  skip_empty_lines(in);
  

				   // then there are again a whole lot
				   // of fields of which I have no
				   // clue what they mean. skip them
				   // all and leave the interpretation
				   // to other implementors...
  while (getline(in,line), ((line.find("End")==string::npos) && (in)));
				   // ok, so we are not at the end of
				   // the file, that's it, mostly

  
				   // check that no forbidden arrays are used
  Assert (subcelldata.check_consistency(dim), ExcInternalError());

  AssertThrow (in, ExcIO());

				   // do some clean-up on vertices
  delete_unused_vertices (vertices, cells, subcelldata);

  ofstream x1("x1");
  debug_output_grid (cells, vertices, x1);
  GridReordering<dim>::reorder_cells (cells);
  ofstream x2("x2");
  debug_output_grid (cells, vertices, x2);
  tria->create_triangulation (vertices, cells, subcelldata);
};



template <int dim>
void GridIn<dim>::skip_empty_lines (istream &in)
{    
  string line;
  while (in) 
    {
				       // get line
      getline (in, line);
				       // eat all spaces from the back
      while ((line.length()>0) && (line[line.length()-1]==' '))
	line.erase (line.length()-1, 1);
				       // if still non-null, then this
				       // is a non-empty line. put
				       // back all info and leave
      if (line.length() > 0)
	{
	  in.putback ('\n');
	  for (int i=line.length()-1; i>=0; --i)
	    in.putback (line[i]);
	  return;
	};

				       // else: go on with next line
    };
};



template <int dim>
void GridIn<dim>::skip_comment_lines (istream    &in,
				      const char  comment_start)
{    
  char c;
  while (in.get(c), c==comment_start) 
    {
      char line[256];
      in.get (line, 255, '\n'); // ignore rest of line, at most 256 chars
      in.get (c);         // ignore '\n' at end of line.
    };
  
				   // put back first character of
				   // first non-comment line
  in.putback (c);
};



template <int dim>
void
GridIn<dim>::delete_unused_vertices (vector<Point<dim> >    &vertices,
				     vector<CellData<dim> > &cells,
				     SubCellData            &subcelldata)
{
				   // first check which vertices are
				   // actually used
  vector<bool> vertex_used (vertices.size(), false);
  for (unsigned int c=0; c<cells.size(); ++c)
    for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
      vertex_used[cells[c].vertices[v]] = true;

				   // then renumber the vertices that
				   // are actually used in the same
				   // order as they were beforehand
  const unsigned int invalid_vertex = static_cast<unsigned int>(-1);
  vector<unsigned int> new_vertex_numbers (vertices.size(), invalid_vertex);
  unsigned int next_free_number = 0;
  for (unsigned int i=0; i<vertices.size(); ++i)
    if (vertex_used[i] == true)
      {
	new_vertex_numbers[i] = next_free_number;
	++next_free_number;
      };

				   // next replace old vertex numbers
				   // by the new ones
  for (unsigned int c=0; c<cells.size(); ++c)
    for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
      cells[c].vertices[v] = new_vertex_numbers[cells[c].vertices[v]];

				   // same for boundary data
  for (unsigned int c=0; c<subcelldata.boundary_lines.size(); ++c)
    for (unsigned int v=0; v<GeometryInfo<1>::vertices_per_cell; ++v)
      subcelldata.boundary_lines[c].vertices[v]
	= new_vertex_numbers[subcelldata.boundary_lines[c].vertices[v]];
  for (unsigned int c=0; c<subcelldata.boundary_quads.size(); ++c)
    for (unsigned int v=0; v<GeometryInfo<2>::vertices_per_cell; ++v)
      subcelldata.boundary_quads[c].vertices[v]
	= new_vertex_numbers[subcelldata.boundary_quads[c].vertices[v]];

				   // finally copy over the vertices
				   // which we really need to a new
				   // array and replace the old one by
				   // the new one
  vector<Point<dim> > tmp;
  tmp.reserve (count(vertex_used.begin(), vertex_used.end(), true));
  for (unsigned int v=0; v<vertices.size(); ++v)
    if (vertex_used[v] == true)
      tmp.push_back (vertices[v]);
  swap (vertices, tmp);
};



template <int dim>
void GridIn<dim>::debug_output_grid (const vector<CellData<dim> > &/*cells*/,
				     const vector<Point<dim> >    &/*vertices*/,
				     ostream                      &/*out*/)
{
  Assert (false, ExcNotImplemented());
};


#if deal_II_dimension == 2

template <>
void
GridIn<2>::debug_output_grid (const vector<CellData<2> > &cells,
			      const vector<Point<2> >    &vertices,
			      ostream                    &out)
{
  double min_x = vertices[cells[0].vertices[0]](0),
	 max_x = vertices[cells[0].vertices[0]](0),
	 min_y = vertices[cells[0].vertices[0]](1),
	 max_y = vertices[cells[0].vertices[0]](1);
  
  for (unsigned int i=0; i<cells.size(); ++i)
    {
      for (unsigned int v=0; v<4; ++v)
	{
	  const Point<2> &p = vertices[cells[i].vertices[v]];
	  
	  if (p(0) < min_x)
	    min_x = p(0);
	  if (p(0) > max_x)
	    max_x = p(0);
	  if (p(1) < min_y)
	    min_y = p(1);
	  if (p(1) > max_y)
	    max_y = p(1);
	};

      out << "# cell " << i << endl;
      Point<2> center;
      for (unsigned int f=0; f<4; ++f)
	center += vertices[cells[i].vertices[f]];
      center /= 4;

      out << "set label \"" << i << "\" at "
	  << center(0) << ',' << center(1)
	  << " center"
	  << endl;

				       // first two line right direction
      for (unsigned int f=0; f<2; ++f)
	out << "set arrow from "
	    << vertices[cells[i].vertices[f]](0) << ',' 
	    << vertices[cells[i].vertices[f]](1)
	    << " to "
	    << vertices[cells[i].vertices[(f+1)%4]](0) << ',' 
	    << vertices[cells[i].vertices[(f+1)%4]](1)
	    << endl;
				       // other two lines reverse direction
      for (unsigned int f=2; f<4; ++f)
	out << "set arrow from "
	    << vertices[cells[i].vertices[(f+1)%4]](0) << ',' 
	    << vertices[cells[i].vertices[(f+1)%4]](1)
	    << " to "
	    << vertices[cells[i].vertices[f]](0) << ',' 
	    << vertices[cells[i].vertices[f]](1)
	    << endl;
      out << endl;
    };
  

  out << endl
      << "set nokey" << endl
      << "pl [" << min_x << ':' << max_x << "]["
      << min_y << ':' << max_y <<  "] "
      << min_y << endl
      << "pause -1" << endl;
};

#endif


//explicit instantiations
template class GridIn<deal_II_dimension>;
