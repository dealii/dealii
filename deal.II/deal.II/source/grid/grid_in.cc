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

#include <map>


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
		Assert (false, ExcInvalidVertexIndex(cell, cells.back().vertices[i]));
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
		Assert (false,
			ExcInvalidVertexIndex(cell,
					      subcelldata.boundary_lines.back().vertices[i]));
		subcelldata.boundary_lines.back().vertices[i] = -1;
	      };
	  }
	else
					   // cannot read this
	  Assert (false, ExcUnknownIdentifier(cell_type));
    };

				   // check that no forbidden arrays are used
  Assert (subcelldata.check_consistency(dim), ExcInternalError());

  AssertThrow (in, ExcIO());

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
  int dummy;
  
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
	  --cells.back().vertices[i];
	};

				       // read and discard Ref phi_i
      in >> dummy;
    };

  skip_empty_lines(in);
  

				   // then there are again a whole lot
				   // of fields of which I have no
				   // clue what they mean. skip them
				   // all and leave the interpretation
				   // to other implementors...
  while (getline(in,line), line.find("End")==string::npos);
				   // ok, so we are not at the end of
				   // the file, that's it, mostly

  
				   // check that no forbidden arrays are used
  Assert (subcelldata.check_consistency(dim), ExcInternalError());

  AssertThrow (in, ExcIO());

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




//explicit instantiations
template class GridIn<deal_II_dimension>;
