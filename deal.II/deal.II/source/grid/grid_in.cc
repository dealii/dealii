//----------------------------  grid_in.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2005 by the deal.II authors
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
		tria(0), default_format(ucd)
{}


template <int dim>
void GridIn<dim>::attach_triangulation (Triangulation<dim> &t)
{
  tria = &t;
}


template <int dim>
void GridIn<dim>::read_ucd (std::istream &in)
{
  Assert (tria != 0, ExcNoTriangulationSelected());
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
  std::vector<Point<dim> >     vertices (n_vertices);
				   // set up mapping between numbering
				   // in ucd-file (key) and in the
				   // vertices vector
  std::map<int,int> vertex_indices;
  
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
  std::vector<CellData<dim> > cells;
  SubCellData                 subcelldata;

  for (unsigned int cell=0; cell<n_cells; ++cell) 
    {
				       // note that since in the input
				       // file we found the number of
				       // cells at the top, there
				       // should still be input here,
				       // so check this:
      AssertThrow (in, ExcIO());
      
      std::string cell_type;
      int material_id;
      
      in >> dummy          // cell number
	 >> material_id;
      in >> cell_type;

      if (((cell_type == "line") && (dim == 1)) ||
	  ((cell_type == "quad") && (dim == 2)) ||
	  ((cell_type == "hex" ) && (dim == 3)))
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
	if ((cell_type == "line") && ((dim == 2) || (dim == 3)))
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
	  if ((cell_type == "quad") && (dim == 3))
					     // boundary info
	    {
 	      subcelldata.boundary_quads.push_back (CellData<2>());
 	      in >> subcelldata.boundary_quads.back().vertices[0]
 	         >> subcelldata.boundary_quads.back().vertices[1]
 		 >> subcelldata.boundary_quads.back().vertices[2]
 		 >> subcelldata.boundary_quads.back().vertices[3];
	      
 	      subcelldata.boundary_quads.back().material_id = material_id;
	      
					       // transform from ucd to
					       // consecutive numbering
 	      for (unsigned int i=0; i<4; ++i)
 	        if (vertex_indices.find (subcelldata.boundary_quads.back().vertices[i]) !=
 		    vertex_indices.end())
 		                                   // vertex with this index exists
		  subcelldata.boundary_quads.back().vertices[i]
		    = vertex_indices[subcelldata.boundary_quads.back().vertices[i]];
 	        else
 	          {
						     // no such vertex index
 	            Assert (false,
 		            ExcInvalidVertexIndex(cell,
 			                          subcelldata.boundary_quads.back().vertices[i]));
 		    subcelldata.boundary_quads.back().vertices[i] = -1;
 		  };
	      
	    }
	  else
					     // cannot read this
	    AssertThrow (false, ExcUnknownIdentifier(cell_type));
    };

  
				   // check that no forbidden arrays are used
  Assert (subcelldata.check_consistency(dim), ExcInternalError());

  AssertThrow (in, ExcIO());

				   // do some clean-up on vertices...
  delete_unused_vertices (vertices, cells, subcelldata);
				   // ... and cells
  GridReordering<dim>::reorder_cells (cells);
  tria->create_triangulation (vertices, cells, subcelldata);
}



template <int dim>
void GridIn<dim>::read_dbmesh (std::istream &in)
{
  Assert (tria != 0, ExcNoTriangulationSelected());
  Assert (dim==2, ExcNotImplemented());

  AssertThrow (in, ExcIO());
  
				   // skip comments at start of file
  skip_comment_lines (in, '#');


				   // first read in identifier string
  std::string line;
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
  while (getline(in,line), line.find("# END")==std::string::npos);
  skip_empty_lines (in);


				   // now read vertices
  getline (in, line);
  AssertThrow (line=="Vertices", ExcInvalidDBMESHInput(line));
  
  unsigned int n_vertices;
  double dummy;
  
  in >> n_vertices;
  std::vector<Point<dim> >     vertices (n_vertices);
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

  std::vector<CellData<dim> > cells;
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
  while (getline(in,line), ((line.find("End")==std::string::npos) && (in)));
				   // ok, so we are not at the end of
				   // the file, that's it, mostly

  
				   // check that no forbidden arrays are used
  Assert (subcelldata.check_consistency(dim), ExcInternalError());

  AssertThrow (in, ExcIO());

				   // do some clean-up on vertices...
  delete_unused_vertices (vertices, cells, subcelldata);
				   // ...and cells
  GridReordering<dim>::reorder_cells (cells);
  tria->create_triangulation (vertices, cells, subcelldata);
}



template <int dim>
void GridIn<dim>::read_xda (std::istream &)
{
  Assert (false, ExcNotImplemented());
}



// 2D XDA meshes
#if deal_II_dimension == 2

template <>
void GridIn<2>::read_xda (std::istream &in)
{
  Assert (tria != 0, ExcNoTriangulationSelected());
  AssertThrow (in, ExcIO());

  std::string line;
				   // skip comments at start of file
  getline (in, line);


  unsigned int n_vertices;
  unsigned int n_cells;

				   // read cells, throw away rest of line
  in >> n_cells;
  getline (in, line);

  in >> n_vertices;
  getline (in, line);

				   // ignore following 8 lines
  for (unsigned int i=0; i<8; ++i)
    getline (in, line);

				   // set up array of cells
  std::vector<CellData<2> > cells (n_cells);
  SubCellData subcelldata;

  for (unsigned int cell=0; cell<n_cells; ++cell) 
    {
				       // note that since in the input
				       // file we found the number of
				       // cells at the top, there
				       // should still be input here,
				       // so check this:
      AssertThrow (in, ExcIO());
      Assert (GeometryInfo<2>::vertices_per_cell == 4,
	      ExcInternalError());
      
      for (unsigned int i=0; i<4; ++i)
	in >> cells[cell].vertices[i];
    };


   
				   // set up array of vertices
  std::vector<Point<2> > vertices (n_vertices);
  for (unsigned int vertex=0; vertex<n_vertices; ++vertex) 
    {
      double x[3];

				       // read vertex
      in >> x[0] >> x[1] >> x[2];

				       // store vertex
      for (unsigned int d=0; d<2; ++d)
	vertices[vertex](d) = x[d];
    };
  AssertThrow (in, ExcIO());

				   // do some clean-up on vertices...
  delete_unused_vertices (vertices, cells, subcelldata);
				   // ... and cells
  GridReordering<2>::reorder_cells (cells);
  tria->create_triangulation (vertices, cells, subcelldata);
}

#endif // #if deal_II_dimension == 2



// 3-D XDA meshes
#if deal_II_dimension == 3

template <>
void GridIn<3>::read_xda (std::istream &in)
{
  Assert (tria != 0, ExcNoTriangulationSelected());
  AssertThrow (in, ExcIO());

  static const unsigned int xda_to_dealII_map[] = {0,1,5,4,3,2,6,7};
  
  std::string line;
				   // skip comments at start of file
  getline (in, line);


  unsigned int n_vertices;
  unsigned int n_cells;

				   // read cells, throw away rest of line
  in >> n_cells;
  getline (in, line);

  in >> n_vertices;
  getline (in, line);

				   // ignore following 8 lines
  for (unsigned int i=0; i<8; ++i)
    getline (in, line);

				   // set up array of cells
  std::vector<CellData<3> > cells (n_cells);
  SubCellData subcelldata;

  for (unsigned int cell=0; cell<n_cells; ++cell) 
    {
				       // note that since in the input
				       // file we found the number of
				       // cells at the top, there
				       // should still be input here,
				       // so check this:
      AssertThrow (in, ExcIO());
      Assert(GeometryInfo<3>::vertices_per_cell == 8,
	     ExcInternalError());
      
      unsigned int xda_ordered_nodes[8];
      
      for (unsigned int i=0; i<8; ++i)
	in >> xda_ordered_nodes[i];

      for (unsigned int i=0; i<8; i++)
	cells[cell].vertices[i] = xda_ordered_nodes[xda_to_dealII_map[i]];
    };


  
				   // set up array of vertices
  std::vector<Point<3> > vertices (n_vertices);
  for (unsigned int vertex=0; vertex<n_vertices; ++vertex) 
    {
      double x[3];

				       // read vertex
      in >> x[0] >> x[1] >> x[2];

				       // store vertex
      for (unsigned int d=0; d<3; ++d)
	vertices[vertex](d) = x[d];
    };
  AssertThrow (in, ExcIO());

				   // do some clean-up on vertices...
  delete_unused_vertices (vertices, cells, subcelldata);
				   // ... and cells
  GridReordering<3>::reorder_cells (cells);
  tria->create_triangulation (vertices, cells, subcelldata);
}

#endif // #if deal_II_dimension == 3



template <int dim>
void GridIn<dim>::read_msh (std::istream &in)
{
  Assert (tria != 0, ExcNoTriangulationSelected());
  AssertThrow (in, ExcIO());

  unsigned int n_vertices;
  unsigned int n_cells;
  unsigned int dummy;
  std::string line;
  
  getline(in, line);

  AssertThrow (line=="$NOD",
	       ExcInvalidGMSHInput(line));

  in >> n_vertices;

  std::vector<Point<dim> >     vertices (n_vertices);
				   // set up mapping between numbering
				   // in msh-file (nod) and in the
				   // vertices vector
  std::map<int,int> vertex_indices;
  
  for (unsigned int vertex=0; vertex<n_vertices; ++vertex) 
    {
      int vertex_number;
      double x[3];

				       // read vertex
      in >> vertex_number
	 >> x[0] >> x[1] >> x[2];
     
      for (unsigned int d=0; d<dim; ++d)
	vertices[vertex](d) = x[d];
				       // store mapping;
      vertex_indices[vertex_number] = vertex;
    };
  
				   // This is needed to flush the last
				   // new line
  getline (in, line);
  
				   // Now read in next bit
  getline (in, line);
  AssertThrow (line=="$ENDNOD",
	       ExcInvalidGMSHInput(line));

  
  getline (in, line);
  AssertThrow (line=="$ELM",
	       ExcInvalidGMSHInput(line));

  in >> n_cells;
				   // set up array of cells
  std::vector<CellData<dim> > cells;
  SubCellData                 subcelldata;

  for (unsigned int cell=0; cell<n_cells; ++cell) 
    {
				       // note that since in the input
				       // file we found the number of
				       // cells at the top, there
				       // should still be input here,
				       // so check this:
      AssertThrow (in, ExcIO());
      
/*  
    $ENDNOD
    $ELM
    NUMBER-OF-ELEMENTS
    ELM-NUMBER ELM-TYPE REG-PHYS REG-ELEM NUMBER-OF-NODES NODE-NUMBER-LIST
    ...
    $ENDELM
*/
      
      unsigned int cell_type;
      unsigned int material_id;
      unsigned int nod_num;
      
      in >> dummy          // ELM-NUMBER
	 >> cell_type	   // ELM-TYPE
	 >> material_id	   // REG-PHYS
	 >> dummy	   // reg_elm
	 >> nod_num;
      
/*       `ELM-TYPE'
	 defines the geometrical type of the N-th element:
	 `1'
	 Line (2 nodes, 1 edge).
                                                                                                
	 `3'
	 Quadrangle (4 nodes, 4 edges).
                                                                                                
	 `5'
	 Hexahedron (8 nodes, 12 edges, 6 faces).
                                                                                    
	 `15'
	 Point (1 node).
*/
      
      if (((cell_type == 1) && (dim == 1)) ||
	  ((cell_type == 3) && (dim == 2)) ||
	  ((cell_type == 5) && (dim == 3)))
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
	if ((cell_type == 1) && ((dim == 2) || (dim == 3)))
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
	  if ((cell_type == 3) && (dim == 3))
					     // boundary info
	    {
 	      subcelldata.boundary_quads.push_back (CellData<2>());
 	      in >> subcelldata.boundary_quads.back().vertices[0]
 	         >> subcelldata.boundary_quads.back().vertices[1]
 		 >> subcelldata.boundary_quads.back().vertices[2]
 		 >> subcelldata.boundary_quads.back().vertices[3];
	      
 	      subcelldata.boundary_quads.back().material_id = material_id;
	      
					       // transform from gmsh to
					       // consecutive numbering
 	      for (unsigned int i=0; i<4; ++i)
 	        if (vertex_indices.find (subcelldata.boundary_quads.back().vertices[i]) !=
 		    vertex_indices.end())
 		                                   // vertex with this index exists
		  subcelldata.boundary_quads.back().vertices[i]
		    = vertex_indices[subcelldata.boundary_quads.back().vertices[i]];
 	        else
 	          {
						     // no such vertex index
 	            Assert (false,
 		            ExcInvalidVertexIndex(cell,
 			                          subcelldata.boundary_quads.back().vertices[i]));
 		    subcelldata.boundary_quads.back().vertices[i] = -1;
 		  };
	      
	    }
	  else
					     // cannot read this
	    AssertThrow (false, ExcGmshUnsupportedGeometry(cell_type));
    };

  
				   // check that no forbidden arrays are used
  Assert (subcelldata.check_consistency(dim), ExcInternalError());

  AssertThrow (in, ExcIO());

				   // do some clean-up on vertices...
  delete_unused_vertices (vertices, cells, subcelldata);
				   // ... and cells
  GridReordering<dim>::reorder_cells (cells);
  tria->create_triangulation (vertices, cells, subcelldata);
}




template <int dim>
void GridIn<dim>::skip_empty_lines (std::istream &in)
{    
  std::string line;
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
}



template <int dim>
void GridIn<dim>::skip_comment_lines (std::istream &in,
				      const char    comment_start)
{
  char c;
				   // loop over the following comment
				   // lines
  while ((c=in.get()) == comment_start)
				     // loop over the characters after
				     // the comment starter
    while (in.get() != '\n');
  
				   // put back first character of
				   // first non-comment line
  in.putback (c);
}



template <int dim>
void
GridIn<dim>::delete_unused_vertices (std::vector<Point<dim> >    &vertices,
				     std::vector<CellData<dim> > &cells,
				     SubCellData                          &subcelldata)
{
				   // first check which vertices are
				   // actually used
  std::vector<bool> vertex_used (vertices.size(), false);
  for (unsigned int c=0; c<cells.size(); ++c)
    for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
      vertex_used[cells[c].vertices[v]] = true;

				   // then renumber the vertices that
				   // are actually used in the same
				   // order as they were beforehand
  const unsigned int invalid_vertex = deal_II_numbers::invalid_unsigned_int;
  std::vector<unsigned int> new_vertex_numbers (vertices.size(), invalid_vertex);
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
  std::vector<Point<dim> > tmp;
  tmp.reserve (count(vertex_used.begin(), vertex_used.end(), true));
  for (unsigned int v=0; v<vertices.size(); ++v)
    if (vertex_used[v] == true)
      tmp.push_back (vertices[v]);
  swap (vertices, tmp);
}



template <int dim>
void GridIn<dim>::debug_output_grid (const std::vector<CellData<dim> > &/*cells*/,
				     const std::vector<Point<dim> >    &/*vertices*/,
				     std::ostream                      &/*out*/)
{
  Assert (false, ExcNotImplemented());
}


#if deal_II_dimension == 2

template <>
void
GridIn<2>::debug_output_grid (const std::vector<CellData<2> > &cells,
			      const std::vector<Point<2> >    &vertices,
			      std::ostream                    &out)
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

      out << "# cell " << i << std::endl;
      Point<2> center;
      for (unsigned int f=0; f<4; ++f)
	center += vertices[cells[i].vertices[f]];
      center /= 4;

      out << "set label \"" << i << "\" at "
	  << center(0) << ',' << center(1)
	  << " center"
	  << std::endl;

				       // first two line right direction
      for (unsigned int f=0; f<2; ++f)
	out << "set arrow from "
	    << vertices[cells[i].vertices[f]](0) << ',' 
	    << vertices[cells[i].vertices[f]](1)
	    << " to "
	    << vertices[cells[i].vertices[(f+1)%4]](0) << ',' 
	    << vertices[cells[i].vertices[(f+1)%4]](1)
	    << std::endl;
				       // other two lines reverse direction
      for (unsigned int f=2; f<4; ++f)
	out << "set arrow from "
	    << vertices[cells[i].vertices[(f+1)%4]](0) << ',' 
	    << vertices[cells[i].vertices[(f+1)%4]](1)
	    << " to "
	    << vertices[cells[i].vertices[f]](0) << ',' 
	    << vertices[cells[i].vertices[f]](1)
	    << std::endl;
      out << std::endl;
    };
  

  out << std::endl
      << "set nokey" << std::endl
      << "pl [" << min_x << ':' << max_x << "]["
      << min_y << ':' << max_y <<  "] "
      << min_y << std::endl
      << "pause -1" << std::endl;
}

#endif


#if deal_II_dimension == 3

template <>
void
GridIn<3>::debug_output_grid (const std::vector<CellData<3> > &cells,
			      const std::vector<Point<3> >    &vertices,
			      std::ostream                    &out)
{
  for (unsigned int cell=0; cell<cells.size(); ++cell)
    {
				       // line 0
      out << vertices[cells[cell].vertices[0]]
	  << std::endl
	  << vertices[cells[cell].vertices[1]]
	  << std::endl << std::endl << std::endl;
				       // line 1
      out << vertices[cells[cell].vertices[1]]
	  << std::endl
	  << vertices[cells[cell].vertices[2]]
	  << std::endl << std::endl << std::endl;
				       // line 2
      out << vertices[cells[cell].vertices[3]]
	  << std::endl
	  << vertices[cells[cell].vertices[2]]
	  << std::endl << std::endl << std::endl;
				       // line 3
      out << vertices[cells[cell].vertices[0]]
	  << std::endl
	  << vertices[cells[cell].vertices[3]]
	  << std::endl << std::endl << std::endl;
				       // line 4
      out << vertices[cells[cell].vertices[4]]
	  << std::endl
	  << vertices[cells[cell].vertices[5]]
	  << std::endl << std::endl << std::endl;
				       // line 5
      out << vertices[cells[cell].vertices[5]]
	  << std::endl
	  << vertices[cells[cell].vertices[6]]
	  << std::endl << std::endl << std::endl;
				       // line 6
      out << vertices[cells[cell].vertices[7]]
	  << std::endl
	  << vertices[cells[cell].vertices[6]]
	  << std::endl << std::endl << std::endl;
				       // line 7
      out << vertices[cells[cell].vertices[4]]
	  << std::endl
	  << vertices[cells[cell].vertices[7]]
	  << std::endl << std::endl << std::endl;
				       // line 8
      out << vertices[cells[cell].vertices[0]]
	  << std::endl
	  << vertices[cells[cell].vertices[4]]
	  << std::endl << std::endl << std::endl;
				       // line 9
      out << vertices[cells[cell].vertices[1]]
	  << std::endl
	  << vertices[cells[cell].vertices[5]]
	  << std::endl << std::endl << std::endl;
				       // line 10
      out << vertices[cells[cell].vertices[2]]
	  << std::endl
	  << vertices[cells[cell].vertices[6]]
	  << std::endl << std::endl << std::endl;
				       // line 11
      out << vertices[cells[cell].vertices[3]]
	  << std::endl
	  << vertices[cells[cell].vertices[7]]
	  << std::endl << std::endl << std::endl;
    };
}

#endif

template <int dim>
void GridIn<dim>::read (const std::string& filename,
			Format format)
{
  std::ifstream in(filename.c_str());
  
  Assert (in.is_open(), ExcFileNotOpen(filename.c_str()));

//TODO:[GK] Do we really need the following, which does not compile with gcc 2.95?  
//   if (!in)
//     {
//       std::ios_base::iostate state = in.rdstate();
//       std::cerr << "File open, but error " << state << std::endl;
//       exit(1);
//     }
  
  read(in, format);
  in.close();
}


template <int dim>
void GridIn<dim>::read (std::istream& in,
			Format format)
{
  if (format == Default)
    format = default_format;
  
  switch (format)
    {
      case dbmesh:
	read_dbmesh (in);
	return;
	
      case msh:
	read_msh (in);
	return;
	
      case ucd:
	read_ucd (in);
	return;
	
      case xda:
	read_xda (in);
	return;

      case Default:
	break;
   }
  Assert (false, ExcInternalError());
}



template <int dim>
std::string
GridIn<dim>::default_suffix (const Format format) 
{
  switch (format) 
    {
      case dbmesh:
        return ".dbmesh";
      case msh: 
	return ".msh";
      case ucd: 
	return ".inp";
      case xda:
	return ".xda";
      default: 
	Assert (false, ExcNotImplemented()); 
	return ".unknown_format";
    }
}



template <int dim>
typename GridIn<dim>::Format
GridIn<dim>::parse_format (const std::string &format_name)
{
  if (format_name == "dbmesh")
    return dbmesh;

  if (format_name == "msh")
    return msh;
  
  if (format_name == "ucd")
    return ucd;

  if (format_name == "xda")
    return xda;

  AssertThrow (false, ExcInvalidState ());
				   // return something weird
  return Format(Default);
}



template <int dim>
std::string GridIn<dim>::get_format_names () 
{
  return "dbmesh|msh|ucd|xda";
}



//explicit instantiations
template class GridIn<deal_II_dimension>;
