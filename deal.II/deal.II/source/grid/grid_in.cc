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
  char c;
  while (in.get(c), c=='#') 
    {
      char line[256];
      in.get (line, 255, '\n'); // ignore rest of line, at most 256 chars
      in.get (c);         // ignore '\n' at end of line.
    };
  
				   // put back first character of
				   // first non-comment line
  in.putback (c);


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


//explicit instantiations
template class GridIn<deal_II_dimension>;
