/* $Id$ */

#include <basic/data_io.h>
#include <grid/dof.h>
#include <grid/dof_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria.h>
#include <grid/geometry_info.h>
#include <fe/fe.h>

#include <map>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <ctime>




bool SubCellData::check_consistency (const unsigned int dim) const {
  switch (dim) 
    {
      case 1:
	    return ((boundary_lines.size() == 0) &&
		    (boundary_quads.size() == 0));
      case 2:
	    return (boundary_quads.size() == 0);
    };
  return true;
};

		    



template <int dim>
DataIn<dim>::DataIn () :
		tria(0) {};



template <int dim>
void DataIn<dim>::attach_triangulation (Triangulation<dim> *t) {
  tria = t;
};



template <int dim>
void DataIn<dim>::read_ucd (istream &in) {
  Assert (tria != 0, ExcNoTriangulationSelected());
  Assert ((1<=dim) && (dim<=2), ExcNotImplemented());


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
	  Assert (false, ExcNotImplemented());
    };

				   // check that no forbidden arrays are used
  Assert (subcelldata.check_consistency(dim), ExcInternalError());
  
  tria->create_triangulation (vertices, cells, subcelldata);
};







template <int dim>
DataOut<dim>::DataEntry::DataEntry () :
		data(0), name(""), units("") {};



template <int dim>
DataOut<dim>::DataEntry::DataEntry (const dVector *data,
				    const string name,
				    const string units) :
			data(data), name(name), units(units) {};




template <int dim>
DataOut<dim>::DataOut () :
		dofs(0) {};



template <int dim>
void DataOut<dim>::attach_dof_handler (DoFHandler<dim> &d) {
  dofs = &d;
};



template <int dim>
void DataOut<dim>::add_data_vector (const dVector &vec,
				    const string  &name,
				    const string  &units) {
  DataEntry new_entry (&vec, name, units);
  data.push_back (new_entry);
};



template <int dim>
void DataOut<dim>::clear_data_vectors () {
  data.erase (data.begin(), data.end());
};



template <int dim>
void DataOut<dim>::write_ucd (ostream &out) const {
  Assert (dofs != 0, ExcNoDoFHandlerSelected());
  Assert (dofs->get_selected_fe().dofs_per_vertex==1,
	  ExcIncorrectDofsPerVertex());
  Assert ((1<=dim) && (dim<=3), ExcNotImplemented());
  
  DoFHandler<dim>::active_cell_iterator cell,
					endc = dofs->end();
  unsigned int n_vertex_dofs;

				   // first loop over all cells to
				   // find out how many degrees of
				   // freedom there are located on
				   // vertices
  if (true) 
    {
				       // block this to have local
				       // variables destroyed after
				       // use
      vector<bool> is_vertex_dof (dofs->n_dofs(), false);
      for (cell=dofs->begin_active(); cell!=endc; ++cell)
	for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell;
	     ++vertex) 
	  is_vertex_dof[cell->vertex_dof_index(vertex,0)] = true;

      n_vertex_dofs = 0;
      for (unsigned i=0; i!=is_vertex_dof.size(); ++i)
	if (is_vertex_dof[i] == true)
	  ++n_vertex_dofs;
    };
  

				   // write preamble
  if (true) 
    {
				       // block this to have local
				       // variables destroyed after
				       // use
      time_t  time1= time (0);
      tm     *time = localtime(&time1); 
      out << "# This file was generated by the deal.II library." << endl
	  << "# Date =  "
	  << time->tm_year+1900 << "/"
	  << time->tm_mon+1 << "/"
	  << time->tm_mday << endl
	  << "# Time =  "
	  << time->tm_hour << ":"
	  << setw(2) << time->tm_min << ":"
	  << setw(2) << time->tm_sec << endl
	  << "#" << endl
	  << "# For a description of the UCD format see the AVS Developer's guide."
	  << endl
	  << "#" << endl;
    };

				   // start with ucd data
  out << n_vertex_dofs << ' '
      << dofs->get_tria().n_active_cells() + n_boundary_faces() << ' '
      << data.size() << ' '
      << 0 << ' '                  // no cell data
      << 0                         // no model data
      << endl;
  
				   // write used nodes
  if (true) 
    {
				       // note if a given vertex was
				       // already written
      vector<bool> already_written (dofs->n_dofs(), false);
				       // write vertices
      for (cell=dofs->begin_active(); cell!=endc; ++cell)
	for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell;
	     ++vertex) 
	  if (already_written[cell->vertex_dof_index(vertex,0)]==false)
	    {
	      out << cell->vertex_dof_index(vertex,0) // vertex index
		  << "   "; 
	      for (unsigned int d=0; d<dim; ++d)      // vertex coordinates
		out << cell->vertex(vertex)(d) << ' ';
	      for (unsigned int d=dim; d<3; ++d)
		out << 0 << ' ';
	      out << endl;

	      already_written[cell->vertex_dof_index(vertex,0)] = true;
	    };
    };

				   // write cells. Enumerate cells
				   // consecutively (doesn't matter since
				   // we do not use cell data)
  if (true)
    {
      unsigned int index;
      for (cell=dofs->begin_active(), index=0; cell!=endc; ++cell, ++index)
	{
	  out << index << ' '
	      << static_cast<unsigned int>(cell->material_id())
	      << " ";
	  switch (dim) 
	    {
	      case 1:  out << "line    "; break;
	      case 2:  out << "quad    "; break;
	      case 3:  out << "hex     "; break;
	      default:
		    Assert (false, ExcNotImplemented());
	    };
	  for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell;
	       ++vertex)
	    out << cell->vertex_dof_index(vertex,0) << ' ';
	  out << endl;
	};

      write_ucd_faces (out, index);
    };

				   // if data given: write data, else
				   // only write grid
  if (data.size() != 0)
    {      
      out << data.size() << "    ";    // number of vectors
      for (unsigned int i=0; i!=data.size(); ++i)
	out << 1 << ' ';               // number of components;
				       // only 1 supported presently
      out << endl;
      
      for (unsigned int i=0; i<data.size(); ++i)
	out << data[i].name << ',' << data[i].units << endl;

      				       // AVS requires that the dof values
				       // be given in exactly the same order
				       // as in the vertex section and only
				       // once.

				       // note if a given vertex was
				       // already written
      vector<bool> already_written (dofs->n_dofs(), false);
				       // write vertices
      for (cell=dofs->begin_active(); cell!=endc; ++cell)
	for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell;
	     ++vertex) 
	  if (already_written[cell->vertex_dof_index(vertex,0)]==false)
	    {
	      out << cell->vertex_dof_index(vertex,0) // vertex index
		  << "   "; 
	      for (unsigned int i=0; i!=data.size(); ++i)
		out << (*data[i].data)(cell->vertex_dof_index(vertex,0)) << ' ';
	      out << endl;

	      already_written[cell->vertex_dof_index(vertex,0)] = true;
	    };
    };
				   // no cell data
				   // no model data
};



template <>
unsigned int DataOut<1>::n_boundary_faces () const {
  return 0;
};



template <int dim>
unsigned int DataOut<dim>::n_boundary_faces () const {
  typename DoFHandler<dim>::active_face_iterator face, endf;
  unsigned long int n_faces = 0;

  for (face=dofs->begin_active_face(), endf=dofs->end_face();
       face != endf; ++face)
    if ((face->boundary_indicator() != 0) &&
	(face->boundary_indicator() != 255))
      n_faces++;

  return n_faces;
};



template <>
void DataOut<1>::write_ucd_faces (ostream &, const unsigned int) const {
  return;
};



template <int dim>
void DataOut<dim>::write_ucd_faces (ostream &out,
				    const unsigned int starting_index) const {
  typename DoFHandler<dim>::active_face_iterator face, endf;
  unsigned int index=starting_index;

  for (face=dofs->begin_active_face(), endf=dofs->end_face();
       face != endf; ++face)
    if ((face->boundary_indicator() != 0) &&
	(face->boundary_indicator() != 255)) 
      {
	out << index << "  "
	    << static_cast<unsigned int>(face->boundary_indicator())
	    << "  ";
	switch (dim) 
	  {
	    case 2: out << "line    ";  break;
	    case 3: out << "quad    ";  break;
	    default:
		  Assert (false, ExcNotImplemented());
	  };
	for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_face; ++vertex)
	  out << face->vertex_dof_index(vertex,0) << ' ';
	out << endl;

	++index;
      };	  
};

      



template <int dim>
void DataOut<dim>::write_gnuplot (ostream &out) const {
  Assert (dofs != 0, ExcNoDoFHandlerSelected());
  Assert ((1<=dim) && (dim<=3), ExcNotImplemented());

  if (dim==1)
    Assert (dofs->get_selected_fe().dofs_per_vertex==1,
	    ExcIncorrectDofsPerVertex());

  				   // write preamble
  if (true) 
    {
				       // block this to have local
				       // variables destroyed after
				       // use
      time_t  time1= time (0);
      tm     *time = localtime(&time1); 
      out << "# This file was generated by the deal.II library." << endl
	  << "# Date =  "
	  << time->tm_year+1900 << "/"
	  << time->tm_mon+1 << "/"
	  << time->tm_mday << endl
	  << "# Time =  "
	  << time->tm_hour << ":"
	  << setw(2) << time->tm_min << ":"
	  << setw(2) << time->tm_sec << endl
	  << "#" << endl
	  << "# For a description of the UCD format see the AVS Developers guide."
	  << endl
	  << "#" << endl
	  << "# ";
      switch (dim) 
	{
	  case 1:
		out << "<x> ";
		break;
	  case 2:
		out << "<x> <y> ";
		break;
	  default:
		Assert (false, ExcNotImplemented());
	};
      for (unsigned int i=0; i!=data.size(); ++i)
	out << '<'
	    << data[i].name << ':' << data[i].units
	    << "> ";
      out << endl;
      
    };


  DoFHandler<dim>::active_cell_iterator cell;
  DoFHandler<dim>::active_cell_iterator endc = dofs->end();

  for (cell=dofs->begin_active(); cell!=endc; ++cell) 
    {
      switch (dim) 
	{
	  case 1:
						 // one dimension: write
						 // left vertex, right vertex
						 // and data values
		for (unsigned int vertex=0; vertex<2; ++vertex) 
		  {
		    out << cell->vertex(vertex)
			<< "   ";
		    for (unsigned int i=0; i!=data.size(); ++i)
		      out << (*data[i].data)(cell->vertex_dof_index(vertex,0))
			  << ' ';
		    out << endl;
		  };
		
		break;

	  case 2:
						 // two dimension: output
						 // grid and data as a sequence
						 // of lines in 3d
		for (unsigned int vertex=0; vertex<4; ++vertex) 
		  {
		    out << cell->vertex(vertex) << "   ";
		    for (unsigned int i=0; i!=data.size(); ++i)
		      out << (*data[i].data)(cell->vertex_dof_index(vertex,0))
			  << ' ';
		    out << endl;
		  };
						 // first vertex again
		out << cell->vertex(0) << "   ";
		for (unsigned int i=0; i!=data.size(); ++i)
		  out << (*data[i].data)(cell->vertex_dof_index(0,0))
		      << ' ';
		out << endl
		    << endl
		    << endl;      // end of cell; two newlines, since this
						 // stops continuous drawing
						 // of lines

		break;

	  default:
		Assert (false, ExcNotImplemented());
	};
    };
};

      

template <int dim>
void DataOut<dim>::write (ostream &out,
			  const OutputFormat output_format) const {
  switch (output_format) 
    {
      case ucd:
	    write_ucd (out);
	    break;
      case gnuplot:
	    write_gnuplot (out);
	    break;
      default:
	    Assert (false, ExcNotImplemented());
    };
};



template <int dim>
string DataOut<dim>::default_suffix (const OutputFormat output_format) {
  switch (output_format) 
    {
      case ucd:
	    return ".inp";
      case gnuplot:
	    return ".gnuplot";
      default:
	    Assert (false, ExcNotImplemented());
	    return "";
    };
};
  



//explicite instantiations
template class DataIn<1>;
template class DataIn<2>;
template class DataOut<1>;
template class DataOut<2>;
