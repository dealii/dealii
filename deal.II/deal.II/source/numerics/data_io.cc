/* $Id$ */
/* Copyright W. Bangerth, Guido Kanschat, Stefan Nauber  */
/* University of Heidelberg, 1998, 1999                  */

#include <basic/data_io.h>
#include <grid/dof.h>
#include <grid/dof_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria.h>
#include <grid/geometry_info.h>
#include <base/quadrature_lib.h>
#include <fe/fe_values.h>
#include <fe/fe.h>
#include <lac/vector.h>

#include <map>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <ctime>





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







template <int dim>
DataOut<dim>::DataEntry::DataEntry () :
		data(0), name(""), units("") {};



template <int dim>
DataOut<dim>::DataEntry::DataEntry (const Vector<double> *data,
				    const string name,
				    const string units) :
			data(data), name(name), units(units) {};




template <int dim>
DataOut<dim>::DataOut () :
		dofs(0) {};



template <int dim>
void DataOut<dim>::attach_dof_handler (const DoFHandler<dim> &d) {
  dofs = &d;
};



template <int dim>
void DataOut<dim>::add_data_vector (const Vector<double> &vec,
				    const string  &name,
				    const string  &units) {
  Assert (dofs != 0, ExcNoDoFHandlerSelected ());
  Assert (name.find_first_not_of("abcdefghijklmnopqrstuvwxyz"
				 "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
				 "0123456789_<>()") == string::npos,
	  ExcInvalidCharacter (name));
  Assert (units.find_first_not_of("abcdefghijklmnopqrstuvwxyz"
				  "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
				  "0123456789_<>()") == string::npos,
	  ExcInvalidCharacter (units));
  
  DataEntry new_entry (&vec, name, units);
  if (vec.size() == dofs->n_dofs())
    dof_data.push_back (new_entry);
  else
    if (vec.size() == dofs->get_tria().n_active_cells())
      cell_data.push_back (new_entry);
    else
      Assert (false,
	      ExcInvalidVectorSize (vec.size(),
				    dofs->n_dofs(),
				    dofs->get_tria().n_active_cells()));
};



template <int dim>
void DataOut<dim>::clear_data_vectors () {
  dof_data.erase (dof_data.begin(), dof_data.end());
  cell_data.erase (cell_data.begin(), cell_data.end());
};



template <int dim>
void DataOut<dim>::write_ucd (ostream &out) const {
  Assert (dofs != 0, ExcNoDoFHandlerSelected());
  Assert (dofs->get_fe().dofs_per_vertex==1,
	  ExcIncorrectDofsPerVertex());
  Assert ((1<=dim) && (dim<=3), ExcNotImplemented());
  
  DoFHandler<dim>::active_cell_iterator cell,
					endc = dofs->end();
  unsigned int n_vertex_dofs = 0;

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

      for (unsigned i=0; i!=is_vertex_dof.size(); ++i)
	if (is_vertex_dof[i] == true)
	  ++n_vertex_dofs;
    };
  

				   // write preamble
  if (true)
    {
/*      
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
*/
    };

				   // start with ucd data
  out << n_vertex_dofs << ' '
      << dofs->get_tria().n_active_cells() + n_boundary_faces() << ' '
      << dof_data.size() << ' '
      << cell_data.size() << ' '
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

					   // it follows a list of the
					   // vertices of each cell. in 1d
					   // this is simply a list of the
					   // two vertices, in 2d its counter
					   // clockwise, as usual in this
					   // library. in 3d, the same applies
					   // (special thanks to AVS for
					   // numbering their vertices in a
					   // way compatible to deal.II!)
					   //
					   // technical reference:
					   // AVS Developer's Guide, Release 4,
					   // May, 1992, p. E6
	  for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell;
	       ++vertex)
	    out << cell->vertex_dof_index(vertex,0) << ' ';
	  out << endl;
	};

      write_ucd_faces (out, index);
    };

				   // if data given: write data, else
				   // only write grid
  if (dof_data.size() != 0)
    {      
      out << dof_data.size() << "    ";    // number of vectors
      for (unsigned int i=0; i!=dof_data.size(); ++i)
	out << 1 << ' ';               // number of components;
				       // only 1 supported presently
      out << endl;
      
      for (unsigned int i=0; i<dof_data.size(); ++i)
	out << dof_data[i].name << ',' << dof_data[i].units << endl;

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
	      for (unsigned int i=0; i!=dof_data.size(); ++i)
		out << (*dof_data[i].data)(cell->vertex_dof_index(vertex,0)) << ' ';
	      out << endl;

	      already_written[cell->vertex_dof_index(vertex,0)] = true;
	    };
    };

  if (cell_data.size() != 0)
    {
      out << cell_data.size() << "    ";    // number of vectors
      for (unsigned int i=0; i!=cell_data.size(); ++i)
	out << 1 << ' ';               // number of components;
				       // only 1 supported presently
      out << endl;

      for (unsigned int i=0; i<cell_data.size(); ++i)
	out << cell_data[i].name << ',' << cell_data[i].units << endl;

      unsigned int index;
      for (cell=dofs->begin_active(), index=0; cell!=endc; ++cell, ++index)
	{
	  out << index << "  ";
	  for (unsigned int i=0; i!=cell_data.size(); ++i)
	    out << (*cell_data[i].data)(index) << ' ';
	  out << endl;
	};

				       // strange enough, but true: the ucd
				       // format requires that the number of
				       // cell data entries be the same as
				       // there were cells. cells however
				       // include those boundary faces with
				       // an indicator other than zero, so
				       // we may have printed faces as well
				       // in the above list of cells. we
				       // have to give respective values
				       // here as well. since faces have no
				       // natural value when cell data is
				       // concerned, we assign a zero.
      for (unsigned int i=0; i<n_boundary_faces(); ++i, ++index)
	{
	  out << index << "  ";
	  for (unsigned int j=0; j!=cell_data.size(); ++j)
	    out << "0 ";
	  out << endl;
	};
      
    };
  
				   // no model data

				   // assert the stream is still ok
  AssertThrow (out, ExcIO());
};



#if deal_II_dimension == 1

template <>
unsigned int DataOut<1>::n_boundary_faces () const {
  return 0;
};

#endif



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




#if deal_II_dimension == 1

template <>
void DataOut<1>::write_ucd_faces (ostream &, const unsigned int) const {
  return;
};

#endif


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

  AssertThrow (out, ExcIO());
};

      



template <int dim>
void DataOut<dim>::write_gnuplot (ostream &out, unsigned int accuracy) const
{
  Assert (dofs != 0, ExcNoDoFHandlerSelected());
  Assert ((1<=dim) && (dim<=3), ExcNotImplemented());

  if (dim==1)
    Assert (dofs->get_fe().dofs_per_vertex==1,
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
	  << "# For a description of the GNUPLOT format see the GNUPLOT manual."
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
	  case 3:
		out << "<x> <y> <z>";
		break;
		
	  default:
		Assert (false, ExcNotImplemented());
	};
      for (unsigned int i=0; i!=dof_data.size(); ++i)
	out << '<'
	    << dof_data[i].name << ':' << dof_data[i].units
	    << "> ";
      for (unsigned int i=0; i!=cell_data.size(); ++i)
	out << '<'
	    << cell_data[i].name << ':' << cell_data[i].units
	    << "> ";
      out << endl;
      
    };


  DoFHandler<dim>::active_cell_iterator cell;
  DoFHandler<dim>::active_cell_iterator endc = dofs->end();

  QIteratedTrapez<dim> points(accuracy);
  
  FEValues<dim> fe(dofs->get_fe(), points, UpdateFlags(update_q_points),
		   dofs->get_tria().get_boundary());
  vector< vector <vector<double> > >
    values (dof_data.size(),
	    vector< vector<double> >(points.n_quadrature_points,
				     vector<double>(dofs->get_fe().n_components
				     )));

  unsigned int cell_index=0;
  for (cell=dofs->begin_active(); cell!=endc; ++cell, ++cell_index) 
    {
      fe.reinit(cell);

      for (unsigned i=0; i<dof_data.size(); ++i)
	fe.get_function_values(*dof_data[i].data, values[i]);
      
      unsigned supp_pt;
      
      switch (dim) 
	{
	  case 1:
						 // one dimension: write
						 // each vertex, right vertex
						 // and data values
		for (supp_pt = 0; supp_pt<points.n_quadrature_points; ++supp_pt) 
		  {
		    Point<dim> pt = fe.quadrature_point(supp_pt);
		    out << pt << "  ";
		    for (unsigned int i=0; i!=dof_data.size(); ++i)
		      for (unsigned int j=0; j < dofs->get_fe().n_components; ++j)
			out << values[i][supp_pt][j]
			    << ' ';
		    for (unsigned int i=0; i<cell_data.size(); ++i)
		      out << (*cell_data[i].data)(cell_index)
			  << ' ';
		    
		    out << endl;
		  };
		
		break;

	  case 2:
						 // two dimension: output
						 // grid and data as a sequence
						 // of points in 3d forming
						 // a tensor mesh on each cell
		for (unsigned xpt = 0, supp_pt = 0;
		     xpt <= accuracy; ++xpt)
		{
		  for(unsigned ypt = 0; ypt <= accuracy; ++ypt, ++supp_pt)
		  {
		    Point<dim> pt = fe.quadrature_point(supp_pt);
		    out << pt << "  ";
		    
		    for (unsigned int i=0; i!=dof_data.size(); ++i)
		      for (unsigned int j=0; j < dofs->get_fe().n_components; ++j)
			out << values[i][supp_pt][j]
			    << ' ';
		    for (unsigned int i=0; i<cell_data.size(); ++i)
		      out << (*cell_data[i].data)(cell_index)
			  << ' ';
		    out << endl;
		  }
		  out << endl;
		}
		
		break;

	  case 3:
						 // three dimension: output
						 // grid and data as a tensor
						 // mesh in 4d. we very
						 // obviously have a problem
						 // here since humans can watch
						 // 4d lines properly. but one
						 // can at least view cuts
						 // through this hypercube,
						 // which is possible this
						 // way
		
		for (unsigned xpt = 0, supp_pt = 0;
		     xpt <= accuracy; ++xpt)
		  {
		    for(unsigned ypt = 0; ypt <= accuracy; ++ypt)
		      {
			for(unsigned zpt = 0; zpt <= accuracy; ++zpt, ++supp_pt)
			  {
			    Point<dim> pt = fe.quadrature_point(supp_pt);
			    out << pt << "  ";
			    
			    for (unsigned int i=0; i!=dof_data.size(); ++i)
			      for (unsigned int j=0; j < dofs->get_fe().n_components; ++j)
				out << values[i][supp_pt][j]
				    << ' ';
			    for (unsigned int i=0; i<cell_data.size(); ++i)
			      out << (*cell_data[i].data)(cell_index)
				  << ' ';
			    out << endl;
			  }
			out << endl;
		      }
		    out << endl;
		  };
		
		
		break;

	  default:
		Assert (false, ExcNotImplemented());
	}
    out << endl;
    }

  AssertThrow (out, ExcIO());
}



template <int dim>
void DataOut<dim>::write_gnuplot_draft (ostream &out) const
{
  Assert (dofs != 0, ExcNoDoFHandlerSelected());
  Assert ((1<=dim) && (dim<=3), ExcNotImplemented());

  if (dim==1)
    Assert (dofs->get_fe().dofs_per_vertex==1,
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
	  << "# For a description of the GNUPLOT format see the GNUPLOT manual."
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
	  case 3:
		out << "<x> <y> <z>";
		break;
	  default:
		Assert (false, ExcNotImplemented());
	};
      for (unsigned int i=0; i!=dof_data.size(); ++i)
	out << '<'
	    << dof_data[i].name << ':' << dof_data[i].units
	    << "> ";
      out << endl;
      
    };


  DoFHandler<dim>::active_cell_iterator cell;
  DoFHandler<dim>::active_cell_iterator endc = dofs->end();

  unsigned int cell_index=0;
  for (cell=dofs->begin_active(); cell!=endc; ++cell, ++cell_index) 
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
		    for (unsigned int i=0; i!=dof_data.size(); ++i)
		      out << (*dof_data[i].data)(cell->vertex_dof_index(vertex,0))
			  << ' ';
		    for (unsigned int i=0; i<cell_data.size(); ++i)
		      out << (*cell_data[i].data)(cell_index)
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
		    for (unsigned int i=0; i!=dof_data.size(); ++i)
		      out << (*dof_data[i].data)(cell->vertex_dof_index(vertex,0))
			  << ' ';
		    for (unsigned int i=0; i<cell_data.size(); ++i)
		      out << (*cell_data[i].data)(cell_index)
			  << ' ';
		    out << endl;
		  };
						 // first vertex again
		out << cell->vertex(0) << "   ";
		for (unsigned int i=0; i!=dof_data.size(); ++i)
		  out << (*dof_data[i].data)(cell->vertex_dof_index(0,0))
		      << ' ';
		for (unsigned int i=0; i<cell_data.size(); ++i)
		  out << (*cell_data[i].data)(cell_index)
		      << ' ';
		out << endl
		    << endl
		    << endl;      // end of cell; two newlines, since this
						 // stops continuous drawing
						 // of lines

		break;

	  case 3:
						 // three dimension: output
						 // grid and data as a sequence
						 // of lines in 4d

						 // first we plot the front face
		for (unsigned int vertex=0; vertex<4; ++vertex) 
		  {
		    out << cell->vertex(vertex) << "   ";
		    for (unsigned int i=0; i!=dof_data.size(); ++i)
		      out << (*dof_data[i].data)(cell->vertex_dof_index(vertex,0))
			  << ' ';
		    for (unsigned int i=0; i<cell_data.size(); ++i)
		      out << (*cell_data[i].data)(cell_index)
			  << ' ';
		    out << endl;
		  };
						 // first vertex again
		out << cell->vertex(0) << "   ";
		for (unsigned int i=0; i!=dof_data.size(); ++i)
		  out << (*dof_data[i].data)(cell->vertex_dof_index(0,0))
		      << ' ';
		for (unsigned int i=0; i<cell_data.size(); ++i)
		  out << (*cell_data[i].data)(cell_index)
		      << ' ';
		out << endl
		    << endl
		    << endl;      // end of front face; two newlines, since this
						 // stops continuous drawing
						 // of lines

						 // first we back the front face
		for (unsigned int vertex=4; vertex<8; ++vertex) 
		  {
		    out << cell->vertex(vertex) << "   ";
		    for (unsigned int i=0; i!=dof_data.size(); ++i)
		      out << (*dof_data[i].data)(cell->vertex_dof_index(vertex,0))
			  << ' ';
		    for (unsigned int i=0; i<cell_data.size(); ++i)
		      out << (*cell_data[i].data)(cell_index)
			  << ' ';
		    out << endl;
		  };
						 // first vertex again
		out << cell->vertex(4) << "   ";
		for (unsigned int i=0; i!=dof_data.size(); ++i)
		  out << (*dof_data[i].data)(cell->vertex_dof_index(4,0))
		      << ' ';
		for (unsigned int i=0; i<cell_data.size(); ++i)
		  out << (*cell_data[i].data)(cell_index)
		      << ' ';
		out << endl
		    << endl
		    << endl;      // end of back face; two newlines, since this
						 // stops continuous drawing
						 // of lines

						 // finally, we plot a
						 // continuous line along
						 // the vertices 0, 4, 7, 3, 2,
						 // 6, 5, 1 to show the other
						 // four lines
		for (unsigned int vertex=0; vertex<8; ++vertex)
		  {
		    static const unsigned int vertex_list[8]
		      = { 0, 4, 7, 3, 2, 6, 5, 1 };
		    
		    out << cell->vertex(vertex_list[vertex]) << "   ";
		    for (unsigned int i=0; i!=dof_data.size(); ++i)
		      out << (*dof_data[i].data)(cell->vertex_dof_index(vertex_list[vertex],0))
			  << ' ';
		    for (unsigned int i=0; i<cell_data.size(); ++i)
		      out << (*cell_data[i].data)(cell_index)
			  << ' ';
		    out << endl;
		  };
						 // again: stop drawing
		out << endl << endl;
		
		
		break;
		
	  default:
		Assert (false, ExcNotImplemented());
	};
    };

  AssertThrow (out, ExcIO());
};


#if deal_II_dimension == 2

template <>
void DataOut<2>::write_povray_mesh (ostream &out) const {
  Assert (dofs != 0, ExcNoDoFHandlerSelected());

  				   // write preamble
  if (true) 
    {
				       // block this to have local
				       // variables destroyed after
				       // use
      time_t  time1= time (0);
      tm     *time = localtime(&time1); 
      out << "/*" << endl
	  << "This file was generated by the deal.II library." << endl
	  << "Date =  "
	  << time->tm_year+1900 << "/"
	  << time->tm_mon+1 << "/"
	  << time->tm_mday << endl
	  << "Time =  "
	  << time->tm_hour << ":"
	  << setw(2) << time->tm_min << ":"
	  << setw(2) << time->tm_sec << endl
	  << endl
	  << "For a description of the POVRAY format see the POVRAY manual."
	  << endl
	  << endl
	  << "*/"
	  << endl << endl;
    };

				   // write list of include files
  out << "#include \"colors.inc\""   << endl
      << endl;

				   // declare standard texture
  out << "#declare default_texture = texture {"         << endl
      << "        pigment { color White }" << endl
      << "        finish  { ambient 0.2 diffuse 0.6 specular 0.5 }" << endl
      << "}" << endl << endl;				

				   // write camera and light sources
  out << "camera {"                  << endl
      << "  location <8, 10, -20>"    << endl
      << "  angle 7" << endl
      << "  look_at  <0., 0., 0.>" << endl
      << "  sky  <0., 0., 1.>" << endl
      << "}"                         << endl
      << endl
      << "light_source {"            << endl
      << "  <20, 20, -20> White"        << endl
      << "}"                         << endl
      << endl;
      
      
				       // write frame object
  out << "mesh {" << endl;
  
  DoFHandler<2>::active_cell_iterator cell;
  DoFHandler<2>::active_cell_iterator endc = dofs->end();

  for (cell=dofs->begin_active(); cell!=endc; ++cell) 
    {
				       // write cell as two triangles
				       // y and z coordinates are switched
      out << "  triangle { ";
      out << '<' << cell->vertex(0)(0) << ','
	  << (*dof_data[0].data)(cell->vertex_dof_index(0,0)) << ','
	  << cell->vertex(0)(1) << '>'
	  << ", "
	  << '<' << cell->vertex(1)(0) << ','
	  << (*dof_data[0].data)(cell->vertex_dof_index(1,0)) << ','
	  << cell->vertex(1)(1) << '>'
	  << ", "
	  << '<' << cell->vertex(2)(0) << ','
	  << (*dof_data[0].data)(cell->vertex_dof_index(2,0)) << ','
	  << cell->vertex(2)(1) << '>'
	  << endl
	  << "              texture { default_texture }"   << endl
	  << "           }"
	  << endl;
      out << "  triangle { ";
      out << '<' << cell->vertex(0)(0) << ','
	  << (*dof_data[0].data)(cell->vertex_dof_index(0,0)) << ','
	  << cell->vertex(0)(1) << '>'
	  << ", "
	  << '<' << cell->vertex(2)(0) << ','
	  << (*dof_data[0].data)(cell->vertex_dof_index(2,0)) << ','
	  << cell->vertex(2)(1) << '>'
	  << ", "
	  << '<' << cell->vertex(3)(0) << ','
	  << (*dof_data[0].data)(cell->vertex_dof_index(3,0)) << ','
	  << cell->vertex(3)(1) << '>'
	  << endl
	  << "              texture { default_texture }"   << endl
	  << "           }"
	  << endl;
    };
  out << "}";     

  AssertThrow (out, ExcIO());
};


#endif



#if deal_II_dimension == 3

template <>
void DataOut<3>::write_povray_mesh (ostream &) const {
  Assert (false, ExcNotImplemented());
};

#endif



#if deal_II_dimension == 2

template <>
void DataOut<2>::write_epsgrid (ostream &out) const {
  Assert (dofs != 0, ExcNoDoFHandlerSelected());
  
  // write preamble
  if (true) 
    {
      // block this to have local
      // variables destroyed after
      // use
      time_t  time1= time (0);
      tm     *time = localtime(&time1); 
      out << "%!PS-Adobe-2.0 EPSF-1.2" << endl
	  << "%%Title: Deal Output" << endl
	  << "%%Creator: the deal.II library" << endl
	  << "%%Creation Date: " 
	  << time->tm_year+1900 << "/"
	  << time->tm_mon+1 << "/"
	  << time->tm_mday << " - "
	  << time->tm_hour << ":"
	  << setw(2) << time->tm_min << ":"
	  << setw(2) << time->tm_sec << endl
	  << "%%BoundingBox: 0 0 310 310" << endl;
    };  

  // Get scaling factors for Bounding Box 310 x 310
  
  DoFHandler<2>::active_cell_iterator cell;
  DoFHandler<2>::active_cell_iterator endc = dofs->end();
  double x, y, xmin=0, xmax=0, ymin=0, ymax=0, scale, xofs, yofs;
  int i;

  for (cell=dofs->begin_active(); cell!=endc; ++cell)
    {
      for (i=0; i<4; i++)
	{
	  x=cell->vertex(i)(0);
	  y=cell->vertex(i)(1);
	  xmin = ( x < xmin ? x : xmin );
	  xmax = ( x > xmax ? x : xmax );
	  ymin = ( y < ymin ? y : ymin );
	  ymax = ( y > ymax ? y : ymax );
	}
    }
  x = xmax - xmin;
  y = ymax - ymin;
  scale = 300 / (x > y ? x : y);
  xofs = -(xmin*scale)+5;
  yofs = -(ymin*scale)+5;

  for (cell=dofs->begin_active(); cell!=endc; ++cell) 
    {
      out << (cell->vertex(0)(0))*scale+xofs << " " 
	  << (cell->vertex(0)(1))*scale+yofs << " moveto "
	  << (cell->vertex(1)(0))*scale+xofs << " " 
	  << (cell->vertex(1)(1))*scale+yofs << " lineto "
	  << (cell->vertex(2)(0))*scale+xofs << " " 
	  << (cell->vertex(2)(1))*scale+yofs << " lineto "
	  << (cell->vertex(3)(0))*scale+xofs << " " 
	  << (cell->vertex(3)(1))*scale+yofs << " lineto "
	  << (cell->vertex(0)(0))*scale+xofs << " " 
	  << (cell->vertex(0)(1))*scale+yofs << " lineto "
	  << " closepath stroke" << endl;
    };
  out << "showpage" << endl;     

  AssertThrow (out, ExcIO());
};



template <>
void DataOut<2>::write_eps (ostream &out) const {
  Assert (dofs != 0, ExcNoDoFHandlerSelected());
  
  // write preamble
  if (true) 
    {
      // block this to have local
      // variables destroyed after
      // use
      time_t  time1= time (0);
      tm     *time = localtime(&time1); 
      out << "%!PS-Adobe-2.0 EPSF-1.2" << endl
	  << "%%Title: Deal Output" << endl
	  << "%%Creator: the deal.II library" << endl
	  << "%%Creation Date: " 
	  << time->tm_year+1900 << "/"
	  << time->tm_mon+1 << "/"
	  << time->tm_mday << " - "
	  << time->tm_hour << ":"
	  << setw(2) << time->tm_min << ":"
	  << setw(2) << time->tm_sec << endl
	  << "%%BoundingBox: -220 -261 220 450" << endl;
    };  

  // Get scaling factors for Bounding Box 310 x 310
  
  DoFHandler<2>::active_cell_iterator cell;
  DoFHandler<2>::active_cell_iterator endc = dofs->end();

  double x, y, z, xmin=0, xmax=0, ymin=0, ymax=0, zmin=0, zmax=0, zscale, scale, xofs, yofs;
  double cx[4], cy[4], cz[4];
  int i;

  for (cell=dofs->begin_active(); cell!=endc; ++cell)
    {
      for (i=0; i<4; i++)
	{
	  x=cell->vertex(i)(0);
	  y=cell->vertex(i)(1);
	  z=(*dof_data[0].data)(cell->vertex_dof_index(i,0));
	  xmin = ( x < xmin ? x : xmin );
	  xmax = ( x > xmax ? x : xmax );
	  ymin = ( y < ymin ? y : ymin );
	  ymax = ( y > ymax ? y : ymax );
	  zmin = ( z < zmin ? z : zmin );
	  zmax = ( z > zmax ? z : zmax );
	}
    }
  x = xmax - xmin;
  y = ymax - ymin;
  z = zmax - zmin;
  zscale = (x>y ? x : y)/z;
  scale = 300 / (x > y ? x : y);
  xofs = -(xmin*scale)+5;
  yofs = -(ymin*scale)+5;


  for (cell=dofs->begin_active(); cell!=endc; ++cell) 
    {
      cx[0]=cell->vertex(0)(0); cy[0]=cell->vertex(0)(1); cz[0]=(*dof_data[0].data)(cell->vertex_dof_index(0,0))*zscale;
      cx[1]=cell->vertex(1)(0); cy[1]=cell->vertex(1)(1); cz[1]=(*dof_data[0].data)(cell->vertex_dof_index(1,0))*zscale;;
      cx[2]=cell->vertex(2)(0); cy[2]=cell->vertex(2)(1); cz[2]=(*dof_data[0].data)(cell->vertex_dof_index(2,0))*zscale;;
      cx[3]=cell->vertex(3)(0); cy[3]=cell->vertex(3)(1); cz[3]=(*dof_data[0].data)(cell->vertex_dof_index(3,0))*zscale;;
      
      // Turn and scale

      for (i=0;i<4;i++)
	{
	  // x =  0.707 * cx[i] - 0.707 * cy[i] + 0.000 * cz[i];
	  // y =  0.354 * cx[i] + 0.354 * cy[i] + 0.866 * cz[i];
	  // z = -0.559 * cx[i] - 0.559 * cy[i] + 0.500 * cz[i];

	  x = 0.707 * (cx[i] - cy[i]);
	  y = 0.354 * (cx[i] + cy[i]) + 0.866 * cz[i];

	  cx[i]=x*scale+xofs;
	  cy[i]=y*scale+yofs;
	}

      out << cx[0] << " " << cy[0] << " moveto "
	  << cx[1] << " " << cy[1] << " lineto "
	  << cx[2] << " " << cy[2] << " lineto "
	  << cx[3] << " " << cy[3] << " lineto "
	  << cx[0] << " " << cy[0] << " lineto "
	  << " closepath stroke" << endl;
    };
  out << "showpage" << endl;     

  AssertThrow (out, ExcIO());
};

#endif


#if deal_II_dimension == 3

template <>
void DataOut<3>::write_epsgrid (ostream &/*out*/) const {
  Assert (false, ExcNotImplemented());
};


template <>
void DataOut<3>::write_eps (ostream &/*out*/) const {
  Assert (false, ExcNotImplemented());
};

#endif


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
    case gnuplot_draft:
      write_gnuplot_draft (out);
      break;
    case povray_mesh:
      write_povray_mesh (out);
      break;
    case eps:
      write_eps(out);
      break;
    case epsgrid:
      write_epsgrid(out);
      break;
    default:
      Assert (false, ExcNotImplemented());
    };
};



template <int dim>
string DataOut<dim>::default_suffix (const OutputFormat output_format) 
{
  switch (output_format) 
    {
    case ucd:
      return ".inp";
    case gnuplot: 
    case gnuplot_draft: 
      return ".gnuplot";
    case povray_mesh: 
      return ".pov";
    case eps: 
      return ".eps";
    case epsgrid: 
      return ".eps";
    default: 
      Assert (false, ExcNotImplemented()); 
      return "";
    };
};
  


template <int dim>
DataOut<dim>::OutputFormat
DataOut<dim>::parse_output_format (const string format_name) {
  if (format_name == "ucd")
    return ucd;

  if (format_name == "gnuplot")
    return gnuplot;

  if (format_name == "gnuplot draft")
    return gnuplot_draft;

  if (format_name == "povray mesh")
    return povray_mesh;

  if (format_name == "eps")
    return eps;

  if (format_name == "epsgrid")
    return epsgrid;

  AssertThrow (false, ExcInvalidState ());
};



template <int dim>
string DataOut<dim>::get_output_format_names () {
  return "ucd|gnuplot|gnuplot draft|povray mesh|eps|epsgrid";
};



//explicite instantiations

template class DataIn<deal_II_dimension>;
template class DataOut<deal_II_dimension>;
