/* $Id$ */
/* Copyright W. Bangerth, Guido Kanschat, Stefan Nauber  */
/* University of Heidelberg, 1998, 1999                  */

#include <numerics/data_io.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
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
#include <set>
#include <cmath>



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
DataOut_Old<dim>::DataEntry::DataEntry () :
		data(0), name(""), units("") {};



template <int dim>
DataOut_Old<dim>::DataEntry::DataEntry (const Vector<double> *data,
				    const string name,
				    const string units) :
			data(data), name(name), units(units) {};




template <int dim>
DataOut_Old<dim>::DataOut_Old () :
		dofs(0) {};



template <int dim>
void DataOut_Old<dim>::attach_dof_handler (const DoFHandler<dim> &d) {
  dofs = &d;
};



template <int dim>
void DataOut_Old<dim>::add_data_vector (const Vector<double> &vec,
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
void DataOut_Old<dim>::clear_data_vectors () {
  dof_data.erase (dof_data.begin(), dof_data.end());
  cell_data.erase (cell_data.begin(), cell_data.end());
};



template <int dim>
void DataOut_Old<dim>::write_ucd (ostream &out) const {
  Assert (dofs != 0, ExcNoDoFHandlerSelected());
  Assert (dofs->get_fe().dofs_per_vertex==1,
	  ExcIncorrectDofsPerVertex());
  Assert ((1<=dim) && (dim<=3), ExcNotImplemented());
  
  DoFHandler<dim>::active_cell_iterator       cell;
  const DoFHandler<dim>::active_cell_iterator endc = dofs->end();
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
				   // consecutively
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
unsigned int DataOut_Old<1>::n_boundary_faces () const {
  return 0;
};

#endif



template <int dim>
unsigned int DataOut_Old<dim>::n_boundary_faces () const {
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
void DataOut_Old<1>::write_ucd_faces (ostream &, const unsigned int) const {
  return;
};

#endif


template <int dim>
void DataOut_Old<dim>::write_ucd_faces (ostream &out,
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
void DataOut_Old<dim>::write_gnuplot (ostream &out, unsigned int accuracy) const
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

  QTrapez<1>     q_trapez;
  QIterated<dim> points (q_trapez, accuracy);
  
  FEValues<dim> fe(dofs->get_fe(), points, UpdateFlags(update_q_points));
  vector< vector <Vector<double> > >
    values (dof_data.size(),
	    vector< Vector<double> >(points.n_quadrature_points,
				     Vector<double>(dofs->get_fe().n_components()
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
		      for (unsigned int j=0; j < dofs->get_fe().n_components(); ++j)
			out << values[i][supp_pt](j)
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
		      for (unsigned int j=0; j < dofs->get_fe().n_components(); ++j)
			out << values[i][supp_pt](j)
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
			      for (unsigned int j=0; j < dofs->get_fe().n_components(); ++j)
				out << values[i][supp_pt](j)
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
void DataOut_Old<dim>::write_gnuplot_draft (ostream &out) const
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
void DataOut_Old<2>::write_povray_mesh (ostream &out) const {
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



template <int dim>
void DataOut_Old<dim>::write_povray_mesh (ostream &) const {
				   // this is for all other dimensions that
				   // are not explicitely specialized
  Assert (false, ExcNotImplemented());
};




#if deal_II_dimension == 2

template <>
void DataOut_Old<2>::write_eps (ostream &out, const EpsOutputData &eod) const {
  Assert (dofs != 0, ExcNoDoFHandlerSelected());

  {
				     // write preamble
				     // block this to have local
				     // variables destroyed after
				     // use
    time_t  time1= time (0);
    tm     *time = localtime(&time1); 
    out << "%!PS-Adobe-2.0 EPSF-1.2" << endl
	<< "%%Title: deal.II Output" << endl
	<< "%%Creator: the deal.II library" << endl
	<< "%%Creation Date: " 
	<< time->tm_year+1900 << "/"
	<< time->tm_mon+1 << "/"
	<< time->tm_mday << " - "
	<< time->tm_hour << ":"
	<< setw(2) << time->tm_min << ":"
	<< setw(2) << time->tm_sec << endl
	<< "%%BoundingBox: -5 -5 305 305" << endl;
  };  

   
				    // Make output values local by
				    // copying them to a multiset.
				    // Perform the necessary turn.
   const DoFHandler<2>::active_cell_iterator endc = dofs->end();
   multiset<DataOut_Old<2>::EpsCellData> cells;
   multiset<DataOut_Old<2>::EpsCellData> cells2;
   
   bool height_data_p = (
                          ((dof_data.size())>0) 
			  && 
			  (
			    ( eod.height_info == EpsOutputData::DefaultHeight)
			    || 
			    ( eod.height_info == EpsOutputData::HeightVector)
			  )
			);

				    // Cells are colored, if there is
				    // cell data and the mode is
				    // ShadingVector or DefaultShading
   bool cell_data_p = (
     ((cell_data.size())>0) 
     && 
     (
       (eod.cell_shading == EpsOutputData::ShadingVector)
       ||
       (eod.cell_shading == EpsOutputData::DefaultShading)
     )
   );

				    // Cells are shaded, i.e. light
				    // shading, if they are not
				    // colored,  there is height
				    // information and the mode is
				    // Default or Light
   bool cell_shade_p = (
     (!cell_data_p)
     &&
     (height_data_p)
     &&
     (
       (eod.cell_shading == EpsOutputData::DefaultShading)
       ||
       (eod.cell_shading == EpsOutputData::LightShaded)
     )
   );

   unsigned cell_index;
   DoFHandler<2>::active_cell_iterator cell;
   for(cell_index=0, cell=dofs->begin_active(); 
       cell!=endc; 
       ++cell, ++cell_index)
     {
       EpsCellData cd;
       for (unsigned int i=0; i<GeometryInfo<2>::vertices_per_cell; ++i)
	 {
	   (cd.vertices[i]).x=cell->vertex(i)(0);
	   (cd.vertices[i]).y=cell->vertex(i)(1);
	   (cd.vertices[i]).z=(height_data_p ? 
			       (*dof_data[eod.height_vector].data)(cell->vertex_dof_index(i,0))
			       : 0); 
	 };
       if (height_data_p)
	 cd.turn(eod.azimuth,eod.elevation);

       if (cell_data_p)
	 cd.red=(*cell_data[eod.cell_vector].data)(cell_index);
       cells.insert(cd);
     };

				    // Now we proceed with the
				    // multiset cells. First we look
				    // for extrema.
   
   double xmin=cells.begin()->vertices[0].x;
   double xmax=xmin;
   double ymin=cells.begin()->vertices[0].y;
   double ymax=ymin;
   float cell_vector_min=cells.begin()->red; 
   float cell_vector_max=cell_vector_min;

   for(multiset<DataOut_Old<2>::EpsCellData>::iterator c=cells.begin();
       c!=cells.end(); ++c, ++cell_index)
     {
       for (unsigned int i=0; i<4; ++i)
	 {
	   double xvv,yvv;
	   xvv=c->vertices[i].x;
	   xmin=(xmin < xvv ? xmin : xvv);
	   xmax=(xmax > xvv ? xmax : xvv);
	   
	   yvv=c->vertices[i].y;
	   ymin=(ymin < yvv ? ymin : yvv);
	   ymax=(ymax > yvv ? ymax : yvv);
	 }
       if (cell_data_p) 
	 {
	   double cvv;
	   cvv = c->red;
	   cell_vector_max = (cell_vector_max > cvv ? cell_vector_max : cvv);
	   cell_vector_min = (cell_vector_min < cvv ? cell_vector_min : cvv);
	}
    };
   cells2.clear();


				    // If we want shaded output we can
				    // do the shading now.
   if (cell_shade_p)
     {
       double spann1[3], spann2[3], normal[3];
       double light_norm, normal_norm;
       float color;

       for (multiset<DataOut_Old<2>::EpsCellData>::iterator c=cells.begin();c!=cells.end();++c)
	 {
	   EpsCellData cd(*c);

	   spann1[0]=spann2[0]=cd.vertices[0].x;
	   spann1[1]=spann2[1]=cd.vertices[0].y;
	   spann1[2]=spann2[2]=cd.vertices[0].z;

	   spann1[0]-=cd.vertices[1].x;
	   spann1[1]-=cd.vertices[1].y;
	   spann1[2]-=cd.vertices[1].z;
		      
	   spann2[0]-=cd.vertices[2].x;
	   spann2[1]-=cd.vertices[2].y;
	   spann2[2]-=cd.vertices[2].z;

	   normal[0] = spann1[1]*spann2[2]-spann1[2]*spann2[1];
	   normal[1] = spann1[2]*spann2[0]-spann1[0]*spann2[2];
	   normal[2] = spann1[0]*spann2[1]-spann1[1]*spann2[0];

	   normal_norm = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
	   light_norm = sqrt(eod.light[0]*eod.light[0]+eod.light[1]*eod.light[1]+eod.light[2]*eod.light[2]);

	   color = eod.light[0]*normal[0]+eod.light[1]*normal[1]+eod.light[2]*normal[2];
	   color /= light_norm * normal_norm;
	   
	   cd.red=color;
	   cd.green=color;
	   cd.blue=color;

	   cells2.insert(cd);
	 };

					// since we don't need #cells# any
					// more, delete it now
       cells.clear ();
     }
   else 
				      // copy #cells# to #cells2#. since
				      // we don't need #cells# any
				      // more, we use a trick for copying
				      // that is significantly faster
     cells2.swap (cells);
   


				    // Next we have to shift and scale
				    // a bit so that everything still
				    // arrives in our bounding box of
				    // 310x310.
				    // If cell_data_p we also scale
				    // this so that it is in the range
				    // between 0 and 1.

   const double scale = 300 / (xmax-xmin > ymax-ymin ? xmax-xmin : ymax-ymin);
   
   for (multiset<DataOut_Old<2>::EpsCellData>::iterator c=cells2.begin();
	c!=cells2.end(); ++c)
     {
       EpsCellData cd (*c);
       for (unsigned int i=0; i<4; ++i)
	 {
	   cd.vertices[i].x=(cd.vertices[i].x-xmin)*scale;
	   cd.vertices[i].y=(cd.vertices[i].y-ymin)*scale;
	 };

       if (cell_data_p)
	 eod.color(cd.red,cell_vector_max,cell_vector_min,cd.red,cd.green,cd.blue);

       cells.insert(cd);
     };


				    //  Now we are ready to output...
   for (multiset<DataOut_Old<2>::EpsCellData>::iterator c=cells.begin();
	c!=cells.end(); ++c)
     {
       if (cell_data_p || cell_shade_p)
	 {
	   out << c->red << " " << c->green << " " << c->blue << " setrgbcolor "
	       << c->vertices[0].x << " " << c->vertices[0].y << " moveto "
	       << c->vertices[1].x << " " << c->vertices[1].y << " lineto "
	       << c->vertices[2].x << " " << c->vertices[2].y << " lineto "
	       << c->vertices[3].x << " " << c->vertices[3].y << " lineto "
	       << " closepath fill" << endl;
	 };

       if (eod.cell_boundary_shading != EpsOutputData::NoBoundary)
	 {
	   switch (eod.cell_boundary_shading)
	     {
	       case EpsOutputData::BlackBoundary:
	       case EpsOutputData::DefaultBoundary:
		     out << "0";
		     break;
	       case EpsOutputData::WhiteBoundary:
		     out << "1";
		     break;
	       case EpsOutputData::NoBoundary:
		     break;
	     };
	   out << " setgray " 
	       << c->vertices[0].x << " " << c->vertices[0].y << " moveto "
	       << c->vertices[1].x << " " << c->vertices[1].y << " lineto "
	       << c->vertices[2].x << " " << c->vertices[2].y << " lineto "
	       << c->vertices[3].x << " " << c->vertices[3].y << " lineto closepath stroke" << endl;
	 };
     };
   out << "showpage" << endl;
};

#endif




template <int dim>
void DataOut_Old<dim>::write_eps (ostream &,
			      const EpsOutputData &) const{
				   // this is for all other dimensions that
				   // are not explicitely specialized
  Assert (false, ExcNotImplemented());
};




template <int dim>
void DataOut_Old<dim>::write_gmv (ostream &out) const
{
				   // this function is mostly copied from
				   // the ucd format function

  Assert (dofs != 0, ExcNoDoFHandlerSelected());
  Assert (dofs->get_fe().dofs_per_vertex==1,
	  ExcIncorrectDofsPerVertex());
  Assert ((1<=dim) && (dim<=3), ExcNotImplemented());

  const unsigned int dofs_per_vertex = dofs->get_fe().dofs_per_vertex;
  
				   ///////////////////////
				   // preamble
  out << "gmvinput ascii"
      << endl
      << endl;


				   ///////////////////////////////
				   // first make up a list of used
				   // vertices along with their
				   // coordinates. since dofs have
				   // no separate numbering in
				   // gmv format, we need to generate
				   // a separate numbering as well
				   // for those dofs on vertices.
				   // note that we can't simply use
				   // the numbers of the vertices,
				   // since this numbering may not
				   // always be continuous if the
				   // triangulation has undergone
				   // coarsening somewhen (well, we
				   // could place all unused vertices
				   // to the origin, but that's not
				   // really elegant...)
				   //
				   // therefore: have a list with
				   // one slot for each dof, where
				   // the value indicates the
				   // respective number of the vertex
				   // if the dof is on a vertex,
				   // and -1 otherwise.
				   //
				   // also create a list of vertices
				   // storing their coordinates in
				   // the order defined by the
				   // above map, which we write to
				   // the file immediately after
				   // creation.
				   //
				   // note: since GMV seems to be
				   // made by Fortranists, the count
				   // their indices from 1 onwards,
				   // so whenever we actually output
				   // a vertex number, we should add
				   // a one.
  DoFHandler<dim>::active_cell_iterator       cell;
  const DoFHandler<dim>::active_cell_iterator endc = dofs->end();
  
  vector<int> dof_to_vertex_map (dofs->n_dofs(), -1);
  unsigned int used_vertices = 0;
  if (true)
    {
      vector<Point<dim> > vertices (dofs->get_tria().n_used_vertices());
      unsigned int next_free_vertex = 0;

      for (cell=dofs->begin_active(); cell!=endc; ++cell)
	for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
					   // check whether we
					   // already treated this vertex
	  if (dof_to_vertex_map[cell->vertex_dof_index(vertex,0)] == -1)
	    {
	      vertices[next_free_vertex] = cell->vertex(vertex);

					       // associate all dofs on
					       // this vertex with the
					       // respective vertex
	      for (unsigned int d=0; d<dofs_per_vertex; ++d)
		dof_to_vertex_map[cell->vertex_dof_index(vertex,d)] = next_free_vertex;

	      ++next_free_vertex;
	    };
      Assert (next_free_vertex == vertices.size(),
	      ExcInternalError());
      used_vertices = vertices.size();

				       // now write out the vertices in
				       // this order
      out << "nodes " << vertices.size() << endl;
      for (unsigned int component=0; component<dim; ++component)
	{
	  for (unsigned int v=0; v<vertices.size(); ++v)
	    out << vertices[v](component) << ' ';
	  out << endl;
	};

				       // and write the missing components
				       // y (for 1d) and z (for 1d and 2d)
				       // by simply setting them to zero
      for (unsigned int d=dim+1; d<=3; ++d)
	{
	  fill_n (ostream_iterator<double>(out, " "), vertices.size(), 0.0);
	  out << endl;
	}
            
      out << endl;
    };
  

				   /////////////////////////////////////
				   // now for the cells. this is simpler
				   // than the above task
  if (true)
    {
      out << "cells " << dofs->get_tria().n_active_cells() << endl;

      const char *cell_description[3] = { "line 2\n  ",
					  "quad 4\n  ",
					  "hex 8\n  "};
      
      for (cell=dofs->begin_active(); cell!=endc; ++cell)
	{
	  out << cell_description[dim-1];
					   // output vertex indices,
					   // counted from 1 onwards
	  for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
	    out << dof_to_vertex_map[cell->vertex_dof_index(v,0)]+1 << ' ';
	  out << endl << endl;
	};

      out << endl;
    };


				   ///////////////////////////////////////
				   // data output.
  out << "variable" << endl;

				   // first go for node data
				   //
				   // since here as with the vertex
				   // coordinates the order is a bit
				   // unpleasant (first all data of
				   // variable 1, then variable 2, etc)
				   // we have to copy them a bit around
				   //
				   // note that we copy vectors when
				   // looping over the cells since we
				   // have to write them one variable
				   // at a time and don't want to use
				   // more than one loop
  if (true)
    {
      vector<vector<double> > data_vectors (dof_data.size(),
					    vector<double> (used_vertices));

				       // loop over all cells and copy
				       // the data into the other vector
				       // note  if a vertex has already
				       // been visited
      vector<bool> vertex_copied (dofs->n_dofs(), false);
      for (cell=dofs->begin_active(); cell!=endc; ++cell)
	for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
	  if (vertex_copied[cell->vertex_dof_index(v,0)] == false)
	    {
					       // global dof number
	      const unsigned int dof_index = cell->vertex_dof_index(v,0);
	      
					       // this is this vertex's
					       // number to GMV
	      const int vertex_index = dof_to_vertex_map[dof_index];
	      Assert (vertex_index >= 0, ExcInternalError());
	      
	      for (unsigned vec=0; vec<dof_data.size(); ++vec)
		data_vectors[vec][vertex_index] = (*dof_data[vec].data)(dof_index);

	      vertex_copied[dof_index] = true;
	    };
      
				       // indicator 0 = node data 
      for (unsigned int vec=0; vec<dof_data.size(); ++vec)
	{
	  out << dof_data[vec].name << " 1" << endl;
	  copy(data_vectors[vec].begin(),
	       data_vectors[vec].end(),
	       ostream_iterator<double>(out, " "));
	  out << endl
	      << endl;
	};
    };


				   // Cell data is the simplest since the order
				   // already is correct
  if (true)
    {
      for (unsigned int vec=0; vec<cell_data.size(); ++vec)
	{
	  out << cell_data[vec].name << " 0" << endl;
	  copy((*cell_data[vec].data).begin(),
	       (*cell_data[vec].data).end(),
	       ostream_iterator<double>(out, " "));
	  out << endl
	      << endl;
	};
    };
  
				   // end of variable section
  out << "endvars" << endl;
  
				   // end of output
  out << "endgmv"
      << endl;
  
				   // assert the stream is still ok
  AssertThrow (out, ExcIO());
};

  



template <int dim>
void DataOut_Old<dim>::write (ostream &out,
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
	    
      case gmv:
	    write_gmv (out);
	    break;
	    
      default:
	    Assert (false, ExcNotImplemented());
    };
};



template <int dim>
string DataOut_Old<dim>::default_suffix (const OutputFormat output_format) 
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

      case gmv:
	    return ".gmv";
	    
      default: 
	    Assert (false, ExcNotImplemented()); 
	    return "";
    };
};
  


template <int dim>
DataOut_Old<dim>::OutputFormat
DataOut_Old<dim>::parse_output_format (const string &format_name) {
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

  if (format_name == "gmv")
    return gmv;
  
  AssertThrow (false, ExcInvalidState ());

				   // return something invalid
  return OutputFormat(-1);
};


template <int dim>
string DataOut_Old<dim>::get_output_format_names () {
  return "ucd|gnuplot|gnuplot draft|povray mesh|eps|gmv";
};



template<int dim>
bool DataOut_Old<dim>::EpsCellData::operator < (const EpsCellData &other) const
{
  double maxz = vertices[0].z, 
         othermaxz = other.vertices[0].z;
  unsigned i;

  for (i=1; i<4; ++i)
    { 
      maxz = (maxz > vertices[i].z ? maxz : vertices[i].z);
      othermaxz = (othermaxz > other.vertices[i].z ? othermaxz : other.vertices[i].z);
    };

  return maxz > othermaxz;
};



template <int dim>
void DataOut_Old<dim>::EpsVertexData::turn(double azi, double ele)
{
  double nx,ny,nz;

  double cx=cos(ele), cz=cos(azi), sx=sin(ele), sz=sin(azi);

  nx = -   cz*x+   sz*y;
  ny = -cx*sz*x-cx*cz*y-sx*z;
  nz = -sx*sz*x-sx*cz*y+cx*z;

  x=nx; z=ny; y=nz;
};

//      ( 1 0    0 )
// Dx = ( 0 cx -sx )
//      ( 0 sx  cx )

//      ( cy 0 sy )
// Dy = (  0 1  0 )
//      (-sy 0 cy )

//      ( cz -sz 0 )
// Dz = ( sz  cz 0 )
//      (  0   0 1 )

//       ( cz -sz 0 )( 1 0    0 )(x)   ( cz*x-sz*(cx*y-sx*z)+0*(sx*y+cx*z) )
// Dxz = ( sz  cz 0 )( 0 cx -sx )(y) = ( sz*x+cz*(cx*y-sx*z)+0*(sx*y+cx*z) )
//	 (  0   0 1 )( 0 sx  cx )(z)   (  0*x+	*(cx*y-sx*z)+1*(sx*y+cx*z) )



template <int dim>
void DataOut_Old<dim>::EpsCellData::turn(double azi, double ele)
{
  for (unsigned i=0; i<4; ++i)
    vertices[i].turn(azi,ele);
};


EpsOutputData::EpsOutputData()
		: height_info(DefaultHeight),
		  cell_shading(DefaultShading),
		  cell_boundary_shading (DefaultBoundary),
                  height_vector(0),
		  cell_vector(0),
		  azimuth(2*3.1415926* (180 - 30)/360),
		  elevation(2*3.1415926* (90-60)/360)
{ 
  light[0]=-1;
  light[1]=-1;
  light[2]=1;
};



void EpsOutputData::color(const float x,
			  const float xmax,
			  const float xmin, 
			  float &r,
			  float &g,
			  float &b) const
{
// A difficult color scale:
//     xmin          = black  (1)
// 3/4*xmin+1/4*xmax = blue   (2)
// 1/2*xmin+1/2*xmax = green  (3)
// 1/4*xmin+3/4*xmax = red    (4)
//              xmax = white  (5)
// Makes the following color functions:
//
// red      green    blue
//       __
//      /      /\  /  /\    /
// ____/    __/  \/  /  \__/

//     { 0                                (1) - (3)
// r = { ( 4*x-2*xmin+2*xmax)/(xmax-xmin) (3) - (4)
//     { 1                                (4) - (5)
//
//     { 0                                (1) - (2)
// g = { ( 4*x-3*xmin-  xmax)/(xmax-xmin) (2) - (3)
//     { (-4*x+  xmin+3*xmax)/(xmax-xmin) (3) - (4)
//     { ( 4*x-  xmin-3*xmax)/(xmax-xmin) (4) - (5)
//
//     { ( 4*x-4*xmin       )/(xmax-xmin) (1) - (2)
// b = { (-4*x+2*xmin+2*xmax)/(xmax-xmin) (2) - (3)
//     { 0                                (3) - (4)
//     { ( 4*x-  xmin-3*xmax)/(xmax-xmin) (4) - (5)

  float sum   =   xmax+  xmin;
  float sum13 =   xmin+3*xmax;
  float sum22 = 2*xmin+2*xmax;
  float sum31 = 3*xmin+  xmax;
  float dif = xmax-xmin;
  float rezdif = 1.0/dif;

  int where;

  if (x<(sum31)/4)
    where = 0;
  else if (x<(sum22)/4)
    where = 1;
  else if (x<(sum13)/4)
    where = 2;
  else
    where = 3;

  if (dif!=0)
    {
      switch (where)
	{
	  case 0:
		r=0;		      g=0;		          b=(x-xmin)*4.*rezdif;
		break;
	  case 1:
		r=0;                  g=(4*x-3*xmin-xmax)*rezdif; b=(sum22-4.*x)*rezdif;
		break;
	  case 2:
		r=(4*x-2*sum)*rezdif; g=(xmin+3*xmax-4*x)*rezdif; b=0;
		break;
	  case 3:
		r=1;                  g=(4*x-xmin-3*xmax)*rezdif; b=(4.*x-sum13)*rezdif;
	  default:
		break;
	};
    }
  else // White 
    {
      r=1;
      g=1;
      b=1;
    };
};



//explicit instantiations

template class DataIn<deal_II_dimension>;
template class DataOut_Old<deal_II_dimension>;
