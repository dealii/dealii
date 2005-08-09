//----------------------------  subcelldata.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  subcelldata.cc  ---------------------------


#include "../tests.h"
#include <grid/tria.h>
#include <grid/grid_out.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>

#include <fstream>


template <int dim>
void test()
{
  Assert(dim==2 || dim==3, ExcNotImplemented());
  
  std::vector<Point<dim> > vertices (GeometryInfo<dim>::vertices_per_cell);
  vertices[0](0)=0;
  vertices[0](1)=0;
  vertices[1](0)=2;
  vertices[1](1)=1;
  vertices[2](0)=3;
  vertices[2](1)=3;
  vertices[3](0)=0;
  vertices[3](1)=1;
  if (dim==3)
    {
				   // for the new numbering
//       for (unsigned int i=0; i<4; ++i)
// 	{
// 	  vertices[i+4]=vertices[i];
// 	  vertices[i+4](2)=1;
// 	}
				       // for the old numbering
      for (unsigned int i=0; i<4; ++i)
	{
	  std::swap(vertices[i](1),vertices[i](2));
	  vertices[i+4]=vertices[i];
	  vertices[i+4](1)=1;
	}
      
    }

  std::vector<CellData<dim> > cells (1);
  for (unsigned int i=0;i<GeometryInfo<dim>::vertices_per_cell;++i)
    cells[0].vertices[i] = i;
  cells[0].material_id = 0;

  SubCellData subcelldata;
  if (dim==2)
    {
      subcelldata.boundary_lines.resize(GeometryInfo<dim>::faces_per_cell);
      for(unsigned int i=0; i<GeometryInfo<dim>::faces_per_cell; ++i)
	{
	  subcelldata.boundary_lines[i].vertices[0]=i;
	  subcelldata.boundary_lines[i].vertices[1]=(i+1)%4;
	  subcelldata.boundary_lines[i].material_id=10*i+1;
	}
    }
  else if (dim==3)
    {
      subcelldata.boundary_quads.resize(GeometryInfo<dim>::faces_per_cell);
      for(unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
	{
	  for(unsigned int i=0; i<GeometryInfo<dim>::vertices_per_face; ++i)
	    subcelldata.boundary_quads[f].vertices[i]=
	      GeometryInfo<dim>::face_to_cell_vertices(f,i);
	  subcelldata.boundary_quads[f].material_id=10*f+1;
	}
    }

  Triangulation<dim> tria;
  tria.create_triangulation (vertices, cells, subcelldata);

  GridOutFlags::Ucd ucd_flags(true,true);
  GridOut grid_out;
  grid_out.set_flags(ucd_flags);
  grid_out.write_ucd(tria, deallog.get_file_stream());

//   std::ofstream gnuplot_file("subcelldata.gnuplot");
//   grid_out.write_gnuplot(tria, gnuplot_file);
//   std::ofstream ucd_file("subcelldata.inp");
//   grid_out.write_ucd(tria, ucd_file);
}


int main ()
{
  std::ofstream logfile("subcelldata.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  
  test<2>();
  test<3>();
  
  return 0;
}
