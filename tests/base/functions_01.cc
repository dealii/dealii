//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

// Plot functions in library and check their consistency

#include "../tests.h"
#include <base/data_out_base.h>
#include <base/logstream.h>
#include <base/quadrature_lib.h>
#include <base/function_lib.h>
#include <base/flow_function.h>
#include <lac/vector.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <string>



template<int dim>
void
check_function(const Function<dim>& f,
	       unsigned int sub,
	       std::ostream& out)
{
				   // Prepare a vector with a single
				   // patch stretching over the cube
				   // [-1,1]^dim
  std::vector<DataOutBase::Patch<dim,dim> > patches(1);
  unsigned int vertex_number = 0;
  for (unsigned int iz=0;iz < ((dim>2) ? 2 : 1) ;++iz)
    for (unsigned int iy=0;iy < ((dim>1) ? 2 : 1) ;++iy)
      for (unsigned int ix=0;ix < 2 ;++ix)
	{
	  if (dim>0) patches[0].vertices[vertex_number](0) = -1. + 2.*ix;
	  if (dim>1) patches[0].vertices[vertex_number](1) = -1. + 2.*iy;
	  if (dim>2) patches[0].vertices[vertex_number](2) = -1. + 2.*iz;
	  ++vertex_number;
	}
  for (unsigned int i=0;i<GeometryInfo<dim>::faces_per_cell;++i)
    patches[0].neighbors[i] = deal_II_numbers::invalid_unsigned_int;
  patches[0].patch_index = 0;
  patches[0].n_subdivisions = sub;
  patches[0].points_are_available = false;

  vertex_number = 1;
  for (unsigned int d=0;d<dim;++d)
    vertex_number *= (sub+1);
  patches[0].data.reinit(f.n_components, vertex_number);

				   // Build the vector of quadrature points;
  std::vector<Point<dim> > points(vertex_number);
  const double h = 2./sub;
  vertex_number = 0;
  for (unsigned int iz=0;iz <= ((dim>2) ? sub : 0) ;++iz)
    for (unsigned int iy=0;iy <= ((dim>1) ? sub : 0) ;++iy)
      for (unsigned int ix=0;ix <= sub ;++ix)
	{
	  if (dim>0) points[vertex_number](0) = -1.+ix*h;
	  if (dim>1) points[vertex_number](1) = -1.+iy*h;
	  if (dim>2) points[vertex_number](2) = -1.+iz*h;
	  ++vertex_number;
	}
  
  std::vector<Vector<double> > values(points.size(), Vector<double>(f.n_components));
  f.vector_value_list(points, values);
  for (unsigned int i=0;i<values.size();++i)
    for (unsigned int j=0;j<values[i].size();++j)
      {
//	deallog << points[i] << '\t' << values[i](j) << std::endl;
	patches[0].data(j, i) = values[i](j);
      }
  
  std::vector<std::string> names(f.n_components);
  for (unsigned int i=0;i<names.size();++i)
    {
      names[i] = std::string("comp");
    }
  
  
  DataOutBase dout;
  DataOutBase::DXFlags flags;
  dout.write_dx(patches, names, flags, out);
}


int main()
{
  std::ofstream logfile("functions_01/output");

  if (true)
    {
      Functions::StokesLSingularity f;
      check_function(f, 4, logfile);
    }
  if (true)
    {
      Functions::PoisseuilleFlow<2> f(.8, 10.);
      check_function(f, 4, logfile);
    }
  if (true)
    {
      Functions::PoisseuilleFlow<3> f(.8, 10.);
      check_function(f, 4, logfile);
    }
}
