//----------------------------  geometry_info_5.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2003, 2004, 2005, 2006, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  geometry_info_5.cc  ---------------------------


// check GeometryInfo::d_linear_shape_function_gradient

#include "../tests.h"
#include <base/logstream.h>
#include <base/geometry_info.h>

#include <fstream>
#include <cstdlib>


template <int dim>
void test ()
{
  deallog << "Checking in " << dim << "d" << std::endl;

				   // check phi_i(v_j) = delta_{ij}
  for (unsigned int i=0;i<GeometryInfo<dim>::vertices_per_cell;++i)
    for (unsigned int v=0;v<GeometryInfo<dim>::vertices_per_cell;++v)
      {
	const Tensor<1,dim>
	  phi_i_grad
	  = GeometryInfo<dim>::d_linear_shape_function_gradient(GeometryInfo<dim>::unit_cell_vertex(v),i);
	
	deallog << phi_i_grad << std::endl;
      }

				   // check that
				   //    sum_i phi_i(x) == const
				   // at all points by verifying that the
				   // gradient of the sum of shape functions
				   // is zero. do so at every vertex, and then
				   // at the center
  for (unsigned int v=0;v<GeometryInfo<dim>::vertices_per_cell;++v)
    {
      Tensor<1,dim> s;
      for (unsigned int i=0;i<GeometryInfo<dim>::vertices_per_cell;++i)
	s += GeometryInfo<dim>::d_linear_shape_function_gradient(GeometryInfo<dim>::unit_cell_vertex(v),i);
      Assert (s.norm() == 0, ExcInternalError());

      deallog << "Sum of shape functions: " << s << std::endl;
    }
  {
    Point<dim> center;
    for (unsigned int i=0; i<dim; ++i)
      center[i] = 0.5;
    
    Tensor<1,dim> s;
    for (unsigned int i=0;i<GeometryInfo<dim>::vertices_per_cell;++i)
      s += GeometryInfo<dim>::d_linear_shape_function_gradient(center,i);
    Assert (s.norm() == 0, ExcInternalError());

    deallog << "Sum of shape functions: " << s << std::endl;
  }
}


int main () 
{
  std::ofstream logfile("geometry_info_5/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1> ();
  test<2> ();
  test<3> ();

  return 0;
}
