//----------------------------  geometry_info_1.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  geometry_info_1.cc  ---------------------------


// check GeometryInfo::cell_to_child and back

#include "../tests.h"
#include <base/logstream.h>
#include <base/geometry_info.h>

#include <fstream>
#include <cstdlib>

double rand_2 ()
{
  return 1.*rand()/RAND_MAX*4-2.;
}


template <int dim>
void test ()
{
				   // Output normal directions for each face
  for (unsigned int f=0;f<GeometryInfo<dim>::faces_per_cell;++f)
    {
      deallog << "Face " << f << ": n = ( ";
      for (unsigned int d=0;d<dim;++d)
	{
	  if (d != 0)
	    deallog << " , ";
	  if (d==GeometryInfo<dim>::unit_normal_direction[f])
	    deallog << GeometryInfo<dim>::unit_normal_orientation[f];
	  else
	    deallog << '0';
	}
      deallog << " )" << std::endl;
    }
  
  
  Point<dim> p;

				   // generate N random points in
				   // [-2:2]^d, and transform them
				   // back and forth between mother
				   // and child cell
  const unsigned int N = 7;
  for (unsigned int i=0; i<N; ++i)
    {
      for (unsigned int d=0; d<dim; ++d)
	p[d] = rand_2();

      deallog << i << ' ' << p << ' '
	      << GeometryInfo<dim>::is_inside_unit_cell (p) << std::endl;
      for (unsigned int c=0; c<GeometryInfo<dim>::children_per_cell; ++c)
	{
	  const Point<dim> q = GeometryInfo<dim>::cell_to_child_coordinates(p,c);
	  const Point<dim> pp = GeometryInfo<dim>::child_to_cell_coordinates(q,c);
	  
	  deallog << "    " << c << " [" << q << "] [" << pp << ']'
		  << std::endl;
	  Assert ((p-pp).square() < 1e-15*1e-15, ExcInternalError());
	  Assert (GeometryInfo<dim>::is_inside_unit_cell (p) ==
		  GeometryInfo<dim>::is_inside_unit_cell (pp),
		  ExcInternalError());
	}
    }
}


int main () 
{
  std::ofstream logfile("geometry_info_1.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1> ();
  test<2> ();
  test<3> ();

  return 0;
}
