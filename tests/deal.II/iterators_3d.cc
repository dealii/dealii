//----------------------------  iterators_3d.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  iterators_3d.cc  ---------------------------


// a test written in an attempt to figure out where I messed up
// converting some iterator functions...

#include "../tests.h"
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <grid/grid_generator.h>
#include <base/logstream.h>

#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>


std::ofstream logfile("iterators_3d/output");


template <int dim>
void test () 
{
  deallog << dim << std::endl;
  
  Triangulation<dim> tria;
  GridGenerator::hyper_cube (tria);
  tria.refine_global (2);

  {
    deallog << "Raw iterators" << std::endl;
    deallog << tria.begin_raw() << std::endl;
    if (dim > 1)
      deallog << tria.begin_raw_face() << std::endl;
    deallog << tria.last_raw() << std::endl;
    if (dim > 1)
      deallog << tria.last_raw_face() << std::endl;
    if (dim > 1)
      deallog << tria.end_raw_face() << std::endl;
  }
  {
    deallog << "Iterators" << std::endl;
    deallog << tria.begin() << std::endl;
    if (dim > 1)
      deallog << tria.begin_face() << std::endl;
    deallog << tria.last() << std::endl;
    if (dim > 1)
      deallog << tria.last_face() << std::endl;
    deallog << tria.end() << std::endl;
    if (dim > 1)
      deallog << tria.end_face() << std::endl;
  }
  {
    deallog << "Active iterators" << std::endl;
    deallog << tria.begin_active() << std::endl;
    if (dim > 1)
      deallog << tria.begin_active_face() << std::endl;
    deallog << tria.last_active() << std::endl;
    if (dim > 1)
      deallog << tria.last_active_face() << std::endl;
    if (dim > 1)
      deallog << tria.end_active_face() << std::endl;
  }

  for (unsigned int l=0; l<tria.n_levels(); ++l)
    {
      {
	deallog << "Raw iterators" << std::endl;
	deallog << tria.begin_raw(l) << std::endl;
	deallog << tria.last_raw(l) << std::endl;
	deallog << tria.end_raw(l) << std::endl;
      }
      {
	deallog << "Iterators" << std::endl;
	deallog << tria.begin(l) << std::endl;
	deallog << tria.last(l) << std::endl;
	deallog << tria.end(l) << std::endl;
      }
      {
	deallog << "Active iterators" << std::endl;
	deallog << tria.begin_active(l) << std::endl;
	deallog << tria.last_active(l) << std::endl;
	deallog << tria.end_active(l) << std::endl;
      }
    }
}


int main ()
{
  deallog << std::setprecision (2);
  logfile << std::setprecision (2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1> ();
  test<2> ();
  test<3> ();
}

