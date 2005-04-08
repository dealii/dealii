//----------------------------  fe_collection_01a.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  fe_collection_01a.cc  ---------------------------

// FECollection would throw an exception when its destructor is called. Check
// that this is fixed

#include "../tests.h"
#include <base/logstream.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/grid_generator.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_handler.h>
#include <fe/fe_collection.h>
#include <fe/fe_q.h>

#include <fstream>
#include <iostream>

template <int dim>
void
check ()
{
  {
    const FE_Q<dim> fe_1(1);
    const FE_Q<dim> fe_2(2);
    
    FECollection<dim> fc;
    fc.add_fe (fe_1);
    fc.add_fe (fe_2);
  }
  deallog << dim << "d: OK" << std::endl;
}

int main()
{
  std::ofstream logfile("fe_collection_01.output");
  logfile.precision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  check<1> ();
  check<2> ();
  check<3> ();
}
