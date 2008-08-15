//----------------------------  refine_flag_set.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2006, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  refine_flag_set.cc  ---------------------------


// after the merge of branch_anisotropic,
// TriaAccessor::refine_flag_set did not just return a boolean, but a
// member of RefinementPossibilities. This broke the possibility of
// doing something like
//     cell->refine_flag_set() == true

#include "../tests.h"
#include <grid/tria_boundary.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/grid_tools.h>
#include <grid/grid_out.h>
#include <base/logstream.h>
#include <cmath>
#include <cstdlib>

#include <fstream>
#include <iomanip>
#include <cstdio>

std::ofstream logfile("refine_flag_set/output");



template <int dim>
void test ()
{
  deallog << dim << "D" << std::endl;
  
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.begin_active()->set_refine_flag();

  deallog << (int)(unsigned char)tria.begin_active()->refine_flag_set() << std::endl;
  
  Assert (tria.begin_active()->refine_flag_set() == true,
	  ExcInternalError());
  
  deallog << "OK" << std::endl;
}


int main ()
{
  deallog << std::setprecision(2);
  logfile << std::setprecision(2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  test<1> ();
  test<2> ();
  test<3> ();
  
  return 0;
}
