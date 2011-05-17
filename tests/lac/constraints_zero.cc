//----------------------------  constraints_zero.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2006, 2007, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  constraints_zero.cc  ---------------------------


// generate the constraints for a case where there are nodes that have
// a constraint x[i]=0, i.e. where the right hand side is a trivial
// linear combination of other degrees of freedom. then print this set
// of constraints.
//
// we used to get this case wrong (we simply forgot to output this
// node).


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/constraint_matrix.h>

#include <fstream>


void test ()
{

  ConstraintMatrix cm;

				   // a "regular" constraint
  cm.add_line (1);
  cm.add_entry (1, 2, 42.);

				   // a "singular" constraint
  cm.add_line (4);
  
  cm.print (deallog.get_file_stream());
}


int main ()
{
  std::ofstream logfile("constraints_zero/output");
  logfile.precision(2);
  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);  

  test ();
  
  deallog << "OK" << std::endl;
}
