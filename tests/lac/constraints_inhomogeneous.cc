//-----------------------  constraints_inhomogeneous.cc  ----------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------  constraints_inhomogeneous.cc  ----------------------


// test: set inhomogeneous constraints and indirectly apply those to other
// constrained nodes.


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/constraint_matrix.h>

#include <fstream>


void test ()
{

  ConstraintMatrix cm;

				   // an inhomogeneous constraint
  cm.add_line (4);
  cm.set_inhomogeneity (4, 3.14159);

				   // a homogeneous constraint that is
				   // constrained to the inhomogeneous one
  cm.add_line (1);
  cm.add_entry (1, 2, 42.);
  cm.add_entry (1, 4, 1.);

				   // and a standard homogeneous constraint
  cm.add_line (17);
  cm.add_entry(17, 6, 2.);
  cm.add_entry(17, 15, 3.);

				   // a "singular" constraint
  cm.add_line (3);

				   // now close the constraint matrix
  cm.close();

  cm.print (deallog.get_file_stream());
}


int main ()
{
  std::ofstream logfile("constraints_inhomogeneous/output");
  logfile.precision(2);
  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);  

  test ();
  
  deallog << "OK" << std::endl;
}
