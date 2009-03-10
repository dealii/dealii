//----------------------------  constraints_zero_merge.cc  ---------------------------
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
//----------------------------  constraints_zero_merge.cc  ---------------------------


// generate the constraints for a case where there are nodes that have
// a constraint x[i]=0, i.e. where the right hand side is a trivial
// linear combination of other degrees of freedom. then merge two such
// constraint matrices


#include "../tests.h"
#include <base/logstream.h>
#include <lac/constraint_matrix.h>

#include <fstream>


void test ()
{

  ConstraintMatrix cm1, cm2;

				   // a "regular" and a singular
				   // constraint to each of the
				   // constraint matrices
  cm1.add_line (1);
  cm1.add_entry (1, 2, 42.);
  cm1.add_line (4);

  cm2.add_line (11);
  cm2.add_entry (11, 22, 42.);
  cm2.add_line (14);

  cm2.merge (cm1);

  deallog << "CM1" << std::endl;
  cm1.print (deallog.get_file_stream());

  deallog << "CM2" << std::endl;
  cm2.print (deallog.get_file_stream());
}


int main ()
{
  std::ofstream logfile("constraints_zero_merge/output");
  logfile.precision(2);
  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);  

  test ();
  
  deallog << "OK" << std::endl;
}
