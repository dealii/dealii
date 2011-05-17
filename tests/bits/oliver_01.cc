//----------------------------  oliver_01.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005 by the deal.II authors and Oliver Kayser-Herold
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  oliver_01.cc  ---------------------------


// Oliver found an example, where sparse_matrix_iterator->value=0 didn't work,
// because the iterator->value expects a double on the right hand side, not an
// integer. If the right hand side is zero, it can also be converted to a
// pointer, which leads to an ambiguity. Fix this by having an additional
// operator= in the iterator/reference class

#include "../tests.h"
#include <deal.II/lac/sparse_matrix.h>
#include <iomanip>
#include <fstream>


int main () 
{
  std::ofstream logfile("oliver_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

                                   // this test only needs to compile, not run
  if (false)
    {
      SparseMatrix<double>::iterator *i;
      (*i)->value () = (int)0;
    }

  deallog << "OK" << std::endl;
  
  return 0;
}
