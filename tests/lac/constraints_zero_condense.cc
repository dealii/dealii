//----------------------------  constraints_zero_condense.cc  ---------------------------
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
//----------------------------  constraints_zero_condense.cc  ---------------------------


// generate the constraints for a case where there are nodes that have
// a constraint x[i]=0, i.e. where the right hand side is a trivial
// linear combination of other degrees of freedom. then check that if
// we condense a matrix with this then the right thing happens


#include "../tests.h"
#include <base/logstream.h>
#include <lac/constraint_matrix.h>
#include <lac/sparse_matrix.h>

#include <fstream>
#include <iomanip>


void test ()
{

				   // constrain each dof to zero. this
				   // should yield a diagonal matrix
  ConstraintMatrix cm;

  for (unsigned int i=0; i<5; ++i)
    cm.add_line (i);
  cm.close ();
  
				   // completely fill a 5x5 matrix
				   // with some values
  SparsityPattern sp (5,5,5);
  for (unsigned int i=0; i<5; ++i)
    for (unsigned int j=0; j<5; ++j)
      sp.add(i,j);
  sp.compress ();
  
  SparseMatrix<double> m(sp);
  for (unsigned int i=0; i<5; ++i)
    for (unsigned int j=0; j<5; ++j)
      m.set(i,j, i+j+2);

				   // now condense it
  cm.condense (m);

				   // and print it
  for (unsigned int i=0; i<5; ++i)
    {
      for (unsigned int j=0; j<5; ++j)
	deallog << m(i,j) << ' ';
      deallog << std::endl;
    }
}


int main ()
{
  std::ofstream logfile("constraints_zero_condense/output");
  deallog << std::setprecision (2);
  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);  

  test ();
  
  deallog << "OK" << std::endl;
}
