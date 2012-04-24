//----------------------------  full_matrix_06.cc,v  ---------------------------
//    $Id$
//
//    Copyright (C) 2011, 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  full_matrix_06.cc,v  ---------------------------

// check method FullMatrix::scatter_matrix_to

#include "../tests.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iomanip>
#include <cstdlib>

#include <deal.II/base/logstream.h>
#include <deal.II/lac/full_matrix.h>

void test ()
{
				   // create a matrix with known
				   // elements
  FullMatrix<double> A(5,6);
  for (unsigned int i=0; i<A.m(); ++i)
    for (unsigned int j=0; j<A.n(); ++j)
      A(i,j) = i+j;

				   // pick every other row and column
  std::vector<unsigned int> rows (A.m());
  for (unsigned int i=0; i<rows.size(); ++i)
    rows[i] = 2*i;

  std::vector<unsigned int> cols (A.n());
  for (unsigned int i=0; i<cols.size(); ++i)
    cols[i] = 2*i;

				   // do the scatter
  FullMatrix<double> X(rows.size()*2, cols.size()*2);
  A.scatter_matrix_to (rows, cols, X);

				   // verify that the elements are
				   // correct
  for (unsigned int i=0; i<X.m(); ++i)
    for (unsigned int j=0; j<X.n(); ++j)
      if ((i % 2 == 0) && (j % 2 == 0))
	Assert (X(i,j) == i/2 + j/2,
		ExcInternalError())
	else
	Assert (X(i,j) == 0,
		ExcInternalError());

  deallog << "OK" << std::endl;
}


int
main ()
{
  const std::string logname = JobIdentifier::base_name(__FILE__) + std::string("/output");
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console(0);

  test();
}
