//----------------------------  full_matrix_05.cc,v  ---------------------------
//    $Id$
//
//    Copyright (C) 2011, 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  full_matrix_05.cc,v  ---------------------------

// check method FullMatrix::extract_submatrix_from

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
  FullMatrix<double> A(10,12);
  for (unsigned int i=0; i<A.m(); ++i)
    for (unsigned int j=0; j<A.n(); ++j)
      A(i,j) = i+j;

				   // pick every other row and column
  std::vector<types::global_dof_index> rows (A.m()/2);
  for (unsigned int i=0; i<rows.size(); ++i)
    rows[i] = 2*i;

  std::vector<types::global_dof_index> cols (A.n()/2);
  for (unsigned int i=0; i<cols.size(); ++i)
    cols[i] = 2*i;

				   // do the extraction
  FullMatrix<double> X(rows.size(), cols.size());
  X.extract_submatrix_from (A, rows, cols);

				   // verify that the elements are
				   // correct
  for (unsigned int i=0; i<X.m(); ++i)
    for (unsigned int j=0; j<X.n(); ++j)
      Assert (X(i,j) == 2*i + 2*j,
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
