//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2006, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

// check FullMatrix::trace


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/full_matrix.h>

#include <fstream>
#include <iomanip>
#include <cmath>


int main()
{
  std::ofstream logfile("trace/output");
  deallog << std::fixed;
  deallog << std::setprecision(0);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  const unsigned int N=20;
  FullMatrix<double> m (N,N);

  double tr = 0;
  for (unsigned int i=0; i<N; ++i)
    for (unsigned int j=0; j<N; ++j)
      {
	m(i,j) = i+j;
	if (i==j)
	  tr += i+j;
      }

  deallog << "Trace=" << m.trace() << std::endl;
  Assert (m.trace() == tr, ExcInternalError());
}
