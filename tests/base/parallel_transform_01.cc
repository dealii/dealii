//-----------------------------------------------------------------------------
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
//-----------------------------------------------------------------------------

// test parallel::transform

#include "../tests.h"
#include <iomanip>
#include <fstream>

#include <base/parallel.h>
#include <lac/vector.h>
#include <boost/lambda/lambda.hpp>



int main()
{
  std::ofstream logfile("parallel_transform_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  const unsigned int N=10000;
  Vector<double> x(N), y(N);

  for (unsigned int i=0; i<N; ++i)
    x(i) = i;

				   // set y=2*x
  parallel::transform (x.begin(), x.end(), y.begin(),
		       (2*boost::lambda::_1),
		       10);

				   // compute y=0 from the previous result
  y -= x;
  y -= x;

  Assert (y.l2_norm() == 0, ExcInternalError());

  deallog << "OK" << std::endl;
}
