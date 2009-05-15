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
  std::ofstream logfile("parallel_transform_02/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  const unsigned int N=10000;
  Vector<double> x(N), y(N), z(N);

  for (unsigned int i=0; i<N; ++i)
    {
      x(i) = 2.*i;
      y(i) = -1.*i;
    }

				   // set z=x+2y, which happens to be zero
  parallel::transform (x.begin(), x.end(),
		       y.begin(),
		       z.begin(),
		       (boost::lambda::_1 + 2*boost::lambda::_2),
		       10);

  Assert (z.l2_norm() == 0, ExcInternalError());

  deallog << "OK" << std::endl;
}
