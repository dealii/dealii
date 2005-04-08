//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

// generate a hierarchical basis and display the values of the shape
// functions at equidistant points. (I needed this output at one point
// in time, so why not make it a testcase -- WB)

#include "../tests.h"
#include <fstream>
#include <cmath>

#include <base/logstream.h>
#include <base/polynomial.h>


using namespace Polynomials;


int main ()
{
  std::ofstream logfile("hierarchical.output");
  logfile.precision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  const std::vector<Polynomial<double> >
    p = Hierarchical::generate_complete_basis (10);

  const unsigned int div=30;
  for (unsigned int i=0; i<=div; ++i)
    {
      const double x = 1.*i/div;
      deallog << x << " ";
      for (unsigned int j=0; j<p.size(); ++j)
        deallog << p[j].value(x) << " ";
      deallog << std::endl;
    }
}

