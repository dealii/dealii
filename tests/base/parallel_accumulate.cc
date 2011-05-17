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

// test parallel::accumulate_from_subranges

#include "../tests.h"
#include <iomanip>
#include <fstream>

#include <deal.II/base/parallel.h>


int sum (const int begin,
	 const int end)
{
  int s=0;
  for (int i=begin; i<end; ++i)
    s += i;
  return s;
}
  

int main()
{
  std::ofstream logfile("parallel_accumulate/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  const int N = 10000;
  const int s = parallel::accumulate_from_subranges<int> (&sum, 0, N, 10);

  Assert (s == N*(N-1)/2, ExcInternalError());

  deallog << s << std::endl;
}
