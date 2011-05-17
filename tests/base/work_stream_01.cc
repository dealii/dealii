//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

// test functions in namespace WorkStream

#include "../tests.h"
#include <iomanip>
#include <iomanip>
#include <fstream>
#include <cmath>

#include <deal.II/base/work_stream.h>


struct ScratchData
{};


struct CopyData
{
    unsigned int computed;
};


struct X
{
    void worker (const std::vector<unsigned int>::iterator &i,
		 ScratchData &,
		 CopyData &ad)
      {
	ad.computed = *i * 2;
      }

    void copier (const CopyData &ad)
      {
	deallog << ad.computed << std::endl;
      }
};


void test () 
{
  std::vector<unsigned int> v;
  for (unsigned int i=0; i<20; ++i)
    v.push_back (i);
  
  X x;
  WorkStream::run (v.begin(), v.end(), x, &X::worker, &X::copier,
		   ScratchData(),
		   CopyData());
}

  
  

int main()
{
  std::ofstream logfile("work_stream_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}
