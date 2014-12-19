// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// test WorkStream with empty functions

#include "../tests.h"
#include <iomanip>
#include <iomanip>
#include <fstream>
#include <cmath>

#include <deal.II/base/work_stream.h>


struct ScratchData
{};


void foo (const std::vector<unsigned int>::iterator,
          ScratchData&,
          unsigned int&)
{
}

void bar (const unsigned int&)
{
}

void test ()
{
  std::vector<unsigned int> v;
  for (unsigned int i=0; i<20; ++i)
    v.push_back (i);

  // first run with only a worker
  WorkStream::run (v.begin(), v.end(),
                   &foo,
                   std_cxx11::function<void(const unsigned int&)>(),
                   ScratchData(),
                   0U);
  
  // next run with only a copier
  WorkStream::run (v.begin(), v.end(),
                   std_cxx11::function<void(const std::vector<unsigned int>::iterator,
                                            ScratchData&,unsigned int&)>(),
                   &bar,
                   ScratchData(),
                   0U);
}




int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}
