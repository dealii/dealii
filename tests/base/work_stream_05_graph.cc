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


// like _05, but with graph coloring

#include "../tests.h"
#include <iomanip>
#include <iomanip>
#include <fstream>
#include <cmath>

#include <deal.II/base/work_stream.h>
#include <deal.II/lac/vector.h>


Vector<double> result(100);


struct ScratchData
{};


struct CopyData
{
  unsigned int computed;
};


void worker (const std::vector<unsigned int>::iterator &i,
	     ScratchData &,
	     CopyData &ad)
{
  ad.computed = *i * 2;
}

void copier (const CopyData &ad)
{
  // write into the five elements of 'result' starting at ad.computed%result.size()
  for (unsigned int j=0; j<5; ++j)
    result((ad.computed+j) % result.size()) += ad.computed;
}


// the function that computes conflicts
std::vector<types::global_dof_index>
conflictor (const std::vector<unsigned int>::iterator &i)
{
  std::vector<types::global_dof_index> conflicts;
  const unsigned int ad_computed = *i * 2;
  for (unsigned int j=0; j<5; ++j)
    conflicts.push_back ((ad_computed+j) % result.size());

  return conflicts;
}



void test ()
{
  std::vector<unsigned int> v;
  for (unsigned int i=0; i<200; ++i)
    v.push_back (i);

  WorkStream::run (GraphColoring::make_graph_coloring (v.begin(), v.end(),
						       std_cxx11::function<std::vector<types::global_dof_index>
									   (const std::vector<unsigned int>::iterator &)>
						       (&conflictor)),
		   &worker, &copier,
                   ScratchData(),
                   CopyData());

  // now simulate what we should have gotten
  Vector<double> comp(result.size());
  for (unsigned int i=0; i<v.size(); ++i)
    {
      const unsigned int ad_computed = v[i] * 2;
      for (unsigned int j=0; j<5; ++j)
	comp((ad_computed+j) % result.size()) += ad_computed;
    }


  // and compare
  for (unsigned int i=0; i<result.size(); ++i)
    Assert (result(i) == comp(i), ExcInternalError());

  for (unsigned int i=0; i<result.size(); ++i)
    deallog << result(i) << std::endl;
}




int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}
