// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
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



#include "../tests.h"
#include <deal.II/algorithms/any_data.h>
#include <deal.II/base/logstream.h>

#include <fstream>

double d1 = 17.;

void fill(AnyData& data)
{
  int i = 7;
  data.add(i, " i  7");
  data.add<double>(d1, " d  17.");
  data.add<double*>(&d1, " d* 17.");
  data.add<const double*>(&d1, "cd* 17.");
  d1 = 18.;
}


void extract(const AnyData& data)
{
  // This set of tests with old functionality. Remove when deprecating
  // index access
  for (unsigned int i=0;i<data.size();++i)
    deallog << i << '\t' << data.name(i) << std::endl;
  deallog << data.name(0)
	  << '\t' << data.entry<int>(0) << std::endl;
  deallog << data.name(1)
	  << '\t' << data.entry<double>(1) << std::endl;
  deallog << data.name(2)
	  << '\t' << *data.entry<double*>(2) << std::endl;
  deallog << data.name(3)
	  << '\t' << *data.entry<const double*>(3) << std::endl;

  // From here on keep after deprecation
  int i1 = data.entry<int>(" i  7");
  const int i2 = data.entry<const int>(" i  7");
  double d = data.entry<double>(" d  17.");
  double* p2 = data.entry<double*>(" d* 17.");
  const double * p3 = data.entry<const double*>("cd* 17.");
  
  deallog << i1 << std::endl
	  << i2 << std::endl
	  << d  << std::endl
	  << *p2 << std::endl
	  << *p3 << std::endl;

  deallog << *data.try_read<double>(" d  17.") << std::endl
	  << data.try_read<char *>(" d  17.") << std::endl
	  << data.try_read<double>("does not exist") << std::endl;
// try
  //   {
  //     double* p3a = data.entry<double*>("cd* 17.");
  //     deallog << p3a;
  //   }
  // // catch(ExceptionBase e)
  // //   {
  // //     deallog << e.what() << std::endl;
  // //   }
  // catch(...)
  //   {
  //     deallog << "Exception duly thrown" << std::endl;
  //   }

}

int main()
{
  deal_II_exceptions::disable_abort_on_exception();
  
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  
  AnyData data;
  fill(data);
  extract(data);
}
