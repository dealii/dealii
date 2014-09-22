// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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



// functionparser: TBB

#include "../tests.h"
#include <fstream>
#include <iomanip>
#include <map>
#include <deal.II/base/logstream.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/work_stream.h>


FunctionParser<2> fp;

struct scratch_data
{
};

struct copy_data
{
    int value;
    copy_data():
		    value(0)
      {}
    
};


void assemble(const std::vector<int>::iterator &it,
              scratch_data &scratch,
              copy_data &data)
{
  double s = *it;
  double value = fp.value(Point<2>(s, 2.5));
  Assert(abs(1.0+s*2.5 - value) < 1e-10, ExcMessage("wrong value"));
  std::cout << data.value  << std::endl;
  
  data.value = (abs(1.0+s*2.5 - value) < 1e-10)?1:0;
}

void copy(int & value, const copy_data &data)
{
  value += data.value;
}


void test2()
{
  std::map<std::string, double> constants;
  constants["c"]=1.0;
  fp.initialize("s,t", "s*t+c", constants);

  std::vector<int> v(10000);
  for (unsigned int i=0;i<v.size();++i)
    v[i] = i;

  int result = 0;
  WorkStream::run(v.begin(),
                  v.end(),
                  &assemble,
                  std_cxx11::bind(&copy,
                                  std_cxx11::ref(result),
                                  std_cxx11::_1),
                  scratch_data(), copy_data());
  std::cout << "result: " << result << std::endl;

  Assert(result == v.size(), ExcMessage("uhuh!"));
}


int main ()
{
  initlog();

  test2();
}




