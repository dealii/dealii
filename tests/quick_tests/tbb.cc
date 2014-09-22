// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
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

// test Threads::new_task and WorkStream::run

#include <deal.II/base/thread_management.h>
#include <deal.II/base/work_stream.h>
#include <tbb/task_scheduler_init.h>
#include <iostream>

using namespace dealii;

void add_one(unsigned int &var)
{
  var += 1;
}

void test1()
{
  unsigned int tmp = 1;
  Threads::Task<> task = Threads::new_task (&add_one,tmp);
  task.join();
  if (tmp!=2)
    exit(1);
}

struct scratch_data
{
};

struct copy_data
{
    int value;
};


void assemble(const std::vector<int>::iterator &it,
	      scratch_data &scratch,
	      copy_data &data)
{
  data.value = (*it);
}

void copy(int & value, const copy_data &data)
{
  value += data.value;
}

void test2()
{
  const int maxi = 10000;
  std::vector<int> v(maxi);
  for (unsigned int i=0;i<v.size();++i)
    v[i] = i+1;
  int result = 0;
  WorkStream::run(v.begin(),
		  v.end(),
		  &assemble,
		  std_cxx11::bind(&copy,
				  std_cxx11::ref(result),
				  std_cxx11::_1),
		  scratch_data(), copy_data());
  std::cout << "result: " << result << std::endl;

  if (result != maxi*(maxi+1)/2)
    exit(2);
}

int main ()
{
  std::cout << "TBB will use " << tbb::task_scheduler_init::default_num_threads() << " threads." << std::endl;

  test1();
  test2();
}
