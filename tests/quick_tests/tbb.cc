// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// test Threads::new_task and WorkStream::run

#include <deal.II/base/thread_management.h>
#include <deal.II/base/work_stream.h>

#include <iostream>

using namespace dealii;

void
add_one(unsigned int &var)
{
  var += 1;
}

void
test1()
{
  unsigned int    tmp  = 1;
  Threads::Task<> task = Threads::new_task(&add_one, tmp);
  task.join();
  if (tmp != 2)
    exit(1);
}

struct scratch_data
{};

struct copy_data
{
  int value;
};


void
assemble(const std::vector<int>::iterator &it,
         scratch_data & /*scratch*/,
         copy_data &data)
{
  data.value = (*it);
}

void
copy(int &value, const copy_data &data)
{
  value += data.value;
}

void
test2()
{
  const int        maxi = 10000;
  std::vector<int> v(maxi);
  for (unsigned int i = 0; i < v.size(); ++i)
    v[i] = i + 1;
  int result = 0;
  WorkStream::run(
    v.begin(),
    v.end(),
    &assemble,
    [&result](const copy_data &data) { copy(result, data); },
    scratch_data(),
    copy_data());
  std::cout << "result: " << result << std::endl;

  if (result != maxi * (maxi + 1) / 2)
    exit(2);
}

int
main()
{
  std::cout << "TBB will use " << MultithreadInfo::n_threads() << " threads."
            << std::endl;

  test1();
  test2();
}
