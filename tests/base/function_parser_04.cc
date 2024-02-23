// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// functionparser: TBB

#include <deal.II/base/function_parser.h>
#include <deal.II/base/point.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/lac/vector.h>

#include <map>

#include "../tests.h"


FunctionParser<2> fp;

struct scratch_data
{};

struct copy_data
{
  int value;
  copy_data()
    : value(0)
  {}
};


void
assemble(const std::vector<int>::iterator &it,
         scratch_data                     &scratch,
         copy_data                        &data)
{
  double s     = *it;
  double value = fp.value(Point<2>(s, 2.5));
  Assert(std::abs(1.0 + s * 2.5 - value) < 1e-10, ExcMessage("wrong value"));
  std::cout << data.value << std::endl;

  data.value = (std::abs(1.0 + s * 2.5 - value) < 1e-10) ? 1 : 0;
}

void
copy(int &value, const copy_data &data)
{
  value += data.value;
}


void
test2()
{
  std::map<std::string, double> constants;
  constants["c"] = 1.0;
  fp.initialize("s,t", "s*t+c", constants);

  std::vector<int> v(10000);
  for (unsigned int i = 0; i < v.size(); ++i)
    v[i] = i;

  int result = 0;
  WorkStream::run(v.begin(),
                  v.end(),
                  &assemble,
                  std::bind(&copy, std::ref(result), std::placeholders::_1),
                  scratch_data(),
                  copy_data());
  std::cout << "result: " << result << std::endl;

  Assert((unsigned int)result == v.size(), ExcMessage("uhuh!"));
}


int
main()
{
  initlog();

  test2();
}
