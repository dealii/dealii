// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check that IncrementalFunction resets the wrapped function time
// to its original state.

#include <deal.II/base/incremental_function.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"

// homogeneous linear in time function
template <int dim>
class MyFunc : public Function<dim>
{
public:
  MyFunc(const double a, const double b)
    : Function<dim>(dim)
    , a(a)
    , b(b)
  {}

  virtual double
  value(const Point<dim> /*p*/, const unsigned int component = 0) const
  {
    (void)component;
    return a * this->get_time() + b;
  }

  virtual void
  vector_value(const Point<dim> & /*p*/, Vector<double> &values) const override
  {
    for (unsigned int i = 0; i < dim; ++i)
      values[i] = a * this->get_time() + b;
  }

private:
  const double a;
  const double b;
};

template <int dim>
void
check()
{
  Vector<double> v1(dim);
  Vector<double> v2(dim);
  Vector<double> v3(dim);

  const Point<dim> point;

  MyFunc<dim>                         func(0.5, 10);
  MyFunc<dim>                         func2(0.5, 10);
  Functions::IncrementalFunction<dim> inc(func2);

  // Original time state:
  func2.set_time(7.0);
  deallog << "Wrapped func time (before): " << func2.get_time() << std::endl;

  std::vector<std::pair<double, double>> time_delta = {{1.0, 0.}, {1.0, 0.2}};

  for (auto &el : time_delta)
    {
      deallog << "time = " << el.first << " delta = " << el.second << std::endl;
      func.set_time(el.first);
      func.vector_value(point, v1);
      v1.print(deallog.get_file_stream());
      deallog << '-' << std::endl;

      func.set_time(el.first - el.second);
      func.vector_value(point, v2);
      v2.print(deallog.get_file_stream());
      deallog << '=' << std::endl;

      inc.set_time(el.first);
      inc.set_decrement(el.second);
      inc.vector_value(point, v3);
      v3.print(deallog.get_file_stream());

      v1.add(-1, v2);
      v1.add(-1, v3);
      AssertThrow(v1.l2_norm() == 0., ExcInternalError());
    }

  // Final time state (should be the same as the initial state:
  deallog << "Wrapped func time (after): " << func2.get_time() << std::endl;

  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();

  check<2>();
}
