// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check the functor interface


#include <deal.II/base/index_set.h>

#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/read_write_vector.templates.h>

#include "../tests.h"


struct Functor
{
  void
  operator()(double &value) const
  {
    value *= 2.;
  }
};

void
test()
{
  const unsigned int                     size = 25;
  LinearAlgebra::ReadWriteVector<double> vector(size);
  for (unsigned int i = 0; i < size; ++i)
    vector[i] = i;

  Functor functor;
  vector.apply(functor);
  for (unsigned int i = 0; i < size; ++i)
    deallog << vector[i] << std::endl;
}

int
main()
{
  initlog();
  test();

  return 0;
}
