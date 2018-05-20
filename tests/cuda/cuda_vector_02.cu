// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2016 by the deal.II authors
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
#include <deal.II/base/utilities.h>
#include <deal.II/lac/cuda_vector.h>
#include <deal.II/lac/read_write_vector.h>
#include <fstream>
#include <iostream>
#include <vector>

// Check LinearAlgebra::CUDAWrappers::Vector add and sadd.

void
test()
{
  const unsigned int                          size = 100;
  LinearAlgebra::CUDAWrappers::Vector<double> a(size);
  LinearAlgebra::CUDAWrappers::Vector<double> b(size);
  LinearAlgebra::CUDAWrappers::Vector<double> c(size);

  LinearAlgebra::ReadWriteVector<double> read_write_1(size);
  LinearAlgebra::ReadWriteVector<double> read_write_2(size);
  LinearAlgebra::ReadWriteVector<double> read_write_3(size);

  for(unsigned int i = 0; i < size; ++i)
    {
      read_write_1[i] = i;
      read_write_2[i] = 5. + i;
    }

  a.import(read_write_1, VectorOperation::insert);
  b.import(read_write_2, VectorOperation::insert);
  c.import(read_write_2, VectorOperation::insert);

  a.add(1.);
  read_write_3.import(a, VectorOperation::insert);
  for(unsigned int i = 0; i < size; ++i)
    AssertThrow(1. + read_write_1[i] == read_write_3[i],
                ExcMessage("Problem in add(scalar)."));

  a.add(2., b);
  read_write_3.import(a, VectorOperation::insert);
  for(unsigned int i = 0; i < size; ++i)
    AssertThrow(1. + read_write_1[i] + 2. * read_write_2[i] == read_write_3[i],
                ExcMessage("Problem in add(scalar,Vector)."));

  LinearAlgebra::CUDAWrappers::Vector<double> d(a);
  a.add(2., b, 3., d);
  read_write_3.import(a, VectorOperation::insert);
  for(unsigned int i = 0; i < size; ++i)
    AssertThrow(4. + 4. * read_write_1[i] + 10. * read_write_2[i]
                  == read_write_3[i],
                ExcMessage("Problem in add(scalar,Vector,scalar,Vector)."));

  a.import(read_write_1, VectorOperation::insert);
  a.sadd(3., 2., c);
  read_write_3.import(a, VectorOperation::insert);
  for(unsigned int i = 0; i < size; ++i)
    AssertThrow(3. * read_write_1[i] + 2. * read_write_2[i] == read_write_3[i],
                ExcMessage("Problem in sadd(scalar,scalar,Vector)."));

  a.import(read_write_1, VectorOperation::insert);
  a.scale(b);
  read_write_3.import(a, VectorOperation::insert);
  for(unsigned int i = 0; i < size; ++i)
    AssertThrow(read_write_1[i] * read_write_2[i] == read_write_3[i],
                ExcMessage("Problem in scale."));

  a.equ(2., c);
  read_write_3.import(a, VectorOperation::insert);
  for(unsigned int i = 0; i < size; ++i)
    AssertThrow(2. * read_write_2[i] == read_write_3[i],
                ExcMessage("Problem in equ."));

  AssertThrow(b.mean_value() == 54.50, ExcMessage("Problem in mean_value."));

  AssertThrow(b.l1_norm() == 5450., ExcMessage("Problem in l1_norm."));

  const double eps = 1e-3;
  AssertThrow(std::fabs(b.l2_norm() - 616.725222) < eps,
              ExcMessage("Problem in l2_norm"));

  AssertThrow(b.linfty_norm() == 104., ExcMessage("Problem in linfty_norm."));

  a.import(read_write_1, VectorOperation::insert);
  const double val = a.add_and_dot(2., a, b);
  AssertThrow(val == 1059300., ExcMessage("Problem in add_and_dot"));
}

int
main(int argc, char** argv)
{
  initlog();
  deallog.depth_console(0);

  test();

  deallog << "OK" << std::endl;

  return 0;
}
