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


// Check LinearAlgebra::CUDAWrappers::Vector assignement and import

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/lac/cuda_vector.h>
#include <deal.II/lac/read_write_vector.h>
#include <fstream>
#include <iostream>

void test()
{
  const unsigned int size = 100;
  LinearAlgebra::CUDAWrappers::Vector<double> a;
  LinearAlgebra::CUDAWrappers::Vector<double> b(size);
  LinearAlgebra::CUDAWrappers::Vector<double> c(b);
  LinearAlgebra::CUDAWrappers::Vector<double> d;
  d.reinit(c);

  AssertThrow(a.size()==0, ExcMessage("Vector has the wrong size."));
  AssertThrow(b.size()==size, ExcMessage("Vector has the wrong size."));
  AssertThrow(c.size()==size, ExcMessage("Vector has the wrong size."));
  AssertThrow(d.size()==size, ExcMessage("Vector has the wrong size."));

  a.reinit(size);
  AssertThrow(a.size()==size, ExcMessage("Vector has the wrong size."));


  LinearAlgebra::ReadWriteVector<double> read_write_1(size);
  LinearAlgebra::ReadWriteVector<double> read_write_2(size);
  LinearAlgebra::ReadWriteVector<double> read_write_3(size);
  for (unsigned int i=0; i<size; ++i)
    {
      read_write_1[i] = i;
      read_write_2[i] = 5.+i;
    }

  a.import(read_write_2, VectorOperation::insert);
  b.import(read_write_1, VectorOperation::insert);
  c.import(read_write_2, VectorOperation::insert);


  read_write_3.import(a, VectorOperation::insert);
  for (unsigned int i=0; i<size; ++i)
    AssertThrow(read_write_2[i] == read_write_3[i],
                ExcMessage("Vector a has been modified."));

  read_write_3.import(b, VectorOperation::insert);
  for (unsigned int i=0; i<size; ++i)
    AssertThrow(read_write_1[i] == read_write_3[i],
                ExcMessage("Vector b has been modified."));

  read_write_3.import(c, VectorOperation::insert);
  for (unsigned int i=0; i<size; ++i)
    AssertThrow(read_write_2[i] == read_write_3[i],
                ExcMessage("Vector c has been modified."));

  a *= 2.;
  read_write_3.import(a, VectorOperation::insert);
  for (unsigned int i=0; i<size; ++i)
    AssertThrow(2.*read_write_2[i]==read_write_3[i],
                ExcMessage("Problem in operator *=."));

  c /= 2.;
  read_write_3.import(c, VectorOperation::insert);
  for (unsigned int i=0; i<size; ++i)
    AssertThrow(0.5*read_write_2[i]==read_write_3[i],
                ExcMessage("Problem in operator /=."));

  b += a;
  read_write_3.import(b, VectorOperation::insert);
  for (unsigned int i=0; i<size; ++i)
    AssertThrow(2.*read_write_2[i]+read_write_1[i]==read_write_3[i],
                ExcMessage("Problem in operator +=."));

  b -= c;
  read_write_3.import(b, VectorOperation::insert);
  for (unsigned int i=0; i<size; ++i)
    AssertThrow(1.5*read_write_2[i]+read_write_1[i]==read_write_3[i],
                ExcMessage("Problem in operator -=."));

  b.import(read_write_1, VectorOperation::insert);
  c.import(read_write_1, VectorOperation::insert);
  const double val = b*c;
  AssertThrow(val==328350., ExcMessage("Problem in operator *."));

  b = 0.;
  read_write_3.import(b, VectorOperation::insert);
  for (unsigned int i=0; i<size; ++i)
    AssertThrow(read_write_3[i] == 0.,ExcMessage("Problem in operator =."));
}

int main(int argc, char **argv)
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test();

  deallog << "OK" <<std::endl;

  return 0;
}
