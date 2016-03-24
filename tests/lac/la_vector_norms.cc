// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2015 by the deal.II authors
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


// checks that the l1 and l2 norm are computed correctly for deal.II vectors
// (check the summation algorithm), including an accuracy test (should not
// lose more than 1 decimal also for 200000 vector entries)

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/la_vector.h>
#include <cmath>
#include <fstream>
#include <iomanip>




template <typename number>
void check_norms ()
{
  const number acc = 1e1*std::numeric_limits<number>::epsilon();
  unsigned int skip = 73;
  for (unsigned int size=1; size<200000; size+=skip)
    {
      // test correctness
      if (size > 10000)
        skip += 17;
      LinearAlgebra::Vector<number> vec(size);
      for (unsigned int i=0; i<size; ++i)
        vec[i] = i+1;

      const number l1_norm = vec.l1_norm();
      AssertThrow (std::abs(l1_norm-0.5*size*(size+1)) < acc*0.5*size*(size+1),
                   ExcInternalError());

      // test accuracy of summation
      const long double value = 3.14159265358979323846;
      for (unsigned int i=0; i<size; ++i)
        vec[i] = (number)value;
      const number l1_norma = vec.l1_norm();
      AssertThrow (std::abs(l1_norma-value*size) < acc*size*value,
                   ExcInternalError());
      const number l2_norma = vec.l2_norm();
      AssertThrow (std::abs(l2_norma-value*std::sqrt((number)size)) < acc*std::sqrt(size)*value,
                   ExcInternalError());
    }
}


template <typename number>
void check_complex_norms ()
{
  const number acc = 1e2*std::numeric_limits<number>::epsilon();
  unsigned int skip = 73;
  for (unsigned int size=1; size<100000; size+=skip)
    {
      // test correctness
      if (size > 10000)
        skip += 17;
      LinearAlgebra::Vector<std::complex<number> > vec(size);
      long double sum = 0.;
      for (unsigned int i=0; i<size; ++i)
        {
          vec(i) = std::complex<number>(i+1,i+2);
          sum += std::sqrt((long double)(i+1)*(i+1) + (long double)(i+2)*(i+2));
        }

      const number l1_norm = vec.l1_norm();
      AssertThrow (std::abs(l1_norm-sum) < acc*sum,
                   ExcInternalError());

      // test accuracy of summation
      const std::complex<long double> value (3.14159265358979323846, 0.1);
      for (unsigned int i=0; i<size; ++i)
        vec[i] = std::complex<number>(value);
      const number l1_norma = vec.l1_norm();
      AssertThrow (std::abs(l1_norma-std::abs(value)*size) <
                   acc*size*std::abs(value),
                   ExcInternalError());
      const number l2_norma = vec.l2_norm();
      AssertThrow (std::abs(l2_norma-std::abs(value)*std::sqrt((number)size)) <
                   acc*std::sqrt((number)size)*std::abs(value),
                   ExcInternalError());
    }
}


int main()
{
  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(2);
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  check_norms<float>();
  check_norms<double>();
  check_norms<long double>();
  check_complex_norms<double>();
  check_complex_norms<float>();
  deallog << "OK" << std::endl;
}


