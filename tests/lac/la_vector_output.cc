// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2017 by the deal.II authors
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


#include <deal.II/lac/la_vector.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <string>

#include "../tests.h"


void
test()
{
  const char *                  filename = "test.txt";
  const unsigned int            size(10);
  LinearAlgebra::Vector<double> vec(size);
  for (unsigned int i = 0; i < size; ++i)
    vec[i] = i;

  // Check block_write and block_read
  std::ofstream file_out(filename);
  // Write the vector in the file
  vec.block_write(file_out);
  file_out.close();
  // Clear the vector
  vec.reinit(0);
  // Load the vector form the file
  std::ifstream file_in(filename);
  vec.block_read(file_in);
  file_in.close();
  // Check the vector
  double eps = 1e-12;
  for (unsigned int i = 0; i < size; ++i)
    AssertThrow(
      std::abs(vec[i] - (double)i) < eps,
      ExcMessage(
        "Value in the vector has been changed by block_write or block_read"));


  // save data to archive
  {
    std::ofstream                 file_out2(filename);
    boost::archive::text_oarchive oa(file_out2);
    oa << vec;
    // archive and stream closed when destructors are called
  }

  // Clear the vector
  vec.reinit(0);
  {
    std::ifstream                 file_in2(filename);
    boost::archive::text_iarchive ia(file_in2);
    ia >> vec;
  }
  // Check the vector
  for (unsigned int i = 0; i < size; ++i)
    AssertThrow(
      std::abs(vec[i] - (double)i) < eps,
      ExcMessage("Value in the vector has been changed by boost archive"));
}


int
main()
{
  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(2);
  deallog.attach(logfile);

  test();

  deallog << "OK" << std::endl;
}
