// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
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



// test memory_consumption() on std_cxx11::array

#include "../tests.h"

#include <deal.II/base/memory_consumption.h>

#include <fstream>


int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);

  typedef std_cxx11::array<std_cxx11::array<int, 3>, 3> IntArray;
  deallog << dealii::MemoryConsumption::memory_consumption(IntArray()) << std::endl;
  deallog << sizeof(IntArray) << std::endl;
}
