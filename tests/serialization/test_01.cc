// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2013 by the deal.II authors
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


// a basic test for some serialization functionality

#include "serialization.h"

#include <typeinfo>

int object_number = 1;

class C
{
public:
  C ()
  {
    object_number = ::object_number++;
    deallog << "Default constructor. Object number "
            << object_number
            << std::endl;
  }

  C (const C &)
  {
    object_number = ::object_number++;
    deallog << "copy constructor. Object number "
            << object_number
            << std::endl;
  }

  template <typename Archive>
  void serialize (Archive &ar, const unsigned int version)
  {
    deallog << "Serializing object number "
            << object_number
            << " via " << typeid(Archive).name()
            << std::endl;
  }

  bool operator == (const C &) const
  {
    return true;
  }

private:
  unsigned int object_number;
};


void test ()
{
  C p1, p2;

  verify (p1, p2);
}


int main ()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();

  deallog << "OK" << std::endl;
}
