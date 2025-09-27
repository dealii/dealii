// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test what happens when serializing a pointer. is a new object created when
// loading into another pointer that is NULL? this is in fact what happens

#include <typeinfo>

#include "serialization.h"

int object_number = 1;

class C
{
public:
  C()
  {
    object_number = ::object_number++;
    deallog << "Default constructor. Object number " << object_number
            << std::endl;
  }

  C(const C &)
  {
    object_number = ::object_number++;
    deallog << "copy constructor. Object number " << object_number << std::endl;
  }

  template <typename Archive>
  void
  serialize(Archive &ar, const unsigned int version)
  {
    deallog << "Serializing object number " << object_number << " via "
            << typeid(Archive).name() << std::endl;
  }

  bool
  operator==(const C &) const
  {
    return true;
  }

private:
  unsigned int object_number;
};


void
test()
{
  C *p1 = new C();
  C *p2;

  verify(p1, p2);

  // p2 should have been created, and should
  // have been different from the address of
  // p1
  AssertThrow(p2 != nullptr, ExcInternalError());
  AssertThrow(p1 != p2, ExcInternalError());

  delete p1;
  delete p2;
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test();

  deallog << "OK" << std::endl;
}
