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
// loading into another pointer that is non-NULL and the old pointer
// destroyed? Or is the old object pointed to being co-opted? the former is in
// fact what happens, and to make things just ever so slightly more awkward,
// the previous object pointed to isn't freed but is left dangling, creating
// the potential for a memory leak

#include <typeinfo>

#include "serialization.h"

int object_number     = 1;
int objects_destroyed = 0;

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

  ~C()
  {
    deallog << "destructor. Object number " << object_number << std::endl;
    ++objects_destroyed;
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
  C *backup;
  {
    C *p1  = new C();
    C *p2  = new C(); // this is the pointer that will be overwritten
    backup = p2;      // but save a pointer to the original object

    verify(p1, p2);

    AssertThrow(p1 != p2, ExcInternalError());

    delete p1;
    delete p2;
  }

  // as mentioned above, p2 is overwritten by
  // a pointer to a new object, leaving the
  // original object pointed to as a memory
  // leak. assert that this behavior persists
  AssertThrow(objects_destroyed == 2, ExcInternalError());

  // we've checked what we wanted to check, so now delete
  // the original object to ensure we don't get undue
  // error messages from memory checkers about leaked memory
  delete backup;
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test();

  deallog << "OK" << std::endl;
}
