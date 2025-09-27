// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test what happens when serializing two objects that have pointers to the
// same common object.
//
// similar to the _03 and _04 tests, but bury the pointer somewhere

#include <boost/serialization/utility.hpp>

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

  ~C()
  {
    deallog << "destructor. Object number " << object_number << std::endl;
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


struct P
{
  C *c;
  template <typename Archive>
  void
  serialize(Archive &ar, const unsigned int)
  {
    ar &c;
  }

  bool
  operator==(const P &p) const
  {
    // there is no useful operation we
    // can do here for the purposes of
    // serialization. do the actual tests
    // below with assertions
    return true;
  }
};



void
test()
{
  {
    C              *c = new C();
    std::pair<P, P> pair_1, pair_2;
    pair_1.first.c  = c;
    pair_1.second.c = c;

    pair_2.first.c  = nullptr;
    pair_2.second.c = nullptr;

    verify(pair_1, pair_2);

    // boost::serialize should have
    // recognized that the two pointers in
    // pair_1 point to the same object and
    // consequently re-create only one object
    // that the two components of the
    // re-created pair point to
    AssertThrow(pair_2.first.c == pair_2.second.c, ExcInternalError());
    AssertThrow(object_number == 3, ExcInternalError());

    // destroy the newly created object. this
    // must succeed and would likely throw
    // some sort of error if the object
    // pointed to was on the stack (like 'c')
    // rather than on the heap
    delete pair_2.first.c;
    delete c;
  }
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test();

  deallog << "OK" << std::endl;
}
