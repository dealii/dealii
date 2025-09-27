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
// same common object. the two objects here are the members of a std::pair

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


template <typename T>
bool
compare(const std::pair<T *, T *> &t1, const std::pair<T *, T *> &t2)
{
  return (*t1.first == *t2.first) && (*t1.second == *t2.second);
}


void
test()
{
  {
    C                  *p = new C();
    std::pair<C *, C *> pair_1(p, p);
    std::pair<C *, C *> pair_2;

    verify(pair_1, pair_2);

    // boost::serialize should have
    // recognized that the two pointers in
    // pair_1 point to the same object and
    // consequently re-create only one object
    // that the two components of the
    // re-created pair point to
    AssertThrow(pair_2.first == pair_2.second, ExcInternalError());
    AssertThrow(object_number == 3, ExcInternalError());

    // destroy the newly created object. this
    // must succeed
    delete pair_2.first;

    delete p;
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
