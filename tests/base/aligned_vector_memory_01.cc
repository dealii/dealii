// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test AlignedVector. This class does some complex data management,
// explicitly allocating objects with placement-new expressions and
// explicitly releasing memory. Make sure this is all works as
// expected.

#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/logstream.h>

#include "../tests.h"


int object_number     = 0;
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

  C(const C &c)
  {
    object_number = ::object_number++;
    deallog << "copy constructor from " << c.object_number << ". Object number "
            << object_number << std::endl;
  }

  C(const C &&c)
  {
    object_number = ::object_number++;
    deallog << "move constructor from " << c.object_number << ". Object number "
            << object_number << std::endl;
  }

  ~C()
  {
    deallog << "destructor. Object number " << object_number << std::endl;
    ++objects_destroyed;
  }

  template <typename Archive>
  void
  serialize(Archive &ar, const unsigned int)
  {
    ar &object_number;
  }


private:
  unsigned int object_number;
};



void
test()
{
  {
    deallog << "Check 1:" << std::endl;
    AlignedVector<C> v(2);
  }

  {
    deallog << "Check 2:" << std::endl;
    AlignedVector<C> v(2);
    v.resize(0);
  }

  {
    deallog << "Check 3:" << std::endl;
    AlignedVector<C> v(2);
    v.resize(1);
    v.resize(3);
  }
}



int
main(int argc, char **argv)
{
  initlog();

  test();

  deallog << "Objects created: " << object_number << std::endl;
  deallog << "Objects destroyed: " << objects_destroyed << std::endl;
}
