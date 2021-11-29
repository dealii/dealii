// --------------------------------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// --------------------------------------------------------------------------------------------


// Check for a bug in memory management in class Table (which uses
// AlignedVector) that had crept in at some point.

#include <deal.II/base/logstream.h>
#include <deal.II/base/table.h>

#include "../tests.h"


int object_number     = 0;
int objects_destroyed = 0;

class C
{
public:
  C()
  {
    object_number = ::object_number++;
    deallog << "default constructor. Object number " << object_number
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

  C &
  operator=(const C &c)
  {
    deallog << "copy operator called for " << object_number << " <- "
            << c.object_number << std::endl;
    return *this;
  }

  C &
  operator=(const C &&c)
  {
    deallog << "move operator called for " << object_number << " <- std::move("
            << c.object_number << ')' << std::endl;
    return *this;
  }


  ~C()
  {
    deallog << "destructor. Object number " << object_number << std::endl;
    ++objects_destroyed;
  }

private:
  unsigned int object_number;
};



void
test()
{
  deallog << "---- Creating outer table" << std::endl;
  Table<1, C> table({1});

  // Copy the object, then destroy the copy again.
  {
    deallog << "---- Cloning outer table" << std::endl;
    Table<1, C> x(table);
    deallog << "---- Destroying the clone" << std::endl;
  }

  deallog << "---- Destroying the source table" << std::endl;
}



int
main(int argc, char **argv)
{
  initlog();

  test();

  deallog << "Objects created: " << object_number << std::endl;
  deallog << "Objects destroyed: " << objects_destroyed << std::endl;
}
