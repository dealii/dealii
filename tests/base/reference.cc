// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1998 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check that ObserverPointers preserve constness etc of the objects they
// point to, through assignment of ObserverPointers to each other and
// other tests.


#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/observer_pointer.h>

#include <iostream>
#include <vector>

#include "../tests.h"

// Provide memory for objects of type T such that access to a deleted
// object does not cause a segmentation fault
std::vector<char> memory(10000);
int               next = 0;

class Test : public EnableObserverPointer
{
  const char *name;

public:
  Test(const char *n)
    : name(n)
  {
    deallog << "Construct " << name << std::endl;
  }
  ~Test()
  {
    deallog << "Destruct " << name << std::endl;
  }
  void
  f()
  {
    deallog << "mutable" << std::endl;
  }
  void
  f() const
  {
    deallog << "const" << std::endl;
  }
};



int
main()
{
  deal_II_exceptions::disable_abort_on_exception();

  initlog();

  Test        a("A");
  const Test &b("B");

  ObserverPointer<Test, Test>       r(&a, "Test R");
  ObserverPointer<const Test, Test> s(&a, "const Test S");
  //  ObserverPointer<Test,Test>       t=&b;    // this one should not work
  ObserverPointer<Test, Test>       t(const_cast<Test *>(&b), "Test T");
  ObserverPointer<const Test, Test> u(&b, "const Test");


  deallog << "a ";
  a.f(); // should print "mutable", since #a# is not const
  deallog << "b ";
  b.f(); // should print "const", since #b# is const
  deallog << "r ";
  r->f(); // should print "mutable", since it points to the non-const #a#
  deallog << "s ";
  s->f(); // should print "const", since it points to the const #b#
  // but we made it const
  deallog << "t ";
  t->f(); // should print "mutable", since #b# is const, but
  // we casted the constness away
  deallog << "u ";
  u->f(); // should print "const" since #b# is const
  // Now try if subscriptor works
  Test c("C");
  r = &c;

  // Test that a dangling pointer is correctly detected.
  try
    {
      {
        Test d("D");
        r = &d;
      }
      const auto dummy = *r;
    }
  catch (const ExceptionBase &exc)
    {
      deallog << exc.get_exc_name() << std::endl;
    }
}
