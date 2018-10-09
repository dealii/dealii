// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2017 by the deal.II authors
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
// ---------------------------------------------------------------------



// check that SmartPointers preserve constness etc of the objects they
// point to, through assignment of SmartPointers to each other and
// other tests.


#include <deal.II/base/smartpointer.h>
#include <deal.II/base/subscriptor.h>

#include <iostream>
#include <vector>

#include "../tests.h"

// Provide memory for objects of type T such that access to a deleted
// object does not cause a segmentation fault
std::vector<char> memory(10000);
int               next = 0;

class Test : public Subscriptor
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


// A sentinel that stores a reference to a SmartPointer and resets the
// Smartpointer in its destructor.

template <typename T, typename P>
class Sentinel
{
public:
  Sentinel(SmartPointer<T, P> &smart_pointer)
    : smart_pointer_(smart_pointer)
  {}

  ~Sentinel()
  {
    // This assumes that the first object stored in SmartPointer is a raw
    // pointer T*
    *reinterpret_cast<T **>(&smart_pointer_) = nullptr;
  }

private:
  SmartPointer<T, P> &smart_pointer_;
};


int
main()
{
  deal_II_exceptions::disable_abort_on_exception();

  initlog();

  Test        a("A");
  const Test &b("B");

  SmartPointer<Test, Test>       r(&a, "Test R");
  SmartPointer<const Test, Test> s(&a, "const Test S");
  //  SmartPointer<Test,Test>       t=&b;    // this one should not work
  SmartPointer<Test, Test>       t(const_cast<Test *>(&b), "Test T");
  SmartPointer<const Test, Test> u(&b, "const Test");


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

  // We have to dance a happy little dance here:
  //
  // In our case the object D defined below will go out of scope before the
  // smartpointer "Test R" does. Under normal circumstances this triggers
  // an ExcNotUsed and aborts the program. BUT, this is a unit test and so,
  // aborting on exception is disabled and we continue execution.
  // Unfortunately, this triggers a "use after scope" memory access error
  // when finally "Test R" goes out of scope and tries to unsubscribe from
  // D that already got destroyed.
  //
  // Work around this issue by creating a sentinel that gets destroyed
  // after D but before "Test R" that simply resets the SmartPointer.
  Sentinel<Test, Test> sentinel(r);

  Test d("D");
  r = &d;
}
