//----------------------------  reference.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  reference.cc  ---------------------------


// check that SmartPointers preserve constness etc of the objects they
// point to, through assignment of SmartPointers to each other and
// other tests.


#include "../tests.h"
#include <fstream>
#include <base/subscriptor.h>
#include <base/smartpointer.h>
#include <base/logstream.h>
#include <iostream>


class Test : public Subscriptor
{
    const char* name;
  public:
    Test(const char* n) : name(n) { deallog << "Construct " << name << std::endl; }
    ~Test()                       { deallog << "Destruct " << name << std::endl;  }
    void f()        { deallog << "mutable" << std::endl; }
    void f() const  { deallog << "const" << std::endl;   }
};


int main()
{
  std::ofstream logfile("reference.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

                                   // we do something rather weird in
                                   // this file: bind the buffer of
                                   // cerr to the log file, so that
                                   // the output of the Assert* macros
                                   // that are written to std::cerr
                                   // end up in the logfiles.
                                   //
                                   // so make sure we store a pointer
                                   // to the old buffer, switch to the
                                   // buffer of the logfile, and at
                                   // the end of main() switch
                                   // back. note that if we don't
                                   // switch back, we get a segfault
                                   // later on in the destruction of
                                   // std::cerr, since it tries to do
                                   // something with its buffer, but
                                   // that is already gone by then
#if __GNUC__ != 2
  std::basic_streambuf<char> *old_cerr_buf = std::cerr.rdbuf();
#else
  streambuf *old_cerr_buf = std::cerr.rdbuf();
#endif
  std::cerr.rdbuf(logfile.rdbuf());

  if (true)
    {
      Test a("A");
      const Test b("B");
      SmartPointer<Test>       r(&a, "Test R");
      SmartPointer<const Test> s(&a, "const Test S");
//  SmartPointer<Test>       t=&b;    // this one should not work
      SmartPointer<Test>       t(const_cast<Test*>(&b), "Test T");
      SmartPointer<const Test> u(&b, "const Test");
      
      
      deallog << "a ";
      a.f();            // should print "mutable", since #a# is not const
      deallog << "b ";
      b.f();            // should print "const", since #b# is const
      deallog << "r ";
      r->f();           // should print "mutable", since it points to the non-const #a#
      deallog << "s ";
      s->f();           // should print "const", since it points to the const #b#
				       // but we made it const
      deallog << "t ";
      t->f();           // should print "mutable", since #b# is const, but
				       // we casted the constness away
      deallog << "u ";
      u->f();           // should print "const" since #b# is const
				       // Now try if subscriptor works
      Test c("C");
      r = &c;
      Test d("D");
      r = &d;
    }
  std::cerr.rdbuf(old_cerr_buf);
}

