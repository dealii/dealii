
#include <iostream>
#include <base/subscriptor.h>



class Test : public Subscriptor
{
public:
  void f()
  {
    cout << "mutable" << endl;
  }
  void f() const
  {
    cout << "const" << endl;
  }
};



main()
{
  Test a;
  const Test b;
  SmartPointer<Test>       r=&a;
  SmartPointer<const Test> s=&a;
  SmartPointer<Test>       t=&b;    // this one should give a warning
  SmartPointer<const Test> u=&b;

  a.f();            // should print "mutable", since #a# is not const
  b.f();            // should print "const", since #b# is const
  r->f();           // should print "mutable", since it points to the non-const #a#
  s->f();           // should print "const", since it points to the non-const #a#
				   // but we made it const
  t->f();           // should print "mutable", since #b# is const, but
				   // we casted the constness away
  u->f();           // should print "const" since #b# is const
}
