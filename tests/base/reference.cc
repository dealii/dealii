
#include <iostream>
#include <base/subscriptor.h>
#include <base/smartpointer.h>


class Test : public Subscriptor
{
  const char* name;
public:
  Test(const char* n) :
		  name(n)
      {
	cout << "Construct " << name << endl;
      }
  ~Test()
      {
	cout << "Destruct " << name << endl;
      }	  
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
  Test a("A");
  const Test b("B");
  SmartPointer<Test>       r=&a;
  SmartPointer<const Test> s=&a;
//  SmartPointer<Test>       t=&b;    // this one should not work
  SmartPointer<Test>       t=const_cast<Test*>(&b);
  SmartPointer<const Test> u=&b;

  
  cout << "a ";
  a.f();            // should print "mutable", since #a# is not const
  cout << "b ";
  b.f();            // should print "const", since #b# is const
  cout << "r ";
  r->f();           // should print "mutable", since it points to the non-const #a#
  cout << "s ";
  s->f();           // should print "const", since it points to the const #b#
				   // but we made it const
  cout << "t ";
  t->f();           // should print "mutable", since #b# is const, but
				   // we casted the constness away
  cout << "u ";
  u->f();           // should print "const" since #b# is const
				   // Now try if subscriptor works
  {
    Test c("C");
    r = &c;
    Test d("D");
    r = &d;
  }
}

void abort()
{}
