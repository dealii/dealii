
#include <iostream>
#include <base/exceptions.h>
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
  SmartPointer<Test> r=&a;
  const SmartPointer<Test> s=&a;
  SmartPointer<Test> t=&b;
  const SmartPointer<Test> u=&b;

  a.f();
  b.f();
  r->f();
  s->f();
  t->f();
  u->f();
}
