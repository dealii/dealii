#include <base/logstream.h>


main()
{
  deallog << "Test" << endl;
  deallog.push("l1");
  deallog << "Test1" << endl;
  deallog.push("l2");
  deallog << "Test2" << "Test3" << endl;
  deallog.push("l3");
  deallog << "Test4";
  deallog.pop();
  deallog << "Test5" << endl;
  deallog.pop();
  deallog << "Test6" << endl;
  deallog.pop();
  deallog << "Test7" << endl;
}
