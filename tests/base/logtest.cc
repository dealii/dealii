//----------------------------  logtest.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  logtest.cc  ---------------------------


#include <base/logstream.h>
#include <fstream>

int main()
{
  ofstream logfile("logtest.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  deallog << "Test" << std::endl;
  deallog.push("l1");
  deallog << "Test1" << std::endl;
  deallog.push("l2");
  deallog << "Test2" << "Test3" << std::endl;
  deallog.push("l3");
  deallog << "Test4";
  deallog.pop();
  deallog << "Test5" << std::endl;
  deallog.pop();
  deallog << "Test6" << std::endl;
  deallog.pop();
  deallog << "Test7" << std::endl;
}
