//----------------------------  get_fe_from_name.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005 by the deal.II authors 
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  get_fe_from_name.cc  ---------------------------
// Test the get_fe_from_name function of the fe_tools. 

#include "../tests.h"
#include <fstream>
#include <base/logstream.h>
#include <fe/fe.h>
#include <fe/fe_tools.h>

int main () 
{ 
  std::ofstream logfile("get_fe_from_name.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  
  FiniteElement<1> * fe1;
  FiniteElement<2> * fe2;
  FiniteElement<3> * fe3;
  
  std::string name;
  name = "FE_Q(1)";
  fe1 = FETools::get_fe_from_name<1>(name);
  fe2 = FETools::get_fe_from_name<2>(name);
  fe3 = FETools::get_fe_from_name<3>(name);
  
  deallog << "Read " << name << std::endl;
  deallog << "Generated :" << std::endl;
  deallog << fe1->get_name()  << std::endl;
  deallog << fe2->get_name()  << std::endl;
  deallog << fe3->get_name()  << std::endl;
  
  delete fe1;
  delete fe2;
  delete fe3;
  
  return 0;
}
    
