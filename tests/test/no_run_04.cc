//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------

// A test program that won't run, but produces results


#include <cstdlib>
#include <fstream>

#include "phony.h"

int main()
{
  std::ofstream out("no_run_04/output");
  out << "My output" << std::endl;
  
  exit(1);
}
