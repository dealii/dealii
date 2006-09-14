//----------------------------------------------------------------------
//    $Id: reference.cc 13395 2006-07-19 12:45:54Z kanschat $
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

// A test program that produces correct results


#include <cstdlib>
#include <fstream>

int main()
{
  std::ofstream out("ok_01/output");
  out << "My output" << std::endl;
}
