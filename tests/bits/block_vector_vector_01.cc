//----------------------------  block_vector_vector_01.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  block_vector_vector_01.cc  ---------------------------


// check existence of
// BlockVector<double>::BlockVector(BlockVector<float>). this conversion
// constructor was disabled previously altogether because of a compiler defect
// that did not honor the 'explicit' keyword on template constructors. this is
// now autoconf'ed.

#include "../tests.h"
#include <lac/vector.h>    
#include <fstream>
#include <iostream>


void test (Vector<double> &v)
{
  for (unsigned int i=0; i<v.size(); ++i)
    v(i) = i+1.;
  Vector<float> w(v);

  Assert (w==v, ExcInternalError());
  
  deallog << "OK" << std::endl;
}



int main () 
{
  std::ofstream logfile("block_vector_vector_01.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  try
    {
      Vector<double> v (100);
      test (v);
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Exception on processing: " << std::endl
		<< exc.what() << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      
      return 1;
    }
  catch (...) 
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    };
}
