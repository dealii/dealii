//----------------------------  petsc_16.cc  ---------------------------
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
//----------------------------  petsc_16.cc  ---------------------------


// check petsc_wrapper::Vector::operator() in set/add-mode alternatingly, but
// writing and overwriting the same elements

#include "../tests.h"
#include <lac/petsc_vector.h>    
#include <fstream>
#include <iostream>
#include <vector>


void test (petsc_wrappers::Vector &v)
{
                                   // set only certain elements of the
                                   // vector. have a bit pattern of where we
                                   // actually wrote elements to
  std::vector<bool> pattern (v.size(), false);
  for (unsigned int i=0; i<v.size(); i+=1+i)
    {
      v(i) = i;
      pattern[i] = true;
    }
  for (unsigned int i=0; i<v.size(); i+=1+i)
    v(i) += i;
  for (unsigned int i=0; i<v.size(); i+=1+i)
    v(i) = i;

  v.compress ();

                                   // check that they are ok, and this time
                                   // all of them
  for (unsigned int i=0; i<v.size(); ++i)
    Assert (((pattern[i] == true) && (v(i) == i)
             ||
             (pattern[i] == false) && (v(i) == 0)),
             ExcInternalError());

  deallog << "OK" << std::endl;
}



int main (int argc,char **argv) 
{
  std::ofstream logfile("petsc_16.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  try
    {
      PetscInitialize(&argc,&argv,0,0);
      {
        petsc_wrappers::Vector v (100);
        test (v);
      }
      PetscFinalize();
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
