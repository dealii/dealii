//----------------------------  petsc_35.cc  ---------------------------
//    petsc_11.cc,v 1.4 2003/07/03 10:31:46 guido Exp
//    Version: 
//
//    Copyright (C) 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  petsc_35.cc  ---------------------------


// check petsc_wrapper::Vector::operator+=(Vector) 

#include "../tests.h"
#include <lac/petsc_vector.h>    
#include <fstream>
#include <iostream>
#include <vector>


void test (petsc_wrappers::Vector &v,
           petsc_wrappers::Vector &w)
{
                                   // set only certain elements of each
                                   // vector
  for (unsigned int i=0; i<v.size(); ++i)
    {
      v(i) = i;
      if (i%3 == 0)
        w(i) = i+1.;
    }
  
  v.compress ();
  w.compress ();

  v += w;

                                   // make sure we get the expected result
  for (unsigned int i=0; i<v.size(); ++i)
    {
      if (i%3 == 0)
        {
          Assert (w(i) == i+1., ExcInternalError());
          Assert (v(i) == i+i+1., ExcInternalError());
        }
      else
        {
          Assert (w(i) == 0, ExcInternalError());
          Assert (v(i) == i, ExcInternalError());
        }
    }


  deallog << "OK" << std::endl;
}



int main (int argc,char **argv) 
{
  std::ofstream logfile("petsc_35.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  try
    {
      PetscInitialize(&argc,&argv,0,0);
      {
        petsc_wrappers::Vector v (100);
        petsc_wrappers::Vector w (100);
        test (v,w);
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
