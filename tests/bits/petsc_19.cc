//----------------------------  petsc_19.cc  ---------------------------
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
//----------------------------  petsc_19.cc  ---------------------------


// check petsc_wrapper::Vector::linfty_norm()

#include "../tests.h"
#include <lac/petsc_vector.h>    
#include <fstream>
#include <iostream>
#include <vector>


void test (petsc_wrappers::Vector &v)
{
                                   // set some elements of the vector
  double norm = 0;
  for (unsigned int i=0; i<v.size(); i+=1+i)
    {
      v(i) = i;
      norm = std::max(norm,fabs(i));
    }
  v.compress ();

                                   // then check the norm
  Assert (v.linfty_norm() == norm, ExcInternalError());

  deallog << "OK" << std::endl;
}



int main (int argc,char **argv) 
{
  std::ofstream logfile("petsc_19.output");
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
