//----------------------------  petsc_50.cc  ---------------------------
//    petsc_11.cc,v 1.4 2003/07/03 10:31:46 guido Exp
//    Version: 
//
//    Copyright (C) 2004 by the deal.II authors and Anna Schneebeli
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  petsc_50.cc  ---------------------------


// check petsc_wrapper::copy_from<T> with T!=PetscScalar

#include "../tests.h"
#include <lac/petsc_vector.h>
#include <lac/vector.h>
#include <fstream>
#include <iostream>
#include <vector>


void test (petsc_wrappers::Vector &v)
{
  Vector<double> w (v.size());
  Vector<float>  x (v.size());

  for (unsigned int i=0; i<w.size(); ++i)
    {
      w(i) = i;
      x(i) = i+1;
    }
  
                                   // first copy from w and make sure we get
                                   // the expected result. then copy from x
                                   // and do the same. in at least one of the
                                   // two cases, the template argument to
                                   // Vector<T> must be different from
                                   // PetscScalar
  v.copy_from (w);
  for (unsigned int i=0; i<v.size(); ++i)
    {
      Assert (w(i) == i, ExcInternalError());
      Assert (v(i) == i, ExcInternalError());
      Assert (x(i) == i+1, ExcInternalError());
    }

  v.copy_from (x);
  for (unsigned int i=0; i<v.size(); ++i)
    {
      Assert (w(i) == i, ExcInternalError());
      Assert (v(i) == i+1, ExcInternalError());
      Assert (x(i) == i+1, ExcInternalError());
    }
  
  deallog << "OK" << std::endl;
}



int main (int argc,char **argv) 
{
  std::ofstream logfile("petsc_50.output");
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
