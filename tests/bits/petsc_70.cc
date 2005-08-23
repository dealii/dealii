//----------------------------  petsc_70.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  petsc_70.cc  ---------------------------


// check PetscScalar

#include "../tests.h"
#include <base/logstream.h>
#include <lac/petsc_vector.h>

#include <fstream>

int main (int argc,char **argv) 
{
  std::ofstream logfile("petsc_70.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  if (typeid(PetscScalar)==typeid(double))
    deallog << "double" << endl;
  else if (typeid(PetscScalar)==typeid(float))
    deallog << "float" << endl;
  else
    Assert(false, ExcNotImplemented());
}
