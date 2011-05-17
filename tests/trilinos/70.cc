//----------------------------  trilinos_70.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2004, 2005, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  trilinos_70.cc  ---------------------------


// check TrilinosScalar

#include "../tests.h" 
#include <deal.II/base/utilities.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/trilinos_vector.h>

#include <fstream>

int main (int argc, char** argv) 
{
  std::ofstream logfile("70/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10); 

  Utilities::System::MPI_InitFinalize mpi_initialization (argc, argv);


  if (typeid(TrilinosScalar)==typeid(double))
    deallog << "double" << std::endl;
  else if (typeid(TrilinosScalar)==typeid(float))
    deallog << "float" << std::endl;
  else
    Assert(false, ExcNotImplemented());
}
