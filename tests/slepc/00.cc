//-------------------------------------------------------
//    $Id: 01.cc 24924 2013-01-28 young $
//    Version: $Name$ 
//
//    Copyright (C) 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-------------------------------------------------------


// just check initialising SLEPc can be done and that it initialises
// PETSc in the way we expect, ie. *a* PETSc object exist.

#include "../tests.h"
#include <deal.II/lac/slepc_solver.h>    
#include <deal.II/base/logstream.h>
#include <deal.II/base/numbers.h>
#include <fstream>
#include <iostream>

std::ofstream logfile ("00/output");

int main (int argc,char **argv) 
{
  deallog.attach (logfile);
  deallog.depth_console (1);

  try
    {

      logfile << "Initializing SLEPc (PETSc): " 
	      << std::flush;

      SlepcInitialize (&argc, &argv, 0, 0);
      {
	logfile << "ok" 
		<< std::endl;

	// Do something simple with PETSc
	logfile << "Using PetscScalar:" 
		<< std::endl;

	const PetscScalar pi  = numbers::PI;
	const PetscScalar two = 2.;

	logfile << "   pi:           " << pi 
		<< std::endl
		<< "   two:          " << two 
		<< std::endl
		<< "   two times pi: " << two*pi 
		<< std::endl;

	
	logfile << "Finalizing SLEPc (PETSc): " 
		<< std::flush;

      }
      SlepcFinalize ();

      logfile << "ok" 
	      << std::endl << std::endl;
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
