//----------------------------  block_matrices.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  block_matrices.cc  ---------------------------


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/block_vector.h>

#include <fstream>
#include <iomanip>
#include <algorithm>
#include <iostream>




void test ()
{
  std::ofstream logfile("is_block_vector/output");
  deallog << std::fixed;
  deallog << std::setprecision(2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog << IsBlockVector<Vector<double> >::value << ' '
	  << IsBlockVector<Vector<float> >::value << ' '
	  << IsBlockVector<BlockVector<double> >::value << ' '
    	  << IsBlockVector<BlockVector<float> >::value << std::endl;
}




int main ()
{
  try
    {
      test ();
    }
  catch (std::exception &e)
    {
      std::cerr << std::endl << std::endl
	   << "----------------------------------------------------"
	   << std::endl;
      std::cerr << "Exception on processing: " << e.what() << std::endl
	   << "Aborting!" << std::endl
	   << "----------------------------------------------------"
	   << std::endl;
				       // abort
      return 2;
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
				       // abort
      return 3;
    };


  return 0;
}
