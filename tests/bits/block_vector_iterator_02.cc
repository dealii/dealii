//----------------------------  block_vector_iterator_02.cc  ---------------------------
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
//----------------------------  block_vector_iterator_02.cc  ---------------------------


// like _01, except that we use operator[] instead of operator*

#include "../tests.h"
#include <lac/block_vector.h>
#include <fstream>
#include <iostream>


void test ()
{
  BlockVector<double> v(2,1);
  v(0) = 1;
  v(1) = 2;

                                   // first check reading through a const
                                   // iterator
  {
    BlockVector<double>::const_iterator i=v.begin();
    Assert (i[0] == 1, ExcInternalError());
    Assert (i[1] == 2, ExcInternalError());
  }

                                   // same, but create iterator in a different
                                   // way
  {
    BlockVector<double>::const_iterator
      i=const_cast<const BlockVector<double>&>(v).begin();
    Assert (i[0] == 1, ExcInternalError());
    Assert (i[1] == 2, ExcInternalError());
  }

                                   // read through a read-write iterator
  {
    BlockVector<double>::iterator i = v.begin();
    Assert (i[0] == 1, ExcInternalError());
    Assert (i[1] == 2, ExcInternalError());
  }

                                   // write through a read-write iterator
  {
    BlockVector<double>::iterator i = v.begin();
    i[0] = 2;
    i[1] = 3;
  }  

                                   // and read again
  {
    BlockVector<double>::iterator i = v.begin();
    Assert (i[0] == 2, ExcInternalError());
    Assert (i[1] == 3, ExcInternalError());
  }
  
  deallog << "OK" << std::endl;
}



int main ()
{
  std::ofstream logfile("block_vector_iterator_02.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  try
    {
      test ();
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
