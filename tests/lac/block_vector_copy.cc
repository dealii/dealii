//--------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2008, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//--------------------------------------------------------------------


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/block_vector.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <list>

void test ()
{
  std::vector<double>		v(9);
  for (unsigned int i = 0; i < v.size(); ++i)
    v[i] = double(i+1);

  std::vector<unsigned int>	partition(3);
  for (unsigned int i = 0; i < partition.size(); ++i)
    partition[i] = 3;

  dealii::BlockVector<double>	b(partition);
  Assert (b.n_blocks() == partition.size(),
	  ExcInternalError());

  unsigned int			size = 0;
  for (unsigned int i = 0; i < b.n_blocks(); ++i)
    {
      Assert (b.block(i).size() == partition[i], ExcInternalError());
      size += b.block(i).size();
    }
  Assert (b.size() == size, ExcInternalError());

  for (unsigned int i = 0; i < b.size(); ++i)
    {
      b(i) = v[i];
      Assert (b(i) == v[i], ExcInternalError());
    }

  dealii::BlockVector<double>	c;
  c = b;
  Assert (c == b, ExcInternalError());
  Assert (c.n_blocks() == b.n_blocks(), ExcInternalError());

  deallog << "OK" << std::endl;
}




int main ()
{
  std::ofstream logfile("block_vector_copy/output");
  deallog << std::fixed;
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

                                   // do the same weird stuff as in
                                   // tests/base/reference.cc
#if __GNUC__ != 2
  std::basic_streambuf<char> *old_cerr_buf = std::cerr.rdbuf();
#else
  streambuf *old_cerr_buf = std::cerr.rdbuf();
#endif
  std::cerr.rdbuf(logfile.rdbuf());

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
      return 0;
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
      return 0;
    };

  std::cerr.rdbuf(old_cerr_buf);
}

