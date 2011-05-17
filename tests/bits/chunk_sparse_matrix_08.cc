//----------------------------  chunk_sparse_matrix_08.cc  ---------------------------
//    chunk_sparse_matrix_08.cc,v 1.4 2004/02/26 17:25:44 wolf Exp
//    Version:  
//
//    Copyright (C) 2004, 2005, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  chunk_sparse_matrix_08.cc  ---------------------------


// check ChunkSparseMatrix::frobenius_norm

#include "../tests.h"
#include <deal.II/lac/chunk_sparse_matrix.h>    
#include <fstream>


void test (const unsigned int chunk_size)
{
  ChunkSparsityPattern sp (5,5,3,chunk_size);
  for (unsigned int i=0; i<5; ++i)
    for (unsigned int j=0; j<5; ++j)
      if ((i+2*j+1) % 3 == 0)
        sp.add (i,j);
  sp.compress ();

  ChunkSparseMatrix<double> m(sp);
  
                                   // first set a few entries. count how many
                                   // entries we have
  double norm = 0;
  for (unsigned int i=0; i<m.m(); ++i)
    for (unsigned int j=0; j<m.n(); ++j)
      if ((i+2*j+1) % 3 == 0)
        {
          m.set (i,j, i*j*.5+.5);
          norm += (i*j*.5+.5)*(i*j*.5+.5);
        }
  norm = std::sqrt(norm);

                                   // compare against the exact value of the
                                   // l2-norm (max row-sum)
  deallog << m.frobenius_norm() << std::endl;
  Assert (std::fabs((m.frobenius_norm() - norm)/norm) < 1e-14, ExcInternalError());
  
  deallog << "OK" << std::endl;
}



int main ()
{
  std::ofstream logfile("chunk_sparse_matrix_08/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      const unsigned int chunk_sizes[] = { 1, 2, 4, 5, 7 };
      for (unsigned int i=0;
	   i<sizeof(chunk_sizes)/sizeof(chunk_sizes[0]);
	   ++i)
	test (chunk_sizes[i]);
    }
  catch (std::exception &exc)
    {
      deallog << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      deallog << "Exception on processing: " << std::endl
		<< exc.what() << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      
      return 1;
    }
  catch (...) 
    {
      deallog << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      deallog << "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    };
}
