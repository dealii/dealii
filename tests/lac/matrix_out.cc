//----------------------------  matrix_out.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  matrix_out.cc  ---------------------------


#include <base/logstream.h>
#include <lac/matrix_out.h>
#include <lac/full_matrix.h>
#include <fstream>

int main () 
{
  std::ofstream logfile("matrix_out.output");
  logfile.setf(std::ios::fixed);
  logfile.precision(2);
  deallog.attach(logfile);
  deallog.depth_console(0);

				   // test for a square full matrix
  if (true)
    {
      FullMatrix<double> full_matrix(4,4);
      for (unsigned int i=0; i<4; ++i)
	full_matrix(i,i) = 1;
      
      MatrixOut matrix_out;
      matrix_out.build_patches (full_matrix, "full_matrix");
      matrix_out.write_gnuplot (logfile);
    };
  
				   // test for a rectangular sparse
				   // matrix
  if (true)
    {
      SparsityPattern sparsity (4,8,7);
      for (unsigned int i=0; i<4; ++i)
	for (unsigned int j=0; j<8; ++j)
	  if (i!=j)
	    sparsity.add (i,j);
      sparsity.compress ();

      SparseMatrix<double> sparse_matrix(sparsity);
      for (unsigned int i=0; i<4; ++i)
	for (unsigned int j=0; j<8; ++j)
	  sparse_matrix.set(i,j, static_cast<signed int>(i-j));
  
      MatrixOut matrix_out;
      matrix_out.build_patches (sparse_matrix, "sparse_matrix",
				MatrixOut::Options (true));
      matrix_out.write_eps (logfile);
    };

				   // test collation of elements
  if (true)
    {
      FullMatrix<double> full_matrix(20,20);
      for (unsigned int i=0; i<20; ++i)
	for (unsigned int j=0; j<20; ++j)
	  full_matrix(i,j) = (1.*i*i/20/20-1.*j*j*j/20/20/20);
  
      MatrixOut matrix_out;
      matrix_out.build_patches (full_matrix, "collated_matrix",
				MatrixOut::Options (false, 4));
      matrix_out.write_gmv (logfile);
    };
};
