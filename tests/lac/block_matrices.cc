//----------------------------  block_matrices.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  block_matrices.cc  ---------------------------


#include <base/logstream.h>
#include <lac/block_sparsity_pattern.h>
#include <lac/block_sparse_matrix.h>
#include <lac/block_vector.h>
#include <fstream>
#include <algorithm>




void test () 
{
  ofstream logfile("block_matrices.output");
  logfile.setf(ios::fixed);
  logfile.precision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);

  BlockSparsityPattern bsp(3,2);
				   // set sizes
  bsp.block(0,0).reinit ( 2, 10, 7);
  bsp.block(0,1).reinit ( 2, 19, 8);
  bsp.block(1,0).reinit ( 7, 10, 6);
  bsp.block(1,1).reinit ( 7, 19, 7);
  bsp.block(2,0).reinit (10, 10, 6);
  bsp.block(2,1).reinit (10, 19, 12);
  bsp.collect_sizes ();

				   // add this pseudo-random sparsity
				   // pattern:
				   // 0   1    1 1 1    1 1 1    1 1 1 
				   // 1    1 1 1    1 1 1    1 1 1    1
				   // 2       1 1 1 1    1 1 1    1 1 1
				   // 3   1 1 1    1 1 1 1    1 1 1    
				   // 4   1    1 1 1    1 1 1 1    1 1 
				   // 5    1 1 1    1 1 1    1 1 1 1   
				   // 6    1    1 1 1    1 1 1    1 1 1
				   // 7   1 1 1 1    1 1 1    1 1 1    
				   // 8   1    1 1 1 1    1 1 1    1 1 
				   // 9   11 1 1    1 1 1 1    1 1 1   
				   // 1    1    1 1 1    1 1 1 1    1 1
				   // 11    1 1 1    1 1 1    1 1 1 1  
				   // 12  1 11   1 1 1    1 1 1    1 1 
				   // 13   1 111 1    1 1 1    1 1 1   
				   // 14   1   11 1 1 1    1 1 1    1 1
				   // 15    1 1 1    1 1 1 1    1 1 1  
				   // 16  1 1    1 1 1    1 1 1 1    1 
				   // 17     1 1 11   1 1 1    1 1 1 1 
				   // 18   1 1    111 1    1 1 1    1 1
  for (unsigned int row=0; row<19; ++row)
    for (unsigned int i=0; i<10; ++i)
      bsp.add (row, (row*5+i*9)%29);
  bsp.compress ();

				   // now check whether the elements
				   // we inserted are indeed those
				   // that are in there. for now, we
				   // only check their number, but
				   // their places are checked later
				   // with the matrix-vector
				   // operations.
  unsigned int total_nonzero_elements = 0;
  for (unsigned int row=0; row<19; ++row)
    {
				       // first count the number of
				       // elements in each row
      vector<bool> t(29, false);
      for (unsigned int i=0; i<10; ++i)
	t[(row*5+i*9)%29] = true;
				       // if we are in the third block
				       // row, then the first matrix
				       // is square, so there may be
				       // an additional element
      if (row>=9)
	t[row-9] = true;

      deallog << "Row " << row << " sparsity:  ";
      for (unsigned int i=0; i<29; ++i)
	deallog << t[i];
      deallog << endl;
      
      const unsigned int c=count(t.begin(), t.end(), true);
      
				       // now see how many elements
				       // there really are:
      unsigned int ac=0;
      for (unsigned int col=0; col<2; ++col)
	if (row<2)
	  ac += bsp.block(0,col).row_length(row-0);
	else
	  if (row<9)
	    ac += bsp.block(1,col).row_length(row-2);
	  else
	    ac += bsp.block(2,col).row_length(row-9);
      deallog << "Row=" << row
	      << ": expected length=" << c
	      << ", actual length=" << ac
	      << endl;
      total_nonzero_elements += ac;
      AssertThrow (c==ac, ExcInternalError());
    };
  deallog << total_nonzero_elements << "=="
	  << bsp.n_nonzero_elements()
	  << endl;
  AssertThrow (total_nonzero_elements == bsp.n_nonzero_elements(),
	       ExcInternalError());



  
				   // now make a matrix with this
				   // sparsity pattern
  BlockSparseMatrix<double,3,2> bsm (bsp);
  deallog << total_nonzero_elements << "=="
	  << bsm.n_nonzero_elements()
	  << endl;
  AssertThrow (total_nonzero_elements == bsm.n_nonzero_elements(),
	       ExcInternalError());
  
				   // try to write something into it,
				   // set entry (i,j) to i*j
  for (unsigned int row=0; row<19; ++row)
    for (unsigned int i=0; i<10; ++i)
      bsm.set (row, (row*5+i*9)%29, row*((row*5+i*9)%29));
				   // and add .5 to each value
  for (unsigned int row=0; row<19; ++row)
    for (unsigned int i=0; i<10; ++i)
      bsm.add (row, (row*5+i*9)%29, 0.5);

				   // now allocate two block vectors
				   // and see what we can get after
				   // vmults:
  BlockVector<double> src;
  vector<unsigned int> src_sizes (2);
  src_sizes[0] = 10;
  src_sizes[1] = 19;
  src.reinit (src_sizes);

  BlockVector<double> dst;
  vector<unsigned int> dst_sizes (3);
  dst_sizes[0] = 2;
  dst_sizes[1] = 7;
  dst_sizes[2] = 10;
  dst.reinit (dst_sizes);

  for (unsigned int i=0; i<29; ++i)
    src(i) = i;

  bsm.vmult (dst, src);
				   // now check what came out
  for (unsigned int row=0; row<19; ++row)
    {
      vector<double> t(29, 0);
				       // first check which elements
				       // in this row exist
      for (unsigned int i=0; i<10; ++i)
	t[(row*5+i*9)%29] = row*((row*5+i*9)%29);
      
      for (unsigned int i=0; i<10; ++i)
	t[(row*5+i*9)%29] += 0.5;

				       // compute the exact result
      double row_sum = 0;
      for (unsigned int i=0; i<29; ++i)
	row_sum += t[i]*i;

				       // compare to vmult result
      Assert (row_sum == dst(row), ExcInternalError());
      deallog << "vmult " << row << ' ' << row_sum << ' ' << dst(row) << endl;
    };


				   // test matrix_scalar_product. note that dst=M*src
  const double msp1 = dst.norm_sqr ();
  const double msp2 = bsm.matrix_scalar_product (dst, src);
  Assert (msp1 == msp2, ExcInternalError());
  deallog << "matrix_scalar_product " << msp1 << ' ' << msp2 << endl;
};




int main ()
{
  try
    {
      test ();
    }
  catch (exception &e)
    {
      cerr << endl << endl
	   << "----------------------------------------------------"
	   << endl;
      cerr << "Exception on processing: " << e.what() << endl
	   << "Aborting!" << endl
	   << "----------------------------------------------------"
	   << endl;
				       // abort
      return 2;
    }
  catch (...) 
    {
      cerr << endl << endl
	   << "----------------------------------------------------"
	   << endl;
      cerr << "Unknown exception!" << endl
	   << "Aborting!" << endl
	   << "----------------------------------------------------"
	   << endl;
				       // abort
      return 3;
    };
  
  
  return 0;
};
