// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// the versions of BlockMatrixBase::set/add that take a pointer to a set of
// values were not thread safe, even if they were writing into disjoint sets
// of elements. this test verifies that the set() function is now indeed
// thread safe under these conditions
//
// we test this by calling add() from a multitude of threads at once. without
// the patch that fixed this, one can elicit a great deal of different memory
// corruption errors and assertions from this test, depending on luck and
// phase of the moon :-)

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <fstream>
#include <iomanip>
#include <algorithm>



void do_set (const bool even_or_odd,
	     BlockSparseMatrix<double> &bsm)
{
  BlockSparseMatrix<double>::size_type col_indices[5];
  for (unsigned int i=0; i<5 ; ++i)
    if (even_or_odd)
      col_indices[i] = 2*i;
    else
      col_indices[i] = 2*i+1;

  BlockSparseMatrix<double>::value_type values[5];
  for (unsigned int i=0; i<5 ; ++i)
    if (even_or_odd)
      values[i] = 1;
    else
      values[i] = 2;

  bsm.set (0, 5, col_indices, values, false);
}


void test ()
{
  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  BlockSparsityPattern bsp(2,2);
  // set sizes
  for (unsigned int i=0; i<2; ++i)
    for (unsigned int j=0; j<2; ++j)
      bsp.block(i,j).reinit ( 5, 5, 5);
  bsp.collect_sizes ();

  // make a full matrix
  for (unsigned int row=0; row<10; ++row)
    for (unsigned int i=0; i<10; ++i)
      bsp.add (row, i);
  bsp.compress ();

  BlockSparseMatrix<double> bsm (bsp);

  Threads::ThreadGroup<> tg;
  for (unsigned int i=0; i<100; ++i)
    {
      tg += Threads::new_thread (&do_set, true, bsm);
      tg += Threads::new_thread (&do_set, false, bsm);
    }
  tg.join_all ();

  bsm.print_formatted (deallog.get_file_stream());
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
