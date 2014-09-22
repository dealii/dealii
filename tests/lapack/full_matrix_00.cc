// ---------------------------------------------------------------------
//
// Copyright (C) 2014 by the deal.II authors
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


// Tests reinitialisation of square and rectangle LAPACKFullMatrix

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/lapack_full_matrix.h>

#include <fstream>
#include <iostream>


void test (const unsigned int size,
	   const bool         reinit_square)
{
  // this test can not currently work with matrices smaller than
  // 1\times2.
  Assert (size>2, ExcInternalError());

  // initialise a first matrix with the standard constructor and fill
  // it with some numbers
  LAPACKFullMatrix<double> M (size, size);
  
  for (unsigned int i=0; i<size; ++i)
    for (unsigned int j=0; j<size; ++j)
      M(i,j) = i+2.*j;

  // initialise a second matrix with the standard constructor and fill
  // it with some numbers
  LAPACKFullMatrix<double> N (size+2, size-2);

  for (unsigned int i=0; i<N.m(); ++i)
    for (unsigned int j=0; j<N.n (); ++j)
      N(i,j) = i+2.*j;
  
  // clearly, this should be the case
  Assert (N.m () != M.m (), ExcInternalError());
  Assert (N.n () != M.n (), ExcInternalError());

  // if reinit_square is true, reinitialise the rectangle matrix to a
  // square matrix (use reinit (const unsigned int))
  if (reinit_square)
    {
      // reinitialise the matrix and fill it with some numbers
      N.reinit (size);
      
      for (unsigned int i=0; i<N.m (); ++i)
	for (unsigned int j=0; j<N.n (); ++j)
	  N(i,j) = i+2.*j;
    }

  // otherwise reinitialise the rectangle matrix to a square one (use
  // reinit (const unsigned int, const unsigned int))
  else
    {
      // reinitialise the matrix and fill it with some numbers
      M.reinit (size+2, size-2);
      
      for (unsigned int i=0; i<M.m (); ++i)
	for (unsigned int j=0; j<M.n (); ++j)
	  M(i,j) = i+2.*j;
    }

  // and now this should be true
  Assert (N.m () == M.m (), ExcInternalError());
  Assert (N.n () == M.n (), ExcInternalError());

  // in fact, this should be true too, so check 
  for (unsigned int i=0; i<M.m (); ++i)
    for (unsigned int j=0; j<M.n (); ++j)
      Assert (M(i,j) == N(i,j), ExcInternalError());

  deallog << "OK" << std::endl;
}


int main()
{
  const std::string logname = "output";
  std::ofstream logfile(logname.c_str());

  logfile.precision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  // Test square matrix initialisation
  test (4, true);
  test (5, true);
  test (6, true);

  // Test rectangle matrix initialisation
  test (4, false);
  test (5, false);
  test (6, false);


}
