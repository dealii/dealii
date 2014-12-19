// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2013 by the deal.II authors
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



// check FullMatrix::left_invert and FullMatrix::right_invert


#include "../tests.h"
#include "full_matrix_common.h"


std::string output_file_name = "output";

template <typename number>
void
calculate(const FullMatrix<number> A,
          const FullMatrix<number> B)
{

  FullMatrix<number> A_r_inv(A.n(),A.m());
  FullMatrix<number> identity(A.m(),A.m());
  A_r_inv.right_invert(A);
  deallog<<"A matrix"<<std::endl;
  display_matrix(A);
  deallog<<"Right inverse"<<std::endl;
  display_matrix(A_r_inv);
  deallog<<"Identity = A*A_r_inv"<<std::endl;
  A.mmult(identity,A_r_inv);
  display_matrix(identity);

  deallog<<std::endl;

  FullMatrix<number> B_l_inv(B.n(),B.m());
  FullMatrix<number> identity2(B.n(),B.n());
  B_l_inv.left_invert(B);
  deallog<<"B matrix"<<std::endl;
  display_matrix(B);
  deallog<<"Left inverse"<<std::endl;
  display_matrix(B_l_inv);
  deallog<<"Identity = B_l_inv*B"<<std::endl;
  B_l_inv.mmult(identity2,B);
  display_matrix(identity2);

}


template <typename number>
void
check ()
{

  FullMatrix<number> A(1,2);
  fill_matrix(A);

  FullMatrix<number> B(2,1);
  fill_matrix(B);

  calculate(A,B);

  FullMatrix<number> A1(2,3);
  fill_matrix(A1);

  FullMatrix<number> B1(3,2);
  fill_matrix(B1);

  calculate(A1,B1);

}

