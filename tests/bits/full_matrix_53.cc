//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------


// check FullMatrix::left_invert and FullMatrix::right_invert


#include "../tests.h"
#include "full_matrix_common.h"


std::string output_file_name = "full_matrix_53/output";

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

