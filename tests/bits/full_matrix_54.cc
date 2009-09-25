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


// another check FullMatrix::Tadd
// proper transposition of rectangular matrices is verified


#include "../tests.h"
#include "full_matrix_common.h"


std::string output_file_name = "full_matrix_54/output";


template <typename number>
void
check ()
{
  FullMatrix<number> A(8,9);
  fill_matrix(A);
  deallog<<"Given matrix"<<std::endl;
  display_matrix(A);
  FullMatrix<number> A_t(A.n(),A.m());
  A_t.Tadd(A,1);
  deallog<<"Transposed matrix"<<std::endl;
  display_matrix(A_t);
}
