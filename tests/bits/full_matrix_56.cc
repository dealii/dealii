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



// check FullMatrix::cholesky for correct functionality

#include "../tests.h"
#include "full_matrix_common.h"


std::string output_file_name = "output";

const double entries2[4] =  { 0.814723686393179 ,  0.516389376684563,
                              0.516389376684563,   0.913375856139019
                            };


const double entries3[9] = {  1.808621732261680,   0.845866689167942,
                              0.667106522517665,   0.845866689167942,
                              1.398873167281503,   0.281706853672865,
                              0.667106522517665,   0.281706853672865,
                              0.741757766593798
                           };



// Create a positive definite random matrix

void random_matrix(FullMatrix<double> &A)
{
  for (unsigned int i=0; i<A.m(); ++i)
    for (unsigned int j=0; j<A.n(); ++j)
      {
        double rnd = Testing::rand();
        rnd /= RAND_MAX;
        A(i,j) = (i==j) ? A.m()+rnd : rnd;
      }
}

template <typename number>
void
check ()
{
  deallog << std::fixed;
  deallog << std::setprecision(4);
//  deallog.attach(logfile);
  deallog.depth_console(0);
//  deallog.threshold_double(1.e-10);
  Testing::srand(3391466);

  FullMatrix<double> T2(2,2,entries2), S2(2,2), R2(2,2);
  FullMatrix<double>T3(3,3,entries3), S3(3,3), R3(3,3);

  // These results for T2, T3, T5 were compared
  // against Matlab.
  deallog << "Original Matrix:" << std::endl;
  deallog << "================" << std::endl;

  display_matrix(T2);
  deallog << std::endl;

  S2.cholesky(T2);
  deallog << "Cholesky Decomposition Matrix:" << std::endl;
  deallog << "==============================" << std::endl;
  display_matrix(S2);
  deallog << std::endl;

  R2.Tadd(1,S2); // R = S^t
  deallog << "Transposed Cholesky Decomposition Matrix:" << std::endl;
  deallog << "=========================================" << std::endl;
  display_matrix(R2);
  deallog << std::endl;

  R2.Tmmult(S2, R2, false);  // S = R^tR
  deallog << "Reconstructed Matrix:" << std::endl;
  deallog << "=====================" << std::endl;
  display_matrix(S2);
  deallog << std::endl;
  S2.add(-1.0,T2);
  if (S2.frobenius_norm()>1.0e-10)
    {
      deallog << "NOT the same to 1e-10 tolerance" << std::endl;
    }
  else
    {
      deallog << "the SAME to 1e-10 tolerance" << std::endl;
    }
  deallog << std::endl;
  deallog << std::endl;

  deallog << "Original Matrix:" << std::endl;
  deallog << "================" << std::endl;
  display_matrix(T3);
  deallog << std::endl;

  S3.cholesky(T3);
  deallog << "Cholesky Decomposition Matrix:" << std::endl;
  deallog << "==============================" << std::endl;
  display_matrix(S3);
  deallog << std::endl;

  R3.Tadd(1,S3); // R = S^t
  deallog << "Transposed Cholesky Decomposition Matrix:" << std::endl;
  deallog << "=========================================" << std::endl;
  display_matrix(R3);
  deallog << std::endl;

  R3.Tmmult(S3, R3, false);  // S = R^tR
  deallog << "Reconstructed Matrix:" << std::endl;
  deallog << "=====================" << std::endl;
  display_matrix(S3);
  deallog << std::endl;
  S3.add(-1.0,T3);
  if (S3.frobenius_norm()>1.0e-10)
    {
      deallog << "NOT the same to 1e-10 tolerance" << std::endl;
    }
  else
    {
      deallog << "the SAME to 1e-10 tolerance" << std::endl;
    }
  deallog << std::endl;
  deallog << std::endl;

  FullMatrix<double> T5(5,5), S5(5,5), R5(5,5);
  random_matrix(T5);
  T5.symmetrize();

  deallog << "Original Matrix:" << std::endl;
  deallog << "================" << std::endl;
  display_matrix(T5);
  deallog << std::endl;

  S5.cholesky(T5);
  deallog << "Cholesky Decomposition Matrix:" << std::endl;
  deallog << "==============================" << std::endl;
  display_matrix(S5);
  deallog << std::endl;

  R5.Tadd(1,S5); // R = S^t
  deallog << "Transposed Cholesky Decomposition Matrix:" << std::endl;
  deallog << "=========================================" << std::endl;
  display_matrix(R5);
  deallog << std::endl;

  R5.Tmmult(S5, R5, false);  // S = R^tR
  deallog << "Reconstructed Matrix:" << std::endl;
  deallog << "=====================" << std::endl;
  display_matrix(S5);

  deallog << std::endl;
  S5.add(-1.0,T5);
  if (S5.frobenius_norm()>1.0e-10)
    {
      deallog << "NOT the same to 1e-10 tolerance" << std::endl;
    }
  else
    {
      deallog << "the SAME to 1e-10 tolerance" << std::endl;
    }
  deallog << std::endl;
  deallog << std::endl;

}



