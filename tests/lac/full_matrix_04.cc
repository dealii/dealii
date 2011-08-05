//----------------------------  full_matrix_03.cc,v  ---------------------------
//    $Id$
//
//    Copyright (C) 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  full_matrix_03.cc,v  ---------------------------

//check method FullMatrix::triple_prduct

#include "../tests.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iomanip>
#include <cstdlib>

#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/eigen.h>

void diff(FullMatrix<double>& M)
{
  const double err = M.frobenius_norm();
  if (err < 1.e-14)
    deallog << "ok" << std::endl;
  else
    deallog << "oops " << err << std::endl;
}


void test (const unsigned int n, const unsigned int m)
{
				   // Create some random matrices
  FullMatrix<double> A(n,n);
  FullMatrix<double> C(m,m);  
  FullMatrix<double> B(n,m);
  FullMatrix<double> D(m,n);
  FullMatrix<double> Bt(m,n);
  FullMatrix<double> Dt(n,m);

  for (unsigned int i=0;i<n;++i)
    {
      for (unsigned int j=0;j<n;++j)
	A(i,j) = std::rand();
      for (unsigned int j=0;j<m;++j)
	{
	  B(i,j) = Bt(j,i) = std::rand();
	  D(j,i) = Dt(i,j) = std::rand();
	}
    }
  for (unsigned int i=0;i<m;++i)
    for (unsigned int j=0;j<m;++j)
      C(i,j) = std::rand();

				   // Compare first Schur complement
				   // with mmult.
  FullMatrix<double> S1(m,m);
  S1.triple_product(A, D, B, false, false);

  FullMatrix<double> aux1(n,m);
  A.mmult(aux1,B);
  FullMatrix<double> aux2(m,m);
  D.mmult(aux2,aux1);
  aux2.add(-1., S1);
  diff(aux2);

				   // Compare second Schur complement
				   // with mmult
  FullMatrix<double> S2(n,n);
  S2.triple_product(C, B, D, false, false);

  FullMatrix<double> aux3(m,n);
  C.mmult(aux3,D);
  FullMatrix<double> aux4(n,n);
  B.mmult(aux4,aux3);
  aux4.add(-1., S2);
  diff(aux4);

				   // Check transpose versions
  aux2 = 0.;
  aux2.triple_product(A, Dt, B, true, false);
  aux2.add(-1., S1);
  diff(aux2);
  aux2 = 0.;
  aux2.triple_product(A, D, Bt, false, true);
  aux2.add(-1., S1);
  diff(aux2);
  aux2 = 0.;
  aux2.triple_product(A, Dt, Bt, true, true);
  aux2.add(-1., S1);
  diff(aux2);
  
  aux4 = 0.;
  aux4.triple_product(C, Bt, D, true, false);
  aux4.add(-1., S2);
  diff(aux4);
  aux4 = 0.;
  aux4.triple_product(C, B, Dt, false, true);
  aux4.add(-1., S2);
  diff(aux4);
  aux4 = 0.;
  aux4.triple_product(C, Bt, Dt, true, true);
  aux4.add(-1., S2);
  diff(aux4);
}


int
main ()
{
  const std::string logname = JobIdentifier::base_name(__FILE__) + std::string("/output");
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console(0);
  std::srand(3391466);

  test(3,4);
  test(4,7);
}    
