//----------------------------  $RCSfile$  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  $RCSfile$  ---------------------------


#include <cmath>
#include <fstream>
#include <iomanip>

#include <base/logstream.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/eigen.h>

int
main ()
{
  FullMatrix<double> A(5,5), C(5,5), D(5,5), H(5,5);
  D(0,0) = 1.;
  D(1,1) = 2.;
  D(2,2) = 3.;
  D(3,3) = 4.;
  D(4,4) = 5.;

  A = D;
  
  for (unsigned int i=0;i<4;++i)
    {
				       // Setup rotation matrix
      C.clear();
      C.diagadd(1.);
      C(i,i) = C(i+1,i+1) = cos(i+1);
      C(i+1,i) = sin(i+1);
      C(i,i+1) = -sin(i+1);

				       // Rotate original matrix
      A.mmult(H,C);
      C.Tmmult(A,H);
    }

  A.print_formatted (cout);

  Vector<double> u(5);
  GrowingVectorMemory<Vector<double> > mem;
  
  SolverControl control (500,1.e-8, false, true);
  
  {
    u = 1.;
    EigenPower<Vector<double> >
      von_Mises(control, mem, 0.);
    double eigen = 0.;
    von_Mises.solve(eigen, A, u);
    deallog << "Eigenvalue: " << eigen << endl;
  }
  {
    u = 1.;
    EigenPower<Vector<double> >
      von_Mises(control, mem, -4.);
    double eigen = 0.;
    von_Mises.solve(eigen, A, u);
    deallog << "Eigenvalue: " << eigen << endl;
  }
}

      
