//--------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//--------------------------------------------------------------------


#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "testmatrix.h"
#include <base/logstream.h>
#include <lac/sparse_matrix.h>
#include <lac/vector.h>
#include <lac/vector_memory.h>
#include <lac/solver_control.h>
#include <lac/eigen.h>
#include <lac/precondition.h>


int main()
{
    std::ofstream logfile("eigen.output");
//  logfile.setf(std::ios::fixed);
  logfile.precision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);

  GrowingVectorMemory<> mem;
  SolverControl control(1000, 1.e-5);
  
  const unsigned int size = 10;
  const unsigned int dim = (size-1)*(size-1);
  
  FDMatrix testproblem(size, size);
  SparsityPattern structure(dim, dim, 5);
  testproblem.five_point_structure(structure);
  structure.compress();
  SparseMatrix<double>  A(structure);
  testproblem.five_point(A);

  Vector<double> u(dim);
  u = 1.;
  double lambda;

  EigenPower<> mises(control, mem);
  mises.solve(lambda, A, u);
  deallog << "Eigenvalue " << lambda << std::endl;

  double lambda2;
  u = 1.;
  
  EigenPower<> mises2(control, mem, -1.5*lambda);
  mises2.solve(lambda2, A, u);
  deallog << "Eigenvalue " << lambda2 << std::endl;

  u = 1.;
  lambda = 0.;
  EigenInverse<> wieland(control, mem);
  wieland.solve(lambda, A, u);
  deallog << "Eigenvalue " << lambda << std::endl;  

  u = 1.;
  lambda = 10.;
  EigenInverse<> wieland2(control, mem);
  wieland2.solve(lambda, A, u);
  deallog << "Eigenvalue " << lambda << std::endl;  
}
