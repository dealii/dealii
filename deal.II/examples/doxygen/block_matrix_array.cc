//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

// See documentation of BlockMatrixArray for documentation of this example

#include <base/logstream.h>
#include <lac/block_matrix_array.h>
#include <lac/full_matrix.h>
#include <lac/vector.h>
#include <lac/block_vector.h>
#include <lac/vector_memory.h>
#include <lac/precondition.h>
#include <lac/solver_cg.h>
#include <lac/solver_gmres.h>

#include <iostream>
#include <fstream>

double Adata[] =
{
      4., .5, .1, 0.,
      .5, 4., .5, .1,
      .1, .5, 4., .5,
      0., .1, .5, 4.
};

double B1data[] =
{
      .5, .1,
      .4, .2,
      .3, .3,
      .2, .4
};

double B2data[] =
{
      .3, 0., -.3, 0.,
      -.3, 0., .3, 0.
};

double Cdata[] =
{
      8., 1.,
      1., 8.
};

int main () 
{
  FullMatrix<float> A(4,4);
  FullMatrix<float> B1(4,2);
  FullMatrix<float> B2(2,4);
  FullMatrix<float> C(2,2);

  A.fill(Adata);
  B1.fill(B1data);
  B2.fill(B2data);
  C.fill(Cdata);
  
  GrowingVectorMemory<Vector<double> > simple_mem;
  
  BlockMatrixArray<FullMatrix<float>, double> matrix(2, 2, simple_mem);
  
  matrix.enter(A,0,0,2.);
  matrix.enter(B1,0,1,-1.);
  matrix.enter(B2,0,1,1., true);
  matrix.enter(B2,1,0,1.);
  matrix.enter(B1,1,0,-1., true);
  matrix.enter(C,1,1);
  matrix.print_latex(deallog);
  
  std::vector<unsigned int> block_sizes(2);
  block_sizes[0] = 4;
  block_sizes[1] = 2;
  
  BlockVector<double> result(block_sizes);
  BlockVector<double> x(block_sizes);
  BlockVector<double> y(block_sizes);
  for (unsigned int i=0;i<result.size();++i)
    result(i) = i;

  matrix.vmult(y, result);

  SolverControl control(100,1.e-10);
  GrowingVectorMemory<BlockVector<double> > mem;
  PreconditionIdentity id;

  SolverCG<BlockVector<double> > cg(control, mem);
  cg.solve(matrix, x, y, id);
  x.add(-1., result);
  deallog << "Error " << x.l2_norm() << std::endl;
  
  deallog << "Error A-norm "
	  << std::sqrt(matrix.matrix_norm_square(x))
	  << std::endl;
  
  FullMatrix<float> Ainv(4,4);
  Ainv.invert(A);
  FullMatrix<float> Cinv(2,2);
  Cinv.invert(C);
  
  BlockTrianglePrecondition<FullMatrix<float>, double>
    precondition(2, simple_mem);
  precondition.enter(Ainv,0,0,.5);
  precondition.enter(Cinv,1,1);

  cg.solve(matrix, x, y, precondition);
  x.add(-1., result);
  deallog << "Error " << x.l2_norm() << std::endl;
  
  precondition.enter(B1,1,0,-1., true);
  precondition.enter(B2,1,0,1.);
  
  SolverGMRES<BlockVector<double> > gmres(control, mem);
  gmres.solve(matrix, x, y, precondition);
  x.add(-1., result);
  deallog << "Error " << x.l2_norm() << std::endl;
  
  return 0;
}
