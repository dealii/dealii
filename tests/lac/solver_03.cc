// ---------------------------------------------------------------------
// $Id$
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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



// test if number of blocks for src and dst can be
// different in GMRES 

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/block_matrix_array.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>


class Preconditioner 
{
  public:
  void vmult (BlockVector<double> &dst, 
      const BlockVector<double> &src) const;
};

void Preconditioner::vmult(BlockVector<double> &dst, 
    const BlockVector<double> &src) const
{
  //check that the sizes of input and
  //output vectors match, however we do not
  //check if the number of blocks is the same.
  Assert (dst.size() == src.size(),
      ExcDimensionMismatch (dst.size(), src.size()));

  for(unsigned int i=0; i<dst.size(); ++i)
    dst(i) = src(i);
}

template<class SOLVER>
void test()
{
  const unsigned int test = 2;
  const unsigned int trial = 3;
  SparsityPattern sparsity(test, trial, trial);
  SparsityPattern row_sp(1, trial, trial);
  for(unsigned int j=0; j<trial; ++j)
  {
    row_sp.add(0,j);
    for(unsigned int i=0; i<test; ++i)
      sparsity.add(i,j);
  }
  sparsity.compress();
  row_sp.compress();
  SparseMatrix<double> mat;
  mat.reinit(sparsity);
  SparseMatrix<double> row;
  row.reinit(row_sp);

  mat.set(0,0, 1);
  mat.set(0,1, 1);
  mat.set(0,2,  2);
  mat.set(1,0,  0);
  mat.set(1,1,  1);
  mat.set(1,2, -3);
  row.set(0,0,  0);
  row.set(0,1,  0);
  row.set(0,2,-19);

  BlockVector<double> rhs;
  rhs.reinit(2);
  rhs.block(0).reinit(test);
  rhs.block(1).reinit(1);
  rhs.collect_sizes();

  rhs(0) =  3;
  rhs(1) =- 4;
  rhs(2) =-19;

  BlockVector<double> solvec;
  solvec.reinit(1);
  solvec.block(0).reinit(trial);
  solvec.collect_sizes();

  BlockMatrixArray<double> bma;
  bma.reinit(2,1);
  bma.enter(mat, 0, 0);
  bma.enter(row, 1, 0);

  SolverGMRES<BlockVector<double> >::AdditionalData right_prec (30, true);
  SolverGMRES<BlockVector<double> >::AdditionalData left_prec (30, false);
  SolverControl solvctrl(1000, 1e-12, true);
  SOLVER solver_rightprec(solvctrl, right_prec);
  SOLVER solver_leftprec(solvctrl, left_prec);

  Preconditioner precond;
  solver_rightprec.solve(bma, solvec, rhs, precond);
  solvec.block(0).print(deallog);
  solvec = 0.;
  solver_leftprec.solve(bma, solvec, rhs, precond);
  solvec.block(0).print(deallog);
}

int main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<SolverGMRES<BlockVector<double> > >();
}

