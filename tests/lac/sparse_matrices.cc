// ---------------------------------------------------------------------
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



//TODO: [GK] Produce some useful output!

#include "../tests.h"
#include "testmatrix.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/sparse_matrix_ez.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/precondition_block.h>

#include <deal.II/lac/sparse_matrix_ez.templates.h>
#include <deal.II/lac/precondition_block.templates.h>

#include <fstream>
#include <cstdio>


#define PREC_CHECK(solver, method, precond) \
  deallog.push(#precond); \
  try { u = 0.; solver.method (A, u, f, precond);  } catch (...) {} \
  deallog.pop(); \
  residuals.push_back(control.last_value())

template<class MATRIX>
void
check_vmult_quadratic(std::vector<double> &residuals,
                      const MATRIX &A,
                      const char *prefix)
{
  deallog.push(prefix);

  Vector<double> u(A.n());
  Vector<double> f(A.m());
  GrowingVectorMemory<> mem;

  SolverControl control(10, 1.e-13, false);
  SolverRichardson<> rich(control, mem, .01);
  SolverRichardson<> prich(control, mem, 1.);

  const types::global_dof_index block_size = (types::global_dof_index) std::sqrt(A.n()+.3);
  const unsigned int n_blocks = A.n()/block_size;

  typename PreconditionBlock<MATRIX, float>::AdditionalData
  data(block_size, 1.2);
  std::vector<types::global_dof_index> perm(A.n());
  std::vector<types::global_dof_index> iperm(A.n());
  for (unsigned int i=0; i<n_blocks; ++i)
    for (unsigned int j=0; j<block_size; ++j)
      {
        perm[block_size*i+j] = block_size*((i+1)%n_blocks)+j;
        iperm[perm[block_size*i+j]] = block_size*i+j;
      }

  PreconditionIdentity identity;
  PreconditionJacobi<MATRIX> jacobi;
  jacobi.initialize(A, .5);
  PreconditionSOR<MATRIX> sor;
  sor.initialize(A, 1.2);
//   PreconditionPSOR<MATRIX> psor;
//   psor.initialize(A, perm, iperm, 1.2);
  PreconditionSSOR<MATRIX> ssor;
  ssor.initialize(A, 1.2);

  PreconditionBlockJacobi<MATRIX, float> block_jacobi;
  block_jacobi.initialize(A, data);
  PreconditionBlockSSOR<MATRIX, float> block_ssor;
  block_ssor.initialize(A, data);
  PreconditionBlockSOR<MATRIX, float> block_sor;
  block_sor.initialize(A, data);
  PreconditionBlockSOR<MATRIX, float> block_psor;
  block_psor.set_permutation(perm, iperm);
  block_psor.initialize(A, data);

  f = 1.;

  PREC_CHECK(rich, solve, identity);
  PREC_CHECK(prich, solve, jacobi);
  PREC_CHECK(prich, solve, ssor);
  PREC_CHECK(prich, solve, sor);
//  PREC_CHECK(prich, solve, psor);
  PREC_CHECK(prich, solve, block_jacobi);
  PREC_CHECK(prich, solve, block_ssor);
  PREC_CHECK(prich, solve, block_sor);
  PREC_CHECK(prich, solve, block_psor);

  deallog << "Transpose" << std::endl;
  PREC_CHECK(rich, Tsolve, identity);
  PREC_CHECK(prich, Tsolve, jacobi);
  PREC_CHECK(prich, Tsolve, ssor);
  PREC_CHECK(prich, Tsolve, sor);
//  PREC_CHECK(prich, Tsolve, psor);
  PREC_CHECK(prich, Tsolve, block_jacobi);
  PREC_CHECK(prich, Tsolve, block_ssor);
  PREC_CHECK(prich, Tsolve, block_sor);
  PREC_CHECK(prich, Tsolve, block_psor);
  deallog.pop();
}



void
check_vmult_quadratic(std::vector<double> &residuals,
                      const BlockSparseMatrix<double> &A,
                      const char *prefix)
{
  deallog.push(prefix);

  Vector<double> u(A.n());
  Vector<double> f(A.m());
  GrowingVectorMemory<> mem;

  SolverControl control(10, 1.e-13, false);
  SolverRichardson<> rich(control, mem, .01);
  SolverRichardson<> prich(control, mem, 1.);
  PreconditionIdentity identity;
  PreconditionJacobi<BlockSparseMatrix<double> > jacobi;
  jacobi.initialize(A, .5);

  PreconditionBlock<BlockSparseMatrix<double>, float>::AdditionalData
  data((unsigned int) std::sqrt(A.n()+.3), 1.2);

  PreconditionBlockJacobi<BlockSparseMatrix<double>, float> block_jacobi;
  block_jacobi.initialize(A, data);

  u = 0.;
  f = 1.;

  PREC_CHECK(rich, solve, identity);
  PREC_CHECK(prich, solve, jacobi);
  u = 0.;
  PREC_CHECK(prich, solve, block_jacobi);

  u = 0.;
  deallog << "Transpose" << std::endl;
  PREC_CHECK(rich, Tsolve, identity);
  PREC_CHECK(prich, Tsolve, jacobi);
  u = 0.;
  PREC_CHECK(prich, Tsolve, block_jacobi);
  deallog.pop();
}


template <class MATRIX>
void
check_iterator (const MATRIX &A)
{
//  deallog.push("it");

  typename MATRIX::const_iterator E = A.end();

  if (A.m() < 10)
    for (unsigned int r=0; r<A.m(); ++r)
      {
        typename MATRIX::const_iterator b = A.begin(r);
        if (b == E)
          deallog << "Final" << std::endl;
        else
          deallog << r
                  << '\t' << b->row()
                  << '\t' << b->column()
                  << '\t' << b->value()
                  << std::endl;
        typename MATRIX::const_iterator e = A.end(r);
        if (e == E)
          deallog << "Final" << std::endl;
        else
          deallog << '\t' << e->row()
                  << std::endl;
        deallog << "cols:";

        for (typename MATRIX::const_iterator i=b; i!=e; ++i)
          deallog << '\t' << ',' << i->column();
        deallog << std::endl;
      }
  for (typename MATRIX::const_iterator i = A.begin(); i!= A.end(); ++i)
    deallog << '\t' << i->row()
            << '\t' << i->column()
            << '\t' << i->value()
            << std::endl;
  deallog << "Repeat row 2" << std::endl;
  for (typename MATRIX::const_iterator i = A.begin(2); i!= A.end(2); ++i)
    deallog << '\t' << i->row()
            << '\t' << i->column()
            << '\t' << i->value()
            << std::endl;

//  deallog.pop();
}


void check_ez_iterator()
{
  SparseMatrixEZ<float> m (6, 6, 0);

  deallog << "Empty matrix" << std::endl;

  check_iterator(m);

  deallog << "Irregular matrix" << std::endl;

  m.set (0, 0, 1.);
  m.set (1, 0, 2.);
  m.set (2, 0, 3.);
  m.set (2, 2, 0.);  // should be ignored
  m.set (0, 1, 1.);
  m.set (0, 2, 2.);
  m.set (0, 3, 3.);
  m.set (4, 3, 17.);

  check_iterator(m);
}


void check_conjugate(std::ostream &out)
{
  SparseMatrixEZ<double> B(3,2);
  SparseMatrixEZ<float>  A1(2,2);
  SparseMatrixEZ<float>  A2(3,3);
  SparseMatrixEZ<double> C1(3,3);
  SparseMatrixEZ<double> C2(2,2);

  B.set(0, 0, 0.);
  B.set(0, 1, 1.);
  B.set(1, 0, 2.);
  B.set(1, 1, 3.);
  B.set(2, 0, 4.);
  B.set(2, 1, 5.);

  A1.set(0, 0, 1.);
  A1.set(0, 1, 2.);
  A1.set(1, 0, 3.);
  A1.set(1, 1, 4.);

  A2.set(0, 0, 1.);
  A2.set(0, 1, 2.);
  A2.set(0, 2, 3.);
  A2.set(1, 0, 4.);
  A2.set(1, 1, 5.);
  A2.set(1, 2, 6.);
  A2.set(2, 0, 7.);
  A2.set(2, 1, 8.);
  A2.set(2, 2, 9.);

  C1.conjugate_add(A1,B);
  C2.conjugate_add(A2,B,true);

  out << "First conjugate" << std::endl;
  C1.print(out);

  out << "Second conjugate" << std::endl;
  C2.print(out);
}


int main()
{
  std::ofstream logfile("output");
  logfile << std::setprecision(3);
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console (0);

  // Switch between regression test
  // and benchmark
//#ifdef DEBUG
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  const unsigned int size = 5;
  const unsigned int row_length = 3;
/*#else
  deallog.depth_console(1000);
  deallog.log_execution_time(true);
  deallog.log_time_differences(true);
  const unsigned int size = 50;
  const unsigned int row_length = 9;
  #endif*/

  check_ez_iterator();
  check_conjugate(logfile);

  FDMatrix testproblem (size, size);
  unsigned int dim = (size-1)*(size-1);

  std::vector<double> A_res;
  std::vector<double> E_res;

  // usual sparse matrix
  SparsityPattern structure(dim, dim, 5);
  SparseMatrix<double>  A;
  {
    deallog << "Structure" << std::endl;
    testproblem.five_point_structure(structure);
    structure.compress();
    A.reinit (structure);
    deallog << "Assemble" << std::endl;
    testproblem.five_point(A, true);
    check_vmult_quadratic(A_res, A, "5-SparseMatrix<double>");
//#ifdef DEBUG
    check_iterator(A);
//#endif
  }

  // block sparse matrix with only
  // one block
  {
    deallog << "Structure" << std::endl;
    BlockSparsityPattern block_structure(1,1);
    block_structure.block(0,0).reinit (dim, dim, 5);
    block_structure.collect_sizes ();
    testproblem.five_point_structure(block_structure.block(0,0));
    block_structure.compress();
    BlockSparseMatrix<double>  block_A(block_structure);
    deallog << "Assemble" << std::endl;
    testproblem.five_point(block_A, true);
    std::vector<double> block_A_res;
    check_vmult_quadratic(block_A_res, block_A, "5-BlockSparseMatrix<double>");
//#ifdef DEBUG
    check_iterator(block_A);
//#endif
  }

  // ez sparse matrix
  SparseMatrixEZ<double> E(dim,dim,row_length,2);
  {
    deallog << "Assemble" << std::endl;
    testproblem.five_point(E, true);
    check_vmult_quadratic(E_res, E, "5-SparseMatrixEZ<double>");
//#ifdef DEBUG
    check_iterator(E);
//#endif
    E.print_statistics(deallog, true);
    E.add(-1., A);
    if (E.l2_norm() < 1.e-14)
      deallog << "Matrices are equal" << std::endl;
    else
      deallog << "Matrices differ!!" << std::endl;
  }

  {
    structure.reinit(dim, dim, 5);
    deallog << "Structure" << std::endl;
    structure.reinit(dim, dim, 9);
    testproblem.nine_point_structure(structure);
    structure.compress();
    A.reinit (structure);
    deallog << "Assemble" << std::endl;
    testproblem.nine_point(A);
    check_vmult_quadratic(A_res, A, "9-SparseMatrix<double>");
  }

  E.clear();
  E.reinit(dim,dim,row_length,2);
  deallog << "Assemble" << std::endl;
  testproblem.nine_point(E);
  check_vmult_quadratic(E_res, E, "9-SparseMatrixEZ<double>");
  E.print_statistics(deallog, true);

  for (unsigned int i=0; i<E_res.size(); ++i)
    if (std::fabs(A_res[i] - E_res[i]) > 1.e-13)
      deallog << "SparseMatrix and SparseMatrixEZ differ!!!"
              << std::endl;
  // dump A into a file, and re-read
  // it, then delete tmp file and
  // check equality
  std::ofstream tmp_write ("sparse_matrices.tmp");
  A.block_write (tmp_write);
  tmp_write.close ();

  std::ifstream tmp_read ("sparse_matrices.tmp");
  SparseMatrix<double> A_tmp;
  A_tmp.reinit (A.get_sparsity_pattern());
  A_tmp.block_read (tmp_read);
  tmp_read.close ();

  std::remove ("sparse_matrices.tmp");

  for (unsigned int i=0; i<A.n_nonzero_elements(); ++i)
    if (std::fabs(A.global_entry(i) - A_tmp.global_entry(i)) <=
        std::fabs(1e-14*A.global_entry(i)))
      deallog << "write/read-error at global position "
              << i << std::endl;
}
