// ---------------------------------------------------------------------
//
// Copyright (C) 2014 by the deal.II authors
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

// check SchurMatrix::prepare_rhs

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_matrix_array.h>
#include <list>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/schur_matrix.h>
#include <vector>

template<class MA_inverse, class MB, class MDt, class MC>
void
checkPrepare_rhs(BlockVector<double> &dst, const BlockVector<double> &src,
                 const MA_inverse &Ainv, const MB &B, const MDt &Dt, const MC &C)
{
  deallog << "prepare_rhs" << std::endl;

  GrowingVectorMemory < BlockVector<double> > mem;

  SchurMatrix<MA_inverse, MB, MDt, MC> S(Ainv, B, Dt, C, mem);

  S.debug_level(1);

  S.prepare_rhs(dst, src);

}

int
main()
{
  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  std::vector < types::global_dof_index > ivector(2);
  ivector[0] = 1;
  ivector[1] = 1;

  double src_array[] =
  { 9, 21 };
  std::list<double> src_l(&src_array[0], &src_array[2]);
  BlockVector<double> src(ivector, src_l.begin(), src_l.end());

  std::vector < types::global_dof_index > ivector2(1);
  ivector2[0] = 1;

  double dst_array[] =
  { 2 };
  std::list<double> dst_l(&dst_array[0], &dst_array[1]);
  BlockVector<double> dst(ivector2, dst_l.begin(), dst_l.end());

  const double Adata1[] =
  { 1 };

  const double Adata2[] =
  { 2 };

  const double Adata3[] =
  { 4 };

  const double Adata4[] =
  { 5 };

  const double Ddata1[] =
  { 3 };

  const double Ddata2[] =
  { 6 };

  const double Bdata1[] =
  { 7 };

  const double Bdata2[] =
  { 8 };

  const double Cdata[] =
  { 9 };

  FullMatrix<double> MA1(1, 1);
  FullMatrix<double> MA2(1, 1);
  FullMatrix<double> MA3(1, 1);
  FullMatrix<double> MA4(1, 1);
  FullMatrix<double> MB1(1, 1);
  FullMatrix<double> MB2(1, 1);
  FullMatrix<double> MC(1, 1);
  FullMatrix<double> MD1(1, 1);
  FullMatrix<double> MD2(1, 1);

  MA1.fill(Adata1);
  MA2.fill(Adata2);
  MA3.fill(Adata3);
  MA4.fill(Adata4);
  MB1.fill(Bdata1);
  MB2.fill(Bdata2);
  MC.fill(Cdata);
  MD1.fill(Ddata1);
  MD2.fill(Ddata2);

  BlockMatrixArray<double> Ainv(2, 2);
  BlockMatrixArray<double> B(1, 2);
  BlockMatrixArray<double> C(1, 1);
  BlockMatrixArray<double> Dt(1, 2);

  Ainv.enter(MA1, 0, 0);
  Ainv.enter(MA2, 0, 1);
  Ainv.enter(MA3, 1, 0);
  Ainv.enter(MA4, 1, 1);
  B.enter(MB1, 0, 0);
  B.enter(MB2, 0, 1);
  C.enter(MC, 0, 0);
  Dt.enter(MD1, 0, 0);
  Dt.enter(MD2, 0, 1);

  checkPrepare_rhs<BlockMatrixArray<double>, BlockMatrixArray<double>,
                   BlockMatrixArray<double>, BlockMatrixArray<double> >(dst, src, Ainv, B,
                       Dt, C);

}
