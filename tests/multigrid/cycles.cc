// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
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


/*
 * Test the general behavior of the multigrid cycles without any
 * numerics. Therefore, all transfer operators are void and we use the
 * same matrix on each level.
 */

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/multigrid/mg_base.h>
#include <deal.II/multigrid/multigrid.h>
#include <deal.II/multigrid/mg_matrix.h>

#include <fstream>

#define N 3
typedef Vector<double> VECTOR;

class MGAll
  :
  public MGSmootherBase<VECTOR>,
  public MGTransferBase<VECTOR>,
  public MGCoarseGridBase<VECTOR>
{
public:
  virtual ~MGAll()
  {}

  virtual void smooth (const unsigned int,
                       VECTOR &, const VECTOR &) const
  {}

  virtual void prolongate (const unsigned int,
                           VECTOR &, const VECTOR &) const
  {}

  virtual void restrict_and_add (const unsigned int,
                                 VECTOR &, const VECTOR &) const
  {}

  virtual void clear ()
  {}

  virtual void operator() (const unsigned int,
                           VECTOR &, const VECTOR &) const
  {}
};

void test_cycles(unsigned int minlevel, unsigned int maxlevel)
{
  MGAll all;
  MGLevelObject<FullMatrix<double> > level_matrices(0, maxlevel);
  for (unsigned int i=0; i<=maxlevel; ++i)
    level_matrices[i].reinit(N, N);
  MGMatrix<FullMatrix<double>, VECTOR> mgmatrix(&level_matrices);

  Multigrid<VECTOR> mg1(minlevel, maxlevel, mgmatrix, all, all, all, all,
                        Multigrid<VECTOR>::v_cycle);
  mg1.set_debug(3);
  for (unsigned int i=minlevel; i<=maxlevel; ++i)
    mg1.defect[i].reinit(N);
  mg1.cycle();
  deallog << std::endl;

  mg1.set_cycle(Multigrid<VECTOR>::w_cycle);
  mg1.cycle();
  deallog << std::endl;

  mg1.set_cycle(Multigrid<VECTOR>::f_cycle);
  mg1.cycle();
}

int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test_cycles (0,4);
  test_cycles (2,5);
}
