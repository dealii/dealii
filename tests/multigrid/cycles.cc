// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2017 by the deal.II authors
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

#include <deal.II/base/mg_level_object.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/multigrid/mg_base.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/multigrid.h>

#include "../tests.h"


#define N 3
typedef Vector<double> VectorType;

class MGAll : public MGSmootherBase<VectorType>,
              public MGTransferBase<VectorType>,
              public MGCoarseGridBase<VectorType>
{
public:
  virtual ~MGAll()
  {}

  virtual void
  smooth(const unsigned int, VectorType &, const VectorType &) const
  {}

  virtual void
  prolongate(const unsigned int, VectorType &, const VectorType &) const
  {}

  virtual void
  restrict_and_add(const unsigned int, VectorType &, const VectorType &) const
  {}

  virtual void
  clear()
  {}

  virtual void
  operator()(const unsigned int, VectorType &, const VectorType &) const
  {}
};

void
test_cycles(unsigned int minlevel, unsigned int maxlevel)
{
  MGAll                             all;
  MGLevelObject<FullMatrix<double>> level_matrices(0, maxlevel);
  for (unsigned int i = 0; i <= maxlevel; ++i)
    level_matrices[i].reinit(N, N);
  mg::Matrix<VectorType> mgmatrix(level_matrices);

  Multigrid<VectorType> mg1(mgmatrix,
                            all,
                            all,
                            all,
                            all,
                            minlevel,
                            maxlevel,
                            Multigrid<VectorType>::v_cycle);
  mg1.set_debug(3);
  for (unsigned int i = minlevel; i <= maxlevel; ++i)
    mg1.defect[i].reinit(N);
  mg1.cycle();
  deallog << std::endl;

  mg1.set_cycle(Multigrid<VectorType>::w_cycle);
  mg1.cycle();
  deallog << std::endl;

  mg1.set_cycle(Multigrid<VectorType>::f_cycle);
  mg1.cycle();
}

int
main()
{
  initlog();

  test_cycles(0, 4);
  test_cycles(2, 5);
}
