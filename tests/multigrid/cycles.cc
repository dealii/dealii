// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
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

  for (unsigned int i = minlevel; i <= maxlevel; ++i)
    mg1.defect[i].reinit(N);

  {
    auto print_coarse_solve = [](const bool start, const unsigned int level) {
      if (start)
        {
          deallog << "V-cycle entering level " << level << std::endl;
          deallog << "V-cycle  Defect norm   0" << std::endl;
          deallog << "Coarse level           " << level << std::endl;
        }
    };

    auto print_restriction = [](const bool start, const unsigned int level) {
      if (start)
        {
          deallog << "Residual on      level " << level << std::endl;
          deallog << "Residual norm          0" << std::endl;
        }
    };

    auto print_prolongation = [](const bool start, const unsigned int level) {
      if (start)
        {
          deallog << "Prolongate norm        0" << std::endl;
          deallog << "V-cycle  Defect norm   0" << std::endl;
        }
    };

    auto print_pre_smoother_step = [](const bool         start,
                                      const unsigned int level) {
      if (start)
        {
          deallog << "V-cycle entering level " << level << std::endl;
          deallog << "V-cycle  Defect norm   0" << std::endl;
          deallog << "Smoothing on     level " << level << std::endl;
          deallog << "Solution norm          0" << std::endl;
        }
    };

    auto print_post_smoother_step = [](const bool         start,
                                       const unsigned int level) {
      if (start)
        {
          deallog << "Smoothing on     level " << level << std::endl;
          deallog << "Solution norm          0" << std::endl;
          deallog << "V-cycle leaving  level " << level << std::endl;
        }
    };

    const auto coarse_connection = mg1.connect_coarse_solve(print_coarse_solve);
    const auto restriction_connection =
      mg1.connect_restriction(print_restriction);
    const auto prolongation_connection =
      mg1.connect_prolongation(print_prolongation);
    const auto pre_smoother_connection =
      mg1.connect_pre_smoother_step(print_pre_smoother_step);
    const auto post_smoother_connection =
      mg1.connect_post_smoother_step(print_post_smoother_step);

    mg1.cycle();
    deallog << std::endl;

    coarse_connection.disconnect();
    restriction_connection.disconnect();
    prolongation_connection.disconnect();
    pre_smoother_connection.disconnect();
    post_smoother_connection.disconnect();
  }

  {
    auto print_coarse_solve_w = [](const bool start, const unsigned int level) {
      if (start)
        {
          deallog << "W-cycle entering level  " << level << std::endl;
          deallog << "W-cycle defect norm     0" << std::endl;
          deallog << "W-cycle coarse level    " << level << std::endl;
        }
    };

    auto print_restriction_w = [](const bool start, const unsigned int level) {
      if (start)
        {
          deallog << "W-cycle residual level  " << level << std::endl;
          deallog << "W-cycle residual norm   0" << std::endl;
        }
    };

    auto print_prolongation_w = [](const bool start, const unsigned int level) {
      if (start)
        {
          deallog << "W-cycle  Defect norm    0" << std::endl;
        }
    };

    auto print_pre_smoother_step_w = [](const bool         start,
                                        const unsigned int level) {
      if (start)
        {
          deallog << "W-cycle entering level  " << level << std::endl;
          deallog << "W-cycle defect norm     0" << std::endl;
          deallog << "W-cycle smoothing level " << level << std::endl;
          deallog << "W-cycle solution norm   0" << std::endl;
        }
    };

    auto print_post_smoother_step_w = [](const bool         start,
                                         const unsigned int level) {
      if (start)
        {
          deallog << "W-cycle smoothing level " << level << std::endl;
          deallog << "W-cycle solution norm   0" << std::endl;
          deallog << "W-cycle leaving level   " << level << std::endl;
        }
    };

    const auto coarse_connection =
      mg1.connect_coarse_solve(print_coarse_solve_w);
    const auto restriction_connection =
      mg1.connect_restriction(print_restriction_w);
    const auto prolongation_connection =
      mg1.connect_prolongation(print_prolongation_w);
    const auto pre_smoother_connection =
      mg1.connect_pre_smoother_step(print_pre_smoother_step_w);
    const auto post_smoother_connection =
      mg1.connect_post_smoother_step(print_post_smoother_step_w);

    mg1.set_cycle(Multigrid<VectorType>::w_cycle);
    mg1.cycle();
    deallog << std::endl;

    coarse_connection.disconnect();
    restriction_connection.disconnect();
    prolongation_connection.disconnect();
    pre_smoother_connection.disconnect();
    post_smoother_connection.disconnect();
  }

  {
    auto print_coarse_solve_f = [](const bool start, const unsigned int level) {
      if (start)
        {
          deallog << "F-cycle entering level  " << level << std::endl;
          deallog << "F-cycle defect norm     0" << std::endl;
          deallog << "F-cycle coarse level    " << level << std::endl;
        }
    };

    auto print_restriction_f = [](const bool start, const unsigned int level) {
      if (start)
        {
          deallog << "F-cycle residual level  " << level << std::endl;
          deallog << "F-cycle residual norm   0" << std::endl;
        }
    };

    auto print_prolongation_f = [](const bool start, const unsigned int level) {
      if (start)
        {
          deallog << "F-cycle  Defect norm    0" << std::endl;
        }
    };

    auto print_pre_smoother_step_f = [](const bool         start,
                                        const unsigned int level) {
      if (start)
        {
          deallog << "F-cycle entering level  " << level << std::endl;
          deallog << "F-cycle defect norm     0" << std::endl;
          deallog << "F-cycle smoothing level " << level << std::endl;
          deallog << "F-cycle solution norm   0" << std::endl;
        }
    };

    auto print_post_smoother_step_f = [](const bool         start,
                                         const unsigned int level) {
      if (start)
        {
          deallog << "F-cycle smoothing level " << level << std::endl;
          deallog << "F-cycle solution norm   0" << std::endl;
          deallog << "F-cycle leaving level   " << level << std::endl;
        }
    };

    mg1.connect_coarse_solve(print_coarse_solve_f);
    mg1.connect_restriction(print_restriction_f);
    mg1.connect_prolongation(print_prolongation_f);
    mg1.connect_pre_smoother_step(print_pre_smoother_step_f);
    mg1.connect_post_smoother_step(print_post_smoother_step_f);

    mg1.set_cycle(Multigrid<VectorType>::f_cycle);
    mg1.cycle();
  }
}

int
main()
{
  initlog();

  test_cycles(0, 4);
  test_cycles(2, 5);
}
