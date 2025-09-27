// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check that Portable::MatrixFree::cell_loop works with BlockVector

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/lac/la_parallel_block_vector.h>

#include <deal.II/matrix_free/portable_fe_evaluation.h>
#include <deal.II/matrix_free/portable_matrix_free.h>

#include <Kokkos_Core.hpp>

#include "../tests.h"

template <int dim,
          int fe_degree,
          int n_q_points_1d = fe_degree + 1,
          typename Number   = double>
class MatrixFreeTest
{
public:
  static const unsigned int n_local_dofs = Utilities::pow(fe_degree + 1, dim);
  static const unsigned int n_q_points   = Utilities::pow(n_q_points_1d, dim);

  MatrixFreeTest(const Portable::MatrixFree<dim, Number> &data_in)
    : data(data_in){};

  DEAL_II_HOST_DEVICE void
  operator()(const typename Portable::MatrixFree<dim, Number>::Data *data,
             const Portable::DeviceVector<Number>                   &src,
             Portable::DeviceVector<Number>                         &dst) const
  {
    Portable::FEEvaluation<dim, fe_degree, n_q_points_1d, 1, Number> fe_eval(
      data);

    fe_eval.read_dof_values(src);
    fe_eval.distribute_local_to_global(dst);
  }

  void
  test() const
  {
    LinearAlgebra::distributed::Vector<Number, MemorySpace::Default> dst;
    LinearAlgebra::distributed::Vector<Number, MemorySpace::Default> src;

    data.initialize_dof_vector(dst);
    data.initialize_dof_vector(src);
    src.add(1.0);

    data.cell_loop(*this, src, dst);

    Kokkos::fence();

    deallog << "OK:" << dst.linfty_norm() << std::endl;
  }

protected:
  const Portable::MatrixFree<dim, Number> &data;
};


template <int dim,
          int fe_degree,
          int n_q_points_1d = fe_degree + 1,
          typename Number   = double>
class MatrixFreeTestBlock
{
public:
  static const unsigned int n_local_dofs = Utilities::pow(fe_degree + 1, dim);
  static const unsigned int n_q_points   = Utilities::pow(n_q_points_1d, dim);

  MatrixFreeTestBlock(const Portable::MatrixFree<dim, Number> &data_in)
    : data(data_in){};

  DEAL_II_HOST_DEVICE void
  operator()(const typename Portable::MatrixFree<dim, Number>::Data *data,
             const Portable::DeviceBlockVector<Number>              &src,
             Portable::DeviceBlockVector<Number>                    &dst) const
  {
    Portable::FEEvaluation<dim, fe_degree, n_q_points_1d, 1, Number> fe_eval(
      data);

    fe_eval.read_dof_values(src.block(0));
    fe_eval.distribute_local_to_global(dst.block(0));
  }

  void
  test() const
  {
    LinearAlgebra::distributed::BlockVector<Number, MemorySpace::Default> dst(
      1);
    LinearAlgebra::distributed::BlockVector<Number, MemorySpace::Default> src(
      1);

    data.initialize_dof_vector(dst.block(0));
    data.initialize_dof_vector(src.block(0));
    src.add(1.0);

    data.cell_loop(*this, src, dst);

    Kokkos::fence();

    deallog << "OK:" << dst.block(0).linfty_norm() << std::endl;
  }

protected:
  const Portable::MatrixFree<dim, Number> &data;
};

template <int dim, int fe_degree, int n_q_points_1d, typename Number>
const unsigned int
  MatrixFreeTest<dim, fe_degree, n_q_points_1d, Number>::n_local_dofs;

template <int dim, int fe_degree, int n_q_points_1d, typename Number>
const unsigned int
  MatrixFreeTest<dim, fe_degree, n_q_points_1d, Number>::n_q_points;

template <int dim, int fe_degree, int n_q_points_1d, typename Number>
const unsigned int
  MatrixFreeTestBlock<dim, fe_degree, n_q_points_1d, Number>::n_local_dofs;

template <int dim, int fe_degree, int n_q_points_1d, typename Number>
const unsigned int
  MatrixFreeTestBlock<dim, fe_degree, n_q_points_1d, Number>::n_q_points;


template <int dim, int fe_degree, typename number>
void
do_test(const DoFHandler<dim>           &dof,
        const AffineConstraints<double> &constraints)
{
  Portable::MatrixFree<dim, number> mf_data;
  {
    const QGauss<1> quad(fe_degree + 1);
    typename Portable::MatrixFree<dim, number>::AdditionalData data;
    data.mapping_update_flags = update_values | update_gradients |
                                update_JxW_values | update_quadrature_points;
    mf_data.reinit(dof, constraints, quad, data);
  }

  deallog << "Testing " << dof.get_fe().get_name() << std::endl;
  MatrixFreeTest<dim, fe_degree, fe_degree + 1, number> mf(mf_data);
  mf.test();
  MatrixFreeTestBlock<dim, fe_degree, fe_degree + 1, number> mfb(mf_data);
  mfb.test();
}


template <int dim, int fe_degree>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);

  // refine first and last cell
  tria.begin(tria.n_levels() - 1)->set_refine_flag();
  tria.last()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);

  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dof, constraints);
  constraints.close();

  do_test<dim, fe_degree, double>(dof, constraints);
}



int
main()
{
  initlog();

  Kokkos::initialize();

  test<2, 1>();

  Kokkos::finalize();

  return 0;
}
