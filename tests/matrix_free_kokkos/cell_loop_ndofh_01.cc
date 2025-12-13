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


// check that Portable::MatrixFree::cell_loop works with >1 DoFHandler
// and different vectors

#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_block_vector.h>

#include <deal.II/matrix_free/portable_fe_evaluation.h>
#include <deal.II/matrix_free/portable_matrix_free.h>

#include <Kokkos_Core.hpp>

#include "../tests.h"

template <int dim, int n_q_points_1d, typename Number = double>
class MatrixFreeTest
{
public:
  static const unsigned int n_q_points = Utilities::pow(n_q_points_1d, dim);

  MatrixFreeTest(const Portable::MatrixFree<dim, Number> &data_in)
    : data(data_in){};

  DEAL_II_HOST_DEVICE void
  operator()(const typename Portable::MatrixFree<dim, Number>::Data *data,
             const Portable::DeviceVector<Number>                   &src,
             Portable::DeviceVector<Number>                         &dst) const
  {
    Portable::FEEvaluation<dim, 1, n_q_points_1d, 1, Number> fe_eval1(data, 0);
    Portable::FEEvaluation<dim, 2, n_q_points_1d, 1, Number> fe_eval2(data, 1);
    Portable::FEEvaluation<dim, 0, n_q_points_1d, 2, Number> fe_eval3(data, 2);

    if (current_index == 0)
      {
        fe_eval1.read_dof_values(src);
        fe_eval1.distribute_local_to_global(dst);
      }
    if (current_index == 1)
      {
        fe_eval2.read_dof_values(src);
        fe_eval2.distribute_local_to_global(dst);
      }
    if (current_index == 2)
      {
        fe_eval3.read_dof_values(src);
        fe_eval3.distribute_local_to_global(dst);
      }
  }

  void
  test()
  {
    LinearAlgebra::distributed::BlockVector<Number, MemorySpace::Default> dst;
    LinearAlgebra::distributed::BlockVector<Number, MemorySpace::Default> src;
    data.initialize_dof_vector(dst);
    data.initialize_dof_vector(src);

    for (unsigned int i = 0; i < 3; ++i)
      src.block(i).add(1.0 + 1.0 * i);

    for (unsigned int i = 0; i < 3; ++i)
      {
        deallog << "*** block = " << i << '\n' << std::endl;

        this->current_index = i;
        data.cell_loop(*this, src.block(i), dst.block(i));
        Kokkos::fence();
        dst.block(i).print(deallog.get_file_stream());
      }
    deallog << "OK." << std::endl;
  }

protected:
  const Portable::MatrixFree<dim, Number> &data;
  unsigned int                             current_index;
};


template <int dim, int n_q_points_1d, typename Number = double>
class MatrixFreeTestBlock
{
public:
  static const unsigned int n_q_points = Utilities::pow(n_q_points_1d, dim);

  MatrixFreeTestBlock(const Portable::MatrixFree<dim, Number> &data_in)
    : data(data_in){};

  DEAL_II_HOST_DEVICE void
  operator()(const typename Portable::MatrixFree<dim, Number>::Data *data,
             const Portable::DeviceBlockVector<Number>              &src,
             Portable::DeviceBlockVector<Number>                    &dst) const
  {
    Portable::FEEvaluation<dim, 1, n_q_points_1d, 1, Number> fe_eval1(data, 0);
    Portable::FEEvaluation<dim, 2, n_q_points_1d, 1, Number> fe_eval2(data, 1);
    Portable::FEEvaluation<dim, 0, n_q_points_1d, 2, Number> fe_eval3(data, 2);

    fe_eval1.read_dof_values(src.block(0));
    fe_eval2.read_dof_values(src.block(1));
    fe_eval3.read_dof_values(src.block(2));
    fe_eval1.distribute_local_to_global(dst.block(0));
    fe_eval2.distribute_local_to_global(dst.block(1));
    fe_eval3.distribute_local_to_global(dst.block(2));
  }

  void
  test() const
  {
    LinearAlgebra::distributed::BlockVector<Number, MemorySpace::Default> dst;
    LinearAlgebra::distributed::BlockVector<Number, MemorySpace::Default> src;

    data.initialize_dof_vector(src);
    data.initialize_dof_vector(dst);

    for (unsigned int i = 0; i < 3; ++i)
      src.block(i).add(1.0 + 1.0 * i);

    data.cell_loop(*this, src, dst);

    Kokkos::fence();

    deallog << "*** all blocks together:\n" << std::endl;

    for (unsigned int i = 0; i < 3; ++i)
      dst.block(i).print(deallog.get_file_stream());
    deallog << "OK." << std::endl;
  }

protected:
  const Portable::MatrixFree<dim, Number> &data;
};


template <int dim, int n_q_points_1d, typename Number>
const unsigned int MatrixFreeTest<dim, n_q_points_1d, Number>::n_q_points;

template <int dim, int n_q_points_1d, typename Number>
const unsigned int MatrixFreeTestBlock<dim, n_q_points_1d, Number>::n_q_points;


template <int dim, typename number>
void
do_test(const std::vector<const DoFHandler<dim> *>           &dof,
        const std::vector<const AffineConstraints<double> *> &constraints)
{
  constexpr int                     q_degree = 2 + 1;
  MappingQ<dim>                     mapping(1);
  Portable::MatrixFree<dim, number> mf_data;
  {
    const QGauss<1>                                            quad(q_degree);
    typename Portable::MatrixFree<dim, number>::AdditionalData data;
    data.mapping_update_flags = update_values | update_gradients |
                                update_JxW_values | update_quadrature_points;
    mf_data.reinit(mapping, dof, constraints, quad, data);
  }

  MatrixFreeTest<dim, q_degree, number> mf(mf_data);
  mf.test();
  MatrixFreeTestBlock<dim, q_degree, number> mfb(mf_data);
  mfb.test();
}


template <int dim>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);

  FE_Q<dim>     fe1(1);
  FE_Q<dim>     fe2(2);
  FESystem<dim> fe3(FE_DGQ<dim>(0), 2);

  FESystem<dim>   fe(fe1, 1, fe2, 1, fe3, 1);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);

  DoFHandler<dim> dof1(tria), dof2(tria), dof3(tria);
  dof1.distribute_dofs(fe1);
  dof2.distribute_dofs(fe2);
  dof3.distribute_dofs(fe3);

  std::vector<const DoFHandler<dim> *> dof_handler_vec = {&dof1, &dof2, &dof3};

  AffineConstraints<double> constraints1, constraints2, constraints3;
  DoFTools::make_hanging_node_constraints(dof1, constraints1);
  DoFTools::make_hanging_node_constraints(dof2, constraints2);
  DoFTools::make_hanging_node_constraints(dof3, constraints3);
  constraints1.close();
  constraints2.close();
  constraints3.close();

  std::vector<const AffineConstraints<double> *> constraint_vec = {
    &constraints1, &constraints2, &constraints3};

  deallog << "Testing " << fe.get_name() << std::endl;

  do_test<dim, double>(dof_handler_vec, constraint_vec);
}



int
main()
{
  initlog();

  Kokkos::initialize();

  test<2>();

  Kokkos::finalize();

  return 0;
}
