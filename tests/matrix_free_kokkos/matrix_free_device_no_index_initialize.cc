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


// check that Portable::FEEvaluation::submit_dof_value/get_dof_value
// works correctly.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/matrix_free/portable_fe_evaluation.h>

#include "../tests.h"

#include "Kokkos_Core.hpp"

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

    // set to unit vector
    auto fe_eval_ptr = &fe_eval;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(data->team_member,
                                                 n_local_dofs),
                         [&](int i) { fe_eval_ptr->submit_dof_value(1., i); });
    data->team_member.team_barrier();
    fe_eval.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);

#ifndef __APPLE__
    Kokkos::parallel_for(Kokkos::TeamThreadRange(data->team_member,
                                                 n_local_dofs),
                         [&](int i) {
                           // values should evaluate to one, derivatives to zero
                           assert(fe_eval_ptr->get_value(i) == 1.);
                           for (unsigned int e = 0; e < dim; ++e)
                             assert(fe_eval_ptr->get_gradient(i)[e] == 0.);
                         });

    fe_eval.integrate(EvaluationFlags::values | EvaluationFlags::gradients);

    Kokkos::parallel_for(
      Kokkos::TeamThreadRange(data->team_member, n_local_dofs),
      KOKKOS_LAMBDA(int i) { assert(fe_eval_ptr->get_dof_value(i) == 1.); });
#endif
  }



  void
  test() const
  {
    LinearAlgebra::distributed::Vector<Number, MemorySpace::Default> dst_dummy;
    LinearAlgebra::distributed::Vector<Number, MemorySpace::Default> src_dummy;

    data.cell_loop(*this, src_dummy, dst_dummy);

    Kokkos::fence();

    deallog << "OK" << std::endl;
  };

protected:
  const Portable::MatrixFree<dim, Number> &data;
};

template <int dim, int fe_degree, int n_q_points_1d, typename Number>
const unsigned int
  MatrixFreeTest<dim, fe_degree, n_q_points_1d, Number>::n_local_dofs;

template <int dim, int fe_degree, int n_q_points_1d, typename Number>
const unsigned int
  MatrixFreeTest<dim, fe_degree, n_q_points_1d, Number>::n_q_points;



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
}


template <int dim, int fe_degree>
void
test()
{
  const SphericalManifold<dim> manifold;
  Triangulation<dim>           tria;
  GridGenerator::hyper_ball(tria);
  for (const auto &cell : tria.active_cell_iterators())
    for (const unsigned int f : GeometryInfo<dim>::face_indices())
      if (cell->at_boundary(f))
        cell->face(f)->set_all_manifold_ids(0);
  tria.set_manifold(0, manifold);

  // refine first and last cell
  tria.begin(tria.n_levels() - 1)->set_refine_flag();
  tria.last()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  tria.refine_global(4 - dim);

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
