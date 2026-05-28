// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2019 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// Solve a Stokes problem with Portable::MatrixFree and a good block
// preconditioner. Otherwise, like stokes_02.cc

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/repartitioning_policy_tools.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_idr.h>

#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/matrix_free/portable_fe_evaluation.h>
#include <deal.II/matrix_free/portable_matrix_free.h>

#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/multigrid.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/vector_tools_integrate_difference.h>

#include "../tests.h"


using namespace dealii;

// Manufactured solution
//
// 2D:
// u = pi * sx^2 * sin(2 * pi * y)
// v = -pi * sin(2 * pi * x) * sy^2
// p = cos(pi * x) * cos(pi * y)
// where sx = sin(pi*x), sy = sin(pi*y), sz = sin(pi*z)

template <int dim>
class VelocityRightHandSide : public Function<dim>
{
public:
  VelocityRightHandSide()
    : Function<dim>(dim)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const override;
};

template <int dim>
double
VelocityRightHandSide<dim>::value(const Point<dim>  &p,
                                  const unsigned int component) const
{
  AssertIndexRange(component, dim);

  const double x   = p[0];
  const double y   = p[1];
  const double pi  = numbers::PI;
  const double pi2 = pi * pi;

  const double sx = std::sin(pi * x);
  const double cx = std::cos(pi * x);
  const double sy = std::sin(pi * y);
  const double cy = std::cos(pi * y);

  if constexpr (dim == 2)
    {
      if (component == 0)
        return pi * cy * (16.0 * pi2 * sx * sx * sy - 4.0 * pi2 * sy - sx);
      else
        return pi * cx * (-16.0 * pi2 * sx * sy * sy + 4.0 * pi2 * sx - sy);
    }
  else
    {
      const double z   = p[2];
      const double sz  = std::sin(pi * z);
      const double cz  = std::cos(pi * z);
      const double sx2 = sx * sx;
      const double sy2 = sy * sy;
      const double sz2 = sz * sz;
      const double cx2 = cx * cx;
      const double cy2 = cy * cy;
      const double cz2 = cz * cz;
      const double pi3 = pi * pi2;

      if (component == 0)
        return -16.0 * pi3 * sx2 * sy2 * sz * cz +
               16.0 * pi3 * sx2 * sy * sz2 * cy -
               4.0 * pi3 * sx2 * sy * cy * cz2 +
               4.0 * pi3 * sx2 * sz * cy2 * cz - pi * sx * cy * cz +
               4.0 * pi3 * sy2 * sz * cx2 * cz -
               4.0 * pi3 * sy * sz2 * cx2 * cy;
      else if (component == 1)
        return 16.0 * pi3 * sx2 * sy2 * sz * cz -
               4.0 * pi3 * sx2 * sz * cy2 * cz -
               16.0 * pi3 * sx * sy2 * sz2 * cx +
               4.0 * pi3 * sx * sy2 * cx * cz2 +
               4.0 * pi3 * sx * sz2 * cx * cy2 -
               4.0 * pi3 * sy2 * sz * cx2 * cz - pi * sy * cx * cz;
      else
        return -16.0 * pi3 * sx2 * sy * sz2 * cy +
               4.0 * pi3 * sx2 * sy * cy * cz2 +
               16.0 * pi3 * sx * sy2 * sz2 * cx -
               4.0 * pi3 * sx * sy2 * cx * cz2 -
               4.0 * pi3 * sx * sz2 * cx * cy2 +
               4.0 * pi3 * sy * sz2 * cx2 * cy - pi * sz * cx * cy;
    }
}

template <int dim>
class VelocitySolution : public Function<dim>
{
public:
  VelocitySolution()
    : Function<dim>(dim)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const override;
};

template <int dim>
double
VelocitySolution<dim>::value(const Point<dim>  &p,
                             const unsigned int component) const
{
  AssertIndexRange(component, dim);

  const double x  = p[0];
  const double y  = p[1];
  const double pi = numbers::PI;

  const double sx  = std::sin(pi * x);
  const double sy  = std::sin(pi * y);
  const double s2x = std::sin(2.0 * pi * x);
  const double s2y = std::sin(2.0 * pi * y);

  if constexpr (dim == 2)
    {
      if (component == 0)
        return pi * sx * sx * s2y;
      else
        return -pi * s2x * sy * sy;
    }
  else
    {
      const double z      = p[2];
      const double sz     = std::sin(pi * z);
      const double s2z    = std::sin(2.0 * pi * z);
      const double sx2    = sx * sx;
      const double sy2    = sy * sy;
      const double sz2    = sz * sz;
      const double dphidx = pi * s2x * sy2 * sz2;
      const double dphidy = pi * s2y * sx2 * sz2;
      const double dphidz = pi * s2z * sx2 * sy2;
      if (component == 0)
        return dphidy - dphidz;
      else if (component == 1)
        return dphidz - dphidx;
      else
        return dphidx - dphidy;
    }
}

template <int dim>
class PressureSolution : public Function<dim>
{
public:
  PressureSolution()
    : Function<dim>(1)
  {}

  virtual double
  value(const Point<dim> &p,
        const unsigned int /*component*/ = 0) const override
  {
    const double pi = numbers::PI;
    if constexpr (dim == 2)
      return std::cos(pi * p[0]) * std::cos(pi * p[1]);
    else
      return std::cos(pi * p[0]) * std::cos(pi * p[1]) * std::cos(pi * p[2]);
  }
};


// Velocity Block

template <int dim,
          int degree_u,
          int degree_p,
          typename Number,
          int n_q_points_1d>
class VelocityCellOperator
{
public:
  static const unsigned int n_q_points =
    dealii::Utilities::pow(n_q_points_1d, dim);

  DEAL_II_HOST_DEVICE void
  operator()(const typename Portable::MatrixFree<dim, Number>::Data *data,
             const Portable::DeviceVector<Number>                   &src,
             Portable::DeviceVector<Number>                         &dst) const;
};

template <int dim, int fe_degree>
class VelocityOperatorQuad
{
public:
  DEAL_II_HOST_DEVICE void
  operator()(
    Portable::FEEvaluation<dim, fe_degree, fe_degree + 1, dim, double> *fe_eval,
    const int q_point) const
  {
    const auto gradient_u = fe_eval->get_gradient(q_point);
    fe_eval->submit_gradient(gradient_u, q_point);
  }
};

template <int dim,
          int degree_u,
          int degree_p,
          typename Number,
          int n_q_points_1d>
DEAL_II_HOST_DEVICE void
VelocityCellOperator<dim, degree_u, degree_p, Number, n_q_points_1d>::
operator()(const typename Portable::MatrixFree<dim, Number>::Data *data,
           const Portable::DeviceVector<Number>                   &src,
           Portable::DeviceVector<Number>                         &dst) const
{
  Portable::FEEvaluation<dim, degree_u, n_q_points_1d, dim, double> fe_u(data,
                                                                         0);

  fe_u.read_dof_values(src);
  fe_u.evaluate(EvaluationFlags::gradients);

  VelocityOperatorQuad<dim, degree_u> quad_operation;

  data->for_each_quad_point(
    [&](const int q_point) { quad_operation(&fe_u, q_point); });

  fe_u.integrate(EvaluationFlags::gradients);
  fe_u.distribute_local_to_global(dst);
}

template <int dim,
          int degree_u,
          int degree_p,
          typename Number = double,
          typename VectorType =
            LinearAlgebra::distributed::Vector<double, MemorySpace::Default>,
          int n_q_points_1d = degree_u + 1>
class PortableMFVelocityOperator : public EnableObserverPointer

{
public:
  PortableMFVelocityOperator()
  {}

  PortableMFVelocityOperator(
    std::shared_ptr<Portable::MatrixFree<dim, Number>> data_in)
    : data(data_in)
  {}

  void
  reinit(std::shared_ptr<Portable::MatrixFree<dim, Number>> data_in)
  {
    data = data_in;
  }

  void
  initialize_dof_vector(VectorType &vec) const
  {
    data->initialize_dof_vector(vec, 0 /* velocity */);
  }


  types::global_dof_index
  m() const
  {
    return data->get_vector_partitioner(0 /* velocity */)->size();
  }


  std::shared_ptr<DiagonalMatrix<
    LinearAlgebra::distributed::Vector<double, MemorySpace::Default>>>
  get_matrix_diagonal_inverse() const
  {
    return inverse_diagonal_entries;
  }

  double
  el(const types::global_dof_index row, const types::global_dof_index col) const
  {
    (void)col;
    Assert(row == col, ExcNotImplemented());
    Assert(inverse_diagonal_entries.get() != nullptr &&
             inverse_diagonal_entries->m() > 0,
           ExcNotInitialized());
    return 1.0 / (*inverse_diagonal_entries)(row, row);
  }

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    dst = static_cast<Number>(0.);
    VelocityCellOperator<dim, degree_u, degree_p, Number, n_q_points_1d>
      velocity_operator;
    data->cell_loop(velocity_operator, src, dst);

    data->copy_constrained_values(src, dst, 0 /* velocity */);
  }

  void
  Tvmult(VectorType &dst, const VectorType &src) const
  {
    AssertThrow(false, ExcNotImplemented());

    (void)dst;
    (void)src;
  }

  void
  compute_diagonal()
  {
    Assert(data.get() != nullptr, ExcNotInitialized());

    this->inverse_diagonal_entries.reset(
      new DiagonalMatrix<
        LinearAlgebra::distributed::Vector<double, MemorySpace::Default>>());
    LinearAlgebra::distributed::Vector<double, MemorySpace::Default>
      &inverse_diagonal = inverse_diagonal_entries->get_vector();
    data->initialize_dof_vector(inverse_diagonal, 0 /* velocity */);
    VelocityOperatorQuad<dim, degree_u> velocity_operator_quad;

    MatrixFreeTools::compute_diagonal<dim, degree_u, degree_u + 1, dim, double>(
      *data.get(),
      inverse_diagonal,
      velocity_operator_quad,
      EvaluationFlags::gradients,
      EvaluationFlags::gradients,
      0 /*velocity*/);

    double *raw_diagonal = inverse_diagonal.get_values();

    Kokkos::parallel_for(
      "invert A diagonal",
      inverse_diagonal.locally_owned_size(),
      KOKKOS_LAMBDA(int i) {
        Assert(raw_diagonal[i] > 0.,
               ExcMessage("No diagonal entry in a positive definite operator "
                          "should be zero"));
        raw_diagonal[i] = 1. / raw_diagonal[i];
      });
  }

private:
  std::shared_ptr<Portable::MatrixFree<dim, Number>> data;
  std::shared_ptr<DiagonalMatrix<
    LinearAlgebra::distributed::Vector<double, MemorySpace::Default>>>
    inverse_diagonal_entries;
};



// Mass Operator


template <int dim,
          int degree_u,
          int degree_p,
          typename Number,
          int n_q_points_1d>
class MassOperatorQuad
{
public:
  DEAL_II_HOST_DEVICE void
  operator()(
    Portable::FEEvaluation<dim, degree_p, n_q_points_1d, 1, Number> *fe_eval,
    const int q_point) const
  {
    fe_eval->submit_value(fe_eval->get_value(q_point), q_point);
  }
};

template <int dim,
          int degree_u,
          int degree_p,
          typename Number,
          int n_q_points_1d>
class MassCellOperator
{
public:
  static const unsigned int n_q_points =
    dealii::Utilities::pow(n_q_points_1d, dim);

  DEAL_II_HOST_DEVICE void
  operator()(const typename Portable::MatrixFree<dim, Number>::Data *data,
             const Portable::DeviceVector<Number>                   &src,
             Portable::DeviceVector<Number>                         &dst) const;
};

template <int dim,
          int degree_u,
          int degree_p,
          typename Number,
          int n_q_points_1d>
DEAL_II_HOST_DEVICE void
MassCellOperator<dim, degree_u, degree_p, Number, n_q_points_1d>::operator()(
  const typename Portable::MatrixFree<dim, Number>::Data *data,
  const Portable::DeviceVector<Number>                   &src,
  Portable::DeviceVector<Number>                         &dst) const
{
  Portable::FEEvaluation<dim, degree_p, n_q_points_1d, 1> fe_p(data, 1);

  fe_p.read_dof_values(src);
  fe_p.evaluate(EvaluationFlags::values);

  MassOperatorQuad<dim, degree_u, degree_p, Number, n_q_points_1d>
    quad_operation;
  data->for_each_quad_point(
    [&](const int &q_point) { quad_operation(&fe_p, q_point); });

  fe_p.integrate(EvaluationFlags::values);
  fe_p.distribute_local_to_global(dst);
}

template <int dim,
          int degree_u,
          int degree_p,
          typename Number = double,
          typename VectorType =
            LinearAlgebra::distributed::Vector<double, MemorySpace::Default>,
          int n_q_points_1d = degree_u + 1>
class PortableMFMassOperator : public EnableObserverPointer
{
public:
  PortableMFMassOperator(const Portable::MatrixFree<dim, double> &data_in)
    : data(data_in)
  {}

  types::global_dof_index
  m() const
  {
    return data.get_vector_partitioner(1 /* pressure */)->size();
  }

  double
  el(const types::global_dof_index row, const types::global_dof_index col) const
  {
    (void)col;
    Assert(row == col, ExcNotImplemented());
    Assert(inverse_diagonal_entries.get() != nullptr &&
             inverse_diagonal_entries->m() > 0,
           ExcNotInitialized());
    return 1.0 / (*inverse_diagonal_entries)(row, row);
  }

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    dst = static_cast<Number>(0.);
    MassCellOperator<dim, degree_u, degree_p, Number, n_q_points_1d>
      mass_operator;
    data.cell_loop(mass_operator, src, dst);

    data.copy_constrained_values(src, dst, 1 /* pressure */);
  }

  std::shared_ptr<DiagonalMatrix<
    LinearAlgebra::distributed::Vector<double, MemorySpace::Default>>>
  get_matrix_diagonal_inverse() const
  {
    return inverse_diagonal_entries;
  }

  void
  compute_diagonal()
  {
    this->inverse_diagonal_entries.reset(
      new DiagonalMatrix<
        LinearAlgebra::distributed::Vector<double, MemorySpace::Default>>());
    LinearAlgebra::distributed::Vector<double, MemorySpace::Default>
      &inverse_diagonal = inverse_diagonal_entries->get_vector();

    MassOperatorQuad<dim, degree_u, degree_p, Number, n_q_points_1d>
      quad_operation;

    MatrixFreeTools::compute_diagonal<dim, degree_p, n_q_points_1d, 1, Number>(
      data,
      inverse_diagonal,
      quad_operation,
      EvaluationFlags::values,
      EvaluationFlags::values,
      1 /* pressure */);

    double *raw_diagonal = inverse_diagonal.get_values();

    Kokkos::parallel_for(
      "invert Mass diagonal",
      inverse_diagonal.locally_owned_size(),
      KOKKOS_LAMBDA(int i) {
        Assert(raw_diagonal[i] > 0.,
               ExcMessage("No diagonal entry in a positive definite operator "
                          "should be zero"));
        raw_diagonal[i] = 1. / raw_diagonal[i];
      });
  }

private:
  const Portable::MatrixFree<dim, Number> &data;
  std::shared_ptr<DiagonalMatrix<
    LinearAlgebra::distributed::Vector<double, MemorySpace::Default>>>
    inverse_diagonal_entries;
};

// Stokes Operator


template <int dim,
          int degree_u,
          int degree_p,
          typename Number,
          int n_q_points_1d>
class StokesCellOperator
{
public:
  static const unsigned int n_q_points =
    dealii::Utilities::pow(n_q_points_1d, dim);

  DEAL_II_HOST_DEVICE void
  operator()(const typename Portable::MatrixFree<dim, Number>::Data *data,
             const Portable::DeviceBlockVector<Number>              &src,
             Portable::DeviceBlockVector<Number>                    &dst) const;
};

template <int dim,
          int degree_u,
          int degree_p,
          typename Number,
          int n_q_points_1d>
DEAL_II_HOST_DEVICE void
StokesCellOperator<dim, degree_u, degree_p, Number, n_q_points_1d>::operator()(
  const typename Portable::MatrixFree<dim, Number>::Data *data,
  const Portable::DeviceBlockVector<Number>              &src,
  Portable::DeviceBlockVector<Number>                    &dst) const
{
  Portable::FEEvaluation<dim, degree_u, n_q_points_1d, dim> fe_u(data, 0);
  Portable::FEEvaluation<dim, degree_p, n_q_points_1d, 1>   fe_p(data, 1);

  fe_u.read_dof_values(src.block(0));
  fe_p.read_dof_values(src.block(1));
  fe_u.evaluate(EvaluationFlags::gradients);
  fe_p.evaluate(EvaluationFlags::values);

  data->for_each_quad_point([&](const int &q_point) {
    const Tensor<2, dim, Number> gradient_u = fe_u.get_gradient(q_point);
    Tensor<2, dim, Number>       vel_term   = gradient_u;
    for (unsigned int d = 0; d < dim; ++d)
      vel_term[d][d] -= fe_p.get_value(q_point);
    fe_u.submit_gradient(vel_term, q_point);

    const Number pressure_term = trace(gradient_u);
    fe_p.submit_value(pressure_term, q_point);
  });

  fe_u.integrate(EvaluationFlags::gradients);
  fe_p.integrate(EvaluationFlags::values);
  fe_u.distribute_local_to_global(dst.block(0));
  fe_p.distribute_local_to_global(dst.block(1));
}

template <
  int dim,
  int degree_u,
  int degree_p,
  typename Number = double,
  typename VectorType =
    LinearAlgebra::distributed::BlockVector<double, MemorySpace::Default>,
  int n_q_points_1d = degree_u + 1>
class PortableMFStokesOperator
{
public:
  PortableMFStokesOperator(const Portable::MatrixFree<dim, double> &data_in)
    : data(data_in)
  {}

  const Portable::MatrixFree<dim, Number> &data;

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    dst = static_cast<Number>(0.);
    StokesCellOperator<dim, degree_u, degree_p, Number, n_q_points_1d>
      stokes_operator;
    data.cell_loop(stokes_operator, src, dst);

    data.copy_constrained_values(src, dst);
  }
};

// BT Operator


template <int dim,
          int degree_u,
          int degree_p,
          typename Number,
          int n_q_points_1d>
class BTCellOperator
{
public:
  static const unsigned int n_q_points =
    dealii::Utilities::pow(n_q_points_1d, dim);

  DEAL_II_HOST_DEVICE void
  operator()(const typename Portable::MatrixFree<dim, Number>::Data *data,
             const Portable::DeviceBlockVector<Number>              &src,
             Portable::DeviceBlockVector<Number>                    &dst) const;
};

template <int dim,
          int degree_u,
          int degree_p,
          typename Number,
          int n_q_points_1d>
DEAL_II_HOST_DEVICE void
BTCellOperator<dim, degree_u, degree_p, Number, n_q_points_1d>::operator()(
  const typename Portable::MatrixFree<dim, Number>::Data *data,
  const Portable::DeviceBlockVector<Number>              &src,
  Portable::DeviceBlockVector<Number>                    &dst) const
{
  Portable::FEEvaluation<dim, degree_u, n_q_points_1d, dim> fe_u(data, 0);
  Portable::FEEvaluation<dim, degree_p, n_q_points_1d, 1>   fe_p(data, 1);

  fe_p.read_dof_values(src.block(1));
  fe_p.evaluate(EvaluationFlags::values);

  data->for_each_quad_point([&](const int &q_point) {
    Tensor<2, dim, Number> vel_term;
    for (unsigned int d = 0; d < dim; ++d)
      vel_term[d][d] = -fe_p.get_value(q_point);
    fe_u.submit_gradient(vel_term, q_point);
  });

  fe_u.integrate(EvaluationFlags::gradients);
  fe_u.distribute_local_to_global(dst.block(0));
}

template <
  int dim,
  int degree_u,
  int degree_p,
  typename Number = double,
  typename VectorType =
    LinearAlgebra::distributed::BlockVector<double, MemorySpace::Default>,
  int n_q_points_1d = degree_u + 1>
class PortableMFBTOperator
{
public:
  PortableMFBTOperator(const Portable::MatrixFree<dim, double> &data_in)
    : data(data_in)
  {}

  const Portable::MatrixFree<dim, Number> &data;

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    dst = static_cast<Number>(0.);
    BTCellOperator<dim, degree_u, degree_p, Number, n_q_points_1d>
      cell_operator;
    data.cell_loop(cell_operator, src, dst);

    data.set_constrained_values(0.0, dst.block(0), 0 /* velocity */);
  }
};



// Preconditioner:


template <class AInvOperator,
          class SInvOperator,
          class BTOperator,
          class VectorType>
class BlockSchurPreconditioner : public
#if DEAL_II_VERSION_GTE(9, 7, 0)
                                 EnableObserverPointer
#else
                                 Subscriptor
#endif

{
public:
  /**
   * @brief Constructor
   * @param A_inverse_operator Approximation of the inverse of the velocity block.
   * @param S_inverse_operator Approximation for the inverse Schur complement.
   * @param BT_operator Operator for the B^T block of the Stokes system.
   */
  BlockSchurPreconditioner(const AInvOperator &A_inverse_operator,
                           const SInvOperator &S_inverse_operator,
                           const BTOperator   &BT_operator);

  /**
   * Matrix vector product with this preconditioner object.
   */
  void
  vmult(VectorType &dst, const VectorType &src) const;

private:
  /**
   * References to the various operators this preconditioner works with.
   */

  mutable VectorType  tmp;
  const AInvOperator &A_inverse_operator;
  const SInvOperator &S_inverse_operator;
  const BTOperator   &BT_operator;
};


template <class AInvOperator,
          class SInvOperator,
          class BTOperator,
          class VectorType>
BlockSchurPreconditioner<AInvOperator, SInvOperator, BTOperator, VectorType>::
  BlockSchurPreconditioner(const AInvOperator &A_inverse_operator,
                           const SInvOperator &S_inverse_operator,
                           const BTOperator   &BT_operator)
  : A_inverse_operator(A_inverse_operator)
  , S_inverse_operator(S_inverse_operator)
  , BT_operator(BT_operator)
{}



template <class AInvOperator,
          class SInvOperator,
          class BTOperator,
          class VectorType>
void
BlockSchurPreconditioner<AInvOperator, SInvOperator, BTOperator, VectorType>::
  vmult(VectorType &dst, const VectorType &src) const
{
  if (tmp.size() == 0)
    tmp.reinit(src);

  // first apply the Schur Complement inverse operator.
  {
    // Zero out output vector, because the Chebychev smoother otherwise
    // incorrectly uses existing values:
    dst.block(1) = 0.0;
    S_inverse_operator.vmult(dst.block(1), src.block(1));
    dst.block(1) *= -1.0;
  }

  // Apply the top right block:
  {
    BT_operator.vmult(tmp, dst);
    tmp.block(0) *= -1.0;
    tmp.block(0) += src.block(0);
  }

  // Finally the velocity block:
  A_inverse_operator.vmult(dst.block(0), tmp.block(0));
}



template <int dim, int degree_p>
void
test(unsigned int n_refinements)
{
  using Number = double;

  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(n_refinements);

  const unsigned int degree_u = degree_p + 1;

  FESystem<dim>   fe_u(FE_Q<dim>(degree_u), dim);
  FE_Q<dim>       fe_p(degree_p);
  DoFHandler<dim> dof_u(tria);
  DoFHandler<dim> dof_p(tria);
  dof_u.distribute_dofs(fe_u);
  dof_p.distribute_dofs(fe_p);

  deallog << "refinement: " << n_refinements
          << ", n_dofs: " << dof_u.n_dofs() + dof_p.n_dofs() << std::endl;

  const IndexSet &owned_set_u = dof_u.locally_owned_dofs();
  const IndexSet  relevant_set_u =
    DoFTools::extract_locally_relevant_dofs(dof_u);
  AffineConstraints<double> constraints_u(owned_set_u, relevant_set_u);
  DoFTools::make_hanging_node_constraints(dof_u, constraints_u);
  VectorTools::interpolate_boundary_values(dof_u,
                                           0,
                                           Functions::ZeroFunction<dim>(dim),
                                           constraints_u);
  constraints_u.close();

  const IndexSet &owned_set_p = dof_p.locally_owned_dofs();
  const IndexSet  relevant_set_p =
    DoFTools::extract_locally_relevant_dofs(dof_p);
  AffineConstraints<double> constraints_p(owned_set_p, relevant_set_p);
  DoFTools::make_hanging_node_constraints(dof_p, constraints_p);
  constraints_p.close();

  std::vector<const DoFHandler<dim> *> dof_handlers          = {&dof_u, &dof_p};
  std::vector<const AffineConstraints<double> *> constraints = {&constraints_u,
                                                                &constraints_p};

  MappingQ<dim>                                      mapping(degree_p);
  std::shared_ptr<Portable::MatrixFree<dim, Number>> mf_data =
    std::make_shared<Portable::MatrixFree<dim, Number>>();
  const QGauss<1>                                            quad(degree_p + 2);
  typename Portable::MatrixFree<dim, Number>::AdditionalData additional_data;
  additional_data.mapping_update_flags = update_values | update_gradients;
  mf_data->reinit(mapping, dof_handlers, constraints, quad, additional_data);

  PortableMFStokesOperator<dim, degree_u, degree_p> stokes_operator(
    *mf_data.get());

  LinearAlgebra::distributed::BlockVector<Number, MemorySpace::Default>
    solution;
  mf_data->initialize_dof_vector(solution);
  LinearAlgebra::distributed::BlockVector<Number, MemorySpace::Default> rhs;
  mf_data->initialize_dof_vector(rhs);

  LinearAlgebra::distributed::BlockVector<Number, MemorySpace::Host>
    solution_host;
  mf_data->initialize_dof_vector(solution_host);

  LinearAlgebra::distributed::BlockVector<Number, MemorySpace::Host> rhs_host;
  mf_data->initialize_dof_vector(rhs_host);

  VectorTools::create_right_hand_side(dof_u,
                                      QGauss<dim>(degree_u + 2),
                                      VelocityRightHandSide<dim>(),
                                      rhs_host.block(0),
                                      constraints_u);

  rhs.block(0).import_elements(rhs_host.block(0), VectorOperation::insert);
  rhs.block(1).import_elements(rhs_host.block(1), VectorOperation::insert);

  SolverControl solver_control(100, 1e-6 * rhs.l2_norm());

  using BlockVectorType =
    LinearAlgebra::distributed::BlockVector<Number, MemorySpace::Default>;
  SolverGMRES<BlockVectorType> solver(
    solver_control,
    typename SolverGMRES<BlockVectorType>::AdditionalData(30, true));

  using VectorType =
    LinearAlgebra::distributed::Vector<Number, MemorySpace::Default>;
  using LevelMatrixType = PortableMFVelocityOperator<dim, degree_u, degree_p>;
  using SmootherPreconditionerType = DiagonalMatrix<VectorType>;
  using SmootherType               = PreconditionChebyshev<LevelMatrixType,
                                             VectorType,
                                             SmootherPreconditionerType>;
  using MGTransferType =
    MGTransferMatrixFree<dim, Number, MemorySpace::Default>;

  const auto coarse_grid_triangulations =
    MGTransferGlobalCoarseningTools::create_geometric_coarsening_sequence(tria);

  const unsigned int max_level = coarse_grid_triangulations.size() - 1;
  // Do not go down to level 0, because this will lead to slower runtime as
  // the problem becomes very small:
  const unsigned int min_level = std::min(2U, max_level - 1);

  MGLevelObject<DoFHandler<dim>> mg_dof_handlers(min_level, max_level);
  MGLevelObject<AffineConstraints<Number>> mg_constraints(min_level, max_level);
  MGLevelObject<LevelMatrixType>           mg_matrices(min_level, max_level);

  MGLevelObject<MGTwoLevelTransferCopyToHost<dim, VectorType>> mg_transfers(
    min_level, max_level);

  std::vector<std::shared_ptr<Portable::MatrixFree<dim, Number>>>
    mf_data_levels;

  // level operators
  for (unsigned int level = min_level; level <= max_level; ++level)
    {
      auto &dof_handler = mg_dof_handlers[level];
      auto &constraint  = mg_constraints[level];

      dof_handler.reinit(*coarse_grid_triangulations[level]);
      dof_handler.distribute_dofs(fe_u);

      constraint.reinit(dof_handler.locally_owned_dofs(),
                        DoFTools::extract_locally_relevant_dofs(dof_handler));

      DoFTools::make_zero_boundary_constraints(dof_handler, constraint);
      constraint.close();

      typename Portable::MatrixFree<dim, Number>::AdditionalData
        additional_data;
      additional_data.mapping_update_flags =
        update_JxW_values | update_gradients;
      if (level == max_level)
        mf_data_levels.emplace_back(mf_data);
      else
        {
          mf_data_levels.emplace_back(
            std::make_shared<Portable::MatrixFree<dim, Number>>());
          mf_data_levels.back()->reinit(
            mapping, dof_handler, constraint, quad, additional_data);
        }

      mg_matrices[level].reinit(mf_data_levels.back());
    }

  mg::Matrix<VectorType> mg_matrix(mg_matrices);

  // transfer operator
  for (unsigned int level = min_level; level < max_level; ++level)
    mg_transfers[level + 1].reinit(mg_dof_handlers[level + 1],
                                   mg_dof_handlers[level],
                                   mg_constraints[level + 1],
                                   mg_constraints[level]);

  MGTransferType mg_transfer(mg_transfers, [&](const auto l, auto &vec) {
    mg_matrices[l].initialize_dof_vector(vec);
  });

  // smoother
  MGLevelObject<typename SmootherType::AdditionalData> smoother_data(min_level,
                                                                     max_level);

  for (unsigned int level = min_level; level <= max_level; ++level)
    {
      mg_matrices[level].compute_diagonal();
      smoother_data[level].preconditioner =
        std::make_shared<SmootherPreconditionerType>(
          *mg_matrices[level].get_matrix_diagonal_inverse());
      smoother_data[level].smoothing_range     = 20;
      smoother_data[level].degree              = 5;
      smoother_data[level].eig_cg_n_iterations = 20;
      smoother_data[level].constraints.copy_from(mg_constraints[level]);
    }

  MGSmootherPrecondition<LevelMatrixType, SmootherType, VectorType> mg_smoother;
  mg_smoother.initialize(mg_matrices, smoother_data);

  for (unsigned int level = min_level; level <= max_level; ++level)
    {
      VectorType vec;
      mg_matrices[level].initialize_dof_vector(vec);
      mg_smoother.smoothers[level].estimate_eigenvalues(vec);
    }

  // coarse-grid solver
  MGCoarseGridApplySmoother<VectorType> mg_coarse;
  mg_coarse.initialize(mg_smoother);

  // put everything together
  Multigrid<VectorType> mg(mg_matrix,
                           mg_coarse,
                           mg_transfer,
                           mg_smoother,
                           mg_smoother,
                           min_level,
                           max_level);

  using APreconditionerType = PreconditionMG<dim, VectorType, MGTransferType>;
  APreconditionerType preconditioner_A(dof_u, mg, mg_transfer);



  PortableMFMassOperator<dim, degree_u, degree_p> mass_operator(*mf_data.get());
  mass_operator.compute_diagonal();

  using SPreconditionerType = PreconditionChebyshev<
    PortableMFMassOperator<dim, degree_u, degree_p>,
    LinearAlgebra::distributed::Vector<double, MemorySpace::Default>>;

  SPreconditionerType preconditioner_schur;
  {
    typename SPreconditionerType::AdditionalData additional_data;
    additional_data.smoothing_range     = 15.;
    additional_data.degree              = 3;
    additional_data.eig_cg_n_iterations = 10;
    additional_data.constraints.copy_from(constraints_p);
    additional_data.preconditioner =
      mass_operator.get_matrix_diagonal_inverse();

    preconditioner_schur.initialize(mass_operator, additional_data);
  }

  using BTOperatorType = PortableMFBTOperator<dim, degree_u, degree_p>;
  BTOperatorType BT_operator(*mf_data.get());


  using BlockVectorType =
    LinearAlgebra::distributed::BlockVector<Number, MemorySpace::Default>;
  BlockSchurPreconditioner<APreconditionerType,
                           SPreconditionerType,
                           BTOperatorType,
                           BlockVectorType>
    preconditioner(preconditioner_A, preconditioner_schur, BT_operator);

  solver.solve(stokes_operator, solution, rhs, preconditioner);

  deallog << "converged in " << solver_control.last_step() << " iterations."
          << std::endl;

  solution_host.block(0).import_elements(solution.block(0),
                                         VectorOperation::insert);
  solution_host.block(1).import_elements(solution.block(1),
                                         VectorOperation::insert);

  constraints_u.distribute(solution_host.block(0));
  constraints_p.distribute(solution_host.block(1));
  solution_host.update_ghost_values();

  const QGauss<dim> quadrature_formula(degree_u + 1);

  Vector<double> cellwise_errors_ul2(tria.n_active_cells());
  Vector<double> cellwise_errors_pl2(tria.n_active_cells());

  VectorTools::integrate_difference(dof_u,
                                    solution_host.block(0),
                                    VelocitySolution<dim>(),
                                    cellwise_errors_ul2,
                                    quadrature_formula,
                                    VectorTools::L2_norm);
  VectorTools::integrate_difference(dof_p,
                                    solution_host.block(1),
                                    PressureSolution<dim>(),
                                    cellwise_errors_pl2,
                                    quadrature_formula,
                                    VectorTools::L2_norm);

  const double u_l2 = VectorTools::compute_global_error(tria,
                                                        cellwise_errors_ul2,
                                                        VectorTools::L2_norm);
  const double p_l2 = VectorTools::compute_global_error(tria,
                                                        cellwise_errors_pl2,
                                                        VectorTools::L2_norm);

  deallog << "velocity error: " << std::setprecision(2) << u_l2
          << " pressure error: " << p_l2 << std::endl;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);

  initlog();

  for (int ref = 1; ref <= 4; ++ref)
    test<2, 1>(ref);

  deallog << "OK" << std::endl;
}
