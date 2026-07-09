/* -----------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
 * Copyright (C) 2026 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Detailed license information governing the source code and contributions
 * can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
 *
 * -----------------------------------------------------------------------------
 */



#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/timer.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_gmres.h>

#include <deal.II/matrix_free/operators.h>
#include <deal.II/matrix_free/portable_fe_evaluation.h>
#include <deal.II/matrix_free/portable_matrix_free.h>

#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/multigrid.h>
#include <deal.II/multigrid/portable_mg_transfer_global_coarsening.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/vector_tools_integrate_difference.h>

namespace Step104
{
  using namespace dealii;

  // The index of the velocity and pressure DoFHandler into the vector
  // of DoFHandlers inside Portable::MatrixFree and the blocks of
  // the solution vector.
  constexpr unsigned int velocity_dof_handler_index = 0;
  constexpr unsigned int pressure_dof_handler_index = 1;

  // @sect3{Problem Definition}
  //
  // We start with the definition of the right-hand side and exact solution of
  // the manufactured solution in 2 and 3 dimensions. In 2D the exact solution
  // is given by
  // @f{align*}{
  // u &= \pi  \sin^2(\pi x)  \sin(2  \pi  y)
  // \\ v &= -\pi  \sin(2 \pi  x)  \sin^2(\pi y)
  // \\ p &= \cos(\pi  x)  \cos(\pi  y)
  // @f}
  // and in 3D by:
  // @f{align*}{
  // u &= \pi  \sin^2(\pi x)(\sin^2(\pi z)\sin(2\pi y)-\sin^2(\pi y)\sin(2\pi
  // z))
  // \\ v &= \pi  \sin^2(\pi y)  (\sin^2(\pi x)  \sin(2\pi z) - \sin^2(\pi
  // z)\sin(2\pi x)) \\ w &= \pi \sin^2(\pi z)  (\sin^2(\pi y)  \sin(2\pi x) -
  // \sin^2(\pi x)\sin(2\pi y))
  // \\ p &= \cos(\pi x)\cos(\pi y)\cos(\pi z)
  // @f}
  //
  // The following classes define these in code.
  template <int dim, typename Number>
  class VelocityRightHandSide : public Function<dim, Number>
  {
  public:
    VelocityRightHandSide()
      : Function<dim, Number>(dim)
    {}

    virtual Number value(const Point<dim>  &p,
                         const unsigned int component = 0) const override;
  };

  template <int dim, typename Number>
  Number
  VelocityRightHandSide<dim, Number>::value(const Point<dim>  &p,
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



  template <int dim, typename Number>
  class VelocitySolution : public Function<dim, Number>
  {
  public:
    VelocitySolution()
      : Function<dim, Number>(dim)
    {}

    virtual Number value(const Point<dim>  &p,
                         const unsigned int component = 0) const override;
  };



  template <int dim, typename Number>
  Number
  VelocitySolution<dim, Number>::value(const Point<dim>  &p,
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



  template <int dim, typename Number>
  class PressureSolution : public Function<dim, Number>
  {
  public:
    PressureSolution()
      : Function<dim, Number>()
    {}

    virtual Number value(const Point<dim> &p,
                         const unsigned int /*component*/ = 0) const override
    {
      const double pi = numbers::PI;
      if constexpr (dim == 2)
        return std::cos(pi * p[0]) * std::cos(pi * p[1]);
      else
        return std::cos(pi * p[0]) * std::cos(pi * p[1]) * std::cos(pi * p[2]);
    }
  };


  // @sect3{The velocity operator}
  //
  // The matrix-free operator for the velocity block $A$
  // given by $(\nabla u,\nabla v)$ is defined by the class
  // PortableMFVelocityOperator. It uses
  // the class VelocityCellOperator, which is evaluated in parallel
  // on each cell. On each cell, we define the action at each
  // quadrature point with the small helper class VelocityOperatorQuad
  // with operator().

  template <int dim, int fe_degree, typename Number>
  class VelocityOperatorQuad
  {
  public:
    DEAL_II_HOST_DEVICE void operator()(
      Portable::FEEvaluation<dim, fe_degree, fe_degree + 1, dim, Number>
               *fe_eval,
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
  class VelocityCellOperator
  {
  public:
    static const unsigned int n_q_points =
      dealii::Utilities::pow(n_q_points_1d, dim);

    DEAL_II_HOST_DEVICE void
    operator()(const typename Portable::MatrixFree<dim, Number>::Data *data,
               const Portable::DeviceVector<Number>                   &src,
               Portable::DeviceVector<Number> &dst) const
    {
      Portable::FEEvaluation<dim, degree_u, n_q_points_1d, dim, Number> fe_u(
        data, velocity_dof_handler_index);

      fe_u.read_dof_values(src);
      fe_u.evaluate(EvaluationFlags::gradients);

      VelocityOperatorQuad<dim, degree_u, Number> quad_operation;

      data->for_each_quad_point(
        [&](const int q_point) { quad_operation(&fe_u, q_point); });

      fe_u.integrate(EvaluationFlags::gradients);
      fe_u.distribute_local_to_global(dst);
    }
  };


  // This class finally provides the matrix-free operator for the velocity
  // block. Note that we also compute the inverse diagonal of the operator,
  // which is used in the Chebyshev smoother when we approximate A^{-1} with a
  // GMG v-cycle. We note that Tvmult() is not implemented because it is not
  // required for the smoother and $A$ is symmetric anyway.
  template <int dim,
            int degree_u,
            int degree_p,
            typename Number = double,
            typename VectorType =
              LinearAlgebra::distributed::Vector<Number, MemorySpace::Default>,
            int n_q_points_1d = degree_u + 1>
  class PortableMFVelocityOperator : public EnableObserverPointer
  {
  public:
    PortableMFVelocityOperator() = default;

    PortableMFVelocityOperator(
      std::shared_ptr<Portable::MatrixFree<dim, Number>> data_in)
      : data(data_in)
    {}

    void reinit(std::shared_ptr<Portable::MatrixFree<dim, Number>> data_in)
    {
      data = data_in;
    }

    void initialize_dof_vector(VectorType &vec) const
    {
      data->initialize_dof_vector(vec, velocity_dof_handler_index);
    }

    types::global_dof_index m() const
    {
      return data->get_vector_partitioner(velocity_dof_handler_index)->size();
    }

    std::shared_ptr<DiagonalMatrix<
      LinearAlgebra::distributed::Vector<Number, MemorySpace::Default>>>
    get_matrix_diagonal_inverse() const
    {
      return inverse_diagonal_entries;
    }

    double el(const types::global_dof_index row,
              const types::global_dof_index col) const
    {
      (void)col;
      Assert(row == col, ExcNotImplemented());
      Assert(inverse_diagonal_entries.get() != nullptr &&
               inverse_diagonal_entries->m() > 0,
             ExcNotInitialized());
      return 1.0 / (*inverse_diagonal_entries)(row, row);
    }

    void vmult(VectorType &dst, const VectorType &src) const
    {
      dst = static_cast<Number>(0.);
      VelocityCellOperator<dim, degree_u, degree_p, Number, n_q_points_1d>
        velocity_operator;
      data->cell_loop(velocity_operator, src, dst);

      data->copy_constrained_values(src, dst, velocity_dof_handler_index);
    }

    void Tvmult(VectorType & /* dst */, const VectorType & /* src */) const
    {
      AssertThrow(false, ExcNotImplemented());
    }

    void compute_diagonal()
    {
      Assert(data.get() != nullptr, ExcNotInitialized());

      this->inverse_diagonal_entries =
        std::make_shared<DiagonalMatrix<VectorType>>();
      VectorType &inverse_diagonal =
        this->inverse_diagonal_entries->get_vector();
      data->initialize_dof_vector(inverse_diagonal, velocity_dof_handler_index);

      VelocityOperatorQuad<dim, degree_u, Number> velocity_operator_quad;

      MatrixFreeTools::
        compute_diagonal<dim, degree_u, degree_u + 1, dim, Number>(
          *data.get(),
          inverse_diagonal,
          velocity_operator_quad,
          EvaluationFlags::gradients,
          EvaluationFlags::gradients,
          velocity_dof_handler_index);

      Number *raw_diagonal = inverse_diagonal.get_values();

      Kokkos::parallel_for(
        "invert A diagonal",
        inverse_diagonal.locally_owned_size(),
        KOKKOS_LAMBDA(int i) {
          Assert(raw_diagonal[i] > 0.,
                 ExcMessage("Diagonal entries of a positive definite operator "
                            "should be positive"));
          raw_diagonal[i] = 1. / raw_diagonal[i];
        });
    }

  private:
    std::shared_ptr<Portable::MatrixFree<dim, Number>> data;
    std::shared_ptr<DiagonalMatrix<VectorType>>        inverse_diagonal_entries;
  };



  //@sect3{The Schur complement operator}

  // The preconditioner requires a Schur complement
  // approximation, which is here given by a mass matrix in the
  // pressure space $(p,q)$. This is implemented in a very similar way
  // to the velocity block above. A notable difference is that
  // we select the pressure by passing a dof_handler_index of 1
  // instead of 0 to the FEEvaluation class.
  template <int dim,
            int degree_u,
            int degree_p,
            typename Number,
            int n_q_points_1d>
  class MassOperatorQuad
  {
  public:
    DEAL_II_HOST_DEVICE void operator()(
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
               Portable::DeviceVector<Number> &dst) const;
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
    Portable::FEEvaluation<dim, degree_p, n_q_points_1d, 1, Number> fe_p(
      data, pressure_dof_handler_index);

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
              LinearAlgebra::distributed::Vector<Number, MemorySpace::Default>,
            int n_q_points_1d = degree_u + 1>
  class PortableMFMassOperator : public EnableObserverPointer
  {
  public:
    PortableMFMassOperator(
      const std::shared_ptr<Portable::MatrixFree<dim, Number>> &data_in)
      : data(data_in)
    {}

    types::global_dof_index m() const
    {
      return data->get_vector_partitioner(pressure_dof_handler_index)->size();
    }

    Number el(const types::global_dof_index row,
              const types::global_dof_index col) const
    {
      (void)col;
      Assert(row == col, ExcNotImplemented());
      Assert(inverse_diagonal_entries.get() != nullptr &&
               inverse_diagonal_entries->m() > 0,
             ExcNotInitialized());
      return 1.0 / (*inverse_diagonal_entries)(row, row);
    }

    void vmult(VectorType &dst, const VectorType &src) const
    {
      dst = static_cast<Number>(0.);
      MassCellOperator<dim, degree_u, degree_p, Number, n_q_points_1d>
        mass_operator;
      data->cell_loop(mass_operator, src, dst);

      data->copy_constrained_values(src, dst, pressure_dof_handler_index);
    }

    std::shared_ptr<DiagonalMatrix<VectorType>>
    get_matrix_diagonal_inverse() const
    {
      return inverse_diagonal_entries;
    }

    void compute_diagonal()
    {
      this->inverse_diagonal_entries =
        std::make_shared<DiagonalMatrix<VectorType>>();
      VectorType &inverse_diagonal =
        this->inverse_diagonal_entries->get_vector();

      MassOperatorQuad<dim, degree_u, degree_p, Number, n_q_points_1d>
        quad_operation;

      MatrixFreeTools::
        compute_diagonal<dim, degree_p, n_q_points_1d, 1, Number>(
          *data.get(),
          inverse_diagonal,
          quad_operation,
          EvaluationFlags::values,
          EvaluationFlags::values,
          pressure_dof_handler_index);

      Number *raw_diagonal = inverse_diagonal.get_values();

      Kokkos::parallel_for(
        "invert Mass diagonal",
        inverse_diagonal.locally_owned_size(),
        KOKKOS_LAMBDA(int i) {
          Assert(raw_diagonal[i] > 0.,
                 ExcMessage("Diagonal entries of a positive definite operator "
                            "should be positive"));
          raw_diagonal[i] = 1. / raw_diagonal[i];
        });
    }

  private:
    std::shared_ptr<Portable::MatrixFree<dim, Number>> data;
    std::shared_ptr<DiagonalMatrix<VectorType>>        inverse_diagonal_entries;
  };



  //@sect3{The Stokes operator}

  // The following set of classes provides the whole Stokes operator
  // @f{eqnarray*}{
  // \begin{bmatrix} A & B^T \\ B & 0 \end{bmatrix}.
  // @f}
  // While structured in a similar way (class PortableMFStokesOperator uses
  // StokesCellOperator and a lambda function for the action at each
  // quadrature point), we now operate on Portable::DeviceBlockVector
  // and use two Portable::FEEvaluation objects, one for the velocity
  // and one for the pressure.
  //
  // Note that we don't need support for computing the diagonal, as this
  // is not needed in the block preconditioner.
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
               Portable::DeviceBlockVector<Number> &dst) const
    {
      Portable::FEEvaluation<dim, degree_u, n_q_points_1d, dim, Number> fe_u(
        data, velocity_dof_handler_index);
      Portable::FEEvaluation<dim, degree_p, n_q_points_1d, 1, Number> fe_p(
        data, pressure_dof_handler_index);

      fe_u.read_dof_values(src.block(0));
      fe_p.read_dof_values(src.block(1));
      fe_u.evaluate(EvaluationFlags::gradients);
      fe_p.evaluate(EvaluationFlags::values);

      data->for_each_quad_point([&](const int &q_point) {
        const Tensor<2, dim, Number> gradient_u = fe_u.get_gradient(q_point);
        const Number                 pressure_value = fe_p.get_value(q_point);

        Tensor<2, dim, Number> velocity_term = gradient_u;
        for (unsigned int d = 0; d < dim; ++d)
          velocity_term[d][d] -= pressure_value;
        fe_u.submit_gradient(velocity_term, q_point);

        const Number pressure_term = trace(gradient_u);
        fe_p.submit_value(pressure_term, q_point);
      });

      fe_u.integrate(EvaluationFlags::gradients);
      fe_p.integrate(EvaluationFlags::values);
      fe_u.distribute_local_to_global(dst.block(0));
      fe_p.distribute_local_to_global(dst.block(1));
    }
  };


  template <
    int dim,
    int degree_u,
    int degree_p,
    typename Number = double,
    typename VectorType =
      LinearAlgebra::distributed::BlockVector<Number, MemorySpace::Default>,
    int n_q_points_1d = degree_u + 1>
  class PortableMFStokesOperator
  {
  public:
    PortableMFStokesOperator(
      const std::shared_ptr<Portable::MatrixFree<dim, Number>> &data_in)
      : data(data_in)
    {}

    void vmult(VectorType &dst, const VectorType &src) const
    {
      dst = static_cast<Number>(0.);
      StokesCellOperator<dim, degree_u, degree_p, Number, n_q_points_1d>
        stokes_operator;
      data->cell_loop(stokes_operator, src, dst);

      data->copy_constrained_values(src, dst);
    }

  private:
    std::shared_ptr<Portable::MatrixFree<dim, Number>> data;
  };



  // @sect3{The BT operator}
  // The last ingredient of the block preconditioner is the action
  // of the block $B^T$ given by $-(p,\nabla \cdot v)$. The operator reads
  // pressure values $p$ and produces
  // a velocity (you can think of it as a rectangular matrix block). Therefore,
  // the implementation is similar to the Stokes operator in that we work with
  // two Portable::FEEvaluation objects and operate on
  // Portable::DeviceBlockVector.
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
               Portable::DeviceBlockVector<Number> &dst) const
    {
      Portable::FEEvaluation<dim, degree_u, n_q_points_1d, dim, Number> fe_u(
        data, velocity_dof_handler_index);
      Portable::FEEvaluation<dim, degree_p, n_q_points_1d, 1, Number> fe_p(
        data, pressure_dof_handler_index);

      fe_p.read_dof_values(src.block(1));
      fe_p.evaluate(EvaluationFlags::values);

      data->for_each_quad_point([&](const int &q_point) {
        const Number pressure_value = fe_p.get_value(q_point);
        fe_u.submit_divergence(-pressure_value, q_point);
      });

      fe_u.integrate(EvaluationFlags::gradients);
      fe_u.distribute_local_to_global(dst.block(0));
    }
  };



  template <
    int dim,
    int degree_u,
    int degree_p,
    typename Number = double,
    typename VectorType =
      LinearAlgebra::distributed::BlockVector<Number, MemorySpace::Default>,
    int n_q_points_1d = degree_u + 1>
  class PortableMFBTOperator
  {
  public:
    PortableMFBTOperator(
      const std::shared_ptr<Portable::MatrixFree<dim, Number>> &data_in)
      : data(data_in)
    {}



    void vmult(VectorType &dst, const VectorType &src) const
    {
      dst = static_cast<Number>(0.);
      BTCellOperator<dim, degree_u, degree_p, Number, n_q_points_1d>
        cell_operator;
      data->cell_loop(cell_operator, src, dst);

      // Instead of copying constrained values, zero them out. The BT operator
      // does not receive an input velocity to copy values from and zeroing out
      // is the correct operation for boundary and hanging nodes for an update
      // of the velocity:
      data->set_constrained_values(0.0,
                                   dst.block(0),
                                   velocity_dof_handler_index);
    }

  private:
    std::shared_ptr<Portable::MatrixFree<dim, Number>> data;
  };



  // @sect3{The Preconditioner <code>BlockSchurPreconditioner</code>}
  //
  // The following class implements the block preconditioner. The class
  // takes the types of the operators for $A^{-1}$, $S^{-1}$, and
  // $B^T$ as template arguments. This is the same preconditioner used
  // in step-32 and step-56.
  //
  // We keep a temporary vector `tmp` inside this class to avoid reallocating
  // memory in every vmult call. It has to be mutable because it is used in the
  // vmult method that is declared as const.
  template <class AInvOperator,
            class SInvOperator,
            class BTOperator,
            class VectorType>
  class BlockSchurPreconditioner : public EnableObserverPointer
  {
  public:
    /**
     * @brief Constructor
     * @param A_inverse_operator Approximation of the inverse of the velocity block.
     * @param S_inverse_operator Approximation of the inverse Schur complement.
     * @param BT_operator Operator for the B^T block of the Stokes system. Note
     * that this operator is exactly applied while @p A_inverse_operator and @p S_inverse_operator
     * are approximations for the purpose of the block preconditioner.
     */
    BlockSchurPreconditioner(const AInvOperator &A_inverse_operator,
                             const SInvOperator &S_inverse_operator,
                             const BTOperator   &BT_operator);

    /**
     * Matrix-vector product with this preconditioner object.
     */
    void vmult(VectorType &dst, const VectorType &src) const;

  private:
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
    // Allocate the temporary vector on first use. Its content doesn't matter,
    // as it will be overwritten in the vmult() call.
    if (tmp.size() == 0)
      tmp.reinit(src);

    // First apply the Schur Complement inverse operator: dst.p = S^-1 * src.p
    {
      S_inverse_operator.vmult(dst.block(1), src.block(1));
      dst.block(1) *= -1.0;
    }

    // Apply the top right block: tmp.u = -B^T * dst.p + src.u
    {
      BT_operator.vmult(tmp, dst);
      tmp.block(0) *= -1.0;
      tmp.block(0) += src.block(0);
    }

    // Finally the velocity block:
    A_inverse_operator.vmult(dst.block(0), tmp.block(0));
  }



  // @sect3{The main class <code>StokesProblem</code>}
  //
  // The remaining part of this tutorial is the StokesProblem class
  // that puts everything together.
  template <int dim, int degree_p, typename Number = double>
  class StokesProblem
  {
  public:
    static constexpr unsigned int degree_u = degree_p + 1;

    StokesProblem();

    void run();

    using VectorType =
      LinearAlgebra::distributed::Vector<Number, MemorySpace::Default>;
    using BlockVectorType =
      LinearAlgebra::distributed::BlockVector<Number, MemorySpace::Default>;

  private:
    void setup_dofs();

    void solve();

    void postprocess();

    parallel::distributed::Triangulation<dim> tria;

    MappingQ<dim> mapping;

    FESystem<dim> fe_u;
    FE_Q<dim>     fe_p;

    DoFHandler<dim> dof_u;
    DoFHandler<dim> dof_p;

    AffineConstraints<Number> constraints_u;
    AffineConstraints<Number> constraints_p;

    std::shared_ptr<Portable::MatrixFree<dim, Number>> mf_data;
    BlockVectorType                                    solution;
    BlockVectorType                                    rhs;
    ConditionalOStream                                 pcout;
  };


  template <int dim, int degree_p, typename Number>
  StokesProblem<dim, degree_p, Number>::StokesProblem()
    : tria(MPI_COMM_WORLD)
    , mapping(1)
    , fe_u(FE_Q<dim>(degree_p + 1), dim)
    , fe_p(degree_p)
    , dof_u(tria)
    , dof_p(tria)
    , pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
  {}

  // The setup_dofs() function distributes the two separate DoFHandlers for the
  // velocity and pressure before initializing the MatrixFree object with an
  // std::vector of both of them. We can later refer to the velocity using
  // DoFHandler index 0 and the pressure using DoFHandler index 1.
  template <int dim, int degree_p, typename Number>
  void StokesProblem<dim, degree_p, Number>::setup_dofs()
  {
    dof_u.distribute_dofs(fe_u);
    dof_p.distribute_dofs(fe_p);

    const IndexSet &owned_set_u = dof_u.locally_owned_dofs();
    const IndexSet  relevant_set_u =
      DoFTools::extract_locally_relevant_dofs(dof_u);
    constraints_u.reinit(owned_set_u, relevant_set_u);
    DoFTools::make_hanging_node_constraints(dof_u, constraints_u);
    VectorTools::interpolate_boundary_values(
      dof_u, 0, Functions::ZeroFunction<dim, Number>(dim), constraints_u);
    constraints_u.close();

    const IndexSet &owned_set_p = dof_p.locally_owned_dofs();
    const IndexSet  relevant_set_p =
      DoFTools::extract_locally_relevant_dofs(dof_p);

    constraints_p.reinit(owned_set_p, relevant_set_p);
    DoFTools::make_hanging_node_constraints(dof_p, constraints_p);
    constraints_p.close();

    std::vector<const DoFHandler<dim> *> dof_handlers = {&dof_u, &dof_p};
    std::vector<const AffineConstraints<Number> *> constraints = {
      &constraints_u, &constraints_p};

    mf_data = std::make_shared<Portable::MatrixFree<dim, Number>>();

    const QGauss<1> quad(degree_p + 2);
    typename Portable::MatrixFree<dim, Number>::AdditionalData additional_data;
    additional_data.mapping_update_flags = update_values | update_gradients;
    mf_data->reinit(mapping, dof_handlers, constraints, quad, additional_data);

    {
      // create the right hand side on the host and move to device:

      LinearAlgebra::distributed::BlockVector<Number, MemorySpace::Host>
        rhs_host;
      mf_data->initialize_dof_vector(rhs_host);

      VectorTools::create_right_hand_side(mapping,
                                          dof_u,
                                          QGauss<dim>(degree_u + 2),
                                          VelocityRightHandSide<dim, Number>(),
                                          rhs_host.block(0),
                                          constraints_u);

      mf_data->initialize_dof_vector(rhs);
      rhs.block(0).import_elements(rhs_host.block(0), VectorOperation::insert);
      rhs.block(1).import_elements(rhs_host.block(1), VectorOperation::insert);
    }
  }


  // In the solve() function we set up the preconditioner and run the GMRES
  // solver. For this, we construct the multigrid
  // hierarchy for the GMG v-cycle with a Chebyshev iteration around the
  // point-Jacobi scheme, i.e., the inverse of the diagonal of $A$, to
  // approximate the action of $A^{-1}$.
  // We approximate the Schur Complement with a Chebyshev iteration
  // applied to the pressure mass matrix (without multigrid).
  template <int dim, int degree_p, typename Number>
  void StokesProblem<dim, degree_p, Number>::solve()
  {
    PortableMFStokesOperator<dim, degree_u, degree_p, Number> stokes_operator(
      mf_data);

    mf_data->initialize_dof_vector(solution);

    {
      dealii::Timer t(tria.get_mpi_communicator());
      stokes_operator.vmult(solution, rhs);
      const double time          = t.wall_time();
      const double dofs_p_second = static_cast<double>(solution.size()) / time;
      pcout << "Stokes operator: " << time << " s, DoFs/s: " << dofs_p_second
            << std::endl;
      solution = 0.0;
    }

    SolverControl solver_control(1000, 1e-8 * rhs.l2_norm());

    SolverGMRES<BlockVectorType> solver(
      solver_control,
      typename SolverGMRES<BlockVectorType>::AdditionalData(50, true));

    using LevelMatrixType =
      PortableMFVelocityOperator<dim, degree_u, degree_p, Number>;
    using SmootherPreconditionerType = DiagonalMatrix<VectorType>;
    using SmootherType               = PreconditionChebyshev<LevelMatrixType,
                                               VectorType,
                                               SmootherPreconditionerType>;
    using MGTransferType =
      MGTransferMatrixFree<dim, Number, MemorySpace::Default>;

    const auto coarse_grid_triangulations =
      MGTransferGlobalCoarseningTools::create_geometric_coarsening_sequence(
        tria);

    const unsigned int max_level = coarse_grid_triangulations.size() - 1;
    // Do not go down to level 0, because this will lead to slower runtime as
    // the problem becomes very small:
    const unsigned int min_level = std::min(3U, max_level - 1);

    MGLevelObject<DoFHandler<dim>> mg_dof_handlers(min_level, max_level);
    MGLevelObject<AffineConstraints<Number>> mg_constraints(min_level,
                                                            max_level);
    MGLevelObject<LevelMatrixType>           mg_matrices(min_level, max_level);

    MGLevelObject<Portable::MGTwoLevelTransfer<dim, VectorType>> mg_transfers(
      min_level, max_level);

    std::vector<std::shared_ptr<Portable::MatrixFree<dim, Number>>>
      mf_data_levels;

    // Prepare the operators and data structures on all levels of the multigrid
    // hierarchy
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
          // On the finest level we can reuse the MatrixFree object from the
          // Stokes operator. This way we can solve significantly larger
          // problems before we run out of device memory.
          mf_data_levels.emplace_back(mf_data);
        else
          {
            const QGauss<1> quad(degree_p + 2);
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
      mg_transfers[level + 1].reinit_geometric_transfer(
        mg_dof_handlers[level + 1],
        mg_dof_handlers[level],
        mg_constraints[level + 1],
        mg_constraints[level]);

    MGTransferType mg_transfer(mg_transfers, [&](const auto l, auto &vec) {
      mg_matrices[l].initialize_dof_vector(vec);
    });

    // smoother
    MGLevelObject<typename SmootherType::AdditionalData> smoother_data(
      min_level, max_level);

    for (unsigned int level = min_level; level <= max_level; ++level)
      {
        mg_matrices[level].compute_diagonal();
        smoother_data[level].preconditioner =
          std::make_shared<SmootherPreconditionerType>(
            *mg_matrices[level].get_matrix_diagonal_inverse());
        smoother_data[level].constraints.copy_from(mg_constraints[level]);

        if (level == min_level)
          {
            // Use the Chebyshev iteration as an (approximate) solver on the
            // coarsest level. In this mode @p smoothing_range is a relative
            // target tolerance and must be strictly less than one; the number
            // of iterations is then chosen automatically by setting
            // @p degree to numbers::invalid_unsigned_int. We also use more
            // CG iterations for the eigenvalue estimate because when
            // @p min_level > 0, the coarse problem can still be reasonably
            // large and badly conditioned.
            smoother_data[level].smoothing_range = 1e-3;
            smoother_data[level].degree = numbers::invalid_unsigned_int;
            smoother_data[level].eig_cg_n_iterations = 40;
          }
        else
          {
            // These values are chosen by experimentation for the problem at
            // hand. We chose the smoothing range first. A good value will allow
            // the smoother to effectively separate large and small scale
            // oscillations in the residual and as such improve the convergence
            // of the Chebyshev iteration and the multigrid method. Finally, the
            // degree is chosen to minimize total runtime (a larger value
            // increases the cost but improves the outer number of GMRES
            // iterations).
            smoother_data[level].smoothing_range     = 5;
            smoother_data[level].degree              = 4;
            smoother_data[level].eig_cg_n_iterations = 20;
          }
      }

    MGSmootherPrecondition<LevelMatrixType, SmootherType, VectorType>
      mg_smoother;
    mg_smoother.initialize(mg_matrices, smoother_data);

    // Estimate and print the eigenvalue spectrum of the velocity block on each
    // level. This spectrum is later used by the Chebyshev iteration.
    pcout << "GMG velocity block smoothers:" << std::endl;
    for (unsigned int level = min_level; level <= max_level; ++level)
      {
        VectorType vec;
        mg_matrices[level].initialize_dof_vector(vec);
        auto eigenvalue_info =
          mg_smoother.smoothers[level].estimate_eigenvalues(vec);
        pcout << "    level: " << level << " n_dofs: " << vec.size()
              << ", eigenvalue spectrum: [ "
              << eigenvalue_info.min_eigenvalue_estimate << ", "
              << eigenvalue_info.max_eigenvalue_estimate << " ]" << std::endl;
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


    dealii::Timer timer_smoother;
    dealii::Timer timer_transfer;
    dealii::Timer timer_coarse;
    dealii::Timer timer_residual;
    {
      timer_smoother.reset();
      timer_transfer.reset();
      timer_coarse.reset();
      timer_residual.reset();

      auto make_timer_lambda = [&](dealii::Timer &timer) {
        return [&](const bool before, const unsigned int /*level*/) {
          if (before)
            timer.start();
          else
            timer.stop();
        };
      };
      mg.connect_pre_smoother_step(make_timer_lambda(timer_smoother));
      mg.connect_post_smoother_step(make_timer_lambda(timer_smoother));
      mg.connect_residual_step(make_timer_lambda(timer_residual));
      mg.connect_restriction(make_timer_lambda(timer_transfer));
      mg.connect_prolongation(make_timer_lambda(timer_transfer));
      mg.connect_coarse_solve(make_timer_lambda(timer_coarse));
    }

    using APreconditionerType = PreconditionMG<dim, VectorType, MGTransferType>;
    APreconditionerType preconditioner_A(dof_u, mg, mg_transfer);


    PortableMFMassOperator<dim, degree_u, degree_p, Number> mass_operator(
      mf_data);
    mass_operator.compute_diagonal();

    using SPreconditionerType = PreconditionChebyshev<
      PortableMFMassOperator<dim, degree_u, degree_p, Number>,
      VectorType>;

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

    using BTOperatorType =
      PortableMFBTOperator<dim, degree_u, degree_p, Number>;
    BTOperatorType BT_operator(mf_data);

    BlockSchurPreconditioner<APreconditionerType,
                             SPreconditionerType,
                             BTOperatorType,
                             BlockVectorType>
      preconditioner(preconditioner_A, preconditioner_schur, BT_operator);

    dealii::Timer t(tria.get_mpi_communicator());
    solver.solve(stokes_operator, solution, rhs, preconditioner);
    t.stop();

    pcout << "Solver converged in " << solver_control.last_step()
          << " iterations in " << t.wall_time() << " seconds" << std::endl;

    pcout << "Velocity block GMG timings:"
          << "\n    smoother: " << timer_smoother.wall_time()
          << " s\n    transfer: " << timer_transfer.wall_time()
          << " s\n    coarse  : " << timer_coarse.wall_time()
          << " s\n    residual: " << timer_residual.wall_time() << " s"
          << std::endl;
  }



  // The postprocess() function moves the solution to host memory
  // and integrates the difference to the manufactured solution to
  // compute errors.
  template <int dim, int degree_p, typename Number>
  void StokesProblem<dim, degree_p, Number>::postprocess()
  {
    LinearAlgebra::distributed::BlockVector<Number, MemorySpace::Host>
      solution_host;
    mf_data->initialize_dof_vector(solution_host);

    solution_host.block(0).import_elements(solution.block(0),
                                           VectorOperation::insert);
    solution_host.block(1).import_elements(solution.block(1),
                                           VectorOperation::insert);

    constraints_u.distribute(solution_host.block(0));
    constraints_p.distribute(solution_host.block(1));
    solution_host.update_ghost_values();
    const double mean_pressure = VectorTools::compute_mean_value(
      dof_p, QGauss<dim>(degree_p + 2), solution_host.block(1), 0);
    solution_host.block(1).add(-mean_pressure);

    const QGauss<dim> quadrature_formula(degree_u + 1);

    Vector<double> cellwise_errors_ul2(tria.n_active_cells());
    Vector<double> cellwise_errors_pl2(tria.n_active_cells());

    VectorTools::integrate_difference(dof_u,
                                      solution_host.block(0),
                                      VelocitySolution<dim, Number>(),
                                      cellwise_errors_ul2,
                                      quadrature_formula,
                                      VectorTools::L2_norm);
    VectorTools::integrate_difference(dof_p,
                                      solution_host.block(1),
                                      PressureSolution<dim, Number>(),
                                      cellwise_errors_pl2,
                                      quadrature_formula,
                                      VectorTools::L2_norm);

    const double u_l2 = VectorTools::compute_global_error(tria,
                                                          cellwise_errors_ul2,
                                                          VectorTools::L2_norm);
    const double p_l2 = VectorTools::compute_global_error(tria,
                                                          cellwise_errors_pl2,
                                                          VectorTools::L2_norm);

    pcout << "velocity error: " << u_l2 << " pressure error: " << p_l2
          << std::endl;
  }



  // The run() function prints some statistics and then performs a familiar
  // refinement loop.
  template <int dim, int degree_p, typename Number>
  void StokesProblem<dim, degree_p, Number>::run()
  {
    pcout << std::setprecision(10);
    pcout << "Running on " << Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)
          << " MPI ranks (with " << MultithreadInfo::n_threads()
          << " threads each) in ";
    if constexpr (running_in_debug_mode())
      pcout << "DEBUG mode";
    else
      pcout << "RELEASE mode";

    pcout << "\nKokkos execution space: "
          << Kokkos::DefaultExecutionSpace::name();
    pcout << '\n'
          << "dim: " << dim << '\n'
          << "Element: Q" << degree_u << "-Q" << degree_p << std::endl;

    unsigned int n_refinements = 10;

    for (unsigned int i = 0; i < n_refinements; ++i)
      {
        if (i == 0)
          {
            GridGenerator::hyper_cube(tria);
            tria.refine_global(2);
          }
        else
          {
            tria.refine_global(1);
          }
        setup_dofs();

        pcout << "\nrefinement: " << i
              << ", n_dofs: " << dof_u.n_dofs() + dof_p.n_dofs() << " = "
              << dof_u.n_dofs() << " + " << dof_p.n_dofs() << std::endl;

        solve();
        postprocess();
      }
  }
} // namespace Step104



// @sect3{The <code>main()</code> function}
//
// The only interesting bits here are the template arguments that
// specify dimension and polynomial degree to be used.
int main(int argc, char **argv)
{
  using namespace Step104;
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);

  const unsigned int                   dim      = 3;
  const unsigned int                   degree_p = 1;
  StokesProblem<dim, degree_p, double> problem;
  problem.run();
}
