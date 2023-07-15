// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2022 by the deal.II authors
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


// test for correctness of matrix free implementation for multigrid stokes


#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>

#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>
#include <deal.II/multigrid/multigrid.h>

#include <deal.II/numerics/vector_tools.h>

#include <cmath>
#include <sstream>

#include "../tests.h"


unsigned int       minlevel        = 0;
const unsigned int velocity_degree = 2;

double pressure_scaling = 1.0;

namespace StokesClass
{
  class QuietException
  {};

  namespace StokesSolver
  {
    /**
     * Implement the block Schur preconditioner for the Stokes system.
     */
    template <class StokesMatrixType,
              class MassMatrixType,
              class PreconditionerA,
              class PreconditionerMp>
    class BlockSchurPreconditioner : public Subscriptor
    {
    public:
      /**
       * brief Constructor
       *
       *  S The entire Stokes matrix
       *  Spre The matrix whose blocks are used in the definition of
       *     the preconditioning of the Stokes matrix, i.e. containing
       *approximations of the A and S blocks. Mppreconditioner Preconditioner
       *object for the Schur complement, typically chosen as the mass matrix.
       *  Apreconditioner Preconditioner object for the matrix A.
       *  do_solve_A A flag indicating whether we should actually solve with
       *     the matrix $A$, or only apply one preconditioner step with it.
       *  A_block_tolerance The tolerance for the CG solver which computes
       *     the inverse of the A block.
       *  S_block_tolerance The tolerance for the CG solver which computes
       *     the inverse of the S block (Schur complement matrix).
       **/
      BlockSchurPreconditioner(const StokesMatrixType &S,
                               const MassMatrixType &  Mass,
                               const PreconditionerMp &Mppreconditioner,
                               const PreconditionerA & Apreconditioner,
                               const bool              do_solve_A,
                               const double            A_block_tolerance,
                               const double            S_block_tolerance);

      /**
       * Matrix vector product with this preconditioner object.
       */
      void
      vmult(LinearAlgebra::distributed::BlockVector<double> &      dst,
            const LinearAlgebra::distributed::BlockVector<double> &src) const;

      unsigned int
      n_iterations_A() const;
      unsigned int
      n_iterations_S() const;

    private:
      /**
       * References to the various matrix object this preconditioner works on.
       */
      const StokesMatrixType &stokes_matrix;
      const MassMatrixType &  mass_matrix;
      const PreconditionerMp &mp_preconditioner;
      const PreconditionerA & a_preconditioner;

      /**
       * Whether to actually invert the $\tilde A$ part of the preconditioner
       *matrix or to just apply a single preconditioner step with it.
       **/
      const bool           do_solve_A;
      mutable unsigned int n_iterations_A_;
      mutable unsigned int n_iterations_S_;
      const double         A_block_tolerance;
      const double         S_block_tolerance;
    };


    template <class StokesMatrixType,
              class MassMatrixType,
              class PreconditionerA,
              class PreconditionerMp>
    BlockSchurPreconditioner<StokesMatrixType,
                             MassMatrixType,
                             PreconditionerA,
                             PreconditionerMp>::
      BlockSchurPreconditioner(const StokesMatrixType &S,
                               const MassMatrixType &  Mass,
                               const PreconditionerMp &Mppreconditioner,
                               const PreconditionerA & Apreconditioner,
                               const bool              do_solve_A,
                               const double            A_block_tolerance,
                               const double            S_block_tolerance)
      : stokes_matrix(S)
      , mass_matrix(Mass)
      , mp_preconditioner(Mppreconditioner)
      , a_preconditioner(Apreconditioner)
      , do_solve_A(do_solve_A)
      , n_iterations_A_(0)
      , n_iterations_S_(0)
      , A_block_tolerance(A_block_tolerance)
      , S_block_tolerance(S_block_tolerance)
    {}

    template <class StokesMatrixType,
              class MassMatrixType,
              class PreconditionerA,
              class PreconditionerMp>
    unsigned int
    BlockSchurPreconditioner<StokesMatrixType,
                             MassMatrixType,
                             PreconditionerA,
                             PreconditionerMp>::n_iterations_A() const
    {
      return n_iterations_A_;
    }

    template <class StokesMatrixType,
              class MassMatrixType,
              class PreconditionerA,
              class PreconditionerMp>
    unsigned int
    BlockSchurPreconditioner<StokesMatrixType,
                             MassMatrixType,
                             PreconditionerA,
                             PreconditionerMp>::n_iterations_S() const
    {
      return n_iterations_S_;
    }

    template <class StokesMatrixType,
              class MassMatrixType,
              class PreconditionerA,
              class PreconditionerMp>
    void
    BlockSchurPreconditioner<StokesMatrixType,
                             MassMatrixType,
                             PreconditionerA,
                             PreconditionerMp>::
      vmult(LinearAlgebra::distributed::BlockVector<double> &      dst,
            const LinearAlgebra::distributed::BlockVector<double> &src) const
    {
      LinearAlgebra::distributed::BlockVector<double> utmp(src);

      // first solve with the bottom left block, which we have built
      // as a mass matrix with the inverse of the viscosity
      {
        SolverControl solver_control(1000,
                                     src.block(1).l2_norm() * S_block_tolerance,
                                     false,
                                     false);

        SolverCG<LinearAlgebra::distributed::Vector<double>> solver(
          solver_control);
        try
          {
            dst.block(1) = 0.0;
            solver.solve(mass_matrix,
                         dst.block(1),
                         src.block(1),
                         mp_preconditioner);
            n_iterations_S_ += solver_control.last_step();
          }
        // if the solver fails, report the error from processor 0 with some
        // additional information about its location, and throw a quiet
        // exception on all other processors
        catch (const std::exception &exc)
          {
            if (Utilities::MPI::this_mpi_process(
                  src.block(0).get_mpi_communicator()) == 0)
              AssertThrow(
                false,
                ExcMessage(
                  std::string(
                    "The iterative (bottom right) solver in BlockSchurPreconditioner::vmult "
                    "did not converge to a tolerance of " +
                    Utilities::to_string(solver_control.tolerance()) +
                    ". It reported the following error:\n\n") +
                  exc.what())) else throw QuietException();
          }
        dst.block(1) *= -1.0;
      }

      // apply the top right block
      {
        LinearAlgebra::distributed::BlockVector<double> dst_tmp(dst);
        dst_tmp.block(0) = 0.0;
        stokes_matrix.vmult(utmp, dst_tmp); // B^T
        utmp.block(0) *= -1.0;
        utmp.block(0) += src.block(0);
      }

      // now either solve with the top left block (if do_solve_A==true)
      // or just apply one preconditioner sweep (for the first few
      // iterations of our two-stage outer GMRES iteration)
      if (do_solve_A == true)
        {
          Assert(false, ExcNotImplemented());
        }
      else
        {
          a_preconditioner.vmult(dst.block(0), utmp.block(0));
          n_iterations_A_ += 1;
        }
    }
  } // namespace StokesSolver


  // Parameters for Sinker example
  double beta  = 10.0;
  double delta = 200.0;
  double omega = 0.1;

  template <int dim>
  struct Sinker
  {
    unsigned int            problem_dim;
    unsigned int            n_sinkers;
    std::vector<Point<dim>> centers;
    double                  DR_mu;
    double                  mu_min;
    double                  mu_max;
  };

  template <int dim>
  class Viscosity
  {
  public:
    Viscosity(const Sinker<dim> &sink);
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const;
    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<double> &          values,
               const unsigned int             component = 0) const;

    Sinker<dim> sinker;
  };
  template <int dim>
  Viscosity<dim>::Viscosity(const Sinker<dim> &sink)
  {
    sinker = sink;
  }
  template <int dim>
  double
  Viscosity<dim>::value(const Point<dim> &p,
                        const unsigned int /*component*/) const
  {
    double Chi = 1.0;
    for (unsigned int s = 0; s < sinker.n_sinkers; ++s)
      {
        double dist = p.distance(sinker.centers[s]);
        double temp =
          1 - std::exp(-delta * std::pow(std::max(0.0, dist - omega / 2.0), 2));
        Chi *= temp;
      }
    return (sinker.mu_max - sinker.mu_min) * (1 - Chi) + sinker.mu_min;
  }
  template <int dim>
  void
  Viscosity<dim>::value_list(const std::vector<Point<dim>> &points,
                             std::vector<double> &          values,
                             const unsigned int             component) const
  {
    Assert(values.size() == points.size(),
           ExcDimensionMismatch(values.size(), points.size()));
    Assert(component == 0, ExcIndexRange(component, 0, 1));
    const unsigned int n_points = points.size();
    for (unsigned int i = 0; i < n_points; ++i)
      values[i] = value(points[i], component);
  }

  template <int dim>
  class RightHandSide
  {
  public:
    RightHandSide(const Sinker<dim> &sink);
    Sinker<dim> sinker;
    virtual void
    vector_value(const Point<dim> &p, Vector<double> &value) const;
  };
  template <int dim>
  RightHandSide<dim>::RightHandSide(const Sinker<dim> &sink)
  {
    sinker = sink;
  }
  template <int dim>
  void
  RightHandSide<dim>::vector_value(const Point<dim> &p,
                                   Vector<double> &  values) const
  {
    double Chi = 1.0;
    for (unsigned int s = 0; s < sinker.n_sinkers; ++s)
      {
        double dist = p.distance(sinker.centers[s]);
        double temp =
          1 - std::exp(-delta * std::pow(std::max(0.0, dist - omega / 2.0), 2));
        Chi *= temp;
      }

    if (sinker.problem_dim == 2)
      {
        values[0] = 0;
        values[1] = beta * (Chi - 1.0);
        values[2] = 0;
      }
    else if (sinker.problem_dim == 3)
      {
        values[0] = 0;
        values[1] = 0;
        values[2] = beta * (Chi - 1.0);
        values[3] = 0;
      }
    return;
  }

  template <int dim>
  class ExactSolution_BoundaryValues : public Function<dim>
  {
  public:
    ExactSolution_BoundaryValues()
      : Function<dim>(dim + 1)
    {}
    virtual void
    vector_value(const Point<dim> &p, Vector<double> &value) const;
  };
  template <int dim>
  void
  ExactSolution_BoundaryValues<dim>::vector_value(const Point<dim> &p,
                                                  Vector<double> &values) const
  {
    (void)p;
    for (unsigned int i = 0; i < values.size(); ++i)
      values(i) = 0.0;
    return;
  }

  template <int dim>
  class ExactSolution_BoundaryValues_u : public Function<dim>
  {
  public:
    ExactSolution_BoundaryValues_u()
      : Function<dim>(dim)
    {}
    virtual void
    vector_value(const Point<dim> &p, Vector<double> &value) const;
  };
  template <int dim>
  void
  ExactSolution_BoundaryValues_u<dim>::vector_value(
    const Point<dim> &p,
    Vector<double> &  values) const
  {
    (void)p;
    for (unsigned int i = 0; i < values.size(); ++i)
      values(i) = 0.0;
    return;
  }



  template <int dim, int degree_v, typename number>
  class StokesOperator
    : public MatrixFreeOperators::
        Base<dim, LinearAlgebra::distributed::BlockVector<number>>
  {
  public:
    StokesOperator()
      : MatrixFreeOperators::
          Base<dim, LinearAlgebra::distributed::BlockVector<number>>()
    {}
    void
    clear();
    void
    evaluate_2_x_viscosity(const Viscosity<dim> &viscosity_function);
    virtual void
    compute_diagonal();

  private:
    virtual void
    apply_add(LinearAlgebra::distributed::BlockVector<number> &      dst,
              const LinearAlgebra::distributed::BlockVector<number> &src) const;

    void
    local_apply(const dealii::MatrixFree<dim, number> &                data,
                LinearAlgebra::distributed::BlockVector<number> &      dst,
                const LinearAlgebra::distributed::BlockVector<number> &src,
                const std::pair<unsigned int, unsigned int> &cell_range) const;

    Table<2, VectorizedArray<number>> viscosity_x_2;
  };
  template <int dim, int degree_v, typename number>
  void
  StokesOperator<dim, degree_v, number>::clear()
  {
    viscosity_x_2.reinit(0, 0);
    MatrixFreeOperators::
      Base<dim, LinearAlgebra::distributed::BlockVector<number>>::clear();
  }
  template <int dim, int degree_v, typename number>
  void
  StokesOperator<dim, degree_v, number>::evaluate_2_x_viscosity(
    const Viscosity<dim> &viscosity_function)
  {
    const unsigned int n_cells = this->data->n_cell_batches();
    FEEvaluation<dim, degree_v, degree_v + 1, dim, number> velocity(*this->data,
                                                                    0);
    viscosity_x_2.reinit(n_cells, velocity.n_q_points);
    for (unsigned int cell = 0; cell < n_cells; ++cell)
      {
        velocity.reinit(cell);
        for (unsigned int q = 0; q < velocity.n_q_points; ++q)
          {
            VectorizedArray<number> return_value =
              make_vectorized_array<number>(1.);
            for (unsigned int i = 0; i < VectorizedArray<number>::size(); ++i)
              {
                Point<dim> p;
                for (unsigned int d = 0; d < dim; ++d)
                  {
                    p(d) = velocity.quadrature_point(q)(d)[i];
                  }
                return_value[i] = 2.0 * viscosity_function.value(p);
              }
            viscosity_x_2(cell, q) = return_value;
          }
      }
  }
  template <int dim, int degree_v, typename number>
  void
  StokesOperator<dim, degree_v, number>::local_apply(
    const dealii::MatrixFree<dim, number> &                data,
    LinearAlgebra::distributed::BlockVector<number> &      dst,
    const LinearAlgebra::distributed::BlockVector<number> &src,
    const std::pair<unsigned int, unsigned int> &          cell_range) const
  {
    using vector_t = VectorizedArray<number>;
    FEEvaluation<dim, degree_v, degree_v + 1, dim, number>   velocity(data, 0);
    FEEvaluation<dim, degree_v - 1, degree_v + 1, 1, number> pressure(data, 1);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        velocity.reinit(cell);
        velocity.read_dof_values(src.block(0));
        velocity.evaluate(EvaluationFlags::gradients);
        pressure.reinit(cell);
        pressure.read_dof_values(src.block(1));
        pressure.evaluate(EvaluationFlags::values);

        for (unsigned int q = 0; q < velocity.n_q_points; ++q)
          {
            SymmetricTensor<2, dim, vector_t> sym_grad_u =
              velocity.get_symmetric_gradient(q);
            vector_t pres = pressure.get_value(q);
            vector_t div  = -trace(sym_grad_u);
            pressure.submit_value(div, q);

            sym_grad_u *= viscosity_x_2(cell, q);
            // subtract p * I
            for (unsigned int d = 0; d < dim; ++d)
              sym_grad_u[d][d] -= pres;

            velocity.submit_symmetric_gradient(sym_grad_u, q);
          }

        velocity.integrate(EvaluationFlags::gradients);
        velocity.distribute_local_to_global(dst.block(0));
        pressure.integrate(EvaluationFlags::values);
        pressure.distribute_local_to_global(dst.block(1));
      }
  }
  template <int dim, int degree_v, typename number>
  void
  StokesOperator<dim, degree_v, number>::apply_add(
    LinearAlgebra::distributed::BlockVector<number> &      dst,
    const LinearAlgebra::distributed::BlockVector<number> &src) const
  {
    MatrixFreeOperators::
      Base<dim, LinearAlgebra::distributed::BlockVector<number>>::data
        ->cell_loop(&StokesOperator::local_apply, this, dst, src);
  }
  template <int dim, int degree_v, typename number>
  void
  StokesOperator<dim, degree_v, number>::compute_diagonal()
  {
    Assert(false, ExcNotImplemented());
  }


  template <int dim, int degree_p, typename number>
  class MassMatrixOperator
    : public MatrixFreeOperators::
        Base<dim, LinearAlgebra::distributed::Vector<number>>
  {
  public:
    MassMatrixOperator()
      : MatrixFreeOperators::Base<dim,
                                  LinearAlgebra::distributed::Vector<number>>()
    {}
    void
    clear();
    void
    evaluate_1_over_viscosity(const Viscosity<dim> &viscosity_function);
    virtual void
    compute_diagonal();

  private:
    virtual void
    apply_add(LinearAlgebra::distributed::Vector<number> &      dst,
              const LinearAlgebra::distributed::Vector<number> &src) const;

    void
    local_apply(const dealii::MatrixFree<dim, number> &           data,
                LinearAlgebra::distributed::Vector<number> &      dst,
                const LinearAlgebra::distributed::Vector<number> &src,
                const std::pair<unsigned int, unsigned int> &cell_range) const;

    void
    local_compute_diagonal(
      const MatrixFree<dim, number> &              data,
      LinearAlgebra::distributed::Vector<number> & dst,
      const unsigned int &                         dummy,
      const std::pair<unsigned int, unsigned int> &cell_range) const;

    Table<2, VectorizedArray<number>> one_over_viscosity;
  };
  template <int dim, int degree_p, typename number>
  void
  MassMatrixOperator<dim, degree_p, number>::clear()
  {
    one_over_viscosity.reinit(0, 0);
    MatrixFreeOperators::Base<dim, LinearAlgebra::distributed::Vector<number>>::
      clear();
  }
  template <int dim, int degree_p, typename number>
  void
  MassMatrixOperator<dim, degree_p, number>::evaluate_1_over_viscosity(
    const Viscosity<dim> &viscosity_function)
  {
    const unsigned int n_cells = this->data->n_cell_batches();
    FEEvaluation<dim, degree_p, degree_p + 2, 1, number> pressure(*this->data,
                                                                  0);
    one_over_viscosity.reinit(n_cells, pressure.n_q_points);
    for (unsigned int cell = 0; cell < n_cells; ++cell)
      {
        pressure.reinit(cell);
        for (unsigned int q = 0; q < pressure.n_q_points; ++q)
          {
            VectorizedArray<number> return_value =
              make_vectorized_array<number>(1.);
            for (unsigned int i = 0; i < VectorizedArray<number>::size(); ++i)
              {
                Point<dim> p;
                for (unsigned int d = 0; d < dim; ++d)
                  p(d) = pressure.quadrature_point(q)(d)[i];
                return_value[i] = 1.0 / viscosity_function.value(p);
              }
            one_over_viscosity(cell, q) = return_value;
          }
      }
  }
  template <int dim, int degree_p, typename number>
  void
  MassMatrixOperator<dim, degree_p, number>::local_apply(
    const dealii::MatrixFree<dim, number> &           data,
    LinearAlgebra::distributed::Vector<number> &      dst,
    const LinearAlgebra::distributed::Vector<number> &src,
    const std::pair<unsigned int, unsigned int> &     cell_range) const
  {
    FEEvaluation<dim, degree_p, degree_p + 2, 1, number> pressure(data);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        AssertDimension(one_over_viscosity.size(0), data.n_cell_batches());
        AssertDimension(one_over_viscosity.size(1), pressure.n_q_points);

        pressure.reinit(cell);
        pressure.read_dof_values(src);
        pressure.evaluate(EvaluationFlags::values);
        for (unsigned int q = 0; q < pressure.n_q_points; ++q)
          pressure.submit_value(one_over_viscosity(cell, q) *
                                  pressure.get_value(q),
                                q);
        pressure.integrate(EvaluationFlags::values);
        pressure.distribute_local_to_global(dst);
      }
  }
  template <int dim, int degree_p, typename number>
  void
  MassMatrixOperator<dim, degree_p, number>::apply_add(
    LinearAlgebra::distributed::Vector<number> &      dst,
    const LinearAlgebra::distributed::Vector<number> &src) const
  {
    MatrixFreeOperators::Base<dim,
                              LinearAlgebra::distributed::Vector<number>>::data
      ->cell_loop(&MassMatrixOperator::local_apply, this, dst, src);
  }
  template <int dim, int degree_p, typename number>
  void
  MassMatrixOperator<dim, degree_p, number>::compute_diagonal()
  {
    this->inverse_diagonal_entries.reset(
      new DiagonalMatrix<LinearAlgebra::distributed::Vector<number>>());
    this->diagonal_entries.reset(
      new DiagonalMatrix<LinearAlgebra::distributed::Vector<number>>());

    LinearAlgebra::distributed::Vector<number> &inverse_diagonal =
      this->inverse_diagonal_entries->get_vector();
    LinearAlgebra::distributed::Vector<number> &diagonal =
      this->diagonal_entries->get_vector();

    unsigned int dummy = 0;
    this->data->initialize_dof_vector(inverse_diagonal);
    this->data->initialize_dof_vector(diagonal);

    this->data->cell_loop(&MassMatrixOperator::local_compute_diagonal,
                          this,
                          diagonal,
                          dummy);

    this->set_constrained_entries_to_one(diagonal);
    inverse_diagonal              = diagonal;
    const unsigned int local_size = inverse_diagonal.locally_owned_size();
    for (unsigned int i = 0; i < local_size; ++i)
      {
        Assert(inverse_diagonal.local_element(i) > 0.,
               ExcMessage("No diagonal entry in a positive definite operator "
                          "should be zero"));
        inverse_diagonal.local_element(i) =
          1. / inverse_diagonal.local_element(i);
      }
  }
  template <int dim, int degree_p, typename number>
  void
  MassMatrixOperator<dim, degree_p, number>::local_compute_diagonal(
    const MatrixFree<dim, number> &             data,
    LinearAlgebra::distributed::Vector<number> &dst,
    const unsigned int &,
    const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, degree_p, degree_p + 2, 1, number> pressure(data, 0);
    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        pressure.reinit(cell);
        AlignedVector<VectorizedArray<number>> diagonal(pressure.dofs_per_cell);
        for (unsigned int i = 0; i < pressure.dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < pressure.dofs_per_cell; ++j)
              pressure.begin_dof_values()[j] = VectorizedArray<number>();
            pressure.begin_dof_values()[i] = make_vectorized_array<number>(1.);

            pressure.evaluate(EvaluationFlags::values);
            for (unsigned int q = 0; q < pressure.n_q_points; ++q)
              pressure.submit_value(one_over_viscosity(cell, q) *
                                      pressure.get_value(q),
                                    q);
            pressure.integrate(EvaluationFlags::values);

            diagonal[i] = pressure.begin_dof_values()[i];
          }

        for (unsigned int i = 0; i < pressure.dofs_per_cell; ++i)
          pressure.begin_dof_values()[i] = diagonal[i];
        pressure.distribute_local_to_global(dst);
      }
  }


  template <int dim, int degree_v, typename number>
  class ABlockOperator : public MatrixFreeOperators::
                           Base<dim, LinearAlgebra::distributed::Vector<number>>
  {
  public:
    ABlockOperator()
      : MatrixFreeOperators::Base<dim,
                                  LinearAlgebra::distributed::Vector<number>>()
    {}
    void
    clear();
    void
    evaluate_2_x_viscosity(const Viscosity<dim> &viscosity_function);
    virtual void
    compute_diagonal();

  private:
    virtual void
    apply_add(LinearAlgebra::distributed::Vector<number> &      dst,
              const LinearAlgebra::distributed::Vector<number> &src) const;

    void
    local_apply(const dealii::MatrixFree<dim, number> &           data,
                LinearAlgebra::distributed::Vector<number> &      dst,
                const LinearAlgebra::distributed::Vector<number> &src,
                const std::pair<unsigned int, unsigned int> &cell_range) const;

    void
    local_compute_diagonal(
      const MatrixFree<dim, number> &              data,
      LinearAlgebra::distributed::Vector<number> & dst,
      const unsigned int &                         dummy,
      const std::pair<unsigned int, unsigned int> &cell_range) const;

    Table<2, VectorizedArray<number>> viscosity_x_2;
  };
  template <int dim, int degree_v, typename number>
  void
  ABlockOperator<dim, degree_v, number>::clear()
  {
    viscosity_x_2.reinit(0, 0);
    MatrixFreeOperators::Base<dim, LinearAlgebra::distributed::Vector<number>>::
      clear();
  }
  template <int dim, int degree_v, typename number>
  void
  ABlockOperator<dim, degree_v, number>::evaluate_2_x_viscosity(
    const Viscosity<dim> &viscosity_function)
  {
    const unsigned int n_cells = this->data->n_cell_batches();
    FEEvaluation<dim, degree_v, degree_v + 1, dim, number> velocity(*this->data,
                                                                    0);
    viscosity_x_2.reinit(n_cells, velocity.n_q_points);
    for (unsigned int cell = 0; cell < n_cells; ++cell)
      {
        velocity.reinit(cell);
        for (unsigned int q = 0; q < velocity.n_q_points; ++q)
          {
            VectorizedArray<number> return_value =
              make_vectorized_array<number>(1.);
            for (unsigned int i = 0; i < VectorizedArray<number>::size(); ++i)
              {
                Point<dim> p;
                for (unsigned int d = 0; d < dim; ++d)
                  {
                    p(d) = velocity.quadrature_point(q)(d)[i];
                  }
                return_value[i] = 2.0 * viscosity_function.value(p);
              }
            viscosity_x_2(cell, q) = return_value;
          }
      }
  }
  template <int dim, int degree_v, typename number>
  void
  ABlockOperator<dim, degree_v, number>::local_apply(
    const dealii::MatrixFree<dim, number> &           data,
    LinearAlgebra::distributed::Vector<number> &      dst,
    const LinearAlgebra::distributed::Vector<number> &src,
    const std::pair<unsigned int, unsigned int> &     cell_range) const
  {
    FEEvaluation<dim, degree_v, degree_v + 1, dim, number> velocity(data);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        AssertDimension(viscosity_x_2.size(0), data.n_cell_batches());
        AssertDimension(viscosity_x_2.size(1), velocity.n_q_points);

        velocity.reinit(cell);
        velocity.read_dof_values(src);
        velocity.evaluate(EvaluationFlags::gradients);
        for (unsigned int q = 0; q < velocity.n_q_points; ++q)
          {
            velocity.submit_symmetric_gradient(
              viscosity_x_2(cell, q) * velocity.get_symmetric_gradient(q), q);
          }
        velocity.integrate(EvaluationFlags::gradients);
        velocity.distribute_local_to_global(dst);
      }
  }
  template <int dim, int degree_v, typename number>
  void
  ABlockOperator<dim, degree_v, number>::apply_add(
    LinearAlgebra::distributed::Vector<number> &      dst,
    const LinearAlgebra::distributed::Vector<number> &src) const
  {
    MatrixFreeOperators::Base<dim,
                              LinearAlgebra::distributed::Vector<number>>::data
      ->cell_loop(&ABlockOperator::local_apply, this, dst, src);
  }
  template <int dim, int degree_v, typename number>
  void
  ABlockOperator<dim, degree_v, number>::compute_diagonal()
  {
    this->inverse_diagonal_entries.reset(
      new DiagonalMatrix<LinearAlgebra::distributed::Vector<number>>());
    LinearAlgebra::distributed::Vector<number> &inverse_diagonal =
      this->inverse_diagonal_entries->get_vector();
    this->data->initialize_dof_vector(inverse_diagonal);
    unsigned int dummy = 0;
    this->data->cell_loop(&ABlockOperator::local_compute_diagonal,
                          this,
                          inverse_diagonal,
                          dummy);

    this->set_constrained_entries_to_one(inverse_diagonal);

    for (unsigned int i = 0; i < inverse_diagonal.locally_owned_size(); ++i)
      {
        Assert(inverse_diagonal.local_element(i) > 0.,
               ExcMessage("No diagonal entry in a positive definite operator "
                          "should be zero"));
        inverse_diagonal.local_element(i) =
          1. / inverse_diagonal.local_element(i);
      }
  }
  template <int dim, int degree_v, typename number>
  void
  ABlockOperator<dim, degree_v, number>::local_compute_diagonal(
    const MatrixFree<dim, number> &             data,
    LinearAlgebra::distributed::Vector<number> &dst,
    const unsigned int &,
    const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, degree_v, degree_v + 1, dim, number> velocity(data, 0);
    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        velocity.reinit(cell);
        AlignedVector<VectorizedArray<number>> diagonal(velocity.dofs_per_cell);
        for (unsigned int i = 0; i < velocity.dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < velocity.dofs_per_cell; ++j)
              velocity.begin_dof_values()[j] = VectorizedArray<number>();
            velocity.begin_dof_values()[i] = make_vectorized_array<number>(1.);

            velocity.evaluate(EvaluationFlags::gradients);
            for (unsigned int q = 0; q < velocity.n_q_points; ++q)
              {
                velocity.submit_symmetric_gradient(
                  viscosity_x_2(cell, q) * velocity.get_symmetric_gradient(q),
                  q);
              }
            velocity.integrate(EvaluationFlags::gradients);

            diagonal[i] = velocity.begin_dof_values()[i];
          }

        for (unsigned int i = 0; i < velocity.dofs_per_cell; ++i)
          velocity.begin_dof_values()[i] = diagonal[i];
        velocity.distribute_local_to_global(dst);
      }
  }



  template <int dim>
  class StokesProblem
  {
  public:
    StokesProblem();

    void
    run();

  private:
    void
    make_grid(const unsigned int ref = 4);
    void
    create_sinker(const unsigned int n_sinkers, const double visc_jump);
    void
    setup_system();
    void
    assemble_rhs();
    void
    solve();

    using vector_t       = LinearAlgebra::distributed::Vector<double>;
    using block_vector_t = LinearAlgebra::distributed::BlockVector<double>;

    using StokesMatrixType = StokesOperator<dim, velocity_degree, double>;
    using MassMatrixType = MassMatrixOperator<dim, velocity_degree - 1, double>;
    using LevelMatrixType = ABlockOperator<dim, velocity_degree, double>;

    unsigned int degree_u;

    FESystem<dim> fe_u;
    FE_Q<dim>     fe_p;

    parallel::distributed::Triangulation<dim> triangulation;
    DoFHandler<dim>                           dof_handler_u;
    DoFHandler<dim>                           dof_handler_p;

    std::vector<IndexSet> owned_partitioning;
    std::vector<IndexSet> relevant_partitioning;

    AffineConstraints<double> constraints_u;
    AffineConstraints<double> constraints_p;

    block_vector_t solution;
    block_vector_t system_rhs;

    StokesMatrixType stokes_matrix;
    MassMatrixType   mass_matrix;

    MGLevelObject<LevelMatrixType> mg_matrices;
    MGConstrainedDoFs              mg_constrained_dofs;

    Sinker<dim> sinker;
  };



  template <int dim>
  StokesProblem<dim>::StokesProblem()
    : degree_u(velocity_degree)
    , fe_u(FE_Q<dim>(degree_u), dim)
    , fe_p(FE_Q<dim>(degree_u - 1))
    , triangulation(MPI_COMM_WORLD,
                    typename Triangulation<dim>::MeshSmoothing(
                      Triangulation<dim>::limit_level_difference_at_vertices |
                      Triangulation<dim>::smoothing_on_refinement |
                      Triangulation<dim>::smoothing_on_coarsening),
                    parallel::distributed::Triangulation<
                      dim>::construct_multigrid_hierarchy)
    , dof_handler_u(triangulation)
    , dof_handler_p(triangulation)
  {}

  template <int dim>
  void
  StokesProblem<dim>::create_sinker(const unsigned int n_sinkers,
                                    const double       visc_jump)
  {
    sinker.problem_dim = dim;
    sinker.n_sinkers   = n_sinkers;
    std::srand(171);
    for (unsigned int s = 0; s < sinker.n_sinkers; ++s)
      {
        std::vector<double> coords(dim);
        if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
          for (unsigned int i = 0; i < dim; ++i)
            coords[i] = std::rand() / (double)RAND_MAX;

        MPI_Bcast(&(coords[0]), dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        Tensor<1, dim, double> coords_tens;
        for (unsigned int i = 0; i < dim; ++i)
          coords_tens[i] = coords[i];
        sinker.centers.push_back(Point<dim>(coords_tens));
      }

    sinker.DR_mu  = visc_jump;
    sinker.mu_max = std::sqrt(sinker.DR_mu);
    sinker.mu_min = 1.0 / std::sqrt(sinker.DR_mu);
  }


  template <int dim>
  void
  StokesProblem<dim>::make_grid(const unsigned int ref)
  {
    GridGenerator::hyper_cube(triangulation, 0, 1);
    triangulation.refine_global(ref);
  }


  template <int dim>
  void
  StokesProblem<dim>::setup_system()
  {
    dof_handler_u.clear();
    dof_handler_u.distribute_dofs(fe_u);
    dof_handler_u.distribute_mg_dofs();

    dof_handler_p.clear();
    dof_handler_p.distribute_dofs(fe_p);

    IndexSet locally_relevant_dofs_u;
    DoFTools::extract_locally_relevant_dofs(dof_handler_u,
                                            locally_relevant_dofs_u);
    constraints_u.reinit(locally_relevant_dofs_u);
    DoFTools::make_hanging_node_constraints(dof_handler_u, constraints_u);

    VectorTools::interpolate_boundary_values(
      dof_handler_u, 0, ExactSolution_BoundaryValues_u<dim>(), constraints_u);
    constraints_u.close();

    IndexSet locally_relevant_dofs_p;
    DoFTools::extract_locally_relevant_dofs(dof_handler_p,
                                            locally_relevant_dofs_p);
    constraints_p.reinit(locally_relevant_dofs_p);
    DoFTools::make_hanging_node_constraints(dof_handler_p, constraints_p);
    constraints_p.close();


    // Stokes matrix stuff...
    typename MatrixFree<dim, double>::AdditionalData additional_data_stokes;
    additional_data_stokes.tasks_parallel_scheme =
      MatrixFree<dim, double>::AdditionalData::none;
    additional_data_stokes.mapping_update_flags =
      (update_values | update_gradients | update_JxW_values |
       update_quadrature_points);

    std::vector<const DoFHandler<dim> *> stokes_dofs;
    stokes_dofs.push_back(&dof_handler_u);
    stokes_dofs.push_back(&dof_handler_p);
    std::vector<const AffineConstraints<double> *> stokes_constraints;
    stokes_constraints.push_back(&constraints_u);
    stokes_constraints.push_back(&constraints_p);

    std::shared_ptr<MatrixFree<dim, double>> stokes_mf_storage(
      new MatrixFree<dim, double>());
    stokes_mf_storage->reinit(MappingQ1<dim>{},
                              stokes_dofs,
                              stokes_constraints,
                              QGauss<1>(degree_u + 1),
                              additional_data_stokes);

    stokes_matrix.initialize(stokes_mf_storage);
    stokes_matrix.evaluate_2_x_viscosity(Viscosity<dim>(sinker));

    // Mass matrix stuff...
    typename MatrixFree<dim, double>::AdditionalData additional_data_mass;
    additional_data_mass.tasks_parallel_scheme =
      MatrixFree<dim, double>::AdditionalData::none;
    additional_data_mass.mapping_update_flags =
      (update_values | update_JxW_values | update_quadrature_points);
    std::shared_ptr<MatrixFree<dim, double>> mass_mf_storage(
      new MatrixFree<dim, double>());
    mass_mf_storage->reinit(MappingQ1<dim>{},
                            dof_handler_p,
                            constraints_p,
                            QGauss<1>(degree_u + 1),
                            additional_data_mass);

    mass_matrix.initialize(mass_mf_storage);
    mass_matrix.evaluate_1_over_viscosity(Viscosity<dim>(sinker));
    mass_matrix.compute_diagonal();

    // GMG stuff...
    const unsigned int n_levels = triangulation.n_global_levels();
    mg_matrices.clear_elements();

    mg_matrices.resize(0, n_levels - 1);

    mg_constrained_dofs.clear();
    const std::set<types::boundary_id> dirichlet_boundary = {0};
    mg_constrained_dofs.initialize(dof_handler_u);
    mg_constrained_dofs.make_zero_boundary_constraints(dof_handler_u,
                                                       dirichlet_boundary);


    for (unsigned int level = 0; level < n_levels; ++level)
      {
        IndexSet relevant_dofs;
        DoFTools::extract_locally_relevant_level_dofs(dof_handler_u,
                                                      level,
                                                      relevant_dofs);
        AffineConstraints<double> level_constraints;
        level_constraints.reinit(relevant_dofs);
        level_constraints.add_lines(
          mg_constrained_dofs.get_boundary_indices(level));
        level_constraints.close();

        typename MatrixFree<dim, double>::AdditionalData additional_data;
        additional_data.tasks_parallel_scheme =
          MatrixFree<dim, double>::AdditionalData::none;
        additional_data.mapping_update_flags =
          (update_gradients | update_JxW_values | update_quadrature_points);
        additional_data.mg_level = level;
        std::shared_ptr<MatrixFree<dim, double>> mg_mf_storage_level(
          new MatrixFree<dim, double>());
        mg_mf_storage_level->reinit(MappingQ1<dim>{},
                                    dof_handler_u,
                                    level_constraints,
                                    QGauss<1>(degree_u + 1),
                                    additional_data);

        mg_matrices[level].initialize(mg_mf_storage_level,
                                      mg_constrained_dofs,
                                      level);
        mg_matrices[level].evaluate_2_x_viscosity(Viscosity<dim>(sinker));
        mg_matrices[level].compute_diagonal();
      }

    solution.reinit(2);
    system_rhs.reinit(2);

    stokes_matrix.initialize_dof_vector(solution);
    stokes_matrix.initialize_dof_vector(system_rhs);

    solution.update_ghost_values();
    solution.collect_sizes();

    system_rhs.update_ghost_values();
    system_rhs.collect_sizes();
  }



  template <int dim>
  void
  StokesProblem<dim>::assemble_rhs()
  {
    system_rhs = 0.0;

    // Create operator with no Dirchlet info
    StokesMatrixType                                 operator_homogeneous;
    typename MatrixFree<dim, double>::AdditionalData data;
    data.tasks_parallel_scheme = MatrixFree<dim, double>::AdditionalData::none;
    data.mapping_update_flags  = (update_values | update_gradients |
                                 update_JxW_values | update_quadrature_points);

    // Create constraints with no Dirchlet info
    AffineConstraints<double> constraints_u_no_dirchlet;
    IndexSet                  locally_relevant_dofs_u;
    DoFTools::extract_locally_relevant_dofs(dof_handler_u,
                                            locally_relevant_dofs_u);
    constraints_u_no_dirchlet.reinit(locally_relevant_dofs_u);
    DoFTools::make_hanging_node_constraints(dof_handler_u,
                                            constraints_u_no_dirchlet);
    constraints_u_no_dirchlet.close();

    std::vector<const AffineConstraints<double> *> constraints_no_dirchlet;
    constraints_no_dirchlet.push_back(&constraints_u_no_dirchlet);
    constraints_no_dirchlet.push_back(&constraints_p);
    std::vector<const DoFHandler<dim> *> dofs;
    dofs.push_back(&dof_handler_u);
    dofs.push_back(&dof_handler_p);

    std::shared_ptr<MatrixFree<dim, double>> matrix_free_homogeneous(
      new MatrixFree<dim, double>());
    matrix_free_homogeneous->reinit(MappingQ1<dim>{},
                                    dofs,
                                    constraints_no_dirchlet,
                                    QGauss<1>(degree_u + 1),
                                    data);

    operator_homogeneous.initialize(matrix_free_homogeneous);
    operator_homogeneous.evaluate_2_x_viscosity(Viscosity<dim>(sinker));
    LinearAlgebra::distributed::BlockVector<double> inhomogeneity(2);
    operator_homogeneous.initialize_dof_vector(inhomogeneity);
    constraints_u.distribute(inhomogeneity.block(0));
    operator_homogeneous.vmult(system_rhs, inhomogeneity);
    system_rhs *= -1.;

    // Normal apply boundary
    RightHandSide<dim> right_hand_side(sinker);

    FEEvaluation<dim, velocity_degree, velocity_degree + 1, dim, double>
      velocity(*stokes_matrix.get_matrix_free(), 0);
    FEEvaluation<dim, velocity_degree - 1, velocity_degree + 1, 1, double>
      pressure(*stokes_matrix.get_matrix_free(), 1);

    for (unsigned int cell = 0;
         cell < stokes_matrix.get_matrix_free()->n_cell_batches();
         ++cell)
      {
        velocity.reinit(cell);
        pressure.reinit(cell);
        for (unsigned int q = 0; q < velocity.n_q_points; ++q)
          {
            Tensor<1, dim, VectorizedArray<double>> rhs_u;
            for (unsigned int d = 0; d < dim; ++d)
              rhs_u[d] = make_vectorized_array<double>(1.0);
            VectorizedArray<double> rhs_p = make_vectorized_array<double>(1.0);
            for (unsigned int i = 0; i < VectorizedArray<double>::size(); ++i)
              {
                Point<dim> p;
                for (unsigned int d = 0; d < dim; ++d)
                  p(d) = velocity.quadrature_point(q)(d)[i];

                Vector<double> rhs_temp(dim + 1);
                right_hand_side.vector_value(p, rhs_temp);

                for (unsigned int d = 0; d < dim; ++d)
                  rhs_u[d][i] = rhs_temp(d);
                rhs_p[i] = rhs_temp(dim);
              }
            velocity.submit_value(rhs_u, q);
            pressure.submit_value(rhs_p, q);
          }
        velocity.integrate(EvaluationFlags::values);
        velocity.distribute_local_to_global(system_rhs.block(0));
        pressure.integrate(EvaluationFlags::values);
        pressure.distribute_local_to_global(system_rhs.block(1));
      }
    system_rhs.compress(VectorOperation::add);
  }

  template <int dim>
  void
  StokesProblem<dim>::solve()
  {
    const double       solver_tolerance = 1e-6 * system_rhs.l2_norm();
    const unsigned int n_cheap_stokes_solver_steps     = 1000;
    const double       linear_solver_A_block_tolerance = 1e-2;
    const double       linear_solver_S_block_tolerance = 1e-6;

    // extract Stokes parts of solution vector, without any ghost elements
    block_vector_t distributed_stokes_solution(solution);

    const unsigned int block_vel = 0;
    const unsigned int block_p   = 1;

    // extract Stokes parts of rhs vector
    block_vector_t distributed_stokes_rhs(system_rhs);

    PrimitiveVectorMemory<block_vector_t> mem;

    SolverControl solver_control_cheap(n_cheap_stokes_solver_steps,
                                       solver_tolerance,
                                       false,
                                       false);

    using Transfer = MGTransferMF<dim, double>;

    Transfer mg_transfer(mg_constrained_dofs);
    mg_transfer.initialize_constraints(mg_constrained_dofs);
    mg_transfer.build(dof_handler_u);

    LevelMatrixType &  coarse_matrix = mg_matrices[0];
    SolverControl      coarse_solver_control(1000, 1e-12, false, false);
    SolverCG<vector_t> coarse_solver(coarse_solver_control);

    PreconditionIdentity coarse_prec_identity;
    MGCoarseGridIterativeSolver<vector_t,
                                SolverCG<vector_t>,
                                LevelMatrixType,
                                PreconditionIdentity>
      mg_coarse;
    mg_coarse.initialize(coarse_solver, coarse_matrix, coarse_prec_identity);

    using SmootherType = PreconditionChebyshev<LevelMatrixType, vector_t>;
    mg::SmootherRelaxation<SmootherType, vector_t>       mg_smoother;
    MGLevelObject<typename SmootherType::AdditionalData> smoother_data;
    smoother_data.resize(0, triangulation.n_global_levels() - 1);
    for (unsigned int level = 0; level < triangulation.n_global_levels();
         ++level)
      {
        if (level > 0)
          {
            smoother_data[level].smoothing_range     = 15.;
            smoother_data[level].degree              = 4;
            smoother_data[level].eig_cg_n_iterations = 10;
          }
        else
          {
            smoother_data[0].smoothing_range = 1e-3;
            smoother_data[0].degree          = numbers::invalid_unsigned_int;
            smoother_data[0].eig_cg_n_iterations = mg_matrices[0].m();
          }
        smoother_data[level].preconditioner =
          mg_matrices[level].get_matrix_diagonal_inverse();
      }
    mg_smoother.initialize(mg_matrices, smoother_data);

    mg::Matrix<vector_t> mg_matrix(mg_matrices);

    MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<LevelMatrixType>>
      mg_interface_matrices;
    mg_interface_matrices.resize(0, triangulation.n_global_levels() - 1);
    for (unsigned int level = 0; level < triangulation.n_global_levels();
         ++level)
      mg_interface_matrices[level].initialize(mg_matrices[level]);
    mg::Matrix<vector_t> mg_interface(mg_interface_matrices);

    Multigrid<vector_t> mg(
      mg_matrix, mg_coarse, mg_transfer, mg_smoother, mg_smoother, 0);
    mg.set_edge_matrices(mg_interface, mg_interface);


    PreconditionMG<dim, vector_t, Transfer> prec_A(dof_handler_u,
                                                   mg,
                                                   mg_transfer);

    using MassPrec = PreconditionChebyshev<MassMatrixType, vector_t>;
    MassPrec                          prec_S;
    typename MassPrec::AdditionalData prec_S_data;
    prec_S_data.smoothing_range     = 1e-3;
    prec_S_data.degree              = numbers::invalid_unsigned_int;
    prec_S_data.eig_cg_n_iterations = mass_matrix.m();
    prec_S_data.preconditioner      = mass_matrix.get_matrix_diagonal_inverse();
    prec_S.initialize(mass_matrix, prec_S_data);

    using A_prec_type = PreconditionMG<dim, vector_t, Transfer>;

    // create a cheap preconditioner that consists of only a single V-cycle
    const StokesSolver::BlockSchurPreconditioner<StokesMatrixType,
                                                 MassMatrixType,
                                                 A_prec_type,
                                                 MassPrec>
      preconditioner_cheap(stokes_matrix,
                           mass_matrix,
                           prec_S,
                           prec_A,
                           false,
                           linear_solver_A_block_tolerance,
                           linear_solver_S_block_tolerance);
    try
      {
        SolverFGMRES<block_vector_t> solver(
          solver_control_cheap,
          mem,
          SolverFGMRES<block_vector_t>::AdditionalData(50));

        solver.solve(stokes_matrix,
                     distributed_stokes_solution,
                     distributed_stokes_rhs,
                     preconditioner_cheap);
      }
    catch (const SolverControl::NoConvergence &)
      {
        deallog
          << "********************************************************************"
          << std::endl
          << "SOLVER DID NOT CONVERGE AFTER " << n_cheap_stokes_solver_steps
          << " ITERATIONS. res=" << solver_control_cheap.last_value()
          << std::endl
          << "********************************************************************"
          << std::endl;
      }

    constraints_u.distribute(distributed_stokes_solution.block(0));

    distributed_stokes_solution.block(block_p) *= pressure_scaling;

    solution.block(block_vel) = distributed_stokes_solution.block(block_vel);
    solution.block(block_p)   = distributed_stokes_solution.block(block_p);

    deallog << "Solved-in "
            << (solver_control_cheap.last_step() !=
                    numbers::invalid_unsigned_int ?
                  solver_control_cheap.last_step() :
                  0)
            << " iterations, final residual: "
            << solver_control_cheap.last_value() << std::endl;
  }


  template <int dim>
  void
  StokesProblem<dim>::run()
  {
    deallog << "Sinker problem in " << dim << "D." << std::endl;

    create_sinker(4, 1000);
    deallog << "n_sinker: " << sinker.n_sinkers
            << "       max/min viscosity ratio: " << sinker.DR_mu << std::endl
            << std::endl;

    unsigned int initial_ref;
    if (dim == 2)
      {
        initial_ref = 5;
      }
    else if (dim == 3)
      {
        initial_ref = 3;
      }

    unsigned int n_cycles = 1;
    if (dim == 2)
      n_cycles = 2;
    for (unsigned int cycle = 0; cycle < n_cycles; ++cycle)
      {
        if (cycle == 0)
          make_grid(initial_ref);
        else
          triangulation.refine_global();

        setup_system();
        deallog << "Number of active cells: "
                << triangulation.n_global_active_cells() << " (on "
                << triangulation.n_global_levels() << " levels)" << std::endl;
        deallog << "Number of degrees of freedom: "
                << dof_handler_u.n_dofs() + dof_handler_p.n_dofs() << " ("
                << dof_handler_u.n_dofs() << '+' << dof_handler_p.n_dofs()
                << ')' << std::endl;

        assemble_rhs();
        solve();

        deallog << std::endl;
      }
  }
} // namespace StokesClass


int
main(int argc, char *argv[])
{
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  mpi_initlog();
  try
    {
      {
        deallog.push("2d");
        StokesClass::StokesProblem<2> problem;
        problem.run();
        deallog.pop();
      }
      {
        deallog.push("3d");
        StokesClass::StokesProblem<3> problem;
        problem.run();
        deallog.pop();
      }
    }
  catch (const std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      throw;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      throw;
    }

  return 0;
}
