// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2018 by the deal.II authors
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

#ifndef dealii_solver_cg_h
#define dealii_solver_cg_h

#include <cmath>
#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/lac/solver.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/tridiagonal_matrix.h>

DEAL_II_NAMESPACE_OPEN

// forward declaration
class PreconditionIdentity;

/*!@addtogroup Solvers */
/*@{*/

/**
 * This class implements the preconditioned Conjugate Gradients (CG)
 * method that can be used to solve linear systems with a symmetric positive
 * definite matrix. This
 * class is used first in step-3 and step-4, but is used in many other
 * tutorial programs as well. Like all other solver classes, it can work on
 * any kind of vector and matrix as long as they satisfy certain requirements
 * (for the requirements on matrices and vectors in order to work with this
 * class, see the documentation of the Solver base class). The type of the
 * solution vector must be passed as template argument, and defaults to
 * dealii::Vector<double>.
 *
 * @note This version of CG is taken from D. Braess's book "Finite Elements".
 * It requires a symmetric preconditioner (i.e., for example, SOR is not a
 * possible choice).
 *
 *
 * <h3>Eigenvalue computation</h3>
 *
 * The cg-method performs an orthogonal projection of the original
 * preconditioned linear system to another system of smaller dimension.
 * Furthermore, the projected matrix @p T is tri-diagonal. Since the
 * projection is orthogonal, the eigenvalues of @p T approximate those of the
 * original preconditioned matrix @p PA. In fact, after @p n steps, where @p n
 * is the dimension of the original system, the eigenvalues of both matrices
 * are equal. But, even for small numbers of iteration steps, the condition
 * number of @p T is a good estimate for the one of @p PA.
 *
 * After @p m steps the matrix T_m can be written in terms of the coefficients
 * @p alpha and @p beta as the tri-diagonal matrix with diagonal elements
 * <tt>1/alpha_0</tt>, <tt>1/alpha_1 + beta_0/alpha_0</tt>, ...,
 * <tt>1/alpha_{m-1</tt>+beta_{m-2}/alpha_{m-2}} and off-diagonal elements
 * <tt>sqrt(beta_0)/alpha_0</tt>, ..., <tt>sqrt(beta_{m-2</tt>)/alpha_{m-2}}.
 * The eigenvalues of this matrix can be computed by postprocessing.
 *
 * @see Y. Saad: "Iterative methods for Sparse Linear Systems", section 6.7.3
 * for details.
 *
 * The coefficients, eigenvalues and condition number (computed as the ratio
 * of the largest over smallest eigenvalue) can be obtained by connecting a
 * function as a slot to the solver using one of the functions @p
 * connect_coefficients_slot, @p connect_eigenvalues_slot and @p
 * connect_condition_number_slot. These slots will then be called from the
 * solver with the estimates as argument.
 *
 * <h3>Observing the progress of linear solver iterations</h3>
 *
 * The solve() function of this class uses the mechanism described in the
 * Solver base class to determine convergence. This mechanism can also be used
 * to observe the progress of the iteration.
 *
 *
 * @author W. Bangerth, G. Kanschat, R. Becker and F.-T. Suttmeier
 */
template <typename VectorType = Vector<double>>
class SolverCG : public Solver<VectorType>
{
public:
  /**
   * Declare type for container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * Standardized data struct to pipe additional data to the solver.
   * Here, it doesn't store anything but just exists for consistency
   * with the other solver classes.
   */
  struct AdditionalData
  {};

  /**
   * Constructor.
   */
  SolverCG(SolverControl&            cn,
           VectorMemory<VectorType>& mem,
           const AdditionalData&     data = AdditionalData());

  /**
   * Constructor. Use an object of type GrowingVectorMemory as a default to
   * allocate memory.
   */
  SolverCG(SolverControl& cn, const AdditionalData& data = AdditionalData());

  /**
   * Virtual destructor.
   */
  virtual ~SolverCG() override = default;

  /**
   * Solve the linear system $Ax=b$ for x.
   */
  template <typename MatrixType, typename PreconditionerType>
  void
  solve(const MatrixType&         A,
        VectorType&               x,
        const VectorType&         b,
        const PreconditionerType& preconditioner);

  /**
   * Connect a slot to retrieve the CG coefficients. The slot will be called
   * with alpha as the first argument and with beta as the second argument,
   * where alpha and beta follow the notation in Y. Saad: "Iterative methods
   * for Sparse Linear Systems", section 6.7. Called once per iteration
   */
  boost::signals2::connection
  connect_coefficients_slot(const std::function<void(double, double)>& slot);

  /**
   * Connect a slot to retrieve the estimated condition number. Called on each
   * iteration if every_iteration=true, otherwise called once when iterations
   * are ended (i.e., either because convergence has been achieved, or because
   * divergence has been detected).
   */
  boost::signals2::connection
  connect_condition_number_slot(const std::function<void(double)>& slot,
                                const bool every_iteration = false);

  /**
   * Connect a slot to retrieve the estimated eigenvalues. Called on each
   * iteration if every_iteration=true, otherwise called once when iterations
   * are ended (i.e., either because convergence has been achieved, or because
   * divergence has been detected).
   */
  boost::signals2::connection
  connect_eigenvalues_slot(
    const std::function<void(const std::vector<double>&)>& slot,
    const bool every_iteration = false);

protected:
  /**
   * Interface for derived class. This function gets the current iteration
   * vector, the residual and the update vector in each step. It can be used
   * for graphical output of the convergence history.
   */
  virtual void
  print_vectors(const unsigned int step,
                const VectorType&  x,
                const VectorType&  r,
                const VectorType&  d) const;

  /**
   * Estimates the eigenvalues from diagonal and offdiagonal. Uses these
   * estimate to compute the condition number. Calls the signals
   * eigenvalues_signal and cond_signal with these estimates as arguments.
   */
  static void
  compute_eigs_and_cond(
    const std::vector<double>& diagonal,
    const std::vector<double>& offdiagonal,
    const boost::signals2::signal<void(const std::vector<double>&)>&
                                                 eigenvalues_signal,
    const boost::signals2::signal<void(double)>& cond_signal);

  /**
   * Additional parameters.
   */
  AdditionalData additional_data;

  /**
   * Signal used to retrieve the CG coefficients. Called on each iteration.
   */
  boost::signals2::signal<void(double, double)> coefficients_signal;

  /**
   * Signal used to retrieve the estimated condition number. Called once when
   * all iterations are ended.
   */
  boost::signals2::signal<void(double)> condition_number_signal;

  /**
   * Signal used to retrieve the estimated condition numbers. Called on each
   * iteration.
   */
  boost::signals2::signal<void(double)> all_condition_numbers_signal;

  /**
   * Signal used to retrieve the estimated eigenvalues. Called once when all
   * iterations are ended.
   */
  boost::signals2::signal<void(const std::vector<double>&)> eigenvalues_signal;

  /**
   * Signal used to retrieve the estimated eigenvalues. Called on each
   * iteration.
   */
  boost::signals2::signal<void(const std::vector<double>&)>
    all_eigenvalues_signal;
};

/*@}*/

/*------------------------- Implementation ----------------------------*/

#ifndef DOXYGEN

template <typename VectorType>
SolverCG<VectorType>::SolverCG(SolverControl&            cn,
                               VectorMemory<VectorType>& mem,
                               const AdditionalData&     data)
  : Solver<VectorType>(cn, mem), additional_data(data)
{}

template <typename VectorType>
SolverCG<VectorType>::SolverCG(SolverControl& cn, const AdditionalData& data)
  : Solver<VectorType>(cn), additional_data(data)
{}

template <typename VectorType>
void
SolverCG<VectorType>::print_vectors(const unsigned int,
                                    const VectorType&,
                                    const VectorType&,
                                    const VectorType&) const
{}

template <typename VectorType>
inline void
SolverCG<VectorType>::compute_eigs_and_cond(
  const std::vector<double>& diagonal,
  const std::vector<double>& offdiagonal,
  const boost::signals2::signal<void(const std::vector<double>&)>&
                                               eigenvalues_signal,
  const boost::signals2::signal<void(double)>& cond_signal)
{
  //Avoid computing eigenvalues unless they are needed.
  if(!cond_signal.empty() || !eigenvalues_signal.empty())
    {
      TridiagonalMatrix<double> T(diagonal.size(), true);
      for(size_type i = 0; i < diagonal.size(); ++i)
        {
          T(i, i) = diagonal[i];
          if(i < diagonal.size() - 1)
            T(i, i + 1) = offdiagonal[i];
        }
      T.compute_eigenvalues();
      //Need two eigenvalues to estimate the condition number.
      if(diagonal.size() > 1)
        {
          double condition_number = T.eigenvalue(T.n() - 1) / T.eigenvalue(0);
          cond_signal(condition_number);
        }
      //Avoid copying the eigenvalues of T to a vector unless a signal is
      //connected.
      if(!eigenvalues_signal.empty())
        {
          std::vector<double> eigenvalues(T.n());
          for(unsigned int j = 0; j < T.n(); ++j)
            {
              eigenvalues.at(j) = T.eigenvalue(j);
            }
          eigenvalues_signal(eigenvalues);
        }
    }
}

template <typename VectorType>
template <typename MatrixType, typename PreconditionerType>
void
SolverCG<VectorType>::solve(const MatrixType&         A,
                            VectorType&               x,
                            const VectorType&         b,
                            const PreconditionerType& preconditioner)
{
  SolverControl::State conv = SolverControl::iterate;

  LogStream::Prefix prefix("cg");

  // Memory allocation
  typename VectorMemory<VectorType>::Pointer g_pointer(this->memory);
  typename VectorMemory<VectorType>::Pointer d_pointer(this->memory);
  typename VectorMemory<VectorType>::Pointer h_pointer(this->memory);

  // define some aliases for simpler access
  VectorType& g = *g_pointer;
  VectorType& d = *d_pointer;
  VectorType& h = *h_pointer;

  // Should we build the matrix for eigenvalue computations?
  const bool do_eigenvalues
    = !condition_number_signal.empty() || !all_condition_numbers_signal.empty()
      || !eigenvalues_signal.empty() || !all_eigenvalues_signal.empty();

  // vectors used for eigenvalue
  // computations
  std::vector<double> diagonal;
  std::vector<double> offdiagonal;

  int    it  = 0;
  double res = -std::numeric_limits<double>::max();

  double eigen_beta_alpha = 0;

  // resize the vectors, but do not set
  // the values since they'd be overwritten
  // soon anyway.
  g.reinit(x, true);
  d.reinit(x, true);
  h.reinit(x, true);

  double gh, beta;

  // compute residual. if vector is
  // zero, then short-circuit the
  // full computation
  if(!x.all_zero())
    {
      A.vmult(g, x);
      g.add(-1., b);
    }
  else
    g.equ(-1., b);
  res = g.l2_norm();

  conv = this->iteration_status(0, res, x);
  if(conv != SolverControl::iterate)
    return;

  if(std::is_same<PreconditionerType, PreconditionIdentity>::value == false)
    {
      preconditioner.vmult(h, g);

      d.equ(-1., h);

      gh = g * h;
    }
  else
    {
      d.equ(-1., g);
      gh = res * res;
    }

  while(conv == SolverControl::iterate)
    {
      it++;
      A.vmult(h, d);

      double alpha = d * h;
      Assert(alpha != 0., ExcDivideByZero());
      alpha = gh / alpha;

      x.add(alpha, d);
      res = std::sqrt(g.add_and_dot(alpha, h, g));

      print_vectors(it, x, g, d);

      conv = this->iteration_status(it, res, x);
      if(conv != SolverControl::iterate)
        break;

      if(std::is_same<PreconditionerType, PreconditionIdentity>::value == false)
        {
          preconditioner.vmult(h, g);

          beta = gh;
          Assert(beta != 0., ExcDivideByZero());
          gh   = g * h;
          beta = gh / beta;
          d.sadd(beta, -1., h);
        }
      else
        {
          beta = gh;
          gh   = res * res;
          beta = gh / beta;
          d.sadd(beta, -1., g);
        }

      this->coefficients_signal(alpha, beta);
      // set up the vectors
      // containing the diagonal
      // and the off diagonal of
      // the projected matrix.
      if(do_eigenvalues)
        {
          diagonal.push_back(1. / alpha + eigen_beta_alpha);
          eigen_beta_alpha = beta / alpha;
          offdiagonal.push_back(std::sqrt(beta) / alpha);
        }
      compute_eigs_and_cond(diagonal,
                            offdiagonal,
                            all_eigenvalues_signal,
                            all_condition_numbers_signal);
    }

  compute_eigs_and_cond(
    diagonal, offdiagonal, eigenvalues_signal, condition_number_signal);

  // in case of failure: throw exception
  if(conv != SolverControl::success)
    AssertThrow(false, SolverControl::NoConvergence(it, res));
  // otherwise exit as normal
}

template <typename VectorType>
boost::signals2::connection
SolverCG<VectorType>::connect_coefficients_slot(
  const std::function<void(double, double)>& slot)
{
  return coefficients_signal.connect(slot);
}

template <typename VectorType>
boost::signals2::connection
SolverCG<VectorType>::connect_condition_number_slot(
  const std::function<void(double)>& slot,
  const bool                         every_iteration)
{
  if(every_iteration)
    {
      return all_condition_numbers_signal.connect(slot);
    }
  else
    {
      return condition_number_signal.connect(slot);
    }
}

template <typename VectorType>
boost::signals2::connection
SolverCG<VectorType>::connect_eigenvalues_slot(
  const std::function<void(const std::vector<double>&)>& slot,
  const bool                                             every_iteration)
{
  if(every_iteration)
    {
      return all_eigenvalues_signal.connect(slot);
    }
  else
    {
      return eigenvalues_signal.connect(slot);
    }
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
