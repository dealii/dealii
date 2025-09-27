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

#ifndef dealii_solver_idr_h
#define dealii_solver_idr_h


#include <deal.II/base/config.h>

#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/template_constraints.h>

#include <deal.II/lac/block_vector_base.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver.h>
#include <deal.II/lac/solver_control.h>

#include <cmath>
#include <random>

DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup Solvers
 * @{
 */

namespace internal
{
  /**
   * A namespace for a helper class to the IDR(s) solver.
   */
  namespace SolverIDRImplementation
  {
    /**
     * Class to hold temporary vectors whose size depends on
     * the solver parameter s.
     */
    template <typename VectorType>
    class TmpVectors
    {
    public:
      /**
       * Constructor. Prepares an array of @p VectorType of length @p s_param.
       */
      TmpVectors(const unsigned int s_param, VectorMemory<VectorType> &vmem);

      /**
       * Destructor. Delete all allocated vectors.
       */
      ~TmpVectors() = default;

      /**
       * Get vector number @p i. If this vector was unused before, an error
       * occurs.
       */
      VectorType &
      operator[](const unsigned int i) const;

      /**
       * Get vector number @p i. Allocate it if necessary.
       *
       * If a vector must be allocated, @p temp is used to reinit it to the
       * proper dimensions.
       */
      VectorType &
      operator()(const unsigned int i, const VectorType &temp);

    private:
      /**
       * Pool where vectors are obtained from.
       */
      VectorMemory<VectorType> &mem;

      /**
       * Field for storing the vectors.
       */
      std::vector<typename VectorMemory<VectorType>::Pointer> data;
    };
  } // namespace SolverIDRImplementation
} // namespace internal

/**
 * This class implements the IDR(s) method used for solving nonsymmetric,
 * indefinite linear systems, developed in <a
 * href="https://epubs.siam.org/doi/abs/10.1137/070685804">
 * IDR(s): A Family of Simple and Fast Algorithms for Solving Large
 * Nonsymmetric Systems of Linear Equations by Martin B. van Gijzen and Peter
 * Sonneveld </a>. The implementation here is the preconditioned version from <a
 * href="https://dl.acm.org/citation.cfm?id=2049667">
 * Algorithm 913: An Elegant IDR(s) Variant that Efficiently Exploits
 * Biorthogonality Properties
 * by Martin B. van Gijzen and Peter Sonneveld</a>. The local structure
 * @p AdditionalData takes the value for the parameter s which can be any
 * integer greater than or equal to 1. For <code>s=1</code>, this method has
 * similar convergence to BiCGStab.
 *
 * @note Each iteration of IDR(s) requires <code>s+1</code> preconditioning steps
 * and matrix-vector products. In this implementation the residual is updated
 * and convergence is checked after each of these inner steps inside the outer
 * iteration. If the user enables the history data, the residual at each of
 * these steps is stored and therefore there will be multiple values per
 * iteration.
 */
template <typename VectorType = Vector<double>>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
class SolverIDR : public SolverBase<VectorType>
{
public:
  /**
   * Structure for storing additional data needed by the solver.
   */
  struct AdditionalData
  {
    /**
     * Constructor. By default, an IDR(2) method is used.
     */
    explicit AdditionalData(const unsigned int s = 2)
      : s(s)
    {}

    const unsigned int s;
  };

  /**
   * Constructor.
   */
  SolverIDR(SolverControl            &cn,
            VectorMemory<VectorType> &mem,
            const AdditionalData     &data = AdditionalData());

  /**
   * Constructor. Use an object of type GrowingVectorMemory as a default to
   * allocate memory.
   */
  explicit SolverIDR(SolverControl        &cn,
                     const AdditionalData &data = AdditionalData());

  /**
   * Virtual destructor.
   */
  virtual ~SolverIDR() override = default;

  /**
   * Solve the linear system <code>Ax=b</code> for x.
   */
  template <typename MatrixType, typename PreconditionerType>
  DEAL_II_CXX20_REQUIRES(
    (concepts::is_linear_operator_on<MatrixType, VectorType> &&
     concepts::is_linear_operator_on<PreconditionerType, VectorType>))
  void solve(const MatrixType         &A,
             VectorType               &x,
             const VectorType         &b,
             const PreconditionerType &preconditioner);

protected:
  /**
   * Interface for derived class. This function gets the current iteration
   * vector, the residual and the update vector in each step. It can be used
   * for graphical output of the convergence history.
   */
  virtual void
  print_vectors(const unsigned int step,
                const VectorType  &x,
                const VectorType  &r,
                const VectorType  &d) const;

private:
  /**
   * Additional solver parameters.
   */
  AdditionalData additional_data;
};

/** @} */
/*------------------------- Implementation ----------------------------*/

#ifndef DOXYGEN


namespace internal
{
  namespace SolverIDRImplementation
  {
    template <typename VectorType>
    inline TmpVectors<VectorType>::TmpVectors(const unsigned int        s_param,
                                              VectorMemory<VectorType> &vmem)
      : mem(vmem)
      , data(s_param)
    {}



    template <typename VectorType>
    inline VectorType &
    TmpVectors<VectorType>::operator[](const unsigned int i) const
    {
      AssertIndexRange(i, data.size());

      Assert(data[i] != nullptr, ExcNotInitialized());
      return *data[i];
    }



    template <typename VectorType>
    inline VectorType &
    TmpVectors<VectorType>::operator()(const unsigned int i,
                                       const VectorType  &temp)
    {
      AssertIndexRange(i, data.size());
      if (data[i] == nullptr)
        {
          data[i] = std::move(typename VectorMemory<VectorType>::Pointer(mem));
          data[i]->reinit(temp, true);
        }
      return *data[i];
    }



    template <typename VectorType,
              std::enable_if_t<!IsBlockVector<VectorType>::value, VectorType>
                * = nullptr>
    unsigned int
    n_blocks(const VectorType &)
    {
      return 1;
    }



    template <typename VectorType,
              std::enable_if_t<IsBlockVector<VectorType>::value, VectorType> * =
                nullptr>
    unsigned int
    n_blocks(const VectorType &vector)
    {
      return vector.n_blocks();
    }



    template <typename VectorType,
              std::enable_if_t<!IsBlockVector<VectorType>::value, VectorType>
                * = nullptr>
    VectorType &
    block(VectorType &vector, const unsigned int b)
    {
      AssertDimension(b, 0);
      return vector;
    }



    template <typename VectorType,
              std::enable_if_t<IsBlockVector<VectorType>::value, VectorType> * =
                nullptr>
    typename VectorType::BlockType &
    block(VectorType &vector, const unsigned int b)
    {
      return vector.block(b);
    }

  } // namespace SolverIDRImplementation
} // namespace internal



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
SolverIDR<VectorType>::SolverIDR(SolverControl            &cn,
                                 VectorMemory<VectorType> &mem,
                                 const AdditionalData     &data)
  : SolverBase<VectorType>(cn, mem)
  , additional_data(data)
{}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
SolverIDR<VectorType>::SolverIDR(SolverControl &cn, const AdditionalData &data)
  : SolverBase<VectorType>(cn)
  , additional_data(data)
{}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
void SolverIDR<VectorType>::print_vectors(const unsigned int,
                                          const VectorType &,
                                          const VectorType &,
                                          const VectorType &) const
{}



template <typename VectorType>
DEAL_II_CXX20_REQUIRES(concepts::is_vector_space_vector<VectorType>)
template <typename MatrixType, typename PreconditionerType>
DEAL_II_CXX20_REQUIRES(
  (concepts::is_linear_operator_on<MatrixType, VectorType> &&
   concepts::is_linear_operator_on<PreconditionerType, VectorType>))
void SolverIDR<VectorType>::solve(const MatrixType         &A,
                                  VectorType               &x,
                                  const VectorType         &b,
                                  const PreconditionerType &preconditioner)
{
  LogStream::Prefix prefix("IDR(s)");

  SolverControl::State iteration_state = SolverControl::iterate;
  unsigned int         step            = 0;

  const unsigned int s = additional_data.s;

  // Define temporary vectors which do not do not depend on s
  typename VectorMemory<VectorType>::Pointer r_pointer(this->memory);
  typename VectorMemory<VectorType>::Pointer v_pointer(this->memory);
  typename VectorMemory<VectorType>::Pointer uhat_pointer(this->memory);

  VectorType &r    = *r_pointer;
  VectorType &v    = *v_pointer;
  VectorType &uhat = *uhat_pointer;

  r.reinit(x, true);
  v.reinit(x, true);
  uhat.reinit(x, true);

  // Initial residual
  A.vmult(r, x);
  r.sadd(-1.0, 1.0, b);

  using value_type = typename VectorType::value_type;
  using real_type  = typename numbers::NumberTraits<value_type>::real_type;

  // Check for convergent initial guess
  real_type res   = r.l2_norm();
  iteration_state = this->iteration_status(step, res, x);
  if (iteration_state == SolverControl::success)
    return;

  // Initialize sets of vectors/matrices whose size dependent on s
  internal::SolverIDRImplementation::TmpVectors<VectorType> G(s, this->memory);
  internal::SolverIDRImplementation::TmpVectors<VectorType> U(s, this->memory);
  internal::SolverIDRImplementation::TmpVectors<VectorType> Q(s, this->memory);
  FullMatrix<value_type>                                    M(s, s);

  // Random number generator for vector entries of
  // Q (normal distribution, mean=0 sigma=1)
  std::mt19937               rng;
  std::normal_distribution<> normal_distribution(0.0, 1.0);
  for (unsigned int i = 0; i < s; ++i)
    {
      // Initialize vectors
      G(i, x);
      U(i, x);

      // Compute random set of s orthonormalized vectors Q
      // Note: the first vector is chosen to be the initial
      // residual to match BiCGStab (as is done in comparisons
      // with BiCGStab in the papers listed in the documentation
      // of this function)
      VectorType &tmp_q = Q(i, x);
      if (i != 0)
        {
          for (unsigned int b = 0;
               b < internal::SolverIDRImplementation::n_blocks(tmp_q);
               ++b)
            for (auto indx : internal::SolverIDRImplementation::block(tmp_q, b)
                               .locally_owned_elements())
              internal::SolverIDRImplementation::block(tmp_q, b)(indx) =
                normal_distribution(rng);
          tmp_q.compress(VectorOperation::insert);
        }
      else
        tmp_q = r;

      for (unsigned int j = 0; j < i; ++j)
        {
          v = Q[j];
          v *= (v * tmp_q) / (tmp_q * tmp_q);
          tmp_q.add(-1.0, v);
        }

      if (i != 0)
        tmp_q *= 1.0 / tmp_q.l2_norm();

      M(i, i) = 1.;
    }

  value_type omega = 1.;

  bool early_exit = false;

  // Outer iteration
  while (iteration_state == SolverControl::iterate)
    {
      ++step;

      // Compute phi
      Vector<value_type> phi(s);
      for (unsigned int i = 0; i < s; ++i)
        phi(i) = Q[i] * r;

      // Inner iteration over s
      for (unsigned int k = 0; k < s; ++k)
        {
          // Solve M(k:s)*gamma = phi(k:s)
          Vector<value_type> gamma(s - k);
          {
            Vector<value_type>        phik(s - k);
            FullMatrix<value_type>    Mk(s - k, s - k);
            std::vector<unsigned int> indices;
            unsigned int              j = 0;
            for (unsigned int i = k; i < s; ++i, ++j)
              {
                indices.push_back(i);
                phik(j) = phi(i);
              }
            Mk.extract_submatrix_from(M, indices, indices);

            FullMatrix<value_type> Mk_inv(s - k, s - k);
            Mk_inv.invert(Mk);
            Mk_inv.vmult(gamma, phik);
          }

          v = r;

          if (step > 1)
            {
              for (unsigned int i = k, j = 0; i < s; ++i, ++j)
                v.add(-gamma(j), G[i]);
            }

          preconditioner.vmult(uhat, v);

          if (step > 1)
            {
              uhat.sadd(omega, gamma(0), U[k]);
              for (unsigned int i = k + 1, j = 1; i < s; ++i, ++j)
                uhat.add(gamma(j), U[i]);
            }
          else
            uhat *= omega;

          A.vmult(G[k], uhat);

          // Update G and U
          // Orthogonalize G[k] to Q0,..,Q_{k-1} and update uhat
          if (k > 0)
            {
              value_type alpha = Q[0] * G[k] / M(0, 0);
              for (unsigned int i = 1; i < k; ++i)
                {
                  const value_type alpha_old = alpha;
                  alpha = G[k].add_and_dot(-alpha, G[i - 1], Q[i]) / M(i, i);

                  // update uhat every other iteration to reduce vector access
                  if (i % 2 == 1)
                    uhat.add(-alpha_old, U[i - 1], -alpha, U[i]);
                }
              M(k, k) = G[k].add_and_dot(-alpha, G[k - 1], Q[k]);
              if (k % 2 == 1)
                uhat.add(-alpha, U[k - 1]);
            }
          else
            M(k, k) = G[k] * Q[k];

          U[k].swap(uhat);

          // Update kth column of M
          for (unsigned int i = k + 1; i < s; ++i)
            M(i, k) = Q[i] * G[k];

          // Orthogonalize r to Q0,...,Qk, update x
          {
            const value_type beta = phi(k) / M(k, k);
            r.add(-beta, G[k]);
            x.add(beta, U[k]);

            print_vectors(step, x, r, U[k]);

            // Check for early convergence. If so, store
            // information in early_exit so that outer iteration
            // is broken before recomputing the residual
            res             = r.l2_norm();
            iteration_state = this->iteration_status(step, res, x);
            if (iteration_state != SolverControl::iterate)
              {
                early_exit = true;
                break;
              }

            // Update phi
            if (k + 1 < s)
              {
                for (unsigned int i = 0; i < k + 1; ++i)
                  phi(i) = 0.0;
                for (unsigned int i = k + 1; i < s; ++i)
                  phi(i) -= beta * M(i, k);
              }
          }
        }
      if (early_exit == true)
        break;

      // Update r and x
      preconditioner.vmult(uhat, r);
      A.vmult(v, uhat);

      omega = (v * r) / (v * v);

      res = std::sqrt(r.add_and_dot(-1.0 * omega, v, r));
      x.add(omega, uhat);

      print_vectors(step, x, r, uhat);

      // Check for convergence
      iteration_state = this->iteration_status(step, res, x);
      if (iteration_state != SolverControl::iterate)
        break;
    }

  if (iteration_state != SolverControl::success)
    AssertThrow(false, SolverControl::NoConvergence(step, res));
}


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
