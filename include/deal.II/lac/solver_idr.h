// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2020 by the deal.II authors
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

#ifndef dealii_solver_idr_h
#define dealii_solver_idr_h


#include <deal.II/base/config.h>

#include <deal.II/base/logstream.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver.h>
#include <deal.II/lac/solver_control.h>

#include <cmath>
#include <random>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup Solvers */
/*@{*/

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
      VectorType &operator[](const unsigned int i) const;

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
 *
 * @author Conrad Clevenger, 2019
 */
template <class VectorType = Vector<double>>
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
  SolverIDR(SolverControl &           cn,
            VectorMemory<VectorType> &mem,
            const AdditionalData &    data = AdditionalData());

  /**
   * Constructor. Use an object of type GrowingVectorMemory as a default to
   * allocate memory.
   */
  explicit SolverIDR(SolverControl &       cn,
                     const AdditionalData &data = AdditionalData());

  /**
   * Virtual destructor.
   */
  virtual ~SolverIDR() override = default;

  /**
   * Solve the linear system <code>Ax=b</code> for x.
   */
  template <typename MatrixType, typename PreconditionerType>
  void
  solve(const MatrixType &        A,
        VectorType &              x,
        const VectorType &        b,
        const PreconditionerType &preconditioner);

protected:
  /**
   * Interface for derived class. This function gets the current iteration
   * vector, the residual and the update vector in each step. It can be used
   * for graphical output of the convergence history.
   */
  virtual void
  print_vectors(const unsigned int step,
                const VectorType & x,
                const VectorType & r,
                const VectorType & d) const;

private:
  /**
   * Additional solver parameters.
   */
  AdditionalData additional_data;
};

/*@}*/
/*------------------------- Implementation ----------------------------*/

#ifndef DOXYGEN


namespace internal
{
  namespace SolverIDRImplementation
  {
    template <class VectorType>
    inline TmpVectors<VectorType>::TmpVectors(const unsigned int        s_param,
                                              VectorMemory<VectorType> &vmem)
      : mem(vmem)
      , data(s_param)
    {}



    template <class VectorType>
    inline VectorType &TmpVectors<VectorType>::
                       operator[](const unsigned int i) const
    {
      AssertIndexRange(i, data.size());

      Assert(data[i] != nullptr, ExcNotInitialized());
      return *data[i];
    }



    template <class VectorType>
    inline VectorType &
    TmpVectors<VectorType>::operator()(const unsigned int i,
                                       const VectorType & temp)
    {
      AssertIndexRange(i, data.size());
      if (data[i] == nullptr)
        {
          data[i] = std::move(typename VectorMemory<VectorType>::Pointer(mem));
          data[i]->reinit(temp);
        }
      return *data[i];
    }
  } // namespace SolverIDRImplementation
} // namespace internal



template <class VectorType>
SolverIDR<VectorType>::SolverIDR(SolverControl &           cn,
                                 VectorMemory<VectorType> &mem,
                                 const AdditionalData &    data)
  : SolverBase<VectorType>(cn, mem)
  , additional_data(data)
{}



template <class VectorType>
SolverIDR<VectorType>::SolverIDR(SolverControl &cn, const AdditionalData &data)
  : SolverBase<VectorType>(cn)
  , additional_data(data)
{}



template <class VectorType>
void
SolverIDR<VectorType>::print_vectors(const unsigned int,
                                     const VectorType &,
                                     const VectorType &,
                                     const VectorType &) const
{}



template <class VectorType>
template <typename MatrixType, typename PreconditionerType>
void
SolverIDR<VectorType>::solve(const MatrixType &        A,
                             VectorType &              x,
                             const VectorType &        b,
                             const PreconditionerType &preconditioner)
{
  LogStream::Prefix prefix("IDR(s)");

  SolverControl::State iteration_state = SolverControl::iterate;
  unsigned int         step            = 0;

  const unsigned int s = additional_data.s;

  // Define temporary vectors which do not do not depend on s
  typename VectorMemory<VectorType>::Pointer r_pointer(this->memory);
  typename VectorMemory<VectorType>::Pointer v_pointer(this->memory);
  typename VectorMemory<VectorType>::Pointer vhat_pointer(this->memory);
  typename VectorMemory<VectorType>::Pointer uhat_pointer(this->memory);
  typename VectorMemory<VectorType>::Pointer ghat_pointer(this->memory);

  VectorType &r    = *r_pointer;
  VectorType &v    = *v_pointer;
  VectorType &vhat = *vhat_pointer;
  VectorType &uhat = *uhat_pointer;
  VectorType &ghat = *ghat_pointer;

  r.reinit(x, true);
  v.reinit(x, true);
  vhat.reinit(x, true);
  uhat.reinit(x, true);
  ghat.reinit(x, true);

  // Initial residual
  A.vmult(r, x);
  r.sadd(-1.0, 1.0, b);

  // Check for convergent initial guess
  double res      = r.l2_norm();
  iteration_state = this->iteration_status(step, res, x);
  if (iteration_state == SolverControl::success)
    return;

  // Initialize sets of vectors/matrices whose size dependent on s
  internal::SolverIDRImplementation::TmpVectors<VectorType> G(s, this->memory);
  internal::SolverIDRImplementation::TmpVectors<VectorType> U(s, this->memory);
  internal::SolverIDRImplementation::TmpVectors<VectorType> Q(s, this->memory);
  FullMatrix<double>                                        M(s, s);

  // Random number generator for vector entries of
  // Q (normal distribution, mean=0 sigma=1)
  std::mt19937               rng;
  std::normal_distribution<> normal_distribution(0.0, 1.0);
  for (unsigned int i = 0; i < s; ++i)
    {
      VectorType &tmp_g = G(i, x);
      VectorType &tmp_u = U(i, x);
      tmp_g             = 0;
      tmp_u             = 0;

      // Compute random set of s orthonormalized vectors Q
      // Note: the first vector is chosen to be the initial
      // residual to match BiCGStab (as is done in comparisons
      // with BiCGStab in the papers listed in the documentation
      // of this function)
      VectorType &tmp_q = Q(i, x);
      if (i != 0)
        {
          for (auto indx : tmp_q.locally_owned_elements())
            tmp_q(indx) = normal_distribution(rng);
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

  double omega = 1.;

  bool early_exit = false;

  // Outer iteration
  while (iteration_state == SolverControl::iterate)
    {
      ++step;

      // Compute phi
      Vector<double> phi(s);
      for (unsigned int i = 0; i < s; ++i)
        phi(i) = Q[i] * r;

      // Inner iteration over s
      for (unsigned int k = 0; k < s; ++k)
        {
          // Solve M(k:s)*gamma = phi(k:s)
          Vector<double> gamma(s - k);
          {
            Vector<double>            phik(s - k);
            FullMatrix<double>        Mk(s - k, s - k);
            std::vector<unsigned int> indices;
            unsigned int              j = 0;
            for (unsigned int i = k; i < s; ++i, ++j)
              {
                indices.push_back(i);
                phik(j) = phi(i);
              }
            Mk.extract_submatrix_from(M, indices, indices);

            FullMatrix<double> Mk_inv(s - k, s - k);
            Mk_inv.invert(Mk);
            Mk_inv.vmult(gamma, phik);
          }

          {
            v = r;

            unsigned int j = 0;
            for (unsigned int i = k; i < s; ++i, ++j)
              v.add(-1.0 * gamma(j), G[i]);
            preconditioner.vmult(vhat, v);

            uhat = vhat;
            uhat *= omega;
            j = 0;
            for (unsigned int i = k; i < s; ++i, ++j)
              uhat.add(gamma(j), U[i]);
            A.vmult(ghat, uhat);
          }

          // Update G and U
          // Orthogonalize ghat to Q0,..,Q_{k-1}
          // and update uhat
          for (unsigned int i = 0; i < k; ++i)
            {
              double alpha = (Q[i] * ghat) / M(i, i);
              ghat.add(-alpha, G[i]);
              uhat.add(-alpha, U[i]);
            }
          G[k] = ghat;
          U[k] = uhat;

          // Update kth column of M
          for (unsigned int i = k; i < s; ++i)
            M(i, k) = Q[i] * G[k];

          // Orthogonalize r to Q0,...,Qk,
          // update x
          {
            double beta = phi(k) / M(k, k);
            r.add(-1.0 * beta, G[k]);
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
      preconditioner.vmult(vhat, r);
      A.vmult(v, vhat);

      omega = (v * r) / (v * v);

      r.add(-1.0 * omega, v);
      x.add(omega, vhat);

      print_vectors(step, x, r, vhat);

      // Check for convergence
      res             = r.l2_norm();
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
