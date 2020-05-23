// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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

#ifndef dealii_lac_utilities_h
#define dealii_lac_utilities_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/signaling_nan.h>

#include <deal.II/lac/lapack_support.h>
#include <deal.II/lac/vector_memory.h>

#include <array>
#include <complex>

DEAL_II_NAMESPACE_OPEN

namespace Utilities
{
  /**
   * A collection of linear-algebra utilities.
   */
  namespace LinearAlgebra
  {
    /**
     * Return the elements of a continuous Givens rotation matrix and
     * the norm of the input vector.
     *
     * That is for a given
     * pair @p x and @p y, return $c$ , $s$ and $\sqrt{x^2+y^2}$ such that
     * \f[
     * \begin{bmatrix}
     * c  & s \\
     * -s & c
     * \end{bmatrix}
     * \begin{bmatrix}
     * x \\
     * y
     * \end{bmatrix}
     * =
     * \begin{bmatrix}
     * \sqrt{x^2+y^2} \\
     * 0
     * \end{bmatrix}
     * \f]
     *
     * @note The function is implemented for real valued numbers only.
     *
     * @author Denis Davydov, 2017
     */
    template <typename NumberType>
    std::array<NumberType, 3>
    givens_rotation(const NumberType &x, const NumberType &y);

    /**
     * Return the elements of a hyperbolic rotation matrix.
     *
     * That is for a given
     * pair @p x and @p y, return $c$ , $s$ and $r$ such that
     * \f[
     * \begin{bmatrix}
     * c  & -s \\
     * -s & c
     * \end{bmatrix}
     * \begin{bmatrix}
     * x \\
     * y
     * \end{bmatrix}
     * =
     * \begin{bmatrix}
     * r \\
     * 0
     * \end{bmatrix}
     * \f]
     *
     * Real valued solution only exists if $|x|>|g|$, the function will
     * throw an error otherwise.
     *
     * @note The function is implemented for real valued numbers only.
     *
     * @author Denis Davydov, 2017
     */
    template <typename NumberType>
    std::array<NumberType, 3>
    hyperbolic_rotation(const NumberType &x, const NumberType &y);

    /**
     * Estimate an upper bound for the largest eigenvalue of @p H by a @p k -step
     * Lanczos process starting from the initial vector @p v0. Typical
     * values of @p k are below 10. This estimator computes a k-step Lanczos
     * decomposition $H V_k=V_k T_k+f_k e_k^T$ where $V_k$ contains k Lanczos
     * basis, $V_k^TV_k=I_k$, $T_k$ is the tridiagonal Lanczos matrix, $f_k$ is
     * a residual vector $f_k^TV_k=0$, and $e_k$ is the k-th canonical basis of
     * $R^k$. The returned value is $ ||T_k||_2 + ||f_k||_2$.
     * If @p eigenvalues is not <code>nullptr</code>, the eigenvalues of $T_k$ will be written there.
     *
     * @p vector_memory is used to allocate memory for temporary vectors.
     * OperatorType has to provide <code>vmult</code> operation with
     * VectorType.
     *
     * This function implements the algorithm from
     * @code{.bib}
     * @article{Zhou2006,
     *   Title   = {Self-consistent-field Calculations Using Chebyshev-filtered
     *              Subspace Iteration},
     *   Author  = {Zhou, Yunkai and Saad, Yousef and Tiago, Murilo L. and
     *              Chelikowsky, James R.},
     *   Journal = {Journal of Computational Physics},
     *   Year    = {2006},
     *   Volume  = {219},
     *   Pages   = {172--184},
     * }
     * @endcode
     *
     * @note This function uses Lapack routines to compute the largest
     * eigenvalue of $T_k$.
     *
     * @note This function provides an alternate estimate to that obtained from
     * several steps of SolverCG with
     * SolverCG<VectorType>::connect_eigenvalues_slot().
     *
     * @author Denis Davydov, 2017
     */
    template <typename OperatorType, typename VectorType>
    double
    lanczos_largest_eigenvalue(const OperatorType &      H,
                               const VectorType &        v0,
                               const unsigned int        k,
                               VectorMemory<VectorType> &vector_memory,
                               std::vector<double> *     eigenvalues = nullptr);

    /**
     * Apply Chebyshev polynomial of the operator @p H to @p x. For a
     * non-defective operator $H$ with a complete set of eigenpairs
     * $H \psi_i = \lambda_i \psi_i$, the action of a polynomial filter $p$ is
     * given by $p(H)x =\sum_i a_i p(\lambda_i) \psi_i$, where $x=: \sum_i a_i
     * \psi_i$. Thus by appropriately choosing the polynomial filter, one can
     * alter the eigenmodes contained in $x$.
     *
     * This function uses Chebyshev polynomials of first kind. Below is an
     * example of polynomial $T_n(x)$ of degree $n=8$ normalized to unity at
     * $-1.2$. <table> <tr> <td align="center">
     *       @image html chebyshev8.png
     *     </td>
     *   </tr>
     * </table>
     * By introducing a linear mapping $L$ from @p unwanted_spectrum to
     * $[-1,1]$, we can dump the corresponding modes in @p x. The higher
     * the polynomial degree $n$, the more rapid it grows outside of the
     * $[-1,1]$. In order to avoid numerical overflow, we normalize
     * polynomial filter to unity at @p tau. Thus, the filtered operator
     * is $p(H) = T_n(L(H))/T_n(L(\tau))$.
     *
     * The action of the Chebyshev filter only requires
     * evaluation of <code>vmult()</code> of @p H and is based on the
     * recursion equation for Chebyshev polynomial of degree $n$:
     * $T_{n}(x) = 2x T_{n-1}(x) - T_{n-2}(x)$ with $T_0(x)=1$ and $T_1(x)=x$.
     *
     * @p vector_memory is used to allocate memory for temporary objects.
     *
     * This function implements the algorithm (with a minor fix of sign of
     * $\sigma_1$) from
     * @code{.bib}
     * @article{Zhou2014,
     *   Title   = {Chebyshev-filtered subspace iteration method free of sparse
     *              diagonalization for solving the Kohn--Sham equation},
     *   Author  = {Zhou, Yunkai and Chelikowsky, James R and Saad, Yousef},
     *   Journal = {Journal of Computational Physics},
     *   Year    = {2014},
     *   Volume  = {274},
     *   Pages   = {770--782},
     * }
     * @endcode
     *
     * @note If @p tau is equal to
     * <code>std::numeric_limits<double>::infinity()</code>, no normalization
     * will be performed.
     *
     * @author Denis Davydov, 2017
     */
    template <typename OperatorType, typename VectorType>
    void
    chebyshev_filter(VectorType &                    x,
                     const OperatorType &            H,
                     const unsigned int              n,
                     const std::pair<double, double> unwanted_spectrum,
                     const double                    tau,
                     VectorMemory<VectorType> &      vector_memory);

  } // namespace LinearAlgebra

} // namespace Utilities


/*------------------------- Implementation ----------------------------*/

#ifndef DOXYGEN

namespace internal
{
  namespace UtilitiesImplementation
  {
    // We want to avoid including our own LAPACK wrapper header in any external
    // headers to avoid possible conflicts with other packages that may define
    // their own such header. At the same time we want to be able to call some
    // LAPACK functions from the template functions below. To resolve both
    // problems define some extra wrappers here that can be in the header:
    template <typename Number>
    void
    call_stev(const char            jobz,
              const types::blas_int n,
              Number *              d,
              Number *              e,
              Number *              z,
              const types::blas_int ldz,
              Number *              work,
              types::blas_int *     info);
  } // namespace UtilitiesImplementation
} // namespace internal

namespace Utilities
{
  namespace LinearAlgebra
  {
    template <typename NumberType>
    std::array<std::complex<NumberType>, 3>
    hyperbolic_rotation(const std::complex<NumberType> & /*f*/,
                        const std::complex<NumberType> & /*g*/)
    {
      AssertThrow(false, ExcNotImplemented());
      std::array<NumberType, 3> res;
      return res;
    }



    template <typename NumberType>
    std::array<NumberType, 3>
    hyperbolic_rotation(const NumberType &f, const NumberType &g)
    {
      Assert(f != 0, ExcDivideByZero());
      const NumberType tau = g / f;
      AssertThrow(std::abs(tau) < 1.,
                  ExcMessage(
                    "real-valued Hyperbolic rotation does not exist for (" +
                    std::to_string(f) + "," + std::to_string(g) + ")"));
      const NumberType u =
        std::copysign(std::sqrt((1. - tau) * (1. + tau)),
                      f); // <-- more stable than std::sqrt(1.-tau*tau)
      std::array<NumberType, 3> csr;
      csr[0] = 1. / u;       // c
      csr[1] = csr[0] * tau; // s
      csr[2] = f * u;        // r
      return csr;
    }



    template <typename NumberType>
    std::array<std::complex<NumberType>, 3>
    givens_rotation(const std::complex<NumberType> & /*f*/,
                    const std::complex<NumberType> & /*g*/)
    {
      AssertThrow(false, ExcNotImplemented());
      std::array<NumberType, 3> res;
      return res;
    }



    template <typename NumberType>
    std::array<NumberType, 3>
    givens_rotation(const NumberType &f, const NumberType &g)
    {
      std::array<NumberType, 3> res;
      // naive calculation for "r" may overflow or underflow:
      // c =  x / \sqrt{x^2+y^2}
      // s = -y / \sqrt{x^2+y^2}

      // See Golub 2013, Matrix computations, Chapter 5.1.8
      // Algorithm 5.1.3
      // and
      // Anderson (2000),
      // Discontinuous Plane Rotations and the Symmetric Eigenvalue Problem.
      // LAPACK Working Note 150, University of Tennessee, UT-CS-00-454,
      // December 4, 2000.
      // Algorithm 4
      // We implement the latter below:
      if (g == NumberType())
        {
          res[0] = std::copysign(1., f);
          res[1] = NumberType();
          res[2] = std::abs(f);
        }
      else if (f == NumberType())
        {
          res[0] = NumberType();
          res[1] = std::copysign(1., g);
          res[2] = std::abs(g);
        }
      else if (std::abs(f) > std::abs(g))
        {
          const NumberType tau = g / f;
          const NumberType u   = std::copysign(std::sqrt(1. + tau * tau), f);
          res[0]               = 1. / u;       // c
          res[1]               = res[0] * tau; // s
          res[2]               = f * u;        // r
        }
      else
        {
          const NumberType tau = f / g;
          const NumberType u   = std::copysign(std::sqrt(1. + tau * tau), g);
          res[1]               = 1. / u;       // s
          res[0]               = res[1] * tau; // c
          res[2]               = g * u;        // r
        }

      return res;
    }



    template <typename OperatorType, typename VectorType>
    double
    lanczos_largest_eigenvalue(const OperatorType &      H,
                               const VectorType &        v0_,
                               const unsigned int        k,
                               VectorMemory<VectorType> &vector_memory,
                               std::vector<double> *     eigenvalues)
    {
      // Do k-step Lanczos:

      typename VectorMemory<VectorType>::Pointer v(vector_memory);
      typename VectorMemory<VectorType>::Pointer v0(vector_memory);
      typename VectorMemory<VectorType>::Pointer f(vector_memory);

      v->reinit(v0_);
      v0->reinit(v0_);
      f->reinit(v0_);

      // two vectors to store diagonal and subdiagonal of the Lanczos
      // matrix
      std::vector<double> diagonal;
      std::vector<double> subdiagonal;

      // 1. Normalize input vector
      (*v)     = v0_;
      double a = v->l2_norm();
      Assert(a != 0, ExcDivideByZero());
      (*v) *= 1. / a;

      // 2. Compute f = Hv; a = f*v; f <- f - av; T(0,0)=a;
      H.vmult(*f, *v);
      a = (*f) * (*v);
      f->add(-a, *v);
      diagonal.push_back(a);

      // 3. Loop over steps
      for (unsigned int i = 1; i < k; ++i)
        {
          // 4. L2 norm of f
          const double b = f->l2_norm();
          Assert(b != 0, ExcDivideByZero());
          // 5. v0 <- v; v <- f/b
          *v0 = *v;
          *v  = *f;
          (*v) *= 1. / b;
          // 6. f = Hv; f <- f - b v0;
          H.vmult(*f, *v);
          f->add(-b, *v0);
          // 7. a = f*v; f <- f - a v;
          a = (*f) * (*v);
          f->add(-a, *v);
          // 8. T(i,i-1) = T(i-1,i) = b;  T(i,i) = a;
          diagonal.push_back(a);
          subdiagonal.push_back(b);
        }

      Assert(diagonal.size() == k, ExcInternalError());
      Assert(subdiagonal.size() == k - 1, ExcInternalError());

      // Use Lapack dstev to get ||T||_2 norm, i.e. the largest eigenvalue
      // of T
      const types::blas_int n = k;
      std::vector<double>   Z;       // unused for eigenvalues-only ("N") job
      const types::blas_int ldz = 1; // ^^   (>=1)
      std::vector<double>   work;    // ^^
      types::blas_int       info;
      // call lapack_templates.h wrapper:
      internal::UtilitiesImplementation::call_stev('N',
                                                   n,
                                                   diagonal.data(),
                                                   subdiagonal.data(),
                                                   Z.data(),
                                                   ldz,
                                                   work.data(),
                                                   &info);

      Assert(info == 0, LAPACKSupport::ExcErrorCode("dstev", info));

      if (eigenvalues != nullptr)
        {
          eigenvalues->resize(diagonal.size());
          std::copy(diagonal.begin(), diagonal.end(), eigenvalues->begin());
        }

      // note that the largest eigenvalue of T is below the largest
      // eigenvalue of the operator.
      // return ||T||_2 + ||f||_2, although it is not guaranteed to be an upper
      // bound.
      return diagonal[k - 1] + f->l2_norm();
    }


    template <typename OperatorType, typename VectorType>
    void
    chebyshev_filter(VectorType &                    x,
                     const OperatorType &            op,
                     const unsigned int              degree,
                     const std::pair<double, double> unwanted_spectrum,
                     const double                    a_L,
                     VectorMemory<VectorType> &      vector_memory)
    {
      const double a = unwanted_spectrum.first;
      const double b = unwanted_spectrum.second;
      Assert(degree > 0, ExcMessage("Only positive degrees make sense."));

      const bool scale = (a_L < std::numeric_limits<double>::infinity());
      Assert(
        a < b,
        ExcMessage(
          "Lower bound of the unwanted spectrum should be smaller than the upper bound."));

      Assert(a_L <= a || a_L >= b || !scale,
             ExcMessage(
               "Scaling point should be outside of the unwanted spectrum."));

      // Setup auxiliary vectors:
      typename VectorMemory<VectorType>::Pointer p_y(vector_memory);
      typename VectorMemory<VectorType>::Pointer p_yn(vector_memory);

      p_y->reinit(x);
      p_yn->reinit(x);

      // convenience to avoid pointers
      VectorType &y  = *p_y;
      VectorType &yn = *p_yn;

      // Below is an implementation of
      // Algorithm 3.2 in Zhou et al, Journal of Computational Physics 274
      // (2014) 770-782 with **a bugfix for sigma1**. Here is the original
      // algorithm verbatim:
      //
      // [Y]=chebyshev_filter_scaled(X, m, a, b, aL).
      // e=(b-a)/2; c=(a+b)/2; σ=e/(c-aL); τ=2/σ;
      // Y=(H∗X-c∗X)∗(σ/e);
      // for i=2 to m do
      //   σnew =1/(τ - σ);
      //   Yt =(H∗Y - c∗Y)∗(2∗σnew/e)-(σ∗σnew)∗X;
      //   X =Y; Y =Yt; σ =σnew;

      const double e     = (b - a) / 2.;
      const double c     = (a + b) / 2.;
      const double alpha = 1. / e;
      const double beta  = -c / e;

      const double sigma1 =
        e / (a_L - c); // BUGFIX which is relevant for odd degrees
      double       sigma = scale ? sigma1 : 1.;
      const double tau   = 2. / sigma;
      op.vmult(y, x);
      y.sadd(alpha * sigma, beta * sigma, x);

      for (unsigned int i = 2; i <= degree; ++i)
        {
          const double sigma_new = scale ? 1. / (tau - sigma) : 1.;
          op.vmult(yn, y);
          yn.sadd(2. * alpha * sigma_new, 2. * beta * sigma_new, y);
          yn.add(-sigma * sigma_new, x);
          x.swap(y);
          y.swap(yn);
          sigma = sigma_new;
        }

      x.swap(y);
    }

  } // namespace LinearAlgebra
} // namespace Utilities

#endif



DEAL_II_NAMESPACE_CLOSE


#endif
