// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2021 by the deal.II authors
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

#ifndef dealii_solver_gmres_h
#define dealii_solver_gmres_h



#include <deal.II/base/config.h>

#include <deal.II/base/logstream.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/lac/block_vector_base.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/householder.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/solver.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/vector.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <vector>

DEAL_II_NAMESPACE_OPEN

// forward declarations
#ifndef DOXYGEN
namespace LinearAlgebra
{
  namespace distributed
  {
    template <typename, typename>
    class Vector;
  } // namespace distributed
} // namespace LinearAlgebra
#endif

/**
 * @addtogroup Solvers
 * @{
 */

namespace internal
{
  /**
   * A namespace for a helper class to the GMRES solver.
   */
  namespace SolverGMRESImplementation
  {
    /**
     * Class to hold temporary vectors.  This class automatically allocates a
     * new vector, once it is needed.
     *
     * A future version should also be able to shift through vectors
     * automatically, avoiding restart.
     */

    template <typename VectorType>
    class TmpVectors
    {
    public:
      /**
       * Constructor. Prepares an array of @p VectorType of length @p
       * max_size.
       */
      TmpVectors(const unsigned int max_size, VectorMemory<VectorType> &vmem);

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

      /**
       * Return size of data vector. It is used in the solver to store
       * the Arnoldi vectors.
       */
      unsigned int
      size() const;


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
  } // namespace SolverGMRESImplementation
} // namespace internal

/**
 * Implementation of the Restarted Preconditioned Direct Generalized Minimal
 * Residual Method. The stopping criterion is the norm of the residual.
 *
 * The AdditionalData structure contains the number of temporary vectors used.
 * The size of the Arnoldi basis is this number minus three. Additionally, it
 * allows you to choose between right or left preconditioning. The default is
 * left preconditioning. Finally it includes a flag indicating whether or not
 * the default residual is used as stopping criterion.
 *
 *
 * <h3>Left versus right preconditioning</h3>
 *
 * @p AdditionalData allows you to choose between left and right
 * preconditioning. As expected, this switches between solving for the systems
 * <i>P<sup>-1</sup>A</i> and <i>AP<sup>-1</sup></i>, respectively.
 *
 * A second consequence is the type of residual which is used to measure
 * convergence. With left preconditioning, this is the <b>preconditioned</b>
 * residual, while with right preconditioning, it is the residual of the
 * unpreconditioned system.
 *
 * Optionally, this behavior can be overridden by using the flag
 * AdditionalData::use_default_residual. A <tt>true</tt> value refers to the
 * behavior described in the previous paragraph, while <tt>false</tt> reverts
 * it. Be aware though that additional residuals have to be computed in this
 * case, impeding the overall performance of the solver.
 *
 *
 * <h3>The size of the Arnoldi basis</h3>
 *
 * The maximal basis size is controlled by AdditionalData::max_n_tmp_vectors,
 * and it is this number minus 2. If the number of iteration steps exceeds
 * this number, all basis vectors are discarded and the iteration starts anew
 * from the approximation obtained so far.
 *
 * Note that the minimizing property of GMRes only pertains to the Krylov
 * space spanned by the Arnoldi basis. Therefore, restarted GMRes is
 * <b>not</b> minimizing anymore. The choice of the basis length is a trade-
 * off between memory consumption and convergence speed, since a longer basis
 * means minimization over a larger space.
 *
 * For the requirements on matrices and vectors in order to work with this
 * class, see the documentation of the Solver base class.
 *
 *
 * <h3>Observing the progress of linear solver iterations</h3>
 *
 * The solve() function of this class uses the mechanism described in the
 * Solver base class to determine convergence. This mechanism can also be used
 * to observe the progress of the iteration.
 *
 *
 * <h3>Eigenvalue and condition number estimates</h3>
 *
 * This class can estimate eigenvalues and condition number during the
 * solution process. This is done by creating the Hessenberg matrix during the
 * inner iterations. The eigenvalues are estimated as the eigenvalues of the
 * Hessenberg matrix and the condition number is estimated as the ratio of the
 * largest and smallest singular value of the Hessenberg matrix. The estimates
 * can be obtained by connecting a function as a slot using @p
 * connect_condition_number_slot and @p connect_eigenvalues_slot. These slots
 * will then be called from the solver with the estimates as argument.
 */
template <class VectorType = Vector<double>>
class SolverGMRES : public SolverBase<VectorType>
{
public:
  /**
   * Standardized data struct to pipe additional data to the solver.
   */
  struct AdditionalData
  {
    enum class OrthogonalizationStrategy
    {
      /**
       * Use modified Gram-Schmidt algorithm.
       */
      modified_gram_schmidt,
      /**
       * Use classical Gram-Schmidt algorithm. Since this approach works on
       * multi-vectors and performs a global reduction only once, it is
       * more efficient than the modified Gram-Schmidt algorithm.
       * However, it might be numerically unstable.
       */
      classical_gram_schmidt
    };

    /**
     * Constructor. By default, set the number of temporary vectors to 30,
     * i.e. do a restart every 28 iterations. Also set preconditioning from
     * left, the residual of the stopping criterion to the default residual,
     * and re-orthogonalization only if necessary. Also, the batched mode with
     * reduced functionality to track information is disabled by default.
     */
    explicit AdditionalData(
      const unsigned int              max_n_tmp_vectors          = 30,
      const bool                      right_preconditioning      = false,
      const bool                      use_default_residual       = true,
      const bool                      force_re_orthogonalization = false,
      const bool                      batched_mode               = false,
      const OrthogonalizationStrategy orthogonalization_strategy =
        OrthogonalizationStrategy::modified_gram_schmidt);

    /**
     * Maximum number of temporary vectors. This parameter controls the size
     * of the Arnoldi basis, which for historical reasons is
     * #max_n_tmp_vectors-2. SolverGMRES assumes that there are at least three
     * temporary vectors, so this value must be greater than or equal to three.
     */
    unsigned int max_n_tmp_vectors;

    /**
     * Flag for right preconditioning.
     *
     * @note Change between left and right preconditioning will also change
     * the way residuals are evaluated. See the corresponding section in the
     * SolverGMRES.
     */
    bool right_preconditioning;

    /**
     * Flag for the default residual that is used to measure convergence.
     */
    bool use_default_residual;

    /**
     * Flag to force re-orthogonalization of orthonormal basis in every step.
     * If set to false, the solver automatically checks for loss of
     * orthogonality every 5 iterations and enables re-orthogonalization only
     * if necessary.
     */
    bool force_re_orthogonalization;

    /**
     * Flag to control whether a reduced mode of the solver should be
     * run. This is necessary when running (several) SolverGMRES instances
     * involving very small and cheap linear systems where the feedback from
     * all signals, eigenvalue computations, and log stream are disabled.
     */
    bool batched_mode;

    /**
     * Strategy to orthogonalize vectors.
     */
    OrthogonalizationStrategy orthogonalization_strategy;
  };

  /**
   * Constructor.
   */
  SolverGMRES(SolverControl &           cn,
              VectorMemory<VectorType> &mem,
              const AdditionalData &    data = AdditionalData());

  /**
   * Constructor. Use an object of type GrowingVectorMemory as a default to
   * allocate memory.
   */
  SolverGMRES(SolverControl &cn, const AdditionalData &data = AdditionalData());

  /**
   * The copy constructor is deleted.
   */
  SolverGMRES(const SolverGMRES<VectorType> &) = delete;

  /**
   * Solve the linear system $Ax=b$ for x.
   */
  template <typename MatrixType, typename PreconditionerType>
  void
  solve(const MatrixType &        A,
        VectorType &              x,
        const VectorType &        b,
        const PreconditionerType &preconditioner);

  /**
   * Connect a slot to retrieve the estimated condition number. Called on each
   * outer iteration if every_iteration=true, otherwise called once when
   * iterations are ended (i.e., either because convergence has been achieved,
   * or because divergence has been detected).
   */
  boost::signals2::connection
  connect_condition_number_slot(const std::function<void(double)> &slot,
                                const bool every_iteration = false);

  /**
   * Connect a slot to retrieve the estimated eigenvalues. Called on each
   * outer iteration if every_iteration=true, otherwise called once when
   * iterations are ended (i.e., either because convergence has been achieved,
   * or because divergence has been detected).
   */
  boost::signals2::connection
  connect_eigenvalues_slot(
    const std::function<void(const std::vector<std::complex<double>> &)> &slot,
    const bool every_iteration = false);

  /**
   * Connect a slot to retrieve the Hessenberg matrix obtained by the
   * projection of the initial matrix on the Krylov basis. Called on each
   * outer iteration if every_iteration=true, otherwise called once when
   * iterations are ended (i.e., either because convergence has been achieved,
   * or because divergence has been detected).
   */
  boost::signals2::connection
  connect_hessenberg_slot(
    const std::function<void(const FullMatrix<double> &)> &slot,
    const bool every_iteration = true);

  /**
   * Connect a slot to retrieve the basis vectors of the Krylov space
   * generated by the Arnoldi algorithm. Called at once when iterations
   * are completed (i.e., either because convergence has been achieved,
   * or because divergence has been detected).
   */
  boost::signals2::connection
  connect_krylov_space_slot(
    const std::function<
      void(const internal::SolverGMRESImplementation::TmpVectors<VectorType> &)>
      &slot);


  /**
   * Connect a slot to retrieve a notification when the vectors are
   * re-orthogonalized.
   */
  boost::signals2::connection
  connect_re_orthogonalization_slot(const std::function<void(int)> &slot);


  DeclException1(ExcTooFewTmpVectors,
                 int,
                 << "The number of temporary vectors you gave (" << arg1
                 << ") is too small. It should be at least 10 for "
                 << "any results, and much more for reasonable ones.");

protected:
  /**
   * Includes the maximum number of tmp vectors.
   */
  AdditionalData additional_data;

  /**
   * Signal used to retrieve the estimated condition number. Called once when
   * all iterations are ended.
   */
  boost::signals2::signal<void(double)> condition_number_signal;

  /**
   * Signal used to retrieve the estimated condition numbers. Called on each
   * outer iteration.
   */
  boost::signals2::signal<void(double)> all_condition_numbers_signal;

  /**
   * Signal used to retrieve the estimated eigenvalues. Called once when all
   * iterations are ended.
   */
  boost::signals2::signal<void(const std::vector<std::complex<double>> &)>
    eigenvalues_signal;

  /**
   * Signal used to retrieve the estimated eigenvalues. Called on each outer
   * iteration.
   */
  boost::signals2::signal<void(const std::vector<std::complex<double>> &)>
    all_eigenvalues_signal;

  /**
   * Signal used to retrieve the Hessenberg matrix. Called once when
   * all iterations are ended.
   */
  boost::signals2::signal<void(const FullMatrix<double> &)> hessenberg_signal;

  /**
   * Signal used to retrieve the Hessenberg matrix. Called on each outer
   * iteration.
   */
  boost::signals2::signal<void(const FullMatrix<double> &)>
    all_hessenberg_signal;

  /**
   * Signal used to retrieve the Krylov space basis vectors. Called once
   * when all iterations are ended.
   */
  boost::signals2::signal<void(
    const internal::SolverGMRESImplementation::TmpVectors<VectorType> &)>
    krylov_space_signal;

  /**
   * Signal used to retrieve a notification
   * when the vectors are re-orthogonalized.
   */
  boost::signals2::signal<void(int)> re_orthogonalize_signal;

  /**
   * A reference to the underlying SolverControl object. In the regular case,
   * this is not needed, as the signal from the base class is used, but the
   * batched variant cannot use those mechanisms due to the high costs.
   */
  SolverControl &solver_control;

  /**
   * Implementation of the computation of the norm of the residual.
   */
  virtual double
  criterion();

  /**
   * Transformation of an upper Hessenberg matrix into tridiagonal structure
   * by givens rotation of the last column
   */
  void
  givens_rotation(Vector<double> &h,
                  Vector<double> &b,
                  Vector<double> &ci,
                  Vector<double> &si,
                  int             col) const;

  /**
   * Estimates the eigenvalues from the Hessenberg matrix, H_orig, generated
   * during the inner iterations. Uses these estimate to compute the condition
   * number. Calls the signals eigenvalues_signal and cond_signal with these
   * estimates as arguments.
   */
  static void
  compute_eigs_and_cond(
    const FullMatrix<double> &H_orig,
    const unsigned int        dim,
    const boost::signals2::signal<
      void(const std::vector<std::complex<double>> &)> &eigenvalues_signal,
    const boost::signals2::signal<void(const FullMatrix<double> &)>
      &                                          hessenberg_signal,
    const boost::signals2::signal<void(double)> &cond_signal);

  /**
   * Projected system matrix
   */
  FullMatrix<double> H;

  /**
   * Auxiliary vector for orthogonalization
   */
  Vector<double> gamma;

  /**
   * Auxiliary vector for orthogonalization
   */
  Vector<double> ci;

  /**
   * Auxiliary vector for orthogonalization
   */
  Vector<double> si;

  /**
   * Auxiliary vector for orthogonalization
   */
  Vector<double> h;
};



/**
 * Implementation of the Generalized minimal residual method with flexible
 * preconditioning (flexible GMRES or FGMRES).
 *
 * This flexible version of the GMRES method allows for the use of a different
 * preconditioner in each iteration step. Therefore, it is also more robust
 * with respect to inaccurate evaluation of the preconditioner. An important
 * application is the use of a Krylov space method inside the
 * preconditioner. As opposed to SolverGMRES which allows one to choose
 * between left and right preconditioning, this solver always applies the
 * preconditioner from the right.
 *
 * FGMRES needs two vectors in each iteration steps yielding a total of
 * <tt>2*SolverFGMRES::%AdditionalData::%max_basis_size+1</tt> auxiliary
 * vectors. Otherwise, FGMRES requires roughly the same number of operations
 * per iteration compared to GMRES, except one application of the
 * preconditioner less at each restart and at the end of solve().
 *
 * For more details see @cite Saad1991.
 */
template <class VectorType = Vector<double>>
class SolverFGMRES : public SolverBase<VectorType>
{
public:
  /**
   * Standardized data struct to pipe additional data to the solver.
   */
  struct AdditionalData
  {
    /**
     * Constructor. By default, set the maximum basis size to 30.
     */
    explicit AdditionalData(const unsigned int max_basis_size = 30)
      : max_basis_size(max_basis_size)
    {}

    /**
     * Maximum basis size.
     */
    unsigned int max_basis_size;
  };

  /**
   * Constructor.
   */
  SolverFGMRES(SolverControl &           cn,
               VectorMemory<VectorType> &mem,
               const AdditionalData &    data = AdditionalData());

  /**
   * Constructor. Use an object of type GrowingVectorMemory as a default to
   * allocate memory.
   */
  SolverFGMRES(SolverControl &       cn,
               const AdditionalData &data = AdditionalData());

  /**
   * Solve the linear system $Ax=b$ for x.
   */
  template <typename MatrixType, typename PreconditionerType>
  void
  solve(const MatrixType &        A,
        VectorType &              x,
        const VectorType &        b,
        const PreconditionerType &preconditioner);

private:
  /**
   * Additional flags.
   */
  AdditionalData additional_data;

  /**
   * Projected system matrix
   */
  FullMatrix<double> H;

  /**
   * Auxiliary matrix for inverting @p H
   */
  FullMatrix<double> H1;
};

/** @} */
/* --------------------- Inline and template functions ------------------- */


#ifndef DOXYGEN
namespace internal
{
  namespace SolverGMRESImplementation
  {
    template <class VectorType>
    inline TmpVectors<VectorType>::TmpVectors(const unsigned int max_size,
                                              VectorMemory<VectorType> &vmem)
      : mem(vmem)
      , data(max_size)
    {}



    template <class VectorType>
    inline VectorType &
    TmpVectors<VectorType>::operator[](const unsigned int i) const
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
          data[i]->reinit(temp, true);
        }
      return *data[i];
    }



    template <class VectorType>
    unsigned int
    TmpVectors<VectorType>::size() const
    {
      return (data.size() > 0 ? data.size() - 1 : 0);
    }



    // A comparator for better printing eigenvalues
    inline bool
    complex_less_pred(const std::complex<double> &x,
                      const std::complex<double> &y)
    {
      return x.real() < y.real() ||
             (x.real() == y.real() && x.imag() < y.imag());
    }

    // A function to solve the (upper) triangular system after Givens
    // rotations on a matrix that has possibly unused rows and columns
    inline void
    solve_triangular(const unsigned int        dim,
                     const FullMatrix<double> &H,
                     const Vector<double> &    rhs,
                     Vector<double> &          solution)
    {
      for (int i = dim - 1; i >= 0; --i)
        {
          double s = rhs(i);
          for (unsigned int j = i + 1; j < dim; ++j)
            s -= solution(j) * H(i, j);
          solution(i) = s / H(i, i);
          AssertIsFinite(solution(i));
        }
    }
  } // namespace SolverGMRESImplementation
} // namespace internal



template <class VectorType>
inline SolverGMRES<VectorType>::AdditionalData::AdditionalData(
  const unsigned int              max_n_tmp_vectors,
  const bool                      right_preconditioning,
  const bool                      use_default_residual,
  const bool                      force_re_orthogonalization,
  const bool                      batched_mode,
  const OrthogonalizationStrategy orthogonalization_strategy)
  : max_n_tmp_vectors(max_n_tmp_vectors)
  , right_preconditioning(right_preconditioning)
  , use_default_residual(use_default_residual)
  , force_re_orthogonalization(force_re_orthogonalization)
  , batched_mode(batched_mode)
  , orthogonalization_strategy(orthogonalization_strategy)
{
  Assert(3 <= max_n_tmp_vectors,
         ExcMessage("SolverGMRES needs at least three "
                    "temporary vectors."));
}



template <class VectorType>
SolverGMRES<VectorType>::SolverGMRES(SolverControl &           cn,
                                     VectorMemory<VectorType> &mem,
                                     const AdditionalData &    data)
  : SolverBase<VectorType>(cn, mem)
  , additional_data(data)
  , solver_control(cn)
{}



template <class VectorType>
SolverGMRES<VectorType>::SolverGMRES(SolverControl &       cn,
                                     const AdditionalData &data)
  : SolverBase<VectorType>(cn)
  , additional_data(data)
  , solver_control(cn)
{}



template <class VectorType>
inline void
SolverGMRES<VectorType>::givens_rotation(Vector<double> &h,
                                         Vector<double> &b,
                                         Vector<double> &ci,
                                         Vector<double> &si,
                                         int             col) const
{
  for (int i = 0; i < col; ++i)
    {
      const double s     = si(i);
      const double c     = ci(i);
      const double dummy = h(i);
      h(i)               = c * dummy + s * h(i + 1);
      h(i + 1)           = -s * dummy + c * h(i + 1);
    };

  const double r = 1. / std::sqrt(h(col) * h(col) + h(col + 1) * h(col + 1));
  si(col)        = h(col + 1) * r;
  ci(col)        = h(col) * r;
  h(col)         = ci(col) * h(col) + si(col) * h(col + 1);
  b(col + 1)     = -si(col) * b(col);
  b(col) *= ci(col);
}



namespace internal
{
  namespace SolverGMRESImplementation
  {
    template <typename VectorType, typename Enable = void>
    struct is_dealii_compatible_distributed_vector;

    template <typename VectorType>
    struct is_dealii_compatible_distributed_vector<
      VectorType,
      typename std::enable_if<!internal::is_block_vector<VectorType>>::type>
    {
      static constexpr bool value = std::is_same<
        VectorType,
        LinearAlgebra::distributed::Vector<typename VectorType::value_type,
                                           MemorySpace::Host>>::value;
    };



    template <typename VectorType>
    struct is_dealii_compatible_distributed_vector<
      VectorType,
      typename std::enable_if<internal::is_block_vector<VectorType>>::type>
    {
      static constexpr bool value = std::is_same<
        typename VectorType::BlockType,
        LinearAlgebra::distributed::Vector<typename VectorType::value_type,
                                           MemorySpace::Host>>::value;
    };



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
      (void)b;
      return vector;
    }



    template <typename VectorType,
              std::enable_if_t<!IsBlockVector<VectorType>::value, VectorType>
                * = nullptr>
    const VectorType &
    block(const VectorType &vector, const unsigned int b)
    {
      AssertDimension(b, 0);
      (void)b;
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



    template <typename VectorType,
              std::enable_if_t<IsBlockVector<VectorType>::value, VectorType> * =
                nullptr>
    const typename VectorType::BlockType &
    block(const VectorType &vector, const unsigned int b)
    {
      return vector.block(b);
    }



    template <class VectorType,
              std::enable_if_t<
                !is_dealii_compatible_distributed_vector<VectorType>::value,
                VectorType> * = nullptr>
    void
    Tvmult_add(const unsigned int dim,
               const VectorType & vv,
               const internal::SolverGMRESImplementation::TmpVectors<VectorType>
                 &             orthogonal_vectors,
               Vector<double> &h)
    {
      for (unsigned int i = 0; i < dim; ++i)
        h[i] += vv * orthogonal_vectors[i];
    }



    template <class VectorType,
              std::enable_if_t<
                is_dealii_compatible_distributed_vector<VectorType>::value,
                VectorType> * = nullptr>
    void
    Tvmult_add(const unsigned int dim,
               const VectorType & vv,
               const internal::SolverGMRESImplementation::TmpVectors<VectorType>
                 &             orthogonal_vectors,
               Vector<double> &h)
    {
      for (unsigned int b = 0; b < n_blocks(vv); ++b)
        {
          unsigned int j = 0;

          if (dim <= 128)
            {
              // optimized path
              static constexpr unsigned int n_lanes =
                VectorizedArray<double>::size();

              VectorizedArray<double> hs[128];
              for (unsigned int d = 0; d < dim; ++d)
                hs[d] = 0.0;

              unsigned int c = 0;

              for (; c < block(vv, b).locally_owned_size() / n_lanes / 4;
                   ++c, j += n_lanes * 4)
                for (unsigned int i = 0; i < dim; ++i)
                  {
                    VectorizedArray<double> vvec[4];
                    for (unsigned int k = 0; k < 4; ++k)
                      vvec[k].load(block(vv, b).begin() + j + k * n_lanes);

                    for (unsigned int k = 0; k < 4; ++k)
                      {
                        VectorizedArray<double> temp;
                        temp.load(block(orthogonal_vectors[i], b).begin() + j +
                                  k * n_lanes);
                        hs[i] += temp * vvec[k];
                      }
                  }

              c *= 4;
              for (; c < block(vv, b).locally_owned_size() / n_lanes;
                   ++c, j += n_lanes)
                for (unsigned int i = 0; i < dim; ++i)
                  {
                    VectorizedArray<double> vvec, temp;
                    vvec.load(block(vv, b).begin() + j);
                    temp.load(block(orthogonal_vectors[i], b).begin() + j);
                    hs[i] += temp * vvec;
                  }

              for (unsigned int i = 0; i < dim; ++i)
                for (unsigned int v = 0; v < n_lanes; ++v)
                  h(i) += hs[i][v];
            }

          // remainder loop of optimized path or non-optimized path (if
          // dim>128)
          for (; j < block(vv, b).locally_owned_size(); ++j)
            for (unsigned int i = 0; i < dim; ++i)
              h(i) += block(orthogonal_vectors[i], b).local_element(j) *
                      block(vv, b).local_element(j);
        }

      Utilities::MPI::sum(h, MPI_COMM_WORLD, h);
    }



    template <class VectorType,
              std::enable_if_t<
                !is_dealii_compatible_distributed_vector<VectorType>::value,
                VectorType> * = nullptr>
    double
    substract_and_norm(
      const unsigned int dim,
      const internal::SolverGMRESImplementation::TmpVectors<VectorType>
        &                   orthogonal_vectors,
      const Vector<double> &h,
      VectorType &          vv)
    {
      Assert(dim > 0, ExcInternalError());

      for (unsigned int i = 0; i < dim; ++i)
        vv.add(-h(i), orthogonal_vectors[i]);

      return std::sqrt(vv.add_and_dot(-h(dim), orthogonal_vectors[dim], vv));
    }



    template <class VectorType,
              std::enable_if_t<
                is_dealii_compatible_distributed_vector<VectorType>::value,
                VectorType> * = nullptr>
    double
    substract_and_norm(
      const unsigned int dim,
      const internal::SolverGMRESImplementation::TmpVectors<VectorType>
        &                   orthogonal_vectors,
      const Vector<double> &h,
      VectorType &          vv)
    {
      static constexpr unsigned int n_lanes = VectorizedArray<double>::size();

      double norm_vv_temp = 0.0;

      for (unsigned int b = 0; b < n_blocks(vv); ++b)
        {
          VectorizedArray<double> norm_vv_temp_vectorized = 0.0;

          unsigned int j = 0;
          unsigned int c = 0;
          for (; c < block(vv, b).locally_owned_size() / n_lanes / 4;
               ++c, j += n_lanes * 4)
            {
              VectorizedArray<double> temp[4];

              for (unsigned int k = 0; k < 4; ++k)
                temp[k].load(block(vv, b).begin() + j + k * n_lanes);

              for (unsigned int i = 0; i < dim; ++i)
                {
                  const double factor = h(i);
                  for (unsigned int k = 0; k < 4; ++k)
                    {
                      VectorizedArray<double> vec;
                      vec.load(block(orthogonal_vectors[i], b).begin() + j +
                               k * n_lanes);
                      temp[k] -= factor * vec;
                    }
                }

              for (unsigned int k = 0; k < 4; ++k)
                temp[k].store(block(vv, b).begin() + j + k * n_lanes);

              norm_vv_temp_vectorized +=
                (temp[0] * temp[0] + temp[1] * temp[1]) +
                (temp[2] * temp[2] + temp[3] * temp[3]);
            }

          c *= 4;
          for (; c < block(vv, b).locally_owned_size() / n_lanes;
               ++c, j += n_lanes)
            {
              VectorizedArray<double> temp;
              temp.load(block(vv, b).begin() + j);

              for (unsigned int i = 0; i < dim; ++i)
                {
                  VectorizedArray<double> vec;
                  vec.load(block(orthogonal_vectors[i], b).begin() + j);
                  temp -= h(i) * vec;
                }

              temp.store(block(vv, b).begin() + j);

              norm_vv_temp_vectorized += temp * temp;
            }

          for (unsigned int v = 0; v < n_lanes; ++v)
            norm_vv_temp += norm_vv_temp_vectorized[v];

          for (; j < block(vv, b).locally_owned_size(); ++j)
            {
              double temp = block(vv, b).local_element(j);
              for (unsigned int i = 0; i < dim; ++i)
                temp -= h(i) * block(orthogonal_vectors[i], b).local_element(j);
              block(vv, b).local_element(j) = temp;

              norm_vv_temp += temp * temp;
            }
        }

      return std::sqrt(Utilities::MPI::sum(norm_vv_temp, MPI_COMM_WORLD));
    }


    template <class VectorType,
              std::enable_if_t<
                !is_dealii_compatible_distributed_vector<VectorType>::value,
                VectorType> * = nullptr>
    double
    sadd_and_norm(VectorType &      v,
                  const double      factor_a,
                  const VectorType &b,
                  const double      factor_b)
    {
      v.sadd(factor_a, factor_b, b);
      return v.l2_norm();
    }


    template <class VectorType,
              std::enable_if_t<
                is_dealii_compatible_distributed_vector<VectorType>::value,
                VectorType> * = nullptr>
    double
    sadd_and_norm(VectorType &      v,
                  const double      factor_a,
                  const VectorType &w,
                  const double      factor_b)
    {
      double norm = 0;

      for (unsigned int b = 0; b < n_blocks(v); ++b)
        for (unsigned int j = 0; j < block(v, b).locally_owned_size(); ++j)
          {
            const double temp = block(v, b).local_element(j) * factor_a +
                                block(w, b).local_element(j) * factor_b;

            block(v, b).local_element(j) = temp;

            norm += temp * temp;
          }

      return std::sqrt(Utilities::MPI::sum(norm, MPI_COMM_WORLD));
    }



    template <class VectorType,
              std::enable_if_t<
                !is_dealii_compatible_distributed_vector<VectorType>::value,
                VectorType> * = nullptr>
    void
    add(VectorType &          p,
        const unsigned int    dim,
        const Vector<double> &h,
        const internal::SolverGMRESImplementation::TmpVectors<VectorType>
          &        tmp_vectors,
        const bool zero_out)
    {
      if (zero_out)
        p.equ(h(0), tmp_vectors[0]);
      else
        p.add(h(0), tmp_vectors[0]);

      for (unsigned int i = 1; i < dim; ++i)
        p.add(h(i), tmp_vectors[i]);
    }



    template <class VectorType,
              std::enable_if_t<
                is_dealii_compatible_distributed_vector<VectorType>::value,
                VectorType> * = nullptr>
    void
    add(VectorType &          p,
        const unsigned int    dim,
        const Vector<double> &h,
        const internal::SolverGMRESImplementation::TmpVectors<VectorType>
          &        tmp_vectors,
        const bool zero_out)
    {
      for (unsigned int b = 0; b < n_blocks(p); ++b)
        for (unsigned int j = 0; j < block(p, b).locally_owned_size(); ++j)
          {
            double temp = zero_out ? 0 : block(p, b).local_element(j);
            for (unsigned int i = 0; i < dim; ++i)
              temp += block(tmp_vectors[i], b).local_element(j) * h(i);
            block(p, b).local_element(j) = temp;
          }
    }



    /**
     * Orthogonalize the vector @p vv against the @p dim (orthogonal) vectors
     * given by @p orthogonal_vectors using the modified or classical
     * Gram-Schmidt algorithm.
     * The factors used for orthogonalization are stored in @p h. The boolean @p
     * re_orthogonalize specifies whether the Gram-Schmidt algorithm
     * should be applied twice. The algorithm checks loss of orthogonality in
     * the procedure every fifth step and sets the flag to true in that case.
     * All subsequent iterations use re-orthogonalization.
     * Calls the signal re_orthogonalize_signal if it is connected.
     */
    template <class VectorType>
    inline double
    iterated_gram_schmidt(
      const typename SolverGMRES<VectorType>::AdditionalData::
        OrthogonalizationStrategy orthogonalization_strategy,
      const internal::SolverGMRESImplementation::TmpVectors<VectorType>
        &                                       orthogonal_vectors,
      const unsigned int                        dim,
      const unsigned int                        accumulated_iterations,
      VectorType &                              vv,
      Vector<double> &                          h,
      bool &                                    reorthogonalize,
      const boost::signals2::signal<void(int)> &reorthogonalize_signal =
        boost::signals2::signal<void(int)>())
    {
      Assert(dim > 0, ExcInternalError());
      const unsigned int inner_iteration = dim - 1;

      // need initial norm for detection of re-orthogonalization, see below
      double     norm_vv_start = 0;
      const bool consider_reorthogonalize =
        (reorthogonalize == false) && (inner_iteration % 5 == 4);
      if (consider_reorthogonalize)
        norm_vv_start = vv.l2_norm();

      for (unsigned int i = 0; i < dim; ++i)
        h[i] = 0;

      for (unsigned int c = 0; c < 2;
           ++c) // 0: orthogonalize, 1: reorthogonalize
        {
          // Orthogonalization
          double norm_vv = 0.0;

          if (orthogonalization_strategy ==
              SolverGMRES<VectorType>::AdditionalData::
                OrthogonalizationStrategy::modified_gram_schmidt)
            {
              double htmp = vv * orthogonal_vectors[0];
              h(0) += htmp;
              for (unsigned int i = 1; i < dim; ++i)
                {
                  htmp = vv.add_and_dot(-htmp,
                                        orthogonal_vectors[i - 1],
                                        orthogonal_vectors[i]);
                  h(i) += htmp;
                }

              norm_vv = std::sqrt(
                vv.add_and_dot(-htmp, orthogonal_vectors[dim - 1], vv));
            }
          else if (orthogonalization_strategy ==
                   SolverGMRES<VectorType>::AdditionalData::
                     OrthogonalizationStrategy::classical_gram_schmidt)
            {
              Tvmult_add(dim, vv, orthogonal_vectors, h);
              norm_vv = substract_and_norm(dim, orthogonal_vectors, h, vv);
            }
          else
            {
              AssertThrow(false, ExcNotImplemented());
            }

          if (c == 1)
            return norm_vv; // reorthogonalization already performed -> finished

          // Re-orthogonalization if loss of orthogonality detected. For the
          // test, use a strategy discussed in C. T. Kelley, Iterative Methods
          // for Linear and Nonlinear Equations, SIAM, Philadelphia, 1995:
          // Compare the norm of vv after orthogonalization with its norm when
          // starting the orthogonalization. If vv became very small (here: less
          // than the square root of the machine precision times 10), it is
          // almost in the span of the previous vectors, which indicates loss of
          // precision.
          if (consider_reorthogonalize)
            {
              if (norm_vv >
                  10. * norm_vv_start *
                    std::sqrt(std::numeric_limits<
                              typename VectorType::value_type>::epsilon()))
                return norm_vv;

              else
                {
                  reorthogonalize = true;
                  if (!reorthogonalize_signal.empty())
                    reorthogonalize_signal(accumulated_iterations);
                }
            }

          if (reorthogonalize == false)
            return norm_vv; // no reorthogonalization needed -> finished
        }

      AssertThrow(false, ExcInternalError());

      return 0.0;
    }
  } // namespace SolverGMRESImplementation
} // namespace internal



template <class VectorType>
inline void
SolverGMRES<VectorType>::compute_eigs_and_cond(
  const FullMatrix<double> &H_orig,
  const unsigned int        dim,
  const boost::signals2::signal<void(const std::vector<std::complex<double>> &)>
    &eigenvalues_signal,
  const boost::signals2::signal<void(const FullMatrix<double> &)>
    &                                          hessenberg_signal,
  const boost::signals2::signal<void(double)> &cond_signal)
{
  // Avoid copying the Hessenberg matrix if it isn't needed.
  if ((!eigenvalues_signal.empty() || !hessenberg_signal.empty() ||
       !cond_signal.empty()) &&
      dim > 0)
    {
      LAPACKFullMatrix<double> mat(dim, dim);
      for (unsigned int i = 0; i < dim; ++i)
        for (unsigned int j = 0; j < dim; ++j)
          mat(i, j) = H_orig(i, j);
      hessenberg_signal(H_orig);
      // Avoid computing eigenvalues if they are not needed.
      if (!eigenvalues_signal.empty())
        {
          // Copy mat so that we can compute svd below. Necessary since
          // compute_eigenvalues will leave mat in state
          // LAPACKSupport::unusable.
          LAPACKFullMatrix<double> mat_eig(mat);
          mat_eig.compute_eigenvalues();
          std::vector<std::complex<double>> eigenvalues(dim);
          for (unsigned int i = 0; i < mat_eig.n(); ++i)
            eigenvalues[i] = mat_eig.eigenvalue(i);
          // Sort eigenvalues for nicer output.
          std::sort(eigenvalues.begin(),
                    eigenvalues.end(),
                    internal::SolverGMRESImplementation::complex_less_pred);
          eigenvalues_signal(eigenvalues);
        }
      // Calculate condition number, avoid calculating the svd if a slot
      // isn't connected. Need at least a 2-by-2 matrix to do the estimate.
      if (!cond_signal.empty() && (mat.n() > 1))
        {
          mat.compute_svd();
          double condition_number =
            mat.singular_value(0) / mat.singular_value(mat.n() - 1);
          cond_signal(condition_number);
        }
    }
}



template <class VectorType>
template <typename MatrixType, typename PreconditionerType>
void
SolverGMRES<VectorType>::solve(const MatrixType &        A,
                               VectorType &              x,
                               const VectorType &        b,
                               const PreconditionerType &preconditioner)
{
  // TODO:[?] Check, why there are two different start residuals.
  // TODO:[GK] Make sure the parameter in the constructor means maximum basis
  // size

  std::unique_ptr<LogStream::Prefix> prefix;
  if (!additional_data.batched_mode)
    prefix = std::make_unique<LogStream::Prefix>("GMRES");

  // extra call to std::max to placate static analyzers: coverity rightfully
  // complains that data.max_n_tmp_vectors - 2 may overflow
  const unsigned int n_tmp_vectors =
    std::max(additional_data.max_n_tmp_vectors, 3u);

  // Generate an object where basis vectors are stored.
  internal::SolverGMRESImplementation::TmpVectors<VectorType> tmp_vectors(
    n_tmp_vectors, this->memory);

  // number of the present iteration; this
  // number is not reset to zero upon a
  // restart
  unsigned int accumulated_iterations = 0;

  const bool do_eigenvalues =
    !additional_data.batched_mode &&
    (!condition_number_signal.empty() ||
     !all_condition_numbers_signal.empty() || !eigenvalues_signal.empty() ||
     !all_eigenvalues_signal.empty() || !hessenberg_signal.empty() ||
     !all_hessenberg_signal.empty());
  // for eigenvalue computation, need to collect the Hessenberg matrix (before
  // applying Givens rotations)
  FullMatrix<double> H_orig;
  if (do_eigenvalues)
    H_orig.reinit(n_tmp_vectors, n_tmp_vectors - 1);

  // matrix used for the orthogonalization process later
  H.reinit(n_tmp_vectors, n_tmp_vectors - 1, /* omit_initialization */ true);

  // some additional vectors, also used in the orthogonalization
  gamma.reinit(n_tmp_vectors);
  ci.reinit(n_tmp_vectors - 1);
  si.reinit(n_tmp_vectors - 1);
  h.reinit(n_tmp_vectors - 1);

  unsigned int dim = 0;

  SolverControl::State iteration_state = SolverControl::iterate;
  double               last_res        = std::numeric_limits<double>::lowest();

  // switch to determine whether we want a left or a right preconditioner. at
  // present, left is default, but both ways are implemented
  const bool left_precondition = !additional_data.right_preconditioning;

  // Per default the left preconditioned GMRes uses the preconditioned
  // residual and the right preconditioned GMRes uses the unpreconditioned
  // residual as stopping criterion.
  const bool use_default_residual = additional_data.use_default_residual;

  // define two aliases
  VectorType &v = tmp_vectors(0, x);
  VectorType &p = tmp_vectors(n_tmp_vectors - 1, x);

  // Following vectors are needed when we are not using the default residuals
  // as stopping criterion
  typename VectorMemory<VectorType>::Pointer r;
  typename VectorMemory<VectorType>::Pointer x_;
  std::unique_ptr<dealii::Vector<double>>    gamma_;
  if (!use_default_residual)
    {
      r  = std::move(typename VectorMemory<VectorType>::Pointer(this->memory));
      x_ = std::move(typename VectorMemory<VectorType>::Pointer(this->memory));
      r->reinit(x);
      x_->reinit(x);

      gamma_ = std::make_unique<dealii::Vector<double>>(gamma.size());
    }

  bool re_orthogonalize = additional_data.force_re_orthogonalization;

  ///////////////////////////////////////////////////////////////////////////
  // outer iteration: loop until we either reach convergence or the maximum
  // number of iterations is exceeded. each cycle of this loop amounts to one
  // restart
  do
    {
      // reset this vector to the right size
      h.reinit(n_tmp_vectors - 1);

      double rho = 0.0;

      if (left_precondition)
        {
          A.vmult(p, x);
          p.sadd(-1., 1., b);
          preconditioner.vmult(v, p);
          rho = v.l2_norm();
        }
      else
        {
          A.vmult(v, x);
          rho = dealii::internal::SolverGMRESImplementation::sadd_and_norm(v,
                                                                           -1,
                                                                           b,
                                                                           1.0);
        }

      // check the residual here as well since it may be that we got the exact
      // (or an almost exact) solution vector at the outset. if we wouldn't
      // check here, the next scaling operation would produce garbage
      if (use_default_residual)
        {
          last_res = rho;
          if (additional_data.batched_mode)
            iteration_state = solver_control.check(accumulated_iterations, rho);
          else
            iteration_state =
              this->iteration_status(accumulated_iterations, rho, x);

          if (iteration_state != SolverControl::iterate)
            break;
        }
      else
        {
          deallog << "default_res=" << rho << std::endl;

          if (left_precondition)
            {
              A.vmult(*r, x);
              r->sadd(-1., 1., b);
            }
          else
            preconditioner.vmult(*r, v);

          double res = r->l2_norm();
          last_res   = res;
          if (additional_data.batched_mode)
            iteration_state = solver_control.check(accumulated_iterations, rho);
          else
            iteration_state =
              this->iteration_status(accumulated_iterations, res, x);

          if (iteration_state != SolverControl::iterate)
            break;
        }

      gamma(0) = rho;

      v *= 1. / rho;

      // inner iteration doing at most as many steps as there are temporary
      // vectors. the number of steps actually been done is propagated outside
      // through the @p dim variable
      for (unsigned int inner_iteration = 0;
           ((inner_iteration < n_tmp_vectors - 2) &&
            (iteration_state == SolverControl::iterate));
           ++inner_iteration)
        {
          ++accumulated_iterations;
          // yet another alias
          VectorType &vv = tmp_vectors(inner_iteration + 1, x);

          if (left_precondition)
            {
              A.vmult(p, tmp_vectors[inner_iteration]);
              preconditioner.vmult(vv, p);
            }
          else
            {
              preconditioner.vmult(p, tmp_vectors[inner_iteration]);
              A.vmult(vv, p);
            }

          dim = inner_iteration + 1;

          const double s =
            internal::SolverGMRESImplementation::iterated_gram_schmidt(
              additional_data.orthogonalization_strategy,
              tmp_vectors,
              dim,
              accumulated_iterations,
              vv,
              h,
              re_orthogonalize,
              re_orthogonalize_signal);
          h(inner_iteration + 1) = s;

          // s=0 is a lucky breakdown, the solver will reach convergence,
          // but we must not divide by zero here.
          if (s != 0)
            vv *= 1. / s;

          // for eigenvalues, get the resulting coefficients from the
          // orthogonalization process
          if (do_eigenvalues)
            for (unsigned int i = 0; i < dim + 1; ++i)
              H_orig(i, inner_iteration) = h(i);

          //  Transformation into tridiagonal structure
          givens_rotation(h, gamma, ci, si, inner_iteration);

          //  append vector on matrix
          for (unsigned int i = 0; i < dim; ++i)
            H(i, inner_iteration) = h(i);

          //  default residual
          rho = std::fabs(gamma(dim));

          if (use_default_residual)
            {
              last_res = rho;
              if (additional_data.batched_mode)
                iteration_state =
                  solver_control.check(accumulated_iterations, rho);
              else
                iteration_state =
                  this->iteration_status(accumulated_iterations, rho, x);
            }
          else
            {
              if (!additional_data.batched_mode)
                deallog << "default_res=" << rho << std::endl;

              *x_     = x;
              *gamma_ = gamma;
              internal::SolverGMRESImplementation::solve_triangular(dim,
                                                                    H,
                                                                    *gamma_,
                                                                    h);

              if (left_precondition)
                for (unsigned int i = 0; i < dim; ++i)
                  x_->add(h(i), tmp_vectors[i]);
              else
                {
                  p = 0.;
                  for (unsigned int i = 0; i < dim; ++i)
                    p.add(h(i), tmp_vectors[i]);
                  preconditioner.vmult(*r, p);
                  x_->add(1., *r);
                };
              A.vmult(*r, *x_);
              r->sadd(-1., 1., b);
              // Now *r contains the unpreconditioned residual!!
              if (left_precondition)
                {
                  const double res = r->l2_norm();
                  last_res         = res;

                  iteration_state =
                    this->iteration_status(accumulated_iterations, res, x);
                }
              else
                {
                  preconditioner.vmult(*x_, *r);
                  const double preconditioned_res = x_->l2_norm();
                  last_res                        = preconditioned_res;

                  if (additional_data.batched_mode)
                    iteration_state =
                      solver_control.check(accumulated_iterations, rho);
                  else
                    iteration_state =
                      this->iteration_status(accumulated_iterations,
                                             preconditioned_res,
                                             x);
                }
            }
        }

      // end of inner iteration. now calculate the solution from the temporary
      // vectors
      internal::SolverGMRESImplementation::solve_triangular(dim, H, gamma, h);

      if (do_eigenvalues)
        compute_eigs_and_cond(H_orig,
                              dim,
                              all_eigenvalues_signal,
                              all_hessenberg_signal,
                              condition_number_signal);

      if (left_precondition)
        dealii::internal::SolverGMRESImplementation::add(
          x, dim, h, tmp_vectors, false);
      else
        {
          dealii::internal::SolverGMRESImplementation::add(
            p, dim, h, tmp_vectors, true);
          preconditioner.vmult(v, p);
          x.add(1., v);
        };
      // end of outer iteration. restart if no convergence and the number of
      // iterations is not exceeded
    }
  while (iteration_state == SolverControl::iterate);

  if (do_eigenvalues)
    compute_eigs_and_cond(H_orig,
                          dim,
                          eigenvalues_signal,
                          hessenberg_signal,
                          condition_number_signal);

  if (!additional_data.batched_mode && !krylov_space_signal.empty())
    krylov_space_signal(tmp_vectors);

  // in case of failure: throw exception
  AssertThrow(iteration_state == SolverControl::success,
              SolverControl::NoConvergence(accumulated_iterations, last_res));
}



template <class VectorType>
boost::signals2::connection
SolverGMRES<VectorType>::connect_condition_number_slot(
  const std::function<void(double)> &slot,
  const bool                         every_iteration)
{
  if (every_iteration)
    {
      return all_condition_numbers_signal.connect(slot);
    }
  else
    {
      return condition_number_signal.connect(slot);
    }
}



template <class VectorType>
boost::signals2::connection
SolverGMRES<VectorType>::connect_eigenvalues_slot(
  const std::function<void(const std::vector<std::complex<double>> &)> &slot,
  const bool every_iteration)
{
  if (every_iteration)
    {
      return all_eigenvalues_signal.connect(slot);
    }
  else
    {
      return eigenvalues_signal.connect(slot);
    }
}



template <class VectorType>
boost::signals2::connection
SolverGMRES<VectorType>::connect_hessenberg_slot(
  const std::function<void(const FullMatrix<double> &)> &slot,
  const bool                                             every_iteration)
{
  if (every_iteration)
    {
      return all_hessenberg_signal.connect(slot);
    }
  else
    {
      return hessenberg_signal.connect(slot);
    }
}



template <class VectorType>
boost::signals2::connection
SolverGMRES<VectorType>::connect_krylov_space_slot(
  const std::function<void(
    const internal::SolverGMRESImplementation::TmpVectors<VectorType> &)> &slot)
{
  return krylov_space_signal.connect(slot);
}



template <class VectorType>
boost::signals2::connection
SolverGMRES<VectorType>::connect_re_orthogonalization_slot(
  const std::function<void(int)> &slot)
{
  return re_orthogonalize_signal.connect(slot);
}



template <class VectorType>
double
SolverGMRES<VectorType>::criterion()
{
  // dummy implementation. this function is not needed for the present
  // implementation of gmres
  Assert(false, ExcInternalError());
  return 0;
}


//----------------------------------------------------------------------//

template <class VectorType>
SolverFGMRES<VectorType>::SolverFGMRES(SolverControl &           cn,
                                       VectorMemory<VectorType> &mem,
                                       const AdditionalData &    data)
  : SolverBase<VectorType>(cn, mem)
  , additional_data(data)
{}



template <class VectorType>
SolverFGMRES<VectorType>::SolverFGMRES(SolverControl &       cn,
                                       const AdditionalData &data)
  : SolverBase<VectorType>(cn)
  , additional_data(data)
{}



template <class VectorType>
template <typename MatrixType, typename PreconditionerType>
void
SolverFGMRES<VectorType>::solve(const MatrixType &        A,
                                VectorType &              x,
                                const VectorType &        b,
                                const PreconditionerType &preconditioner)
{
  LogStream::Prefix prefix("FGMRES");

  SolverControl::State iteration_state = SolverControl::iterate;

  const unsigned int basis_size = additional_data.max_basis_size;

  // Generate an object where basis vectors are stored.
  typename internal::SolverGMRESImplementation::TmpVectors<VectorType> v(
    basis_size, this->memory);
  typename internal::SolverGMRESImplementation::TmpVectors<VectorType> z(
    basis_size, this->memory);

  // number of the present iteration; this number is not reset to zero upon a
  // restart
  unsigned int accumulated_iterations = 0;

  // matrix used for the orthogonalization process later
  H.reinit(basis_size + 1, basis_size);

  // Vectors for projected system
  Vector<double> projected_rhs;
  Vector<double> y;

  // Iteration starts here
  double res = std::numeric_limits<double>::lowest();

  typename VectorMemory<VectorType>::Pointer aux(this->memory);
  aux->reinit(x);
  do
    {
      A.vmult(*aux, x);
      aux->sadd(-1., 1., b);

      double beta     = aux->l2_norm();
      res             = beta;
      iteration_state = this->iteration_status(accumulated_iterations, res, x);
      if (iteration_state == SolverControl::success)
        break;

      H.reinit(basis_size + 1, basis_size);
      double a = beta;

      for (unsigned int j = 0; j < basis_size; ++j)
        {
          if (a != 0) // treat lucky breakdown
            v(j, x).equ(1. / a, *aux);
          else
            v(j, x) = 0.;


          preconditioner.vmult(z(j, x), v[j]);
          A.vmult(*aux, z[j]);

          // Gram-Schmidt
          H(0, j) = *aux * v[0];
          for (unsigned int i = 1; i <= j; ++i)
            H(i, j) = aux->add_and_dot(-H(i - 1, j), v[i - 1], v[i]);
          H(j + 1, j) = a = std::sqrt(aux->add_and_dot(-H(j, j), v[j], *aux));

          // Compute projected solution

          if (j > 0)
            {
              H1.reinit(j + 1, j);
              projected_rhs.reinit(j + 1);
              y.reinit(j);
              projected_rhs(0) = beta;
              H1.fill(H);

              // check convergence. note that the vector 'x' we pass to the
              // criterion is not the final solution we compute if we
              // decide to jump out of the iteration (we update 'x' again
              // right after the current loop)
              Householder<double> house(H1);
              res = house.least_squares(y, projected_rhs);
              iteration_state =
                this->iteration_status(++accumulated_iterations, res, x);
              if (iteration_state != SolverControl::iterate)
                break;
            }
        }

      // Update solution vector
      for (unsigned int j = 0; j < y.size(); ++j)
        x.add(y(j), z[j]);
    }
  while (iteration_state == SolverControl::iterate);

  // in case of failure: throw exception
  if (iteration_state != SolverControl::success)
    AssertThrow(false,
                SolverControl::NoConvergence(accumulated_iterations, res));
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
