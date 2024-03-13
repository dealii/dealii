// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1999 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_solver_gmres_h
#define dealii_solver_gmres_h



#include <deal.II/base/config.h>

#include <deal.II/base/logstream.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/lac/block_vector_base.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/orthogonalization.h>
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
 * The AdditionalData structure allows to control the size of the Arnoldi
 * basis used for orthogonalization (default: 30 vectors). It is related to
 * the number of temporary vectors used, which is the basis size plus
 * two. Additionally, it allows you to choose between right or left
 * preconditioning (default: left preconditioning). Furthermore, it
 * includes a flag indicating whether or not the default residual is used as
 * stopping criterion and an option for the orthogonalization algorithm, see
 * LinearAlgebra::OrthogonalizationStrategy for available options.
 *
 *
 * <h3>Left versus right preconditioning</h3>
 *
 * @p AdditionalData allows you to choose between left and right
 * preconditioning. As expected, this switches between solving for the systems
 * <i>P<sup>-1</sup>A</i> and <i>AP<sup>-1</sup></i>, respectively.
 *
 * A second consequence is the type of residual used to measure
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
 * The maximal basis size is controlled by AdditionalData::max_basis_size. If
 * the number of iteration steps exceeds this number, all basis vectors are
 * discarded and the iteration starts anew from the approximation obtained so
 * far. This algorithm strategy is typically called restarted GMRES method.
 *
 * Note that the minimizing property of GMRES only pertains to the Krylov
 * space spanned by the Arnoldi basis. Therefore, restarted GMRES is
 * <b>not</b> minimizing anymore. The choice of the basis length is a
 * trade-off between memory consumption and convergence speed, since a longer
 * basis means minimization over a larger space.
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
template <typename VectorType = Vector<double>>
class SolverGMRES : public SolverBase<VectorType>
{
public:
  /**
   * Standardized data struct to pipe additional data to the solver.
   */
  struct AdditionalData
  {
    /**
     * Constructor. By default, set the size of the Arnoldi basis to 30. Also
     * set preconditioning from left, the residual of the stopping criterion
     * to the default residual, and re-orthogonalization only if
     * necessary. Also, the batched mode with reduced functionality to track
     * information is disabled by default. Finally, the default
     * orthogonalization algorithm is the modified Gram-Schmidt method.
     */
    explicit AdditionalData(
      const unsigned int max_basis_size             = 30,
      const bool         right_preconditioning      = false,
      const bool         use_default_residual       = true,
      const bool         force_re_orthogonalization = false,
      const bool         batched_mode               = false,
      const LinearAlgebra::OrthogonalizationStrategy
        orthogonalization_strategy =
          LinearAlgebra::OrthogonalizationStrategy::modified_gram_schmidt);

    /**
     * Maximum number of temporary vectors. Together with max_basis_size, this
     * parameter controls the size of the Arnoldi basis, which corresponds to
     * max_n_tmp_vectors-2 as used in previous versions of the deal.II
     * library. SolverGMRES assumes that there are at least three temporary
     * vectors, so this value must be greater than or equal to three. If both
     * this variable and max_basis_size are set to a non-zero value, the
     * choice in max_basis_size takes precedence.
     *
     * @deprecated Use max_basis_size instead.
     */
    unsigned int max_n_tmp_vectors;

    /**
     * Maximum size of the Arnoldi basis. SolverGMRES assumes that there is at
     * least one vector in the Arnoldi basis, so this value must be greater
     * than or equal to one. Note that whenever this variable is set to a
     * non-zero value, including the value set by the default constructor,
     * this variable takes precedence over max_n_tmp_vectors.
     */
    unsigned int max_basis_size;

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
    LinearAlgebra::OrthogonalizationStrategy orthogonalization_strategy;
  };

  /**
   * Constructor.
   */
  SolverGMRES(SolverControl            &cn,
              VectorMemory<VectorType> &mem,
              const AdditionalData     &data = AdditionalData());

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
  solve(const MatrixType         &A,
        VectorType               &x,
        const VectorType         &b,
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
                                                &hessenberg_signal,
    const boost::signals2::signal<void(double)> &cond_signal);

  /**
   * Projected system matrix
   */
  FullMatrix<double> H;

  /**
   * Auxiliary vector for orthogonalization
   */
  Vector<double> projected_rhs;

  /**
   * Auxiliary vector for orthogonalization
   */
  std::vector<std::pair<double, double>> givens_rotations;

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
template <typename VectorType = Vector<double>>
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
    explicit AdditionalData(
      const unsigned int max_basis_size = 30,
      const LinearAlgebra::OrthogonalizationStrategy
        orthogonalization_strategy =
          LinearAlgebra::OrthogonalizationStrategy::modified_gram_schmidt)
      : max_basis_size(max_basis_size)
      , orthogonalization_strategy(orthogonalization_strategy)
    {}

    /**
     * Maximum basis size.
     */
    unsigned int max_basis_size;

    /**
     * Strategy to orthogonalize vectors.
     */
    LinearAlgebra::OrthogonalizationStrategy orthogonalization_strategy;
  };

  /**
   * Constructor.
   */
  SolverFGMRES(SolverControl            &cn,
               VectorMemory<VectorType> &mem,
               const AdditionalData     &data = AdditionalData());

  /**
   * Constructor. Use an object of type GrowingVectorMemory as a default to
   * allocate memory.
   */
  SolverFGMRES(SolverControl        &cn,
               const AdditionalData &data = AdditionalData());

  /**
   * Solve the linear system $Ax=b$ for x.
   */
  template <typename MatrixType, typename PreconditionerType>
  void
  solve(const MatrixType         &A,
        VectorType               &x,
        const VectorType         &b,
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

template <typename VectorType>
inline SolverGMRES<VectorType>::AdditionalData::AdditionalData(
  const unsigned int                             max_basis_size,
  const bool                                     right_preconditioning,
  const bool                                     use_default_residual,
  const bool                                     force_re_orthogonalization,
  const bool                                     batched_mode,
  const LinearAlgebra::OrthogonalizationStrategy orthogonalization_strategy)
  : max_n_tmp_vectors(0)
  , max_basis_size(max_basis_size)
  , right_preconditioning(right_preconditioning)
  , use_default_residual(use_default_residual)
  , force_re_orthogonalization(force_re_orthogonalization)
  , batched_mode(batched_mode)
  , orthogonalization_strategy(orthogonalization_strategy)
{
  Assert(max_basis_size >= 1,
         ExcMessage("SolverGMRES needs at least one vector in the "
                    "Arnoldi basis."));
}



template <typename VectorType>
SolverGMRES<VectorType>::SolverGMRES(SolverControl            &cn,
                                     VectorMemory<VectorType> &mem,
                                     const AdditionalData     &data)
  : SolverBase<VectorType>(cn, mem)
  , additional_data(data)
  , solver_control(cn)
{}



template <typename VectorType>
SolverGMRES<VectorType>::SolverGMRES(SolverControl        &cn,
                                     const AdditionalData &data)
  : SolverBase<VectorType>(cn)
  , additional_data(data)
  , solver_control(cn)
{}



namespace internal
{
  namespace SolverGMRESImplementation
  {
    template <typename VectorType>
    inline TmpVectors<VectorType>::TmpVectors(const unsigned int max_size,
                                              VectorMemory<VectorType> &vmem)
      : mem(vmem)
      , data(max_size)
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



    template <typename VectorType>
    unsigned int
    TmpVectors<VectorType>::size() const
    {
      return (data.size() > 0 ? data.size() - 1 : 0);
    }



    template <typename VectorType, typename Enable = void>
    struct is_dealii_compatible_distributed_vector;

    template <typename VectorType>
    struct is_dealii_compatible_distributed_vector<
      VectorType,
      std::enable_if_t<!internal::is_block_vector<VectorType>>>
    {
      static constexpr bool value = std::is_same_v<
        VectorType,
        LinearAlgebra::distributed::Vector<typename VectorType::value_type,
                                           MemorySpace::Host>>;
    };



    template <typename VectorType>
    struct is_dealii_compatible_distributed_vector<
      VectorType,
      std::enable_if_t<internal::is_block_vector<VectorType>>>
    {
      static constexpr bool value = std::is_same_v<
        typename VectorType::BlockType,
        LinearAlgebra::distributed::Vector<typename VectorType::value_type,
                                           MemorySpace::Host>>;
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



    template <bool delayed_reorthogonalization,
              typename VectorType,
              std::enable_if_t<
                !is_dealii_compatible_distributed_vector<VectorType>::value,
                VectorType> * = nullptr>
    void
    Tvmult_add(const unsigned int dim,
               const VectorType  &vv,
               const internal::SolverGMRESImplementation::TmpVectors<VectorType>
                              &orthogonal_vectors,
               Vector<double> &h)
    {
      for (unsigned int i = 0; i < dim; ++i)
        {
          h(i) += vv * orthogonal_vectors[i];
          if (delayed_reorthogonalization)
            h(dim + i) += orthogonal_vectors[i] * orthogonal_vectors[dim - 1];
        }
      if (delayed_reorthogonalization)
        h(dim + dim) += vv * vv;
    }



    template <bool delayed_reorthogonalization,
              typename VectorType,
              std::enable_if_t<
                is_dealii_compatible_distributed_vector<VectorType>::value,
                VectorType> * = nullptr>
    void
    Tvmult_add(const unsigned int dim,
               const VectorType  &vv,
               const internal::SolverGMRESImplementation::TmpVectors<VectorType>
                              &orthogonal_vectors,
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
              for (unsigned int i = 0; i < dim; ++i)
                hs[i] = 0.0;
              VectorizedArray<double>
                correct[delayed_reorthogonalization ? 129 : 1];
              if (delayed_reorthogonalization)
                for (unsigned int i = 0; i < dim + 1; ++i)
                  correct[i] = 0.0;

              unsigned int c = 0;

              constexpr unsigned int inner_batch_size =
                delayed_reorthogonalization ? 4 : 8;

              for (; c < block(vv, b).locally_owned_size() / n_lanes /
                           inner_batch_size;
                   ++c, j += n_lanes * inner_batch_size)
                {
                  VectorizedArray<double> vvec[inner_batch_size];
                  for (unsigned int k = 0; k < inner_batch_size; ++k)
                    vvec[k].load(block(vv, b).begin() + j + k * n_lanes);
                  VectorizedArray<double> last_vector[inner_batch_size];
                  for (unsigned int k = 0; k < inner_batch_size; ++k)
                    last_vector[k].load(
                      block(orthogonal_vectors[dim - 1], b).begin() + j +
                      k * n_lanes);

                  {
                    VectorizedArray<double> local_sum_0 =
                      last_vector[0] * vvec[0];
                    VectorizedArray<double> local_sum_1 =
                      last_vector[0] * last_vector[0];
                    VectorizedArray<double> local_sum_2 = vvec[0] * vvec[0];
                    for (unsigned int k = 1; k < inner_batch_size; ++k)
                      {
                        local_sum_0 += last_vector[k] * vvec[k];
                        if (delayed_reorthogonalization)
                          {
                            local_sum_1 += last_vector[k] * last_vector[k];
                            local_sum_2 += vvec[k] * vvec[k];
                          }
                      }
                    hs[dim - 1] += local_sum_0;
                    if (delayed_reorthogonalization)
                      {
                        correct[dim - 1] += local_sum_1;
                        correct[dim] += local_sum_2;
                      }
                  }

                  for (unsigned int i = 0; i < dim - 1; ++i)
                    {
                      // break the dependency chain into the field hs[i] for
                      // small sizes i by first accumulating 4 or 8 results
                      // into a local variable
                      VectorizedArray<double> temp;
                      temp.load(block(orthogonal_vectors[i], b).begin() + j);
                      VectorizedArray<double> local_sum_0 = temp * vvec[0];
                      VectorizedArray<double> local_sum_1 =
                        delayed_reorthogonalization ? temp * last_vector[0] :
                                                      0.;
                      for (unsigned int k = 1; k < inner_batch_size; ++k)
                        {
                          temp.load(block(orthogonal_vectors[i], b).begin() +
                                    j + k * n_lanes);
                          local_sum_0 += temp * vvec[k];
                          if (delayed_reorthogonalization)
                            local_sum_1 += temp * last_vector[k];
                        }
                      hs[i] += local_sum_0;
                      if (delayed_reorthogonalization)
                        correct[i] += local_sum_1;
                    }
                }

              c *= inner_batch_size;
              for (; c < block(vv, b).locally_owned_size() / n_lanes;
                   ++c, j += n_lanes)
                {
                  VectorizedArray<double> vvec, last_vector;
                  vvec.load(block(vv, b).begin() + j);
                  last_vector.load(
                    block(orthogonal_vectors[dim - 1], b).begin() + j);
                  hs[dim - 1] += last_vector * vvec;
                  if (delayed_reorthogonalization)
                    {
                      correct[dim - 1] += last_vector * last_vector;
                      correct[dim] += vvec * vvec;
                    }

                  for (unsigned int i = 0; i < dim - 1; ++i)
                    {
                      VectorizedArray<double> temp;
                      temp.load(block(orthogonal_vectors[i], b).begin() + j);
                      hs[i] += temp * vvec;
                      if (delayed_reorthogonalization)
                        correct[i] += temp * last_vector;
                    }
                }

              for (unsigned int i = 0; i < dim; ++i)
                {
                  h(i) += hs[i].sum();
                  if (delayed_reorthogonalization)
                    h(i + dim) += correct[i].sum();
                }
              if (delayed_reorthogonalization)
                h(dim + dim) += correct[dim].sum();
            }

          // remainder loop of optimized path or non-optimized path (if
          // dim>128)
          for (; j < block(vv, b).locally_owned_size(); ++j)
            {
              const double vvec = block(vv, b).local_element(j);
              const double last_vector =
                block(orthogonal_vectors[dim - 1], b).local_element(j);
              h(dim - 1) += last_vector * vvec;
              if (delayed_reorthogonalization)
                {
                  h(dim + dim - 1) += last_vector * last_vector;
                  h(dim + dim) += vvec * vvec;
                }
              for (unsigned int i = 0; i < dim - 1; ++i)
                {
                  const double temp =
                    block(orthogonal_vectors[i], b).local_element(j);
                  h(i) += temp * vvec;
                  if (delayed_reorthogonalization)
                    h(dim + i) += temp * last_vector;
                }
            }
        }

      Utilities::MPI::sum(h, block(vv, 0).get_mpi_communicator(), h);
    }



    template <bool delayed_reorthogonalization,
              typename VectorType,
              std::enable_if_t<
                !is_dealii_compatible_distributed_vector<VectorType>::value,
                VectorType> * = nullptr>
    double
    subtract_and_norm(
      const unsigned int dim,
      const internal::SolverGMRESImplementation::TmpVectors<VectorType>
                           &orthogonal_vectors,
      const Vector<double> &h,
      VectorType           &vv)
    {
      Assert(dim > 0, ExcInternalError());

      VectorType &last_vector =
        const_cast<VectorType &>(orthogonal_vectors[dim - 1]);
      for (unsigned int i = 0; i < dim - 1; ++i)
        {
          if (delayed_reorthogonalization && i + 2 < dim)
            last_vector.add(-h(dim + i), orthogonal_vectors[i]);
          vv.add(-h(i), orthogonal_vectors[i]);
        }

      if (delayed_reorthogonalization)
        {
          if (dim > 1)
            last_vector.sadd(1. / h(dim + dim - 1),
                             -h(dim + dim - 2) / h(dim + dim - 1),
                             orthogonal_vectors[dim - 2]);

          // h(dim + dim) is lucky breakdown
          const double scaling_factor_vv =
            h(dim + dim) > 0.0 ? 1. / (h(dim + dim - 1) * h(dim + dim)) :
                                 1. / (h(dim + dim - 1) * h(dim + dim - 1));
          vv.sadd(scaling_factor_vv,
                  -h(dim - 1) * scaling_factor_vv,
                  last_vector);
          return vv.l2_norm();
        }
      else
        return std::sqrt(
          vv.add_and_dot(-h(dim - 1), orthogonal_vectors[dim - 1], vv));
    }



    template <bool delayed_reorthogonalization,
              typename VectorType,
              std::enable_if_t<
                is_dealii_compatible_distributed_vector<VectorType>::value,
                VectorType> * = nullptr>
    double
    subtract_and_norm(
      const unsigned int dim,
      const internal::SolverGMRESImplementation::TmpVectors<VectorType>
                           &orthogonal_vectors,
      const Vector<double> &h,
      VectorType           &vv)
    {
      static constexpr unsigned int n_lanes = VectorizedArray<double>::size();

      double      norm_vv_temp = 0.0;
      VectorType &last_vector =
        const_cast<VectorType &>(orthogonal_vectors[dim - 1]);
      const double inverse_norm_previous =
        delayed_reorthogonalization ? 1. / h(dim + dim - 1) : 0.;
      const double scaling_factor_vv =
        delayed_reorthogonalization ?
          (h(dim + dim) > 0.0 ? inverse_norm_previous / h(dim + dim) :
                                inverse_norm_previous / h(dim + dim - 1)) :
          0.;

      for (unsigned int b = 0; b < n_blocks(vv); ++b)
        {
          VectorizedArray<double> norm_vv_temp_vectorized = 0.0;

          constexpr unsigned int inner_batch_size =
            delayed_reorthogonalization ? 4 : 8;

          unsigned int j = 0;
          unsigned int c = 0;
          for (; c <
                 block(vv, b).locally_owned_size() / n_lanes / inner_batch_size;
               ++c, j += n_lanes * inner_batch_size)
            {
              VectorizedArray<double> temp[inner_batch_size];
              VectorizedArray<double> last_vec[inner_batch_size];

              const double last_factor = h(dim - 1);
              for (unsigned int k = 0; k < inner_batch_size; ++k)
                {
                  temp[k].load(block(vv, b).begin() + j + k * n_lanes);
                  last_vec[k].load(block(last_vector, b).begin() + j +
                                   k * n_lanes);
                  if (!delayed_reorthogonalization)
                    temp[k] -= last_factor * last_vec[k];
                }

              for (unsigned int i = 0; i < dim - 1; ++i)
                {
                  const double factor = h(i);
                  const double correction_factor =
                    (delayed_reorthogonalization ? h(dim + i) : 0.0);
                  for (unsigned int k = 0; k < inner_batch_size; ++k)
                    {
                      VectorizedArray<double> vec;
                      vec.load(block(orthogonal_vectors[i], b).begin() + j +
                               k * n_lanes);
                      temp[k] -= factor * vec;
                      if (delayed_reorthogonalization)
                        last_vec[k] -= correction_factor * vec;
                    }
                }

              if (delayed_reorthogonalization)
                for (unsigned int k = 0; k < inner_batch_size; ++k)
                  {
                    last_vec[k] = last_vec[k] * inverse_norm_previous;
                    last_vec[k].store(block(last_vector, b).begin() + j +
                                      k * n_lanes);
                    temp[k] -= last_factor * last_vec[k];
                    temp[k] = temp[k] * scaling_factor_vv;
                    temp[k].store(block(vv, b).begin() + j + k * n_lanes);
                  }
              else
                for (unsigned int k = 0; k < inner_batch_size; ++k)
                  {
                    temp[k].store(block(vv, b).begin() + j + k * n_lanes);
                    norm_vv_temp_vectorized += temp[k] * temp[k];
                  }
            }

          c *= inner_batch_size;
          for (; c < block(vv, b).locally_owned_size() / n_lanes;
               ++c, j += n_lanes)
            {
              VectorizedArray<double> temp, last_vec;
              temp.load(block(vv, b).begin() + j);
              last_vec.load(block(last_vector, b).begin() + j);
              if (!delayed_reorthogonalization)
                temp -= h(dim - 1) * last_vec;

              for (unsigned int i = 0; i < dim - 1; ++i)
                {
                  VectorizedArray<double> vec;
                  vec.load(block(orthogonal_vectors[i], b).begin() + j);
                  temp -= h(i) * vec;
                  if (delayed_reorthogonalization)
                    last_vec -= h(dim + i) * vec;
                }

              if (delayed_reorthogonalization)
                {
                  last_vec = last_vec * inverse_norm_previous;
                  last_vec.store(block(last_vector, b).begin() + j);
                  temp -= h(dim - 1) * last_vec;
                  temp = temp * scaling_factor_vv;
                  temp.store(block(vv, b).begin() + j);
                }
              else
                {
                  temp.store(block(vv, b).begin() + j);
                  norm_vv_temp_vectorized += temp * temp;
                }
            }

          if (!delayed_reorthogonalization)
            norm_vv_temp += norm_vv_temp_vectorized.sum();

          for (; j < block(vv, b).locally_owned_size(); ++j)
            {
              double temp     = block(vv, b).local_element(j);
              double last_vec = block(last_vector, b).local_element(j);
              if (delayed_reorthogonalization)
                {
                  for (unsigned int i = 0; i < dim - 1; ++i)
                    {
                      const double vec =
                        block(orthogonal_vectors[i], b).local_element(j);
                      temp -= h(i) * vec;
                      last_vec -= h(dim + i) * vec;
                    }
                  last_vec *= inverse_norm_previous;
                  block(last_vector, b).local_element(j) = last_vec;
                  temp -= h(dim - 1) * last_vec;
                  temp *= scaling_factor_vv;
                }
              else
                {
                  temp -= h(dim - 1) * last_vec;
                  for (unsigned int i = 0; i < dim - 1; ++i)
                    temp -=
                      h(i) * block(orthogonal_vectors[i], b).local_element(j);
                  norm_vv_temp += temp * temp;
                }
              block(vv, b).local_element(j) = temp;
            }
        }

      return std::sqrt(
        Utilities::MPI::sum(norm_vv_temp, block(vv, 0).get_mpi_communicator()));
    }


    template <typename VectorType,
              std::enable_if_t<
                !is_dealii_compatible_distributed_vector<VectorType>::value,
                VectorType> * = nullptr>
    double
    sadd_and_norm(VectorType       &v,
                  const double      factor_a,
                  const VectorType &b,
                  const double      factor_b)
    {
      v.sadd(factor_a, factor_b, b);
      return v.l2_norm();
    }


    template <typename VectorType,
              std::enable_if_t<
                is_dealii_compatible_distributed_vector<VectorType>::value,
                VectorType> * = nullptr>
    double
    sadd_and_norm(VectorType       &v,
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

      return std::sqrt(
        Utilities::MPI::sum(norm, block(v, 0).get_mpi_communicator()));
    }



    template <typename VectorType,
              std::enable_if_t<
                !is_dealii_compatible_distributed_vector<VectorType>::value,
                VectorType> * = nullptr>
    void
    add(VectorType           &p,
        const unsigned int    dim,
        const Vector<double> &h,
        const internal::SolverGMRESImplementation::TmpVectors<VectorType>
                  &tmp_vectors,
        const bool zero_out)
    {
      if (zero_out)
        p.equ(h(0), tmp_vectors[0]);
      else
        p.add(h(0), tmp_vectors[0]);

      for (unsigned int i = 1; i < dim; ++i)
        p.add(h(i), tmp_vectors[i]);
    }



    template <typename VectorType,
              std::enable_if_t<
                is_dealii_compatible_distributed_vector<VectorType>::value,
                VectorType> * = nullptr>
    void
    add(VectorType           &p,
        const unsigned int    dim,
        const Vector<double> &h,
        const internal::SolverGMRESImplementation::TmpVectors<VectorType>
                  &tmp_vectors,
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
    template <typename VectorType>
    inline void
    iterated_gram_schmidt(
      const LinearAlgebra::OrthogonalizationStrategy orthogonalization_strategy,
      const TmpVectors<VectorType>                  &orthogonal_vectors,
      const unsigned int                             dim,
      const unsigned int                             accumulated_iterations,
      VectorType                                    &vv,
      Vector<double>                                &h,
      FullMatrix<double>                            &H,
      FullMatrix<double>                            &H_orig,
      bool                                          &reorthogonalize,
      const boost::signals2::signal<void(int)>      &reorthogonalize_signal =
        boost::signals2::signal<void(int)>())
    {
      Assert(dim > 0, ExcInternalError());
      if (orthogonalization_strategy ==
          LinearAlgebra::OrthogonalizationStrategy::
            delayed_classical_gram_schmidt)
        {
          const double scaling_norm_previous = dim > 0 ? h(dim + dim - 2) : 1.;

          for (unsigned int i = 0; i < dim + dim + 1; ++i)
            h(i) = 0;

          // This is algorithm 4 of Bielich et al. (2022)
          Tvmult_add<true>(dim, vv, orthogonal_vectors, h);

          // delayed correction terms
          double tmp = 0;
          for (unsigned int i = 0; i < dim - 1; ++i)
            tmp += h(dim + i) * h(dim + i);
          const double alpha_j = h(dim + dim - 1) > tmp ?
                                   std::sqrt(h(dim + dim - 1) - tmp) :
                                   std::sqrt(h(dim + dim - 1));
          h(dim + dim - 1)     = alpha_j;

          tmp = 0;
          for (unsigned int i = 0; i < dim - 1; ++i)
            tmp += h(i) * h(dim + i);
          h(dim - 1) = (h(dim - 1) - tmp) / alpha_j;

          // representation of H(j-1)
          if (dim > 1)
            {
              for (unsigned int i = 0; i < dim - 1; ++i)
                H(i, dim - 2) += h(dim + i) * scaling_norm_previous;
              H(dim - 1, dim - 2) = alpha_j * scaling_norm_previous;

              // correct H_orig according to H
              for (unsigned int i = 0; i < dim; ++i)
                H_orig(i, dim - 2) = H(i, dim - 2);
            }
          for (unsigned int i = 0; i < dim; ++i)
            {
              double sum = 0;
              for (unsigned int j = (i == 0 ? 0 : i - 1); j < dim - 1; ++j)
                sum += H_orig(i, j) * h(dim + j);
              H(i, dim - 1) = (h(i) - sum) / alpha_j;
            }

          // Compute estimate norm for approximate convergence criterion (to
          // be corrected in next iteration)
          double sum = 0;
          for (unsigned int i = 0; i < dim - 1; ++i)
            sum += h(i) * h(i);
          sum += (2. - 1.) * h(dim - 1) * h(dim - 1);
          H(dim, dim - 1) = std::sqrt(std::abs(h(dim + dim) - sum)) / alpha_j;

          // projection and delayed reorthogonalization. We scale the vector
          // vv here by the preliminary norm to avoid working with too large
          // values and correct to the actual norm in high precision in the
          // next iteration.
          h(dim + dim) = H(dim, dim - 1);
          subtract_and_norm<true>(dim, orthogonal_vectors, h, vv);
        }
      else
        {
          const unsigned int inner_iteration = dim - 1;

          // need initial norm for detection of re-orthogonalization, see below
          double     norm_vv       = 0.0;
          double     norm_vv_start = 0;
          const bool consider_reorthogonalize =
            (reorthogonalize == false) && (inner_iteration % 5 == 4);
          if (consider_reorthogonalize)
            norm_vv_start = vv.l2_norm();

          for (unsigned int i = 0; i < dim; ++i)
            h(i) = 0;

          // run two loops with index 0: orthogonalize, 1: reorthogonalize
          for (unsigned int c = 0; c < 2; ++c)
            {
              // Orthogonalization
              if (orthogonalization_strategy ==
                  LinearAlgebra::OrthogonalizationStrategy::
                    modified_gram_schmidt)
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
                       LinearAlgebra::OrthogonalizationStrategy::
                         classical_gram_schmidt)
                {
                  Tvmult_add<false>(dim, vv, orthogonal_vectors, h);
                  norm_vv =
                    subtract_and_norm<false>(dim, orthogonal_vectors, h, vv);
                }
              else
                {
                  AssertThrow(false, ExcNotImplemented());
                }

              if (c == 1)
                break; // reorthogonalization already performed -> finished

              // Re-orthogonalization if loss of orthogonality detected. For the
              // test, use a strategy discussed in C. T. Kelley, Iterative
              // Methods for Linear and Nonlinear Equations, SIAM, Philadelphia,
              // 1995: Compare the norm of vv after orthogonalization with its
              // norm when starting the orthogonalization. If vv became very
              // small (here: less than the square root of the machine precision
              // times 10), it is almost in the span of the previous vectors,
              // which indicates loss of precision.
              if (consider_reorthogonalize)
                {
                  if (norm_vv >
                      10. * norm_vv_start *
                        std::sqrt(std::numeric_limits<
                                  typename VectorType::value_type>::epsilon()))
                    break;

                  else
                    {
                      reorthogonalize = true;
                      if (!reorthogonalize_signal.empty())
                        reorthogonalize_signal(accumulated_iterations);
                    }
                }

              if (reorthogonalize == false)
                break; // no reorthogonalization needed -> finished
            }

          for (unsigned int i = 0; i < dim; ++i)
            H(i, dim - 1) = h(i);
          H(dim, dim - 1) = norm_vv;

          // norm_vv is a lucky breakdown, the solver will reach convergence,
          // but we must not divide by zero here.
          if (norm_vv != 0)
            vv *= 1. / H(dim, inner_iteration);
        }
    }



    // A comparator for better printing eigenvalues
    inline bool
    complex_less_pred(const std::complex<double> &x,
                      const std::complex<double> &y)
    {
      return x.real() < y.real() ||
             (x.real() == y.real() && x.imag() < y.imag());
    }



    // A function to compute the Givens rotation for the QR factorization of
    // the Hessenberg matrix involved in the Arnoldi process, transforming it
    // into an upper triangular matrix.
    inline void
    givens_rotation(FullMatrix<double>                     &H,
                    Vector<double>                         &b,
                    std::vector<std::pair<double, double>> &rotations,
                    const int                               col)
    {
      for (int i = 0; i < col; ++i)
        {
          const double c   = rotations[i].first;
          const double s   = rotations[i].second;
          const double tmp = H(i, col);
          H(i, col)        = c * tmp + s * H(i + 1, col);
          H(i + 1, col)    = -s * tmp + c * H(i + 1, col);
        }

      const double H_col1   = H(col + 1, col);
      double      &H_col    = H(col, col);
      const double r        = 1. / std::sqrt(H_col * H_col + H_col1 * H_col1);
      rotations[col].second = H_col1 * r;
      rotations[col].first  = H_col * r;
      H_col = rotations[col].first * H_col + rotations[col].second * H_col1;
      b(col + 1) = -rotations[col].second * b(col);
      b(col) *= rotations[col].first;
    }



    // Function that determines factor for givens rotation in the right hand
    // side, without actually performing the elimination in the matrix. This
    // function is necessary to get a residual estimate for the classical
    // Gram-Schmidt algorithm with delayed reorthogonalization, which
    // maintains an accurate Hessenberg matrix that lags behind by one
    // iteration compared to the residual we want to estimate. For how the
    // code is derive, compare with the other function above and how itwould
    // compute b(col + 1), removing all unnecessary computations.
    inline double
    compute_givens_rotation_rhs(
      const FullMatrix<double>                     &H,
      const Vector<double>                         &b,
      const std::vector<std::pair<double, double>> &rotations,
      const int                                     col)
    {
      double H_col = H(0, col);
      for (int i = 0; i < col; ++i)
        {
          const double c = rotations[i].first;
          const double s = rotations[i].second;
          H_col          = -s * H_col + c * H(i + 1, col);
        }

      const double H_col1 = H(col + 1, col);
      const double r      = 1. / std::sqrt(H_col * H_col + H_col1 * H_col1);
      return -H_col1 * r * b(col);
    }



    // A function to solve the (upper) triangular system after Givens
    // rotations on a matrix that has possibly unused rows and columns
    inline void
    solve_triangular(const unsigned int        dim,
                     const FullMatrix<double> &H,
                     const Vector<double>     &rhs,
                     Vector<double>           &solution)
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



template <typename VectorType>
inline void
SolverGMRES<VectorType>::compute_eigs_and_cond(
  const FullMatrix<double> &H_orig,
  const unsigned int        dim,
  const boost::signals2::signal<void(const std::vector<std::complex<double>> &)>
    &eigenvalues_signal,
  const boost::signals2::signal<void(const FullMatrix<double> &)>
                                              &hessenberg_signal,
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



template <typename VectorType>
template <typename MatrixType, typename PreconditionerType>
void
SolverGMRES<VectorType>::solve(const MatrixType         &A,
                               VectorType               &x,
                               const VectorType         &b,
                               const PreconditionerType &preconditioner)
{
  std::unique_ptr<LogStream::Prefix> prefix;
  if (!additional_data.batched_mode)
    prefix = std::make_unique<LogStream::Prefix>("GMRES");

  // extra call to std::max to placate static analyzers: coverity rightfully
  // complains that data.max_n_tmp_vectors - 2 may overflow
  const unsigned int basis_size =
    (additional_data.max_basis_size > 0 ?
       additional_data.max_basis_size :
       std::max(additional_data.max_n_tmp_vectors, 3u) - 2);

  // Generate an object where basis vectors are stored.
  internal::SolverGMRESImplementation::TmpVectors<VectorType> basis_vectors(
    basis_size + 2, this->memory);
  const bool delayed_reorthogonalization =
    additional_data.orthogonalization_strategy ==
    LinearAlgebra::OrthogonalizationStrategy::delayed_classical_gram_schmidt;

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
  if (do_eigenvalues || delayed_reorthogonalization)
    H_orig.reinit(basis_size + 1, basis_size);

  // matrix used for the orthogonalization process later
  H.reinit(basis_size + 1, basis_size, /* omit_initialization */ true);

  // some additional vectors, also used in the orthogonalization
  projected_rhs.reinit(basis_size + 1);
  givens_rotations.resize(basis_size);
  if (delayed_reorthogonalization)
    h.reinit(2 * basis_size + 3);
  else
    h.reinit(basis_size + 1);

  SolverControl::State iteration_state = SolverControl::iterate;
  double               res             = std::numeric_limits<double>::lowest();

  // switch to determine whether we want a left or a right preconditioner. at
  // present, left is default, but both ways are implemented
  const bool left_precondition = !additional_data.right_preconditioning;

  // Per default the left preconditioned GMRES uses the preconditioned
  // residual and the right preconditioned GMRES uses the unpreconditioned
  // residual as stopping criterion.
  const bool use_default_residual = additional_data.use_default_residual;

  // define an alias
  VectorType &p = basis_vectors(basis_size + 1, x);

  // Following vectors are needed when we are not using the default residuals
  // as stopping criterion
  typename VectorMemory<VectorType>::Pointer r;
  typename VectorMemory<VectorType>::Pointer x_;
  std::unique_ptr<dealii::Vector<double>>    gamma;
  if (!use_default_residual)
    {
      r  = std::move(typename VectorMemory<VectorType>::Pointer(this->memory));
      x_ = std::move(typename VectorMemory<VectorType>::Pointer(this->memory));
      r->reinit(x);
      x_->reinit(x);

      gamma = std::make_unique<dealii::Vector<double>>(projected_rhs.size());
    }

  bool re_orthogonalize = additional_data.force_re_orthogonalization;

  ///////////////////////////////////////////////////////////////////////////
  // outer iteration: loop until we either reach convergence or the maximum
  // number of iterations is exceeded. each cycle of this loop amounts to one
  // restart
  do
    {
      VectorType &v      = basis_vectors(0, x);
      double      norm_v = 0.;

      if (left_precondition)
        {
          A.vmult(p, x);
          p.sadd(-1., 1., b);
          preconditioner.vmult(v, p);
          norm_v = v.l2_norm();
        }
      else
        {
          A.vmult(v, x);
          norm_v = dealii::internal::SolverGMRESImplementation::sadd_and_norm(
            v, -1, b, 1.0);
        }

      projected_rhs(0) = norm_v;
      if (norm_v != 0)
        v /= norm_v;

      // check the residual here as well since it may be that we got the exact
      // (or an almost exact) solution vector at the outset. if we wouldn't
      // check here, the next scaling operation would produce garbage
      if (use_default_residual)
        {
          res = norm_v;
          if (additional_data.batched_mode)
            iteration_state = solver_control.check(accumulated_iterations, res);
          else
            iteration_state =
              this->iteration_status(accumulated_iterations, res, x);

          if (iteration_state != SolverControl::iterate)
            break;
        }
      else
        {
          deallog << "default_res=" << norm_v << std::endl;

          if (left_precondition)
            {
              A.vmult(*r, x);
              r->sadd(-1., 1., b);
            }
          else
            preconditioner.vmult(*r, v);

          res = r->l2_norm();
          if (additional_data.batched_mode)
            iteration_state = solver_control.check(accumulated_iterations, res);
          else
            iteration_state =
              this->iteration_status(accumulated_iterations, res, x);

          if (iteration_state != SolverControl::iterate)
            break;
        }

      // inner iteration doing at most as many steps as the size of the
      // Arnoldi basis
      unsigned int inner_iteration = 0;
      for (; (inner_iteration < basis_size &&
              iteration_state == SolverControl::iterate);
           ++inner_iteration)
        {
          ++accumulated_iterations;
          // yet another alias
          VectorType &vv = basis_vectors(inner_iteration + 1, x);

          if (left_precondition)
            {
              A.vmult(p, basis_vectors[inner_iteration]);
              preconditioner.vmult(vv, p);
            }
          else
            {
              preconditioner.vmult(p, basis_vectors[inner_iteration]);
              A.vmult(vv, p);
            }

          internal::SolverGMRESImplementation::iterated_gram_schmidt(
            additional_data.orthogonalization_strategy,
            basis_vectors,
            inner_iteration + 1,
            accumulated_iterations,
            vv,
            h,
            H,
            H_orig,
            re_orthogonalize,
            re_orthogonalize_signal);

          // for eigenvalues, get the resulting coefficients from the
          // orthogonalization process
          if (do_eigenvalues)
            for (unsigned int i = 0; i < inner_iteration + 2; ++i)
              H_orig(i, inner_iteration) = H(i, inner_iteration);

          //  Transformation into upper triangular structure
          if (delayed_reorthogonalization)
            {
              if (inner_iteration > 0)
                internal::SolverGMRESImplementation::givens_rotation(
                  H, projected_rhs, givens_rotations, inner_iteration - 1);
              res = std::fabs(internal::SolverGMRESImplementation::
                                compute_givens_rotation_rhs(H,
                                                            projected_rhs,
                                                            givens_rotations,
                                                            inner_iteration));
            }
          else
            {
              internal::SolverGMRESImplementation::givens_rotation(
                H, projected_rhs, givens_rotations, inner_iteration);

              //  default residual
              res = std::fabs(projected_rhs(inner_iteration + 1));
            }

          if (use_default_residual)
            {
              if (additional_data.batched_mode)
                iteration_state =
                  solver_control.check(accumulated_iterations, res);
              else
                iteration_state =
                  this->iteration_status(accumulated_iterations, res, x);
            }
          else
            {
              if (!additional_data.batched_mode)
                deallog << "default_res=" << res << std::endl;

              *x_    = x;
              *gamma = projected_rhs;
              internal::SolverGMRESImplementation::solve_triangular(
                inner_iteration + 1, H, *gamma, h);

              if (left_precondition)
                for (unsigned int i = 0; i < inner_iteration + 1; ++i)
                  x_->add(h(i), basis_vectors[i]);
              else
                {
                  p = 0.;
                  for (unsigned int i = 0; i < inner_iteration + 1; ++i)
                    p.add(h(i), basis_vectors[i]);
                  preconditioner.vmult(*r, p);
                  x_->add(1., *r);
                };
              A.vmult(*r, *x_);
              r->sadd(-1., 1., b);
              // Now *r contains the unpreconditioned residual!!
              if (left_precondition)
                {
                  res = r->l2_norm();
                  iteration_state =
                    this->iteration_status(accumulated_iterations, res, x);
                }
              else
                {
                  preconditioner.vmult(*x_, *r);
                  res = x_->l2_norm();

                  if (additional_data.batched_mode)
                    iteration_state =
                      solver_control.check(accumulated_iterations, res);
                  else
                    iteration_state =
                      this->iteration_status(accumulated_iterations, res, x);
                }
            }
        }

      // end of inner iteration. now calculate the solution from the temporary
      // vectors. do the last orthogonalization step (delayed by the algorithm
      // design) without reorthogonalization when solving the triangular
      // system
      if (delayed_reorthogonalization)
        {
          internal::SolverGMRESImplementation::givens_rotation(
            H, projected_rhs, givens_rotations, inner_iteration - 1);
        }
      internal::SolverGMRESImplementation::solve_triangular(inner_iteration,
                                                            H,
                                                            projected_rhs,
                                                            h);

      if (do_eigenvalues)
        compute_eigs_and_cond(H_orig,
                              inner_iteration,
                              all_eigenvalues_signal,
                              all_hessenberg_signal,
                              condition_number_signal);

      if (left_precondition)
        dealii::internal::SolverGMRESImplementation::add(
          x, inner_iteration, h, basis_vectors, false);
      else
        {
          dealii::internal::SolverGMRESImplementation::add(
            p, inner_iteration, h, basis_vectors, true);
          preconditioner.vmult(v, p);
          x.add(1., v);
        }

      // in the last round, print the eigenvalues from the last Arnoldi step
      if (iteration_state != SolverControl::iterate)
        {
          if (do_eigenvalues)
            compute_eigs_and_cond(H_orig,
                                  inner_iteration,
                                  eigenvalues_signal,
                                  hessenberg_signal,
                                  condition_number_signal);

          if (!additional_data.batched_mode && !krylov_space_signal.empty())
            {
              // Must normalize the last vector
              if (delayed_reorthogonalization &&
                  H(inner_iteration, inner_iteration - 1) != 0.0)
                basis_vectors[inner_iteration] /=
                  H(inner_iteration, inner_iteration - 1);

              krylov_space_signal(basis_vectors);
            }

          // end of outer iteration. restart if no convergence and the number of
          // iterations is not exceeded
        }
    }
  while (iteration_state == SolverControl::iterate);

  // in case of failure: throw exception
  AssertThrow(iteration_state == SolverControl::success,
              SolverControl::NoConvergence(accumulated_iterations, res));
}



template <typename VectorType>
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



template <typename VectorType>
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



template <typename VectorType>
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



template <typename VectorType>
boost::signals2::connection
SolverGMRES<VectorType>::connect_krylov_space_slot(
  const std::function<void(
    const internal::SolverGMRESImplementation::TmpVectors<VectorType> &)> &slot)
{
  return krylov_space_signal.connect(slot);
}



template <typename VectorType>
boost::signals2::connection
SolverGMRES<VectorType>::connect_re_orthogonalization_slot(
  const std::function<void(int)> &slot)
{
  return re_orthogonalize_signal.connect(slot);
}



template <typename VectorType>
double
SolverGMRES<VectorType>::criterion()
{
  // dummy implementation. this function is not needed for the present
  // implementation of gmres
  DEAL_II_ASSERT_UNREACHABLE();
  return 0;
}


//----------------------------------------------------------------------//

template <typename VectorType>
SolverFGMRES<VectorType>::SolverFGMRES(SolverControl            &cn,
                                       VectorMemory<VectorType> &mem,
                                       const AdditionalData     &data)
  : SolverBase<VectorType>(cn, mem)
  , additional_data(data)
{}



template <typename VectorType>
SolverFGMRES<VectorType>::SolverFGMRES(SolverControl        &cn,
                                       const AdditionalData &data)
  : SolverBase<VectorType>(cn)
  , additional_data(data)
{}



template <typename VectorType>
template <typename MatrixType, typename PreconditionerType>
void
SolverFGMRES<VectorType>::solve(const MatrixType         &A,
                                VectorType               &x,
                                const VectorType         &b,
                                const PreconditionerType &preconditioner)
{
  LogStream::Prefix prefix("FGMRES");

  SolverControl::State iteration_state = SolverControl::iterate;

  const unsigned int basis_size = additional_data.max_basis_size;

  // Generate an object where basis vectors are stored.
  typename internal::SolverGMRESImplementation::TmpVectors<VectorType> v(
    basis_size + 1, this->memory);
  typename internal::SolverGMRESImplementation::TmpVectors<VectorType> z(
    basis_size, this->memory);

  const bool delayed_reorthogonalization =
    additional_data.orthogonalization_strategy ==
    LinearAlgebra::OrthogonalizationStrategy::delayed_classical_gram_schmidt;

  // number of the present iteration; this number is not reset to zero upon a
  // restart
  unsigned int accumulated_iterations = 0;

  // matrix used for the orthogonalization process later
  H.reinit(basis_size + 1, basis_size);
  FullMatrix<double>                     H_orig(H);
  std::vector<std::pair<double, double>> givens_rotations(basis_size);
  Vector<double> h(delayed_reorthogonalization ? 2 * basis_size + 3 :
                                                 basis_size + 1);

  // Vectors for projected system
  Vector<double> projected_rhs(basis_size + 1);
  Vector<double> y(basis_size);

  // Iteration starts here
  double res = std::numeric_limits<double>::lowest();

  do
    {
      A.vmult(v(0, x), x);
      v[0].sadd(-1., 1., b);

      double norm_v   = v[0].l2_norm();
      res             = norm_v;
      iteration_state = this->iteration_status(accumulated_iterations, res, x);
      if (iteration_state == SolverControl::success)
        break;

      projected_rhs(0) = norm_v;
      if (norm_v != 0)
        v[0] /= norm_v;

      unsigned int inner_iteration = 0;
      for (; (inner_iteration < basis_size &&
              iteration_state == SolverControl::iterate);
           ++inner_iteration)
        {
          preconditioner.vmult(z(inner_iteration, x), v[inner_iteration]);
          A.vmult(v(inner_iteration + 1, x), z[inner_iteration]);

          // Gram-Schmidt
          bool re_orthogonalize = false;
          internal::SolverGMRESImplementation::iterated_gram_schmidt<
            VectorType>(additional_data.orthogonalization_strategy,
                        v,
                        inner_iteration + 1,
                        accumulated_iterations,
                        v[inner_iteration + 1],
                        h,
                        H,
                        H_orig,
                        re_orthogonalize);

          // Compute projected solution
          if (delayed_reorthogonalization)
            {
              if (inner_iteration > 0)
                internal::SolverGMRESImplementation::givens_rotation(
                  H, projected_rhs, givens_rotations, inner_iteration - 1);
              res = std::fabs(internal::SolverGMRESImplementation::
                                compute_givens_rotation_rhs(H,
                                                            projected_rhs,
                                                            givens_rotations,
                                                            inner_iteration));
            }
          else
            {
              internal::SolverGMRESImplementation::givens_rotation(
                H, projected_rhs, givens_rotations, inner_iteration);

              //  default residual
              res = std::fabs(projected_rhs(inner_iteration + 1));
            }

          // check convergence. note that the vector 'x' we pass to the
          // criterion is not the final solution we compute if we
          // decide to jump out of the iteration (we update 'x' again
          // right after the current loop)
          iteration_state =
            this->iteration_status(++accumulated_iterations, res, x);
        }

      // Solve triangular system with projected quantities and update solution
      // vector
      if (delayed_reorthogonalization)
        internal::SolverGMRESImplementation::givens_rotation(
          H, projected_rhs, givens_rotations, inner_iteration - 1);
      internal::SolverGMRESImplementation::solve_triangular(inner_iteration,
                                                            H,
                                                            projected_rhs,
                                                            y);
      dealii::internal::SolverGMRESImplementation::add(
        x, inner_iteration, y, z, false);
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
