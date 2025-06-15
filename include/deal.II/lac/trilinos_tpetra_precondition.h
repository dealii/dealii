// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_trilinos_tpetra_precondition_h
#define dealii_trilinos_tpetra_precondition_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/memory_space.h>

#include <deal.II/lac/vector.h>

#ifdef DEAL_II_TRILINOS_WITH_TPETRA

#  include <deal.II/base/enable_observer_pointer.h>

#  include <deal.II/lac/la_parallel_vector.h>
#  include <deal.II/lac/trilinos_tpetra_sparse_matrix.h>

#  include <Teuchos_BLAS_types.hpp>
#  include <Teuchos_ConfigDefs.hpp>
#  include <Teuchos_ParameterList.hpp>
#  include <Teuchos_RCPDecl.hpp>
#  include <Tpetra_Operator.hpp>


DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  namespace TpetraWrappers
  {

    /**
     * The base class for all preconditioners based on Tpetra sparse matrices.
     *
     * @ingroup TpetraWrappers
     * @ingroup Preconditioners
     */
    template <typename Number, typename MemorySpace = dealii::MemorySpace::Host>
    class PreconditionBase : public EnableObserverPointer
    {
    public:
      /**
       * @brief Standardized data struct to pipe additional flags to the
       * preconditioner.
       *
       */
      struct AdditionalData
      {};

      /**
       * @brief Constructor. Does not do anything. The <tt>initialize</tt> function of
       * the derived classes will have to create the preconditioner from a given
       * sparse matrix.
       *
       */
      PreconditionBase() = default;

      /**
       * @brief Destructor.
       * Destroys the preconditioner, leaving an object like just after having
       * called the constructor.
       */
      void
      clear();

      /**
       * @brief Apply the preconditioner.
       *
       * @param dst Input vector to apply the preconditioner to
       * @param src Result vector
       */
      virtual void
      vmult(Vector<Number, MemorySpace>       &dst,
            const Vector<Number, MemorySpace> &src) const;

      /**
       * @brief Apply the transpose preconditioner
       *
       * @param dst Input vector to apply the preconditioner to
       * @param src Result vector
       */
      virtual void
      Tvmult(Vector<Number, MemorySpace>       &dst,
             const Vector<Number, MemorySpace> &src) const;

      /**
       * @brief Apply the preconditioner
       *
       * @param dst Input vector to apply the preconditioner to
       * @param src Result vector
       */
      virtual void
      vmult(dealii::Vector<Number> &dst, dealii::Vector<Number> &src) const;

      /**
       * @brief Apply the transpose preconditioner
       */
      virtual void
      Tvmult(dealii::Vector<Number> &dst, dealii::Vector<Number> &src) const;

      /**
       * @brief Access to underlying Trilinos data
       *
       * Calling this function from an uninitialized object will cause an
       * exception.
       */
      const TpetraTypes::LinearOperator<Number, MemorySpace> &
      trilinos_operator() const;

      /**
       * @brief Access to underlying Trilinos data
       *
       * Calling this function from an uninitialized object will cause an
       * exception.
       */
      Teuchos::RCP<TpetraTypes::LinearOperator<Number, MemorySpace>>
      trilinos_rcp() const;

      /**
       * @name Partitioners
       */
      /** @{ */

      /**
       * Return the partitioning of the domain space of this matrix, i.e., the
       * partitioning of the vectors this matrix has to be multiplied with.
       *
       * @return IndexSet of the domain of this preconditioner
       */
      IndexSet
      locally_owned_domain_indices() const;

      /**
       * Return the partitioning of the range space of this matrix, i.e., the
       * partitioning of the vectors that are result from matrix-vector
       * products.
       *
       * @return IndexSet of the range of this preconditioner
       */
      IndexSet
      locally_owned_range_indices() const;

      /** @} */

      /**
       * @addtogroup Exceptions
       */
      /** @{ */
      /**
       * The maps of the underlying matrix and one of the vectors
       * do not match.
       */
      DeclException1(
        ExcNonMatchingMaps,
        std::string,
        << "The sparse matrix the preconditioner is based on "
        << "uses a map that is not compatible to the one in vector " << arg1
        << ". Check preconditioner and matrix setup.");

      /**
       * The chosen preconditioner does not support a transposed apply.
       */
      DeclExceptionMsg(
        ExcTransposeNotSupported,
        "The chosen preconditioner does not support transposing the matrix.");
      /** @} */

    protected:
      /**
       * This is a pointer to the preconditioner object that is used when
       * applying the preconditioner.
       */
      Teuchos::RCP<TpetraTypes::LinearOperator<Number, MemorySpace>>
        preconditioner;

      /**
       * @brief The list of preconditioner parameters.
       *
       * This structure is Trilinos counterpart to the AdditionalData structures
       * in deal.II. Therefore any initialize will at some point pass this to
       * the preconditioner. Most derived classes will handle building this list
       * based on an AdditionalData object that exposes and defaults the most
       * sensible parameters. But some classes will also offer full
       * customization for experienced Trilinos users.
       *
       */
      Teuchos::ParameterList parameter_list;
    };



#  ifdef DEAL_II_TRILINOS_WITH_IFPACK2
    /**
     * @brief Wrapper class for the IdentitySolver preconditioner of Ifpack2.
     *
     * @ingroup TpetraWrappers
     * @ingroup Preconditioners
     *
     * @tparam Number scalar type of the preconditioner
     * @tparam MemorySpace
     */
    template <typename Number, typename MemorySpace = dealii::MemorySpace::Host>
    class PreconditionIdentity : public PreconditionBase<Number, MemorySpace>
    {
    public:
      /**
       * @brief Construct identity preconditioner.
       *
       */
      PreconditionIdentity() = default;

      /**
       * @brief Initializes the preconditioner for the matrix <tt>A</tt>.
       *
       * Note that this only needs the matrix for information on the parallel
       * distribution of the matrix and will return the original vector on
       * <tt>vmult</tt> or <tt>Tvmult</tt>.
       *
       * @param A Matrix to base the preconditioner on.
       */
      void
      initialize(const SparseMatrix<Number, MemorySpace> &A);
    };



    /**
     * @brief The base class for all Ifpack2 preconditioners which are handled through its Factory.
     *
     * @ingroup TpetraWrappers
     * @ingroup Preconditioners
     */
    template <typename Number, typename MemorySpace = dealii::MemorySpace::Host>
    class PreconditionIfpackBase : public PreconditionBase<Number, MemorySpace>
    {
    public:
      /**
       * @brief Constructor.
       *
       * This constructor will only save the type name given.
       * The actual preconditioner has to be constructed through the
       * initialize method of one of the underlying classes.
       *
       * @param preconditioner_type
       */
      PreconditionIfpackBase(const std::string &preconditioner_type);

      /**
       * Initializes the preconditioner for the matrix <tt>A</tt> based on
       * the <tt>parameter_set</tt>.
       *
       */
      void
      initialize(const SparseMatrix<Number, MemorySpace> &A);

      /**
       * The chosen preconditioner is not supported or configured with Ifpack2.
       */
      DeclException1(
        ExcTrilinosIpack2PreconditionerUnsupported,
        std::string,
        << "You tried to select the preconditioner type <" << arg1 << ">\n"
        << "but this preconditioner is not supported by Trilinos/Ifpack22\n"
        << "due to one of the following reasons:\n"
        << "* This preconditioner does not exist\n"
        << "* This preconditioner has a specialized constructor not supported by the Ifpack2 Factory.\n"
        << "* This preconditioner is not (yet) supported by Trilinos/Ifpack2\n"
        << "* Trilinos/Ifpack2 was not configured for its use.");

    protected:
      /**
       * The set preconditioner type to be handed to
       * <tt>Ifpack2::Factory::create()</tt>
       */
      std::string preconditioner_type;
    };



    /**
     * @brief Wrapper to create custom Ifpack2 preconditioners.
     *
     * This general purpose class supports any preconditioner of Ifpack2
     * that can be constructed through its <tt>Ifpack2::Factory::create()</tt>
     * function.
     *
     * @note This is a class for users familiar with Trilinos.
     *
     *
     * @ingroup TpetraWrappers
     * @ingroup Preconditioners
     */
    template <typename Number, typename MemorySpace = dealii::MemorySpace::Host>
    class PreconditionIfpack
      : public PreconditionIfpackBase<Number, MemorySpace>
    {
    public:
      /**
       * @brief Construct a new custom Ifpack2 preconditioner.
       *
       * The currently tested options for the type are
       * <ul>
       * <li> "RELAXATION" </li>
       * <li> "CHEBYSHEV" </li>
       * <li> "RILUK" </li>
       * <li> "ILUT" </li>
       * <li> "SCHWARZ" </li>
       * </ul>
       *
       * But please refer to the User's Guide of Ifpack2 for more information.
       *
       * @param preconditioner_type the type based on Ifpack2 notation
       */
      PreconditionIfpack(const std::string &preconditioner_type);

      /**
       * @brief Set the parameter list for the preconditioner.
       *
       * This list will be passed to the Ifpack2 preconditioner during
       * initialization. For parameter options based on the chosen
       * preconditioner type see the User's Guide of Ifpack2.
       *
       * @param parameter_list
       */
      void
      set_parameter_list(Teuchos::ParameterList &parameter_list);
    };



    /**
     * @brief The classical Jacobi preconditioner within Ifpack2.
     *
     * @ingroup TpetraWrappers
     * @ingroup Preconditioners
     *
     */
    template <typename Number, typename MemorySpace = dealii::MemorySpace::Host>
    class PreconditionJacobi
      : public PreconditionIfpackBase<Number, MemorySpace>
    {
    public:
      /**
       * @brief The set of additional parameters to tune the preconditioner.
       *
       */
      struct AdditionalData
      {
        /**
         * @brief Constructor.
         *
         * Set the parameters to be used in the preconditioner.
         * The defaults yield a textbook Jacobi preconditioner.
         *
         * @param omega The damping factor for the relaxation.
         * @param fix_diagonal Fix small diagonal entries for the inversion?
         * @param min_diagonal Threshold for fix_diagonal.
         * @param n_sweeps Number of relaxation sweeps per call to vmult or Tvmult.
         */
        AdditionalData(const double omega        = 1.,
                       const bool   fix_diagonal = false,
                       const double min_diagonal = 0.,
                       const int    n_sweeps     = 1);

        /**
         * @brief Relaxation damping factor.
         *
         */
        double omega;
        /**
         * @brief Whether or not to enlarge entries below a threshold
         *
         * If the matrix diagonal contains entries close to zero, their
         * inversion will be close to division by zero. This can be corrected by
         * setting this value to <tt>true</tt>. This will result in additional
         * computational work during initialization.
         */
        bool fix_diagonal;
        /**
         * @brief Threshold below which entries will be fixed.
         *
         * If the threshold is zero (default) only entries which are exactly
         * zero will be replaced will small nonzero values.
         *
         */
        double min_diagonal;
        /**
         * @brief Set how often the preconditioner should be applied during vmult() or Tvmult().
         *
         */
        int n_sweeps;
      };
      /**
       * @brief Constructor.
       *
       */
      PreconditionJacobi();

      /**
       * @brief Compute the preconditioner based on the given matrix and parameters.
       *
       * @param A The matrix to base the preconditioner on.
       * @param additional_data The set of parameters to tune the preconditioner.
       */
      void
      initialize(const SparseMatrix<Number, MemorySpace> &A,
                 const AdditionalData &additional_data = AdditionalData());
    };



    /**
     * @brief l1 variant of the Jacobi preconditioner.
     *
     * This variant adds the l1-Norm of all off-processor entries of a local row
     * i to its diagonal entry d_ii to improve coupling in MPI-parallel
     * applications, as described and introduced in @cite BFKY2011.
     *
     * @ingroup TpetraWrappers
     * @ingroup Preconditioners
     */
    template <typename Number, typename MemorySpace = dealii::MemorySpace::Host>
    class PreconditionL1Jacobi
      : public PreconditionIfpackBase<Number, MemorySpace>
    {
    public:
      /**
       * @brief The set of additional parameters to tune the preconditioner.
       *
       */
      struct AdditionalData
      {
        /**
         * @brief Constructor.
         *
         * Set the parameters to be used in the preconditioner.
         * The defaults yield the l1 Jacobi preconditioner specified in the
         * original publication.
         *
         * @param omega The damping factor for the relaxation.
         * @param fix_diagonal Fix small diagonal entries for the inversion?
         * @param min_diagonal Threshold for fix_diagonal.
         * @param n_sweeps Number of relaxation sweeps per call to vmult or Tvmult.
         */
        AdditionalData(const double omega        = 1.,
                       const bool   fix_diagonal = false,
                       const double min_diagonal = 0.,
                       const int    n_sweeps     = 1);

        /**
         * @brief Relaxation damping factor.
         *
         */
        double omega;
        /**
         * @brief Whether or not to enlarge entries below a threshold
         *
         * If the matrix diagonal contains entries close to zero, their
         * inversion will be close to division by zero. This can be corrected by
         * setting this value to <tt>true</tt>. This will result in additional
         * computational work during initialization.
         */
        bool fix_diagonal;
        /**
         * @brief Threshold below which entries will be fixed.
         *
         * If the threshold is zero (default) only entries which are exactly
         * zero will be replaced will small nonzero values.
         *
         */
        double min_diagonal;
        /**
         * @brief Set how often the preconditioner should be applied during vmult() or Tvmult().
         *
         */
        int n_sweeps;
      };

      /**
       * @brief Constructor.
       *
       */
      PreconditionL1Jacobi();

      /**
       * @brief Compute the preconditioner based on the given matrix and parameters.
       *
       * @param A The matrix to base the preconditioner on.
       * @param additional_data The set of parameters to tune the preconditioner.
       */
      void
      initialize(const SparseMatrix<Number, MemorySpace> &A,
                 const AdditionalData &additional_data = AdditionalData());
    };



    /**
     * @brief l1 variant of the Gauss-Seidel preconditioner.
     *
     * This variant adds the l1-Norm of all off-processor entries of a single
     * row i of the upper triangular part of A to its diagonal entry d_ii to
     * improve coupling in MPI-parallel applications, , as described and
     * introduced in @cite BFKY2011.
     *
     * @ingroup TpetraWrappers
     * @ingroup Preconditioners
     */
    template <typename Number, typename MemorySpace = dealii::MemorySpace::Host>
    class PreconditionL1GaussSeidel
      : public PreconditionIfpackBase<Number, MemorySpace>
    {
    public:
      /**
       * @brief The set of additional parameters to tune the preconditioner.
       *
       */
      struct AdditionalData
      {
        /**
         * @brief Constructor.
         *
         * Set the parameters to be used in the preconditioner.
         * The defaults yield the l1 Gauss-Seidel preconditioner specified in
         * the original publication.
         *
         * @param omega The damping factor for the relaxation.
         * @param eta Threshold parameter for diagonal correction.
         * @param fix_diagonal Fix small diagonal entries for the inversion?
         * @param min_diagonal Threshold for fix_diagonal.
         * @param n_sweeps Number of relaxation sweeps per call to vmult or Tvmult.
         */
        AdditionalData(const double omega        = 1,
                       const double eta          = 1.5,
                       const bool   fix_diagonal = false,
                       const double min_diagonal = 0,
                       const int    n_sweeps     = 1);

        /**
         * @brief Relaxation damping factor.
         *
         */
        double omega;

        /**
         * @brief Threshold parameter for diagonal correction.
         *
         * The l1 correction for a row i will only be added if it is more than
         * eta times larger than the original diagonal entry.
         *
         */
        double eta;
        /**
         * @brief Whether or not to enlarge entries below a threshold
         *
         * If the matrix diagonal contains entries close to zero, their
         * inversion will be close to division by zero. This can be corrected by
         * setting this value to <tt>true</tt>. This will result in additional
         * computational work during initialization.
         */
        bool fix_diagonal;
        /**
         * @brief Threshold below which entries will be fixed.
         *
         * If the threshold is zero (default) only entries which are exactly
         * zero will be replaced will small nonzero values.
         *
         */
        double min_diagonal;
        /**
         * @brief Set how often the preconditioner should be applied during vmult() or Tvmult().
         *
         */
        int n_sweeps;
      };

      /**
       * @brief Constructor.
       *
       */
      PreconditionL1GaussSeidel();

      /**
       * @brief Compute the preconditioner based on the given matrix and parameters.
       *
       * @param A The matrix to base the preconditioner on.
       * @param additional_data The set of parameters to tune the preconditioner.
       */
      void
      initialize(const SparseMatrix<Number, MemorySpace> &A,
                 const AdditionalData &additional_data = AdditionalData());
    };



    /**
     * @brief The class for the SOR preconditioner within Ifpack2.
     *
     * If the code is executed MPI-parallel the individual processors will be
     * connected with an additive Schwarz method.
     *
     * @ingroup TpetraWrappers
     * @ingroup Preconditioners
     */
    template <typename Number, typename MemorySpace = dealii::MemorySpace::Host>
    class PreconditionSOR : public PreconditionIfpackBase<Number, MemorySpace>
    {
    public:
      /**
       * @brief The set of additional parameters to tune the preconditioner.
       *
       */
      struct AdditionalData
      {
        /**
         * @brief Constructor.
         *
         * Set the parameters to be used in the preconditioner.
         * The defaults yield a textbook Jacobi preconditioner.
         *
         * @param omega The damping factor for the relaxation.
         * @param overlap Overlap between processor local matrices.
         * @param fix_diagonal Fix small diagonal entries for the inversion?
         * @param min_diagonal Threshold for fix_diagonal.
         * @param n_sweeps Number of relaxation sweeps per call to vmult or Tvmult.
         */
        AdditionalData(const double omega        = 1.,
                       const int    overlap      = 0,
                       const bool   fix_diagonal = false,
                       const double min_diagonal = 0.,
                       const int    n_sweeps     = 1);

        /**
         * @brief Relaxation damping factor.
         *
         */
        double omega;
        /**
         * @brief Overlap between processor local matrices.
         *
         * The amount of overlap between the processor local matrix blocks
         * used in the underlying additive Schwarz method.
         *
         * The default will yield a block Jacobi preconditioner with each
         * processor forming its own local block.
         */
        int overlap;
        /**
         * @brief Whether or not to enlarge entries below a threshold
         *
         * If the matrix diagonal contains entries close to zero, their
         * inversion will be close to division by zero. This can be corrected by
         * setting this value to <tt>true</tt>. This will result in additional
         * computational work during initialization.
         */
        bool fix_diagonal;
        /**
         * @brief Threshold below which entries will be fixed.
         *
         * If the threshold is zero (default) only entries which are exactly
         * zero will be replaced will small nonzero values.
         *
         */
        double min_diagonal;
        /**
         * @brief Set how often the preconditioner should be applied during vmult() or Tvmult().
         *
         */
        int n_sweeps;
      };

      /**
       * @brief Constructor.
       *
       */
      PreconditionSOR();

      /**
       * @brief Compute the preconditioner based on the given matrix and parameters.
       *
       * @param A The matrix to base the preconditioner on.
       * @param additional_data The set of parameters to tune the preconditioner.
       */
      void
      initialize(const SparseMatrix<Number, MemorySpace> &A,
                 const AdditionalData &additional_data = AdditionalData());
    };



    /**
     * The class for the classical SSOR preconditioner within Ifpack2.
     *
     * @ingroup TpetraWrappers
     * @ingroup Preconditioners
     */
    template <typename Number, typename MemorySpace = dealii::MemorySpace::Host>
    class PreconditionSSOR : public PreconditionIfpackBase<Number, MemorySpace>
    {
    public:
      /**
       * @brief The set of additional parameters to tune the preconditioner.
       *
       */
      struct AdditionalData
      {
        /**
         * @brief Constructor.
         *
         * Set the parameters to be used in the preconditioner.
         * The defaults yield a textbook Jacobi preconditioner.
         *
         * @param omega The damping factor for the relaxation.
         * @param overlap Overlap between processor local matrices.
         * @param fix_diagonal Fix small diagonal entries for the inversion?
         * @param min_diagonal Threshold for fix_diagonal.
         * @param n_sweeps Number of relaxation sweeps per call to vmult or Tvmult.
         */
        AdditionalData(const double omega        = 1.,
                       const int    overlap      = 0,
                       const bool   fix_diagonal = false,
                       const double min_diagonal = 0.,
                       const int    n_sweeps     = 1);

        /**
         * @brief Relaxation damping factor.
         *
         */
        double omega;
        /**
         * @brief Overlap between processor local matrices.
         *
         * The amount of overlap between the processor local matrix blocks
         * used in the underlying additive Schwarz method.
         *
         * The default will yield a block Jacobi preconditioner with each
         * processor forming its own local block.
         */
        int overlap;
        /**
         * @brief Whether or not to enlarge entries below a threshold
         *
         * If the matrix diagonal contains entries close to zero, their
         * inversion will be close to division by zero. This can be corrected by
         * setting this value to <tt>true</tt>. This will result in additional
         * computational work during initialization.
         */
        bool fix_diagonal;
        /**
         * @brief Threshold below which entries will be fixed.
         *
         * If the threshold is zero (default) only entries which are exactly
         * zero will be replaced will small nonzero values.
         *
         */
        double min_diagonal;
        /**
         * @brief Set how often the preconditioner should be applied during vmult() or Tvmult().
         *
         */
        int n_sweeps;
      };


      /**
       * @brief Constructor.
       *
       */
      PreconditionSSOR();

      /**
       * @brief Compute the preconditioner based on the given matrix and parameters.
       *
       * @param A The matrix to base the preconditioner on.
       * @param additional_data The set of parameters to tune the preconditioner.
       */
      void
      initialize(const SparseMatrix<Number, MemorySpace> &A,
                 const AdditionalData &additional_data = AdditionalData());
    };



    /**
     * The class for the Chebyshev preconditioner within Ifpack2.
     *
     * @ingroup TpetraWrappers
     * @ingroup Preconditioners
     */
    template <typename Number, typename MemorySpace = dealii::MemorySpace::Host>
    class PreconditionChebyshev
      : public PreconditionIfpackBase<Number, MemorySpace>
    {
    public:
      /**
       * @brief The set of additional parameters to tune the preconditioner.
       *
       */
      struct AdditionalData
      {
        /**
         * @brief Constructor.
         *
         * @param degree Degree of the Chebyshev polynomial.
         * @param max_eigenvalue Upper bound for the maximum eigenvalue of the matrix.
         * @param min_eigenvalue Lower bound for the minimum eigenvalue of the matrix.
         * @param eigenvalue_ratio Ratio between maximum and minimum eigenvalue.
         * @param min_diagonal Threshold for increasing diagonal entries.
         * @param nonzero_starting Do not zero starting entries of solution.
         */
        AdditionalData(const int    degree           = 1,
                       const double max_eigenvalue   = 10.,
                       const double min_eigenvalue   = 1.,
                       const double eigenvalue_ratio = 30.,
                       const double min_diagonal     = 1e-12,
                       const bool   nonzero_starting = false);
        /**
         * @brief Degree of the Chebyshev polynomial.
         *
         * The degree directly corresponds to the number of matrix-vector
         * products that have to be performed during a single application of the
         * vmult() and Tvmult() operations.
         *
         */
        int degree;
        /*
         * @brief Upper bound for the maximum eigenvalue of the matrix.
         *
         * This needs to be set properly for the appropriate performance of this
         * preconditioner.
         */
        double max_eigenvalue;
        /*
         * @brief Lower bound for the minimum eigenvalue of the matrix.
         *
         */
        double min_eigenvalue;
        /**
         * @brief Estimated ratio between maximum and minimum eigenvalue.
         *
         */
        double eigenvalue_ratio;
        /**
         * @brief Threshold below which entries will be fixed.
         *
         * If the threshold is zero (default) only entries which are exactly
         * zero will be replaced will small nonzero values.
         *
         */
        double min_diagonal;
        /**
         * @brief Do not zero starting entries of solution vector.
         *
         * The default (false) zeroes out the entries of dst during vmult() and
         * Tvmult() which is the recommended setting.
         *
         * However, in some situations (e.g. high-frequency error smoothing) it
         * can be useful to append previous data to the Chebyshev corrections.
         * The user should really know what they are doing when setting this
         * flag to true.
         *
         */
        bool nonzero_starting;
      };

      /**
       * @brief Constructor.
       *
       */
      PreconditionChebyshev();

      /**
       * @brief Compute the preconditioner based on the given matrix and parameters.
       *
       * @param A The matrix to base the preconditioner on.
       * @param additional_data The set of parameters to tune the preconditioner.
       */
      void
      initialize(const SparseMatrix<Number, MemorySpace> &A,
                 const AdditionalData &additional_data = AdditionalData());
    };



    /**
     * @brief The ILU/ILU(K)/RILU(K) preconditioner.
     *
     * The class for the modified incomplete LU factorization preconditioner
     * RILUK within Ifpack2. This preconditioner works both in serial and
     * parallel, depending on the matrix it is based on.
     *
     * In general, an incomplete factorization only has values on the sparsity
     * pattern of the matrix, but one can gradually increase <tt>ilu_fill</tt>
     * to eventually obtain a full LU factorization.
     *
     * For parallel applications this preconditioner implements an additive
     * Schwarz preconditioner with RILUK as local smoothing on each processor.
     *
     * @ingroup TpetraWrappers
     * @ingroup Preconditioners
     */
    template <typename Number, typename MemorySpace = dealii::MemorySpace::Host>
    class PreconditionILU : public PreconditionIfpackBase<Number, MemorySpace>
    {
    public:
      /**
       * @brief The set of additional parameters to tune the preconditioner.
       *
       */
      struct AdditionalData
      {
        /**
         * @brief Constructor.
         *
         * Before factorization each diagonal entry will be modified by the
         * following formula \f[a_{ii}^{new} = \alpha\sign(a_{ii})+\beta a_{ii}
         * \f] with \f[\alpha\f] given by ilu_atol and \f[\beta\f] given by
         * ilu_rtol.
         *
         * @param ilu_fill Amount of additional fill-in.
         * @param ilu_atol Constant to be to each diagonal entry before factorization.
         * @param ilu_rtol Factor to scale all diagonal entries by before factorization.
         * @param overlapOverlap between processor local matrices.
         */
        AdditionalData(const int    ilu_fill = 0,
                       const double ilu_atol = 0.,
                       const double ilu_rtol = 1.,
                       const int    overlap  = 0);

        /**
         * @brief Amount of additional fill-in.
         *
         * Level-of-fill to increase the sparsity pattern of the preconditioner.
         * A large enough value will result in a complete LU factorization.
         *
         */
        int ilu_fill;
        /**
         * @brief Constant to be added to each diagonal entry before factorization.
         *
         */
        double ilu_atol;
        /**
         * @brief Factor to scale all diagonal entries by before factorization.
         *
         */
        double ilu_rtol;
        /**
         * @brief Overlap between processor local matrices.
         *
         * The amount of overlap between the processor local matrix blocks
         * used in the underlying additive Schwarz method.
         *
         * The default will yield a block Jacobi preconditioner with each
         * processor forming its own local block.
         */
        int overlap;
      };

      /**
       * @brief Constructor.
       *
       */
      PreconditionILU();

      /**
       * @brief Compute the preconditioner based on the given matrix and parameters.
       *
       * @param A The matrix to base the preconditioner on.
       * @param additional_data The set of parameters to tune the preconditioner.
       */
      void
      initialize(const SparseMatrix<Number, MemorySpace> &A,
                 const AdditionalData &additional_data = AdditionalData());
    };



    /**
     * The class for the thresholded incomplete LU (ILUT) preconditioner within
     * Ifpack2.
     *
     * @ingroup TpetraWrappers
     * @ingroup Preconditioners
     */
    template <typename Number, typename MemorySpace = dealii::MemorySpace::Host>
    class PreconditionILUT : public PreconditionIfpackBase<Number, MemorySpace>
    {
    public:
      /**
       * @brief The set of additional parameters to tune the preconditioner.
       *
       */
      struct AdditionalData
      {
        /**
         * @brief Constructor.
         *
         * Before factorization each diagonal entry will be modified by the
         * following formula \f[a_{ii}^{new} = \alpha\sign(a_{ii})+\beta a_{ii}
         * \f] with \f[\alpha\f] given by ilut_atol and \f[\beta\f] given by
         * ilut_rtol.
         *
         * @param ilut_drop Threshold for dropping entries.
         * @param ilut_fill Amount of additional fill-in.
         * @param ilut_atol Constant to be to each diagonal entry before factorization.
         * @param ilut_rtol Factor to scale all diagonal entries by before factorization.
         * @param overlapOverlap between processor local matrices.
         */
        AdditionalData(const double ilut_drop = 0.,
                       const double ilut_fill = 0.,
                       const double ilut_atol = 0.,
                       const double ilut_rtol = 1.,
                       const int    overlap   = 0);

        /**
         * @brief Threshold for dropping entries.
         *
         * Together with <tt>ilut_fill</tt> this controls the amount of fill-in
         * and the actual values to be used or dropped.
         */
        double ilut_drop;
        /**
         * @brief Amount of additional fill-in.
         *
         * Level-of-fill to increase the sparsity pattern of the preconditioner.
         * A large enough value will result in a complete LU factorization.
         * Note, however, that this will drastically increase the memory
         * requirement, especially for 3d simulations.
         *
         */
        double ilut_fill;
        /**
         * @brief Constant to be added to each diagonal entry before factorization.
         *
         */
        double ilut_atol;
        /**
         * @brief Factor to scale all diagonal entries by before factorization.
         *
         */
        double ilut_rtol;
        /**
         * @brief Overlap between processor local matrices.
         *
         * The amount of overlap between the processor local matrix blocks
         * used in the underlying additive Schwarz method.
         *
         * The default will yield a block Jacobi preconditioner with each
         * processor forming its own local block.
         */
        int overlap;
      };



      /**
       * @brief Constructor.
       *
       */
      PreconditionILUT();

      /**
       * @brief Compute the preconditioner based on the given matrix and parameters.
       *
       * @param A The matrix to base the preconditioner on.
       * @param additional_data The set of parameters to tune the preconditioner.
       */
      void
      initialize(const SparseMatrix<Number, MemorySpace> &A,
                 const AdditionalData &additional_data = AdditionalData());
    };



    /**
     * The class for the BlockJacobi preconditioner within Ifpack2.
     *
     * @note This preconditioner always uses linear partitioning.
     * There are other partitioners available that need additional user provided
     * data, so we suggest using PreconditionIfpack with a custom
     * Teuchos::ParameterList. Details on the parameters can be found in the
     * User's Guide of Ifpack2.
     *
     *
     * @ingroup TpetraWrappers
     * @ingroup Preconditioners
     */
    template <typename Number, typename MemorySpace = dealii::MemorySpace::Host>
    class PreconditionBlockJacobi
      : public PreconditionIfpackBase<Number, MemorySpace>
    {
    public:
      /**
       * @brief The set of additional parameters to tune the preconditioner.
       *
       */
      struct AdditionalData
      {
        /**
         * @brief Constructor.
         *
         * Set the parameters to be used in the preconditioner.
         * The defaults yield a textbook Jacobi preconditioner.
         *
         * @param n_local_parts The number of blocks per processor.
         * @param omega The damping factor for the relaxation.
         * @param block_overlap Amount of overlap between blocks
         * @param n_sweeps Number of relaxation sweeps per call to vmult or Tvmult
         */
        AdditionalData(const int    n_local_parts = 1,
                       const double omega         = 1.,
                       const int    block_overlap = 0,
                       const int    n_sweeps      = 1);

        /**
         * @brief Number of blocks per processor.
         *
         * The default of 1 results in a single block per processor
         * which corresponds to Preconditioner from PreconditionJacobi
         */
        int n_local_parts;
        /**
         * @brief Relaxation damping factor.
         *
         */
        double omega;
        /**
         * @brief Amount of overlap between blocks
         *
         */
        int block_overlap;
        /**
         * @brief Set how often the preconditioner should be applied during vmult() or Tvmult().
         *
         */
        int n_sweeps;
      };


      /**
       * @brief Constructor.
       *
       */
      PreconditionBlockJacobi();

      /**
       * @brief Compute the preconditioner based on the given matrix and parameters.
       *
       * @param A The matrix to base the preconditioner on.
       * @param additional_data The set of parameters to tune the preconditioner.
       */
      void
      initialize(const SparseMatrix<Number, MemorySpace> &A,
                 const AdditionalData &additional_data = AdditionalData());
    };



    /**
     * @brief The class for the Block SOR preconditioner within Ifpack2.
     *
     * If the code is executed MPI-parallel the individual processors will be
     * connected with an additive Schwarz method.
     *
     * @ingroup TpetraWrappers
     * @ingroup Preconditioners
     */
    template <typename Number, typename MemorySpace = dealii::MemorySpace::Host>
    class PreconditionBlockSOR
      : public PreconditionIfpackBase<Number, MemorySpace>
    {
    public:
      /**
       * @brief The set of additional parameters to tune the preconditioner.
       *
       */
      struct AdditionalData
      {
        /**
         * @brief Constructor.
         *
         * Set the parameters to be used in the preconditioner.
         *
         * @param n_local_parts The number of blocks per processor.
         * @param omega The damping factor for the relaxation.
         * @param overlap Overlap between processors
         * @param n_sweeps Number of relaxation sweeps per call to vmult or Tvmult
         */
        AdditionalData(const int    n_local_parts = 1,
                       const double omega         = 1,
                       const int    overlap       = 0,
                       const int    n_sweeps      = 1);

        /**
         * @brief Number of blocks per processor.
         *
         * The default of 1 results in a single block per processor
         * which corresponds to Preconditioner from PreconditionJacobi
         */
        int n_local_parts;
        /**
         * @brief Relaxation damping factor.
         *
         */
        double omega;
        /**
         * @brief Overlap between processor local matrices.
         *
         * The amount of overlap between the processor local matrix blocks
         * used in the underlying additive Schwarz method.
         *
         * The default will yield a block Jacobi preconditioner with each
         * processor forming its own local block.
         */
        int overlap;
        /**
         * @brief Set how often the preconditioner should be applied during vmult() or Tvmult().
         *
         */
        int n_sweeps;
      };

      PreconditionBlockSOR();

      void
      initialize(const SparseMatrix<Number, MemorySpace> &A,
                 const AdditionalData &additional_data = AdditionalData());
    };



    /**
     * @brief The class for the Block SSOR preconditioner within Ifpack2.
     *
     * If the code is executed MPI-parallel the individual processors will be
     * connected with an additive Schwarz method.
     *
     * @ingroup TpetraWrappers
     * @ingroup Preconditioners
     */
    template <typename Number, typename MemorySpace = dealii::MemorySpace::Host>
    class PreconditionBlockSSOR
      : public PreconditionIfpackBase<Number, MemorySpace>
    {
    public:
      /**
       * @brief The set of additional parameters to tune the preconditioner.
       *
       */
      struct AdditionalData
      {
        /**
         * @brief Constructor.
         *
         * Set the parameters to be used in the preconditioner.
         *
         * @param n_local_parts The number of blocks per processor.
         * @param omega The damping factor for the relaxation.
         * @param overlap Overlap between processors
         * @param n_sweeps Number of relaxation sweeps per call to vmult or Tvmult
         */
        AdditionalData(const int    n_local_parts = 1,
                       const double omega         = 1,
                       const int    overlap       = 1,
                       const int    n_sweeps      = 1);
        /**
         * @brief Number of blocks per processor.
         *
         * The default of 1 results in a single block per processor
         * which corresponds to Preconditioner from PreconditionJacobi
         */
        int n_local_parts;
        /**
         * @brief Relaxation damping factor.
         *
         */
        double omega;
        /**
         * @brief Overlap between processor local matrices.
         *
         * The amount of overlap between the processor local matrix blocks
         * used in the underlying additive Schwarz method.
         *
         * The default will yield a block Jacobi preconditioner with each
         * processor forming its own local block.
         */
        int overlap;
        /**
         * @brief Set how often the preconditioner should be applied during vmult() or Tvmult().
         *
         */
        int n_sweeps;
      };

      PreconditionBlockSSOR();

      void
      initialize(const SparseMatrix<Number, MemorySpace> &A,
                 const AdditionalData &additional_data = AdditionalData());
    };
#    ifdef DEAL_II_TRILINOS_WITH_TPETRA_MUELU

    /**
     * @brief The class for the algebraic multigrid preconditioner MueLu.
     *
     * @ingroup TpetraWrappers
     * @ingroup Preconditioners
     */
    template <typename Number, typename MemorySpace = dealii::MemorySpace::Host>
    class PreconditionAMGMueLu : public PreconditionBase<Number, MemorySpace>
    {
    public:
      /**
       * The chosen smoother is not supported or configured with Ifpack2.
       */
      DeclException1(
        ExcTrilinosMueLuSmootherUnsupported,
        std::string,
        << "You tried to select the smoother type <" << arg1 << ">\n"
        << "but this is not supported by our interface\n"
        << "due to one of the following reasons:\n"
        << "* This smoother does not exist\n"
        << "* This smoother is not (yet) supported by our interface\n"
        << "* This smoother is not (yet) supported by MueLu.\n");
      DeclException1(
        ExcTrilinosMueLuCoarseSolverUnsupported,
        std::string,
        << "You tried to select the coarse solver type <" << arg1 << ">\n"
        << "but this is not supported by our interface\n"
        << "due to one of the following reasons:\n"
        << "* This coarse solver does not exist\n"
        << "* This coarse solver is not (yet) supported by our interface\n"
        << "* This coarse solver is not (yet) supported by MueLu.\n");

      /**
       * @brief The set of additional parameters to tune the preconditioner.
       *
       */
      struct AdditionalData
      {
        /**
         * @brief Constructor.
         *
         * By default, we pretend to work on elliptic problems with linear
         * finite elements on a scalar equation.
         *
         * @param elliptic Optimize MueLu for elliptic problems?
         * @param symmetric Assume for A to be symmetric?
         * @param w_cycle Use W-cycle instead of V-cycle?
         * @param aggregation_threshold Threshold for coarsening.
         * @param smoother_sweeps Number of times to apply the smoother.
         * @param smoother_overlap Overlap of smoother if run in parallel.
         * @param output_details Print additional info to screen?
         * @param smoother_type Determine the smoother to use.
         * @param coarse_type Determine the coarse solver.
         */
        AdditionalData(const bool         elliptic              = true,
                       const bool         symmetric             = true,
                       const bool         w_cycle               = false,
                       const double       aggregation_threshold = 1e-4,
                       const int          smoother_sweeps       = 2,
                       const int          smoother_overlap      = 0,
                       const bool         output_details        = false,
                       const std::string &smoother_type         = "Chebyshev",
                       const std::string &coarse_type           = "KLU2");
        /**
         * @brief Optimize for elliptic problems?
         *
         * Determines what underlying multigrid algorithm should be used.
         * For elliptic problems this would be smoothed aggregation
         * and unsmoothed for non-elliptic problems.
         *
         */
        bool elliptic;
        /**
         * @brief Assume A is symmetric?
         *
         * If A is symmetric the prolongation and restriction operators are
         * transposes of each other.
         */
        bool symmetric;
        /**
         * @brief Use W cycle instead of V cycle.
         *
         * Use W-cycle instead of the default V-cycle.
         */
        bool w_cycle;
        /**
         * @brief Threshold for coarsening.
         *
         * This threshold tells the AMG setup how the coarsening should be
         * performed. In MueLu all points that strongly couple with the
         * tentative coarse-level point form one aggregate.
         * The strong coupling is controlled by aggregation_threshold,
         * meaning that all elements that are larger than aggregation_threshold
         * times the diagonal element couple strongly.
         */
        double aggregation_threshold;
        /**
         * @brief Number of times pre- and post-smoothing is applied.
         *
         * For a Chebyshev smoother this determines the polynomial degree.
         * Otherwise it determines how often the matrix-vector product of the
         * smoother is applied.
         */
        int smoother_sweeps;
        /**
         * @brief Determine the overlap in the smoother when run in parallel.
         *
         */
        int smoother_overlap;
        /**
         * @brief Print internal information details of MueLu to screen.
         *
         * This can be useful for debugging.
         *
         */
        bool output_details;
        /**
         * @brief Determines which smoother to use for the AMG cycle
         * The following possibilities are implemented in this interface.
         * From Ifpack2:
         * <ul>
         * <li> "Jacobi" </li>
         * <li> "l1 Jacobi" </li>
         * <li> "Gauss Seidel" </li>
         * <li> "l1 Gauss Seidel" </li>
         * <li> "Symmetric Gauss Seidel" </li>
         * <li> "Chebyshev" </li>
         * <li> "ILU" </li>
         * <li> "ILUT" </li>
         * </ul>
         * From Amesos2:
         * <ul>
         * <li> "KLU2" </li>
         * <li> "SuperLU" </li>
         * <li> "SuperLU_dist" </li>
         * </ul>
         *
         */
        std::string smoother_type;
        /**
         * @brief Type of smoother to use on the finer levels.
         *
         * The same settings as for the smoother type are possible.
         */
        std::string coarse_type;
      };

      PreconditionAMGMueLu() = default;

      /**
       * @brief Full control interface for constructing a MueLu preconditioner.
       *
       * MueLu allows the construction of very complex Multigrid hierarchies,
       * through a complex Teuchos::ParameterList.
       * This function allows anyone familiar with MueLu to directly pass such a
       * list to MueLu.
       *
       * If you want to learn more about constructing such hierarchies both
       * the MueLu User Guide and Tutorial are highly recommended.
       *
       */
      void
      initialize(SparseMatrix<Number, MemorySpace> &A,
                 Teuchos::ParameterList            &parameters);

      /**
       * @brief Simplified interface for constructing a MueLu preconditioner
       *
       * MueLu allows the construction of very complex Multigrid hierarchies,
       * through a complex Teuchos::ParameterList.
       *
       * This function constructs such a list internally based on the given
       * additional_data. However, if you want full control over the hierarchy,
       * you should use the other initialize function.
       */
      void
      initialize(SparseMatrix<Number, MemorySpace> &A,
                 const AdditionalData &additional_data = AdditionalData());
    };
#    endif
#  endif // DEAL_II_TRILINOS_WITH_IFPACK2
  }      // namespace TpetraWrappers
} // namespace LinearAlgebra


DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_TRILINOS_WITH_TPETRA

#endif
