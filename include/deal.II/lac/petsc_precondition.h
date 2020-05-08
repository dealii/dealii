// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2020 by the deal.II authors
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

#ifndef dealii_petsc_precondition_h
#  define dealii_petsc_precondition_h


#  include <deal.II/base/config.h>

#  ifdef DEAL_II_WITH_PETSC

#    include <deal.II/lac/exceptions.h>

#    include <petscpc.h>

DEAL_II_NAMESPACE_OPEN



namespace PETScWrappers
{
  // forward declarations
#    ifndef DOXYGEN
  class MatrixBase;
  class VectorBase;
  class SolverBase;
#    endif

  /**
   * Base class for preconditioner classes using the PETSc functionality. The
   * classes in this hierarchy don't do a whole lot, except for providing a
   * function that sets the preconditioner and certain parameters on the
   * preconditioning context of the solver. These classes are basically here
   * only to allow a similar interface as already used for the deal.II solver
   * and preconditioner classes.
   *
   * Note that derived classes only provide interfaces to the relevant
   * functionality of PETSc. PETSc does not implement all preconditioners for
   * all matrix types. In particular, some preconditioners are not going to
   * work for parallel jobs, such as for example the ILU preconditioner.
   *
   * @ingroup PETScWrappers
   * @author Wolfgang Bangerth, Timo Heister, 2004, 2011
   */
  class PreconditionerBase
  {
  public:
    /**
     * Constructor.
     */
    PreconditionerBase();

    /**
     * Destructor.
     */
    virtual ~PreconditionerBase();

    /**
     * Destroys the preconditioner, leaving an object like just after having
     * called the constructor.
     */
    void
    clear();

    /**
     * Apply the preconditioner once to the given src vector.
     */
    void
    vmult(VectorBase &dst, const VectorBase &src) const;

    /**
     * Apply the transpose preconditioner once to the given src vector.
     */
    void
    Tvmult(VectorBase &dst, const VectorBase &src) const;


    /**
     * Give access to the underlying PETSc object.
     */
    const PC &
    get_pc() const;

  protected:
    /**
     * the PETSc preconditioner object
     */
    PC pc;

    /**
     * A pointer to the matrix that acts as a preconditioner.
     */
    Mat matrix;

    /**
     * Internal function to create the PETSc preconditioner object. Fails if
     * called twice.
     */
    void
    create_pc();

    /**
     * Conversion operator to get a representation of the matrix that
     * represents this preconditioner. We use this inside the actual solver,
     * where we need to pass this matrix to the PETSc solvers.
     */
    operator Mat() const;

    // Make the solver class a friend, since it needs to call the conversion
    // operator.
    friend class SolverBase;
  };



  /**
   * A class that implements the interface to use the PETSc Jacobi
   * preconditioner.
   *
   * See the comment in the base class
   * @ref PreconditionerBase
   * for when this preconditioner may or may not work.
   *
   * @ingroup PETScWrappers
   * @author Wolfgang Bangerth, Timo Heister, 2004, 2011
   */
  class PreconditionJacobi : public PreconditionerBase
  {
  public:
    /**
     * Standardized data struct to pipe additional flags to the
     * preconditioner.
     */
    struct AdditionalData
    {};

    /**
     * Empty Constructor. You need to call initialize() before using this
     * object.
     */
    PreconditionJacobi() = default;


    /**
     * Constructor. Take the matrix which is used to form the preconditioner,
     * and additional flags if there are any.
     */
    PreconditionJacobi(
      const MatrixBase &    matrix,
      const AdditionalData &additional_data = AdditionalData());

    /**
     * Same as above but without setting a matrix to form the preconditioner.
     * Intended to be used with SLEPc objects.
     */
    PreconditionJacobi(
      const MPI_Comm        communicator,
      const AdditionalData &additional_data = AdditionalData());

    /**
     * Initialize the preconditioner object and calculate all data that is
     * necessary for applying it in a solver. This function is automatically
     * called when calling the constructor with the same arguments and is only
     * used if you create the preconditioner without arguments.
     */
    void
    initialize(const MatrixBase &    matrix,
               const AdditionalData &additional_data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular preconditioner.
     */
    AdditionalData additional_data;

    /**
     * Initialize the preconditioner object without knowing a particular
     * matrix. This function sets up appropriate parameters to the underlying
     * PETSc object after it has been created.
     */
    void
    initialize();
  };



  /**
   * A class that implements the interface to use the PETSc Block Jacobi
   * preconditioner. PETSc defines the term "block Jacobi" as a preconditioner
   * in which it looks at a number of diagonal blocks of the matrix and then
   * defines a preconditioner in which the preconditioner matrix has the same
   * block structure as only these diagonal blocks, and each diagonal block
   * of the preconditioner is an approximation of the inverse of the
   * corresponding block of the original matrix.
   * The blocking structure of the matrix is determined by the
   * association of degrees of freedom to the individual processors in an
   * MPI-parallel job. If you use this preconditioner on a sequential job (or an
   * MPI job with only one process) then the entire matrix is the only block.
   *
   * By default, PETSc uses an ILU(0) decomposition of each diagonal block of
   * the matrix for preconditioning. This can be changed, as is explained in
   * the relevant section of the PETSc manual, but is not implemented here.
   *
   * See the comment in the base class
   * @ref PreconditionerBase
   * for when this preconditioner may or may not work.
   *
   * @ingroup PETScWrappers
   * @author Wolfgang Bangerth, Timo Heister, 2004, 2011
   */
  class PreconditionBlockJacobi : public PreconditionerBase
  {
  public:
    /**
     * Standardized data struct to pipe additional flags to the
     * preconditioner.
     */
    struct AdditionalData
    {};

    /**
     * Empty Constructor. You need to call initialize() before using this
     * object.
     */
    PreconditionBlockJacobi() = default;

    /**
     * Constructor. Take the matrix which is used to form the preconditioner,
     * and additional flags if there are any.
     */
    PreconditionBlockJacobi(
      const MatrixBase &    matrix,
      const AdditionalData &additional_data = AdditionalData());

    /**
     * Same as above but without setting a matrix to form the preconditioner.
     * Intended to be used with SLEPc objects.
     */
    PreconditionBlockJacobi(
      const MPI_Comm        communicator,
      const AdditionalData &additional_data = AdditionalData());


    /**
     * Initialize the preconditioner object and calculate all data that is
     * necessary for applying it in a solver. This function is automatically
     * called when calling the constructor with the same arguments and is only
     * used if you create the preconditioner without arguments.
     */
    void
    initialize(const MatrixBase &    matrix,
               const AdditionalData &additional_data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular preconditioner.
     */
    AdditionalData additional_data;

    /**
     * Initialize the preconditioner object without knowing a particular
     * matrix. This function sets up appropriate parameters to the underlying
     * PETSc object after it has been created.
     */
    void
    initialize();
  };



  /**
   * A class that implements the interface to use the PETSc SOR
   * preconditioner.
   *
   * @note Only works in serial with a PETScWrappers::SparseMatrix.
   *
   * @ingroup PETScWrappers
   * @author Wolfgang Bangerth, Timo Heister, 2004, 2011
   */
  class PreconditionSOR : public PreconditionerBase
  {
  public:
    /**
     * Standardized data struct to pipe additional flags to the
     * preconditioner.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default, set the damping parameter to one.
       */
      AdditionalData(const double omega = 1);

      /**
       * Relaxation parameter.
       */
      double omega;
    };

    /**
     * Empty Constructor. You need to call initialize() before using this
     * object.
     */
    PreconditionSOR() = default;

    /**
     * Constructor. Take the matrix which is used to form the preconditioner,
     * and additional flags if there are any.
     */
    PreconditionSOR(const MatrixBase &    matrix,
                    const AdditionalData &additional_data = AdditionalData());

    /**
     * Initialize the preconditioner object and calculate all data that is
     * necessary for applying it in a solver. This function is automatically
     * called when calling the constructor with the same arguments and is only
     * used if you create the preconditioner without arguments.
     */
    void
    initialize(const MatrixBase &    matrix,
               const AdditionalData &additional_data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular preconditioner.
     */
    AdditionalData additional_data;
  };



  /**
   * A class that implements the interface to use the PETSc SSOR
   * preconditioner.
   *
   * @note Only works in serial with a PETScWrappers::SparseMatrix.
   *
   * @ingroup PETScWrappers
   * @author Wolfgang Bangerth, Timo Heister, 2004, 2011
   */
  class PreconditionSSOR : public PreconditionerBase
  {
  public:
    /**
     * Standardized data struct to pipe additional flags to the
     * preconditioner.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default, set the damping parameter to one.
       */
      AdditionalData(const double omega = 1);

      /**
       * Relaxation parameter.
       */
      double omega;
    };

    /**
     * Empty Constructor. You need to call initialize() before using this
     * object.
     */
    PreconditionSSOR() = default;

    /**
     * Constructor. Take the matrix which is used to form the preconditioner,
     * and additional flags if there are any.
     */
    PreconditionSSOR(const MatrixBase &    matrix,
                     const AdditionalData &additional_data = AdditionalData());

    /**
     * Initialize the preconditioner object and calculate all data that is
     * necessary for applying it in a solver. This function is automatically
     * called when calling the constructor with the same arguments and is only
     * used if you create the preconditioner without arguments.
     */
    void
    initialize(const MatrixBase &    matrix,
               const AdditionalData &additional_data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular preconditioner.
     */
    AdditionalData additional_data;
  };



  /**
   * A class that implements the interface to use the PETSc Eisenstat
   * preconditioner, which implements SSOR on the diagonal block owned by
   * each processor.
   *
   * See the comment in the base class
   * @ref PreconditionerBase
   * for when this preconditioner may or may not work.
   *
   * @ingroup PETScWrappers
   * @author Wolfgang Bangerth, Timo Heister, 2004, 2011
   */
  class PreconditionEisenstat : public PreconditionerBase
  {
  public:
    /**
     * Standardized data struct to pipe additional flags to the
     * preconditioner.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default, set the damping parameter to one.
       */
      AdditionalData(const double omega = 1);

      /**
       * Relaxation parameter.
       */
      double omega;
    };

    /**
     * Empty Constructor. You need to call initialize() before using this
     * object.
     */
    PreconditionEisenstat() = default;

    /**
     * Constructor. Take the matrix which is used to form the preconditioner,
     * and additional flags if there are any.
     */
    PreconditionEisenstat(
      const MatrixBase &    matrix,
      const AdditionalData &additional_data = AdditionalData());

    /**
     * Initialize the preconditioner object and calculate all data that is
     * necessary for applying it in a solver. This function is automatically
     * called when calling the constructor with the same arguments and is only
     * used if you create the preconditioner without arguments.
     */
    void
    initialize(const MatrixBase &    matrix,
               const AdditionalData &additional_data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular preconditioner.
     */
    AdditionalData additional_data;
  };



  /**
   * A class that implements the interface to use the PETSc Incomplete
   * Cholesky preconditioner.
   *
   * @note Only works in serial with a PETScWrappers::SparseMatrix.
   *
   * @ingroup PETScWrappers
   * @author Wolfgang Bangerth, Timo Heister, 2004, 2011
   */
  class PreconditionICC : public PreconditionerBase
  {
  public:
    /**
     * Standardized data struct to pipe additional flags to the
     * preconditioner.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default, set the fill-in parameter to zero.
       */
      AdditionalData(const unsigned int levels = 0);

      /**
       * Fill-in parameter.
       */
      unsigned int levels;
    };

    /**
     * Empty Constructor. You need to call initialize() before using this
     * object.
     */
    PreconditionICC() = default;

    /**
     * Constructor. Take the matrix which is used to form the preconditioner,
     * and additional flags if there are any.
     */
    PreconditionICC(const MatrixBase &    matrix,
                    const AdditionalData &additional_data = AdditionalData());

    /**
     * Initialize the preconditioner object and calculate all data that is
     * necessary for applying it in a solver. This function is automatically
     * called when calling the constructor with the same arguments and is only
     * used if you create the preconditioner without arguments.
     */
    void
    initialize(const MatrixBase &    matrix,
               const AdditionalData &additional_data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular preconditioner.
     */
    AdditionalData additional_data;
  };



  /**
   * A class that implements the interface to use the PETSc ILU
   * preconditioner.
   *
   * @note Only works in serial with a PETScWrappers::SparseMatrix.
   *
   * @ingroup PETScWrappers
   * @author Wolfgang Bangerth, Timo Heister, 2004, 2011
   */
  class PreconditionILU : public PreconditionerBase
  {
  public:
    /**
     * Standardized data struct to pipe additional flags to the
     * preconditioner.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default, set the fill-in parameter to zero.
       */
      AdditionalData(const unsigned int levels = 0);

      /**
       * Fill-in parameter.
       */
      unsigned int levels;
    };

    /**
     * Empty Constructor. You need to call initialize() before using this
     * object.
     */
    PreconditionILU() = default;

    /**
     * Constructor. Take the matrix which is used to form the preconditioner,
     * and additional flags if there are any.
     */
    PreconditionILU(const MatrixBase &    matrix,
                    const AdditionalData &additional_data = AdditionalData());

    /**
     * Initialize the preconditioner object and calculate all data that is
     * necessary for applying it in a solver. This function is automatically
     * called when calling the constructor with the same arguments and is only
     * used if you create the preconditioner without arguments.
     */
    void
    initialize(const MatrixBase &    matrix,
               const AdditionalData &additional_data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular preconditioner.
     */
    AdditionalData additional_data;
  };



  /**
   * A class that implements the interface to use the PETSc LU preconditioner
   * (@p PCLU). Unlike classes like PreconditionILU, this class usually
   * (depending on the settings) performs an exact factorization of the
   * matrix, so it is not necessary to wrap it in an iterative solver. This
   * class is typically used with SolverPreOnly to get a direct
   * solver. Alternatively, you can use PreconditionerBase::vmult() directly.
   *
   * @note This is not a parallel preconditioner so it only works in serial
   * computations with a single processor.
   *
   * @ingroup PETScWrappers
   * @author Oliver Kayser-Herold, 2004
   */
  class PreconditionLU : public PreconditionerBase
  {
  public:
    /**
     * Standardized data struct to pipe additional flags to the
     * preconditioner.
     */
    struct AdditionalData
    {
      /**
       * Constructor. (Default values taken from function PCCreate_LU of the
       * PETSc lib.)
       */
      AdditionalData(const double pivoting   = 1.e-6,
                     const double zero_pivot = 1.e-12,
                     const double damping    = 0.0);

      /**
       * Determines, when Pivoting is done during LU decomposition. 0.0
       * indicates no pivoting, and 1.0 complete pivoting. Confer PETSc manual
       * for more details.
       */
      double pivoting;

      /**
       * Size at which smaller pivots are declared to be zero. Confer PETSc
       * manual for more details.
       */
      double zero_pivot;

      /**
       * This quantity is added to the diagonal of the matrix during
       * factorization.
       */
      double damping;
    };

    /**
     * Empty Constructor. You need to call initialize() before using this
     * object.
     */
    PreconditionLU() = default;

    /**
     * Constructor. Take the matrix which is used to form the preconditioner,
     * and additional flags if there are any.
     */
    PreconditionLU(const MatrixBase &    matrix,
                   const AdditionalData &additional_data = AdditionalData());

    /**
     * Initialize the preconditioner object and calculate all data that is
     * necessary for applying it in a solver. This function is automatically
     * called when calling the constructor with the same arguments and is only
     * used if you create the preconditioner without arguments.
     */
    void
    initialize(const MatrixBase &    matrix,
               const AdditionalData &additional_data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular preconditioner.
     */
    AdditionalData additional_data;
  };



  /**
   * A class that implements the interface to use the BoomerAMG algebraic
   * multigrid preconditioner from the HYPRE suite. Note that PETSc has to be
   * configured with HYPRE (e.g. with \--download-hypre=1).
   *
   * The preconditioner does support parallel distributed computations. See
   * step-40 for an example.
   *
   * @ingroup PETScWrappers
   * @author Timo Heister, 2010
   */
  class PreconditionBoomerAMG : public PreconditionerBase
  {
  public:
    /**
     * Standardized data struct to pipe additional flags to the
     * preconditioner.
     */
    struct AdditionalData
    {
      /**
       * Constructor. Note that BoomerAMG offers a lot more options to set
       * than what is exposed here.
       */
      AdditionalData(const bool         symmetric_operator = false,
                     const double       strong_threshold   = 0.25,
                     const double       max_row_sum        = 0.9,
                     const unsigned int aggressive_coarsening_num_levels = 0,
                     const bool         output_details = false);

      /**
       * Set this flag to true if you have a symmetric system matrix and you
       * want to use a solver which assumes a symmetric preconditioner like
       * CG. The relaxation is done with SSOR/Jacobi when set to true and with
       * SOR/Jacobi otherwise.
       */
      bool symmetric_operator;

      /**
       * Threshold of when nodes are considered strongly connected. See
       * HYPRE_BoomerAMGSetStrongThreshold(). Recommended values are 0.25 for
       * 2d and 0.5 for 3d problems, but it is problem dependent.
       */
      double strong_threshold;

      /**
       * If set to a value smaller than 1.0 then diagonally dominant parts of
       * the matrix are treated as having no strongly connected nodes. If the
       * row sum weighted by the diagonal entry is bigger than the given
       * value, it is considered diagonally dominant. This feature is turned
       * of by setting the value to 1.0. This is the default as some matrices
       * can result in having only diagonally dominant entries and thus no
       * multigrid levels are constructed. The default in BoomerAMG for this
       * is 0.9. When you try this, check for a reasonable number of levels
       * created.
       */
      double max_row_sum;

      /**
       * Number of levels of aggressive coarsening. Increasing this value
       * reduces the construction time and memory requirements but may
       * decrease effectiveness.
       */
      unsigned int aggressive_coarsening_num_levels;

      /**
       * Setting this flag to true produces debug output from HYPRE, when the
       * preconditioner is constructed.
       */
      bool output_details;
    };

    /**
     * Empty Constructor. You need to call initialize() before using this
     * object.
     */
    PreconditionBoomerAMG() = default;

    /**
     * Constructor. Take the matrix which is used to form the preconditioner,
     * and additional flags if there are any.
     */
    PreconditionBoomerAMG(
      const MatrixBase &    matrix,
      const AdditionalData &additional_data = AdditionalData());

    /**
     * Same as above but without setting a matrix to form the preconditioner.
     * Intended to be used with SLEPc objects.
     */
    PreconditionBoomerAMG(
      const MPI_Comm        communicator,
      const AdditionalData &additional_data = AdditionalData());


    /**
     * Initialize the preconditioner object and calculate all data that is
     * necessary for applying it in a solver. This function is automatically
     * called when calling the constructor with the same arguments and is only
     * used if you create the preconditioner without arguments.
     */
    void
    initialize(const MatrixBase &    matrix,
               const AdditionalData &additional_data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular preconditioner.
     */
    AdditionalData additional_data;

    /**
     * Initialize the preconditioner object without knowing a particular
     * matrix. This function sets up appropriate parameters to the underlying
     * PETSc object after it has been created.
     */
    void
    initialize();
  };



  /**
   * A class that implements the interface to use the ParaSails sparse
   * approximate inverse preconditioner from the HYPRE suite. Note that PETSc
   * has to be configured with HYPRE (e.g. with \--download-hypre=1).
   *
   * ParaSails uses least-squares minimization to compute a sparse approximate
   * inverse. The sparsity pattern used is the pattern of a power of a
   * sparsified matrix. ParaSails also uses a post-filtering technique to
   * reduce the cost of applying the preconditioner.
   *
   * ParaSails solves symmetric positive definite (SPD) problems using a
   * factorized SPD preconditioner and can also solve general (nonsymmetric
   * and/or indefinite) problems with a nonfactorized preconditioner. The
   * problem type has to be set in @p AdditionalData.
   *
   * The preconditioner does support parallel distributed computations.
   *
   * @ingroup PETScWrappers
   * @author Martin Steigemann, 2012
   */
  class PreconditionParaSails : public PreconditionerBase
  {
  public:
    /**
     * Standardized data struct to pipe additional flags to the
     * preconditioner.
     */
    struct AdditionalData
    {
      /**
       * Constructor.
       */
      AdditionalData(const unsigned int symmetric      = 1,
                     const unsigned int n_levels       = 1,
                     const double       threshold      = 0.1,
                     const double       filter         = 0.05,
                     const bool         output_details = false);

      /**
       * This parameter specifies the type of problem to solve:
       * <ul>
       * <li> @p 0: nonsymmetric and/or indefinite problem, and nonsymmetric
       * preconditioner
       * <li> @p 1: SPD problem, and SPD (factored) preconditioner
       * <li> @p 2: nonsymmetric, definite problem, and SPD (factored)
       * preconditioner
       * </ul>
       * Default is <tt>symmetric = 1</tt>.
       */
      unsigned int symmetric;

      /**
       * The sparsity pattern used for the approximate inverse is the pattern
       * of a power <tt>B^m</tt> where <tt>B</tt> has been sparsified from the
       * given matrix <tt>A</tt>, <tt>n_level</tt> is equal to <tt>m+1</tt>.
       * Default value is <tt>n_levels = 1</tt>.
       */
      unsigned int n_levels;

      /**
       * Sparsification is performed by dropping nonzeros which are smaller
       * than <tt>thresh</tt> in magnitude. Lower values of <tt>thresh</tt>
       * lead to more accurate, but also more expensive preconditioners.
       * Default value is <tt>thresh = 0.1</tt>. Setting <tt>thresh < 0</tt> a
       * threshold is selected automatically, such that <tt>-thresh</tt>
       * represents the fraction of nonzero elements that are dropped. For
       * example, if <tt>thresh = -0.9</tt>, then <tt>B</tt> will contain
       * about ten percent of the nonzeros of the given matrix <tt>A</tt>.
       */
      double threshold;

      /**
       * Filtering is a post-processing procedure, <tt>filter</tt> represents
       * a fraction of nonzero elements that are dropped after creating the
       * approximate inverse sparsity pattern. Default value is <tt>filter =
       * 0.05</tt>. Setting <tt>filter < 0</tt> a value is selected
       * automatically, such that <tt>-filter</tt> represents the fraction of
       * nonzero elements that are dropped. For example, if <tt>thresh =
       * -0.9</tt>, then about 90 percent of the entries in the computed
       * approximate inverse are dropped.
       */
      double filter;

      /**
       * Setting this flag to true produces output from HYPRE, when the
       * preconditioner is constructed.
       */
      bool output_details;
    };



    /**
     * Empty Constructor. You need to call initialize() before using this
     * object.
     */
    PreconditionParaSails() = default;

    /**
     * Constructor. Take the matrix which is used to form the preconditioner,
     * and additional flags if there are any.
     */
    PreconditionParaSails(
      const MatrixBase &    matrix,
      const AdditionalData &additional_data = AdditionalData());

    /**
     * Initialize the preconditioner object and calculate all data that is
     * necessary for applying it in a solver. This function is automatically
     * called when calling the constructor with the same arguments and is only
     * used if you create the preconditioner without arguments.
     */
    void
    initialize(const MatrixBase &    matrix,
               const AdditionalData &additional_data = AdditionalData());

  private:
    /**
     * Store a copy of the flags for this particular preconditioner.
     */
    AdditionalData additional_data;
  };



  /**
   * A class that implements a non-preconditioned method.
   *
   * @ingroup PETScWrappers
   * @author Martin Steigemann, 2012
   */
  class PreconditionNone : public PreconditionerBase
  {
  public:
    /**
     * Standardized data struct to pipe additional flags to the
     * preconditioner.
     */
    struct AdditionalData
    {};

    /**
     * Empty Constructor. You need to call initialize() before using this
     * object.
     */
    PreconditionNone() = default;

    /**
     * Constructor. Take the matrix which is used to form the preconditioner,
     * and additional flags if there are any. The matrix is completely ignored
     * in computations.
     */
    PreconditionNone(const MatrixBase &    matrix,
                     const AdditionalData &additional_data = AdditionalData());

    /**
     * Initialize the preconditioner object and calculate all data that is
     * necessary for applying it in a solver. This function is automatically
     * called when calling the constructor with the same arguments and is only
     * used if you create the preconditioner without arguments. The matrix is
     * completely ignored in computations.
     */
    void
    initialize(const MatrixBase &    matrix,
               const AdditionalData &additional_data = AdditionalData());

  private:
    /**
     * Store a copy of the flags for this particular preconditioner.
     */
    AdditionalData additional_data;
  };
} // namespace PETScWrappers



DEAL_II_NAMESPACE_CLOSE


#  endif // DEAL_II_WITH_PETSC

#endif
/*--------------------------- petsc_precondition.h --------------------------*/
