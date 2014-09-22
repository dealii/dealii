// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2014 by the deal.II authors
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

#ifndef __deal2__trilinos_precondition_h
#define __deal2__trilinos_precondition_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_TRILINOS

#  include <deal.II/base/subscriptor.h>
#  include <deal.II/base/std_cxx11/shared_ptr.h>

#  include <deal.II/lac/trilinos_vector_base.h>
#  include <deal.II/lac/parallel_vector.h>

#  ifdef DEAL_II_WITH_MPI
#    include <Epetra_MpiComm.h>
#  else
#    include <Epetra_SerialComm.h>
#  endif
#  include <Epetra_Map.h>

#  include <Teuchos_ParameterList.hpp>
#  include <Epetra_RowMatrix.h>
#  include <Epetra_Vector.h>

// forward declarations
class Ifpack_Preconditioner;
class Ifpack_Chebyshev;
namespace ML_Epetra
{
  class MultiLevelPreconditioner;
}


DEAL_II_NAMESPACE_OPEN

// forward declarations
template <typename number> class SparseMatrix;
template <typename number> class Vector;
class SparsityPattern;

/*! @addtogroup TrilinosWrappers
 *@{
 */

namespace TrilinosWrappers
{
  // forward declarations
  class SparseMatrix;
  class BlockSparseMatrix;
  class SolverBase;

  /**
   * The base class for all preconditioners based on Trilinos sparse
   * matrices.
   *
   * @ingroup TrilinosWrappers
   * @ingroup Preconditioners
   * @author Martin Kronbichler, 2008
   */
  class PreconditionBase : public Subscriptor
  {
  public:
    /**
     * Declare the type for container size.
     */
    typedef dealii::types::global_dof_index size_type;

    /**
     * Standardized data struct to pipe additional flags to the
     * preconditioner.
     */
    struct AdditionalData
    {};

    /**
     * Constructor. Does not do anything. The <tt>initialize</tt> function of
     * the derived classes will have to create the preconditioner from a given
     * sparse matrix.
     */
    PreconditionBase ();

    /**
     * Copy constructor.
     */
    PreconditionBase (const PreconditionBase &);

    /**
     * Destructor.
     */
    ~PreconditionBase ();

    /**
     * Destroys the preconditioner, leaving an object like just after having
     * called the constructor.
     */
    void clear ();

    /**
     * Apply the preconditioner.
     */
    virtual void vmult (VectorBase       &dst,
                        const VectorBase &src) const;

    /**
     * Apply the transpose preconditioner.
     */
    virtual void Tvmult (VectorBase       &dst,
                         const VectorBase &src) const;

    /**
     * Apply the preconditioner on deal.II data structures instead of the ones
     * provided in the Trilinos wrapper class.
     */
    virtual void vmult (dealii::Vector<double>       &dst,
                        const dealii::Vector<double> &src) const;

    /**
     * Apply the transpose preconditioner on deal.II data structures instead
     * of the ones provided in the Trilinos wrapper class.
     */
    virtual void Tvmult (dealii::Vector<double>       &dst,
                         const dealii::Vector<double> &src) const;

    /**
     * Apply the preconditioner on deal.II parallel data structures instead of
     * the ones provided in the Trilinos wrapper class.
     */
    virtual void vmult (dealii::parallel::distributed::Vector<double>       &dst,
                        const dealii::parallel::distributed::Vector<double> &src) const;

    /**
     * Apply the transpose preconditioner on deal.II parallel data structures
     * instead of the ones provided in the Trilinos wrapper class.
     */
    virtual void Tvmult (dealii::parallel::distributed::Vector<double>       &dst,
                         const dealii::parallel::distributed::Vector<double> &src) const;

    /**
     * Exception.
     */
    DeclException1 (ExcNonMatchingMaps,
                    std::string,
                    << "The sparse matrix the preconditioner is based on "
                    << "uses a map that is not compatible to the one in vector "
                    << arg1
                    << ". Check preconditioner and matrix setup.");

    friend class SolverBase;
    friend class PreconditionStokes;

  protected:
    /**
     * This is a pointer to the preconditioner object that is used when
     * applying the preconditioner.
     */
    std_cxx11::shared_ptr<Epetra_Operator> preconditioner;

    /**
     * Internal communication pattern in case the matrix needs to be copied
     * from deal.II format.
     */
#ifdef DEAL_II_WITH_MPI
    Epetra_MpiComm     communicator;
#else
    Epetra_SerialComm  communicator;
#endif

    /**
     * Internal Trilinos map in case the matrix needs to be copied from
     * deal.II format.
     */
    std_cxx11::shared_ptr<Epetra_Map>   vector_distributor;
  };


  /**
   * A wrapper class for a (pointwise) Jacobi preconditioner for
   * Trilinos matrices. This preconditioner works both in serial and in
   * parallel, depending on the matrix it is based on.
   *
   * The AdditionalData data structure allows to set preconditioner
   * options. For the Jacobi preconditioner, these options are the
   * damping parameter <tt>omega</tt> and a <tt>min_diagonal</tt>
   * argument that can be used to make the preconditioner work even if
   * the matrix contains some zero elements on the diagonal. The default
   * settings are 1 for the damping parameter and zero for the diagonal
   * augmentation.
   *
   * @ingroup TrilinosWrappers
   * @ingroup Preconditioners
   * @author Martin Kronbichler, 2008
   */
  class PreconditionJacobi : public PreconditionBase
  {
  public:

    /**
     * Standardized data struct to pipe additional flags to the
     * preconditioner. The parameter <tt>omega</tt> specifies the relaxation
     * parameter in the Jacobi preconditioner. The parameter
     * <tt>min_diagonal</tt> can be used to make the application of the
     * preconditioner also possible when some diagonal elements are zero. In a
     * default application this would mean that we divide by zero, so by
     * setting the parameter <tt>min_diagonal</tt> to a small nonzero value
     * the SOR will work on a matrix that is not too far away from the one we
     * want to treat.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default, set the damping parameter to one, and do not
       * modify the diagonal.
       */
      AdditionalData (const double       omega = 1,
                      const double       min_diagonal = 0,
                      const unsigned int n_sweeps = 1);

      /**
       * This specifies the relaxation parameter in the Jacobi preconditioner.
       */
      double omega;

      /**
       * This specifies the minimum value the diagonal elements should
       * have. This might be necessary when the Jacobi preconditioner is used
       * on matrices with zero diagonal elements. In that case, a
       * straight-forward application would not be possible since we would
       * divide by zero.
       */
      double min_diagonal;

      /**
       * Sets how many times the given operation should be applied during the
       * vmult() operation.
       */
      unsigned int n_sweeps;
    };

    /**
     * Take the sparse matrix the preconditioner object should be built of,
     * and additional flags (damping parameter, etc.) if there are any.
     */
    void initialize (const SparseMatrix   &matrix,
                     const AdditionalData &additional_data = AdditionalData());
  };




  /**
   * A wrapper class for a (pointwise) SSOR preconditioner for Trilinos
   * matrices. This preconditioner works both in serial and in parallel,
   * depending on the matrix it is based on.
   *
   * The AdditionalData data structure allows to set preconditioner
   * options. For the SSOR preconditioner, these options are the
   * damping/relaxation parameter <tt>omega</tt>, a
   * <tt>min_diagonal</tt> argument that can be used to make the
   * preconditioner work even if the matrix contains some zero elements
   * on the diagonal, and a parameter <tt>overlap</tt> that determines
   * if and how much overlap there should be between the matrix
   * partitions on the various MPI processes. The default settings are 1
   * for the relaxation parameter, 0 for the diagonal augmentation and 0
   * for the overlap.
   *
   * Note that a parallel application of the SSOR preconditioner is
   * actually a block-Jacobi preconditioner with block size equal to the
   * local matrix size. Spoken more technically, this parallel operation
   * is an <a
   * href="http://en.wikipedia.org/wiki/Additive_Schwarz_method">additive
   * Schwarz method</a> with an SSOR <em>approximate solve</em> as inner
   * solver, based on the outer parallel partitioning.
   *
   * @ingroup TrilinosWrappers
   * @ingroup Preconditioners
   * @author Wolfgang Bangerth, 2008
   */
  class PreconditionSSOR : public PreconditionBase
  {
  public:

    /**
     * Standardized data struct to pipe additional flags to the
     * preconditioner. The parameter <tt>omega</tt> specifies the relaxation
     * parameter in the SSOR preconditioner. The parameter
     * <tt>min_diagonal</tt> can be used to make the application of the
     * preconditioner also possible when some diagonal elements are zero. In a
     * default application this would mean that we divide by zero, so by
     * setting the parameter <tt>min_diagonal</tt> to a small nonzero value
     * the SOR will work on a matrix that is not too far away from the one we
     * want to treat. Finally, <tt>overlap</tt> governs the overlap of the
     * partitions when the preconditioner runs in parallel, forming a
     * so-called additive Schwarz preconditioner.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default, set the damping parameter to one, we do not
       * modify the diagonal, and there is no overlap (i.e. in parallel, we
       * run a BlockJacobi preconditioner, where each block is inverted
       * approximately by an SSOR).
       */
      AdditionalData (const double       omega = 1,
                      const double       min_diagonal = 0,
                      const unsigned int overlap = 0,
                      const unsigned int n_sweeps = 1);

      /**
       * This specifies the (over-) relaxation parameter in the SSOR
       * preconditioner.
       */
      double omega;

      /**
       * This specifies the minimum value the diagonal elements should
       * have. This might be necessary when the SSOR preconditioner is used on
       * matrices with zero diagonal elements. In that case, a
       * straight-forward application would not be possible since we divide by
       * the diagonal element.
       */
      double min_diagonal;

      /**
       * This determines how large the overlap of the local matrix portions on
       * each processor in a parallel application should be.
       */
      unsigned int overlap;

      /**
       * Sets how many times the given operation should be applied during the
       * vmult() operation.
       */
      unsigned int n_sweeps;
    };

    /**
     * Take the sparse matrix the preconditioner object should be built of,
     * and additional flags (damping parameter, overlap in parallel
     * computations, etc.) if there are any.
     */
    void initialize (const SparseMatrix   &matrix,
                     const AdditionalData &additional_data = AdditionalData());
  };




  /**
   * A wrapper class for a (pointwise) SOR preconditioner for Trilinos
   * matrices. This preconditioner works both in serial and in parallel,
   * depending on the matrix it is based on.
   *
   * The AdditionalData data structure allows to set preconditioner
   * options. For the SOR preconditioner, these options are the
   * damping/relaxation parameter <tt>omega</tt>, a
   * <tt>min_diagonal</tt> argument that can be used to make the
   * preconditioner work even if the matrix contains some zero elements
   * on the diagonal, and a parameter <tt>overlap</tt> that determines
   * if and how much overlap there should be between the matrix
   * partitions on the various MPI processes. The default settings are 1
   * for the relaxation parameter, 0 for the diagonal augmentation and 0
   * for the overlap.
   *
   * Note that a parallel application of the SOR preconditioner is
   * actually a block-Jacobi preconditioner with block size equal to the
   * local matrix size. Spoken more technically, this parallel operation
   * is an <a
   * href="http://en.wikipedia.org/wiki/Additive_Schwarz_method">additive
   * Schwarz method</a> with an SOR <em>approximate solve</em> as inner
   * solver, based on the outer parallel partitioning.
   *
   * @ingroup TrilinosWrappers
   * @ingroup Preconditioners
   * @author Martin Kronbichler, 2008
   */
  class PreconditionSOR : public PreconditionBase
  {
  public:

    /**
     * Standardized data struct to pipe additional flags to the
     * preconditioner. The parameter <tt>omega</tt> specifies the relaxation
     * parameter in the SOR preconditioner. The parameter
     * <tt>min_diagonal</tt> can be used to make the application of the
     * preconditioner also possible when some diagonal elements are zero. In a
     * default application this would mean that we divide by zero, so by
     * setting the parameter <tt>min_diagonal</tt> to a small nonzero value
     * the SOR will work on a matrix that is not too far away from the one we
     * want to treat. Finally, <tt>overlap</tt> governs the overlap of the
     * partitions when the preconditioner runs in parallel, forming a
     * so-called additive Schwarz preconditioner.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default, set the damping parameter to one, we do not
       * modify the diagonal, and there is no overlap (i.e. in parallel, we
       * run a BlockJacobi preconditioner, where each block is inverted
       * approximately by an SOR.
       */
      AdditionalData (const double       omega = 1,
                      const double       min_diagonal = 0,
                      const unsigned int overlap = 0,
                      const unsigned int n_sweeps = 1);

      /**
       * This specifies the (over-) relaxation parameter in the SOR
       * preconditioner.
       */
      double omega;

      /**
       * This specifies the minimum value the diagonal elements should
       * have. This might be necessary when the SOR preconditioner is used on
       * matrices with zero diagonal elements. In that case, a
       * straight-forward application would not be possible since we divide by
       * the diagonal element.
       */
      double min_diagonal;

      /**
       * This determines how large the overlap of the local matrix portions on
       * each processor in a parallel application should be.
       */
      unsigned int overlap;

      /**
       * Sets how many times the given operation should be applied during the
       * vmult() operation.
       */
      unsigned int n_sweeps;
    };

    /**
     * Take the sparse matrix the preconditioner object should be built of,
     * and additional flags (damping parameter, overlap in parallel
     * computations etc.) if there are any.
     */
    void initialize (const SparseMatrix   &matrix,
                     const AdditionalData &additional_data = AdditionalData());
  };



  /**
   * A wrapper class for a block Jacobi preconditioner for Trilinos
   * matrices. As opposed to PreconditionSOR where each row is treated
   * separately, this scheme collects block of a given size and inverts a full
   * matrix for all these rows simultaneously. Trilinos allows to select
   * several strategies for selecting which rows form a block, including
   * "linear" (i.e., divide the local range of the matrix in slices of the
   * block size), "greedy" or "metis". Note that the term <em>block
   * Jacobi</em> does not relate to possible blocks in the MPI setting, but
   * small blocks of dense matrices extracted from the sparse matrix local to
   * each processor.
   *
   * The AdditionalData data structure allows to set preconditioner
   * options.
   *
   * @ingroup TrilinosWrappers
   * @ingroup Preconditioners
   * @author Martin Kronbichler, 2014
   */
  class PreconditionBlockJacobi : public PreconditionBase
  {
  public:

    /**
     * Standardized data struct to pipe additional flags to the
     * preconditioner. The parameter <tt>block_size</tt> sets the size of
     * small blocks. It is recommended to choose this parameter not too large
     * (a few hundreds at most) since this implementation uses a dense matrix
     * for the block. The parameter <tt>block_creation_type</tt> allows to
     * pass the strategy for finding the blocks to Ifpack. The parameter
     * <tt>omega</tt> specifies the relaxation parameter in the SOR
     * preconditioner. The parameter <tt>min_diagonal</tt> can be used to make
     * the application of the preconditioner also possible when some diagonal
     * elements are zero. In a default application this would mean that we
     * divide by zero, so by setting the parameter <tt>min_diagonal</tt> to a
     * small nonzero value the SOR will work on a matrix that is not too far
     * away from the one we want to treat. Finally, <tt>overlap</tt> governs
     * the overlap of the partitions when the preconditioner runs in parallel,
     * forming a so-called additive Schwarz preconditioner.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default, use a block size of 1, use linear
       * subdivision of the rows, set the damping parameter to one, and do not
       * modify the diagonal.
       */
      AdditionalData (const unsigned int block_size = 1,
                      const std::string  block_creation_type = "linear",
                      const double       omega = 1,
                      const double       min_diagonal = 0,
                      const unsigned int n_sweeps = 1);

      /**
       * This specifies the size of blocks.
       */
      unsigned int block_size;

      /**
       * Strategy for creation of blocks passed on to Ifpack block relaxation
       * (variable 'partitioner: type') with this string as the given
       * value. Available types in Ifpack include "linear" (i.e., divide the
       * local range of the matrix in slices of the block size), "greedy"
       * "metis". For a full list, see the documentation of Ifpack.
       */
      std::string block_creation_type;

      /**
       * This specifies the (over-) relaxation parameter in the Jacobi
       * preconditioner.
       */
      double omega;

      /**
       * This specifies the minimum value the diagonal elements should
       * have. This might be necessary when the block Jacobi preconditioner is
       * used on matrices with zero diagonal elements. In that case, a
       * straight-forward application would not be possible since we divide by
       * the diagonal element.
       */
      double min_diagonal;

      /**
       * Sets how many times the given operation should be applied during the
       * vmult() operation.
       */
      unsigned int n_sweeps;
    };

    /**
     * Take the sparse matrix the preconditioner object should be built of,
     * and additional flags (damping parameter, etc.) if there are any.
     */
    void initialize (const SparseMatrix   &matrix,
                     const AdditionalData &additional_data = AdditionalData());
  };




  /**
   * A wrapper class for a block SSOR preconditioner for Trilinos matrices. As
   * opposed to PreconditionSSOR where each row is treated separately
   * (point-wise), this scheme collects block of a given size and inverts a
   * full matrix for all these rows simultaneously. Trilinos allows to select
   * several strategies for selecting which rows form a block, including
   * "linear" (i.e., divide the local range of the matrix in slices of the
   * block size), "greedy" or "metis".
   *
   * The AdditionalData data structure allows to set preconditioner
   * options.
   *
   * Note that a parallel application of this preconditioner is actually a
   * block-Jacobi preconditioner with (outer) block size equal to the local
   * matrix size. Spoken more technically, this parallel operation is an <a
   * href="http://en.wikipedia.org/wiki/Additive_Schwarz_method">additive
   * Schwarz method</a> with a block SSOR <em>approximate solve</em> as inner
   * solver, based on the outer parallel partitioning.
   *
   * @ingroup TrilinosWrappers
   * @ingroup Preconditioners
   * @author Martin Kronbichler, 2014
   */
  class PreconditionBlockSSOR : public PreconditionBase
  {
  public:

    /**
     * Standardized data struct to pipe additional flags to the
     * preconditioner. The parameter <tt>block_size</tt> sets the size of
     * small blocks. It is recommended to choose this parameter not too large
     * (a few hundreds at most) since this implementation uses a dense matrix
     * for the block. The parameter <tt>block_creation_type</tt> allows to
     * pass the strategy for finding the blocks to Ifpack. The parameter
     * <tt>omega</tt> specifies the relaxation parameter in the SSOR
     * preconditioner. The parameter <tt>min_diagonal</tt> can be used to make
     * the application of the preconditioner also possible when some diagonal
     * elements are zero. In a default application this would mean that we
     * divide by zero, so by setting the parameter <tt>min_diagonal</tt> to a
     * small nonzero value the SOR will work on a matrix that is not too far
     * away from the one we want to treat. Finally, <tt>overlap</tt> governs
     * the overlap of the partitions when the preconditioner runs in parallel,
     * forming a so-called additive Schwarz preconditioner.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default, use a block size of 1, use linear
       * subdivision of the rows, set the damping parameter to one, we do not
       * modify the diagonal, and there is no overlap (i.e. in parallel, we
       * run a BlockJacobi preconditioner, where each block is inverted
       * approximately by a block SOR).
       */
      AdditionalData (const unsigned int block_size = 1,
                      const std::string  block_creation_type = "linear",
                      const double       omega = 1,
                      const double       min_diagonal = 0,
                      const unsigned int overlap = 0,
                      const unsigned int n_sweeps = 1);

      /**
       * This specifies the size of blocks.
       */
      unsigned int block_size;

      /**
       * Strategy for creation of blocks passed on to Ifpack block relaxation
       * (variable 'partitioner: type') with this string as the given
       * value. Available types in Ifpack include "linear" (i.e., divide the
       * local range of the matrix in slices of the block size), "greedy"
       * "metis". For a full list, see the documentation of Ifpack.
       */
      std::string block_creation_type;

      /**
       * This specifies the (over-) relaxation parameter in the SOR
       * preconditioner.
       */
      double omega;

      /**
       * This specifies the minimum value the diagonal elements should
       * have. This might be necessary when the SSOR preconditioner is used on
       * matrices with zero diagonal elements. In that case, a
       * straight-forward application would not be possible since we divide by
       * the diagonal element.
       */
      double min_diagonal;

      /**
       * This determines how large the overlap of the local matrix portions on
       * each processor in a parallel application should be.
       */
      unsigned int overlap;

      /**
       * Sets how many times the given operation should be applied during the
       * vmult() operation.
       */
      unsigned int n_sweeps;
    };

    /**
     * Take the sparse matrix the preconditioner object should be built of,
     * and additional flags (damping parameter, overlap in parallel
     * computations, etc.) if there are any.
     */
    void initialize (const SparseMatrix   &matrix,
                     const AdditionalData &additional_data = AdditionalData());
  };




  /**
   * A wrapper class for a block SOR preconditioner for Trilinos matrices. As
   * opposed to PreconditionSOR where each row is treated separately, this
   * scheme collects block of a given size and inverts a full matrix for all
   * these rows simultaneously. Trilinos allows to select several strategies
   * for selecting which rows form a block, including "linear" (i.e., divide
   * the local range of the matrix in slices of the block size), "greedy" or
   * "metis".
   *
   * The AdditionalData data structure allows to set preconditioner
   * options.
   *
   * Note that a parallel application of this preconditioner is actually a
   * block-Jacobi preconditioner with (outer) block size equal to the local
   * matrix size. Spoken more technically, this parallel operation is an <a
   * href="http://en.wikipedia.org/wiki/Additive_Schwarz_method">additive
   * Schwarz method</a> with a block SOR <em>approximate solve</em> as inner
   * solver, based on the outer parallel partitioning.
   *
   * @ingroup TrilinosWrappers
   * @ingroup Preconditioners
   * @author Martin Kronbichler, 2014
   */
  class PreconditionBlockSOR : public PreconditionBase
  {
  public:

    /**
     * Standardized data struct to pipe additional flags to the
     * preconditioner. The parameter <tt>block_size</tt> sets the size of
     * small blocks. It is recommended to choose this parameter not too large
     * (a few hundreds at most) since this implementation uses a dense matrix
     * for the block. The parameter <tt>block_creation_type</tt> allows to
     * pass the strategy for finding the blocks to Ifpack. The parameter
     * <tt>omega</tt> specifies the relaxation parameter in the SOR
     * preconditioner. The parameter <tt>min_diagonal</tt> can be used to make
     * the application of the preconditioner also possible when some diagonal
     * elements are zero. In a default application this would mean that we
     * divide by zero, so by setting the parameter <tt>min_diagonal</tt> to a
     * small nonzero value the SOR will work on a matrix that is not too far
     * away from the one we want to treat. Finally, <tt>overlap</tt> governs
     * the overlap of the partitions when the preconditioner runs in parallel,
     * forming a so-called additive Schwarz preconditioner.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default, use a block size of 1, use linear
       * subdivision of the rows, set the damping parameter to one, we do not
       * modify the diagonal, and there is no overlap (i.e. in parallel, we
       * run a BlockJacobi preconditioner, where each block is inverted
       * approximately by a block SOR).
       */
      AdditionalData (const unsigned int block_size = 1,
                      const std::string  block_creation_type = "linear",
                      const double       omega = 1,
                      const double       min_diagonal = 0,
                      const unsigned int overlap = 0,
                      const unsigned int n_sweeps = 1);

      /**
       * This specifies the size of blocks.
       */
      unsigned int block_size;

      /**
       * Strategy for creation of blocks passed on to Ifpack block relaxation
       * (variable 'partitioner: type') with this string as the given
       * value. Available types in Ifpack include "linear" (i.e., divide the
       * local range of the matrix in slices of the block size), "greedy"
       * "metis". For a full list, see the documentation of Ifpack.
       */
      std::string block_creation_type;

      /**
       * This specifies the (over-) relaxation parameter in the SOR
       * preconditioner.
       */
      double omega;

      /**
       * This specifies the minimum value the diagonal elements should
       * have. This might be necessary when the SOR preconditioner is used on
       * matrices with zero diagonal elements. In that case, a
       * straight-forward application would not be possible since we divide by
       * the diagonal element.
       */
      double min_diagonal;

      /**
       * This determines how large the overlap of the local matrix portions on
       * each processor in a parallel application should be.
       */
      unsigned int overlap;

      /**
       * Sets how many times the given operation should be applied during the
       * vmult() operation.
       */
      unsigned int n_sweeps;
    };

    /**
     * Take the sparse matrix the preconditioner object should be built of,
     * and additional flags (damping parameter, overlap in parallel
     * computations etc.) if there are any.
     */
    void initialize (const SparseMatrix   &matrix,
                     const AdditionalData &additional_data = AdditionalData());
  };



  /**
   * A wrapper class for an incomplete Cholesky factorization (IC)
   * preconditioner for @em symmetric Trilinos matrices. This
   * preconditioner works both in serial and in parallel, depending on
   * the matrix it is based on. In general, an incomplete factorization
   * does not take all fill-in elements that would appear in a full
   * factorization (that is the basis for a direct solve). Trilinos
   * allows to set the amount of fill-in elements, governed by the
   * additional data argument <tt>ic_fill</tt>, so one can gradually
   * choose between a factorization on the sparse matrix structure only
   * (<tt>ic_fill=0</tt>) to a full factorization (<tt>ic_fill</tt> in
   * the range of 10 to 50, depending on the spatial dimension of the
   * PDE problem and the degree of the finite element basis functions;
   * generally, more required fill-in elements require this parameter to
   * be set to a higher integer value).
   *
   * The AdditionalData data structure allows to set preconditioner
   * options. Besides the fill-in argument, these options are some
   * options for perturbations (see the documentation of the
   * AdditionalData structure for details), and a parameter
   * <tt>overlap</tt> that determines if and how much overlap there
   * should be between the matrix partitions on the various MPI
   * processes.  The default settings are 0 for the additional fill-in, 0
   * for the absolute augmentation tolerance, 1 for the relative
   * augmentation tolerance, 0 for the overlap.
   *
   * Note that a parallel application of the IC preconditioner is
   * actually a block-Jacobi preconditioner with block size equal to the
   * local matrix size. Spoken more technically, this parallel operation
   * is an <a
   * href="http://en.wikipedia.org/wiki/Additive_Schwarz_method">additive
   * Schwarz method</a> with an IC <em>approximate solve</em> as inner
   * solver, based on the (outer) parallel partitioning.
   *
   * @ingroup TrilinosWrappers
   * @ingroup Preconditioners
   * @author Martin Kronbichler, 2008
   */
  class PreconditionIC : public PreconditionBase
  {
  public:
    /**
     * Standardized data struct to pipe additional parameters to the
     * preconditioner. The Trilinos IC decomposition allows for some fill-in,
     * so it actually is a threshold incomplete Cholesky factorization. The
     * amount of fill-in, and hence, the amount of memory used by this
     * preconditioner, is controlled by the parameter <tt>ic_fill</tt>, which
     * specifies this as a double. When forming the preconditioner, for
     * certain problems bad conditioning (or just bad luck) can cause the
     * preconditioner to be very poorly conditioned. Hence it can help to add
     * diagonal perturbations to the original matrix and form the
     * preconditioner for this slightly better matrix. <tt>ic_atol</tt> is an
     * absolute perturbation that is added to the diagonal before forming the
     * prec, and <tt>ic_rtol</tt> is a scaling factor $rtol \geq 1$. The last
     * parameter specifies the overlap of the partitions when the
     * preconditioner runs in parallel.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default, set the drop tolerance to 0, the level of
       * extra fill-ins is set to be zero (just use the matrix structure, do
       * not generate any additional fill-in), the tolerance level are 0 and
       * 1, respectively, and the overlap in case of a parallel execution is
       * zero. This overlap in a block-application of the IC in the parallel
       * case makes the preconditioner a so-called additive Schwarz
       * preconditioner.
       */
      AdditionalData (const unsigned int ic_fill = 0,
                      const double       ic_atol = 0.,
                      const double       ic_rtol = 1.,
                      const unsigned int overlap = 0);

      /**
       * This specifies the amount of additional fill-in elements besides the
       * sparse matrix structure. When <tt>ic_fill</tt> is large, this means
       * that many fill-ins will be added, so that the IC preconditioner comes
       * closer to a direct sparse Cholesky decomposition. Note, however, that
       * this will drastically increase the memory requirement, especially
       * when the preconditioner is used in 3D.
       */
      unsigned int ic_fill;

      /**
       * This specifies the amount of an absolute perturbation that will be
       * added to the diagonal of the matrix, which sometimes can help to get
       * better preconditioners.
       */
      double ic_atol;

      /**
       * This specifies the factor by which the diagonal of the matrix will be
       * scaled, which sometimes can help to get better preconditioners.
       */
      double ic_rtol;

      /**
       * This determines how large the overlap of the local matrix portions on
       * each processor in a parallel application should be.
       */
      unsigned int overlap;
    };

    /**
     * Initialize function. Takes the matrix the preconditioner should be
     * computed of, and additional flags if there are any.
     */
    void initialize (const SparseMatrix   &matrix,
                     const AdditionalData &additional_data = AdditionalData());
  };



  /**
   * A wrapper class for an incomplete LU factorization (ILU)
   * preconditioner for Trilinos matrices. This preconditioner works
   * both in serial and in parallel, depending on the matrix it is based
   * on. In general, an incomplete factorization does not take all
   * fill-in elements that would appear in a full factorization (that is
   * the basis for a direct solve). Trilinos allows to set the amount of
   * fill-in elements, governed by the additional data argument
   * <tt>ilu_fill</tt>, so one can gradually choose between a
   * factorization on the sparse matrix structure only
   * (<tt>ilu_fill=0</tt>) to a full factorization (<tt>ilu_fill</tt> in
   * the range of 10 to 50, depending on the spatial dimension of the
   * PDE problem and the degree of the finite element basis functions;
   * generally, more required fill-in elements require this parameter to
   * be set to a higher integer value).
   *
   * The AdditionalData data structure allows to set preconditioner
   * options. Besides the fill-in argument, these options are some
   * options for perturbations (see the documentation of the
   * AdditionalData structure for details), and a parameter
   * <tt>overlap</tt> that determines if and how much overlap there
   * should be between the matrix partitions on the various MPI
   * processes. The default settings are 0 for the additional fill-in, 0
   * for the absolute augmentation tolerance, 1 for the relative
   * augmentation tolerance, 0 for the overlap.
   *
   * Note that a parallel application of the ILU preconditioner is
   * actually a block-Jacobi preconditioner with block size equal to the
   * local matrix size. Spoken more technically, this parallel operation
   * is an <a
   * href="http://en.wikipedia.org/wiki/Additive_Schwarz_method">additive
   * Schwarz method</a> with an ILU <em>approximate solve</em> as inner
   * solver, based on the (outer) parallel partitioning.
   *
   * @ingroup TrilinosWrappers
   * @ingroup Preconditioners
   * @author Martin Kronbichler, 2008
   */
  class PreconditionILU : public PreconditionBase
  {
  public:
    /**
     * Standardized data struct to pipe additional parameters to the
     * preconditioner. The Trilinos ILU decomposition allows for some fill-in,
     * so it actually is a threshold incomplete LU factorization. The amount
     * of fill-in, and hence, the amount of memory used by this
     * preconditioner, is controlled by the parameter <tt>ilu_fill</tt>, which
     * specifies this as a double. When forming the preconditioner, for
     * certain problems bad conditioning (or just bad luck) can cause the
     * preconditioner to be very poorly conditioned. Hence it can help to add
     * diagonal perturbations to the original matrix and form the
     * preconditioner for this slightly better matrix. <tt>ilu_atol</tt> is an
     * absolute perturbation that is added to the diagonal before forming the
     * prec, and <tt>ilu_rtol</tt> is a scaling factor $rtol \geq 1$. The last
     * parameter specifies the overlap of the partitions when the
     * preconditioner runs in parallel.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default, the level of extra fill-ins is set to be
       * zero (just use the matrix structure, do not generate any additional
       * fill-in), the tolerance level are 0 and 1, respectively, and the
       * overlap in case of a parallel execution is zero. This overlap in a
       * block-application of the ILU in the parallel case makes the
       * preconditioner a so-called additive Schwarz preconditioner.
       */
      AdditionalData (const unsigned int ilu_fill = 0,
                      const double       ilu_atol = 0.,
                      const double       ilu_rtol = 1.,
                      const unsigned int overlap  = 0);

      /**
       * This specifies the amount of additional fill-in elements besides the
       * sparse matrix structure. When <tt>ilu_fill</tt> is large, this means
       * that many fill-ins will be added, so that the ILU preconditioner
       * comes closer to a (direct) sparse LU decomposition. Note, however,
       * that this will drastically increase the memory requirement,
       * especially when the preconditioner is used in 3D.
       */
      unsigned int ilu_fill;

      /**
       * This specifies the amount of an absolute perturbation that will be
       * added to the diagonal of the matrix, which sometimes can help to get
       * better preconditioners.
       */
      double ilu_atol;

      /**
       * This specifies the factor by which the diagonal of the matrix will be
       * scaled, which sometimes can help to get better preconditioners.
       */
      double ilu_rtol;

      /**
       * This determines how large the overlap of the local matrix portions on
       * each processor in a parallel application should be.
       */
      unsigned int overlap;
    };

    /**
     * Initialize function. Takes the matrix which is used to form the
     * preconditioner, and additional flags if there are any.
     */
    void initialize (const SparseMatrix   &matrix,
                     const AdditionalData &additional_data = AdditionalData());
  };






  /**
   * A wrapper class for a thresholded incomplete LU factorization (ILU-T)
   * preconditioner for Trilinos matrices. This preconditioner works both in
   * serial and in parallel, depending on the matrix it is based on. In
   * general, an incomplete factorization does not take all fill-in elements
   * that would appear in a full factorization (that is the basis for a direct
   * solve). For the ILU-T precondtioner, the parameter <tt>ilut_drop</tt>
   * lets the user specify which elements should be dropped (i.e., should not
   * be part of the incomplete decomposition). Trilinos calculates first the
   * complete factorization for one row, and then skips those elements that
   * are lower than the threshold. This is the main difference to the
   * non-thresholded ILU preconditioner, where the parameter
   * <tt>ilut_fill</tt> governs the incomplete factorization structure. This
   * parameter is available here as well, but provides only some extra
   * information here.
   *
   * The AdditionalData data structure allows to set preconditioner
   * options. Besides the fill-in arguments, these options are some options
   * for perturbations (see the documentation of the AdditionalData structure
   * for details), and a parameter <tt>overlap</tt> that determines if and how
   * much overlap there should be between the matrix partitions on the various
   * MPI processes. The default settings are 0 for the additional fill-in, 0
   * for the absolute augmentation tolerance, 1 for the relative augmentation
   * tolerance, 0 for the overlap.
   *
   * Note that a parallel application of the ILU-T preconditioner is
   * actually a block-Jacobi preconditioner with block size equal to the
   * local matrix size. Spoken more technically, this parallel operation
   * is an <a
   * href="http://en.wikipedia.org/wiki/Additive_Schwarz_method">additive
   * Schwarz method</a> with an ILU <em>approximate solve</em> as inner
   * solver, based on the (outer) parallel partitioning.
   *
   * @ingroup TrilinosWrappers
   * @ingroup Preconditioners
   * @author Martin Kronbichler, 2009
   */
  class PreconditionILUT : public PreconditionBase
  {
  public:
    /**
     * Standardized data struct to pipe additional parameters to the
     * preconditioner. The Trilinos ILU-T decomposition allows for some
     * fill-in, so it actually is a threshold incomplete LU factorization. The
     * amount of fill-in, and hence, the amount of memory used by this
     * preconditioner, is controlled by the parameters <tt>ilut_drop</tt> and
     * <tt>ilut_fill</tt>, which specifies a threshold about which values
     * should form the incomplete factorization and the level of additional
     * fill-in. When forming the preconditioner, for certain problems bad
     * conditioning (or just bad luck) can cause the preconditioner to be very
     * poorly conditioned. Hence it can help to add diagonal perturbations to
     * the original matrix and form the preconditioner for this slightly
     * better matrix. <tt>ilut_atol</tt> is an absolute perturbation that is
     * added to the diagonal before forming the prec, and <tt>ilu_rtol</tt> is
     * a scaling factor $rtol \geq 1$. The last parameter specifies the
     * overlap of the partitions when the preconditioner runs in parallel.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default, no element will be dropped, the level of
       * extra fill-ins is set to be zero (just use the matrix structure, do
       * not generate any additional fill-in except the one that results from
       * non-dropping large elements), the tolerance level are 0 and 1,
       * respectively, and the overlap in case of a parallel execution is
       * zero. This overlap in a block-application of the ILU in the parallel
       * case makes the preconditioner a so-called additive Schwarz
       * preconditioner.
       */
      AdditionalData (const double       ilut_drop = 0.,
                      const unsigned int ilut_fill = 0,
                      const double       ilut_atol = 0.,
                      const double       ilut_rtol = 1.,
                      const unsigned int overlap  = 0);

      /**
       * This specifies the relative size of elements which should be dropped
       * when forming an incomplete LU decomposition with threshold.
       */
      double ilut_drop;

      /**
       * This specifies the amount of additional fill-in elements besides the
       * sparse matrix structure. When <tt>ilu_fill</tt> is large, this means
       * that many fill-ins will be added, so that the ILU preconditioner
       * comes closer to a (direct) sparse LU decomposition. Note, however,
       * that this will drastically increase the memory requirement,
       * especially when the preconditioner is used in 3D.
       */
      unsigned int ilut_fill;

      /**
       * This specifies the amount of an absolute perturbation that will be
       * added to the diagonal of the matrix, which sometimes can help to get
       * better preconditioners.
       */
      double ilut_atol;

      /**
       * This specifies the factor by which the diagonal of the matrix will be
       * scaled, which sometimes can help to get better preconditioners.
       */
      double ilut_rtol;

      /**
       * This determines how large the overlap of the local matrix portions on
       * each processor in a parallel application should be.
       */
      unsigned int overlap;
    };

    /**
     * Initialize function. Takes the matrix which is used to form the
     * preconditioner, and additional flags if there are any.
     */
    void initialize (const SparseMatrix   &matrix,
                     const AdditionalData &additional_data = AdditionalData());
  };



  /**
   * A wrapper class for a sparse direct LU decomposition on parallel
   * blocks for Trilinos matrices. When run in serial, this corresponds
   * to a direct solve on the matrix.
   *
   * The AdditionalData data structure allows to set preconditioner
   * options.
   *
   * Note that a parallel application of the block direct solve
   * preconditioner is actually a block-Jacobi preconditioner with block
   * size equal to the local matrix size. Spoken more technically, this
   * parallel operation is an <a
   * href="http://en.wikipedia.org/wiki/Additive_Schwarz_method">additive
   * Schwarz method</a> with an <em>exact solve</em> as inner solver,
   * based on the (outer) parallel partitioning.
   *
   * @ingroup TrilinosWrappers
   * @ingroup Preconditioners
   * @author Martin Kronbichler, 2008
   */
  class PreconditionBlockwiseDirect : public PreconditionBase
  {
  public:
    /**
     * Standardized data struct to pipe additional parameters to the
     * preconditioner.
     */
    struct AdditionalData
    {
      /**
       * Constructor.
       */
      AdditionalData (const unsigned int overlap  = 0);


      /**
       * This determines how large the overlap of the local matrix portions on
       * each processor in a parallel application should be.
       */
      unsigned int overlap;
    };

    /**
     * Initialize function. Takes the matrix which is used to form the
     * preconditioner, and additional flags if there are any.
     */
    void initialize (const SparseMatrix   &matrix,
                     const AdditionalData &additional_data = AdditionalData());
  };






  /**
   * A wrapper class for a Chebyshev preconditioner for Trilinos matrices.
   *
   * The AdditionalData data structure allows to set preconditioner
   * options.
   *
   * @ingroup TrilinosWrappers
   * @ingroup Preconditioners
   * @author Martin Kronbichler, 2008
   */
  class PreconditionChebyshev : public PreconditionBase
  {
  public:
    /**
     * Standardized data struct to pipe additional parameters to the
     * preconditioner.
     */
    struct AdditionalData
    {
      /**
       * Constructor.
       */
      AdditionalData (const unsigned int degree           = 1,
                      const double       max_eigenvalue   = 10.,
                      const double       eigenvalue_ratio = 30.,
                      const double       min_eigenvalue   = 1.,
                      const double       min_diagonal     = 1e-12,
                      const bool         nonzero_starting = false);

      /**
       * This determines the degree of the Chebyshev polynomial. The degree of
       * the polynomial gives the number of matrix-vector products to be
       * performed for one application of the vmult() operation.
       */
      unsigned int degree;

      /**
       * This sets the maximum eigenvalue of the matrix, which needs to be set
       * properly for appropriate performance of the Chebyshev preconditioner.
       */
      double max_eigenvalue;

      /**
       * This sets the ratio between the maximum and the minimum eigenvalue.
       */
      double eigenvalue_ratio;

      /**
       * This sets the minimum eigenvalue, which is an optional parameter only
       * used internally for checking whether we use an identity matrix.
       */
      double min_eigenvalue;

      /**
       * This sets a threshold below which the diagonal element will not be
       * inverted in the Chebyshev algorithm.
       */
      double min_diagonal;

      /**
       * When this flag is set to <tt>true</tt>, it enables the method
       * <tt>vmult(dst, src)</tt> to use non-zero data in the vector
       * <tt>dst</tt>, appending to it the Chebyshev corrections. This can be
       * useful in some situations (e.g. when used for high-frequency error
       * smoothing), but not the way the solver classes expect a
       * preconditioner to work (where one ignores the content in <tt>dst</tt>
       * for the preconditioner application). The user should really know what
       * she is doing when touching this flag.
       */
      bool nonzero_starting;
    };

    /**
     * Initialize function. Takes the matrix which is used to form the
     * preconditioner, and additional flags if there are any.
     */
    void initialize (const SparseMatrix   &matrix,
                     const AdditionalData &additional_data = AdditionalData());
  };



  /**
   * This class implements an algebraic multigrid (AMG) preconditioner based
   * on the Trilinos ML implementation, which is a black-box preconditioner
   * that works well for many PDE-based linear problems.  What this class does
   * is twofold.  When the initialize() function is invoked, a ML
   * preconditioner object is created based on the matrix that we want the
   * preconditioner to be based on. A call of the respective
   * <code>vmult</code> function does call the respective operation in the
   * Trilinos package, where it is called <code>ApplyInverse</code>. Use of
   * this class is explained in the step-31 tutorial program.
   *
   * Since the Trilinos objects we want to use are heavily dependent on Epetra
   * objects, we recommend using this class in conjunction with Trilinos
   * (Epetra) sparse matrices and vectors. There is support for use with
   * matrices of the deal.II::SparseMatrix class and corresponding vectors,
   * too, but this requires generating a copy of the matrix, which is slower
   * and takes (much) more memory. When doing such a copy operation, we can
   * still profit from the fact that some of the entries in the preconditioner
   * matrix are zero and hence can be neglected.
   *
   * The implementation is able to distinguish between matrices from elliptic
   * problems and convection dominated problems. We use the standard options
   * provided by Trilinos ML for elliptic problems, except that we use a
   * Chebyshev smoother instead of a symmetric Gauss-Seidel smoother.  For
   * most elliptic problems, Chebyshev provides a better damping of high
   * frequencies (in the algebraic sense) than Gauss-Seidel (SSOR), and is
   * faster (Chebyshev requires only some matrix-vector products, whereas SSOR
   * requires substitutions which are more expensive). Moreover, Chebyshev is
   * perfectly parallel in the sense that it does not degenerate when used on
   * many processors. SSOR, on the other hand, gets more Jacobi-like on many
   * processors.
   *
   * For proper functionality of this class we recommend using Trilinos v9.0
   * and higher. Older versions may have problems with generating the
   * coarse-matrix structure when using matrices with many nonzero entries per
   * row (i.e., matrices stemming from higher order finite element
   * discretizations).
   *
   * @ingroup TrilinosWrappers
   * @ingroup Preconditioners
   * @author Martin Kronbichler, 2008
   */
  class PreconditionAMG : public PreconditionBase
  {
  public:

    /**
     * A data structure that is used to control details of how the algebraic
     * multigrid is set up. The flags detailed in here are then passed to
     * the Trilinos ML implementation. A structure of the current type are
     * passed to the constructor of PreconditionAMG.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default, we pretend to work on elliptic problems with
       * linear finite elements on a scalar equation.
       */
      AdditionalData (const bool                             elliptic = true,
                      const bool                             higher_order_elements = false,
                      const unsigned int                     n_cycles = 1,
                      const bool                             w_cyle = false,
                      const double                           aggregation_threshold = 1e-4,
                      const std::vector<std::vector<bool> > &constant_modes = std::vector<std::vector<bool> > (0),
                      const unsigned int                     smoother_sweeps = 2,
                      const unsigned int                     smoother_overlap = 0,
                      const bool                             output_details = false,
                      const char                            *smoother_type = "Chebyshev",
                      const char                            *coarse_type = "Amesos-KLU");

      /**
       * Determines whether the AMG preconditioner should be optimized for
       * elliptic problems (ML option smoothed aggregation SA, using a
       * Chebyshev smoother) or for non-elliptic problems (ML option
       * non-symmetric smoothed aggregation NSSA, smoother is SSOR with
       * underrelaxation).
       */
      bool elliptic;

      /**
       * Determines whether the matrix that the preconditioner is built upon
       * is generated from linear or higher-order elements.
       */
      bool higher_order_elements;

      /**
       * Defines how many multigrid cycles should be performed by the
       * preconditioner.
       */
      unsigned int n_cycles;

      /**
       * Defines whether a w-cycle should be used instead of the standard
       * setting of a v-cycle.
       */
      bool w_cycle;

      /**
       * This threshold tells the AMG setup how the coarsening should be
       * performed. In the AMG used by ML, all points that strongly couple
       * with the tentative coarse-level point form one aggregate. The term
       * <em>strong coupling</em> is controlled by the variable
       * <tt>aggregation_threshold</tt>, meaning that all elements that are
       * not smaller than <tt>aggregation_threshold</tt> times the diagonal
       * element do couple strongly.
       */
      double aggregation_threshold;

      /**
       * Specifies the constant modes (near null space) of the matrix. This
       * parameter tells AMG whether we work on a scalar equation (where the
       * near null space only consists of ones) or on a vector-valued
       * equation.
       */
      std::vector<std::vector<bool> > constant_modes;

      /**
       * Determines how many sweeps of the smoother should be performed. When
       * the flag <tt>elliptic</tt> is set to <tt>true</tt>, i.e., for
       * elliptic or almost elliptic problems, the polynomial degree of the
       * Chebyshev smoother is set to <tt>smoother_sweeps</tt>. The term
       * sweeps refers to the number of matrix-vector products performed in
       * the Chebyshev case. In the non-elliptic case,
       * <tt>smoother_sweeps</tt> sets the number of SSOR relaxation sweeps
       * for post-smoothing to be performed.
       */
      unsigned int smoother_sweeps;

      /**
       * Determines the overlap in the SSOR/Chebyshev error smoother when run
       * in parallel.
       */
      unsigned int smoother_overlap;

      /**
       * If this flag is set to <tt>true</tt>, then internal information from
       * the ML preconditioner is printed to screen. This can be useful when
       * debugging the preconditioner.
       */
      bool output_details;

      /**
       * Determines which smoother to use for the AMG cycle. Possibilities
       * for smoother_type are the following:
       * <ul>
       *   <li>  "Aztec" </li>
       *   <li>  "IFPACK" </li>
       *   <li>  "Jacobi" </li>
       *   <li>  "ML symmetric Gauss-Seidel" </li>
       *   <li>  "symmetric Gauss-Seidel" </li>
       *   <li>  "ML Gauss-Seidel" </li>
       *   <li>  "Gauss-Seidel" </li>
       *   <li>  "block Gauss-Seidel" </li>
       *   <li>  "symmetric block Gauss-Seidel" </li>
       *   <li>  "Chebyshev" </li>
       *   <li>  "MLS" </li>
       *   <li>  "Hiptmair" </li>
       *   <li>  "Amesos-KLU" </li>
       *   <li>  "Amesos-Superlu" </li>
       *   <li>  "Amesos-UMFPACK" </li>
       *   <li>  "Amesos-Superludist" </li>
       *   <li>  "Amesos-MUMPS" </li>
       *   <li>  "user-defined" </li>
       *   <li>  "SuperLU" </li>
       *   <li>  "IFPACK-Chebyshev" </li>
       *   <li>  "self" </li>
       *   <li>  "do-nothing" </li>
       *   <li>  "IC" </li>
       *   <li>  "ICT" </li>
       *   <li>  "ILU" </li>
       *   <li>  "ILUT" </li>
       *   <li>  "Block Chebyshev" </li>
       *   <li>  "IFPACK-Block Chebyshev" </li>
       * </ul>
       */
      const char *smoother_type;

      /**
       * Determines which solver to use on the coarsest level. The same
       * settings as for the smoother type are possible.
       */
      const char *coarse_type;
    };

    /**
     * Destructor.
     */
    ~PreconditionAMG();


    /**
     * Let Trilinos compute a multilevel hierarchy for the solution of a
     * linear system with the given matrix. The function uses the matrix
     * format specified in TrilinosWrappers::SparseMatrix.
     */
    void initialize (const SparseMatrix   &matrix,
                     const AdditionalData &additional_data = AdditionalData());

    /**
     * Let Trilinos compute a multilevel hierarchy for the solution of a
     * linear system with the given matrix. As opposed to the other initialize
     * function above, this function uses an abstract interface to an object
     * of type Epetra_RowMatrix which allows a user to pass quite general
     * objects to the ML preconditioner.
     *
     * This initialization routine is useful in cases where the operator to be
     * preconditioned is not a TrilinosWrappers::SparseMatrix object but still
     * allows to get a copy of the entries in each of the locally owned matrix
     * rows (method ExtractMyRowCopy) and implements a matrix-vector product
     * (methods Multiply or Apply). An example are operators which provide
     * faster matrix-vector multiplications than possible with matrix entries
     * (matrix-free methods). These implementations can be beneficially
     * combined with Chebyshev smoothers that only perform matrix-vector
     * products. The interface class Epetra_RowMatrix is very flexible to
     * enable this kind of implementation.
     */
    void initialize (const Epetra_RowMatrix &matrix,
                     const AdditionalData   &additional_data = AdditionalData());

    /**
     * Let Trilinos compute a multilevel hierarchy for the solution of a
     * linear system with the given matrix. The function uses the matrix
     * format specified in TrilinosWrappers::SparseMatrix.
     *
     * This function is similar to the one above, but allows the user to set
     * all the options of the Trilinos ML preconditioner. In order to find out
     * about all the options for ML, we refer to the <a
     * href=http://trilinos.sandia.gov/packages/ml/mlguide5.pdf>ML user's
     * guide</a>. In particular, users need to follow the ML instructions in
     * case a vector-valued problem ought to be solved.
     */
    void initialize (const SparseMatrix           &matrix,
                     const Teuchos::ParameterList &ml_parameters);

    /**
     * Let Trilinos compute a multilevel hierarchy for the solution of a
     * linear system with the given matrix. As opposed to the other initialize
     * function above, this function uses an abstract interface to an object
     * of type Epetra_RowMatrix which allows a user to pass quite general
     * objects to the ML preconditioner.
     */
    void initialize (const Epetra_RowMatrix       &matrix,
                     const Teuchos::ParameterList &ml_parameters);

    /**
     * Let Trilinos compute a multilevel hierarchy for the solution of a
     * linear system with the given matrix. This function takes a deal.ii
     * matrix and copies the content into a Trilinos matrix, so the function
     * can be considered rather inefficient.
     */
    template <typename number>
    void initialize (const ::dealii::SparseMatrix<number> &deal_ii_sparse_matrix,
                     const AdditionalData                 &additional_data = AdditionalData(),
                     const double                          drop_tolerance = 1e-13,
                     const ::dealii::SparsityPattern      *use_this_sparsity = 0);

    /**
     * This function can be used for a faster recalculation of the
     * preconditioner construction when the matrix entries underlying the
     * preconditioner have changed, but the matrix sparsity pattern has
     * remained the same. What this function does is taking the already
     * generated coarsening structure, computing the AMG prolongation and
     * restriction according to a smoothed aggregation strategy and then
     * building the whole multilevel hiearchy. This function can be
     * considerably faster than the initialize function, since the coarsening
     * pattern is usually the most difficult thing to do when setting up the
     * AMG ML preconditioner.
     */
    void reinit ();

    /**
     * Destroys the preconditioner, leaving an object like just after having
     * called the constructor.
     */
    void clear ();

    /**
     * Prints an estimate of the memory consumption of this class.
     */
    size_type memory_consumption () const;

  private:
    /**
     * A copy of the deal.II matrix into Trilinos format.
     */
    std_cxx11::shared_ptr<SparseMatrix> trilinos_matrix;
  };



  /**
   * A wrapper class for an identity preconditioner for Trilinos matrices.
   *
   * @ingroup TrilinosWrappers
   * @ingroup Preconditioners
   * @author Bruno Turcksin, 2013
   */
  class PreconditionIdentity : public PreconditionBase
  {
  public:

    /**
     * Apply the preconditioner, i.e., dst = src.
     */
    void vmult (VectorBase       &dst,
                const VectorBase &src) const;

    /**
     * Apply the transport conditioner, i.e., dst = src.
     */
    void Tvmult (VectorBase       &dst,
                 const VectorBase &src) const;

    /**
     * Apply the preconditioner on deal.II data structures
     * instead of the ones provided in the Trilinos wrapper class,
     * i.e., dst = src.
     */
    void vmult (dealii::Vector<double>       &dst,
                const dealii::Vector<double> &src) const;

    /**
     * Apply the transpose preconditioner on deal.II data structures
     * instead of the ones provided in the Trilinos wrapper class,
     * i.e. dst = src.
     */
    void Tvmult (dealii::Vector<double>       &dst,
                 const dealii::Vector<double> &src) const;

    /**
     * Apply the preconditioner on deal.II parallel data structures
     * instead of the ones provided in the Trilinos wrapper class,
     * i.e., dst = src.
     */
    void vmult (parallel::distributed::Vector<double>       &dst,
                const dealii::parallel::distributed::Vector<double> &src) const;

    /**
     * Apply the transpose preconditioner on deal.II parallel data structures
     * instead of the ones provided in the Trilinos wrapper class,
     * i.e., dst = src.
     */
    void Tvmult (parallel::distributed::Vector<double>       &dst,
                 const dealii::parallel::distributed::Vector<double> &src) const;
  };



// -------------------------- inline and template functions ----------------------


#ifndef DOXYGEN

  inline
  void
  PreconditionBase::vmult (VectorBase       &dst,
                           const VectorBase &src) const
  {
    Assert (dst.vector_partitioner().SameAs(preconditioner->OperatorRangeMap()),
            ExcNonMatchingMaps("dst"));
    Assert (src.vector_partitioner().SameAs(preconditioner->OperatorDomainMap()),
            ExcNonMatchingMaps("src"));

    const int ierr = preconditioner->ApplyInverse (src.trilinos_vector(),
                                                   dst.trilinos_vector());
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }

  inline
  void
  PreconditionBase::Tvmult (VectorBase       &dst,
                            const VectorBase &src) const
  {
    Assert (dst.vector_partitioner().SameAs(preconditioner->OperatorRangeMap()),
            ExcNonMatchingMaps("dst"));
    Assert (src.vector_partitioner().SameAs(preconditioner->OperatorDomainMap()),
            ExcNonMatchingMaps("src"));

    preconditioner->SetUseTranspose(true);
    const int ierr = preconditioner->ApplyInverse (src.trilinos_vector(),
                                                   dst.trilinos_vector());
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
    preconditioner->SetUseTranspose(false);
  }


  // For the implementation of the <code>vmult</code> function with deal.II
  // data structures we note that invoking a call of the Trilinos
  // preconditioner requires us to use Epetra vectors as well. We do this by
  // providing a view, i.e., feed Trilinos with a pointer to the data, so we
  // avoid copying the content of the vectors during the iteration (this
  // function is only useful when used in serial anyway). In the declaration
  // of the right hand side, we need to cast the source vector (that is
  // <code>const</code> in all deal.II calls) to non-constant value, as this
  // is the way Trilinos wants to have them.
  inline
  void PreconditionBase::vmult (dealii::Vector<double>       &dst,
                                const dealii::Vector<double> &src) const
  {
    AssertDimension (static_cast<TrilinosWrappers::types::int_type>(dst.size()),
                     preconditioner->OperatorDomainMap().NumMyElements());
    AssertDimension (static_cast<TrilinosWrappers::types::int_type>(src.size()),
                     preconditioner->OperatorRangeMap().NumMyElements());
    Epetra_Vector tril_dst (View, preconditioner->OperatorDomainMap(),
                            dst.begin());
    Epetra_Vector tril_src (View, preconditioner->OperatorRangeMap(),
                            const_cast<double *>(src.begin()));

    const int ierr = preconditioner->ApplyInverse (tril_src, tril_dst);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }


  inline
  void PreconditionBase::Tvmult (dealii::Vector<double>       &dst,
                                 const dealii::Vector<double> &src) const
  {
    AssertDimension (static_cast<TrilinosWrappers::types::int_type>(dst.size()),
                     preconditioner->OperatorDomainMap().NumMyElements());
    AssertDimension (static_cast<TrilinosWrappers::types::int_type>(src.size()),
                     preconditioner->OperatorRangeMap().NumMyElements());
    Epetra_Vector tril_dst (View, preconditioner->OperatorDomainMap(),
                            dst.begin());
    Epetra_Vector tril_src (View, preconditioner->OperatorRangeMap(),
                            const_cast<double *>(src.begin()));

    preconditioner->SetUseTranspose(true);
    const int ierr = preconditioner->ApplyInverse (tril_src, tril_dst);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
    preconditioner->SetUseTranspose(false);
  }



  inline
  void
  PreconditionBase::vmult (parallel::distributed::Vector<double>       &dst,
                           const parallel::distributed::Vector<double> &src) const
  {
    AssertDimension (static_cast<TrilinosWrappers::types::int_type>(dst.local_size()),
                     preconditioner->OperatorDomainMap().NumMyElements());
    AssertDimension (static_cast<TrilinosWrappers::types::int_type>(src.local_size()),
                     preconditioner->OperatorRangeMap().NumMyElements());
    Epetra_Vector tril_dst (View, preconditioner->OperatorDomainMap(),
                            dst.begin());
    Epetra_Vector tril_src (View, preconditioner->OperatorRangeMap(),
                            const_cast<double *>(src.begin()));

    const int ierr = preconditioner->ApplyInverse (tril_src, tril_dst);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }

  inline
  void
  PreconditionBase::Tvmult (parallel::distributed::Vector<double>       &dst,
                            const parallel::distributed::Vector<double> &src) const
  {
    AssertDimension (static_cast<TrilinosWrappers::types::int_type>(dst.local_size()),
                     preconditioner->OperatorDomainMap().NumMyElements());
    AssertDimension (static_cast<TrilinosWrappers::types::int_type>(src.local_size()),
                     preconditioner->OperatorRangeMap().NumMyElements());
    Epetra_Vector tril_dst (View, preconditioner->OperatorDomainMap(),
                            dst.begin());
    Epetra_Vector tril_src (View, preconditioner->OperatorRangeMap(),
                            const_cast<double *>(src.begin()));

    preconditioner->SetUseTranspose(true);
    const int ierr = preconditioner->ApplyInverse (tril_src, tril_dst);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
    preconditioner->SetUseTranspose(false);
  }

#endif

}


/*@}*/


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_TRILINOS

/*----------------------------   trilinos_precondition.h     ---------------------------*/

#endif
/*----------------------------   trilinos_precondition.h     ---------------------------*/
