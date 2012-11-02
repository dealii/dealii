//---------------------------------------------------------------------------
//    $Id: linear_algebra.h 27260 2012-10-31 14:38:43Z heister $
//
//    Copyright (C) 2008, 2009, 2010, 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__parallel_linear_algebra_h
#define __deal2__parallel_linear_algebra_h

#include <deal.II/base/config.h>


#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebraDealII
{
  typedef Vector<double> Vector;
  typedef BlockVector<double> BlockVector;

  typedef SparseMatrix<double> SparseMatrix;

  typedef PreconditionSSOR<SparseMatrix > PreconditionSSOR;
  
}


DEAL_II_NAMESPACE_CLOSE


#ifdef DEAL_II_USE_PETSC

#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_block_sparse_matrix.h>
#include <deal.II/lac/petsc_precondition.h>

DEAL_II_NAMESPACE_OPEN

  namespace LinearAlgebraPETSc
  {
    using namespace dealii;

  typedef PETScWrappers::Vector Vector;
  typedef PETScWrappers::BlockVector BlockVector;

  typedef PETScWrappers::SparseMatrix SparseMatrix;
    typedef PETScWrappers::PreconditionSSOR PreconditionSSOR;
    

    namespace MPI
    {
      
    /**
     * Typedef for the vector type used.
     */
    typedef PETScWrappers::MPI::Vector Vector;

    /**
     * Typedef for the type used to describe vectors that
     * consist of multiple blocks.
     */
    typedef PETScWrappers::MPI::BlockVector BlockVector;

    /**
     * Typedef for the sparse matrix type used.
     */
    typedef PETScWrappers::MPI::SparseMatrix SparseMatrix;

    /**
     * Typedef for the type used to describe sparse matrices that
     * consist of multiple blocks.
     */
    typedef PETScWrappers::MPI::BlockSparseMatrix BlockSparseMatrix;

    /**
     * Typedef for the AMG preconditioner type used for the
     * top left block of the Stokes matrix.
     */
    typedef PETScWrappers::PreconditionBoomerAMG PreconditionAMG;

    /**
     * Typedef for the Incomplete Cholesky preconditioner used
     * for other blocks of the system matrix.
     */
    typedef PETScWrappers::PreconditionICC PreconditionIC;

    /**
     * Typedef for the Incomplete LU decomposition preconditioner used
     * for other blocks of the system matrix.
     */
    typedef PETScWrappers::PreconditionILU PreconditionILU;
    }
    
  }
DEAL_II_NAMESPACE_CLOSE


#endif // DEAL_II_USE_PETSC

#ifdef DEAL_II_USE_TRILINOS

#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_precondition.h>

DEAL_II_NAMESPACE_OPEN

  namespace LinearAlgebraTrilinos
  {
    using namespace dealii;
    namespace MPI
    {
      
    /**
     * Typedef for the vector type used.
     */
    typedef TrilinosWrappers::MPI::Vector Vector;

    /**
     * Typedef for the type used to describe vectors that
     * consist of multiple blocks.
     */
    typedef TrilinosWrappers::MPI::BlockVector BlockVector;

    /**
     * Typedef for the sparse matrix type used.
     */
    typedef TrilinosWrappers::SparseMatrix SparseMatrix;

    /**
     * Typedef for the type used to describe sparse matrices that
     * consist of multiple blocks.
     */
    typedef TrilinosWrappers::BlockSparseMatrix BlockSparseMatrix;

    /**
     * Typedef for the AMG preconditioner type used for the
     * top left block of the Stokes matrix.
     */
    typedef TrilinosWrappers::PreconditionAMG PreconditionAMG;

    /**
     * Typedef for the Incomplete Cholesky preconditioner used
     * for other blocks of the system matrix.
     */
    typedef TrilinosWrappers::PreconditionIC PreconditionIC;

    /**
     * Typedef for the Incomplete LU decomposition preconditioner used
     * for other blocks of the system matrix.
     */
    typedef TrilinosWrappers::PreconditionILU PreconditionILU;
    }
    
  }

DEAL_II_NAMESPACE_CLOSE


#endif // DEAL_II_USE_TRILINOS



#endif
