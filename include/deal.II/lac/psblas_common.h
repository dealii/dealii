// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2019 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#ifndef dealii_psblas_common_h
#define dealii_psblas_common_h

#include <deal.II/base/config.h>

#include "deal.II/base/exception_macros.h"
#include <deal.II/base/exceptions.h>


#ifdef DEAL_II_WITH_PSBLAS

#  include <psb_base_cbind.h>
#  include <psb_c_base.h>
#  include <psb_c_dbase.h>

DEAL_II_NAMESPACE_OPEN

namespace PSCToolkitWrappers
{

  namespace internal
  {
    /*
     * Custom deleter for PSBLAS descriptor. This will be invoked automatically
     * by the shared pointer.
     */
    struct DescriptorDeleter
    {
      void
      operator()(psb_c_descriptor *p) const
      {
        if (p != nullptr)
          psb_c_cdfree(p);
      }
    };

    /**
     * Enum to indicate the state of the vector (building or assembled).
     */
    enum State
    {
      /**
       * State entered after the default constructor, before any allocation.
       * In this state, no operations are possible.
       */
      Default,

      /**
       * State entered after the first allocation, and before the first
       * assembly; in this state it is possible to add communication
       * requirements among different processes.
       */
      Build,

      /**
       * State entered after the assembly; computations such as matrix-vector
       * products, are only possible in this state.
       */
      Assembled
    };

  } // namespace internal

  /**
   * An error occurred code during the initialization of a PSBLAS distributed
   * vector.
   */
  DeclException1(ExcInitializePSBLASVector,
                 int,
                 << "An error with error number " << arg1
                 << " occurred while initializing a PSBLAS vector.");

  /**
   * An error occurred while freeing (deallocating) a PSBLAS distributed vector.
   */
  DeclException1(ExcFreePSBLASVector,
                 int,
                 << "An error with error number " << arg1
                 << " occurred while freeing a PSBLAS vector.");

  /**
   * An error occurred while initializing a communication descriptor.
   */
  DeclException1(ExcInitializePSBLASDescriptor,
                 int,
                 << "An error with error number " << arg1
                 << " occurred while initializing a PSBLAS descriptor.");

  /**
   * An error occurred while assembling (finalizing) a PSBLAS distributed
   * vector.
   */
  DeclException1(ExcAssemblePSBLASVector,
                 int,
                 << "An error with error number " << arg1
                 << " occurred while assembling a PSBLAS vector.");

  /**
   * An error occurred while assembling (finalizing) a PSBLAS communication
   * descriptor.
   */
  DeclException1(ExcAssemblePSBLASDescriptor,
                 int,
                 << "An error with error number " << arg1
                 << " occurred while assembling a PSBLAS descriptor.");


  /**
   * The object is in a state not suitable for this operation.
   */
  DeclException1(ExcInvalidState,
                 int,
                 << "Vector's state is invalid. It is in "
                 << (arg1 == 0 ? "Default" :
                     arg1 == 1 ? "Build" :
                                 "Assembled")
                 << " state. Did you forget to call reinit() or compress()?");

  /**
   * This exception is thrown when an operation that requires the object to be
   * in the Assembled state is called while the object is in a different state
   * (Default or Build). The argument encodes the actual state: 0 = Default,
   * 1 = Build.
   */
  DeclException1(ExcInvalidStateAssembled,
                 int,
                 << "Vector's state is invalid. It should be in "
                 << "Assembled state but it is in "
                 << (arg1 == 0 ? "Default" : "Build")
                 << " state. Did you forget to call compress()?");

  /**
   * This exception is thrown when an operation that requires the object to be
   * in the Build or Assembled state is called while the object is still in
   * the Default state, i.e., before reinit() has been called.
   */
  DeclExceptionMsg(
    ExcInvalidDefault,
    "Vector's state is invalid. It should be in "
    "Build or Assembled state but it is in Default state. Did you forget"
    " to call reinit()?");

  /**
   * This exception is thrown when a generic PSBLAS C interface function
   * returns a nonzero error code. The first argument is the error code and
   * the second argument is the name of the PSBLAS function that failed.
   */
  DeclException2(ExcCallingPSBLASFunction,
                 int,
                 std::string,
                 << "An error with error number " << arg1
                 << " occurred while calling a PSBLAS function: " << arg2
                 << std::endl);

  /**
   * An error occurred while performing an AXPBY (scaled vector addition)
   * operation.
   */
  DeclException1(ExcAXPBY,
                 int,
                 << "An error with error number " << arg1
                 << " occurred while performing the AXPBY operation.");

  /**
   * An error occurred while inserting or setting entries in a distributed
   * vector.
   */
  DeclException1(ExcInsertionInPSBLASVector,
                 int,
                 << "An error with error number " << arg1
                 << " occurred while inserting values into a PSBLAS vector.");

  /**
   * This exception is thrown when the caller attempts a 'set' or 'add'
   * operation on a vector that is currently in the opposite mode. Both
   * modes cannot be mixed without an intervening call to compress().
   * The first argument is the requested mode and the second is the current
   * mode (1 = set, 2 = add).
   */
  DeclException2(ExcWrongMode,
                 int,
                 int,
                 << "You tried to do a "
                 << (arg1 == 1 ? "'set'" : (arg1 == 2 ? "'add'" : "???"))
                 << " operation but the vector is currently in "
                 << (arg2 == 1 ? "'set'" : (arg2 == 2 ? "'add'" : "???"))
                 << " mode. You first have to call 'compress()'.");

  /**
   * Error while accessing an element of a distributed vector that is not stored
   * locally.
   */
  DeclException3(
    ExcAccessToNonlocalElement,
    int,
    int,
    int,
    << "You tried to access element " << arg1
    << " of a distributed vector, but only elements in range [" << arg2 << ','
    << arg3 << "] are stored locally and can be accessed."
    << "\n\n"
    << "A common source for this kind of problem is that you "
    << "are passing a 'fully distributed' vector into a function "
    << "that needs read access to vector elements that correspond "
    << "to degrees of freedom on ghost cells (or at least to "
    << "'locally active' degrees of freedom that are not also "
    << "'locally owned'). You need to pass a vector that has these "
    << "elements as ghost entries.");

  /**
   * An error occurred while allocating a distributed sparse matrix.
   */
  DeclException1(ExcAllocationPSBLASMatrix,
                 int,
                 << "An error with error number " << arg1
                 << " occurred while allocating a PSBLAS sparse matrix.");

  /**
   * An error occurred while assembling (finalizing) a distributed sparse
   * matrix, i.e., during the transition from the Build state to the Assembled
   * state.
   */
  DeclException1(ExcAssemblePSBLASMatrix,
                 int,
                 << "An error with error number " << arg1
                 << " occurred while assembling a PSBLAS sparse matrix.");

  /**
   * An error occurred while inserting or setting entries in a distributed
   * sparse matrix.
   */
  DeclException1(
    ExcInsertionInPSBLASMatrix,
    int,
    << "An error with error number " << arg1
    << " occurred while inserting values into a PSBLAS sparse matrix.");

  /**
   * An error occurred while performing a sparse matrix-vector product.
   */
  DeclException1(
    ExcMatVecPSBLAS,
    int,
    << "An error with error number " << arg1
    << " occurred while performing a matrix-vector operation with a PSBLAS sparse matrix.");

  /**
   * The source and destination objects of an operation are the same object,
   * but the operation requires them to be distinct (e.g., a vector copy or a
   * matrix-vector product into itself).
   */
  DeclExceptionMsg(ExcSourceEqualsDestination,
                   "You are attempting an operation on two vectors that "
                   "are the same object, but the operation requires that the "
                   "two objects are in fact different.");


} // namespace PSCToolkitWrappers

DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PSBLAS
#endif
