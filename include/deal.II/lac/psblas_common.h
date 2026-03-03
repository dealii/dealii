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
     * Custom deleter for PSBLAS descriptor.
     */
    struct PSBLASDescriptorDeleter
    {
      void
      operator()(psb_c_descriptor *p) const
      {
        if (p)
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
      /*
       * State entered after the assembly; computations such as matrix-vector
       * products, are only possible in this state.
       */
      Assembled
    };

  } // namespace internal

  /**
   * Exception
   */
  DeclException1(ExcInitializePSBLASVector,
                 int,
                 << "An error with error number " << arg1
                 << " occurred while initializing a PSBLAS vector.");

  /**
   * Exception
   */
  DeclException1(ExcFreePSBLASVector,
                 int,
                 << "An error with error number " << arg1
                 << " occurred while freeing a PSBLAS vector.");

  /**
   * Exception
   */
  DeclException1(ExcInitializePSBLASDescriptor,
                 int,
                 << "An error with error number " << arg1
                 << " occurred while initializing a PSBLAS descriptor.");

  /**
   * Exception
   */
  DeclException1(ExcAssemblePSBLASVector,
                 int,
                 << "An error with error number " << arg1
                 << " occurred while assembling a PSBLAS vector.");

  /**
   * Exception
   */
  DeclException1(ExcAssemblePSBLASDescriptor,
                 int,
                 << "An error with error number " << arg1
                 << " occurred while assembling a PSBLAS descriptor.");


  /**
   * Exception
   */
  DeclException1(ExcInvalidState,
                 int,
                 << "Vector's state is invalid. It is in "
                 << (arg1 == 0 ? "Default" :
                     arg1 == 1 ? "Build" :
                                 "Assembled")
                 << " state. Did you forget to call reinit() or compress()?");

  /**
   * Exception
   */
  DeclException1(ExcInvalidStateBuild,
                 int,
                 << "Vector's state is invalid. It should be in "
                 << "Build state but it is in "
                 << (arg1 == 0 ? "Default" : "Assembled")
                 << " state. Did you forget to call reinit()?");

  /**
   * Exception
   */
  DeclException1(ExcInvalidStateAssembled,
                 int,
                 << "Vector's state is invalid. It should be in "
                 << "Assembled state but it is in "
                 << (arg1 == 0 ? "Default" : "Build")
                 << " state. Did you forget to call compress()?");

  /**
   * Exception
   */
  DeclExceptionMsg(
    ExcInvalidDefault,
    "Vector's state is invalid. It should be in "
    "Build or Assembled state but it is in Default state. Did you forget"
    " to reinit it()?");

  /**
   * Exception
   */
  DeclException2(ExcCallingPSBLASFunction,
                 int,
                 std::string,
                 << "An error with error number " << arg1
                 << " occurred while calling a PSBLAS function: " << arg2
                 << std::endl);

  /**
   * Exception
   */
  DeclException1(ExcAXPBY,
                 int,
                 << "An error with error number " << arg1
                 << " occurred while performing the AXPBY operation.");

  /**
   * Exception
   */
  DeclException1(ExcInsertionInPSBLASVector,
                 int,
                 << "An error with error number " << arg1
                 << " occurred while inserting values into a PSBLAS vector.");

  /**
   * Exception
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
   * Exception
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

  // Matrix-related exceptions

  DeclException1(ExcAllocationPSBLASMatrix,
                 int,
                 << "An error with error number " << arg1
                 << " occurred while allocating a PSBLAS sparse matrix.");

  DeclException1(ExcAssemblePSBLASMatrix,
                 int,
                 << "An error with error number " << arg1
                 << " occurred while assembling a PSBLAS sparse matrix.");

  DeclException1(
    ExcInsertionInPSBLASMatrix,
    int,
    << "An error with error number " << arg1
    << " occurred while inserting values into a PSBLAS sparse matrix.");

  DeclException1(
    ExcMatVecPSBLAS,
    int,
    << "An error with error number " << arg1
    << " occurred while performing a matrix-vector operation with a PSBLAS sparse matrix.");


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
