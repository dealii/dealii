// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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

#ifndef dealii_cuda_h
#define dealii_cuda_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_CUDA

#include <cusparse.h>

DEAL_II_NAMESPACE_OPEN
namespace Utilities
{
  /**
   * A namespace for utility structures for CUDA.
   */
  namespace CUDA
  {
    /**
     * This structure creates, stores, and destroys the handles of the different
     * CUDA libraries used inside deal.II.
     */
    struct Handle
    {
      /**
       * Constructor. Create the handles for the different libraries.
       */
      Handle();

      /**
       * Copy constructor is deleted.
       */
      Handle(Handle const &) = delete;

      /**
       * Destructor. Destroy the handles and free all the memory allocated by
       * GrowingVectorMemory.
       */
      ~Handle();

      /**
       * Handle to the cuSPARSE library.
       */
      cusparseHandle_t cusparse_handle;
    };
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif

#endif
