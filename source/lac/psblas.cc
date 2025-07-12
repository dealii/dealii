// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/logstream.h>
#include <deal.II/lac/psctoolkit.h>


#ifdef DEAL_II_WITH_PSBLAS

DEAL_II_NAMESPACE_OPEN

namespace PSCToolkit
{

    namespace Communicator
    {
    
      /**
       * Initialize the PSBLAS communicator and return a pointer to the context.
       * This function creates the context over the MPI_COMM_WORLD communicator
       *
       * @return Pointer to the initialized PSBLAS context.
       */
      psb_c_ctxt *Init()
      {
        psb_c_ctxt *cctxt = psb_c_new_ctxt();
        psb_c_init(cctxt);
        return cctxt;
      }

      /**
       * Exit the PSBLAS communicator and free the context.
       *
       * @param cctxt Pointer to the PSBLAS context to be exited and freed.
       */
      void Exit(psb_c_ctxt *cctxt)
      {
        psb_c_exit(*cctxt);
        free(cctxt);
      }
      /**
       * Get information about the current process in the PSBLAS communicator.
       *
       * @param cctxt Pointer to the PSBLAS context.
       * @param iam Pointer to store the rank of the current process.
       * @param nproc Pointer to store the total number of processes.
       */
      void Info(psb_c_ctxt *cctxt, int *iam, int* nproc)
      {
        psb_c_info(*cctxt, iam, nproc);
      }

      /**
       * Synchronize all processes in the PSBLAS communicator.
       *
       * @param cctxt Pointer to the PSBLAS context.
       */
      void Barrier(psb_c_ctxt *cctxt)
      {
        psb_c_barrier(*cctxt);
      }

    }

} // namespace PSCToolkit

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PSBLAS
