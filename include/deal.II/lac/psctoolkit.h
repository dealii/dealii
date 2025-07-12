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

#ifndef dealii_psctoolkit_h
#define dealii_psctoolkit_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_PSBLAS

#include "psb_base_cbind.h"
#include "psb_config.h"

DEAL_II_NAMESPACE_OPEN

namespace PSCToolkit
{
    /**
     * Namespace for PSBLAS communicator functions.
     */
    namespace Communicator
    {
        psb_c_ctxt *Init();
        void Exit(psb_c_ctxt *cctxt);
        void Info(psb_c_ctxt *cctxt, int *iam, int *nproc);
        void Barrier(psb_c_ctxt *cctxt);
    }   // namespace Communicator

} // namespace PSCToolkit


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PSBLAS
#endif // dealii_psctoolkit_h