// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Try to compile a file that includes PSBLAS headers

extern "C"
{
#include <psb_base_cbind.h>
}

#include "../tests.h"

int
main()
{
  psb_i_t     iam, np;
  psb_c_ctxt *cctxt;
  cctxt = psb_c_new_ctxt();
  psb_c_init(cctxt);
  psb_c_info(*cctxt, &iam, &np);

  mpi_initlog();

  deallog << "iam = " << iam << ", np = " << np << std::endl;

  psb_c_exit(*cctxt);
  free(cctxt);
  return 0;
}
