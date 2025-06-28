// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Try to compile a file that includes PSBLAS headers, and do an MPI hello world
// program using psblas interfaces.

#include <fstream>

#include "../tests.h"

#ifdef DEAL_II_WITH_64BIT_INDICES
#  define IPK4
#  define LPK8
#else
#  define IPK4
#  define LPK4
#endif

#include "psb_base_cbind.h"

int
main(int argc, char **argv)
{
  psb_i_t     iam, np;
  psb_c_ctxt *cctxt;
  cctxt = psb_c_new_ctxt();
  psb_c_init(cctxt);
  psb_c_info(*cctxt, &iam, &np);

  // We don't handle mpi ourselves, so we output stuff manually
  std::string   ofname = "output_" + std::to_string(iam);
  std::ofstream output(ofname);

  output << "iam = " << iam << ", np = " << np << std::endl;
  output.close();

  // Add MPI barrier
  psb_c_barrier(*cctxt);

  // Concatenate output_i to the file "output"
  if (iam == 0)
    {
      std::ofstream final_output("output");
      for (int i = 0; i < np; i++)
        {
          std::string   ofname = "output_" + std::to_string(i);
          std::ifstream input(ofname);
          final_output << input.rdbuf();
          input.close();
          std::remove(ofname.c_str());
        }
    }

  psb_c_exit(*cctxt);
  free(cctxt);
}
