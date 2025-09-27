// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test assemble_flags.h

#include <deal.II/meshworker/assemble_flags.h>

#include <fstream>

#include "../tests.h"

using namespace MeshWorker;

int
main()
{
  initlog();
  AssembleFlags flag = assemble_own_cells | assemble_boundary_faces;

  deallog.get_file_stream() << flag << std::endl;
}
