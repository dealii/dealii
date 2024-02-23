// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


int
main()
{
  initlog();
  Triangulation<2, 3> tria;
  GridIn<2, 3>        gridin;
  gridin.attach_triangulation(tria);
  gridin.read_assimp(SOURCE_DIR "/grids/torus.obj", -1, true, 1e-3);

  GridOut go;
  go.write_msh(tria, deallog.get_file_stream());
  return 0;
}
