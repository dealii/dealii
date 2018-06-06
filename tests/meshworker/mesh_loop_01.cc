// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



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
