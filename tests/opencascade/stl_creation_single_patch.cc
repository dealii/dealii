// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2020 by the deal.II authors
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


// Create an STL starting from a single patch IGES file.

#include <deal.II/opencascade/utilities.h>

#include <Standard_Stream.hxx>
#include <TopTools.hxx>
#include <TopoDS_Shape.hxx>

#include "../tests.h"

using namespace OpenCASCADE;

int
main()
{
  initlog();

  const std::string open_iges =
    SOURCE_DIR "/iges_files/DTMB-5415_bulbous_bow.iges";
  const std::string output_open_stl = "output.stl";
  auto              shape_bow       = OpenCASCADE::read_IGES(open_iges, 1e-3);
  write_STL(shape_bow, output_open_stl, 0.5, false, 1e-6, false, 0.5, false);
  deallog << "OK" << std::endl;
}
