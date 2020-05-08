// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2020 by the deal.II authors
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


#include <deal.II/base/data_out_base.h>

#include <string>
#include <vector>

#include "../tests.h"

#include "patches.h"

// Check that the compression levels supported by VtkFlags are valid (that is,
// they are successfully passed to zlib).

template <int dim, int spacedim>
void
check(DataOutBase::VtkFlags flags, std::ostream &out)
{
  const unsigned int np = 4;

  std::vector<DataOutBase::Patch<dim, spacedim>> patches(np);

  create_patches(patches);

  std::vector<std::string> names(5);
  names[0] = "x1";
  names[1] = "x2";
  names[2] = "x3";
  names[3] = "x4";
  names[4] = "i";
  std::vector<
    std::tuple<unsigned int,
               unsigned int,
               std::string,
               DataComponentInterpretation::DataComponentInterpretation>>
    vectors;
  DataOutBase::write_vtu(patches, names, vectors, flags, out);
}


template <int dim, int spacedim>
void
check_all(std::ostream &log)
{
  for (unsigned int i = 0; i < 4; ++i)
    {
      DataOutBase::VtkFlags flags;
      switch (i)
        {
          case (0):
            flags.compression_level = DataOutBase::VtkFlags::no_compression;
            break;
          case (1):
            flags.compression_level = DataOutBase::VtkFlags::best_speed;
            break;
          case (2):
            flags.compression_level = DataOutBase::VtkFlags::best_compression;
            break;
          case (3):
            flags.compression_level =
              DataOutBase::VtkFlags::default_compression;
            break;
          default:
            Assert(false, ExcInternalError());
        }
      log << "==============================\n"
          << dim << spacedim << "-" << flags.compression_level << ".vtu"
          << "\n==============================\n";
      check<dim, spacedim>(flags, log);
    }
}

int
main()
{
  std::ofstream logfile("output");
  check_all<1, 1>(logfile);
  check_all<2, 2>(logfile);
  check_all<3, 3>(logfile);
}
