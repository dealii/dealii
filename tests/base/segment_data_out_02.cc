// ------------------------------------------------------------------------
//
// Copyright (C) 2025 by the deal.II Authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test SegmentDataOut with line segments and attached data.

#include <deal.II/base/segment_data_out.h>

#include "../tests.h"


template <int dim>
void
test()
{
  const unsigned int                             N = 3;
  std::vector<std::pair<Point<dim>, Point<dim>>> segments(N);
  std::vector<std::vector<double>> datasets(N, std::vector<double>(2));

  // Create some test line segments with associated data
  for (unsigned int i = 0; i < N; ++i)
    {
      Point<dim> start, end;
      for (unsigned int d = 0; d < dim; ++d)
        {
          start[d] = static_cast<double>(i);
          end[d]   = static_cast<double>(i + 1);
        }
      segments[i] = std::make_pair(start, end);

      // Attach some data: segment id and length
      datasets[i][0] = static_cast<double>(i);
      datasets[i][1] = start.distance(end);
    }

  std::string fname = "segments_with_data_" + std::to_string(dim) + ".vtk";
  {
    std::ofstream            ofile(fname);
    SegmentDataOut<dim>      data_out;
    DataOutBase::VtkFlags    flags;
    std::vector<std::string> names = {"segment_id", "length"};
    flags.print_date_and_time      = false;
    data_out.set_flags(flags);
    data_out.build_patches(segments);
    data_out.add_datasets(datasets, names);
    data_out.write_vtk(ofile);
  }
  cat_file(fname.c_str());
}


int
main()
{
  initlog();

  test<1>();
  test<2>();
  test<3>();
}
