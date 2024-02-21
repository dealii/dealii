// ------------------------------------------------------------------------
//
// Copyright (C) 2020 - 2024 by the deal.II Authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test BoundingBoxDataOut with an rtree of std::pair<BoundingBox<dim>, int>,
// and output also as data the second entry of the pair.

#include <deal.II/base/bounding_box_data_out.h>
#include <deal.II/base/geometry_info.h>

#include <deal.II/numerics/rtree.h>

#include "../tests.h"


template <int dim>
void
test()
{
  const unsigned int N = 3;
  // Use a non-trivial object, convertible to BoundingBox by boost:
  std::vector<std::pair<BoundingBox<dim>, unsigned int>> boxes(N);
  auto unit = create_unit_bounding_box<dim>();

  Point<dim> ones;
  for (unsigned int i = 0; i < dim; ++i)
    ones[i] = 1;

  unsigned int i = 0;
  for (auto &box : boxes)
    {
      const auto c = random_point<dim>();
      const auto d = random_value();
      box          = std::make_pair(BoundingBox<dim>({Point<dim>(c - d * ones),
                                                      Point<dim>(c + d * ones)}),
                           i++);
    }

  const auto tree = pack_rtree(boxes);

  std::vector<std::vector<double>> datasets(boxes.size(),
                                            std::vector<double>(1));

  {
    unsigned i = 0;
    for (auto p : tree)
      datasets[i++][0] = p.second;
  }

  std::string fname = "boxes_" + std::to_string(dim) + ".vtk";
  {
    std::ofstream            ofile(fname);
    BoundingBoxDataOut<dim>  data_out;
    DataOutBase::VtkFlags    flags;
    std::vector<std::string> names = {"id"};
    flags.print_date_and_time      = false;
    data_out.set_flags(flags);
    data_out.build_patches(tree);
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
