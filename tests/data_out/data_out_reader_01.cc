// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2006 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


#include <deal.II/base/data_out_base.h>

#include <string>
#include <vector>

#include "../tests.h"

#include "patches.h"

// test DataOutReader::merge

template <int dim, int spacedim>
void
check()
{
  const unsigned int np = 1;

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

  std::ostringstream old_data;
  DataOutBase::write_deal_II_intermediate(
    patches,
    names,
    vectors,
    DataOutBase::Deal_II_IntermediateFlags(),
    old_data);

  DataOutReader<dim, spacedim> data;
  {
    std::istringstream input(old_data.str());
    data.read(input);
  }
  DataOutReader<dim, spacedim> additional_data;
  {
    std::istringstream input(old_data.str());
    additional_data.read(input);
  }

  data.merge(additional_data);

  {
    std::ofstream out2("outfile");
    data.write_deal_II_intermediate(out2);
  }

  cat_file("outfile");
  std::remove("outfile");
}


int
main()
{
  initlog();

  check<1, 1>();
  check<1, 2>();
  check<2, 2>();
  check<2, 3>();
  check<3, 3>();
}
