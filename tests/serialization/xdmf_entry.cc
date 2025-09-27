// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check serialization for XDMFEntry

#include <deal.II/base/data_out_base.h>

#include "serialization.h"

void
test()
{
  const std::string mesh_filename     = "mesh.h5";
  const std::string solution_filename = "solution.h5";

  const double       time     = 0.5;
  const unsigned int nodes    = 128;
  const unsigned int cells    = 16;
  const unsigned int dim      = 2;
  const unsigned int spacedim = 3;

  XDMFEntry entry1(mesh_filename,
                   solution_filename,
                   time,
                   nodes,
                   cells,
                   dim,
                   spacedim,
                   ReferenceCells::Quadrilateral);
  XDMFEntry entry2;

  // save data to archive
  std::ostringstream oss;
  {
    boost::archive::text_oarchive oa(oss, boost::archive::no_header);
    oa << entry1;
    // archive and stream closed when destructors are called
  }

  // verify correctness of the serialization
  {
    std::istringstream            iss(oss.str());
    boost::archive::text_iarchive ia(iss, boost::archive::no_header);
    ia >> entry2;
  }

  deallog << "XDMFEntry before serialization: " << std::endl
          << std::endl
          << entry1.get_xdmf_content(0) << std::endl;

  deallog << "XDMFEntry after de-serialization: " << std::endl
          << std::endl
          << entry2.get_xdmf_content(0) << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test();

  deallog << "OK" << std::endl;
}
