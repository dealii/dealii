// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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

  XDMFEntry entry1(
    mesh_filename, solution_filename, time, nodes, cells, dim, spacedim);
  XDMFEntry entry2;

  // save data to archive
  std::ostringstream oss;
  {
    boost::archive::text_oarchive oa(oss, boost::archive::no_header);
    oa << entry1;
    // archive and stream closed when destructors are called
  }
  deallog << oss.str() << std::endl;

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
