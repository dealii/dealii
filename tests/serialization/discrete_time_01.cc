// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check and illustrate the serialization process for DiscreteTime class

#include <deal.II/base/discrete_time.h>

#include "serialization.h"



void
test()
{
  // save data to archive
  std::ostringstream oss;
  {
    DiscreteTime time(0.0, 1.0, 0.6);
    time.advance_time();

    deallog << "Before serialization start time = " << time.get_start_time()
            << ", end time = " << time.get_end_time()
            << ", current time = " << time.get_current_time()
            << ", next time = " << time.get_next_time()
            << ", previous time = " << time.get_previous_time()
            << ", next step size = " << time.get_next_step_size()
            << ", previous step size = " << time.get_previous_step_size()
            << ", step number = " << time.get_step_number() << std::endl;

    boost::archive::text_oarchive oa(oss, boost::archive::no_header);
    oa << time;

    // archive and stream closed when
    // destructors are called
  }

  // verify correctness of the serialization.
  {
    std::istringstream            iss(oss.str());
    boost::archive::text_iarchive ia(iss, boost::archive::no_header);

    DiscreteTime time(1.0, 2.0, 0.5); // Initialize with garbage values.
    ia >> time;

    deallog << "After serialization start time = " << time.get_start_time()
            << ", end time = " << time.get_end_time()
            << ", current time = " << time.get_current_time()
            << ", next time = " << time.get_next_time()
            << ", previous time = " << time.get_previous_time()
            << ", next step size = " << time.get_next_step_size()
            << ", previous step size = " << time.get_previous_step_size()
            << ", step number = " << time.get_step_number() << std::endl;
  }

  deallog << "OK" << std::endl << std::endl;
}


int
main(int argc, char *argv[])
{
  initlog();
  deallog << std::setprecision(4);

  test();
}
