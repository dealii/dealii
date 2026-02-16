// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2009 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// it used to happen that if we destroyed logstream (and presumably
// all objects of the same type) that whatever we had put into with
// operator<< after the last use of std::endl was lost. make sure that
// that isn't the case anymore: logstream should flush whatever it has
// left over when it is destroyed


#include <limits>

#include "../tests.h"


int
main()
{
  initlog();

  {
    LogStream log;

    log.attach(deallog.get_file_stream());
    log.log_thread_id(false);

    log << "This should be printed!";
  }
  deallog << "OK" << std::endl;
}
