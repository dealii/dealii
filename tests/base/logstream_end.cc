// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2018 by the deal.II authors
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
