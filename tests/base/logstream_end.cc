//----------------------------  logstream_end.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2007, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  logstream_end.cc  ---------------------------


// it used to happen that if we destroyed logstream (and presumably
// all objects of the same type) that whatever we had put into with
// operator<< after the last use of std::endl was lost. make sure that
// that isn't the case anymore: logstream should flush whatever it has
// left over when it is destroyed


#include "../tests.h"
#include <base/logstream.h>
#include <fstream>
#include <iomanip>
#include <limits>


int main()
{
  std::ofstream logfile("logstream_end/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  {
    LogStream log;

    log.attach(logfile);
    log.depth_console(0);
    log.threshold_double(1.e-10);
    log.log_thread_id (false);

    log << "This should be printed!";
  }
  deallog << "OK" << std::endl;
}
