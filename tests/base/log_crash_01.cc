//----------------------------  logtest.cc  ---------------------------
//    $Id: logtest.cc 23710 2011-05-17 04:50:10Z bangerth $
//    Version: $Name$ 
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  logtest.cc  ---------------------------


// document crash in deallog related to missing newline


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <fstream>
#include <iomanip>
#include <limits>


int main()
{
  std::ofstream logfile("log_crash_01/output");
  deallog.attach(logfile);
  deallog << "no newline here!";
}
