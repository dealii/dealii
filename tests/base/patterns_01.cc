//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

#include "../tests.h"
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/logstream.h>
#include <memory>

int main()
{
  std::ofstream logfile("patterns_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

				   // create a pattern and let it
				   // output its description
  Patterns::Integer pattern(-1,42);
  const std::string desc = pattern.description();

				   // now let the same class re-create
				   // a pattern object from the
				   // description and verify that the
				   // result is the same as what we
				   // started out with
  Patterns::Integer * pattern2 = Patterns::Integer::create (desc);

  Assert (pattern2 != 0, ExcInternalError());
  Assert (desc == pattern2->description(), ExcInternalError());

  deallog << desc << std::endl;

  delete pattern2;
}
