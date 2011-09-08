//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2008, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

// test ThreadLocalStorage::operator T&

#include "../tests.h"
#include <iomanip>
#include <fstream>

#include <deal.II/base/thread_management.h>
#include <deal.II/base/thread_local_storage.h>


struct X
{
    Threads::ThreadLocalStorage<int> tls_data;

    X ()
		    :
		    tls_data (42)
      {}

    int f ()
      {
					 // access TLS data and have it
					 // converted to the right data type
					 // without the need to call
					 // tls_data.get()
	return tls_data;
      }
};


void test ()
{
  X x;
  Threads::Thread<int> t;
  t = Threads::new_thread (&X::f, x);

  Assert (t.return_value() == 42,
	  ExcInternalError());
}




int main()
{
  std::ofstream logfile("thread_local_storage_03/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
  deallog << "OK" << std::endl;
}
