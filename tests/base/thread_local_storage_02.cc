// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// verify that thread local storage works as advertised. like _01 but using
// the initialization with an exemplar

#include "../tests.h"
#include <iomanip>
#include <fstream>

#include <deal.II/base/thread_management.h>
#include <deal.II/base/thread_local_storage.h>


struct X
{
  X ()
  {
    Assert (false, ExcInternalError());
  };
  X (int n)
  {
    deallog << "Creating" << std::endl;
    Assert (n==42, ExcInternalError());
  };
  X (const X &)
  {
    deallog << "Copying" << std::endl;
  };
  ~X ()
  {
    deallog << "Destroying " << std::endl;
  };
  int i;
};

Threads::ThreadLocalStorage<X> *tls_data;

volatile int counter = 0;

void execute (int i)
{
  tls_data->get().i = i;

  // indicate that the TLS object has been
  // accessed
  static Threads::Mutex m;
  {
    Threads::Mutex::ScopedLock l(m);
    ++counter;
  }

  // wait in order to make sure that the
  // thread lives longer than the TLS object
  sleep (5);
}


void test ()
{
  // create a thread local storage object
  X exemplar(42);
  tls_data = new Threads::ThreadLocalStorage<X>(exemplar);

  // create 5 threads and wait for their
  // return. the OS will create 5 individual
  // thread ids, which means that we will
  // create 5 individual thread specific
  // storage locations
  Threads::ThreadGroup<> tg;
  for (unsigned int i=10; i<15; ++i)
    tg += Threads::new_thread (execute, i);

  // spin lock until all threads have created
  // their objects
  while (counter != 5);

  // delete the TLS object. this should also
  // destroy all the objects created so far,
  // namely 5 of them, plus the one copy the
  // TLS object stores; while that happens,
  // we should get output into logstream
  delete tls_data;

  // write something into the output
  // file. this also makes sure that the
  // output file records the order in which
  // this output is generated and in which
  // the TLS objects were destroyed. we want
  // that these objects are destroyed when
  // the TLS object is destroyed, not when
  // the threads exit
  deallog << "Done." << std::endl;

  // now make sure the threads all finish
  tg.join_all ();

  // at this point, the seventh object will
  // be destroyed, which is the exemplar
  // object local to this function
}




int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}
