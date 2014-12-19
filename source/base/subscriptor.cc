// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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

#include <deal.II/base/thread_management.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/logstream.h>

#include <typeinfo>
#include <string>
#include <iostream>

DEAL_II_NAMESPACE_OPEN


namespace
{
// create a lock that might be used to control subscription to and
// unsubscription from objects, as that might happen in parallel.
// since it should happen rather seldom that several threads try to
// operate on different objects at the same time (the usual case is
// that they subscribe to the same object right after thread
// creation), a global lock should be sufficient, rather than one that
// operates on a per-object base (in which case we would have to
// include the huge <thread_management.h> file into the
// <subscriptor.h> file).
  Threads::Mutex subscription_lock;
}


static const char *unknown_subscriber = "unknown subscriber";


Subscriptor::Subscriptor ()
  :
  counter (0),
  object_info (0)
{}


Subscriptor::Subscriptor (const Subscriptor &)
  :
  counter (0),
  object_info (0)
{}


Subscriptor::~Subscriptor ()
{
  // check whether there are still
  // subscriptions to this object. if
  // so, output the actual name of
  // the class to which this object
  // belongs, i.e. the most derived
  // class. note that the name may be
  // mangled, so it need not be the
  // clear-text class name. however,
  // you can obtain the latter by
  // running the c++filt program over
  // the output.
#ifdef DEBUG

  // if there are still active pointers, show
  // a message and kill the program. However,
  // under some circumstances, this is not so
  // desirable. For example, in code like this
  //
  // Triangulation tria;
  // DoFHandler *dh = new DoFHandler(tria);
  // ...some function that throws an exception
  //
  // the exception will lead to the
  // destruction of the triangulation, but
  // since the dof_handler is on the heap it
  // will not be destroyed. This will trigger
  // an assertion in the triangulation. If we
  // kill the program at this point, we will
  // never be able to learn what caused the
  // problem. In this situation, just display
  // a message and continue the program.
  if (counter != 0)
    {
      if (std::uncaught_exception() == false)
        {
          std::string infostring;
          for (map_iterator it = counter_map.begin(); it != counter_map.end(); ++it)
            {
              if (it->second > 0)
                infostring += std::string("\n  from Subscriber ")
                              + std::string(it->first);
            }

          if (infostring == "")
            infostring = "<none>";

          AssertNothrow (counter == 0,
                         ExcInUse (counter, object_info->name(), infostring));
        }
      else
        {
          std::cerr << "---------------------------------------------------------"
                    << std::endl
                    << "An object pointed to by a SmartPointer is being destroyed."
                    << std::endl
                    << "Under normal circumstances, this would abort the program."
                    << std::endl
                    << "However, another exception is being processed at the"
                    << std::endl
                    << "moment, so the program will continue to run to allow"
                    << std::endl
                    << "this exception to be processed."
                    << std::endl
                    << "---------------------------------------------------------"
                    << std::endl;
        }
    }
  // In case we do not abort
  // on error, this will tell
  // do_unsubscribe below that the
  // object is unused now.
  counter = 0;
#endif
}



Subscriptor &Subscriptor::operator = (const Subscriptor &s)
{
  object_info = s.object_info;
  return *this;
}



void
Subscriptor::subscribe(const char *id) const
{
#ifdef DEBUG
  if (object_info == 0)
    object_info = &typeid(*this);
  Threads::Mutex::ScopedLock lock (subscription_lock);
  ++counter;

#  ifndef DEAL_II_WITH_THREADS
  const char *const name = (id != 0) ? id : unknown_subscriber;

  map_iterator it = counter_map.find(name);
  if (it == counter_map.end())
    counter_map.insert(map_value_type(name, 1U));

  else
    it->second++;
#  endif
#else
  (void)id;
#endif
}


void
Subscriptor::unsubscribe(const char *id) const
{
#ifdef DEBUG
  const char *name = (id != 0) ? id : unknown_subscriber;
  AssertNothrow (counter>0, ExcNoSubscriber(object_info->name(), name));
  // This is for the case that we do
  // not abort after the exception
  if (counter == 0)
    return;

  Threads::Mutex::ScopedLock lock (subscription_lock);
  --counter;

#  ifndef DEAL_II_WITH_THREADS
  map_iterator it = counter_map.find(name);
  AssertNothrow (it != counter_map.end(), ExcNoSubscriber(object_info->name(), name));
  AssertNothrow (it->second > 0, ExcNoSubscriber(object_info->name(), name));

  it->second--;
#  endif
#else
  (void)id;
#endif
}



unsigned int Subscriptor::n_subscriptions () const
{
  return counter;
}



void Subscriptor::list_subscribers () const
{
#ifndef DEAL_II_WITH_THREADS
  for (map_iterator it = counter_map.begin();
       it != counter_map.end(); ++it)
    deallog << it->second << '/'
            << counter << " subscriptions from \""
            << it->first << '\"' << std::endl;
#else
  deallog << "No subscriber listing with multithreading" << std::endl;
#endif
}

DEAL_II_NAMESPACE_CLOSE
