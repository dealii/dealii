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

#include <deal.II/base/logstream.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/thread_management.h>

#include <iostream>
#include <string>
#include <typeinfo>

DEAL_II_NAMESPACE_OPEN


static const char *unknown_subscriber = "unknown subscriber";


std::mutex Subscriptor::mutex;


Subscriptor::Subscriptor()
  : counter(0)
  , object_info(nullptr)
{
  // this has to go somewhere to avoid an extra warning.
  (void)unknown_subscriber;
}



Subscriptor::Subscriptor(const Subscriptor &)
  : counter(0)
  , object_info(nullptr)
{}



Subscriptor::Subscriptor(Subscriptor &&subscriptor) noexcept
  : counter(0)
  , object_info(subscriptor.object_info)
{
  subscriptor.check_no_subscribers();
}



Subscriptor::~Subscriptor()
{
  check_no_subscribers();
  object_info = nullptr;
}


void
Subscriptor::check_no_subscribers() const noexcept
{
  // Check whether there are still subscriptions to this object. If so, output
  // the actual name of the class to which this object belongs, i.e. the most
  // derived class. Note that the name may be mangled, so it need not be the
  // clear-text class name. However, you can obtain the latter by running the
  // c++filt program over the output.
#ifdef DEBUG

  // If there are still active pointers, show a message and kill the program.
  // However, under some circumstances, this is not so desirable. For example,
  // in code like this:
  //
  //     Triangulation tria;
  //     DoFHandler *dh = new DoFHandler(tria);
  //     ...some function that throws an exception
  //
  // the exception will lead to the destruction of the triangulation, but since
  // the dof_handler is on the heap it will not be destroyed. This will trigger
  // an assertion in the triangulation. If we kill the program at this point, we
  // will never be able to learn what caused the problem. In this situation,
  // just display a message and continue the program.
  if (counter != 0)
    {
      if (std::uncaught_exception() == false)
        {
          std::string infostring;
          for (map_iterator it = counter_map.begin(); it != counter_map.end();
               ++it)
            {
              if (it->second > 0)
                infostring +=
                  std::string("\n  from Subscriber ") + std::string(it->first);
            }

          if (infostring == "")
            infostring = "<none>";

          DEAL_II_AssertNothrow(counter == 0,
                                ExcInUse(counter.load(),
                                         object_info->name(),
                                         infostring));
        }
      else
        {
          std::cerr
            << "---------------------------------------------------------"
            << std::endl
            << "An object pointed to by a SmartPointer is being destroyed."
            << std::endl
            << "Under normal circumstances, this would abort the program."
            << std::endl
            << "However, another exception is being processed at the"
            << std::endl
            << "moment, so the program will continue to run to allow"
            << std::endl
            << "this exception to be processed." << std::endl
            << "---------------------------------------------------------"
            << std::endl;
        }
    }
#endif
}



Subscriptor &
Subscriptor::operator=(const Subscriptor &s)
{
  check_no_subscribers();
  object_info = s.object_info;
  return *this;
}



Subscriptor &
Subscriptor::operator=(Subscriptor &&s) noexcept
{
  check_no_subscribers();
  s.check_no_subscribers();
  object_info = s.object_info;
  return *this;
}



void
Subscriptor::subscribe(const char *id) const
{
  std::lock_guard<std::mutex> lock(mutex);

  if (object_info == nullptr)
    object_info = &typeid(*this);
  ++counter;

  const char *const name = (id != 0) ? id : unknown_subscriber;

  map_iterator it = counter_map.find(name);
  if (it == counter_map.end())
    counter_map.insert(map_value_type(name, 1U));

  else
    it->second++;
}


void
Subscriptor::unsubscribe(const char *id) const
{
  const char *name = (id != nullptr) ? id : unknown_subscriber;
  DEAL_II_AssertNothrow(counter > 0,
                        ExcNoSubscriber(object_info->name(), name));
  // This is for the case that we do
  // not abort after the exception
  if (counter == 0)
    return;

  std::lock_guard<std::mutex> lock(mutex);

  map_iterator it = counter_map.find(name);
  if (it == counter_map.end() || it->second == 0)
    {
      DEAL_II_AssertNothrow(it != counter_map.end(),
                            ExcNoSubscriber(object_info->name(), name));
      DEAL_II_AssertNothrow(it->second > 0,
                            ExcNoSubscriber(object_info->name(), name));
    }
  else
    {
      --counter;
      it->second--;
    }
}



unsigned int
Subscriptor::n_subscriptions() const
{
  return counter;
}



void
Subscriptor::list_subscribers() const
{
  list_subscribers(deallog);
}

DEAL_II_NAMESPACE_CLOSE
