// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1998 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/logstream.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <typeinfo>

DEAL_II_NAMESPACE_OPEN


static const char *unknown_subscriber = "unknown subscriber";


std::mutex EnableObserverPointer::mutex;


EnableObserverPointer::EnableObserverPointer(
  EnableObserverPointer &&subscriptor) noexcept
  : counter(0)
  , object_info(subscriptor.object_info)
{
  for (auto *const validity_ptr : subscriptor.validity_pointers)
    *validity_ptr = false;
  subscriptor.validity_pointers.clear();
}



EnableObserverPointer::~EnableObserverPointer()
{
  for (auto *const validity_ptr : validity_pointers)
    *validity_ptr = false;
  object_info = nullptr;
}


void
EnableObserverPointer::check_no_subscribers() const noexcept
{
  // Check whether there are still subscriptions to this object. If so, output
  // the actual name of the class to which this object belongs, i.e. the most
  // derived class. Note that the name may be mangled, so it need not be the
  // clear-text class name. However, you can obtain the latter by running the
  // c++filt program over the output.
  if constexpr (running_in_debug_mode())
    {
      // If there are still active pointers, show a message and kill the
      // program. However, under some circumstances, this is not so desirable.
      // For example, in code like this:
      //
      //     Triangulation tria;
      //     DoFHandler *dh = new DoFHandler(tria);
      //     ...some function that throws an exception
      //
      // the exception will lead to the destruction of the triangulation, but
      // since the dof_handler is on the heap it will not be destroyed. This
      // will trigger an assertion in the triangulation. If we kill the program
      // at this point, we will never be able to learn what caused the problem.
      // In this situation, just display a message and continue the program.
      if (counter != 0)
        {
          if (std::uncaught_exceptions() == 0)
            {
              std::string infostring;
              for (const auto &map_entry : counter_map)
                {
                  if (map_entry.second > 0)
                    infostring +=
                      "\n  from Subscriber " + std::string(map_entry.first);
                }

              if (infostring.empty())
                infostring = "<none>";

              AssertNothrow(counter == 0,
                            ExcInUse(counter.load(),
                                     object_info->name(),
                                     infostring));
            }
          else
            {
              std::cerr
                << "---------------------------------------------------------"
                << std::endl
                << "An object pointed to by a ObserverPointer is being destroyed."
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
    }
}



EnableObserverPointer &
EnableObserverPointer::operator=(EnableObserverPointer &&s) noexcept
{
  for (auto *const validity_ptr : s.validity_pointers)
    *validity_ptr = false;
  s.validity_pointers.clear();
  object_info = s.object_info;
  return *this;
}



void
EnableObserverPointer::subscribe(std::atomic<bool> *const validity,
                                 const std::string       &id) const
{
  std::lock_guard<std::mutex> lock(mutex);

  if (object_info == nullptr)
    object_info = &typeid(*this);
  ++counter;

  const std::string &name = id.empty() ? unknown_subscriber : id;

  ++counter_map[name];

  *validity = true;
  validity_pointers.push_back(validity);
}



void
EnableObserverPointer::unsubscribe(std::atomic<bool> *const validity,
                                   const std::string       &id) const
{
  const std::string &name = id.empty() ? unknown_subscriber : id;

  if (counter == 0)
    {
      AssertNothrow(counter > 0, ExcNoSubscriber(object_info->name(), name));
      // This is for the case that we do not abort after the exception
      return;
    }

  std::lock_guard<std::mutex> lock(mutex);

  map_iterator it = counter_map.find(name);
  if (it == counter_map.end())
    {
      AssertNothrow(it != counter_map.end(),
                    ExcNoSubscriber(object_info->name(), name));
      // This is for the case that we do not abort after the exception
      return;
    }
  if (it->second == 0)
    {
      AssertNothrow(it->second > 0, ExcNoSubscriber(object_info->name(), name));
      // This is for the case that we do not abort after the exception
      return;
    }

  auto validity_ptr_it =
    std::find(validity_pointers.begin(), validity_pointers.end(), validity);
  if (validity_ptr_it == validity_pointers.end())
    {
      AssertNothrow(
        validity_ptr_it != validity_pointers.end(),
        ExcMessage(
          "This EnableObserverPointer object does not know anything about the supplied pointer!"));
      return;
    }

  validity_pointers.erase(validity_ptr_it);
  --counter;
  --it->second;
}



void
EnableObserverPointer::list_subscribers() const
{
  list_subscribers(deallog);
}

DEAL_II_NAMESPACE_CLOSE
