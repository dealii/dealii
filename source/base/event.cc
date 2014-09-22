// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2013 by the deal.II authors
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


#include <deal.II/base/event.h>

DEAL_II_NAMESPACE_OPEN

//TODO: Thread safety

namespace Algorithms
{
  std::vector<std::string> Event::names;

  Event
  Event::assign(const char *name)
  {
    unsigned int index = names.size();
    names.push_back(name);

    Event result;
    // The constructor generated an
    // object with all flags equal
    // zero. Now we set the new one.
    result.flags[index] = true;

    return result;
  }


  Event::Event ()
    :
    all_true(false),
    flags(names.size(), false)
  {}


  void
  Event::clear ()
  {
    all_true = false;
    std::fill(flags.begin(), flags.end(), false);
  }


  void
  Event::all ()
  {
    all_true = true;
  }

  namespace Events
  {
    const Event initial = Event::assign("Initial");
    const Event remesh = Event::assign("Remesh");
    const Event bad_derivative = Event::assign("Bad Derivative");
    const Event new_time = Event::assign("New Time");
    const Event new_timestep_size = Event::assign("New Time Step Size");
  }
}

DEAL_II_NAMESPACE_CLOSE
