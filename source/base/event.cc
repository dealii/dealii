// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/event.h>

#include <string>
#include <vector>


DEAL_II_NAMESPACE_OPEN

// TODO: Thread safety

namespace Algorithms
{
  std::vector<std::string> Event::names;

  Event
  Event::assign(const std::string &name)
  {
    unsigned int index = names.size();
    names.emplace_back(name);

    Event result;
    // The constructor generated an
    // object with all flags equal
    // zero. Now we set the new one.
    result.flags[index] = true;

    return result;
  }


  Event::Event()
    : all_true(false)
    , flags(names.size(), false)
  {}


  void
  Event::clear()
  {
    all_true = false;
    std::fill(flags.begin(), flags.end(), false);
  }


  void
  Event::all()
  {
    all_true = true;
  }

  namespace Events
  {
    const Event initial           = Event::assign("Initial");
    const Event remesh            = Event::assign("Remesh");
    const Event bad_derivative    = Event::assign("Bad Derivative");
    const Event new_time          = Event::assign("New Time");
    const Event new_timestep_size = Event::assign("New Time Step Size");
  } // namespace Events
} // namespace Algorithms

DEAL_II_NAMESPACE_CLOSE
