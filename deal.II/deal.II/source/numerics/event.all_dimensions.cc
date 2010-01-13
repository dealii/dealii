//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <numerics/event.h>

DEAL_II_NAMESPACE_OPEN

//TODO: Thread safety

namespace Algorithms
{
  std::vector<std::string> Event::names;

  Event
  Event::assign(const char* name)
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
    const Event bad_derivative = Event::assign("Bad Derivative");
    const Event new_time = Event::assign("New Time");
    const Event new_time_step_size = Event::assign("New Time Step Size");
  }
}

DEAL_II_NAMESPACE_CLOSE
