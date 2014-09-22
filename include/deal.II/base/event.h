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


#ifndef __deal2__event_h
#define __deal2__event_h

#include <deal.II/base/config.h>

#include <vector>
#include <string>
#include <iostream>

DEAL_II_NAMESPACE_OPEN

namespace Algorithms
{
  /**
   * Objects of this kind are used to notify interior applications of
   * changes provoked by an outer loop. They are handed to the
   * application through Operator::notify() and it is up to the
   * actual application how to handle them.
   *
   * Event is organized as an extensible binary enumerator. Every class
   * can add its own events using assign(). A typical code example is
   *
   * @code
   * class A
   * {
   *   static Event event;
   * };
   *
   * Event A::event = Event::assign("Event for A");
   * @endcode
   */
  class Event
  {
  public:
    /**
     * This function registers a
     * new event type and assigns a
     * unique identifier to it. The
     * result of this function
     * should be stored for later
     * use.
     */
    static Event assign (const char *name);

    /**
     * If you forgot to store the
     * result of assign, here is
     * how to retrieve it knowing
     * the name.
     */
//      static Event find(const std::string& name);

    /**
     * Constructor, generating a
     * clear Event.
     */
    Event ();

    /**
     * Clear all flags
     */
    void clear();

    /**
     * Set all flags
     */
    void all();

    /**
     * Add the flags of the other event
     */
    Event &operator += (const Event &event);

    /**
     * Clear the flags of the other event
     */
    Event &operator -= (const Event &event);

    /**
     * Test whether all the flags
     * set in the other Event are
     * also set in this one.
     */
    bool test (const Event &event) const;

    /**
     * Return <tt>true</tt> if any
     * event is set.
     */
    bool any () const;

    /**
     * List the flags to a stream.
     */
    template <class OS>
    void print (OS &os) const;

    /**
     * List all assigned events.
     */
    template <class OS>
    static void print_assigned (OS &os);

  private:
    /**
     * Sometimes, actions have to
     * be taken by all
     * means. Therefore, if this
     * value is true, test() always
     * returns true.
     */
    bool all_true;

    /**
     * The actual list of events
     */
    std::vector<bool> flags;

    /**
     * The names of registered events
     */
//TODO: This static field must be guarded by a mutex to be thread-safe!
    static std::vector<std::string> names;
  };

  /**
   * Events used by library operators
   */
  namespace Events
  {
    /**
     * The program has just started
     * and everything should be
     * new.
     */
    extern const Event initial;

    /**
     * The mesh has changed.
     */
    extern const Event remesh;

    /**
     * The current derivative leads
     * to slow convergence of
     * Newton's method.
     */
    extern const Event bad_derivative;

    /**
     * The time stepping scheme
     * starts a new time step.
     */
    extern const Event new_time;

    /**
     * The time stepping scheme has changed the time step size.
     */
    extern const Event new_timestep_size;
  }


//----------------------------------------------------------------------//


  inline
  bool
  Event::any () const
  {
    if (all_true) return true;
    for (std::vector<bool>::const_iterator i=flags.begin();
         i != flags.end(); ++i)
      if (*i) return true;
    return false;
  }


  inline
  bool
  Event::test (const Event &event) const
  {

    // First, test all_true in this
    if (all_true) return true;

    const unsigned int n = flags.size();
    const unsigned int m = event.flags.size();
    const unsigned int n_min = (n<m)?n:m;

    // Now, if all_true set in the
    // other, then all must be true
    // in this
    if (event.all_true)
      {
        // Non existing flags are
        // always assumed false
        if (m > n)
          return false;

        // Test all flags separately
        // and return false if one is
        // not set
        for (std::vector<bool>::const_iterator i=flags.begin();
             i != flags.end(); ++i)
          if (!*i) return false;
        // All flags are set
        return true;
      }

    // Finally, compare each flag
    // separately
    for (unsigned int i=0; i<n_min; ++i)
      if (event.flags[i] && !flags[i])
        return false;
    for (unsigned int i=n_min; i<m; ++i)
      if (event.flags[i])
        return false;
    return true;
  }



  inline
  Event &Event::operator += (const Event &event)
  {
    all_true |= event.all_true;
    if (all_true) return *this;

    if (flags.size() < event.flags.size())
      flags.resize(event.flags.size());
    for (unsigned int i=0; i<event.flags.size(); ++i)
      flags[i] = flags[i] || event.flags[i];

    return *this;
  }


  inline
  Event &Event::operator -= (const Event &event)
  {
    if (!event.any()) return *this;

    all_true = false;
    if (event.all_true)
      {
        for (std::vector<bool>::iterator i=flags.begin();
             i != flags.end(); ++i)
          *i = false;
        return *this;
      }

    if (flags.size() < event.flags.size())
      flags.resize(event.flags.size());
    for (unsigned int i=0; i<event.flags.size(); ++i)
      if (event.flags[i]) flags[i] = false;

    return *this;
  }


  template <class OS>
  inline
  void
  Event::print (OS &os) const
  {
    if (all_true)
      os << " ALL";

    for (unsigned int i=0; i<flags.size(); ++i)
      if (flags[i])
        os << ' ' << names[i];
  }


  template <class OS>
  inline
  void
  Event::print_assigned (OS &os)
  {
    for (unsigned int i=0; i<names.size(); ++i)
      os << i << '\t' << names[i] << std::endl;
  }


  /**
   * Output shift operator for
   * events. Calls Event::print().
   *
   * @relates Event
   */
  template <class OS>
  OS &operator << (OS &o, const Event &e)
  {
    e.print(o);
    return o;
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
