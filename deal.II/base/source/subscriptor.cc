//----------------------------  subscriptor.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2005 by the deal authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  subscriptor.cc  ---------------------------


#include <base/thread_management.h>
#include <base/subscriptor.h>

#include <typeinfo>
#include <string>
#include <iostream>

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
  Threads::ThreadMutex subscription_lock;
}


static const char* unknown_subscriber = "unknown subscriber";


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
  std::string infostring;
  for (map_iterator it = counter_map.begin(); it != counter_map.end(); ++it)
    {
      if (it->second > 0)
	infostring += std::string("\n  from Subscriber ")
		      + std::string(it->first);
    }
  
  Assert (counter == 0, ExcInUse(counter,
				 object_info->name(),
				 infostring));
#endif
}



Subscriptor & Subscriptor::operator = (const Subscriptor &s)
{
  object_info = s.object_info;
  return *this;
}

// These are the implementations for debug mode. The optimized
// versions are inlined in the header file.

#ifdef DEBUG
void Subscriptor::subscribe (const char* id) const
{
  if (object_info == 0)
    object_info = &typeid(*this);
  Threads::ThreadMutex::ScopedLock lock (subscription_lock);
  ++counter;
  
#if DEAL_USE_MT == 0
  const char* const name = (id != 0) ? id : unknown_subscriber;
  
  map_iterator it = counter_map.find(name);
  if (it == counter_map.end())
    counter_map.insert(map_value_type(name, 1U));
  
  else
    it->second++;
#endif
}


void Subscriptor::unsubscribe (const char* id) const
{
  Assert (counter>0, ExcNotUsed());
  Threads::ThreadMutex::ScopedLock lock (subscription_lock);
  --counter;
  
#if DEAL_USE_MT == 0
  const char* name = (id != 0) ? id : unknown_subscriber;

  map_iterator it = counter_map.find(name);
  Assert (it != counter_map.end(), ExcNoSubscriber(object_info->name(), name));
  Assert (it->second > 0, ExcNoSubscriber(object_info->name(), name));
  
  it->second--;
#endif
}
#endif


unsigned int Subscriptor::n_subscriptions () const
{
  return counter;
}
