//----------------------------  subscriptor.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal authors
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
};



/*
#include <set>

template <class Class>
class ActiveObjectMonitor
{
  public:
    ~ActiveObjectMonitor ();
    
    void register_object (const Class *p);
    void deregister_object (const Class *p);
  private:
    std::set<const Class*> registered_objects;
};

ActiveObjectMonitor<Subscriptor> active_object_monitor;


ActiveObjectMonitor::~ActiveObjectMonitor ()
{
  if (registered_objects.size() > 0)
    {
      for (std::set<const Subscriptor*>::const_iterator i=registered_objects.begin();
	   i!=registered_objects.end(); ++i)
	std::cout << "Object still exists of type "
		  << typeid(**i).name()
		  << std::endl;
      Assert (false, ExcInternalError());
    };
};


void ActiveObjectMonitor::register_object (const Subscriptor *p)
{
  Assert (registered_objects.find(p) == registered_objects.end(),
	  ExcInternalError());
  registered_objects.insert (p);
};


void
ActiveObjectMonitor::deregister_object (const Subscriptor *p)
{
  Assert (registered_objects.find(p) != registered_objects.end(),
	  ExcInternalError());
  registered_objects.erase (registered_objects.find(p));
};
*/




Subscriptor::Subscriptor () :
		counter (0),
		object_info (0)
{};


Subscriptor::Subscriptor (const Subscriptor &) :
		counter (0),
		object_info (0)
{};


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
  Assert (counter == 0, ExcInUse(counter, object_info->name()));
}



Subscriptor & Subscriptor::operator = (const Subscriptor &s)
{
  object_info = s.object_info;
  return *this;
};



void Subscriptor::subscribe () const
{
#ifdef DEBUG
  if (object_info == 0)
    object_info = &typeid(*this);
#endif

  subscription_lock.acquire();
  ++counter;
  subscription_lock.release();
};


void Subscriptor::unsubscribe () const {
  Assert (counter>0, ExcNotUsed());
  subscription_lock.acquire();
  --counter;
  subscription_lock.release();
};


unsigned int Subscriptor::n_subscriptions () const 
{
  return counter;
};
