//----------------------------  object_monitor.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  object_monitor.cc  ---------------------------


#include <base/object_monitor.h>
#include <base/subscriptor.h>


ObjectMonitor *object_monitor;



ObjectMonitor::ObjectMonitor ():
		activated (false)
{};


ObjectMonitor::~ObjectMonitor ()
{};



void ObjectMonitor::activate ()
{
  activated=true;
};



void ObjectMonitor::deactivate ()
{
  activated = false;
  if (registered_objects.size() > 0)
    {
      for (std::set<const Subscriptor*>::const_iterator i=registered_objects.begin();
	   i!=registered_objects.end(); ++i)
	std::cout << "Object still exists of type "
		  << typeid(**i).name()
		  << std::endl;
      abort ();
    };
};


void ObjectMonitor::register_object (const Subscriptor *p)
{
  if (activated)
    {
      cout << typeid(*p).name() << endl;
      if (registered_objects.find(p) != registered_objects.end())
	abort();
      registered_objects.insert (p);
    };
};


void
ObjectMonitor::deregister_object (const Subscriptor *p)
{
  if (activated)
    {
      if (registered_objects.find(p) == registered_objects.end())
	abort();
      registered_objects.erase (registered_objects.find(p));
    };
};
