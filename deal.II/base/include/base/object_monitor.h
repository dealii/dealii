//----------------------------  object_monitor.h  ---------------------------
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
//----------------------------  object_monitor.h  ---------------------------
#ifndef __deal2__object_monitor_h
#define __deal2__object_monitor_h



#include <set>

class Subscriptor;



class ObjectMonitor
{
  public:
    ObjectMonitor ();
    ~ObjectMonitor ();

    void activate ();
    void deactivate ();
    
    void register_object (const Subscriptor *p);
    void deregister_object (const Subscriptor *p);
  private:
    std::set<const Subscriptor*> registered_objects;
    bool activated;
};


extern ObjectMonitor *object_monitor;



#endif



