//----------------------------  object_monitor_activator.cc  ---------------------------
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
//----------------------------  object_monitor_activator.cc  ---------------------------


#include <base/object_monitor.h>
#include <base/object_monitor_activator.h>
#include <base/subscriptor.h>


ObjectMonitorActivator::ObjectMonitorActivator ()
{
  if (true)
    Subscriptor b;
  object_monitor->activate();
};


ObjectMonitorActivator::~ObjectMonitorActivator ()
{
  object_monitor->deactivate();
};
