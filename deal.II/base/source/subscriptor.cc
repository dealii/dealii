//----------------------------  subscriptor.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  subscriptor.cc  ---------------------------


#include <base/subscriptor.h>


Subscriptor::Subscriptor () :
		counter (0)
{};


Subscriptor::Subscriptor (const Subscriptor &) :
		counter (0)
{};


Subscriptor::~Subscriptor ()
{
#ifndef QUIET_SUBSCRIPTOR
  Assert (counter == 0, ExcInUse(counter, classname ));
#else
  Assert (counter == 0, ExcInUse(counter));
#endif
}


Subscriptor & Subscriptor::operator = (const Subscriptor &)
{
  return *this;
};


void Subscriptor::subscribe () const
{
#ifdef DEBUG
#ifndef QUIET_SUBSCRIPTOR
  if(classname.size() == 0)
    classname = string(typeid(*this).name());
#endif
#endif
  ++counter;
};


void Subscriptor::unsubscribe () const {
  Assert (counter>0, ExcNotUsed());
  --counter;
};


unsigned int Subscriptor::n_subscriptions () const 
{
  return counter;
};
