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


#include <base/subscriptor.h>
#include <typeinfo>



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
