/*      $Id$                 */

#include <base/subscriptor.h>



Subscriptor::Subscriptor () :
		counter (0)
{};



Subscriptor::Subscriptor (const Subscriptor &) :
		counter (0)
{};



Subscriptor::~Subscriptor () {
  Assert (counter == 0, ExcInUse(counter));
};



Subscriptor & Subscriptor::operator = (const Subscriptor &) {
  return *this;
};



void Subscriptor::subscribe () const {
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
