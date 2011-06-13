//----------------------------  pointer_05.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  pointer_05.cc  ---------------------------

// test what happens when serializing two objects that have pointers to the
// same common object.
//
// similar to the _03 and _04 tests, but bury the pointer somewhere

#include "serialization.h"

#include <boost/serialization/utility.hpp>
#include <typeinfo>

int object_number = 1;

class C
{
  public:
    C ()
      {
	object_number = ::object_number++;
	deallog << "Default constructor. Object number "
		<< object_number
		<< std::endl;
      }

    C (const C&)
      {
	object_number = ::object_number++;
	deallog << "copy constructor. Object number "
		<< object_number
		<< std::endl;
      }

    ~C ()
      {
	deallog << "destructor. Object number "
		<< object_number
		<< std::endl;
      }
    
    template <typename Archive>
    void serialize (Archive &ar, const unsigned int version)
      {
	deallog << "Serializing object number "
		<< object_number
		<< " via " << typeid(Archive).name()
		<< std::endl;
      }

    bool operator == (const C &) const
      {
	return true;
      }

  private:
    unsigned int object_number;
};


struct P
{
    C *c;
    template <typename Archive>
    void serialize (Archive &ar, const unsigned int)
      {
	ar & c;
      }

    bool operator == (const P &p) const
      {
					 // there is no useful operation we
					 // can do here for the purposes of
					 // serialization. do the actual tests
					 // below with assertions
	return true;
      }
};



void test ()
{
  {
    C *c = new C();
    std::pair<P,P> pair_1, pair_2;
    pair_1.first.c = c;
    pair_1.second.c = c;

    pair_2.first.c = 0;
    pair_2.second.c = 0;
  
    verify (pair_1, pair_2);

				     // boost::serialize should have
				     // recognized that the two pointers in
				     // pair_1 point to the same object and
				     // consequently re-create only one object
				     // that the two components of the
				     // re-created pair point to
    Assert (pair_2.first.c == pair_2.second.c, ExcInternalError());
    Assert (object_number == 3, ExcInternalError());

				     // destroy the newly created object. this
				     // must succeed and would likely throw
				     // some sort of error if the object
				     // pointed to was on the stack (like 'c')
				     // rather than on the heap
    delete pair_2.first.c;
    delete c;
  }
}


int main ()
{
  std::ofstream logfile("pointer_05/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();

  deallog << "OK" << std::endl;
}
