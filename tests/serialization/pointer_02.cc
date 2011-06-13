//----------------------------  pointer_02.cc  ---------------------------
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
//----------------------------  pointer_02.cc  ---------------------------

// test what happens when serializing a pointer. is a new object created when
// loading into another pointer that is non-NULL and the old pointer
// destroyed? Or is the old object pointed to being co-opted? the former is in
// fact what happens, and to make things just ever so slightly more awkward,
// the previous object pointed to isn't freed but is left dangling, creating
// the potential for a memory leak

#include "serialization.h"

#include <typeinfo>

int object_number = 1;
int objects_destroyed = 0;

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
	++objects_destroyed;
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


void test ()
{
  {
    C *p1 = new C();
    C *p2 = new C();
  
    verify (p1, p2);

    Assert (p1 != p2, ExcInternalError());

    delete p1;
    delete p2;
  }

				   // as mentioned above, p2 is overwritten by
				   // a pointer to a new object, leaving the
				   // original object pointed to as a memory
				   // leak. assert that this behavior persists
  Assert (objects_destroyed == 2, ExcInternalError());
}


int main ()
{
  std::ofstream logfile("pointer_02/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();

  deallog << "OK" << std::endl;
}
