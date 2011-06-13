//----------------------------  test_01.cc  ---------------------------
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
//----------------------------  test_01.cc  ---------------------------

// a basic test for some serialization functionality

#include "serialization.h"

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
  C p1, p2;
  
  verify (p1, p2);
}


int main ()
{
  std::ofstream logfile("test_01/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();

  deallog << "OK" << std::endl;
}
