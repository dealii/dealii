/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// test_dll_exported.cpp

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// should pass compilation and execution

// This is an overly complex test.  The purpose of this test is to
// demonstrate and test the ability to serialize a hierarchy of class
// through a base class pointer even though that class might be
// implemented in different dlls and use different extended type info
// systems.
//
// polymorphic_ base is locally declared and defined.  It use the
// "no_rtti" extended type info system.

// polymorphic_derived1 is locally declared and defined.  It uses
// the default "type_id" extended type info system

// polymorphic_derived2 is declared in polymorphic_derived.hpp
// and defined in dll_polymorphic_derived2.  It uses the typeid
// system.

#include <cstddef> // NULL
#include <fstream>

#include <boost/config.hpp>

#include <boost/archive/archive_exception.hpp>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "polymorphic_base.hpp"

#include "polymorphic_derived2.hpp"

// save exported polymorphic class
void save_exported(const char *testfile)
{
    std::ofstream os(testfile);
    boost::archive::text_oarchive oa(os);

    polymorphic_base *rb2 = new polymorphic_derived2;
    polymorphic_derived2 *rd21 = new polymorphic_derived2;

    // export will permit correct serialization
    // through a pointer to a base class
    oa << BOOST_SERIALIZATION_NVP(rb2);
    oa << BOOST_SERIALIZATION_NVP(rd21);

    delete rb2;
    delete rd21;
}

// save exported polymorphic class
void load_exported(const char *testfile)
{
    std::ifstream is(testfile);
    boost::archive::text_iarchive ia(is);

    polymorphic_base *rb2 = NULL;
    polymorphic_derived2 *rd21 = NULL;

    // export will permit correct serialization
    // through a pointer to a base class
    ia >> BOOST_SERIALIZATION_NVP(rb2);
    assert(
        boost::serialization::type_info_implementation<polymorphic_derived2>
            ::type::get_const_instance()
        ==
        * boost::serialization::type_info_implementation<polymorphic_base>
            ::type::get_const_instance().get_derived_extended_type_info(*rb2));
    ia >> BOOST_SERIALIZATION_NVP(rd21);
    assert(
        boost::serialization::type_info_implementation<polymorphic_derived2>
            ::type::get_const_instance()
        ==
        * boost::serialization::type_info_implementation<polymorphic_derived2>
            ::type::get_const_instance().get_derived_extended_type_info(*rd21));
    delete rb2;
    delete rd21;
}

int main( int /* argc */, char* /* argv */[] )
{
    save_exported("testfile");
    load_exported("testfile");
    return EXIT_SUCCESS;
}

// EOF
