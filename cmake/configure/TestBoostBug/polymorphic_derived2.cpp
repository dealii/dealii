#include <boost/serialization/type_info_implementation.hpp>
#include <boost/serialization/extended_type_info_no_rtti.hpp>
#include <boost/serialization/export.hpp>

#include "polymorphic_derived2.hpp"

template<class Archive>
void polymorphic_derived2::serialize(
    Archive &ar, 
    const unsigned int /* file_version */
){
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(polymorphic_base);
}

// instantiate code for text archives
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

template
void polymorphic_derived2::serialize(
    boost::archive::text_oarchive & ar,
    const unsigned int version
);
template
void polymorphic_derived2::serialize(
    boost::archive::text_iarchive & ar,
    const unsigned int version
);

// instantiate code for polymorphic archives
#include <boost/archive/polymorphic_iarchive.hpp>
#include <boost/archive/polymorphic_oarchive.hpp>

template
void polymorphic_derived2::serialize(
    boost::archive::polymorphic_oarchive & ar,
    const unsigned int version
);
template
/*POLYMORPHIC_DERIVED2_DLL_DECL*/
void polymorphic_derived2::serialize(
    boost::archive::polymorphic_iarchive & ar,
    const unsigned int version
);

// note: export has to be AFTER #includes for all archive classes
BOOST_CLASS_EXPORT_IMPLEMENT(polymorphic_derived2)
