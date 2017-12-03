#ifndef POLYMORPHIC_BASE_HPP
#define POLYMORPHIC_BASE_HPP

#include <boost/config.hpp>

#include <boost/serialization/access.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/type_info_implementation.hpp>
#include <boost/serialization/extended_type_info_no_rtti.hpp>

class BOOST_SYMBOL_VISIBLE polymorphic_base
{
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(
        Archive & /* ar */, 
        const unsigned int /* file_version */
    ){}
public:
    // note that since this class uses the "no_rtti"
    // extended_type_info implementation, it MUST
    // implement this function
    virtual const char * get_key() const = 0;
    virtual ~polymorphic_base(){};
};

BOOST_SERIALIZATION_ASSUME_ABSTRACT(polymorphic_base)

// the no_rtti system requires this !!!
BOOST_CLASS_EXPORT_KEY(polymorphic_base)

BOOST_CLASS_TYPE_INFO(
    polymorphic_base,
    boost::serialization::extended_type_info_no_rtti<polymorphic_base>
)

#endif // POLYMORPHIC_BASE_HPP
