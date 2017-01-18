/*
 Copyright (c) 2014 Glen Joseph Fernandes
 glenfe at live dot com

 Distributed under the Boost Software License,
 Version 1.0. (See accompanying file LICENSE_1_0.txt
 or copy at http://boost.org/LICENSE_1_0.txt)
*/
#ifndef BOOST_ALIGN_DETAIL_PADDED_HPP
#define BOOST_ALIGN_DETAIL_PADDED_HPP

namespace boost {
    namespace alignment {
        namespace detail {
            template<class T>
            struct padded {
                char unit;
                T object;
            };
        }
    }
}

#endif
