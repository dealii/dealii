// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// BOOST defines a number of macros that are, of course, not exported
// by the module partition that we use to wrap BOOST. Moreover, BOOST
// does not have individual files in which these are defined, and so
// we cannot selectively just get the macros -- we'd have to get the
// whole shebang, which of course defeats the purpose of wrapping all
// of BOOST into a module partition to begin with. Rather -- perhaps
// imprudently -- we have to repeat these macros here.

#ifndef dealii_boost_macros_h
#define dealii_boost_macros_h


// Taken from boost/serialization/split_member.hpp:
#ifndef BOOST_SERIALIZATION_SPLIT_MEMBER
#  define BOOST_SERIALIZATION_SPLIT_MEMBER()                       \
    template <class Archive>                                       \
    void serialize(Archive &ar, const unsigned int file_version)   \
    {                                                              \
      boost::serialization::split_member(ar, *this, file_version); \
    }
#endif

#endif
