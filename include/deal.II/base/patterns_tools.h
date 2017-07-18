// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef dealii__patterns_tools_h
#define dealii__patterns_tools_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/point.h>
#include <deal.II/base/std_cxx14/memory.h>

#include <boost/core/demangle.hpp>
#include <boost/any.hpp>


#include <map>
#include <vector>
#include <string>
#include <memory>
#include <type_traits>
#include <sstream>

DEAL_II_NAMESPACE_OPEN


/**
 * Namespace for a few class and functions that act on values and patterns.
 *
 * @ingroup input
 */
namespace PatternsTools
{
  using namespace Patterns;
  /**
   * Converter class. This class is used to generate strings and Patterns associated to
   * the given type, and to convert from a string to the given type and viceversa.
   */
  template<class T, class Enable=void>
  struct Convert
  {

    /**
     * Return a std::unique_ptr to a Pattern that can be used to interpret a
     * string as the type of the template argument, and the other way around.
     */
    static std::unique_ptr<PatternBase> to_pattern() = delete;


    /**
     * Return a string containing a textual version of the variable s. Use the pattern
     * passed to perform the conversion, or create and use a default one.
     */
    static std::string to_string(const T &s,
                                 std::unique_ptr<PatternBase> p = Convert<T>::to_pattern()) = delete;


    /**
     * Convert a string to a value, using the given pattern, or a default one.
     */
    static T to_value(const std::string &s,
                      std::unique_ptr<PatternBase> p = Convert<T>::to_pattern()) = delete;
  };


}

// ---------------------- inline and template functions --------------------

namespace PatternsTools
{
  template<class T>
  struct Convert<T, typename std::enable_if<std::is_integral<T>::value>::type>
  {

    static std::unique_ptr<PatternBase> to_pattern()
    {
      return std_cxx14::make_unique<Integer>(std::numeric_limits<T>::min(), std::numeric_limits<T>::max());
    }

    static std::string to_string(const T &value, std::unique_ptr<PatternBase> p = Convert<T>::to_pattern())
    {
      std::stringstream str;
      str << value;
      AssertThrow(p->match(str.str()), ExcMessage("No match"));
      return str.str();
    }

    static T to_value(const std::string &s,
                      std::unique_ptr<PatternBase> p = Convert<T>::to_pattern())
    {
      AssertThrow(p->match(s), ExcMessage("No match"));
      std::istringstream is(s);
      T i;
      is >> i;
      AssertThrow(!is.fail(), ExcMessage("Failed to convert"));
      return i;
    }
  };
}

DEAL_II_NAMESPACE_CLOSE

#endif
