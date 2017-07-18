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
   * Store information about the rank types of the given class.
   *
   * A class has Rank equal to the number of different separators
   * that are required to specify its element(s) in a string.
   *
   * This class is used to detect wether the class T is compatible
   * with a Patterns::List pattern or with a Patterns::Map pattern, and
   * to choose among a default list of separators for complex types,
   * like vectors of vectors.
   *
   * Objects like Point<dim> or std::complex<double> are vector-likes, and
   * have vector_rank 1. Elementary types, like `int`, `unsigned int`, `double`, etc.
   * have vector_rank 0. `std::vector`, `std::list` and in general containers have rank
   * equal to 1 + vector_rank of the contained type. Similarly for map types.
   *
   * A class with vector_rank_type::value = 0 is either elementary or a
   * map. A class with map_rank_type::value = 0 is either a List compatible
   * class, or an elementary type.
   *
   * Elementary types are not compatible with Patterns::List, but
   * non elementary types, like Point<dim>(), or std::complex<double>, are
   * compatible with the List type. Adding more compatible types is a matter
   * of adding a specialization of this struct for the given type.
   */
  template<class T>
  struct RankInfo
  {
    typedef std::integral_constant<int, 0>::type vector_rank_type;
  };


  /**
   * Return the default list separator for an object with the given vector rank.
   *
   * Objects like Point<dim> or std::complex<double> are vector-likes, and
   * have rank 1. Elementary types, like `int`, `unsigned int`, `double`, etc.
   * have rank 0. `std::vector`, `std::list` and in general containers have rank
   * equal to 1 + rank of the contained type.
   *
   * This function helps in constructing patterns for non elementary types,
   * like for example std::vector<std::vector<std::complex<double>>>.
   */
  std::string default_list_separator(unsigned int);



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

  /**
   * @addtogroup Exceptions
   * @{
   */

  /**
   * Exception.
   */
  DeclException2 (ExcNoMatch,
                  std::string, PatternBase &,
                  << "The string " << arg1 << " does not match the pattern \""
                  << arg2.description() << "\"");
  //@}

}

// ---------------------- inline and template functions --------------------

namespace PatternsTools
{
  template<class T>
  struct Convert<T, typename std::enable_if<std::is_arithmetic<T>::value>::type>
  {

    static std::unique_ptr<PatternBase> to_pattern()
    {
      if(std::is_integral<T>::value)
        return std_cxx14::make_unique<Integer>(std::numeric_limits<T>::min(), std::numeric_limits<T>::max());
      else if(std::is_floating_point<T>::value)
        return std_cxx14::make_unique<Double>(std::numeric_limits<T>::min(), std::numeric_limits<T>::max());
      else {
        AssertThrow(false, ExcNotImplemented());
        return std::unique_ptr<PatternBase>();
        }
    }

    static std::string to_string(const T &value, std::unique_ptr<PatternBase> p = Convert<T>::to_pattern())
    {
      std::stringstream str;
      str << value;
      AssertThrow(p->match(str.str()), ExcNoMatch(str.str(), *p));
      return str.str();
    }

    static T to_value(const std::string &s,
                      std::unique_ptr<PatternBase> p = Convert<T>::to_pattern())
    {
      AssertThrow(p->match(s), ExcNoMatch(s, *p));
      std::istringstream is(s);
      T i;
      is >> i;
      AssertThrow(!is.fail(), ExcMessage("Failed to convert from \"" + s + "\" to the type \""
                                         +boost::core::demangle(typeid(T).name()) + "\""));
      return i;
    }
  };

  // Rank of vector types
  template<template <class T1, class A1> class Container, class T, class Allocator>
  struct RankInfo<Container<T,Allocator>>
  {
    typedef typename std::integral_constant<int, RankInfo<T>::vector_rank_type::value+1>::type vector_rank_type;
  };





  template<template <class T1, class A1> class Container, class T, class Allocator>
  struct Convert<Container<T,Allocator>>
  {
    static std::unique_ptr<PatternBase> to_pattern()
    {
      return std_cxx14::make_unique<List>(*Convert<Container<T,Allocator>>::to_pattern(),
                                          0, std::numeric_limits<int>::max,
                                          default_list_separator(RankInfo<Container<T,Allocator>>::vector_rank_type::value));
    }

    static std::string to_string(const Container<T,Allocator> &t,
                                 std::unique_ptr<PatternBase> pattern = Convert<Container<T,Allocator>>::to_pattern())
    {

      auto p = dynamic_cast<const Patterns::List *>(pattern.get());
      AssertThrow(p, ExcMessage("I need a List pattern to convert a string to a List type."));
      auto base_p = p->get_base_pattern().clone();
      std::vector<std::string> vec(t.size());

      unsigned int i=0;
      for (auto &ti : t)
        vec[i++] = Convert<T>::to_string(ti, base_p);

      std::string s;
      if (vec.size() > 0)
        s = vec[0];
      for (unsigned int i=1; i<vec.size(); ++i)
        s += p->get_separator() + " " + vec[i];

      Assert(p->match(s), ExcMessage("No match for " + s +
                                     " with pattern " + p->description()));
      return s;
    }

    /**
     * Convert a string to a value, using the given pattern, or a default one.
     */
    static Container<T,Allocator>  to_value(const std::string &s,
                                            std::unique_ptr<PatternBase> pattern = Convert<Container<T,Allocator>>::to_pattern())
    {

      AssertThrow(pattern->match(s), ExcMessage("No match for " + s +
                                                " using pattern " + pattern->description()));

      auto p = dynamic_cast<const Patterns::List *>(pattern.get());
      AssertThrow(p, ExcMessage("I need a List pattern to convert a string to a List type."));

      auto base_p = p->get_base_pattern().clone();
      Container<T,Allocator> t;

      auto v = Utilities::split_string_list(s,p->get_separator());
      for (auto str : v)
        t.insert(t.end(), Convert<T>::to_value(str, base_p));

      return t;
    }
  };



}

DEAL_II_NAMESPACE_CLOSE

#endif
