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
#include <deal.II/base/utilities.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/point.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/std_cxx14/memory.h>

#include <boost/core/demangle.hpp>

#include <map>
#include <vector>
#include <string>
#include <memory>
#include <sstream>
#include <deque>
#include <forward_list>
#include <list>
#include <map>
#include <queue>
#include <set>
#include <stack>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <type_traits>


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
  template<class T, class Enable=void>
  struct RankInfo
  {
    typedef std::integral_constant<int, 0>::type vector_rank_type;
    typedef std::integral_constant<int, 0>::type map_rank_type;
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
   * Return the default map separator for an object with the given map rank.
   *
   * This function helps in constructing patterns for non elementary types,
   * like for example std::map<unsigned int, std::map<unsigned int, double>>
   */
  std::string default_map_separator(unsigned int);


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
                                 const std::unique_ptr<PatternBase> &p = Convert<T>::to_pattern()) = delete;


    /**
     * Convert a string to a value, using the given pattern, or a default one.
     */
    static T to_value(const std::string &s,
                      const std::unique_ptr<PatternBase> &p = Convert<T>::to_pattern()) = delete;
  };

  /**
   * Declare a new entry in @p prm with name @p entry, set its default value
   * to the content of the variable @p parameter, and create an action
   * that will fill @p parameter with updated values when a file is parsed,
   * or the entry is set to a new value.
   *
   * By default, the pattern to use is obtained by calling the function
   * PatternsTools::Convert<T>::to_pattern(), but a custom one can be used.
   */
  template <class ParameterType>
  void add_parameter(const std::string           &entry,
                     ParameterType               &parameter,
                     ParameterHandler            &prm,
                     const std::string           &documentation = std::string(),
                     const Patterns::PatternBase &pattern = *Convert<ParameterType>::to_pattern());

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
      if (std::is_integral<T>::value)
        return std_cxx14::make_unique<Integer>(std::numeric_limits<T>::min(), std::numeric_limits<T>::max());
      else if (std::is_floating_point<T>::value)
        return std_cxx14::make_unique<Double>(std::numeric_limits<T>::min(), std::numeric_limits<T>::max());
      else
        {
          AssertThrow(false, ExcNotImplemented());
          return std::unique_ptr<PatternBase>();
        }
    }

    static std::string to_string(const T &value, const std::unique_ptr<PatternBase> &p = Convert<T>::to_pattern())
    {
      std::stringstream str;
      str << value;
      AssertThrow(p->match(str.str()), ExcNoMatch(str.str(), *p));
      return str.str();
    }

    static T to_value(const std::string &s,
                      const std::unique_ptr<PatternBase> &p = Convert<T>::to_pattern())
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

  //specialize a type for all of the STL containers and maps
  namespace internal
  {
    template <typename T>       struct is_stl_container:std::false_type {};
    template <typename T, std::size_t N> struct is_stl_container<std::array    <T,N>>    :std::true_type {};
    template <typename... Args> struct is_stl_container<std::vector            <Args...>>:std::true_type {};
    template <typename... Args> struct is_stl_container<std::deque             <Args...>>:std::true_type {};
    template <typename... Args> struct is_stl_container<std::list              <Args...>>:std::true_type {};
    template <typename... Args> struct is_stl_container<std::forward_list      <Args...>>:std::true_type {};
    template <typename... Args> struct is_stl_container<std::set               <Args...>>:std::true_type {};
    template <typename... Args> struct is_stl_container<std::multiset          <Args...>>:std::true_type {};
    template <typename... Args> struct is_stl_container<std::unordered_set     <Args...>>:std::true_type {};
    template <typename... Args> struct is_stl_container<std::unordered_multiset<Args...>>:std::true_type {};
    template <typename... Args> struct is_stl_container<std::stack             <Args...>>:std::true_type {};
    template <typename... Args> struct is_stl_container<std::queue             <Args...>>:std::true_type {};
    template <typename... Args> struct is_stl_container<std::priority_queue    <Args...>>:std::true_type {};

    template <typename T>       struct is_stl_map:std::false_type {};
    template <typename... Args> struct is_stl_map<std::map               <Args...>>:std::true_type {};
    template <typename... Args> struct is_stl_map<std::multimap          <Args...>>:std::true_type {};
    template <typename... Args> struct is_stl_map<std::unordered_map     <Args...>>:std::true_type {};
    template <typename... Args> struct is_stl_map<std::unordered_multimap<Args...>>:std::true_type {};
  }

  //type trait to use the implementation type traits as well as decay the type
  template <typename T> struct is_stl_container
  {
    static constexpr bool const value = internal::is_stl_container<std::decay_t<T>>::value;
  };


  template <typename T> struct is_stl_map
  {
    static constexpr bool const value = internal::is_stl_map<std::decay_t<T>>::value;
  };


  // Rank of vector types
  template<class T>
  struct RankInfo<T, typename std::enable_if<is_stl_container<T>::value>::type >
  {
    typedef typename std::integral_constant<int, RankInfo<typename T::value_type>::vector_rank_type::value+1>::type vector_rank_type;
    typedef typename std::integral_constant<int, 0>::type map_rank_type;
  };


  // Rank of vector types
  template<class T>
  struct RankInfo<T, typename std::enable_if<is_stl_map<T>::value>::type >
  {
    typedef typename std::integral_constant<int, std::max(RankInfo<typename T::key_type>::vector_rank_type::value,
                                                          RankInfo<typename T::value_type>::vector_rank_type::value)+1>::type vector_rank_type;
    typedef typename std::integral_constant<int, std::max(RankInfo<typename T::key_type>::map_rank_type::value,
                                                          RankInfo<typename T::value_type>::map_rank_type::value)+1>::type map_rank_type;
  };



  template<class T>
  struct Convert<T, typename std::enable_if<is_stl_container<T>::value>::type>
  {
    static std::unique_ptr<PatternBase> to_pattern()
    {
      return std_cxx14::make_unique<List>(*Convert<typename T::value_type>::to_pattern(),
                                          0, std::numeric_limits<unsigned int>::max(),
                                          default_list_separator(RankInfo<T>::vector_rank_type::value));
    }

    static std::string to_string(const T &t,
                                 const std::unique_ptr<PatternBase> &pattern = Convert<T>::to_pattern())
    {

      auto p = dynamic_cast<const Patterns::List *>(pattern.get());
      AssertThrow(p, ExcMessage("I need a List pattern to convert a string to a List type."));
      auto base_p = p->get_base_pattern().clone();
      std::vector<std::string> vec(t.size());

      unsigned int i=0;
      for (auto &ti : t)
        vec[i++] = Convert<typename T::value_type>::to_string(ti, base_p);

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
    static T to_value(const std::string &s,
                      const std::unique_ptr<PatternBase> &pattern = Convert<T>::to_pattern())
    {

      AssertThrow(pattern->match(s), ExcMessage("No match for " + s +
                                                " using pattern " + pattern->description()));

      auto p = dynamic_cast<const Patterns::List *>(pattern.get());
      AssertThrow(p, ExcMessage("I need a List pattern to convert a string to a List type."));

      auto base_p = p->get_base_pattern().clone();
      T t;

      auto v = Utilities::split_string_list(s,p->get_separator());
      for (auto str : v)
        t.insert(t.end(), Convert<typename T::value_type>::to_value(str, base_p));

      return t;
    }
  };


  template <class T>
  struct Convert<T, typename std::enable_if<is_stl_map<T>::value>::type>
  {
    static std::unique_ptr<PatternBase> to_pattern()
    {
      return std_cxx14::make_unique<Map>(*Convert<typename T::key_type>::to_pattern(),
                                         *Convert<typename T::mapped_type>::to_pattern(),
                                         0, std::numeric_limits<unsigned int>::max(),
                                         default_list_separator(RankInfo<T>::vector_rank_type::value),
                                         default_map_separator(RankInfo<T>::map_rank_type::value)
                                        );
    }

    static std::string to_string(const T &t,
                                 const std::unique_ptr<PatternBase> &pattern = Convert<T>::to_pattern())
    {

      auto p = dynamic_cast<const Patterns::Map *>(pattern.get());
      AssertThrow(p, ExcMessage("I need a Map pattern to convert a string to a List type."));
      auto key_p = p->get_key_pattern().clone();
      auto val_p = p->get_value_pattern().clone();
      std::vector<std::string> vec(t.size());

      unsigned int i=0;
      for (auto &ti : t)
        vec[i++] =
          Convert<typename T::key_type>::to_string(ti.first, key_p) +
          p->get_key_value_separator()+
          Convert<typename T::mapped_type>::to_string(ti.second, val_p);

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
    static T  to_value(const std::string &s,
                       const std::unique_ptr<PatternBase> &pattern = Convert<T>::to_pattern())
    {

      AssertThrow(pattern->match(s), ExcMessage("No match for " + s +
                                                " using pattern " + pattern->description()));

      auto p = dynamic_cast<const Patterns::Map *>(pattern.get());
      AssertThrow(p, ExcMessage("I need a Map pattern to convert a string to a List type."));

      auto key_p = p->get_key_pattern().clone();
      auto val_p = p->get_value_pattern().clone();
      T t;

      auto v = Utilities::split_string_list(s,p->get_separator());
      for (auto str : v)
        {
          auto key_val = Utilities::split_string_list(str, p->get_key_value_separator());
          AssertDimension(key_val.size(), 2);
          t[Convert<typename T::key_type>::to_value(key_val[0], key_p)] =
            Convert<typename T::mapped_type>::to_value(key_val[1]);
        }

      return t;
    }
  };

}

DEAL_II_NAMESPACE_CLOSE

#endif
