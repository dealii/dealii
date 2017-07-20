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

#ifndef dealii__pattern_tools_h
#define dealii__pattern_tools_h

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/base/std_cxx14/memory.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/utilities.h>

#include <boost/core/demangle.hpp>

#include <array>
#include <deque>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * Namespace for a few classes and functions that act on values and patterns,
 * and allow to convert from non elementary types to strings and vice versa.
 *
 * A typical usage of these tools is in the following example:
 *
 * @code
 * typedef std::vector<unsigned int> T;
 *
 * T vec(3);
 * vec[0] = 1;
 * vec[1] = 3;
 * vec[2] = 5;
 *
 * auto pattern = PatternTools::Convert<T>::to_pattern();
 *
 * std::cout << pattern->description() << std::endl;
 * // [List of <[Integer]> of length 0...4294967295 (inclusive)]
 *
 * auto s = PatternTools::Convert<T>::to_string(vec);
 * std::cout << s << std::endl;
 * // 1, 2, 3
 *
 * auto vec = PatternTools::Convert<T>::to_value("2,3,4,5");
 * // now vec has size 4, and contains the elements 2,3,4,5
 *
 * std::cout << internal::RankInfo<T>::list_rank << std::endl; // Outputs 1
 * std::cout << internal::RankInfo<T>::map_rank  << std::endl; // Outputs 0
 * @endcode
 *
 * Convert<T> is used by the function PatternTools::add_parameter() in this
 * namespace. Internally it uses the internal::RankInfo<T> class to decide how
 * many different separators are required to convert the given type to a string.
 *
 * For example, to write vectors of vectors, the default is to use "," for the
 * first (inner) separator, and ";" for the second (outer) separator, i.e.
 *
 * @code
 * std::vector<std::vector<unsigned int>> vec;
 * vec = Convert<decltype(vec)>::to_value("1,2,3 ; 4,5,6");
 *
 * s = convert<decltype(vec[0])>::to_string(vec[0]);
 * // s now contains the string "1,2,3"
 * @endcode
 *
 * Separators for Patterns::List and Patterns::Map compatible types are
 * selected according to the
 * rank of the list and map objects, using the arrays
 * PatternTools::internal::default_list_separator and
 * PatternTools::internal::default_map_separator.
 *
 * They are currently set to:
 *
 * @code
 * default_list_separator{{","  ,  ";"  ,  "|"  ,   "%"}};
 * default_map_separator {{":"  ,  "="  ,  "@"  ,   "#"}};
 * @endcode
 *
 * When one needs a mixture of Patterns::List and Patterns::Map types, their
 * RankInfo is computed by taking the maximum of the vector_rank of the Key and
 * of the Value type, so that, for example, it is possible to have the following
 * @code
 * ... // Build compare class
 * std::map<std::vector<unsigned int>, std::vector<double>, compare> map;
 *
 * map = convert<decltype(map)>::to_value("1,2,3 : 5.0,6.0,7.0  ; 8,9,10 :
 * 11.0,12.0,13.0");
 *
 * @endcode
 *
 * Some non elementary types are supported, like Point(), or
 * std::complex<double>. If you wish to support more types, you have to
 * specialize the Convert struct as well as the RankInfo struct.
 *
 * @ingroup input
 * @author Luca Heltai, 2017
 */
namespace PatternTools
{
  /**
   * Converter class. This class is used to generate strings and Patterns
   * associated to the given type, and to convert from a string to the given
   * type and viceversa.
   *
   * The second template parameter is used internally to allow for advanced
   * SFINAE (substitution failure is not an error) tricks used to specialise
   * this class for arbitrary STL containers and maps.
   *
   * @author Luca Heltai, 2017
   */
  template <class T, class Enable = void>
  struct Convert
  {

    /**
     * Return a std::unique_ptr to a Pattern that can be used to interpret a
     * string as the type of the template argument, and the other way around.
     *
     * While the current function (in the general Convert template) is deleted,
     * it is implemented and available in the specializations of the Convert
     * class template for particular kinds of template arguments @p T.
     */
    static std::unique_ptr<Patterns::PatternBase> to_pattern() = delete;

    /**
     * Return a string containing a textual version of the variable s. Use the
     * pattern passed to perform the conversion, or create and use a default
     * one.
     *
     * While the current function (in the general Convert template) is deleted,
     * it is implemented and available in the specializations of the Convert
     * class template for particular kinds of template arguments @p T.
     */
    static std::string to_string(const T &s,
                                 const std::unique_ptr<Patterns::PatternBase>
                                 &p = Convert<T>::to_pattern()) = delete;

    /**
     * Convert a string to a value, using the given pattern. Use the pattern
     * passed to perform the conversion, or create and use a default one.
     *
     * While the current function (in the general Convert template) is deleted,
     * it is implemented and available in the specializations of the Convert
     * class template for particular kinds of template arguments @p T.
     */
    static T to_value(const std::string &s,
                      const std::unique_ptr<Patterns::PatternBase> &p =
                        Convert<T>::to_pattern()) = delete;
  };

  /**
   * Declare a new entry in @p prm with name @p entry, set its default value
   * to the content of the variable @p parameter, and create an action
   * that will fill @p parameter with updated values when a file is parsed,
   * or the entry is set to a new value.
   *
   * By default, the pattern to use is obtained by calling the function
   * PatternTools::Convert<T>::to_pattern(), but a custom one can be used.
   */
  template <class ParameterType>
  void add_parameter(const std::string &entry,
                     ParameterType &parameter,
                     ParameterHandler &prm,
                     const std::string &documentation = std::string(),
                     const Patterns::PatternBase &pattern =
                       *Convert<ParameterType>::to_pattern());

  /**
   * @addtogroup Exceptions
   * @{
   */

  /**
   * Exception.
   */
  DeclException2(ExcNoMatch,
                 std::string,
                 Patterns::PatternBase &,
                 << "The string " << arg1 << " does not match the pattern \""
                 << arg2.description()
                 << "\"");
  //@}
}

// ---------------------- inline and template functions --------------------
namespace PatternTools
{
  namespace internal
  {
    /**
     * Store information about the rank types of the given class.
     *
     * A class has Rank equal to the number of different separators
     * that are required to uniquely identify its element(s) in a string.
     *
     * This class is used to detect wether the class T is compatible
     * with a Patterns::List pattern or with a Patterns::Map pattern.
     *
     * Objects like Point() or std::complex<double> are vector-likes, and
     * have vector_rank 1. Elementary types, like `int`, `unsigned int`,
     * `double`, etc. have vector_rank 0. `std::vector`, `std::list` and in
     * general containers have rank equal to 1 + vector_rank of the contained
     * type. Similarly for map types.
     *
     * A class with list_rank::value = 0 is either elementary or a
     * map. A class with map_rank::value = 0 is either a List compatible
     * class, or an elementary type.
     *
     * Elementary types are not compatible with Patterns::List, but non
     * elementary types, like Point(), or std::complex<double>, are compatible
     * with the List type. Adding more compatible types is a matter of adding a
     * specialization of this struct for the given type.
     *
     * @author Luca Heltai, 2017
     */
    template <class T, class Enable = void>
    struct RankInfo
    {
      static constexpr int list_rank = 0;
      static constexpr int map_rank = 0;
    };
  }

  // Arithmetic types
  template <class T>
  struct Convert<T, typename std::enable_if<std::is_arithmetic<T>::value>::type>
  {

    static std::unique_ptr<Patterns::PatternBase> to_pattern()
    {
      if (std::is_integral<T>::value)
        return std_cxx14::make_unique<Patterns::Integer>(
                 std::numeric_limits<T>::min(), std::numeric_limits<T>::max());
      else if (std::is_floating_point<T>::value)
        return std_cxx14::make_unique<Patterns::Double>(
                 -std::numeric_limits<T>::max(), std::numeric_limits<T>::max());
    }

    static std::string to_string(const T &value,
                                 const std::unique_ptr<Patterns::PatternBase>
                                 &p = Convert<T>::to_pattern())
    {
      std::string str;
      if (std::is_same<T, unsigned char>() || std::is_same<T, char>())
        str = std::to_string((int)value);
      else
        str = std::to_string(value);
      AssertThrow(p->match(str), ExcNoMatch(str, *p));
      return str;
    }

    static T to_value(const std::string &s,
                      const std::unique_ptr<Patterns::PatternBase> &p =
                        Convert<T>::to_pattern())
    {
      AssertThrow(p->match(s), ExcNoMatch(s, *p));
      std::istringstream is(s);
      T value;
      if (std::is_same<T, unsigned char>::value || std::is_same<T, char>::value)
        {
          int i;
          is >> i;
          value = i;
        }
      else
        is >> value;

      // If someone passes "123 abc" to the function, the method yelds an
      // integer 123 alright, but the space terminates the read from the string
      // although there is more to come. This case, however, is checked for in
      // the call p->match(s) at the beginning of this function, and would
      // throw earlier. Here it is safe to assume that if we didn't fail the
      // conversion with the operator >>, then we are good to go.
      AssertThrow(!is.fail(),
                  ExcMessage("Failed to convert from \"" + s +
                             "\" to the type \"" +
                             boost::core::demangle(typeid(T).name()) + "\""));
      return value;
    }
  };

  namespace internal
  {
    const std::array<std::string, 4> default_list_separator {{","  ,  ";"  ,  "|"  ,   "%"}};
    const std::array<std::string, 4> default_map_separator {{":"  ,  "="  ,  "@"  ,   "#"}};

    //specialize a type for all of the STL containers and maps
    template <typename T>       struct is_stl_container : std::false_type {};
    template <typename T, std::size_t N> struct is_stl_container<std::array    <T,N>>     : std::true_type {};
    template <typename... Args> struct is_stl_container<std::vector            <Args...>> : std::true_type {};
    template <typename... Args> struct is_stl_container<std::deque             <Args...>> : std::true_type {};
    template <typename... Args> struct is_stl_container<std::list              <Args...>> : std::true_type {};
    template <typename... Args> struct is_stl_container<std::set               <Args...>> : std::true_type {};
    template <typename... Args> struct is_stl_container<std::multiset          <Args...>> : std::true_type {};
    template <typename... Args> struct is_stl_container<std::unordered_set     <Args...>> : std::true_type {};
    template <typename... Args> struct is_stl_container<std::unordered_multiset<Args...>> : std::true_type {};

    template <typename T>       struct is_stl_map : std::false_type {};
    template <typename... Args> struct is_stl_map<std::map               <Args...>> : std::true_type {};
    template <typename... Args> struct is_stl_map<std::multimap          <Args...>> : std::true_type {};
    template <typename... Args> struct is_stl_map<std::unordered_map     <Args...>> : std::true_type {};
    template <typename... Args> struct is_stl_map<std::unordered_multimap<Args...>> : std::true_type {};
  }

  // type trait to use the implementation type traits as well as decay the type
  template <typename T>
  struct is_stl_container
  {
    static constexpr bool const value =
      internal::is_stl_container<typename std::decay<T>::type>::value;
  };

  template <typename T>
  struct is_stl_map
  {
    static constexpr bool const value =
      internal::is_stl_map<typename std::decay<T>::type>::value;
  };

  namespace internal
  {
    // Rank of vector types
    template <class T>
    struct RankInfo<T,
      typename std::enable_if<is_stl_container<T>::value>::type>
    {
      static constexpr int list_rank =
        RankInfo<typename T::value_type>::list_rank + 1;
      static constexpr int map_rank =
        RankInfo<typename T::value_type>::map_rank;
    };

    // Rank of map types
    template <class T>
    struct RankInfo<T, typename std::enable_if<is_stl_map<T>::value>::type>
    {
      static constexpr int list_rank =
        std::max(internal::RankInfo<typename T::key_type>::list_rank,
                 RankInfo<typename T::mapped_type>::list_rank) +
        1;
      static constexpr int map_rank =
        std::max(internal::RankInfo<typename T::key_type>::map_rank,
                 RankInfo<typename T::mapped_type>::map_rank) +
        1;
    };

    // Rank of Tensor types
    template <int rank, int dim, class Number>
    struct RankInfo<Tensor<rank, dim, Number>>
    {
      static constexpr int list_rank = rank + RankInfo<Number>::list_rank;
      static constexpr int map_rank = RankInfo<Number>::map_rank;
    };

    template <int dim, class Number>
    struct RankInfo<Point<dim, Number>> : RankInfo<Tensor<1, dim, Number>>
    {
    };

    // Rank of complex types
    template <class Number>
    struct RankInfo<std::complex<Number>>
    {
      static constexpr int list_rank = RankInfo<Number>::list_rank + 1;
      static constexpr int map_rank = RankInfo<Number>::map_rank;
    };

    template <class Key, class Value>
    struct RankInfo<std::pair<Key,Value>>
    {
      static constexpr int list_rank = std::max(RankInfo<Key>::list_rank, RankInfo<Value>::list_rank);
      static constexpr int map_rank = std::max(RankInfo<Key>::map_rank, RankInfo<Value>::map_rank)+1;
    };
  }

  // stl containers
  template <class T>
  struct Convert<T, typename std::enable_if<is_stl_container<T>::value>::type>
  {
    static std::unique_ptr<Patterns::PatternBase> to_pattern()
    {
      static_assert(internal::RankInfo<T>::list_rank > 0,
                    "Cannot use this class for non List-compatible types.");
      return std_cxx14::make_unique<Patterns::List>(
               *Convert<typename T::value_type>::to_pattern(),
               0,
               std::numeric_limits<unsigned int>::max(),
               internal::default_list_separator[internal::RankInfo<T>::list_rank - 1]);
    }

    static std::string to_string(const T &t,
                                 const std::unique_ptr<Patterns::PatternBase>
                                 &pattern = Convert<T>::to_pattern())
    {

      auto p = dynamic_cast<const Patterns::List *>(pattern.get());
      AssertThrow(
        p,
        ExcMessage(
          "I need a List pattern to convert a string to a List type."));
      auto base_p = p->get_base_pattern().clone();
      std::vector<std::string> vec(t.size());

      unsigned int i = 0;
      for (const auto &ti : t)
        vec[i++] = Convert<typename T::value_type>::to_string(ti, base_p);

      std::string s;
      if (vec.size() > 0)
        s = vec[0];
      for (unsigned int i = 1; i < vec.size(); ++i)
        s += p->get_separator() + " " + vec[i];

      AssertThrow(pattern->match(s), ExcNoMatch(s, *p));
      return s;
    }

    static T to_value(const std::string &s,
                      const std::unique_ptr<Patterns::PatternBase> &pattern =
                        Convert<T>::to_pattern())
    {

      AssertThrow(pattern->match(s), ExcNoMatch(s, *pattern));

      auto p = dynamic_cast<const Patterns::List *>(pattern.get());
      AssertThrow(
        p,
        ExcMessage(
          "I need a List pattern to convert a string to a List type."));

      auto base_p = p->get_base_pattern().clone();
      T t;

      auto v = Utilities::split_string_list(s, p->get_separator());
      for (const auto &str : v)
        t.insert(t.end(),
                 Convert<typename T::value_type>::to_value(str, base_p));

      return t;
    }
  };

  // stl maps
  template <class T>
  struct Convert<T, typename std::enable_if<is_stl_map<T>::value>::type>
  {
    static std::unique_ptr<Patterns::PatternBase> to_pattern()
    {
      static_assert(internal::RankInfo<T>::list_rank > 0,
                    "Cannot use this class for non List-compatible types.");
      static_assert(internal::RankInfo<T>::map_rank > 0,
                    "Cannot use this class for non Map-compatible types.");
      return std_cxx14::make_unique<Patterns::Map>(
               *Convert<typename T::key_type>::to_pattern(),
               *Convert<typename T::mapped_type>::to_pattern(),
               0,
               std::numeric_limits<unsigned int>::max(),
               internal::default_list_separator[internal::RankInfo<T>::list_rank - 1],
               internal::default_map_separator[internal::RankInfo<T>::map_rank - 1]);
    }

    static std::string to_string(const T &t,
                                 const std::unique_ptr<Patterns::PatternBase>
                                 &pattern = Convert<T>::to_pattern())
    {
      auto p = dynamic_cast<const Patterns::Map *>(pattern.get());
      AssertThrow(
        p,
        ExcMessage("I need a Map pattern to convert a string to a List type."));
      auto key_p = p->get_key_pattern().clone();
      auto val_p = p->get_value_pattern().clone();
      std::vector<std::string> vec(t.size());

      unsigned int i = 0;
      for (const auto &ti : t)
        vec[i++] =
          Convert<typename T::key_type>::to_string(ti.first, key_p) +
          p->get_key_value_separator() +
          Convert<typename T::mapped_type>::to_string(ti.second, val_p);

      std::string s;
      if (vec.size() > 0)
        s = vec[0];
      for (unsigned int i = 1; i < vec.size(); ++i)
        s += p->get_separator() + " " + vec[i];

      AssertThrow(p->match(s), ExcNoMatch(s, *p));
      return s;
    }

    static T to_value(const std::string &s,
                      const std::unique_ptr<Patterns::PatternBase> &pattern =
                        Convert<T>::to_pattern())
    {

      AssertThrow(pattern->match(s), ExcNoMatch(s, *pattern));

      auto p = dynamic_cast<const Patterns::Map *>(pattern.get());
      AssertThrow(
        p,
        ExcMessage("I need a Map pattern to convert a string to a List type."));

      auto key_p = p->get_key_pattern().clone();
      auto val_p = p->get_value_pattern().clone();
      T t;

      auto v = Utilities::split_string_list(s, p->get_separator());
      for (const auto &str : v)
        {
          auto key_val =
            Utilities::split_string_list(str, p->get_key_value_separator());
          AssertDimension(key_val.size(), 2);
          t.insert(std::make_pair(
                     Convert<typename T::key_type>::to_value(key_val[0], key_p),
                     Convert<typename T::mapped_type>::to_value(key_val[1])));
        }

      return t;
    }
  };

  // Tensors
  template <int rank, int dim, class Number>
  struct Convert<Tensor<rank, dim, Number>>
  {
    typedef Tensor<rank, dim, Number> T;
    static std::unique_ptr<Patterns::PatternBase> to_pattern()
    {
      static_assert(internal::RankInfo<T>::list_rank > 0,
                    "Cannot use this class for non List-compatible types.");
      return std_cxx14::make_unique<Patterns::List>(
               *Convert<typename T::value_type>::to_pattern(),
               dim,
               dim,
               internal::default_list_separator[internal::RankInfo<T>::list_rank - 1]);
    }

    static std::string to_string(const T &t,
                                 const std::unique_ptr<Patterns::PatternBase>
                                 &pattern = Convert<T>::to_pattern())
    {

      auto p = dynamic_cast<const Patterns::List *>(pattern.get());
      AssertThrow(
        p,
        ExcMessage(
          "I need a List pattern to convert a string to a List type."));
      auto base_p = p->get_base_pattern().clone();
      std::vector<std::string> vec(dim);

      for (unsigned int i = 0; i < dim; ++i)
        vec[i] = Convert<typename T::value_type>::to_string(t[i], base_p);

      std::string s;
      if (vec.size() > 0)
        s = vec[0];
      for (unsigned int i = 1; i < vec.size(); ++i)
        s += p->get_separator() + " " + vec[i];

      AssertThrow(p->match(s), ExcNoMatch(s, *p));
      return s;
    }

    static T to_value(const std::string &s,
                      const std::unique_ptr<Patterns::PatternBase> &pattern =
                        Convert<T>::to_pattern())
    {

      AssertThrow(pattern->match(s), ExcNoMatch(s, *pattern));

      auto p = dynamic_cast<const Patterns::List *>(pattern.get());
      AssertThrow(
        p,
        ExcMessage(
          "I need a List pattern to convert a string to a List type."));

      auto base_p = p->get_base_pattern().clone();
      T t;

      auto v = Utilities::split_string_list(s, p->get_separator());
      unsigned int i = 0;
      for (const auto &str : v)
        t[i++] = Convert<typename T::value_type>::to_value(str, base_p);

      return t;
    }
  };

  // Points
  template <int dim, class Number>
  struct Convert<Point<dim, Number>>
  {

    typedef Point<dim, Number> T;

    static std::unique_ptr<Patterns::PatternBase> to_pattern()
    {
      return Convert<Tensor<1, dim, Number>>::to_pattern();
    }

    static std::string to_string(const T &t,
                                 const std::unique_ptr<Patterns::PatternBase>
                                 &pattern = Convert<T>::to_pattern())
    {
      return Convert<Tensor<1, dim, Number>>::to_string(
               Tensor<1, dim, Number>(t), pattern);
    }

    static T to_value(const std::string &s,
                      const std::unique_ptr<Patterns::PatternBase> &pattern =
                        Convert<T>::to_pattern())
    {
      return T(Convert<Tensor<1, dim, Number>>::to_value(s, pattern));
    }
  };

  // Complex numbers
  template <class Number>
  struct Convert<std::complex<Number>>
  {
    typedef std::complex<Number> T;

    static std::unique_ptr<Patterns::PatternBase> to_pattern()
    {
      static_assert(internal::RankInfo<T>::list_rank > 0,
                    "Cannot use this class for non List-compatible types.");
      return std_cxx14::make_unique<Patterns::List>(
               *Convert<typename T::value_type>::to_pattern(),
               2,
               2,
               internal::default_list_separator[internal::RankInfo<T>::list_rank - 1]);
    }

    static std::string to_string(const T &t,
                                 const std::unique_ptr<Patterns::PatternBase>
                                 &pattern = Convert<T>::to_pattern())
    {

      auto p = dynamic_cast<const Patterns::List *>(pattern.get());
      AssertThrow(
        p,
        ExcMessage(
          "I need a List pattern to convert a string to a List type."));
      auto base_p = p->get_base_pattern().clone();
      std::string s =
        Convert<typename T::value_type>::to_string(t.real(), base_p) +
        p->get_separator() + " " +
        Convert<typename T::value_type>::to_string(t.imag(), base_p);

      AssertThrow(pattern->match(s), ExcNoMatch(s, *p));
      return s;
    }

    /**
     * Convert a string to a value, using the given pattern, or a default one.
     */
    static T to_value(const std::string &s,
                      const std::unique_ptr<Patterns::PatternBase> &pattern =
                        Convert<T>::to_pattern())
    {

      AssertThrow(pattern->match(s), ExcNoMatch(s, *pattern));

      auto p = dynamic_cast<const Patterns::List *>(pattern.get());
      AssertThrow(
        p,
        ExcMessage(
          "I need a List pattern to convert a string to a List type."));

      auto base_p = p->get_base_pattern().clone();

      auto v = Utilities::split_string_list(s, p->get_separator());
      AssertDimension(v.size(), 2);
      T t(Convert<typename T::value_type>::to_value(v[0], base_p),
          Convert<typename T::value_type>::to_value(v[1], base_p));
      return t;
    }
  };

  // Strings
  template <>
  struct Convert<std::string>
  {
    typedef std::string T;

    static std::unique_ptr<Patterns::PatternBase> to_pattern()
    {
      return std_cxx14::make_unique<Patterns::Anything>();
    }

    static std::string to_string(const T &t,
                                 const std::unique_ptr<Patterns::PatternBase>
                                 &pattern = Convert<T>::to_pattern())
    {
      AssertThrow(pattern->match(t), ExcNoMatch(t, *pattern));
      return t;
    }

    static T to_value(const std::string &s,
                      const std::unique_ptr<Patterns::PatternBase> &pattern =
                        Convert<T>::to_pattern())
    {
      AssertThrow(pattern->match(s), ExcNoMatch(s, *pattern));
      return s;
    }
  };


  // Pairs
  template <class Key, class Value>
  struct Convert<std::pair<Key,Value>>
  {
    typedef std::pair<Key,Value> T;

    static std::unique_ptr<Patterns::PatternBase> to_pattern()
    {
      static_assert(internal::RankInfo<T>::map_rank > 0,
                    "Cannot use this class for non Map-compatible types.");
      return std_cxx14::make_unique<Patterns::Map>(
               *Convert<Key>::to_pattern(),
               *Convert<Value>::to_pattern(),
               1, 1,
               // We keep the same list separator of the previous level, as this is
               // a map with only 1 possible entry
               internal::default_list_separator[internal::RankInfo<T>::list_rank],
               internal::default_map_separator[internal::RankInfo<T>::map_rank - 1]);
    }

    static std::string to_string(const T &t,
                                 const std::unique_ptr<Patterns::PatternBase>
                                 &pattern = Convert<T>::to_pattern())
    {
      std::unordered_map<Key, Value> m;
      m.insert(t);
      std:: string s = Convert<decltype(m)>::to_string(m, pattern);
      AssertThrow(pattern->match(s), ExcNoMatch(s, *pattern));
      return s;
    }

    static T to_value(const std::string &s,
                      const std::unique_ptr<Patterns::PatternBase> &pattern =
                        Convert<T>::to_pattern())
    {
      std::unordered_map<Key, Value> m;
      m = Convert<decltype(m)>::to_value(s, pattern);
      return *m.begin();
    }
  };

  template <class ParameterType>
  void add_parameter(const std::string &entry,
                     ParameterType &parameter,
                     ParameterHandler &prm,
                     const std::string &documentation,
                     const Patterns::PatternBase &pattern)
  {

    static_assert(std::is_const<ParameterType>::value == false,
                  "You tried to add a parameter using a type "
                  "that is const. Use a non-const type.");

    prm.declare_entry(entry,
                      PatternTools::Convert<ParameterType>::to_string(
                        parameter, pattern.clone()),
                      pattern,
                      documentation);

    auto action = [&](const std::string &val)
    {
      parameter =
        PatternTools::Convert<ParameterType>::to_value(val, pattern.clone());
    };
    prm.add_action(entry, action);
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
