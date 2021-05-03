// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


#include <deal.II/base/logstream.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/path_search.h>
#include <deal.II/base/patterns.h>
#include <deal.II/base/utilities.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include <boost/io/ios_state.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#undef BOOST_BIND_GLOBAL_PLACEHOLDERS
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>


DEAL_II_NAMESPACE_OPEN



// TODO[WB]: various functions here could be simplified by using namespace
// Utilities

namespace Patterns
{
  namespace internal
  {
    std::string
    escape(const std::string &input, const PatternBase::OutputStyle style)
    {
      switch (style)
        {
          case PatternBase::Machine:
          case PatternBase::Text:
            return input;
          case PatternBase::LaTeX:
            {
              std::string u;
              u.reserve(input.size());
              for (const auto c : input)
                {
                  switch (c)
                    {
                      case '#':
                      case '$':
                      case '%':
                      case '&':
                      case '_':
                      case '{':
                      case '}':
                        // simple escaping:
                        u.push_back('\\');
                        u.push_back(c);
                        break;

                      case '\\':
                        u.append("\\textbackslash{}");
                        break;

                      case '^':
                        u.append("\\^{}");
                        break;

                      case '~':
                        u.append("\\~{}");
                        break;

                      default:
                        // all other chars are just copied:
                        u.push_back(c);
                    }
                }
              return u;
            }
          default:
            Assert(false, ExcNotImplemented());
        }
      return "";
    }

  } // namespace internal


  namespace
  {
    /**
     * Read to the end of the stream and
     * return whether all there is is
     * whitespace or whether there are other
     * characters as well.
     */
    bool
    has_only_whitespace(std::istream &in)
    {
      while (in)
        {
          char c;

          // skip if we've reached the end of
          // the line
          if (!(in >> c))
            break;

          if ((c != ' ') && (c != '\t'))
            return false;
        }
      return true;
    }
  } // namespace



  std::unique_ptr<PatternBase>
  pattern_factory(const std::string &description)
  {
    std::unique_ptr<PatternBase> p;

    p = Integer::create(description);
    if (p != nullptr)
      return p;

    p = Double::create(description);
    if (p != nullptr)
      return p;

    p = Selection::create(description);
    if (p != nullptr)
      return p;

    p = List::create(description);
    if (p != nullptr)
      return p;

    p = Map::create(description);
    if (p != nullptr)
      return p;

    p = MultipleSelection::create(description);
    if (p != nullptr)
      return p;

    p = Bool::create(description);
    if (p != nullptr)
      return p;

    p = Anything::create(description);
    if (p != nullptr)
      return p;

    p = FileName::create(description);
    if (p != nullptr)
      return p;

    p = DirectoryName::create(description);
    if (p != nullptr)
      return p;

    Assert(false, ExcNotImplemented());

    return p;
  }



  std::size_t
  PatternBase::memory_consumption() const
  {
    if (dynamic_cast<const Integer *>(this) != nullptr)
      return sizeof(Integer);
    else if (dynamic_cast<const Double *>(this) != nullptr)
      return sizeof(Double);
    else if (dynamic_cast<const Bool *>(this) != nullptr)
      return sizeof(Bool);
    else if (dynamic_cast<const Anything *>(this) != nullptr)
      return sizeof(Anything);
    else
      return sizeof(*this) + 32;
  }



  const int Integer::min_int_value = std::numeric_limits<int>::min();
  const int Integer::max_int_value = std::numeric_limits<int>::max();

  const char *Integer::description_init = "[Integer";

  Integer::Integer(const int lower_bound, const int upper_bound)
    : lower_bound(lower_bound)
    , upper_bound(upper_bound)
  {}



  bool
  Integer::match(const std::string &test_string) const
  {
    std::istringstream str(test_string);

    int i;
    if (!(str >> i))
      return false;

    if (!has_only_whitespace(str))
      return false;
    // check whether valid bounds
    // were specified, and if so
    // enforce their values
    if (lower_bound <= upper_bound)
      return ((lower_bound <= i) && (upper_bound >= i));
    else
      return true;
  }



  std::string
  Integer::description(const OutputStyle style) const
  {
    switch (style)
      {
        case Machine:
          {
            // check whether valid bounds
            // were specified, and if so
            // output their values
            if (lower_bound <= upper_bound)
              {
                std::ostringstream description;

                description << description_init << " range " << lower_bound
                            << "..." << upper_bound << " (inclusive)]";
                return description.str();
              }
            else
              // if no bounds were given, then
              // return generic string
              return "[Integer]";
          }
        case Text:
          {
            if (lower_bound <= upper_bound)
              {
                std::ostringstream description;

                description << "An integer n such that " << lower_bound
                            << " <= n <= " << upper_bound;

                return description.str();
              }
            else
              return "An integer";
          }
        case LaTeX:
          {
            if (lower_bound <= upper_bound)
              {
                std::ostringstream description;

                description << "An integer $n$ such that $" << lower_bound
                            << "\\leq n \\leq " << upper_bound << "$";

                return description.str();
              }
            else
              return "An integer";
          }
        default:
          AssertThrow(false, ExcNotImplemented());
      }
    // Should never occur without an exception, but prevent compiler from
    // complaining
    return "";
  }



  std::unique_ptr<PatternBase>
  Integer::clone() const
  {
    return std::unique_ptr<PatternBase>(new Integer(lower_bound, upper_bound));
  }



  std::unique_ptr<Integer>
  Integer::create(const std::string &description)
  {
    if (description.compare(0,
                            std::strlen(description_init),
                            description_init) == 0)
      {
        std::istringstream is(description);

        if (is.str().size() > strlen(description_init) + 1)
          {
            // TODO: verify that description matches the pattern "^\[Integer
            // range \d+\.\.\.\d+\]$"
            int lower_bound, upper_bound;

            is.ignore(strlen(description_init) + strlen(" range "));

            if (!(is >> lower_bound))
              return std::make_unique<Integer>();

            is.ignore(strlen("..."));

            if (!(is >> upper_bound))
              return std::make_unique<Integer>();

            return std::make_unique<Integer>(lower_bound, upper_bound);
          }
        else
          return std::make_unique<Integer>();
      }
    else
      return std::unique_ptr<Integer>();
  }



  const double Double::min_double_value = -std::numeric_limits<double>::max();
  const double Double::max_double_value = std::numeric_limits<double>::max();

  const char *Double::description_init = "[Double";

  Double::Double(const double lower_bound, const double upper_bound)
    : lower_bound(lower_bound)
    , upper_bound(upper_bound)
  {}



  bool
  Double::match(const std::string &test_string) const
  {
    std::istringstream str(test_string);

    double d;
    str >> d;
    if (str.fail())
      return false;

    if (!has_only_whitespace(str))
      return false;
    // check whether valid bounds
    // were specified, and if so
    // enforce their values
    if (lower_bound <= upper_bound)
      return ((lower_bound <= d) && (upper_bound >= d));
    else
      return true;
  }



  std::string
  Double::description(const OutputStyle style) const
  {
    switch (style)
      {
        case Machine:
          {
            std::ostringstream description;

            if (lower_bound <= upper_bound)
              {
                // bounds are valid
                description << description_init << " ";
                // We really want to compare with ==, but -Wfloat-equal would
                // create a warning here, so work around it.
                if (0 == std::memcmp(&lower_bound,
                                     &min_double_value,
                                     sizeof(lower_bound)))
                  description << "-MAX_DOUBLE";
                else
                  description << lower_bound;
                description << "...";
                if (0 == std::memcmp(&upper_bound,
                                     &max_double_value,
                                     sizeof(upper_bound)))
                  description << "MAX_DOUBLE";
                else
                  description << upper_bound;
                description << " (inclusive)]";
                return description.str();
              }
            else
              {
                // invalid bounds, assume unbounded double:
                description << description_init << "]";
                return description.str();
              }
          }
        case Text:
          {
            if (lower_bound <= upper_bound)
              {
                std::ostringstream description;

                description << "A floating point number v such that ";
                if (0 == std::memcmp(&lower_bound,
                                     &min_double_value,
                                     sizeof(lower_bound)))
                  description << "-MAX_DOUBLE";
                else
                  description << lower_bound;
                description << " <= v <= ";
                if (0 == std::memcmp(&upper_bound,
                                     &max_double_value,
                                     sizeof(upper_bound)))
                  description << "MAX_DOUBLE";
                else
                  description << upper_bound;

                return description.str();
              }
            else
              return "A floating point number";
          }
        case LaTeX:
          {
            if (lower_bound <= upper_bound)
              {
                std::ostringstream description;

                description << "A floating point number $v$ such that $";
                if (0 == std::memcmp(&lower_bound,
                                     &min_double_value,
                                     sizeof(lower_bound)))
                  description << "-\\text{MAX\\_DOUBLE}";
                else
                  description << lower_bound;
                description << " \\leq v \\leq ";
                if (0 == std::memcmp(&upper_bound,
                                     &max_double_value,
                                     sizeof(upper_bound)))
                  description << "\\text{MAX\\_DOUBLE}";
                else
                  description << upper_bound;
                description << "$";

                return description.str();
              }
            else
              return "A floating point number";
          }
        default:
          AssertThrow(false, ExcNotImplemented());
      }
    // Should never occur without an exception, but prevent compiler from
    // complaining
    return "";
  }


  std::unique_ptr<PatternBase>
  Double::clone() const
  {
    return std::unique_ptr<PatternBase>(new Double(lower_bound, upper_bound));
  }



  std::unique_ptr<Double>
  Double::create(const std::string &description)
  {
    const std::string description_init_str = description_init;
    if (description.compare(0,
                            description_init_str.size(),
                            description_init_str) != 0)
      return std::unique_ptr<Double>();
    if (*description.rbegin() != ']')
      return std::unique_ptr<Double>();

    std::string temp = description.substr(description_init_str.size());
    if (temp == "]")
      return std::make_unique<Double>(1.0,
                                      -1.0); // return an invalid range

    if (temp.find("...") != std::string::npos)
      temp.replace(temp.find("..."), 3, " ");

    double lower_bound = min_double_value, upper_bound = max_double_value;

    std::istringstream is(temp);
    if (0 == temp.compare(0, std::strlen(" -MAX_DOUBLE"), " -MAX_DOUBLE"))
      is.ignore(std::strlen(" -MAX_DOUBLE"));
    else
      {
        // parse lower bound and give up if not a double
        if (!(is >> lower_bound))
          return std::unique_ptr<Double>();
      }

    // ignore failure here and assume we got MAX_DOUBLE as upper bound:
    is >> upper_bound;
    if (is.fail())
      upper_bound = max_double_value;

    return std::make_unique<Double>(lower_bound, upper_bound);
  }



  const char *Selection::description_init = "[Selection";


  Selection::Selection(const std::string &seq)
    : sequence(seq)
  {
    while (sequence.find(" |") != std::string::npos)
      sequence.replace(sequence.find(" |"), 2, "|");
    while (sequence.find("| ") != std::string::npos)
      sequence.replace(sequence.find("| "), 2, "|");
  }



  bool
  Selection::match(const std::string &test_string) const
  {
    std::string tmp(sequence);

    // remove whitespace at beginning
    while ((tmp.length() != 0) && (std::isspace(tmp[0])))
      tmp.erase(0, 1);

    // check the different possibilities
    while (tmp.find('|') != std::string::npos)
      {
        if (test_string == std::string(tmp, 0, tmp.find('|')))
          return true;

        tmp.erase(0, tmp.find('|') + 1);
      }

    // remove whitespace at the end
    while ((tmp.length() != 0) && (std::isspace(*(tmp.end() - 1))))
      tmp.erase(tmp.end() - 1);

    // check last choice, not finished by |
    if (test_string == tmp)
      return true;

    // not found
    return false;
  }



  std::string
  Selection::description(const OutputStyle style) const
  {
    switch (style)
      {
        case Machine:
          {
            std::ostringstream description;

            description << description_init << " " << sequence << " ]";

            return description.str();
          }
        case Text:
        case LaTeX:
          {
            std::ostringstream description;

            description << "Any one of "
                        << internal::escape(
                             Utilities::replace_in_string(sequence, "|", ", "),
                             style);

            return description.str();
          }
        default:
          AssertThrow(false, ExcNotImplemented());
      }
    // Should never occur without an exception, but prevent compiler from
    // complaining
    return "";
  }



  std::unique_ptr<PatternBase>
  Selection::clone() const
  {
    return std::unique_ptr<PatternBase>(new Selection(sequence));
  }


  std::size_t
  Selection::memory_consumption() const
  {
    return (sizeof(PatternBase) +
            MemoryConsumption::memory_consumption(sequence));
  }



  std::unique_ptr<Selection>
  Selection::create(const std::string &description)
  {
    if (description.compare(0,
                            std::strlen(description_init),
                            description_init) == 0)
      {
        std::string sequence(description);

        sequence.erase(0, std::strlen(description_init) + 1);
        sequence.erase(sequence.length() - 2, 2);

        return std::make_unique<Selection>(sequence);
      }
    else
      return std::unique_ptr<Selection>();
  }



  const unsigned int List::max_int_value =
    std::numeric_limits<unsigned int>::max();

  const char *List::description_init = "[List";


  List::List(const PatternBase &p,
             const unsigned int min_elements,
             const unsigned int max_elements,
             const std::string &separator)
    : pattern(p.clone())
    , min_elements(min_elements)
    , max_elements(max_elements)
    , separator(separator)
  {
    Assert(min_elements <= max_elements,
           ExcInvalidRange(min_elements, max_elements));
    Assert(separator.size() > 0,
           ExcMessage("The separator must have a non-zero length."));
  }



  List::List(const List &other)
    : pattern(other.pattern->clone())
    , min_elements(other.min_elements)
    , max_elements(other.max_elements)
    , separator(other.separator)
  {}


  const std::string &
  List::get_separator() const
  {
    return separator;
  }



  const PatternBase &
  List::get_base_pattern() const
  {
    return *pattern;
  }



  bool
  List::match(const std::string &test_string_list) const
  {
    const std::vector<std::string> split_list =
      Utilities::split_string_list(test_string_list, separator);

    if ((split_list.size() < min_elements) ||
        (split_list.size() > max_elements))
      return false;

    // check the different possibilities
    for (const std::string &string : split_list)
      if (pattern->match(string) == false)
        return false;

    return true;
  }



  std::string
  List::description(const OutputStyle style) const
  {
    switch (style)
      {
        case Machine:
          {
            std::ostringstream description;

            description << description_init << " of <"
                        << pattern->description(style) << ">"
                        << " of length " << min_elements << "..."
                        << max_elements << " (inclusive)";
            if (separator != ",")
              description << " separated by <" << separator << ">";
            description << "]";

            return description.str();
          }
        case Text:
        case LaTeX:
          {
            std::ostringstream description;

            description << "A list of " << min_elements << " to "
                        << max_elements << " elements ";
            if (separator != ",")
              description << "separated by <"
                          << internal::escape(separator, style) << "> ";
            description << "where each element is ["
                        << pattern->description(style) << "]";

            return description.str();
          }
        default:
          AssertThrow(false, ExcNotImplemented());
      }
    // Should never occur without an exception, but prevent compiler from
    // complaining
    return "";
  }



  std::unique_ptr<PatternBase>
  List::clone() const
  {
    return std::unique_ptr<PatternBase>(
      new List(*pattern, min_elements, max_elements, separator));
  }


  std::size_t
  List::memory_consumption() const
  {
    return (sizeof(*this) + MemoryConsumption::memory_consumption(*pattern) +
            MemoryConsumption::memory_consumption(separator));
  }



  std::unique_ptr<List>
  List::create(const std::string &description)
  {
    if (description.compare(0,
                            std::strlen(description_init),
                            description_init) == 0)
      {
        unsigned int min_elements = 0, max_elements = 0;

        std::istringstream is(description);
        is.ignore(strlen(description_init) + strlen(" of <"));

        std::string str;
        std::getline(is, str, '>');

        std::unique_ptr<PatternBase> base_pattern(pattern_factory(str));

        is.ignore(strlen(" of length "));
        if (!(is >> min_elements))
          return std::make_unique<List>(*base_pattern);

        is.ignore(strlen("..."));
        if (!(is >> max_elements))
          return std::make_unique<List>(*base_pattern, min_elements);

        is.ignore(strlen(" (inclusive) separated by <"));
        std::string separator;
        if (!is.eof())
          std::getline(is, separator, '>');
        else
          separator = ",";

        return std::make_unique<List>(*base_pattern,
                                      min_elements,
                                      max_elements,
                                      separator);
      }
    else
      return std::unique_ptr<List>();
  }



  const unsigned int Map::max_int_value =
    std::numeric_limits<unsigned int>::max();

  const char *Map::description_init = "[Map";


  Map::Map(const PatternBase &p_key,
           const PatternBase &p_value,
           const unsigned int min_elements,
           const unsigned int max_elements,
           const std::string &separator,
           const std::string &key_value_separator)
    : key_pattern(p_key.clone())
    , value_pattern(p_value.clone())
    , min_elements(min_elements)
    , max_elements(max_elements)
    , separator(separator)
    , key_value_separator(key_value_separator)
  {
    Assert(min_elements <= max_elements,
           ExcInvalidRange(min_elements, max_elements));
    Assert(separator.size() > 0,
           ExcMessage("The separator must have a non-zero length."));
    Assert(key_value_separator.size() > 0,
           ExcMessage("The key_value_separator must have a non-zero length."));
    Assert(separator != key_value_separator,
           ExcMessage(
             "The separator can not be the same of the key_value_separator "
             "since that is used as the separator between the two elements "
             "of <key:value> pairs"));
  }



  Map::Map(const Map &other)
    : key_pattern(other.key_pattern->clone())
    , value_pattern(other.value_pattern->clone())
    , min_elements(other.min_elements)
    , max_elements(other.max_elements)
    , separator(other.separator)
    , key_value_separator(other.key_value_separator)
  {}



  bool
  Map::match(const std::string &test_string_list) const
  {
    std::vector<std::string> split_list =
      Utilities::split_string_list(test_string_list, separator);
    if ((split_list.size() < min_elements) ||
        (split_list.size() > max_elements))
      return false;

    for (const auto &key_value_pair : split_list)
      {
        std::vector<std::string> pair =
          Utilities::split_string_list(key_value_pair, key_value_separator);

        // Check that we have in fact two matches
        if (pair.size() != 2)
          return false;

        // then verify that the patterns are satisfied
        if (key_pattern->match(pair[0]) == false)
          return false;
        if (value_pattern->match(pair[1]) == false)
          return false;
      }

    return true;
  }



  std::string
  Map::description(const OutputStyle style) const
  {
    switch (style)
      {
        case Machine:
          {
            std::ostringstream description;

            description << description_init << " of <"
                        << key_pattern->description(style) << ">"
                        << key_value_separator << "<"
                        << value_pattern->description(style) << ">"
                        << " of length " << min_elements << "..."
                        << max_elements << " (inclusive)";
            if (separator != ",")
              description << " separated by <" << separator << ">";
            description << "]";

            return description.str();
          }
        case Text:
        case LaTeX:
          {
            std::ostringstream description;

            description << "A key"
                        << internal::escape(key_value_separator, style)
                        << "value map of " << min_elements << " to "
                        << max_elements << " elements ";
            if (separator != ",")
              description << " separated by <"
                          << internal::escape(separator, style) << "> ";
            description << " where each key is ["
                        << key_pattern->description(style) << "]"
                        << " and each value is ["
                        << value_pattern->description(style) << "]";

            return description.str();
          }
        default:
          AssertThrow(false, ExcNotImplemented());
      }
    // Should never occur without an exception, but prevent compiler from
    // complaining
    return "";
  }



  std::unique_ptr<PatternBase>
  Map::clone() const
  {
    return std::unique_ptr<PatternBase>(new Map(*key_pattern,
                                                *value_pattern,
                                                min_elements,
                                                max_elements,
                                                separator,
                                                key_value_separator));
  }


  std::size_t
  Map::memory_consumption() const
  {
    return (sizeof(*this) +
            MemoryConsumption::memory_consumption(*key_pattern) +
            MemoryConsumption::memory_consumption(*value_pattern) +
            MemoryConsumption::memory_consumption(separator) +
            MemoryConsumption::memory_consumption(key_value_separator));
  }



  std::unique_ptr<Map>
  Map::create(const std::string &description)
  {
    if (description.compare(0,
                            std::strlen(description_init),
                            description_init) == 0)
      {
        unsigned int min_elements = 0, max_elements = 0;

        std::istringstream is(description);
        is.ignore(strlen(description_init) + strlen(" of <"));

        std::string key;
        std::getline(is, key, '>');

        std::string key_value_separator;
        std::getline(is, key_value_separator, '<');

        // split 'str' into key and value
        std::string value;
        std::getline(is, value, '>');

        std::unique_ptr<PatternBase> key_pattern(pattern_factory(key));
        std::unique_ptr<PatternBase> value_pattern(pattern_factory(value));

        is.ignore(strlen(" of length "));
        if (!(is >> min_elements))
          return std::make_unique<Map>(*key_pattern, *value_pattern);

        is.ignore(strlen("..."));
        if (!(is >> max_elements))
          return std::make_unique<Map>(*key_pattern,
                                       *value_pattern,
                                       min_elements);

        is.ignore(strlen(" (inclusive) separated by <"));
        std::string separator;
        if (!is.eof())
          std::getline(is, separator, '>');
        else
          separator = ",";

        return std::make_unique<Map>(*key_pattern,
                                     *value_pattern,
                                     min_elements,
                                     max_elements,
                                     separator,
                                     key_value_separator);
      }
    else
      return std::unique_ptr<Map>();
  }



  const PatternBase &
  Map::get_key_pattern() const
  {
    return *key_pattern;
  }



  const PatternBase &
  Map::get_value_pattern() const
  {
    return *value_pattern;
  }



  const std::string &
  Map::get_separator() const
  {
    return separator;
  }


  const std::string &
  Map::get_key_value_separator() const
  {
    return key_value_separator;
  }



  const char *Tuple::description_init = "[Tuple";


  Tuple::Tuple(const std::vector<std::unique_ptr<PatternBase>> &ps,
               const std::string &                              separator)
    : separator(separator)
  {
    Assert(ps.size() > 0,
           ExcMessage("The Patterns list must have a non-zero length."));
    Assert(separator.size() > 0,
           ExcMessage("The separator must have a non-zero length."));
    patterns.resize(ps.size());
    for (unsigned int i = 0; i < ps.size(); ++i)
      patterns[i] = ps[i]->clone();
  }



  Tuple::Tuple(const std::vector<std::unique_ptr<PatternBase>> &ps,
               const char *                                     separator)
    : Tuple(ps, std::string(separator))
  {}



  Tuple::Tuple(const Tuple &other)
    : separator(other.separator)
  {
    patterns.resize(other.patterns.size());
    for (unsigned int i = 0; i < other.patterns.size(); ++i)
      patterns[i] = other.patterns[i]->clone();
  }



  bool
  Tuple::match(const std::string &test_string_list) const
  {
    std::vector<std::string> split_list =
      Utilities::split_string_list(test_string_list, separator);
    if (split_list.size() != patterns.size())
      return false;

    for (unsigned int i = 0; i < patterns.size(); ++i)
      {
        if (patterns[i]->match(split_list[i]) == false)
          return false;
      }

    return true;
  }



  std::string
  Tuple::description(const OutputStyle style) const
  {
    switch (style)
      {
        case Machine:
          {
            std::ostringstream description;

            description << description_init << " of <" << patterns.size()
                        << "> elements <" << patterns[0]->description(style)
                        << ">";
            for (unsigned int i = 1; i < patterns.size(); ++i)
              description << ", <" << patterns[i]->description(style) << ">";

            if (separator != ":")
              description << " separated by <" << separator << ">";
            description << "]";

            return description.str();
          }
        case Text:
        case LaTeX:
          {
            std::ostringstream description;

            description << "A Tuple of " << patterns.size() << " elements ";
            if (separator != ":")
              description << " separated by <"
                          << internal::escape(separator, style) << "> ";
            description << " where each element is ["
                        << patterns[0]->description(style) << "]";
            for (unsigned int i = 1; i < patterns.size(); ++i)
              {
                description << internal::escape(separator, style) << "["
                            << patterns[i]->description(style) << "]";
              }
            return description.str();
          }

        default:
          AssertThrow(false, ExcNotImplemented());
      }
    // Should never occur without an exception, but prevent compiler from
    // complaining
    return "";
  }



  std::unique_ptr<PatternBase>
  Tuple::clone() const
  {
    return std::unique_ptr<PatternBase>(new Tuple(patterns, separator));
  }


  std::size_t
  Tuple::memory_consumption() const
  {
    return (sizeof(*this) + MemoryConsumption::memory_consumption(patterns) +
            MemoryConsumption::memory_consumption(separator));
  }



  std::unique_ptr<Tuple>
  Tuple::create(const std::string &description)
  {
    if (description.compare(0,
                            std::strlen(description_init),
                            description_init) == 0)
      {
        std::vector<std::unique_ptr<PatternBase>> patterns;

        std::istringstream is(description);
        is.ignore(strlen(description_init) + strlen(" of <"));

        std::string len;
        std::getline(is, len, '>');
        const unsigned int n_elements = Utilities::string_to_int(len);
        Assert(n_elements > 0,
               ExcMessage("Provide at least 1 element in the tuple."));
        patterns.resize(n_elements);

        is.ignore(strlen(" elements <"));

        std::string element;
        std::getline(is, element, '>');
        patterns[0] = pattern_factory(element);

        for (unsigned int i = 1; i < n_elements; ++i)
          {
            is.ignore(strlen(", <"));
            std::getline(is, element, '>');
            patterns[i] = pattern_factory(element);
          }

        is.ignore(strlen(" separated by <"));

        std::string separator;
        if (!is.eof())
          std::getline(is, separator, '>');
        else
          separator = ":";

        return std::make_unique<Tuple>(patterns, separator);
      }
    else
      return std::unique_ptr<Tuple>();
  }



  const PatternBase &
  Tuple::get_pattern(const unsigned int i) const
  {
    return *patterns[i];
  }



  const std::string &
  Tuple::get_separator() const
  {
    return separator;
  }



  const char *MultipleSelection::description_init = "[MultipleSelection";


  MultipleSelection::MultipleSelection(const std::string &seq)
  {
    Assert(seq.find(',') == std::string::npos,
           ExcCommasNotAllowed(seq.find(',')));

    sequence = seq;
    while (sequence.find(" |") != std::string::npos)
      sequence.replace(sequence.find(" |"), 2, "|");
    while (sequence.find("| ") != std::string::npos)
      sequence.replace(sequence.find("| "), 2, "|");
  }



  bool
  MultipleSelection::match(const std::string &test_string_list) const
  {
    std::string              tmp = test_string_list;
    std::vector<std::string> split_names;

    // first split the input list
    while (tmp.length() != 0)
      {
        std::string name;
        name = tmp;

        if (name.find(',') != std::string::npos)
          {
            name.erase(name.find(','), std::string::npos);
            tmp.erase(0, tmp.find(',') + 1);
          }
        else
          tmp = "";

        while ((name.length() != 0) && (std::isspace(name[0])))
          name.erase(0, 1);
        while (std::isspace(name[name.length() - 1]))
          name.erase(name.length() - 1, 1);

        split_names.push_back(name);
      }


    // check the different possibilities
    for (const auto &test_string : split_names)
      {
        bool string_found = false;

        tmp = sequence;
        while (tmp.find('|') != std::string::npos)
          {
            if (test_string == std::string(tmp, 0, tmp.find('|')))
              {
                // string found, quit
                // loop. don't change
                // tmp, since we don't
                // need it anymore.
                string_found = true;
                break;
              }

            tmp.erase(0, tmp.find('|') + 1);
          }
        // check last choice, not finished by |
        if (!string_found)
          if (test_string == tmp)
            string_found = true;

        if (!string_found)
          return false;
      }

    return true;
  }



  std::string
  MultipleSelection::description(const OutputStyle style) const
  {
    switch (style)
      {
        case Machine:
          {
            std::ostringstream description;

            description << description_init << " " << sequence << " ]";

            return description.str();
          }
        case Text:
        case LaTeX:
          {
            std::ostringstream description;

            description << "A comma-separated list of any of "
                        << internal::escape(
                             Utilities::replace_in_string(sequence, "|", ", "),
                             style);

            return description.str();
          }
        default:
          AssertThrow(false, ExcNotImplemented());
      }
    // Should never occur without an exception, but prevent compiler from
    // complaining
    return "";
  }



  std::unique_ptr<PatternBase>
  MultipleSelection::clone() const
  {
    return std::unique_ptr<PatternBase>(new MultipleSelection(sequence));
  }


  std::size_t
  MultipleSelection::memory_consumption() const
  {
    return (sizeof(PatternBase) +
            MemoryConsumption::memory_consumption(sequence));
  }



  std::unique_ptr<MultipleSelection>
  MultipleSelection::create(const std::string &description)
  {
    if (description.compare(0,
                            std::strlen(description_init),
                            description_init) == 0)
      {
        std::string sequence(description);

        sequence.erase(0, std::strlen(description_init) + 1);
        sequence.erase(sequence.length() - 2, 2);

        return std::make_unique<MultipleSelection>(sequence);
      }
    else
      return std::unique_ptr<MultipleSelection>();
  }



  const char *Bool::description_init = "[Bool";


  Bool::Bool()
    : Selection("true|false")
  {}



  std::string
  Bool::description(const OutputStyle style) const
  {
    switch (style)
      {
        case Machine:
          {
            std::ostringstream description;

            description << description_init << "]";

            return description.str();
          }
        case Text:
        case LaTeX:
          {
            return "A boolean value (true or false)";
          }
        default:
          AssertThrow(false, ExcNotImplemented());
      }
    // Should never occur without an exception, but prevent compiler from
    // complaining
    return "";
  }



  std::unique_ptr<PatternBase>
  Bool::clone() const
  {
    return std::unique_ptr<PatternBase>(new Bool());
  }



  std::unique_ptr<Bool>
  Bool::create(const std::string &description)
  {
    if (description.compare(0,
                            std::strlen(description_init),
                            description_init) == 0)
      return std::make_unique<Bool>();
    else
      return std::unique_ptr<Bool>();
  }



  const char *Anything::description_init = "[Anything";



  bool
  Anything::match(const std::string &) const
  {
    return true;
  }



  std::string
  Anything::description(const OutputStyle style) const
  {
    switch (style)
      {
        case Machine:
          {
            std::ostringstream description;

            description << description_init << "]";

            return description.str();
          }
        case Text:
        case LaTeX:
          {
            return "Any string";
          }
        default:
          AssertThrow(false, ExcNotImplemented());
      }
    // Should never occur without an exception, but prevent compiler from
    // complaining
    return "";
  }



  std::unique_ptr<PatternBase>
  Anything::clone() const
  {
    return std::unique_ptr<PatternBase>(new Anything());
  }



  std::unique_ptr<Anything>
  Anything::create(const std::string &description)
  {
    if (description.compare(0,
                            std::strlen(description_init),
                            description_init) == 0)
      return std::make_unique<Anything>();
    else
      return std::unique_ptr<Anything>();
  }



  const char *FileName::description_init = "[FileName";


  FileName::FileName(const FileType type)
    : file_type(type)
  {}



  bool
  FileName::match(const std::string &) const
  {
    return true;
  }



  std::string
  FileName::description(const OutputStyle style) const
  {
    switch (style)
      {
        case Machine:
          {
            std::ostringstream description;

            description << description_init;

            if (file_type == input)
              description << " (Type: input)]";
            else
              description << " (Type: output)]";

            return description.str();
          }
        case Text:
        case LaTeX:
          {
            if (file_type == input)
              return "an input filename";
            else
              return "an output filename";
          }
        default:
          AssertThrow(false, ExcNotImplemented());
      }
    // Should never occur without an exception, but prevent compiler from
    // complaining
    return "";
  }



  std::unique_ptr<PatternBase>
  FileName::clone() const
  {
    return std::unique_ptr<PatternBase>(new FileName(file_type));
  }



  std::unique_ptr<FileName>
  FileName::create(const std::string &description)
  {
    if (description.compare(0,
                            std::strlen(description_init),
                            description_init) == 0)
      {
        std::istringstream is(description);
        std::string        file_type;
        FileType           type;

        is.ignore(strlen(description_init) + strlen(" (Type:"));

        is >> file_type;

        if (file_type == "input)]")
          type = input;
        else
          type = output;

        return std::make_unique<FileName>(type);
      }
    else
      return std::unique_ptr<FileName>();
  }



  const char *DirectoryName::description_init = "[DirectoryName";



  bool
  DirectoryName::match(const std::string &) const
  {
    return true;
  }



  std::string
  DirectoryName::description(const OutputStyle style) const
  {
    switch (style)
      {
        case Machine:
          {
            std::ostringstream description;

            description << description_init << "]";

            return description.str();
          }
        case Text:
        case LaTeX:
          {
            return "A directory name";
          }
        default:
          AssertThrow(false, ExcNotImplemented());
      }
    // Should never occur without an exception, but prevent compiler from
    // complaining
    return "";
  }



  std::unique_ptr<PatternBase>
  DirectoryName::clone() const
  {
    return std::unique_ptr<PatternBase>(new DirectoryName());
  }



  std::unique_ptr<DirectoryName>
  DirectoryName::create(const std::string &description)
  {
    if (description.compare(0,
                            std::strlen(description_init),
                            description_init) == 0)
      return std::make_unique<DirectoryName>();
    else
      return std::unique_ptr<DirectoryName>();
  }

} // end namespace Patterns

DEAL_II_NAMESPACE_CLOSE
