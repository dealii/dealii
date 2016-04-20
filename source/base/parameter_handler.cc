// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2016 by the deal.II authors
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


#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/path_search.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/utilities.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <algorithm>
#include <list>
#include <sstream>
#include <cctype>
#include <limits>


DEAL_II_NAMESPACE_OPEN



//TODO[WB]: various functions here could be simplified by using namespace Utilities

namespace Patterns
{

  namespace
  {
    /**
     * Read to the end of the stream and
     * return whether all there is is
     * whitespace or whether there are other
     * characters as well.
     */
    bool has_only_whitespace (std::istream &in)
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
  }



  PatternBase *pattern_factory (const std::string &description)
  {
    PatternBase *p;

    p = Integer::create(description);
    if (p != 0)
      return p;

    p = Double::create(description);
    if (p !=0 )
      return p;

    p = Selection::create(description);
    if (p !=0 )
      return p;

    p = List::create(description);
    if (p !=0 )
      return p;

    p = MultipleSelection::create(description);
    if (p !=0 )
      return p;

    p = Bool::create(description);
    if (p!=0 )
      return p;

    p = Anything::create(description);
    if (p !=0 )
      return p;

    p = FileName::create(description);
    if (p !=0 )
      return p;

    p = DirectoryName::create(description);
    if (p!=0 )
      return p;

    Assert(false, ExcNotImplemented());

    return 0;
  }



  PatternBase::~PatternBase ()
  {}


  std::size_t
  PatternBase::memory_consumption () const
  {
    if (dynamic_cast<const Integer *>(this) != 0)
      return sizeof(Integer);
    else if (dynamic_cast<const Double *>(this) != 0)
      return sizeof(Double);
    else if (dynamic_cast<const Bool *>(this) != 0)
      return sizeof(Bool);
    else if (dynamic_cast<const Anything *>(this) != 0)
      return sizeof(Anything);
    else
      return sizeof(*this) + 32;
  }



  const int Integer::min_int_value = std::numeric_limits<int>::min();
  const int Integer::max_int_value = std::numeric_limits<int>::max();

  const char *Integer::description_init = "[Integer";

  Integer::Integer (const int lower_bound,
                    const int upper_bound)
    :
    lower_bound (lower_bound),
    upper_bound (upper_bound)
  {}



  bool Integer::match (const std::string &test_string) const
  {
    std::istringstream str(test_string);

    int i;
    if (!(str >> i))
      return false;

    if (!has_only_whitespace (str))
      return false;
    // check whether valid bounds
    // were specified, and if so
    // enforce their values
    if (lower_bound <= upper_bound)
      return ((lower_bound <= i) &&
              (upper_bound >= i));
    else
      return true;
  }



  std::string Integer::description () const
  {
    // check whether valid bounds
    // were specified, and if so
    // output their values
    if (lower_bound <= upper_bound)
      {
        std::ostringstream description;

        description << description_init
                    <<" range "
                    << lower_bound << "..." << upper_bound
                    << " (inclusive)]";
        return description.str();
      }
    else
      // if no bounds were given, then
      // return generic string
      return "[Integer]";
  }



  PatternBase *
  Integer::clone () const
  {
    return new Integer(lower_bound, upper_bound);
  }



  Integer *Integer::create (const std::string &description)
  {
    if (description.compare(0, std::strlen(description_init), description_init) == 0)
      {
        std::istringstream is(description);

        if (is.str().size() > strlen(description_init) + 1)
          {
//TODO: verify that description matches the pattern "^\[Integer range \d+\.\.\.\d+\]$"
            int lower_bound, upper_bound;

            is.ignore(strlen(description_init) + strlen(" range "));

            if (!(is >> lower_bound))
              return new Integer();

            is.ignore(strlen("..."));

            if (!(is >> upper_bound))
              return new Integer();

            return new Integer(lower_bound, upper_bound);
          }
        else
          return new Integer();
      }
    else
      return 0;
  }



  const double Double::min_double_value = -std::numeric_limits<double>::max();
  const double Double::max_double_value = std::numeric_limits<double>::max();

  const char *Double::description_init = "[Double";

  Double::Double (const double lower_bound,
                  const double upper_bound)
    :
    lower_bound (lower_bound),
    upper_bound (upper_bound)
  {}



  bool Double::match (const std::string &test_string) const
  {
    std::istringstream str(test_string);

    double d;
    if (!(str >> d))
      return false;

    if (!has_only_whitespace (str))
      return false;
    // check whether valid bounds
    // were specified, and if so
    // enforce their values
    if (lower_bound <= upper_bound)
      return ((lower_bound <= d) &&
              (upper_bound >= d));
    else
      return true;
  }



  std::string Double::description () const
  {
    std::ostringstream description;

    // check whether valid bounds
    // were specified, and if so
    // output their values
    if (lower_bound <= upper_bound)
      {
        description << description_init
                    << " "
                    << lower_bound << "..." << upper_bound
                    << " (inclusive)]";
        return description.str();
      }
    else
      // if no bounds were given, then
      // return generic string
      {
        description << description_init
                    << "]";
        return description.str();
      }
  }


  PatternBase *
  Double::clone () const
  {
    return new Double(lower_bound, upper_bound);
  }



  Double *Double::create (const std::string &description)
  {
    if (description.compare(0, std::strlen(description_init), description_init) == 0)
      {
        std::istringstream is(description);

        if (is.str().size() > strlen(description_init) + 1)
          {
            double lower_bound, upper_bound;

            is.ignore(strlen(description_init) + strlen(" range "));

            if (!(is >> lower_bound))
              return new Double();

            is.ignore(strlen("..."));

            if (!(is >> upper_bound))
              return new Double();

            return new Double(lower_bound, upper_bound);
          }
        else
          return new Double();
      }
    else
      return 0;
  }



  const char *Selection::description_init = "[Selection";


  Selection::Selection (const std::string &seq)
  {
    sequence = seq;

    while (sequence.find(" |") != std::string::npos)
      sequence.replace (sequence.find(" |"), 2, "|");
    while (sequence.find("| ") != std::string::npos)
      sequence.replace (sequence.find("| "), 2, "|");
  }



  bool Selection::match (const std::string &test_string) const
  {
    std::vector<std::string> choices;
    std::string tmp(sequence);
    // check the different possibilities
    while (tmp.find('|') != std::string::npos)
      {
        if (test_string == std::string(tmp, 0, tmp.find('|')))
          return true;

        tmp.erase (0, tmp.find('|')+1);
      };
    // check last choice, not finished by |
    if (test_string == tmp)
      return true;

    // not found
    return false;
  }



  std::string Selection::description () const
  {
    std::ostringstream description;

    description << description_init
                << " "
                << sequence
                << " ]";

    return description.str();
  }



  PatternBase *
  Selection::clone () const
  {
    return new Selection(sequence);
  }


  std::size_t
  Selection::memory_consumption () const
  {
    return (sizeof(PatternBase) +
            MemoryConsumption::memory_consumption(sequence));
  }



  Selection *Selection::create (const std::string &description)
  {
    if (description.compare(0, std::strlen(description_init), description_init) == 0)
      {
        std::string sequence(description);

        sequence.erase(0, std::strlen(description_init) + 1);
        sequence.erase(sequence.length()-2, 2);

        return new Selection(sequence);
      }
    else
      return 0;
  }



  const unsigned int List::max_int_value
    = std::numeric_limits<unsigned int>::max();

  const char *List::description_init = "[List";


  List::List (const PatternBase  &p,
              const unsigned int  min_elements,
              const unsigned int  max_elements,
              const std::string  &separator)
    :
    pattern (p.clone()),
    min_elements (min_elements),
    max_elements (max_elements),
    separator (separator)
  {
    Assert (min_elements <= max_elements,
            ExcInvalidRange (min_elements, max_elements));
    Assert (separator.size() > 0,
            ExcMessage ("The separator must have a non-zero length."));
  }



  List::~List ()
  {
    delete pattern;
    pattern = 0;
  }



  bool List::match (const std::string &test_string_list) const
  {
    std::string tmp = test_string_list;
    std::vector<std::string> split_list;

    // first split the input list
    while (tmp.length() != 0)
      {
        std::string name;
        name = tmp;

        if (name.find(separator) != std::string::npos)
          {
            name.erase (name.find(separator), std::string::npos);
            tmp.erase (0, tmp.find(separator)+separator.size());
          }
        else
          tmp = "";

        while ((name.length() != 0) &&
               (std::isspace (name[0])))
          name.erase (0,1);

        while (std::isspace (name[name.length()-1]))
          name.erase (name.length()-1, 1);

        split_list.push_back (name);
      }

    if ((split_list.size() < min_elements) ||
        (split_list.size() > max_elements))
      return false;

    // check the different possibilities
    for (std::vector<std::string>::const_iterator
         test_string = split_list.begin();
         test_string != split_list.end(); ++test_string)
      if (pattern->match (*test_string) == false)
        return false;

    return true;
  }



  std::string List::description () const
  {
    std::ostringstream description;

    description << description_init
                << " list of <" << pattern->description() << ">"
                << " of length " << min_elements << "..." << max_elements
                << " (inclusive)";
    if (separator != ",")
      description << " separated by <" << separator << ">";
    description << "]";

    return description.str();
  }



  PatternBase *
  List::clone () const
  {
    return new List(*pattern, min_elements, max_elements, separator);
  }


  std::size_t
  List::memory_consumption () const
  {
    return (sizeof(*this) +
            MemoryConsumption::memory_consumption(*pattern) +
            MemoryConsumption::memory_consumption(separator));
  }



  List *List::create (const std::string &description)
  {
    if (description.compare(0, std::strlen(description_init), description_init) == 0)
      {
        int min_elements, max_elements;

        std::istringstream is(description);
        is.ignore(strlen(description_init) + strlen(" list of <"));

        std::string str;
        std::getline(is, str, '>');

        std_cxx11::shared_ptr<PatternBase> base_pattern (pattern_factory(str));

        is.ignore(strlen(" of length "));
        if (!(is >> min_elements))
          return new List(*base_pattern);

        is.ignore(strlen("..."));
        if (!(is >> max_elements))
          return new List(*base_pattern, min_elements);

        is.ignore(strlen(" separated by <"));
        std::string separator;
        if (!is)
          std::getline(is, separator, '>');
        else
          separator = ",";

        return new List(*base_pattern, min_elements, max_elements, separator);
      }
    else
      return 0;
  }



  const unsigned int Map::max_int_value
    = std::numeric_limits<unsigned int>::max();

  const char *Map::description_init = "[Map";


  Map::Map (const PatternBase  &p_key,
            const PatternBase  &p_value,
            const unsigned int  min_elements,
            const unsigned int  max_elements,
            const std::string  &separator)
    :
    key_pattern (p_key.clone()),
    value_pattern (p_value.clone()),
    min_elements (min_elements),
    max_elements (max_elements),
    separator (separator)
  {
    Assert (min_elements <= max_elements,
            ExcInvalidRange (min_elements, max_elements));
    Assert (separator.size() > 0,
            ExcMessage ("The separator must have a non-zero length."));
    Assert (separator != ":",
            ExcMessage ("The separator can not be a colon ':' since that "
                        "is the separator between the two elements of <key:value> pairs"));
  }



  Map::~Map ()
  {
    delete key_pattern;
    key_pattern = 0;

    delete value_pattern;
    value_pattern = 0;
  }



  bool Map::match (const std::string &test_string_list) const
  {
    std::string tmp = test_string_list;
    std::vector<std::string> split_list;

    // first split the input list at comma sites
    while (tmp.length() != 0)
      {
        std::string map_entry;
        map_entry = tmp;

        if (map_entry.find(separator) != std::string::npos)
          {
            map_entry.erase (map_entry.find(separator), std::string::npos);
            tmp.erase (0, tmp.find(separator)+separator.size());
          }
        else
          tmp = "";

        while ((map_entry.length() != 0) &&
               (std::isspace (map_entry[0])))
          map_entry.erase (0,1);

        while (std::isspace (map_entry[map_entry.length()-1]))
          map_entry.erase (map_entry.length()-1, 1);

        split_list.push_back (map_entry);
      }

    if ((split_list.size() < min_elements) ||
        (split_list.size() > max_elements))
      return false;

    // check the different possibilities
    for (std::vector<std::string>::const_iterator
         test_string = split_list.begin();
         test_string != split_list.end(); ++test_string)
      {
        // separate key and value from the test_string
        if (test_string->find(":") == std::string::npos)
          return false;

        // we know now that there is a ':', so split the string there
        // and trim spaces
        std::string key = *test_string;
        key.erase (key.find(":"), std::string::npos);
        while ((key.length() > 0) && (std::isspace (key[key.length()-1])))
          key.erase (key.length()-1, 1);

        std::string value = *test_string;
        value.erase (0, value.find(":")+1);
        while ((value.length() > 0) && (std::isspace (value[0])))
          value.erase (0, 1);

        // then verify that the patterns are satisfied
        if (key_pattern->match (key) == false)
          return false;
        if (value_pattern->match (value) == false)
          return false;
      }

    return true;
  }



  std::string Map::description () const
  {
    std::ostringstream description;

    description << description_init
                << " map of <"
                << key_pattern->description() << ":"
                << value_pattern->description() << ">"
                << " of length " << min_elements << "..." << max_elements
                << " (inclusive)";
    if (separator != ",")
      description << " separated by <" << separator << ">";
    description << "]";

    return description.str();
  }



  PatternBase *
  Map::clone () const
  {
    return new Map(*key_pattern, *value_pattern,
                   min_elements, max_elements,
                   separator);
  }


  std::size_t
  Map::memory_consumption () const
  {
    return (sizeof(*this) +
            MemoryConsumption::memory_consumption (*key_pattern) +
            MemoryConsumption::memory_consumption (*value_pattern) +
            MemoryConsumption::memory_consumption (separator));
  }



  Map *Map::create (const std::string &description)
  {
    if (description.compare(0, std::strlen(description_init), description_init) == 0)
      {
        int min_elements, max_elements;

        std::istringstream is(description);
        is.ignore(strlen(description_init) + strlen(" map of <"));

        std::string str;
        std::getline(is, str, '>');

        // split 'str' into key and value
        std::string key = str;
        key.erase (key.find(":"), std::string::npos);

        std::string value = str;
        value.erase (0, value.find(":")+1);

        std_cxx11::shared_ptr<PatternBase> key_pattern (pattern_factory(key));
        std_cxx11::shared_ptr<PatternBase> value_pattern (pattern_factory(value));

        is.ignore(strlen(" of length "));
        if (!(is >> min_elements))
          return new Map(*key_pattern, *value_pattern);

        is.ignore(strlen("..."));
        if (!(is >> max_elements))
          return new Map(*key_pattern, *value_pattern, min_elements);

        is.ignore(strlen(" separated by <"));
        std::string separator;
        if (!is)
          std::getline(is, separator, '>');
        else
          separator = ",";

        return new Map(*key_pattern, *value_pattern,
                       min_elements, max_elements,
                       separator);
      }
    else
      return 0;
  }



  const char *MultipleSelection::description_init = "[MultipleSelection";


  MultipleSelection::MultipleSelection (const std::string &seq)
  {
    Assert (seq.find (",") == std::string::npos, ExcCommasNotAllowed(seq.find(",")));

    sequence = seq;
    while (sequence.find(" |") != std::string::npos)
      sequence.replace (sequence.find(" |"), 2, "|");
    while (sequence.find("| ") != std::string::npos)
      sequence.replace (sequence.find("| "), 2, "|");
  }



  bool MultipleSelection::match (const std::string &test_string_list) const
  {
    std::string tmp = test_string_list;
    std::list<std::string> split_list;

    // first split the input list
    while (tmp.length() != 0)
      {
        std::string name;
        name = tmp;

        if (name.find(",") != std::string::npos)
          {
            name.erase (name.find(","), std::string::npos);
            tmp.erase (0, tmp.find(",")+1);
          }
        else
          tmp = "";

        while ((name.length() != 0) &&
               (std::isspace (name[0])))
          name.erase (0,1);
        while (std::isspace (name[name.length()-1]))
          name.erase (name.length()-1, 1);

        split_list.push_back (name);
      };


    // check the different possibilities
    for (std::list<std::string>::const_iterator test_string = split_list.begin();
         test_string != split_list.end(); ++test_string)
      {
        bool string_found = false;

        tmp = sequence;
        while (tmp.find('|') != std::string::npos)
          {
            if (*test_string == std::string(tmp, 0, tmp.find('|')))
              {
                // string found, quit
                // loop. don't change
                // tmp, since we don't
                // need it anymore.
                string_found = true;
                break;
              };

            tmp.erase (0, tmp.find('|')+1);
          };
        // check last choice, not finished by |
        if (!string_found)
          if (*test_string == tmp)
            string_found = true;

        if (!string_found)
          return false;
      };

    return true;
  }



  std::string MultipleSelection::description () const
  {
    std::ostringstream description;

    description << description_init
                << " "
                << sequence
                << " ]";

    return description.str();
  }



  PatternBase *
  MultipleSelection::clone () const
  {
    return new MultipleSelection(sequence);
  }


  std::size_t
  MultipleSelection::memory_consumption () const
  {
    return (sizeof(PatternBase) +
            MemoryConsumption::memory_consumption(sequence));
  }



  MultipleSelection *MultipleSelection::create (const std::string &description)
  {
    if (description.compare(0, std::strlen(description_init), description_init) == 0)
      {
        std::string sequence(description);

        sequence.erase(0, std::strlen(description_init) + 1);
        sequence.erase(sequence.length()-2, 2);

        return new MultipleSelection(sequence);
      }
    else
      return 0;
  }



  const char *Bool::description_init = "[Bool";


  Bool::Bool ()
    :
    Selection ("true|false")
  {}



  std::string Bool::description () const
  {
    std::ostringstream description;

    description << description_init
                << "]";

    return description.str();
  }



  PatternBase *
  Bool::clone () const
  {
    return new Bool();
  }



  Bool *Bool::create (const std::string &description)
  {
    if (description.compare(0, std::strlen(description_init), description_init) == 0)
      return new Bool();
    else
      return 0;
  }



  const char *Anything::description_init = "[Anything";


  Anything::Anything ()
  {}



  bool Anything::match (const std::string &) const
  {
    return true;
  }



  std::string Anything::description () const
  {
    std::ostringstream description;

    description << description_init
                << "]";

    return description.str();
  }



  PatternBase *
  Anything::clone () const
  {
    return new Anything();
  }



  Anything *Anything::create (const std::string &description)
  {
    if (description.compare(0, std::strlen(description_init), description_init) == 0)
      return new Anything();
    else
      return 0;
  }



  const char *FileName::description_init = "[FileName";


  FileName::FileName (const FileType type)
    : file_type (type)
  {}



  bool FileName::match (const std::string &) const
  {
    return true;
  }



  std::string FileName::description () const
  {
    std::ostringstream description;

    description << description_init;

    if (file_type == input)
      description << " (Type: input)]";
    else
      description << " (Type: output)]";

    return description.str();
  }



  PatternBase *
  FileName::clone () const
  {
    return new FileName(file_type);
  }



  FileName *FileName::create (const std::string &description)
  {
    if (description.compare(0, std::strlen(description_init), description_init) == 0)
      {
        std::istringstream is(description);
        std::string file_type;
        FileType type;

        is.ignore(strlen(description_init) + strlen(" (Type:"));

        is >> file_type;

        if (file_type == "input)]")
          type = input;
        else
          type = output;

        return new FileName(type);
      }
    else
      return 0;
  }



  const char *DirectoryName::description_init = "[DirectoryName";


  DirectoryName::DirectoryName ()
  {}



  bool DirectoryName::match (const std::string &) const
  {
    return true;
  }



  std::string DirectoryName::description () const
  {
    std::ostringstream description;

    description << description_init << "]";

    return description.str();
  }



  PatternBase *
  DirectoryName::clone () const
  {
    return new DirectoryName();
  }



  DirectoryName *DirectoryName::create (const std::string &description)
  {
    if (description.compare(0, std::strlen(description_init), description_init) == 0)
      return new DirectoryName();
    else
      return 0;
  }

}   // end namespace Patterns



ParameterHandler::ParameterHandler ()
  :
  entries (new boost::property_tree::ptree())
{}



ParameterHandler::~ParameterHandler ()
{}



std::string
ParameterHandler::mangle (const std::string &s)
{
  std::string u;

  // reserve the minimum number of characters we will need. it may
  // be more but this is the least we can do
  u.reserve (s.size());

  // see if the name is special and if so mangle the whole thing
  const bool mangle_whole_string = (s == "value");

  // for all parts of the string, see
  // if it is an allowed character or
  // not
  for (unsigned int i=0; i<s.size(); ++i)
    {
      static const std::string allowed_characters
      ("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789");

      if ((! mangle_whole_string)
          &&
          (allowed_characters.find (s[i]) != std::string::npos))
        u.push_back (s[i]);
      else
        {
          u.push_back ('_');
          static const char hex[16]
            = { '0','1','2','3','4','5','6','7','8','9','a','b','c','d','e','f'};
          u.push_back (hex[static_cast<unsigned char>(s[i])/16]);
          u.push_back (hex[static_cast<unsigned char>(s[i])%16]);
        }
    }

  return u;
}



std::string
ParameterHandler::demangle (const std::string &s)
{
  std::string u;
  u.reserve (s.size());

  for (unsigned int i=0; i<s.size(); ++i)
    if (s[i] != '_')
      u.push_back (s[i]);
    else
      {
        Assert (i+2 < s.size(),
                ExcMessage ("Trying to demangle an invalid string."));

        unsigned char c = 0;
        switch (s[i+1])
          {
          case '0':
            c = 0 * 16;
            break;
          case '1':
            c = 1 * 16;
            break;
          case '2':
            c = 2 * 16;
            break;
          case '3':
            c = 3 * 16;
            break;
          case '4':
            c = 4 * 16;
            break;
          case '5':
            c = 5 * 16;
            break;
          case '6':
            c = 6 * 16;
            break;
          case '7':
            c = 7 * 16;
            break;
          case '8':
            c = 8 * 16;
            break;
          case '9':
            c = 9 * 16;
            break;
          case 'a':
            c = 10 * 16;
            break;
          case 'b':
            c = 11 * 16;
            break;
          case 'c':
            c = 12 * 16;
            break;
          case 'd':
            c = 13 * 16;
            break;
          case 'e':
            c = 14 * 16;
            break;
          case 'f':
            c = 15 * 16;
            break;
          default:
            Assert (false, ExcInternalError());
          }
        switch (s[i+2])
          {
          case '0':
            c += 0;
            break;
          case '1':
            c += 1;
            break;
          case '2':
            c += 2;
            break;
          case '3':
            c += 3;
            break;
          case '4':
            c += 4;
            break;
          case '5':
            c += 5;
            break;
          case '6':
            c += 6;
            break;
          case '7':
            c += 7;
            break;
          case '8':
            c += 8;
            break;
          case '9':
            c += 9;
            break;
          case 'a':
            c += 10;
            break;
          case 'b':
            c += 11;
            break;
          case 'c':
            c += 12;
            break;
          case 'd':
            c += 13;
            break;
          case 'e':
            c += 14;
            break;
          case 'f':
            c += 15;
            break;
          default:
            Assert (false, ExcInternalError());
          }

        u.push_back (static_cast<char>(c));

        // skip the two characters
        i += 2;
      }

  return u;
}



namespace
{
  /**
   * Return whether a given node is a parameter node (as opposed
   * to being a subsection or alias node)
   */
  bool
  is_parameter_node (const boost::property_tree::ptree &p)
  {
    return static_cast<bool>(p.get_optional<std::string>("value"));
  }


  /**
   * Return whether a given node is a alias node (as opposed
   * to being a subsection or parameter node)
   */
  bool
  is_alias_node (const boost::property_tree::ptree &p)
  {
    return static_cast<bool>(p.get_optional<std::string>("alias"));
  }
}


std::string
ParameterHandler::get_current_path () const
{
  if (subsection_path.size() > 0)
    {
      std::string p = mangle(subsection_path[0]);
      for (unsigned int i=1; i<subsection_path.size(); ++i)
        {
          p += path_separator;
          p += mangle(subsection_path[i]);
        }
      return p;
    }
  else
    return "";
}



std::string
ParameterHandler::get_current_full_path (const std::string &name) const
{
  std::string path = get_current_path ();
  if (path.empty() == false)
    path += path_separator;

  path += mangle(name);

  return path;
}



bool ParameterHandler::read_input (std::istream &input,
                                   const std::string &filename,
                                   const std::string &last_line)
{
  AssertThrow (input, ExcIO());

  // store subsections we are currently in
  std::vector<std::string> saved_path = subsection_path;

  std::string input_line;
  std::string fully_concatenated_line;
  bool is_concatenated = false;
  unsigned int current_line_n = 0;
  bool status = true;

  while (std::getline (input, input_line))
    {
      ++current_line_n;
      // Trim the whitespace at the ends of the line here instead of in
      // scan_line. This makes the continuation line logic a lot simpler.
      input_line = Utilities::trim (input_line);

      // If we see the line which is the same as @p last_line ,
      // terminate the parsing.
      if (last_line.length() != 0 &&
          input_line == last_line)
        break;

      // Check whether or not the current line should be joined with the next
      // line before calling scan_line.
      if (input_line.length() != 0 &&
          input_line.find_last_of('\\') == input_line.length() - 1)
        {
          input_line.erase (input_line.length() - 1); // remove the last '\'
          is_concatenated = true;

          fully_concatenated_line += input_line;
        }
      // If the previous line ended in a '\' but the current did not, then we
      // should proceed to scan_line.
      else if (is_concatenated)
        {
          fully_concatenated_line += input_line;
          is_concatenated = false;
        }
      // Finally, if neither the previous nor current lines are continuations,
      // then the current input line is entirely concatenated.
      else
        {
          fully_concatenated_line = input_line;
        }

      if (!is_concatenated)
        {
          status &= scan_line (fully_concatenated_line, filename, current_line_n);
          fully_concatenated_line.clear();
        }
    }

  // While it does not make much sense for anyone to actually do this, allow
  // the last line to end in a backslash.
  if (is_concatenated)
    {
      status &= scan_line (fully_concatenated_line, filename, current_line_n);
    }

  if (status && (saved_path != subsection_path))
    {
      std::cerr << "Unbalanced 'subsection'/'end' in file <" << filename
                << ">." << std::endl;
      if (saved_path.size()>0)
        {
          std::cerr << "Path before loading input:" << std::endl;
          for (unsigned int i=0; i<saved_path.size(); ++i)
            std::cerr << std::setw(i*2+4) << " "
                      << "subsection " << saved_path[i] << std::endl;
        }
      std::cerr << "Current path:" << std::endl;
      for (unsigned int i=0; i<subsection_path.size(); ++i)
        std::cerr << std::setw(i*2+4) << " "
                  << "subsection " << subsection_path[i] << std::endl;

      // restore subsection we started with and return failure:
      subsection_path = saved_path;
      return false;
    }

  return status;
}



bool ParameterHandler::read_input (const std::string &filename,
                                   const bool optional,
                                   const bool write_compact,
                                   const std::string &last_line)
{
  PathSearch search("PARAMETERS");

  try
    {
      std::string openname = search.find(filename);
      std::ifstream file_stream (openname.c_str());
      AssertThrow(file_stream, ExcIO());

      return read_input (file_stream, filename, last_line);
    }
  catch (const PathSearch::ExcFileNotFound &)
    {
      std::cerr << "ParameterHandler::read_input: could not open file <"
                << filename << "> for reading." << std::endl;
      if (!optional)
        {
          std:: cerr << "Trying to make file <"
                     << filename << "> with default values for you." << std::endl;
          std::ofstream output (filename.c_str());
          if (output)
            print_parameters (output, (write_compact ? ShortText : Text));
        }
    }
  return false;
}



bool ParameterHandler::read_input_from_string (const char *s,
                                               const std::string &last_line)
{
  std::istringstream input_stream (s);
  return read_input (input_stream, "input string", last_line);
}



namespace
{
  // Recursively go through the 'source' tree
  // and see if we can find corresponding
  // entries in the 'destination' tree. If
  // not, error out (i.e. we have just read
  // an XML file that has entries that
  // weren't declared in the ParameterHandler
  // object); if so, copy the value of these
  // nodes into the destination object
  bool
  read_xml_recursively (const boost::property_tree::ptree &source,
                        const std::string                 &current_path,
                        const char                         path_separator,
                        const std::vector<std_cxx11::shared_ptr<const Patterns::PatternBase> > &
                        patterns,
                        boost::property_tree::ptree       &destination)
  {
    for (boost::property_tree::ptree::const_iterator p = source.begin();
         p != source.end(); ++p)
      {
        // a sub-tree must either be a
        // parameter node or a subsection
        if (p->second.get_optional<std::string>("value"))
          {
            // make sure we have a
            // corresponding entry in the
            // destination object as well
            const std::string full_path
              = (current_path == ""
                 ?
                 p->first
                 :
                 current_path + path_separator + p->first);
            if (destination.get_optional<std::string> (full_path)
                &&
                destination.get_optional<std::string> (full_path +
                                                       path_separator +
                                                       "value"))
              {
                // first make sure that the
                // new entry actually
                // satisfies its constraints
                const std::string new_value
                  = p->second.get<std::string>("value");

                const unsigned int pattern_index
                  = destination.get<unsigned int> (full_path +
                                                   path_separator +
                                                   "pattern");
                if (patterns[pattern_index]->match(new_value) == false)
                  {
                    std::cerr << "    The entry value" << std::endl
                              << "        " << new_value << std::endl
                              << "    for the entry named" << std::endl
                              << "        " << full_path << std::endl
                              << "    does not match the given pattern" << std::endl
                              << "        " << patterns[pattern_index]->description()
                              << std::endl;
                    return false;
                  }

                // set the found parameter in
                // the destination argument
                destination.put (full_path + path_separator + "value",
                                 new_value);

                // this node might have
                // sub-nodes in addition to
                // "value", such as
                // "default_value",
                // "documentation", etc. we
                // might at some point in the
                // future want to make sure
                // that if they exist that
                // they match the ones in the
                // 'destination' tree
              }
            else
              {
                std::cerr << "The entry <" << full_path
                          << "> with value <"
                          << p->second.get<std::string>("value")
                          << "> has not been declared."
                          << std::endl;
                return false;
              }
          }
        else if (p->second.get_optional<std::string>("alias"))
          {
            // it is an alias node. alias nodes are static and
            // there is nothing to do here (but the same applies as
            // mentioned in the comment above about the static
            // nodes inside parameter nodes
          }
        else
          {
            // it must be a subsection
            const bool result
              = read_xml_recursively (p->second,
                                      (current_path == "" ?
                                       p->first :
                                       current_path + path_separator + p->first),
                                      path_separator,
                                      patterns,
                                      destination);

            // see if the recursive read
            // succeeded. if yes, continue,
            // otherwise exit now
            if (result == false)
              return false;
          }
      }

    return true;
  }
}



bool ParameterHandler::read_input_from_xml (std::istream &in)
{
  // read the XML tree assuming that (as we
  // do in print_parameters(XML) it has only
  // a single top-level node called
  // "ParameterHandler"
  boost::property_tree::ptree single_node_tree;
  try
    {
      read_xml (in, single_node_tree);
    }
  catch (...)
    {
      std::cerr << "This input stream appears not to be valid XML"
                << std::endl;
      return false;
    }

  // make sure there is a top-level element
  // called "ParameterHandler"
  if (!single_node_tree.get_optional<std::string>("ParameterHandler"))
    {
      std::cerr << "There is no top-level XML element called \"ParameterHandler\"."
                << std::endl;
      return false;
    }

  // ensure that there is only a single
  // top-level element
  if (std::distance (single_node_tree.begin(), single_node_tree.end()) != 1)
    {
      std::cerr << "The top-level XML element \"ParameterHandler\" is "
                << "not the only one."
                << std::endl;
      std::cerr << "(There are "
                << std::distance (single_node_tree.begin(),
                                  single_node_tree.end())
                << " top-level elements.)"
                << std::endl;
      return false;
    }

  // read the child elements recursively
  const boost::property_tree::ptree
  &my_entries = single_node_tree.get_child("ParameterHandler");

  return read_xml_recursively (my_entries, "", path_separator, patterns,
                               *entries);
}



void ParameterHandler::clear ()
{
  entries.reset (new boost::property_tree::ptree());
}



void
ParameterHandler::declare_entry (const std::string           &entry,
                                 const std::string           &default_value,
                                 const Patterns::PatternBase &pattern,
                                 const std::string           &documentation)
{
  entries->put (get_current_full_path(entry) + path_separator + "value",
                default_value);
  entries->put (get_current_full_path(entry) + path_separator + "default_value",
                default_value);
  entries->put (get_current_full_path(entry) + path_separator + "documentation",
                documentation);

  // clone the pattern and store its
  // index in the node
  patterns.push_back (std_cxx11::shared_ptr<const Patterns::PatternBase>
                      (pattern.clone()));
  entries->put (get_current_full_path(entry) + path_separator + "pattern",
                static_cast<unsigned int>(patterns.size()-1));
  // also store the description of
  // the pattern. we do so because we
  // may wish to export the whole
  // thing as XML or any other format
  // so that external tools can work
  // on the parameter file; in that
  // case, they will have to be able
  // to re-create the patterns as far
  // as possible
  entries->put (get_current_full_path(entry) + path_separator +
                "pattern_description",
                patterns.back()->description());

  // as documented, do the default value checking at the very end
  AssertThrow (pattern.match (default_value),
               ExcValueDoesNotMatchPattern (default_value, pattern.description()));
}



void
ParameterHandler::declare_alias(const std::string &existing_entry_name,
                                const std::string &alias_name,
                                const bool         alias_is_deprecated)
{
  // see if there is anything to refer to already
  Assert (entries->get_optional<std::string>(get_current_full_path(existing_entry_name)),
          ExcMessage ("You are trying to declare an alias entry <"
                      + alias_name +
                      "> that references an entry <"
                      + existing_entry_name +
                      ">, but the latter does not exist."));
  // then also make sure that what is being referred to is in
  // fact a parameter (not an alias or subsection)
  Assert (entries->get_optional<std::string>(get_current_full_path(existing_entry_name) + path_separator + "value"),
          ExcMessage ("You are trying to declare an alias entry <"
                      + alias_name +
                      "> that references an entry <"
                      + existing_entry_name +
                      ">, but the latter does not seem to be a "
                      "parameter declaration."));


  // now also make sure that if the alias has already been
  // declared, that it is also an alias and refers to the same
  // entry
  if (entries->get_optional<std::string>(get_current_full_path(alias_name)))
    {
      Assert (entries->get_optional<std::string> (get_current_full_path(alias_name) + path_separator + "alias"),
              ExcMessage ("You are trying to declare an alias entry <"
                          + alias_name +
                          "> but a non-alias entry already exists in this "
                          "subsection (i.e., there is either a preexisting "
                          "further subsection, or a parameter entry, with "
                          "the same name as the alias)."));
      Assert (entries->get<std::string> (get_current_full_path(alias_name) + path_separator + "alias")
              ==
              existing_entry_name,
              ExcMessage ("You are trying to declare an alias entry <"
                          + alias_name +
                          "> but an alias entry already exists in this "
                          "subsection and this existing alias references a "
                          "different parameter entry. Specifically, "
                          "you are trying to reference the entry <"
                          + existing_entry_name +
                          "> whereas the existing alias references "
                          "the entry <"
                          + entries->get<std::string> (get_current_full_path(alias_name) + path_separator + "alias") +
                          ">."));
    }

  entries->put (get_current_full_path(alias_name) + path_separator + "alias",
                existing_entry_name);
  entries->put (get_current_full_path(alias_name) + path_separator + "deprecation_status",
                (alias_is_deprecated ? "true" : "false"));
}



void ParameterHandler::enter_subsection (const std::string &subsection)
{
  const std::string current_path = get_current_path ();

  // if necessary create subsection
  if (!entries->get_child_optional (get_current_full_path(subsection)))
    entries->add_child (get_current_full_path(subsection),
                        boost::property_tree::ptree());

  // then enter it
  subsection_path.push_back (subsection);
}



void ParameterHandler::leave_subsection ()
{
  // assert there is a subsection that
  // we may leave
  Assert (subsection_path.size() != 0, ExcAlreadyAtTopLevel());

  if (subsection_path.size() > 0)
    subsection_path.pop_back ();
}



std::string
ParameterHandler::get (const std::string &entry_string) const
{
  // assert that the entry is indeed
  // declared
  if (boost::optional<std::string> value
      = entries->get_optional<std::string> (get_current_full_path(entry_string) + path_separator + "value"))
    return value.get();
  else
    {
      Assert (false, ExcEntryUndeclared(entry_string));
      return "";
    }
}



long int ParameterHandler::get_integer (const std::string &entry_string) const
{
  try
    {
      return Utilities::string_to_int (get (entry_string));
    }
  catch (...)
    {
      AssertThrow (false,
                   ExcMessage("Can't convert the parameter value <"
                              + get(entry_string) +
                              "> for entry <"
                              + entry_string +
                              " to an integer."));
      return 0;
    }
}



double ParameterHandler::get_double (const std::string &entry_string) const
{
  try
    {
      return Utilities::string_to_double (get (entry_string));
    }
  catch (...)
    {
      AssertThrow (false,
                   ExcMessage("Can't convert the parameter value <"
                              + get(entry_string) +
                              "> for entry <"
                              + entry_string +
                              " to a double precision variable."));
      return 0;
    }
}



bool ParameterHandler::get_bool (const std::string &entry_string) const
{
  const std::string s = get(entry_string);

  AssertThrow ((s=="true") || (s=="false") ||
               (s=="yes") || (s=="no"),
               ExcMessage("Can't convert the parameter value <"
                          + get(entry_string) +
                          "> for entry <"
                          + entry_string +
                          " to a boolean."));
  if (s=="true" || s=="yes")
    return true;
  else
    return false;
}



void
ParameterHandler::set (const std::string &entry_string,
                       const std::string &new_value)
{
  // resolve aliases before looking up the correct entry
  std::string path = get_current_full_path(entry_string);
  if (entries->get_optional<std::string>(path + path_separator + "alias"))
    path = get_current_full_path(entries->get<std::string>(path + path_separator + "alias"));

  // assert that the entry is indeed declared
  if (entries->get_optional<std::string>(path + path_separator + "value"))
    {
      const unsigned int pattern_index
        = entries->get<unsigned int> (path + path_separator + "pattern");
      AssertThrow (patterns[pattern_index]->match(new_value),
                   ExcValueDoesNotMatchPattern (new_value,
                                                entries->get<std::string>
                                                (path +
                                                 path_separator +
                                                 "pattern_description")));

      entries->put (path + path_separator + "value",
                    new_value);
    }
  else
    AssertThrow (false, ExcEntryUndeclared(entry_string));
}


void
ParameterHandler::set (const std::string &entry_string,
                       const char        *new_value)
{
  // simply forward
  set (entry_string, std::string(new_value));
}


void
ParameterHandler::set (const std::string &entry_string,
                       const double      &new_value)
{
  std::ostringstream s;
  s << std::setprecision(16);
  s << new_value;

  // hand this off to the function that
  // actually sets the value as a string
  set (entry_string, s.str());
}



void
ParameterHandler::set (const std::string &entry_string,
                       const long int    &new_value)
{
  std::ostringstream s;
  s << new_value;

  // hand this off to the function that
  // actually sets the value as a string
  set (entry_string, s.str());
}



void
ParameterHandler::set (const std::string &entry_string,
                       const bool        &new_value)
{
  // hand this off to the function that
  // actually sets the value as a string
  set (entry_string,
       (new_value ? "true" : "false"));
}



std::ostream &
ParameterHandler::print_parameters (std::ostream     &out,
                                    const OutputStyle style)
{
  AssertThrow (out, ExcIO());

  switch (style)
    {
    case XML:
    {
      // call the writer
      // function and exit as
      // there is nothing
      // further to do down in
      // this function
      //
      // XML has a requirement that
      // there can only be one
      // single top-level entry,
      // but we may have multiple
      // entries and sections.  we
      // work around this by
      // creating a tree just for
      // this purpose with the
      // single top-level node
      // "ParameterHandler" and
      // assign the existing tree
      // under it
      boost::property_tree::ptree single_node_tree;
      single_node_tree.add_child("ParameterHandler",
                                 *entries);

      write_xml (out, single_node_tree);
      return out;
    }


    case JSON:
      // call the writer
      // function and exit as
      // there is nothing
      // further to do down in
      // this function
      write_json (out, *entries);
      return out;

    case Text:
      out << "# Listing of Parameters" << std::endl
          << "# ---------------------" << std::endl;
      break;
    case LaTeX:
      out << "\\subsection{Global parameters}" << std::endl;
      out << "\\label{parameters:global}" << std::endl;
      out << std::endl << std::endl;
      break;
    case Description:
      out << "Listing of Parameters:" << std::endl << std::endl;
      break;
    case ShortText:
      break;
    default:
      Assert (false, ExcNotImplemented());
    };

  // dive recursively into the subsections
  print_parameters_section (out, style, 0);

  switch (style)
    {
    case Text:
    case Description:
    case ShortText:
      break;
    case LaTeX:
      break;
    default:
      Assert (false, ExcNotImplemented());
    };

  return out;
}



// Print a section in the desired style. The styles are separated into
// several verbosity classes depending on the higher bits.
//
// If bit 7 (128) is set, comments are not printed.
// If bit 6 (64) is set, default values after change are not printed.
void
ParameterHandler::print_parameters_section (std::ostream      &out,
                                            const OutputStyle  style,
                                            const unsigned int indent_level,
                                            const bool         include_top_level_elements)
{
  AssertThrow (out, ExcIO());

  const boost::property_tree::ptree &current_section
    = entries->get_child (get_current_path());

  unsigned int overall_indent_level = indent_level;

  switch (style)
    {
    case XML:
    {
      if (include_top_level_elements)
        {
          // call the writer
          // function and exit as
          // there is nothing
          // further to do down in
          // this function
          //
          // XML has a requirement that
          // there can only be one
          // single top-level entry,
          // but a section has multiple
          // entries and sections. we
          // work around this by
          // creating a tree just for
          // this purpose with the
          // single top-level node
          // "ParameterHandler" and
          // assign the full path of
          // down to the current section
          // under it
          boost::property_tree::ptree single_node_tree;

          // if there is no subsection selected,
          // add the whole tree of entries,
          // otherwise add a root element
          // and the selected subsection under it
          if (subsection_path.size() == 0)
            {
              single_node_tree.add_child("ParameterHandler",
                                         *entries);
            }
          else
            {
              std::string  path ("ParameterHandler");

              single_node_tree.add_child(path,
                                         boost::property_tree::ptree());

              path += path_separator + get_current_path ();
              single_node_tree.add_child (path, current_section);
            };

          write_xml (out, single_node_tree);
        }
      else
        Assert (false, ExcNotImplemented());

      break;
    }
    case Text:
    case ShortText:
    {
      // if there are top level elements to print, do it
      if (include_top_level_elements && (subsection_path.size() > 0))
        for (unsigned int i=0; i<subsection_path.size(); ++i)
          {
            out << std::setw(overall_indent_level*2) << ""
                << "subsection " << demangle (subsection_path[i]) << std::endl;
            overall_indent_level += 1;
          };

      // first find out the longest
      // entry name to be able to
      // align the equal signs
      //
      // to do this loop over all
      // nodes of the current tree,
      // select the parameter nodes
      // (and discard sub-tree
      // nodes) and take the
      // maximum of their lengths
      std::size_t longest_name = 0;
      for (boost::property_tree::ptree::const_iterator
           p = current_section.begin();
           p != current_section.end(); ++p)
        if (is_parameter_node (p->second) == true)
          longest_name = std::max (longest_name,
                                   demangle(p->first).length());

      // likewise find the longest
      // actual value string to
      // make sure we can align the
      // default and documentation
      // strings
      std::size_t longest_value = 0;
      for (boost::property_tree::ptree::const_iterator
           p = current_section.begin();
           p != current_section.end(); ++p)
        if (is_parameter_node (p->second) == true)
          longest_value = std::max (longest_value,
                                    p->second.get<std::string>("value").length());


      // print entries one by
      // one. make sure they are
      // sorted by using the
      // appropriate iterators
      bool first_entry = true;
      for (boost::property_tree::ptree::const_assoc_iterator
           p = current_section.ordered_begin();
           p != current_section.not_found(); ++p)
        if (is_parameter_node (p->second) == true)
          {
            const std::string value = p->second.get<std::string>("value");

            // if there is documentation,
            // then add an empty line (unless
            // this is the first entry in a
            // subsection), print the
            // documentation, and then the
            // actual entry; break the
            // documentation into readable
            // chunks such that the whole
            // thing is at most 78 characters
            // wide
            if ((!(style & 128)) &&
                !p->second.get<std::string>("documentation").empty())
              {
                if (first_entry == false)
                  out << std::endl;
                else
                  first_entry = false;

                const std::vector<std::string> doc_lines
                  = Utilities::
                    break_text_into_lines (p->second.get<std::string>("documentation"),
                                           78 - overall_indent_level*2 - 2);

                for (unsigned int i=0; i<doc_lines.size(); ++i)
                  out << std::setw(overall_indent_level*2) << ""
                      << "# "
                      << doc_lines[i]
                      << std::endl;
              }



            // print name and value
            // of this entry
            out << std::setw(overall_indent_level*2) << ""
                << "set "
                << demangle(p->first)
                << std::setw(longest_name-demangle(p->first).length()+1) << " "
                << "= " << value;

            // finally print the
            // default value, but
            // only if it differs
            // from the actual value
            if ((!(style & 64)) && value != p->second.get<std::string>("default_value"))
              {
                out << std::setw(longest_value-value.length()+1) << ' '
                    << "# ";
                out << "default: " << p->second.get<std::string>("default_value");
              }

            out << std::endl;
          }

      break;
    }

    case LaTeX:
    {
      // if there are any parameters in
      // this section then print them as an
      // itemized list
      bool parameters_exist_here = false;
      for (boost::property_tree::ptree::const_assoc_iterator
           p = current_section.ordered_begin();
           p != current_section.not_found(); ++p)
        if ((is_parameter_node (p->second) == true)
            ||
            (is_alias_node (p->second) == true))
          {
            parameters_exist_here = true;
            break;
          }

      if (parameters_exist_here)
        {
          out << "\\begin{itemize}"
              << std::endl;

          // print entries one by
          // one. make sure they are
          // sorted by using the
          // appropriate iterators
          for (boost::property_tree::ptree::const_assoc_iterator
               p = current_section.ordered_begin();
               p != current_section.not_found(); ++p)
            if (is_parameter_node (p->second) == true)
              {
                const std::string value = p->second.get<std::string>("value");

                // print name
                out << "\\item {\\it Parameter name:} {\\tt " << demangle(p->first) << "}\n"
                    << "\\phantomsection\\label{parameters:";
                for (unsigned int i=0; i<subsection_path.size(); ++i)
                  out << subsection_path[i] << "/";
                out << demangle(p->first);
                out << "}\n\n"
                    << std::endl;

                out << "\\index[prmindex]{"
                    << demangle(p->first)
                    << "}\n";
                out << "\\index[prmindexfull]{";
                for (unsigned int i=0; i<subsection_path.size(); ++i)
                  out << subsection_path[i] << "!";
                out << demangle(p->first)
                    << "}\n";

                // finally print value and default
                out << "{\\it Value:} " << value << "\n\n"
                    << std::endl
                    << "{\\it Default:} "
                    << p->second.get<std::string>("default_value") << "\n\n"
                    << std::endl;

                // if there is a
                // documenting string,
                // print it as well
                if (!p->second.get<std::string>("documentation").empty())
                  out << "{\\it Description:} "
                      << p->second.get<std::string>("documentation") << "\n\n"
                      << std::endl;

                // also output possible values
                out << "{\\it Possible values:} "
                    << p->second.get<std::string> ("pattern_description")
                    << std::endl;
              }
            else if (is_alias_node (p->second) == true)
              {
                const std::string alias = p->second.get<std::string>("alias");

                // print name
                out << "\\item {\\it Parameter name:} {\\tt " << demangle(p->first) << "}\n"
                    << "\\phantomsection\\label{parameters:";
                for (unsigned int i=0; i<subsection_path.size(); ++i)
                  out << subsection_path[i] << "/";
                out << demangle(p->first);
                out << "}\n\n"
                    << std::endl;

                out << "\\index[prmindex]{"
                    << demangle(p->first)
                    << "}\n";
                out << "\\index[prmindexfull]{";
                for (unsigned int i=0; i<subsection_path.size(); ++i)
                  out << subsection_path[i] << "!";
                out << demangle(p->first)
                    << "}\n";

                // finally print alias and indicate if it is deprecated
                out << "This parameter is an alias for the parameter ``\\texttt{"
                    << alias << "}''."
                    << (p->second.get<std::string>("deprecation_status") == "true"
                        ?
                        " Its use is deprecated."
                        :
                        "")
                    << "\n\n"
                    << std::endl;
              }
          out << "\\end{itemize}" << std::endl;
        }

      break;
    }

    case Description:
    {
      // if there are top level elements to print, do it
      if (include_top_level_elements && (subsection_path.size() > 0))
        for (unsigned int i=0; i<subsection_path.size(); ++i)
          {
            out << std::setw(overall_indent_level*2) << ""
                << "subsection " << demangle (subsection_path[i]) << std::endl;
            overall_indent_level += 1;
          };

      // first find out the longest
      // entry name to be able to
      // align the equal signs
      std::size_t longest_name = 0;
      for (boost::property_tree::ptree::const_iterator
           p = current_section.begin();
           p != current_section.end(); ++p)
        if (is_parameter_node (p->second) == true)
          longest_name = std::max (longest_name,
                                   demangle(p->first).length());

      // print entries one by
      // one. make sure they are
      // sorted by using the
      // appropriate iterators
      for (boost::property_tree::ptree::const_assoc_iterator
           p = current_section.ordered_begin();
           p != current_section.not_found(); ++p)
        if (is_parameter_node (p->second) == true)
          {
            const std::string value = p->second.get<std::string>("value");

            // print name and value
            out << std::setw(overall_indent_level*2) << ""
                << "set "
                << demangle(p->first)
                << std::setw(longest_name-demangle(p->first).length()+1) << " "
                << " = ";

            // print possible values:
            const std::vector<std::string> description_str
              = Utilities::break_text_into_lines (p->second.get<std::string>
                                                  ("pattern_description"),
                                                  78 - overall_indent_level*2 - 2, '|');
            if (description_str.size() > 1)
              {
                out << std::endl;
                for (unsigned int i=0; i<description_str.size(); ++i)
                  out << std::setw(overall_indent_level*2+6) << ""
                      << description_str[i] << std::endl;
              }
            else if (description_str.empty() == false)
              out << "  " << description_str[0] << std::endl;
            else
              out << std::endl;

            // if there is a
            // documenting string,
            // print it as well
            if (p->second.get<std::string>("documentation").length() != 0)
              out << std::setw(overall_indent_level*2 + longest_name + 10) << ""
                  << "(" << p->second.get<std::string>("documentation") << ")" << std::endl;
          }

      break;
    }

    default:
      Assert (false, ExcNotImplemented());
    }


  // if there was text before and there are
  // sections to come, put two newlines
  // between the last entry and the first
  // subsection
  if (style != XML)
    {
      unsigned int n_parameters = 0;
      unsigned int n_sections   = 0;
      for (boost::property_tree::ptree::const_iterator
           p = current_section.begin();
           p != current_section.end(); ++p)
        if (is_parameter_node (p->second) == true)
          ++n_parameters;
        else if (is_alias_node (p->second) == false)
          ++n_sections;

      if ((style != Description)
          &&
          (!(style & 128))
          &&
          (n_parameters != 0)
          &&
          (n_sections != 0))
        out << std::endl << std::endl;

      // now traverse subsections tree,
      // in alphabetical order
      for (boost::property_tree::ptree::const_assoc_iterator
           p = current_section.ordered_begin();
           p != current_section.not_found(); ++p)
        if ((is_parameter_node (p->second) == false)
            &&
            (is_alias_node (p->second) == false))
          {
            // first print the subsection header
            switch (style)
              {
              case Text:
              case Description:
              case ShortText:
                out << std::setw(overall_indent_level*2) << ""
                    << "subsection " << demangle(p->first) << std::endl;
                break;
              case LaTeX:
              {
                out << std::endl
                    << "\\subsection{Parameters in section \\tt ";

                // find the path to the
                // current section so that we
                // can print it in the
                // \subsection{...} heading
                for (unsigned int i=0; i<subsection_path.size(); ++i)
                  out << subsection_path[i] << "/";
                out << demangle(p->first);

                out << "}" << std::endl;
                out << "\\label{parameters:";
                for (unsigned int i=0; i<subsection_path.size(); ++i)
                  out << mangle(subsection_path[i]) << "/";
                out << p->first << "}";
                out << std::endl;

                out << std::endl;
                break;
              }

              default:
                Assert (false, ExcNotImplemented());
              };

            // then the contents of the
            // subsection
            enter_subsection (demangle(p->first));
            print_parameters_section (out, style, overall_indent_level+1);
            leave_subsection ();
            switch (style)
              {
              case Text:
                // write end of
                // subsection. one
                // blank line after
                // each subsection
                out << std::setw(overall_indent_level*2) << ""
                    << "end" << std::endl
                    << std::endl;

                // if this is a toplevel
                // subsection, then have two
                // newlines
                if (overall_indent_level == 0)
                  out << std::endl;

                break;
              case Description:
                break;
              case ShortText:
                // write end of
                // subsection.
                out << std::setw(overall_indent_level*2) << ""
                    << "end" << std::endl;
                break;
              case LaTeX:
                break;
              default:
                Assert (false, ExcNotImplemented());
              }
          }
    }

  // close top level elements, if there are any
  switch (style)
    {
    case XML:
    case LaTeX:
    case Description:
      break;
    case Text:
    case ShortText:
    {
      if (include_top_level_elements && (subsection_path.size() > 0))
        for (unsigned int i=0; i<subsection_path.size(); ++i)
          {
            overall_indent_level -= 1;
            out << std::setw(overall_indent_level*2) << ""
                << "end" << std::endl;
          };

      break;
    }

    default:
      Assert (false, ExcNotImplemented());
    }

}



void
ParameterHandler::log_parameters (LogStream &out)
{
  out.push("parameters");
  // dive recursively into the
  // subsections
  log_parameters_section (out);

  out.pop();
}



void
ParameterHandler::log_parameters_section (LogStream &out)
{
  const boost::property_tree::ptree &current_section
    = entries->get_child (get_current_path());

  // print entries one by
  // one. make sure they are
  // sorted by using the
  // appropriate iterators
  for (boost::property_tree::ptree::const_assoc_iterator
       p = current_section.ordered_begin();
       p != current_section.not_found(); ++p)
    if (is_parameter_node (p->second) == true)
      out << demangle(p->first) << ": "
          << p->second.get<std::string>("value") << std::endl;

  // now transverse subsections tree
  // now traverse subsections tree,
  // in alphabetical order
  for (boost::property_tree::ptree::const_assoc_iterator
       p = current_section.ordered_begin();
       p != current_section.not_found(); ++p)
    if (is_parameter_node (p->second) == false)
      {
        out.push (demangle(p->first));
        enter_subsection (demangle(p->first));
        log_parameters_section (out);
        leave_subsection ();
        out.pop ();
      }
}



bool
ParameterHandler::scan_line (std::string         line,
                             const std::string  &input_filename,
                             const unsigned int  current_line_n)
{
  // if there is a comment, delete it
  if (line.find('#') != std::string::npos)
    line.erase (line.find("#"), std::string::npos);

  // replace \t by space:
  while (line.find('\t') != std::string::npos)
    line.replace (line.find('\t'), 1, " ");

  //trim start and end:
  line = Utilities::trim(line);

  // if line is now empty: leave
  if (line.length() == 0)
    return true;

  // enter subsection
  if ((line.find ("SUBSECTION ") == 0) ||
      (line.find ("subsection ") == 0))
    {
      // delete this prefix
      line.erase (0, std::string("subsection").length()+1);

      const std::string subsection = Utilities::trim(line);

      // check whether subsection exists
      if (!entries->get_child_optional (get_current_full_path(subsection)))
        {
          std::cerr << "Line <" << current_line_n
                    << "> of file <" << input_filename
                    << ">: There is no such subsection to be entered: "
                    << demangle(get_current_full_path(subsection)) << std::endl;
          for (unsigned int i=0; i<subsection_path.size(); ++i)
            std::cerr << std::setw(i*2+4) << " "
                      << "subsection " << subsection_path[i] << std::endl;
          std::cerr << std::setw(subsection_path.size()*2+4) << " "
                    << "subsection " << subsection << std::endl;
          return false;
        }

      // subsection exists
      subsection_path.push_back (subsection);
      return true;
    }

  // exit subsection
  if ((line.find ("END") == 0) ||
      (line.find ("end") == 0))
    {
      line.erase (0, 3);
      while ((line.size() > 0) && (std::isspace(line[0])))
        line.erase (0, 1);

      if (line.size()>0)
        {
          std::cerr << "Line <" << current_line_n
                    << "> of file <" << input_filename
                    << ">: invalid content after 'end'!" << std::endl;
          return false;
        }

      if (subsection_path.size() == 0)
        {
          std::cerr << "Line <" << current_line_n
                    << "> of file <" << input_filename
                    << ">: There is no subsection to leave here!" << std::endl;
          return false;
        }
      else
        {
          leave_subsection ();
          return true;
        }

    }

  // regular entry
  if ((line.find ("SET ") == 0) ||
      (line.find ("set ") == 0))
    {
      // erase "set" statement
      line.erase (0, 4);

      std::string::size_type pos = line.find("=");
      if (pos == std::string::npos)
        {
          std::cerr << "Line <" << current_line_n
                    << "> of file <" << input_filename
                    << ">: invalid format of set expression!" << std::endl;
          return false;
        }

      // extract entry name and value and trim
      std::string entry_name = Utilities::trim(std::string(line, 0, pos));
      std::string entry_value = Utilities::trim(std::string(line, pos+1, std::string::npos));

      // resolve aliases before we look up the entry. if necessary, print
      // a warning that the alias is deprecated
      std::string path = get_current_full_path(entry_name);
      if (entries->get_optional<std::string>(path + path_separator + "alias"))
        {
          if (entries->get<std::string>(path + path_separator + "deprecation_status") == "true")
            {
              std::cerr << "Warning in line <" << current_line_n
                        << "> of file <" << input_filename
                        << ">: You are using the deprecated spelling <"
                        << entry_name
                        << "> of the parameter <"
                        << entries->get<std::string>(path + path_separator + "alias")
                        << ">." << std::endl;
            }
          path = get_current_full_path(entries->get<std::string>(path + path_separator + "alias"));
        }

      // assert that the entry is indeed declared
      if (entries->get_optional<std::string> (path + path_separator + "value"))
        {
          // if entry was declared:
          // does it match the regex? if not,
          // don't enter it into the database
          // exception: if it contains characters
          // which specify it as a multiple loop
          // entry, then ignore content
          if (entry_value.find ('{') == std::string::npos)
            {
              const unsigned int pattern_index
                = entries->get<unsigned int> (path + path_separator + "pattern");
              if (!patterns[pattern_index]->match(entry_value))
                {
                  std::cerr << "Line <" << current_line_n
                            << "> of file <" << input_filename
                            << ">:" << std::endl
                            << "    The entry value" << std::endl
                            << "        " << entry_value << std::endl
                            << "    for the entry named" << std::endl
                            << "        " << entry_name << std::endl
                            << "    does not match the given pattern" << std::endl
                            << "        " << patterns[pattern_index]->description()
                            << std::endl;
                  return false;
                }
            }

          entries->put (path + path_separator + "value",
                        entry_value);
          return true;
        }
      else
        {
          std::cerr << "Line <" << current_line_n
                    << "> of file <" << input_filename
                    << ">: No such entry was declared:" << std::endl
                    << "    " << entry_name << std::endl
                    << "    <Present subsection:" << std::endl;
          for (unsigned int i=0; i<subsection_path.size(); ++i)
            std::cerr << std::setw(i*2+8) << " "
                      << "subsection " << subsection_path[i] << std::endl;
          std::cerr << "    >" << std::endl;

          return false;
        }
    }

  // an include statement?
  if ((line.find ("INCLUDE ") == 0) ||
      (line.find ("include ") == 0))
    {
      // erase "include " statement and eliminate spaces
      line.erase (0, 7);
      while ((line.size() > 0) && (line[0] == ' '))
        line.erase (0, 1);

      // the remainder must then be a filename
      if (line.size() == 0)
        {
          std::cerr << "Line <" << current_line_n
                    << "> of file <" << input_filename
                    << "> is an include statement but does not name a file!"
                    << std::endl;

          return false;
        }

      std::ifstream input (line.c_str());
      if (!input)
        {
          std::cerr << "Line <" << current_line_n
                    << "> of file <" << input_filename
                    << "> is an include statement but the file <"
                    << line << "> could not be opened!"
                    << std::endl;

          return false;
        }
      else
        return read_input (input);
    }

  // this line matched nothing known
  std::cerr << "Line <" << current_line_n
            << "> of file <" << input_filename
            << ">: This line matched nothing known ('set' or 'subsection' missing!?):" << std::endl
            << "    " << line << std::endl;
  return false;
}



std::size_t
ParameterHandler::memory_consumption () const
{
//TODO: add to this an estimate of the memory in the property_tree
  return (MemoryConsumption::memory_consumption (subsection_path));
}



bool
ParameterHandler::operator == (const ParameterHandler &prm2)  const
{
  if (patterns.size() != prm2.patterns.size())
    return false;

  for (unsigned int j=0; j<patterns.size(); ++j)
    if (patterns[j]->description() != prm2.patterns[j]->description())
      return false;

  // instead of walking through all
  // the nodes of the two trees
  // entries and prm2.entries and
  // comparing them for equality,
  // simply dump the content of the
  // entire structure into a string
  // and compare those for equality
  std::ostringstream o1, o2;
  write_json (o1, *entries);
  write_json (o2, *prm2.entries);
  return (o1.str() == o2.str());
}




MultipleParameterLoop::UserClass::~UserClass ()
{}



MultipleParameterLoop::MultipleParameterLoop()
  :
  n_branches(0)
{}



MultipleParameterLoop::~MultipleParameterLoop ()
{}



bool MultipleParameterLoop::read_input (std::istream &input,
                                        const std::string &filename)
{
  AssertThrow (input, ExcIO());

  bool x = ParameterHandler::read_input (input, filename);
  if (x)
    init_branches ();
  return x;
}



void MultipleParameterLoop::loop (MultipleParameterLoop::UserClass &uc)
{
  for (unsigned int run_no=0; run_no<n_branches; ++run_no)
    {
      // give create_new one-based numbers
      uc.create_new (run_no+1);
      fill_entry_values (run_no);
      uc.run (*this);
    };
}



void MultipleParameterLoop::init_branches ()
{
  multiple_choices.clear ();
  init_branches_current_section ();

  // split up different values
  for (unsigned int i=0; i<multiple_choices.size(); ++i)
    multiple_choices[i].split_different_values ();

  // finally calculate number of branches
  n_branches = 1;
  for (unsigned int i=0; i<multiple_choices.size(); ++i)
    if (multiple_choices[i].type == Entry::variant)
      n_branches *= multiple_choices[i].different_values.size();

  // check whether array entries have the correct
  // number of entries
  for (unsigned int i=0; i<multiple_choices.size(); ++i)
    if (multiple_choices[i].type == Entry::array)
      if (multiple_choices[i].different_values.size() != n_branches)
        std::cerr << "    The entry value" << std::endl
                  << "        " << multiple_choices[i].entry_value << std::endl
                  << "    for the entry named" << std::endl
                  << "        " << multiple_choices[i].entry_name << std::endl
                  << "    does not have the right number of entries for the " << std::endl
                  << "        " << n_branches << " variant runs that will be performed."
                  << std::endl;


  // do a first run on filling the values to
  // check for the conformance with the regexp
  // (later on, this will be lost in the whole
  // other output)
  for (unsigned int i=0; i<n_branches; ++i)
    fill_entry_values (i);
}



void MultipleParameterLoop::init_branches_current_section ()
{
  const boost::property_tree::ptree &current_section
    = entries->get_child (get_current_path());

  // check all entries in the present
  // subsection whether they are
  // multiple entries
  //
  // we loop over entries in sorted
  // order to guarantee backward
  // compatibility to an earlier
  // implementation
  for (boost::property_tree::ptree::const_assoc_iterator
       p = current_section.ordered_begin();
       p != current_section.not_found(); ++p)
    if (is_parameter_node (p->second) == true)
      {
        const std::string value = p->second.get<std::string>("value");
        if (value.find('{') != std::string::npos)
          multiple_choices.push_back (Entry(subsection_path,
                                            demangle(p->first),
                                            value));
      }

  // then loop over all subsections
  for (boost::property_tree::ptree::const_iterator
       p = current_section.begin();
       p != current_section.end(); ++p)
    if (is_parameter_node (p->second) == false)
      {
        enter_subsection (demangle(p->first));
        init_branches_current_section ();
        leave_subsection ();
      }
}




void MultipleParameterLoop::fill_entry_values (const unsigned int run_no)
{
  unsigned int possibilities = 1;

  std::vector<Entry>::iterator choice;
  for (choice = multiple_choices.begin();
       choice != multiple_choices.end();
       ++choice)
    {
      const unsigned int selection
        = (run_no/possibilities) % choice->different_values.size();
      std::string entry_value;
      if (choice->type == Entry::variant)
        entry_value = choice->different_values[selection];
      else
        {
          if (run_no>=choice->different_values.size())
            {
              std::cerr << "The given array for entry <"
                        << choice->entry_name
                        << "> does not contain enough elements! Taking empty string instead."
                        << std::endl;
              entry_value = "";
            }
          else
            entry_value = choice->different_values[run_no];
        }

      // temporarily enter the
      // subsection tree of this
      // multiple entry, set the
      // value, and get out
      // again. the set() operation
      // also tests for the
      // correctness of the value
      // with regard to the pattern
      subsection_path.swap (choice->subsection_path);
      set (choice->entry_name, entry_value);
      subsection_path.swap (choice->subsection_path);

      // move ahead if it was a variant entry
      if (choice->type == Entry::variant)
        possibilities *= choice->different_values.size();
    }
}




std::size_t
MultipleParameterLoop::memory_consumption () const
{
  std::size_t mem = ParameterHandler::memory_consumption ();
  for (unsigned int i=0; i<multiple_choices.size(); ++i)
    mem += multiple_choices[i].memory_consumption ();

  return mem;
}



MultipleParameterLoop::Entry::Entry (const std::vector<std::string> &ssp,
                                     const std::string              &Name,
                                     const std::string              &Value)
  :
  subsection_path (ssp), entry_name(Name), entry_value(Value), type (Entry::array)
{}



void MultipleParameterLoop::Entry::split_different_values ()
{
  // split string into three parts:
  // part before the opening "{",
  // the selection itself, final
  // part after "}"
  std::string prefix  (entry_value, 0, entry_value.find('{'));
  std::string multiple(entry_value, entry_value.find('{')+1,
                       entry_value.rfind('}')-entry_value.find('{')-1);
  std::string postfix (entry_value, entry_value.rfind('}')+1, std::string::npos);
  // if array entry {{..}}: delete inner
  // pair of braces
  if (multiple[0]=='{')
    multiple.erase (0,1);
  if (multiple[multiple.size()-1] == '}')
    multiple.erase (multiple.size()-1, 1);
  // erase leading and trailing spaces
  // in multiple
  while (std::isspace (multiple[0])) multiple.erase (0,1);
  while (std::isspace (multiple[multiple.size()-1])) multiple.erase (multiple.size()-1,1);

  // delete spaces around '|'
  while (multiple.find(" |") != std::string::npos)
    multiple.replace (multiple.find(" |"), 2, "|");
  while (multiple.find("| ") != std::string::npos)
    multiple.replace (multiple.find("| "), 2, "|");

  while (multiple.find('|') != std::string::npos)
    {
      different_values.push_back (prefix +
                                  std::string(multiple, 0, multiple.find('|'))+
                                  postfix);
      multiple.erase (0, multiple.find('|')+1);
    };
  // make up the last selection ("while" broke
  // because there was no '|' any more
  different_values.push_back (prefix+multiple+postfix);
  // finally check whether this was a variant
  // entry ({...}) or an array ({{...}})
  if ((entry_value.find("{{") != std::string::npos) &&
      (entry_value.find("}}") != std::string::npos))
    type = Entry::array;
  else
    type = Entry::variant;
}


std::size_t
MultipleParameterLoop::Entry::memory_consumption () const
{
  return (MemoryConsumption::memory_consumption (subsection_path) +
          MemoryConsumption::memory_consumption (entry_name) +
          MemoryConsumption::memory_consumption (entry_value) +
          MemoryConsumption::memory_consumption (different_values) +
          sizeof (type));
}

DEAL_II_NAMESPACE_CLOSE
