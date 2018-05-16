// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2018 by the deal.II authors
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

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <boost/io/ios_state.hpp>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <algorithm>
#include <sstream>
#include <cctype>
#include <limits>
#include <cstring>


DEAL_II_NAMESPACE_OPEN


ParameterHandler::ParameterHandler ()
  :
  entries (new boost::property_tree::ptree())
{}


namespace
{
  std::string
  mangle (const std::string &s)
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
  demangle (const std::string &s)
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
ParameterHandler::collate_path_string (const std::vector<std::string> &subsection_path)
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
ParameterHandler::get_current_path () const
{
  return collate_path_string (subsection_path);
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



void
ParameterHandler::parse_input (std::istream &input,
                               const std::string &filename,
                               const std::string &last_line)
{
  AssertThrow (input, ExcIO());

  // store subsections we are currently in
  const std::vector<std::string> saved_path = subsection_path;

  std::string input_line;
  std::string fully_concatenated_line;
  bool is_concatenated = false;
  // Maintain both the current line number and the current logical line
  // number, where the latter refers to the line number where (possibly) the
  // current line continuation started.
  unsigned int current_line_n = 0;
  unsigned int current_logical_line_n = 0;

  // define an action that tries to scan a line.
  //
  // if that fails, i.e., if scan_line throws
  // an exception either because a parameter doesn't match its
  // pattern or because an associated action throws an exception,
  // then try to rewind the set of subsections to the same
  // point where we were when the current function was called.
  // this at least allows to read parameters from a predictable
  // state, rather than leave the subsection stack in some
  // unknown state.
  //
  // after unwinding the subsection stack, just re-throw the exception
  auto scan_line_or_cleanup =
    [this, &saved_path](const std::string &line,
                        const std::string &filename,
                        const unsigned int line_number)
  {
    try
      {
        scan_line (line, filename, line_number);
      }
    catch (...)
      {
        while ((saved_path != subsection_path)
               &&
               (subsection_path.size() > 0))
          leave_subsection ();

        throw;
      }
  };


  while (std::getline (input, input_line))
    {
      ++current_line_n;
      if (!is_concatenated)
        current_logical_line_n = current_line_n;
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
          scan_line_or_cleanup (fully_concatenated_line,
                                filename,
                                current_logical_line_n);

          fully_concatenated_line.clear();
        }
    }

  // While it does not make much sense for anyone to actually do this, allow
  // the last line to end in a backslash. to do do, we need to parse
  // whatever was left in the stash of concatenated lines
  if (is_concatenated)
    scan_line_or_cleanup (fully_concatenated_line,
                          filename,
                          current_line_n);

  if (saved_path != subsection_path)
    {
      std::stringstream paths_message;
      if (saved_path.size() > 0)
        {
          paths_message << "Path before loading input:\n";
          for (unsigned int i=0; i < subsection_path.size(); ++i)
            {
              paths_message << std::setw(i*2 + 4)
                            << " "
                            << "subsection "
                            << saved_path[i]
                            << '\n';
            }
          paths_message << "Current path:\n";
          for (unsigned int i=0; i < subsection_path.size(); ++i)
            {
              paths_message << std::setw(i*2 + 4)
                            << " "
                            << "subsection "
                            << subsection_path[i]
                            << (i == subsection_path.size() - 1 ? "" : "\n");
            }
        }
      // restore subsection we started with before throwing the exception:
      subsection_path = saved_path;
      AssertThrow (false, ExcUnbalancedSubsections (filename, paths_message.str()));
    }
}



void
ParameterHandler::parse_input (const std::string &filename,
                               const std::string &last_line)
{
  PathSearch search("PARAMETERS");

  std::string openname = search.find(filename);
  std::ifstream file_stream (openname.c_str());
  parse_input (file_stream, filename, last_line);
}



void
ParameterHandler::parse_input_from_string (const char *s,
                                           const std::string &last_line)
{
  std::istringstream input_stream (s);
  parse_input (input_stream, "input string", last_line);
}



namespace
{
  // Recursively go through the 'source' tree and see if we can find
  // corresponding entries in the 'destination' tree. If not, error out
  // (i.e. we have just read an XML file that has entries that weren't
  // declared in the ParameterHandler object); if so, copy the value of these
  // nodes into the destination object
  void
  read_xml_recursively (const boost::property_tree::ptree &source,
                        const std::string                 &current_path,
                        const char                         path_separator,
                        const std::vector<std::unique_ptr<const Patterns::PatternBase> > &
                        patterns,
                        boost::property_tree::ptree       &destination)
  {
    for (boost::property_tree::ptree::const_iterator p = source.begin();
         p != source.end(); ++p)
      {
        // a sub-tree must either be a parameter node or a subsection
        if (p->second.get_optional<std::string>("value"))
          {
            // make sure we have a corresponding entry in the destination
            // object as well
            const std::string full_path
              = (current_path == ""
                 ?
                 p->first
                 :
                 current_path + path_separator + p->first);

            const std::string new_value
              = p->second.get<std::string>("value");
            AssertThrow (destination.get_optional<std::string> (full_path)
                         &&
                         destination.get_optional<std::string>
                         (full_path + path_separator + "value"),
                         ParameterHandler::ExcEntryUndeclared (p->first));

            // first make sure that the new entry actually satisfies its
            // constraints
            const unsigned int pattern_index
              = destination.get<unsigned int> (full_path +
                                               path_separator +
                                               "pattern");
            AssertThrow (patterns[pattern_index]->match(new_value),
                         ParameterHandler::ExcInvalidEntryForPatternXML
                         // XML entries sometimes have extra surrounding
                         // newlines
                         (Utilities::trim(new_value),
                          p->first,
                          patterns[pattern_index]->description()));

            // set the found parameter in the destination argument
            destination.put (full_path + path_separator + "value",
                             new_value);

            // this node might have sub-nodes in addition to "value", such as
            // "default_value", "documentation", etc. we might at some point
            // in the future want to make sure that if they exist that they
            // match the ones in the 'destination' tree
          }
        else if (p->second.get_optional<std::string>("alias"))
          {
            // it is an alias node. alias nodes are static and there is
            // nothing to do here (but the same applies as mentioned in the
            // comment above about the static nodes inside parameter nodes
          }
        else
          {
            // it must be a subsection
            read_xml_recursively (p->second,
                                  (current_path == "" ?
                                   p->first :
                                   current_path + path_separator + p->first),
                                  path_separator,
                                  patterns,
                                  destination);
          }
      }
  }
}



void
ParameterHandler::parse_input_from_xml (std::istream &in)
{
  AssertThrow(in, ExcIO());
  // read the XML tree assuming that (as we
  // do in print_parameters(XML) it has only
  // a single top-level node called
  // "ParameterHandler"
  boost::property_tree::ptree single_node_tree;
  // This boost function will raise an exception if this is not a valid XML
  // file.
  read_xml (in, single_node_tree);

  // make sure there is a single top-level element
  // called "ParameterHandler"
  AssertThrow (single_node_tree.get_optional<std::string>("ParameterHandler"),
               ExcInvalidXMLParameterFile ("There is no top-level XML element "
                                           "called \"ParameterHandler\"."));

  const std::size_t n_top_level_elements = std::distance (single_node_tree.begin(),
                                                          single_node_tree.end());
  if (n_top_level_elements != 1)
    {
      std::ostringstream top_level_message;
      top_level_message << "The ParameterHandler input parser found "
                        << n_top_level_elements
                        << " top level elements while reading\n "
                        << "    an XML format input file, but there should be"
                        << " exactly one top level element.\n"
                        << "    The top level elements are:\n";

      unsigned int entry_n = 0;
      for (boost::property_tree::ptree::iterator it = single_node_tree.begin();
           it != single_node_tree.end(); ++it, ++entry_n)
        {
          top_level_message << "        "
                            << it->first
                            << (entry_n != n_top_level_elements - 1 ? "\n" : "");
        }

      // repeat assertion condition to make the printed version easier to read
      AssertThrow (n_top_level_elements == 1,
                   ExcInvalidXMLParameterFile (top_level_message.str()));
    }

  // read the child elements recursively
  const boost::property_tree::ptree
  &my_entries = single_node_tree.get_child("ParameterHandler");

  read_xml_recursively (my_entries, "", path_separator, patterns, *entries);
}


void
ParameterHandler::parse_input_from_json (std::istream &in)
{
  AssertThrow(in, ExcIO());

  boost::property_tree::ptree node_tree;
  // This boost function will raise an exception if this is not a valid JSON
  // file.
  read_json (in, node_tree);

  // The xml function is reused to read in the xml into the parameter file.
  // This means that only mangled files can be read.
  read_xml_recursively (node_tree, "", path_separator, patterns, *entries);
}



void
ParameterHandler::clear ()
{
  entries = std_cxx14::make_unique<boost::property_tree::ptree> ();
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

  patterns.reserve (patterns.size() + 1);
  patterns.emplace_back (pattern.clone());
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
ParameterHandler::add_action(const std::string &entry,
                             const std::function<void (const std::string &)> &action)
{
  actions.push_back (action);

  // get the current list of actions, if any
  boost::optional<std::string>
  current_actions
    =
      entries->get_optional<std::string> (get_current_full_path(entry) + path_separator + "actions");

  // if there were actions already associated with this parameter, add
  // the current one to it; otherwise, create a one-item list and use
  // that
  if (current_actions)
    {
      const std::string all_actions = current_actions.get() +
                                      "," +
                                      Utilities::int_to_string (actions.size()-1);
      entries->put (get_current_full_path(entry) + path_separator + "actions",
                    all_actions);
    }
  else
    entries->put (get_current_full_path(entry) + path_separator + "actions",
                  Utilities::int_to_string(actions.size()-1));


  // as documented, run the action on the default value at the very end
  const std::string default_value = entries->get<std::string>(get_current_full_path(entry) +
                                                              path_separator +
                                                              "default_value");
  action (default_value);
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



void
ParameterHandler::enter_subsection (const std::string &subsection)
{
  // if necessary create subsection
  if (!entries->get_child_optional (get_current_full_path(subsection)))
    entries->add_child (get_current_full_path(subsection),
                        boost::property_tree::ptree());

  // then enter it
  subsection_path.push_back (subsection);
}



void
ParameterHandler::leave_subsection ()
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



long int
ParameterHandler::get_integer (const std::string &entry_string) const
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



double
ParameterHandler::get_double (const std::string &entry_string) const
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



bool
ParameterHandler::get_bool (const std::string &entry_string) const
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

  // get the node for the entry. if it doesn't exist, then we end up
  // in the else-branch below, which asserts that the entry is indeed
  // declared
  if (entries->get_optional<std::string>(path + path_separator + "value"))
    {
      // verify that the new value satisfies the provided pattern
      const unsigned int pattern_index
        = entries->get<unsigned int> (path + path_separator + "pattern");
      AssertThrow (patterns[pattern_index]->match(new_value),
                   ExcValueDoesNotMatchPattern (new_value,
                                                entries->get<std::string>
                                                (path +
                                                 path_separator +
                                                 "pattern_description")));

      // then also execute the actions associated with this
      // parameter (if any have been provided)
      const boost::optional<std::string> action_indices_as_string
        = entries->get_optional<std::string> (path + path_separator + "actions");
      if (action_indices_as_string)
        {
          std::vector<int> action_indices
            = Utilities::string_to_int (Utilities::split_string_list (action_indices_as_string.get()));
          for (unsigned int index : action_indices)
            if (actions.size() >= index+1)
              actions[index](new_value);
        }

      // finally write the new value into the database
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
                                    const OutputStyle style) const
{
  AssertThrow (out, ExcIO());

  // we'll have to print some text that is padded with spaces;
  // set the appropriate fill character, but also make sure that
  // we will restore the previous setting (and all other stream
  // flags) when we exit this function
  boost::io::ios_flags_saver restore_flags(out);
  boost::io::basic_ios_fill_saver<char> restore_fill_state(out);
  out.fill(' ');

  // we treat XML and JSON is one step via BOOST, whereas all of the others are
  // done recursively in our own code. take care of the two special formats
  // first
  if (style == XML)
    {
      // call the writer function and exit as there is nothing
      // further to do down in this function
      //
      // XML has a requirement that there can only be one
      // single top-level entry, but we may have multiple
      // entries and sections.  we work around this by
      // creating a tree just for this purpose with the
      // single top-level node "ParameterHandler" and
      // assign the existing tree under it
      boost::property_tree::ptree single_node_tree;
      single_node_tree.add_child("ParameterHandler",
                                 *entries);

      write_xml (out, single_node_tree);
      return out;
    }

  if (style == JSON)
    {
      write_json (out, *entries);
      return out;
    }

  // for all of the other formats, print a preamble:
  switch (style)
    {
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
    }

  // dive recursively into the subsections
  recursively_print_parameters (std::vector<std::string>(), // start at the top level
                                style,
                                0,
                                out);

  return out;
}



void
ParameterHandler::recursively_print_parameters (const std::vector<std::string> &target_subsection_path,
                                                const OutputStyle               style,
                                                const unsigned int              indent_level,
                                                std::ostream                   &out) const
{
  AssertThrow (out, ExcIO());

  // this function should not be necessary for XML or JSON output...
  Assert ((style != XML) && (style != JSON),
          ExcInternalError());

  const boost::property_tree::ptree &current_section
    = entries->get_child (collate_path_string(target_subsection_path));

  unsigned int overall_indent_level = indent_level;

  switch (style)
    {
    case Text:
    case ShortText:
    {
      // first find out the longest entry name to be able to align the
      // equal signs to do this loop over all nodes of the current
      // tree, select the parameter nodes (and discard sub-tree nodes)
      // and take the maximum of their lengths
      //
      // likewise find the longest actual value string to make sure we
      // can align the default and documentation strings
      std::size_t longest_name  = 0;
      std::size_t longest_value = 0;
      for (boost::property_tree::ptree::const_iterator
           p = current_section.begin();
           p != current_section.end(); ++p)
        if (is_parameter_node (p->second) == true)
          {
            longest_name = std::max (longest_name,
                                     demangle(p->first).length());
            longest_value = std::max (longest_value,
                                      p->second.get<std::string>("value").length());
          }

      // print entries one by one. make sure they are sorted by using
      // the appropriate iterators
      bool first_entry = true;
      for (boost::property_tree::ptree::const_assoc_iterator
           p = current_section.ordered_begin();
           p != current_section.not_found(); ++p)
        if (is_parameter_node (p->second) == true)
          {
            const std::string value = p->second.get<std::string>("value");

            // if there is documentation, then add an empty line
            // (unless this is the first entry in a subsection), print
            // the documentation, and then the actual entry; break the
            // documentation into readable chunks such that the whole
            // thing is at most 78 characters wide
            if ((style == Text) &&
                !p->second.get<std::string>("documentation").empty())
              {
                if (first_entry == false)
                  out << '\n';
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
                      << '\n';
              }



            // print name and value of this entry
            out << std::setw(overall_indent_level*2) << ""
                << "set "
                << demangle(p->first)
                << std::setw(longest_name-demangle(p->first).length()+1) << " "
                << "= " << value;

            // finally print the default value, but only if it differs
            // from the actual value
            if ((style == Text) && value != p->second.get<std::string>("default_value"))
              {
                out << std::setw(longest_value-value.length()+1) << ' '
                    << "# ";
                out << "default: " << p->second.get<std::string>("default_value");
              }

            out << '\n';
          }

      break;
    }

    case LaTeX:
    {
      auto escape = [] (const std::string &input)
      {
        return Patterns::internal::escape(input, Patterns::PatternBase::LaTeX);
      };

      // if there are any parameters in this section then print them
      // as an itemized list
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
              << '\n';

          // print entries one by one. make sure they are sorted by
          // using the appropriate iterators
          for (boost::property_tree::ptree::const_assoc_iterator
               p = current_section.ordered_begin();
               p != current_section.not_found(); ++p)
            if (is_parameter_node (p->second) == true)
              {
                const std::string value = p->second.get<std::string>("value");

                // print name
                out << "\\item {\\it Parameter name:} {\\tt "
                    << escape(demangle(p->first)) << "}\n"
                    << "\\phantomsection";
                {
                  // create label: labels are not to be escaped but mangled
                  std::string label = "parameters:";
                  for (unsigned int i=0; i<target_subsection_path.size(); ++i)
                    {
                      label.append(mangle(target_subsection_path[i]));
                      label.append("/");
                    }
                  label.append(p->first);
                  // Backwards-compatibility. Output the label with and without
                  // escaping whitespace:
                  if (label.find("_20")!=std::string::npos)
                    out << "\\label{" << Utilities::replace_in_string(label,"_20", " ") << "}\n";
                  out << "\\label{" << label << "}\n";
                }
                out << "\n\n";

                out << "\\index[prmindex]{"
                    << escape(demangle(p->first))
                    << "}\n";
                out << "\\index[prmindexfull]{";
                for (unsigned int i=0; i<target_subsection_path.size(); ++i)
                  out << escape(target_subsection_path[i]) << "!";
                out << escape(demangle(p->first))
                    << "}\n";

                // finally print value and default
                out << "{\\it Value:} " << escape(value) << "\n\n"
                    << '\n'
                    << "{\\it Default:} "
                    << escape(p->second.get<std::string>("default_value"))
                    << "\n\n"
                    << '\n';

                // if there is a documenting string, print it as well but
                // don't escape to allow formatting/formulas
                if (!p->second.get<std::string>("documentation").empty())
                  out << "{\\it Description:} "
                      << p->second.get<std::string>("documentation") << "\n\n"
                      << '\n';

                // also output possible values, do not escape because the
                // description internally will use LaTeX formatting
                const unsigned int pattern_index = p->second.get<unsigned int> ("pattern");
                const std::string desc_str = patterns[pattern_index]->description (Patterns::PatternBase::LaTeX);
                out << "{\\it Possible values:} "
                    << desc_str
                    << '\n';
              }
            else if (is_alias_node (p->second) == true)
              {
                const std::string alias = p->second.get<std::string>("alias");

                // print name
                out << "\\item {\\it Parameter name:} {\\tt "
                    << escape(demangle(p->first)) << "}\n"
                    << "\\phantomsection";
                {
                  // create label: labels are not to be escaped but mangled
                  std::string label = "parameters:";
                  for (unsigned int i=0; i<target_subsection_path.size(); ++i)
                    {
                      label.append(mangle(target_subsection_path[i]));
                      label.append("/");
                    }
                  label.append(p->first);
                  // Backwards-compatibility. Output the label with and without
                  // escaping whitespace:
                  if (label.find("_20")!=std::string::npos)
                    out << "\\label{" << Utilities::replace_in_string(label,"_20", " ") << "}\n";
                  out << "\\label{" << label << "}\n";
                }
                out << "\n\n";

                out << "\\index[prmindex]{"
                    << escape(demangle(p->first))
                    << "}\n";
                out << "\\index[prmindexfull]{";
                for (unsigned int i=0; i<target_subsection_path.size(); ++i)
                  out << escape(target_subsection_path[i]) << "!";
                out << escape(demangle(p->first))
                    << "}\n";

                // finally print alias and indicate if it is deprecated
                out << "This parameter is an alias for the parameter ``\\texttt{"
                    << escape(alias) << "}''."
                    << (p->second.get<std::string>("deprecation_status") == "true"
                        ?
                        " Its use is deprecated."
                        :
                        "")
                    << "\n\n"
                    << '\n';
              }
          out << "\\end{itemize}"
              << '\n';
        }

      break;
    }

    case Description:
    {
      // first find out the longest entry name to be able to align the
      // equal signs
      std::size_t longest_name = 0;
      for (boost::property_tree::ptree::const_iterator
           p = current_section.begin();
           p != current_section.end(); ++p)
        if (is_parameter_node (p->second) == true)
          longest_name = std::max (longest_name,
                                   demangle(p->first).length());

      // print entries one by one. make sure they are sorted by using
      // the appropriate iterators
      for (boost::property_tree::ptree::const_assoc_iterator
           p = current_section.ordered_begin();
           p != current_section.not_found(); ++p)
        if (is_parameter_node (p->second) == true)
          {
            // print name and value
            out << std::setw(overall_indent_level*2) << ""
                << "set "
                << demangle(p->first)
                << std::setw(longest_name-demangle(p->first).length()+1) << " "
                << " = ";

            // print possible values:
            const unsigned int pattern_index = p->second.get<unsigned int> ("pattern");
            const std::string full_desc_str = patterns[pattern_index]->description (Patterns::PatternBase::Text);
            const std::vector<std::string> description_str
              = Utilities::break_text_into_lines (full_desc_str,
                                                  78 - overall_indent_level*2 - 2, '|');
            if (description_str.size() > 1)
              {
                out << '\n';
                for (unsigned int i=0; i<description_str.size(); ++i)
                  out << std::setw(overall_indent_level*2+6) << ""
                      << description_str[i] << '\n';
              }
            else if (description_str.empty() == false)
              out << "  " << description_str[0] << '\n';
            else
              out << '\n';

            // if there is a documenting string, print it as well
            if (p->second.get<std::string>("documentation").length() != 0)
              out << std::setw(overall_indent_level*2 + longest_name + 10) << ""
                  << "(" << p->second.get<std::string>("documentation") << ")"
                  << '\n';
          }

      break;
    }

    default:
      Assert (false, ExcNotImplemented());
    }


  // if there was text before and there are sections to come, put two
  // newlines between the last entry and the first subsection
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
        (style != ShortText)
        &&
        (n_parameters != 0)
        &&
        (n_sections != 0))
      out << "\n\n";
  }

  // now traverse subsections tree, in alphabetical order
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
                << "subsection " << demangle(p->first)
                << '\n';
            break;

          case LaTeX:
          {
            auto escape = [] (const std::string &input)
            {
              return Patterns::internal::escape(input, Patterns::PatternBase::LaTeX);
            };

            out << '\n'
                << "\\subsection{Parameters in section \\tt ";

            // find the path to the current section so that we can
            // print it in the \subsection{...} heading
            for (unsigned int i=0; i<target_subsection_path.size(); ++i)
              out << escape(target_subsection_path[i]) << "/";
            out << escape(demangle(p->first));

            out << "}" << '\n';
            out << "\\label{parameters:";
            for (unsigned int i=0; i<target_subsection_path.size(); ++i)
              out << mangle(target_subsection_path[i]) << "/";
            out << p->first << "}";
            out << '\n';

            out << '\n';
            break;
          }

          default:
            Assert (false, ExcNotImplemented());
          }

        // then the contents of the subsection
        const std::string subsection = demangle(p->first);
        std::vector<std::string> directory_path = target_subsection_path;
        directory_path.emplace_back (subsection);

        recursively_print_parameters (directory_path, style,
                                      overall_indent_level+1,
                                      out);

        switch (style)
          {
          case Text:
            // write end of subsection. one blank line after each
            // subsection
            out << std::setw(overall_indent_level*2) << ""
                << "end" << '\n'
                << '\n';

            // if this is a toplevel subsection, then have two
            // newlines
            if (overall_indent_level == 0)
              out << '\n';

            break;

          case Description:
            break;

          case ShortText:
            // write end of subsection.
            out << std::setw(overall_indent_level*2) << ""
                << "end" << '\n';
            break;

          case LaTeX:
            break;

          default:
            Assert (false, ExcNotImplemented());
          }
      }
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
          // call the writer function and exit as there is nothing
          // further to do down in this function
          //
          // XML has a requirement that there can only be one single
          // top-level entry, but a section has multiple entries and
          // sections. we work around this by creating a tree just for
          // this purpose with the single top-level node
          // "ParameterHandler" and assign the full path of down to
          // the current section under it
          boost::property_tree::ptree single_node_tree;

          // if there is no subsection selected, add the whole tree of entries,
          // otherwise add a root element and the selected subsection under it
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
            }

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
          }

      // first find out the longest entry name to be able to align the
      // equal signs to do this loop over all nodes of the current
      // tree, select the parameter nodes (and discard sub-tree nodes)
      // and take the maximum of their lengths
      std::size_t longest_name = 0;
      for (boost::property_tree::ptree::const_iterator
           p = current_section.begin();
           p != current_section.end(); ++p)
        if (is_parameter_node (p->second) == true)
          longest_name = std::max (longest_name,
                                   demangle(p->first).length());

      // likewise find the longest actual value string to make sure we
      // can align the default and documentation strings
      std::size_t longest_value = 0;
      for (boost::property_tree::ptree::const_iterator
           p = current_section.begin();
           p != current_section.end(); ++p)
        if (is_parameter_node (p->second) == true)
          longest_value = std::max (longest_value,
                                    p->second.get<std::string>("value").length());


      // print entries one by one. make sure they are sorted by using
      // the appropriate iterators
      bool first_entry = true;
      for (boost::property_tree::ptree::const_assoc_iterator
           p = current_section.ordered_begin();
           p != current_section.not_found(); ++p)
        if (is_parameter_node (p->second) == true)
          {
            const std::string value = p->second.get<std::string>("value");

            // if there is documentation, then add an empty line
            // (unless this is the first entry in a subsection), print
            // the documentation, and then the actual entry; break the
            // documentation into readable chunks such that the whole
            // thing is at most 78 characters wide
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



            // print name and value of this entry
            out << std::setw(overall_indent_level*2) << ""
                << "set "
                << demangle(p->first)
                << std::setw(longest_name-demangle(p->first).length()+1) << " "
                << "= " << value;

            // finally print the default value, but only if it differs
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
      auto escape = [] (const std::string &input)
      {
        return Patterns::internal::escape(input, Patterns::PatternBase::LaTeX);
      };

      // if there are any parameters in this section then print them
      // as an itemized list
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

          // print entries one by one. make sure they are sorted by
          // using the appropriate iterators
          for (boost::property_tree::ptree::const_assoc_iterator
               p = current_section.ordered_begin();
               p != current_section.not_found(); ++p)
            if (is_parameter_node (p->second) == true)
              {
                const std::string value = p->second.get<std::string>("value");

                // print name
                out << "\\item {\\it Parameter name:} {\\tt " << escape(demangle(p->first)) << "}\n"
                    << "\\phantomsection\\label{parameters:";
                for (unsigned int i=0; i<subsection_path.size(); ++i)
                  out << mangle(subsection_path[i]) << "/";
                out << p->first;
                out << "}\n\n"
                    << std::endl;

                out << "\\index[prmindex]{"
                    << escape(demangle(p->first))
                    << "}\n";
                out << "\\index[prmindexfull]{";
                for (unsigned int i=0; i<subsection_path.size(); ++i)
                  out << escape(subsection_path[i]) << "!";
                out << escape(demangle(p->first))
                    << "}\n";

                // finally print value and default
                out << "{\\it Value:} " << escape(value) << "\n\n"
                    << std::endl
                    << "{\\it Default:} "
                    << escape(p->second.get<std::string>("default_value")) << "\n\n"
                    << std::endl;

                // if there is a documenting string, print it as well
                if (!p->second.get<std::string>("documentation").empty())
                  out << "{\\it Description:} "
                      << p->second.get<std::string>("documentation") << "\n\n"
                      << std::endl;

                // also output possible values
                const unsigned int pattern_index = p->second.get<unsigned int> ("pattern");
                const std::string desc_str = patterns[pattern_index]->description (Patterns::PatternBase::LaTeX);
                out << "{\\it Possible values:} "
                    << desc_str
                    << std::endl;
              }
            else if (is_alias_node (p->second) == true)
              {
                const std::string alias = p->second.get<std::string>("alias");

                // print name
                out << "\\item {\\it Parameter name:} {\\tt " << escape(demangle(p->first)) << "}\n"
                    << "\\phantomsection\\label{parameters:";
                for (unsigned int i=0; i<subsection_path.size(); ++i)
                  out << mangle(subsection_path[i]) << "/";
                out << p->first;
                out << "}\n\n"
                    << std::endl;

                out << "\\index[prmindex]{"
                    << escape(demangle(p->first))
                    << "}\n";
                out << "\\index[prmindexfull]{";
                for (unsigned int i=0; i<subsection_path.size(); ++i)
                  out << escape(subsection_path[i]) << "!";
                out << escape(demangle(p->first))
                    << "}\n";

                // finally print alias and indicate if it is deprecated
                out << "This parameter is an alias for the parameter ``\\texttt{"
                    << escape(alias) << "}''."
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

      // first find out the longest entry name to be able to align the
      // equal signs
      std::size_t longest_name = 0;
      for (boost::property_tree::ptree::const_iterator
           p = current_section.begin();
           p != current_section.end(); ++p)
        if (is_parameter_node (p->second) == true)
          longest_name = std::max (longest_name,
                                   demangle(p->first).length());

      // print entries one by one. make sure they are sorted by using
      // the appropriate iterators
      for (boost::property_tree::ptree::const_assoc_iterator
           p = current_section.ordered_begin();
           p != current_section.not_found(); ++p)
        if (is_parameter_node (p->second) == true)
          {
            // print name and value
            out << std::setw(overall_indent_level*2) << ""
                << "set "
                << demangle(p->first)
                << std::setw(longest_name-demangle(p->first).length()+1) << " "
                << " = ";

            // print possible values:
            const unsigned int pattern_index = p->second.get<unsigned int> ("pattern");
            const std::string full_desc_str = patterns[pattern_index]->description (Patterns::PatternBase::Text);
            const std::vector<std::string> description_str
              = Utilities::break_text_into_lines (full_desc_str,
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

            // if there is a documenting string, print it as well
            if (p->second.get<std::string>("documentation").length() != 0)
              out << std::setw(overall_indent_level*2 + longest_name + 10) << ""
                  << "(" << p->second.get<std::string>("documentation") << ")" << std::endl;
          }

      break;
    }

    default:
      Assert (false, ExcNotImplemented());
    }


  // if there was text before and there are sections to come, put two
  // newlines between the last entry and the first subsection
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

      // now traverse subsections tree, in alphabetical order
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
                auto escape = [] (const std::string &input)
                {
                  return Patterns::internal::escape(input, Patterns::PatternBase::LaTeX);
                };

                out << std::endl
                    << "\\subsection{Parameters in section \\tt ";

                // find the path to the current section so that we can
                // print it in the \subsection{...} heading
                for (unsigned int i=0; i<subsection_path.size(); ++i)
                  out << escape(subsection_path[i]) << "/";
                out << escape(demangle(p->first));

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

            // then the contents of the subsection
            enter_subsection (demangle(p->first));
            print_parameters_section (out, style, overall_indent_level+1);
            leave_subsection ();
            switch (style)
              {
              case Text:
                // write end of subsection. one blank line after each
                // subsection
                out << std::setw(overall_indent_level*2) << ""
                    << "end" << std::endl
                    << std::endl;

                // if this is a toplevel subsection, then have two
                // newlines
                if (overall_indent_level == 0)
                  out << std::endl;

                break;
              case Description:
                break;
              case ShortText:
                // write end of subsection.
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
  // dive recursively into the subsections
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



void
ParameterHandler::scan_line (std::string         line,
                             const std::string  &input_filename,
                             const unsigned int  current_line_n)
{
  // save a copy for some error messages
  const std::string original_line = line;

  // if there is a comment, delete it
  if (line.find('#') != std::string::npos)
    line.erase (line.find('#'), std::string::npos);

  // replace \t by space:
  while (line.find('\t') != std::string::npos)
    line.replace (line.find('\t'), 1, " ");

  //trim start and end:
  line = Utilities::trim(line);

  // if line is now empty: leave
  if (line.length() == 0)
    {
      return;
    }
  // enter subsection
  else if (Utilities::match_at_string_start(line, "SUBSECTION ") ||
           Utilities::match_at_string_start(line, "subsection "))
    {
      // delete this prefix
      line.erase (0, std::string("subsection").length()+1);

      const std::string subsection = Utilities::trim(line);

      // check whether subsection exists
      AssertThrow (entries->get_child_optional (get_current_full_path(subsection)),
                   ExcNoSubsection (current_line_n,
                                    input_filename,
                                    demangle(get_current_full_path(subsection))));

      // subsection exists
      subsection_path.push_back (subsection);
    }
  // exit subsection
  else if (Utilities::match_at_string_start(line, "END") ||
           Utilities::match_at_string_start(line, "end"))
    {
      line.erase (0, 3);
      while ((line.size() > 0) && (std::isspace(line[0])))
        line.erase (0, 1);

      AssertThrow (line.size() == 0,
                   ExcCannotParseLine (current_line_n,
                                       input_filename,
                                       "Invalid content after 'end' or 'END' statement."));
      AssertThrow (subsection_path.size() != 0,
                   ExcCannotParseLine (current_line_n,
                                       input_filename,
                                       "There is no subsection to leave here."));
      leave_subsection ();
    }
  // regular entry
  else if (Utilities::match_at_string_start(line, "SET ") ||
           Utilities::match_at_string_start(line, "set "))
    {
      // erase "set" statement
      line.erase (0, 4);

      std::string::size_type pos = line.find('=');
      AssertThrow (pos != std::string::npos,
                   ExcCannotParseLine (current_line_n,
                                       input_filename,
                                       "Invalid format of 'set' or 'SET' statement."));

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

      // get the node for the entry. if it doesn't exist, then we end up
      // in the else-branch below, which asserts that the entry is indeed
      // declared
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
              // verify that the new value satisfies the provided pattern
              const unsigned int pattern_index
                = entries->get<unsigned int> (path + path_separator + "pattern");
              AssertThrow (patterns[pattern_index]->match (entry_value),
                           ExcInvalidEntryForPattern (current_line_n,
                                                      input_filename,
                                                      entry_value,
                                                      entry_name,
                                                      patterns[pattern_index]->description()));

              // then also execute the actions associated with this
              // parameter (if any have been provided)
              const boost::optional<std::string> action_indices_as_string
                = entries->get_optional<std::string> (path + path_separator + "actions");
              if (action_indices_as_string)
                {
                  std::vector<int> action_indices
                    = Utilities::string_to_int (Utilities::split_string_list (action_indices_as_string.get()));
                  for (unsigned int index : action_indices)
                    if (actions.size() >= index+1)
                      actions[index](entry_value);
                }
            }

          // finally write the new value into the database
          entries->put (path + path_separator + "value",
                        entry_value);
        }
      else
        {
          AssertThrow (false,
                       ExcCannotParseLine(current_line_n,
                                          input_filename,
                                          ("No entry with name <" + entry_name +
                                           "> was declared in the current subsection.")));
        }
    }
  // an include statement?
  else if (Utilities::match_at_string_start(line, "include ") ||
           Utilities::match_at_string_start(line, "INCLUDE "))
    {
      // erase "include " statement and eliminate spaces
      line.erase (0, 7);
      while ((line.size() > 0) && (line[0] == ' '))
        line.erase (0, 1);

      // the remainder must then be a filename
      AssertThrow (line.size() !=0,
                   ExcCannotParseLine (current_line_n,
                                       input_filename,
                                       "The current line is an 'include' or "
                                       "'INCLUDE' statement, but it does not "
                                       "name a file for inclusion."));

      std::ifstream input (line.c_str());
      AssertThrow (input, ExcCannotOpenIncludeStatementFile (current_line_n,
                                                             input_filename,
                                                             line));
      parse_input (input);
    }
  else
    {
      AssertThrow (false,
                   ExcCannotParseLine (current_line_n,
                                       input_filename,
                                       "The line\n\n"
                                       "        <"
                                       + original_line
                                       + ">\n\n"
                                       "could not be parsed: please check to "
                                       "make sure that the file is not missing a "
                                       "'set', 'include', 'subsection', or 'end' "
                                       "statement."));
    }
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



MultipleParameterLoop::MultipleParameterLoop()
  :
  n_branches(0)
{}



void
MultipleParameterLoop::parse_input (std::istream &input,
                                    const std::string &filename,
                                    const std::string &last_line)
{
  AssertThrow (input, ExcIO());

  // Note that (to avoid infinite recursion) we have to explicitly call the
  // base class version of parse_input and *not* a wrapper (which may be
  // virtual and lead us back here)
  ParameterHandler::parse_input (input, filename, last_line);
  init_branches ();
}



void
MultipleParameterLoop::loop (MultipleParameterLoop::UserClass &uc)
{
  for (unsigned int run_no=0; run_no<n_branches; ++run_no)
    {
      // give create_new one-based numbers
      uc.create_new (run_no+1);
      fill_entry_values (run_no);
      uc.run (*this);
    };
}



void
MultipleParameterLoop::init_branches ()
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



void
MultipleParameterLoop::init_branches_current_section ()
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
          multiple_choices.emplace_back (subsection_path,
                                         demangle(p->first),
                                         value);
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




void
MultipleParameterLoop::fill_entry_values (const unsigned int run_no)
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



void
MultipleParameterLoop::Entry::split_different_values ()
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
