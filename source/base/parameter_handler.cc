// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1998 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/logstream.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/path_search.h>
#include <deal.II/base/utilities.h>

#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include <boost/algorithm/string.hpp>
#include <boost/io/ios_state.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#undef BOOST_BIND_GLOBAL_PLACEHOLDERS

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>


DEAL_II_NAMESPACE_OPEN


ParameterHandler::ParameterHandler()
  : entries(new boost::property_tree::ptree())
{}


namespace
{
  std::string
  mangle(const std::string &s)
  {
    std::string u;

    // reserve the minimum number of characters we will need. it may
    // be more but this is the least we can do
    u.reserve(s.size());

    // see if the name is special and if so mangle the whole thing
    const bool mangle_whole_string = (s == "value");

    // for all parts of the string, see if it is an allowed character or not
    for (const char c : s)
      {
        static const std::string allowed_characters(
          "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789");

        if ((!mangle_whole_string) &&
            (allowed_characters.find(c) != std::string::npos))
          u.push_back(c);
        else
          {
            u.push_back('_');
            static const char hex[16] = {'0',
                                         '1',
                                         '2',
                                         '3',
                                         '4',
                                         '5',
                                         '6',
                                         '7',
                                         '8',
                                         '9',
                                         'a',
                                         'b',
                                         'c',
                                         'd',
                                         'e',
                                         'f'};
            u.push_back(hex[static_cast<unsigned char>(c) / 16]);
            u.push_back(hex[static_cast<unsigned char>(c) % 16]);
          }
      }

    return u;
  }



  std::string
  demangle(const std::string &s)
  {
    std::string u;
    u.reserve(s.size());

    for (unsigned int i = 0; i < s.size(); ++i)
      if (s[i] != '_')
        u.push_back(s[i]);
      else
        {
          Assert(i + 2 < s.size(),
                 ExcMessage("Trying to demangle an invalid string."));

          unsigned char c = 0;
          switch (s[i + 1])
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
                DEAL_II_ASSERT_UNREACHABLE();
            }
          switch (s[i + 2])
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
                DEAL_II_ASSERT_UNREACHABLE();
            }

          u.push_back(static_cast<char>(c));

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
  is_parameter_node(const boost::property_tree::ptree &p)
  {
    return static_cast<bool>(p.get_optional<std::string>("value"));
  }


  /**
   * Return whether a given node is a alias node (as opposed
   * to being a subsection or parameter node)
   */
  bool
  is_alias_node(const boost::property_tree::ptree &p)
  {
    return static_cast<bool>(p.get_optional<std::string>("alias"));
  }

  /**
   * Given a list of directories and subdirectories that identify
   * a particular place in the tree, return the string that identifies
   * this place in the way the BOOST property tree libraries likes
   * to identify things.
   */
  std::string
  collate_path_string(const char                      separator,
                      const std::vector<std::string> &subsection_path)
  {
    if (subsection_path.size() > 0)
      {
        std::string p = mangle(subsection_path[0]);
        for (unsigned int i = 1; i < subsection_path.size(); ++i)
          {
            p += separator;
            p += mangle(subsection_path[i]);
          }
        return p;
      }
    else
      return "";
  }

  /**
   * Sort all parameters of the subsection given by the
   * @p target_subsection_path argument in alphabetical order,
   * as well as all subsections within it recursively.
   */
  void
  recursively_sort_parameters(
    const char                      separator,
    const std::vector<std::string> &target_subsection_path,
    boost::property_tree::ptree    &tree)
  {
    boost::property_tree::ptree &current_section =
      tree.get_child(collate_path_string(separator, target_subsection_path));

    // Custom comparator to ensure that the order of sorting is:
    // - sorted parameters and aliases;
    // - sorted subsections.
    static auto compare =
      [](const std::pair<std::string, boost::property_tree::ptree> &a,
         const std::pair<std::string, boost::property_tree::ptree> &b) {
        bool a_is_param =
          (is_parameter_node(a.second) || is_alias_node(a.second));

        bool b_is_param =
          (is_parameter_node(b.second) || is_alias_node(b.second));

        // If a is a parameter/alias and b is a subsection,
        // a should go first, and vice-versa.
        if (a_is_param && !b_is_param)
          return true;

        if (!a_is_param && b_is_param)
          return false;

        // Otherwise, compare a and b.
        return a.first < b.first;
      };

    current_section.sort(compare);

    // Now transverse subsections tree recursively.
    for (const auto &p : current_section)
      {
        if ((is_parameter_node(p.second) == false) &&
            (is_alias_node(p.second) == false))
          {
            const std::string subsection = demangle(p.first);

            std::vector<std::string> subsection_path = target_subsection_path;
            subsection_path.emplace_back(subsection);

            recursively_sort_parameters(separator, subsection_path, tree);
          }
      }
  }

  /**
   * Mangle (@p do_mangle = true) or demangle (@p do_mangle = false) all
   * parameters recursively and attach them to @p tree_out.
   * @p is_parameter_or_alias_node indicates, whether a given node
   * @p tree_in is a parameter node or an alias node (as opposed to being
   * a subsection).
   */
  void
  recursively_mangle_or_demangle(const boost::property_tree::ptree &tree_in,
                                 boost::property_tree::ptree       &tree_out,
                                 const bool                         do_mangle,
                                 const bool is_parameter_or_alias_node = false)
  {
    for (const auto &p : tree_in)
      {
        if (is_parameter_or_alias_node)
          {
            tree_out.put_child(p.first, p.second);
          }
        else
          {
            boost::property_tree::ptree temp;

            if (const auto val = p.second.get_value_optional<std::string>())
              temp.put_value<std::string>(*val);

            recursively_mangle_or_demangle(p.second,
                                           temp,
                                           do_mangle,
                                           is_parameter_node(p.second) ||
                                             is_alias_node(p.second));
            tree_out.put_child(do_mangle ? mangle(p.first) : demangle(p.first),
                               temp);
          }
      }
  }

  /**
   * Assert validity of of given output @p style.
   */
  void
  assert_validity_of_output_style(const ParameterHandler::OutputStyle style)
  {
    AssertThrow(
      (((style & ParameterHandler::XML) != 0) +
       ((style & ParameterHandler::JSON) != 0) +
       ((style & ParameterHandler::PRM) != 0) +
       ((style & ParameterHandler::Description) != 0) +
       ((style & ParameterHandler::LaTeX) != 0)) == 1,
      ExcMessage(
        "You have chosen either no or multiple style formats. You can choose "
        "between: PRM, Description, LaTeX, XML, JSON."));
  }

} // namespace



std::string
ParameterHandler::get_current_path() const
{
  return collate_path_string(path_separator, subsection_path);
}



std::string
ParameterHandler::get_current_full_path(const std::string &name) const
{
  std::string path = get_current_path();
  if (path.empty() == false)
    path += path_separator;

  path += mangle(name);

  return path;
}



std::string
ParameterHandler::get_current_full_path(
  const std::vector<std::string> &sub_path,
  const std::string              &name) const
{
  std::string path = get_current_path();
  if (path.empty() == false)
    path += path_separator;

  if (sub_path.empty() == false)
    path += collate_path_string(path_separator, sub_path) + path_separator;

  path += mangle(name);

  return path;
}



void
ParameterHandler::parse_input(std::istream      &input,
                              const std::string &filename,
                              const std::string &last_line,
                              const bool         skip_undefined)
{
  AssertThrow(input.fail() == false, ExcIO());

  // store subsections we are currently in
  const std::vector<std::string> saved_path = subsection_path;

  std::string input_line;
  std::string fully_concatenated_line;
  bool        is_concatenated = false;
  // Maintain both the current line number and the current logical line
  // number, where the latter refers to the line number where (possibly) the
  // current line continuation started.
  unsigned int current_line_n         = 0;
  unsigned int current_logical_line_n = 0;

  // Keep a list of deprecation messages if we encounter deprecated entries.
  std::string deprecation_messages;

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
    [this,
     &skip_undefined,
     &saved_path,
     &deprecation_messages](const std::string &line,
                            const std::string &filename,
                            const unsigned int line_number) {
      try
        {
          scan_line(line, filename, line_number, skip_undefined);
        }
      catch (ExcEntryIsDeprecated &e)
        {
          deprecation_messages += e.what();
        }
      catch (...)
        {
          while ((saved_path != subsection_path) &&
                 (subsection_path.size() > 0))
            leave_subsection();

          throw;
        }
    };


  while (std::getline(input, input_line))
    {
      ++current_line_n;
      if (!is_concatenated)
        current_logical_line_n = current_line_n;
      // Trim the whitespace at the ends of the line here instead of in
      // scan_line. This makes the continuation line logic a lot simpler.
      input_line = Utilities::trim(input_line);

      // If we see the line which is the same as @p last_line ,
      // terminate the parsing.
      if (last_line.size() != 0 && input_line == last_line)
        break;

      // Check whether or not the current line should be joined with the next
      // line before calling scan_line.
      if (input_line.size() != 0 &&
          input_line.find_last_of('\\') == input_line.size() - 1)
        {
          input_line.erase(input_line.size() - 1); // remove the last '\'
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
          scan_line_or_cleanup(fully_concatenated_line,
                               filename,
                               current_logical_line_n);

          fully_concatenated_line.clear();
        }
    }

  // While it does not make much sense for anyone to actually do this, allow
  // the last line to end in a backslash. To do so, we need to parse
  // whatever was left in the stash of concatenated lines
  if (is_concatenated)
    scan_line_or_cleanup(fully_concatenated_line, filename, current_line_n);

  if (saved_path != subsection_path)
    {
      std::stringstream paths_message;
      if (saved_path.size() > 0)
        {
          paths_message << "Path before loading input:\n";
          for (unsigned int i = 0; i < subsection_path.size(); ++i)
            {
              paths_message << std::setw(i * 2 + 4) << " "
                            << "subsection " << saved_path[i] << '\n';
            }
          paths_message << "Current path:\n";
          for (unsigned int i = 0; i < subsection_path.size(); ++i)
            {
              paths_message << std::setw(i * 2 + 4) << " "
                            << "subsection " << subsection_path[i]
                            << (i == subsection_path.size() - 1 ? "" : "\n");
            }
        }
      // restore subsection we started with before throwing the exception:
      subsection_path = saved_path;
      AssertThrow(false,
                  ExcUnbalancedSubsections(filename, paths_message.str()));
    }

  // if we encountered deprecated entries, throw an exception
  // that contains all the deprecation messages
  AssertThrow(deprecation_messages.empty(),
              ExcEncounteredDeprecatedEntries(deprecation_messages));
}



void
ParameterHandler::parse_input(const std::string &filename,
                              const std::string &last_line,
                              const bool         skip_undefined,
                              const bool assert_mandatory_entries_are_found)
{
  std::ifstream is(filename);
  AssertThrow(is, ExcFileNotOpen(filename));

  std::string file_ending = filename.substr(filename.find_last_of('.') + 1);
  boost::algorithm::to_lower(file_ending);
  if (file_ending == "prm")
    parse_input(is, filename, last_line, skip_undefined);
  else if (file_ending == "xml")
    parse_input_from_xml(is, skip_undefined);
  else if (file_ending == "json")
    parse_input_from_json(is, skip_undefined);
  else
    AssertThrow(false,
                ExcMessage("The given input file <" + filename +
                           "> has a file name extension <" + file_ending +
                           "> that is not recognized. Supported types "
                           "are .prm, .xml, and .json."));

  if (assert_mandatory_entries_are_found)
    assert_that_entries_have_been_set();
}



void
ParameterHandler::parse_input_from_string(const std::string &s,
                                          const std::string &last_line,
                                          const bool         skip_undefined)
{
  std::istringstream input_stream(s);
  parse_input(input_stream, "input string", last_line, skip_undefined);
}



namespace
{
  // Recursively go through the 'source' tree and see if we can find
  // corresponding entries in the 'destination' tree. If not, error out
  // (i.e. we have just read an XML file that has entries that weren't
  // declared in the ParameterHandler object); if so, copy the value of these
  // nodes into the destination object
  void
  read_xml_recursively(
    const boost::property_tree::ptree &source,
    const std::string                 &current_path,
    const char                         path_separator,
    const std::vector<std::unique_ptr<const Patterns::PatternBase>> &patterns,
    const bool        skip_undefined,
    ParameterHandler &prm)
  {
    for (const auto &p : source)
      {
        // a sub-tree must either be a parameter node or a subsection
        if (p.second.empty())
          {
            // set the found parameter in the destination argument
            try
              {
                prm.set(demangle(p.first), p.second.data());
              }
            catch (const ParameterHandler::ExcEntryUndeclared &entry_string)
              {
                // ignore undeclared entry assert
                AssertThrow(skip_undefined,
                            ParameterHandler::ExcEntryUndeclared(entry_string));
              }
          }
        else if (p.second.get_optional<std::string>("value"))
          {
            // set the found parameter in the destination argument
            try
              {
                prm.set(demangle(p.first), p.second.get<std::string>("value"));
              }
            catch (const ParameterHandler::ExcEntryUndeclared &entry_string)
              {
                // ignore undeclared entry assert
                AssertThrow(skip_undefined,
                            ParameterHandler::ExcEntryUndeclared(entry_string));
              }

            // this node might have sub-nodes in addition to "value", such as
            // "default_value", "documentation", etc. we might at some point
            // in the future want to make sure that if they exist that they
            // match the ones in the 'destination' tree
          }
        else if (p.second.get_optional<std::string>("alias"))
          {
            // it is an alias node. alias nodes are static and there is
            // nothing to do here (but the same applies as mentioned in the
            // comment above about the static nodes inside parameter nodes
          }
        else
          {
            try
              {
                // it must be a subsection
                prm.enter_subsection(demangle(p.first), !skip_undefined);
                read_xml_recursively(p.second,
                                     (current_path.empty() ?
                                        p.first :
                                        current_path + path_separator +
                                          p.first),
                                     path_separator,
                                     patterns,
                                     skip_undefined,
                                     prm);
                prm.leave_subsection();
              }
            catch (const ParameterHandler::ExcEntryUndeclared &entry_string)
              {
                // ignore undeclared entry assert
                AssertThrow(skip_undefined,
                            ParameterHandler::ExcEntryUndeclared(entry_string));
              }
          }
      }
  }

  // Recursively go through the @p source tree and collapse the nodes of the
  // format:
  //
  //"key":
  //      {
  //          "value"              : "val",
  //          "default_value"      : "...",
  //          "documentation"      : "...",
  //          "pattern"            : "...",
  //          "pattern_description": "..."
  //      },
  //
  // to
  //
  // "key" : "val";
  //
  // As an example a JSON file is shown. However, this function also works for
  // XML since both formats are build around the same BOOST data structures.
  //
  // This function is strongly based on read_xml_recursively().
  void
  recursively_compress_tree(boost::property_tree::ptree &source)
  {
    for (auto &p : source)
      {
        if (p.second.get_optional<std::string>("value"))
          {
            // save the value in a temporal variable
            const auto temp = p.second.get<std::string>("value");
            // clear node (children and value)
            p.second.clear();
            // set the correct value
            p.second.put_value<std::string>(temp);
          }
        else if (p.second.get_optional<std::string>("alias"))
          {
          }
        else
          {
            // it must be a subsection
            recursively_compress_tree(p.second);
          }
      }
  }

  void
  recursively_keep_non_default(const boost::property_tree::ptree &tree_in,
                               boost::property_tree::ptree       &tree_out)
  {
    for (const auto &p : tree_in)
      {
        if (is_parameter_node(p.second))
          {
            const std::string value = p.second.get<std::string>("value");

            if (value != p.second.get<std::string>("default_value"))
              tree_out.put_child(p.first, p.second);
          }
        else if (is_alias_node(p.second))
          {
            // nothing to do
          }
        else
          {
            boost::property_tree::ptree temp;

            if (const auto val = p.second.get_value_optional<std::string>())
              temp.put_value<std::string>(*val);

            recursively_keep_non_default(p.second, temp);

            if (temp.size() > 0)
              tree_out.put_child(p.first, temp);
          }
      }
  }
} // namespace



void
ParameterHandler::parse_input_from_xml(std::istream &in,
                                       const bool    skip_undefined)
{
  AssertThrow(in.fail() == false, ExcIO());
  // read the XML tree assuming that (as we
  // do in print_parameters(XML) it has only
  // a single top-level node called
  // "ParameterHandler"
  boost::property_tree::ptree single_node_tree;
  // This boost function will raise an exception if this is not a valid XML
  // file.
  read_xml(in, single_node_tree);

  // make sure there is a single top-level element
  // called "ParameterHandler"
  AssertThrow(single_node_tree.get_optional<std::string>("ParameterHandler"),
              ExcInvalidXMLParameterFile("There is no top-level XML element "
                                         "called \"ParameterHandler\"."));

  const std::size_t n_top_level_elements =
    std::distance(single_node_tree.begin(), single_node_tree.end());
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
           it != single_node_tree.end();
           ++it, ++entry_n)
        {
          top_level_message
            << "        " << it->first
            << (entry_n != n_top_level_elements - 1 ? "\n" : "");
        }

      // repeat assertion condition to make the printed version easier to read
      AssertThrow(n_top_level_elements == 1,
                  ExcInvalidXMLParameterFile(top_level_message.str()));
    }

  // read the child elements recursively
  const boost::property_tree::ptree &my_entries =
    single_node_tree.get_child("ParameterHandler");

  read_xml_recursively(
    my_entries, "", path_separator, patterns, skip_undefined, *this);
}


void
ParameterHandler::parse_input_from_json(std::istream &in,
                                        const bool    skip_undefined)
{
  AssertThrow(in.fail() == false, ExcIO());

  boost::property_tree::ptree node_tree;
  // This boost function will raise an exception if this is not a valid JSON
  // file.
  try
    {
      read_json(in, node_tree);
    }
  catch (const std::exception &e)
    {
      AssertThrow(
        false,
        ExcMessage(
          "The provided JSON file is not valid. Boost aborted with the "
          "following assert message: \n\n" +
          std::string(e.what())));
    }

  // The xml function is reused to read in the xml into the parameter file.
  // This function can only read mangled files. Therefore, we create a mangled
  // tree first.
  boost::property_tree::ptree node_tree_mangled;
  recursively_mangle_or_demangle(node_tree,
                                 node_tree_mangled,
                                 true /*do_mangle*/);
  read_xml_recursively(
    node_tree_mangled, "", path_separator, patterns, skip_undefined, *this);
}



void
ParameterHandler::clear()
{
  entries = std::make_unique<boost::property_tree::ptree>();
  entries_set_status.clear();
}



void
ParameterHandler::declare_entry(const std::string           &entry,
                                const std::string           &default_value,
                                const Patterns::PatternBase &pattern,
                                const std::string           &documentation,
                                const bool                   has_to_be_set)
{
  entries->put(get_current_full_path(entry) + path_separator + "value",
               default_value);
  entries->put(get_current_full_path(entry) + path_separator + "default_value",
               default_value);
  entries->put(get_current_full_path(entry) + path_separator + "documentation",
               documentation);

  // initialize with false
  const std::pair<bool, bool> set_status =
    std::pair<bool, bool>(has_to_be_set, false);
  entries_set_status.insert(
    std::pair<std::string, std::pair<bool, bool>>(get_current_full_path(entry),
                                                  set_status));

  patterns.reserve(patterns.size() + 1);
  patterns.emplace_back(pattern.clone());
  entries->put(get_current_full_path(entry) + path_separator + "pattern",
               static_cast<unsigned int>(patterns.size() - 1));
  // also store the description of
  // the pattern. we do so because we
  // may wish to export the whole
  // thing as XML or any other format
  // so that external tools can work
  // on the parameter file; in that
  // case, they will have to be able
  // to re-create the patterns as far
  // as possible
  entries->put(get_current_full_path(entry) + path_separator +
                 "pattern_description",
               patterns.back()->description());

  // as documented, do the default value checking at the very end
  AssertThrow(pattern.match(default_value),
              ExcValueDoesNotMatchPattern(default_value,
                                          pattern.description()));
}



void
ParameterHandler::add_action(
  const std::string                              &entry,
  const std::function<void(const std::string &)> &action,
  const bool                                      execute_action)
{
  actions.push_back(action);

  // get the current list of actions, if any
  boost::optional<std::string> current_actions =
    entries->get_optional<std::string>(get_current_full_path(entry) +
                                       path_separator + "actions");

  // if there were actions already associated with this parameter, add
  // the current one to it; otherwise, create a one-item list and use
  // that
  if (current_actions)
    {
      const std::string all_actions =
        current_actions.get() + "," +
        Utilities::int_to_string(actions.size() - 1);
      entries->put(get_current_full_path(entry) + path_separator + "actions",
                   all_actions);
    }
  else
    entries->put(get_current_full_path(entry) + path_separator + "actions",
                 Utilities::int_to_string(actions.size() - 1));

  if (execute_action)
    {
      const std::string default_value = entries->get<std::string>(
        get_current_full_path(entry) + path_separator + "default_value");
      action(default_value);
    }
}



void
ParameterHandler::declare_alias(const std::string &existing_entry_name,
                                const std::string &alias_name,
                                const bool         alias_is_deprecated)
{
  // see if there is anything to refer to already
  Assert(entries->get_optional<std::string>(
           get_current_full_path(existing_entry_name)),
         ExcMessage("You are trying to declare an alias entry <" + alias_name +
                    "> that references an entry <" + existing_entry_name +
                    ">, but the latter does not exist."));
  // then also make sure that what is being referred to is in
  // fact a parameter (not an alias or subsection)
  Assert(entries->get_optional<std::string>(
           get_current_full_path(existing_entry_name) + path_separator +
           "value"),
         ExcMessage("You are trying to declare an alias entry <" + alias_name +
                    "> that references an entry <" + existing_entry_name +
                    ">, but the latter does not seem to be a "
                    "parameter declaration."));


  // now also make sure that if the alias has already been
  // declared, that it is also an alias and refers to the same
  // entry
  if (entries->get_optional<std::string>(get_current_full_path(alias_name)))
    {
      Assert(entries->get_optional<std::string>(
               get_current_full_path(alias_name) + path_separator + "alias"),
             ExcMessage("You are trying to declare an alias entry <" +
                        alias_name +
                        "> but a non-alias entry already exists in this "
                        "subsection (i.e., there is either a preexisting "
                        "further subsection, or a parameter entry, with "
                        "the same name as the alias)."));
      Assert(entries->get<std::string>(get_current_full_path(alias_name) +
                                       path_separator + "alias") ==
               existing_entry_name,
             ExcMessage(
               "You are trying to declare an alias entry <" + alias_name +
               "> but an alias entry already exists in this "
               "subsection and this existing alias references a "
               "different parameter entry. Specifically, "
               "you are trying to reference the entry <" +
               existing_entry_name +
               "> whereas the existing alias references "
               "the entry <" +
               entries->get<std::string>(get_current_full_path(alias_name) +
                                         path_separator + "alias") +
               ">."));
    }

  entries->put(get_current_full_path(alias_name) + path_separator + "alias",
               existing_entry_name);

  if (alias_is_deprecated)
    mark_as_deprecated(alias_name);
  else
    mark_as_deprecated(alias_name, false);
}



void
ParameterHandler::mark_as_deprecated(const std::string &existing_entry_name,
                                     const bool         is_deprecated)
{
  // assert that the entry exists
  Assert(entries->get_optional<std::string>(
           get_current_full_path(existing_entry_name)),
         ExcMessage("You are trying to mark the entry <" + existing_entry_name +
                    "> as deprecated, but the entry does not exist."));

  // then also make sure that what is being referred to is in
  // fact a parameter or alias (not a subsection)
  Assert(
    entries->get_optional<std::string>(
      get_current_full_path(existing_entry_name) + path_separator + "value") ||
      entries->get_optional<std::string>(
        get_current_full_path(existing_entry_name) + path_separator + "alias"),
    ExcMessage("You are trying to mark the entry <" + existing_entry_name +
               "> as deprecated, but the entry does not seem to be a "
               "parameter declaration or a parameter alias."));

  entries->put(get_current_full_path(existing_entry_name) + path_separator +
                 "deprecation_status",
               is_deprecated ? "true" : "false");
}



void
ParameterHandler::enter_subsection(const std::string &subsection,
                                   const bool         create_path_if_needed)
{
  // if necessary create subsection
  const auto path_exists =
    entries->get_child_optional(get_current_full_path(subsection));

  if (!create_path_if_needed)
    {
      AssertThrow(path_exists,
                  ExcEntryUndeclared(get_current_full_path(subsection)));
    }

  if (!path_exists)
    entries->add_child(get_current_full_path(subsection),
                       boost::property_tree::ptree());

  // then enter it
  subsection_path.push_back(subsection);
}



void
ParameterHandler::leave_subsection()
{
  // assert there is a subsection that
  // we may leave
  Assert(subsection_path.size() != 0, ExcAlreadyAtTopLevel());

  if (subsection_path.size() > 0)
    subsection_path.pop_back();
}



bool
ParameterHandler::subsection_path_exists(
  const std::vector<std::string> &sub_path) const
{
  // Get full path to sub_path (i.e. prepend subsection_path to sub_path).
  std::vector<std::string> full_path(subsection_path);
  full_path.insert(full_path.end(), sub_path.begin(), sub_path.end());

  boost::optional<const boost::property_tree::ptree &> subsection(
    entries->get_child_optional(
      collate_path_string(path_separator, full_path)));

  // If subsection is boost::null (i.e. it does not exist)
  // or it exists as a parameter/alias node, return false.
  // Otherwise (i.e. it exists as a subsection node), return true.
  return !(!subsection || is_parameter_node(subsection.get()) ||
           is_alias_node(subsection.get()));
}



std::string
ParameterHandler::get(const std::string &entry_string) const
{
  // assert that the entry is indeed
  // declared
  if (boost::optional<std::string> value = entries->get_optional<std::string>(
        get_current_full_path(entry_string) + path_separator + "value"))
    return value.get();
  else
    {
      Assert(false, ExcEntryUndeclared(entry_string));
      return "";
    }
}



std::string
ParameterHandler::get(const std::vector<std::string> &entry_subsection_path,
                      const std::string              &entry_string) const
{
  // assert that the entry is indeed
  // declared
  if (boost::optional<std::string> value = entries->get_optional<std::string>(
        get_current_full_path(entry_subsection_path, entry_string) +
        path_separator + "value"))
    return value.get();
  else
    {
      Assert(false,
             ExcEntryUndeclared(demangle(
               get_current_full_path(entry_subsection_path, entry_string))));
      return "";
    }
}



long int
ParameterHandler::get_integer(const std::string &entry_string) const
{
  try
    {
      return Utilities::string_to_int(get(entry_string));
    }
  catch (...)
    {
      AssertThrow(false,
                  ExcMessage("Can't convert the parameter value <" +
                             get(entry_string) + "> for entry <" +
                             entry_string + "> to an integer."));
      return 0;
    }
}



long int
ParameterHandler::get_integer(
  const std::vector<std::string> &entry_subsection_path,
  const std::string              &entry_string) const
{
  try
    {
      return Utilities::string_to_int(get(entry_subsection_path, entry_string));
    }
  catch (...)
    {
      AssertThrow(false,
                  ExcMessage(
                    "Can't convert the parameter value <" +
                    get(entry_subsection_path, entry_string) + "> for entry <" +
                    demangle(get_current_full_path(entry_subsection_path,
                                                   entry_string)) +
                    "> to an integer."));
      return 0;
    }
}



double
ParameterHandler::get_double(const std::string &entry_string) const
{
  try
    {
      return Utilities::string_to_double(get(entry_string));
    }
  catch (...)
    {
      AssertThrow(false,
                  ExcMessage("Can't convert the parameter value <" +
                             get(entry_string) + "> for entry <" +
                             entry_string +
                             "> to a double precision variable."));
      return 0;
    }
}



double
ParameterHandler::get_double(
  const std::vector<std::string> &entry_subsection_path,
  const std::string              &entry_string) const
{
  try
    {
      return Utilities::string_to_double(
        get(entry_subsection_path, entry_string));
    }
  catch (...)
    {
      AssertThrow(false,
                  ExcMessage(
                    "Can't convert the parameter value <" +
                    get(entry_subsection_path, entry_string) + "> for entry <" +
                    demangle(get_current_full_path(entry_subsection_path,
                                                   entry_string)) +
                    "> to a double precision variable."));
      return 0;
    }
}



bool
ParameterHandler::get_bool(const std::string &entry_string) const
{
  const std::string s = get(entry_string);

  AssertThrow((s == "true") || (s == "false") || (s == "yes") || (s == "no"),
              ExcMessage("Can't convert the parameter value <" +
                         get(entry_string) + "> for entry <" + entry_string +
                         "> to a boolean."));
  if (s == "true" || s == "yes")
    return true;
  else
    return false;
}



bool
ParameterHandler::get_bool(
  const std::vector<std::string> &entry_subsection_path,
  const std::string              &entry_string) const
{
  const std::string s = get(entry_subsection_path, entry_string);

  AssertThrow((s == "true") || (s == "false") || (s == "yes") || (s == "no"),
              ExcMessage("Can't convert the parameter value <" +
                         get(entry_subsection_path, entry_string) +
                         "> for entry <" +
                         demangle(get_current_full_path(entry_subsection_path,
                                                        entry_string)) +
                         "> to a boolean."));
  if (s == "true" || s == "yes")
    return true;
  else
    return false;
}



void
ParameterHandler::set(const std::string &entry_string,
                      const std::string &new_value)
{
  // resolve aliases before looking up the correct entry
  std::string path = get_current_full_path(entry_string);
  if (entries->get_optional<std::string>(path + path_separator + "alias"))
    path = get_current_full_path(
      entries->get<std::string>(path + path_separator + "alias"));

  // get the node for the entry. if it doesn't exist, then we end up
  // in the else-branch below, which asserts that the entry is indeed
  // declared
  if (entries->get_optional<std::string>(path + path_separator + "value"))
    {
      // verify that the new value satisfies the provided pattern
      const unsigned int pattern_index =
        entries->get<unsigned int>(path + path_separator + "pattern");
      AssertThrow(patterns[pattern_index]->match(new_value),
                  ExcValueDoesNotMatchPattern(new_value,
                                              entries->get<std::string>(
                                                path + path_separator +
                                                "pattern_description")));

      // then also execute the actions associated with this
      // parameter (if any have been provided)
      const boost::optional<std::string> action_indices_as_string =
        entries->get_optional<std::string>(path + path_separator + "actions");
      if (action_indices_as_string)
        {
          std::vector<int> action_indices = Utilities::string_to_int(
            Utilities::split_string_list(action_indices_as_string.get()));
          for (const unsigned int index : action_indices)
            if (actions.size() >= index + 1)
              actions[index](new_value);
        }

      // finally write the new value into the database
      entries->put(path + path_separator + "value", new_value);

      auto map_iter = entries_set_status.find(path);
      if (map_iter != entries_set_status.end())
        map_iter->second = std::pair<bool, bool>(map_iter->second.first, true);
      else
        AssertThrow(false,
                    ExcMessage("Could not find parameter " + path +
                               " in map entries_set_status."));
    }
  else
    AssertThrow(false, ExcEntryUndeclared(entry_string));
}


void
ParameterHandler::set(const std::string &entry_string, const char *new_value)
{
  // simply forward
  set(entry_string, std::string(new_value));
}


void
ParameterHandler::set(const std::string &entry_string, const double new_value)
{
  std::ostringstream s;
  s << std::setprecision(16);
  s << new_value;

  // hand this off to the function that
  // actually sets the value as a string
  set(entry_string, s.str());
}



void
ParameterHandler::set(const std::string &entry_string, const long int new_value)
{
  std::ostringstream s;
  s << new_value;

  // hand this off to the function that
  // actually sets the value as a string
  set(entry_string, s.str());
}



void
ParameterHandler::set(const std::string &entry_string, const bool new_value)
{
  // hand this off to the function that
  // actually sets the value as a string
  set(entry_string, (new_value ? "true" : "false"));
}



std::ostream &
ParameterHandler::print_parameters(std::ostream     &out,
                                   const OutputStyle style) const
{
  AssertThrow(out.fail() == false, ExcIO());

  assert_validity_of_output_style(style);

  // Create entries copy and sort it, if needed.
  // In this way we ensure that the class state is never
  // modified by this function.
  boost::property_tree::ptree current_entries = *entries.get();

  // Sort parameters alphabetically, if needed.
  if ((style & KeepDeclarationOrder) == 0)
    {
      // Dive recursively into the subsections,
      // starting from the top level.
      recursively_sort_parameters(path_separator,
                                  std::vector<std::string>(),
                                  current_entries);
    }

  if ((style & KeepOnlyChanged) != 0)
    {
      boost::property_tree::ptree current_entries_without_default;
      recursively_keep_non_default(current_entries,
                                   current_entries_without_default);
      current_entries = current_entries_without_default;
    }

  // we'll have to print some text that is padded with spaces;
  // set the appropriate fill character, but also make sure that
  // we will restore the previous setting (and all other stream
  // flags) when we exit this function
  boost::io::ios_flags_saver            restore_flags(out);
  boost::io::basic_ios_fill_saver<char> restore_fill_state(out);
  out.fill(' ');

  // we treat XML and JSON is one step via BOOST, whereas all of the others are
  // done recursively in our own code. take care of the two special formats
  // first

  // explicitly compress the tree if requested
  if (((style & Short) != 0) && ((style & (XML | JSON)) != 0))
    {
      // modify the copy of the tree
      recursively_compress_tree(current_entries);
    }

  if ((style & XML) != 0)
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
      single_node_tree.add_child("ParameterHandler", current_entries);

      // set indentation character and indentation length
      boost::property_tree::xml_writer_settings<
        boost::property_tree::ptree::key_type>
        settings(' ', 2);

      write_xml(out, single_node_tree, settings);
      return out;
    }

  if ((style & JSON) != 0)
    {
      boost::property_tree::ptree current_entries_demangled;
      recursively_mangle_or_demangle(current_entries,
                                     current_entries_demangled,
                                     false /*do_mangle*/);
      write_json(out, current_entries_demangled);
      return out;
    }

  // for all of the other formats, print a preamble:
  if (((style & Short) != 0) && ((style & PRM) != 0))
    {
      // nothing to do
    }
  else if ((style & PRM) != 0)
    {
      out << "# Listing of Parameters" << std::endl
          << "# ---------------------" << std::endl;
    }
  else if ((style & LaTeX) != 0)
    {
      out << "\\subsection{Global parameters}" << std::endl;
      out << "\\label{parameters:global}" << std::endl;
      out << std::endl << std::endl;
    }
  else if ((style & Description) != 0)
    {
      out << "Listing of Parameters:" << std::endl << std::endl;
    }
  else
    {
      DEAL_II_NOT_IMPLEMENTED();
    }

  // dive recursively into the subsections
  recursively_print_parameters(
    current_entries,
    std::vector<std::string>(), // start at the top level
    style,
    0,
    out);

  return out;
}



void
ParameterHandler::print_parameters(
  const std::string                  &filename,
  const ParameterHandler::OutputStyle style) const
{
  std::string extension = filename.substr(filename.find_last_of('.') + 1);
  boost::algorithm::to_lower(extension);

  ParameterHandler::OutputStyle output_style = style;
  if (extension == "prm")
    output_style = style | PRM;
  else if (extension == "xml")
    output_style = style | XML;
  else if (extension == "json")
    output_style = style | JSON;
  else if (extension == "tex")
    output_style = style | LaTeX;

  std::ofstream out(filename);
  AssertThrow(out.fail() == false, ExcIO());
  print_parameters(out, output_style);
}



void
ParameterHandler::recursively_print_parameters(
  const boost::property_tree::ptree &tree,
  const std::vector<std::string>    &target_subsection_path,
  const OutputStyle                  style,
  const unsigned int                 indent_level,
  std::ostream                      &out) const
{
  AssertThrow(out.fail() == false, ExcIO());

  // this function should not be necessary for XML or JSON output...
  Assert(!(style & (XML | JSON)), ExcInternalError());

  const boost::property_tree::ptree &current_section =
    tree.get_child(collate_path_string(path_separator, target_subsection_path));

  unsigned int overall_indent_level = indent_level;

  const bool is_short = (style & Short) != 0;

  if ((style & PRM) != 0)
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
      for (const auto &p : current_section)
        if (is_parameter_node(p.second) == true)
          {
            longest_name  = std::max(longest_name, demangle(p.first).size());
            longest_value = std::max(longest_value,
                                     p.second.get<std::string>("value").size());
          }

      // print entries one by one
      bool first_entry = true;
      for (const auto &p : current_section)
        if (is_parameter_node(p.second) == true)
          {
            const std::string value = p.second.get<std::string>("value");

            // if there is documentation, then add an empty line
            // (unless this is the first entry in a subsection), print
            // the documentation, and then the actual entry; break the
            // documentation into readable chunks such that the whole
            // thing is at most 78 characters wide
            if (!is_short &&
                !p.second.get<std::string>("documentation").empty())
              {
                if (first_entry == false)
                  out << '\n';
                else
                  first_entry = false;

                const std::vector<std::string> doc_lines =
                  Utilities::break_text_into_lines(
                    p.second.get<std::string>("documentation"),
                    78 - overall_indent_level * 2 - 2);

                for (const auto &doc_line : doc_lines)
                  {
                    // Start with the comment start ('#'), padded with
                    // overall_indent_level*2 spaces at the front.
                    out << std::setw(overall_indent_level * 2) << "" << '#';

                    if (!doc_line.empty())
                      out << ' ' << doc_line;

                    out << '\n';
                  }
              }

            // Print the name and (if set) value of this entry. Ensure proper
            // padding with overall_indent_level*2 spaces at the front.
            out << std::setw(overall_indent_level * 2) << ""
                << "set " << demangle(p.first)
                << std::setw(longest_name - demangle(p.first).size() + 1) << ' '
                << '=';
            if (!value.empty())
              out << ' ' << value;

            // finally print the default value, but only if it differs
            // from the actual value
            if (!is_short &&
                value != p.second.get<std::string>("default_value"))
              {
                out << std::setw(longest_value - value.size() + 1) << ' '
                    << "# ";
                out << "default: "
                    << p.second.get<std::string>("default_value");
              }

            out << '\n';
          }
    }
  else if ((style & LaTeX) != 0)
    {
      auto escape = [](const std::string &input) {
        return Patterns::internal::escape(input, Patterns::PatternBase::LaTeX);
      };

      // if there are any parameters in this section then print them
      // as an itemized list
      const bool parameters_exist_here =
        std::any_of(current_section.begin(),
                    current_section.end(),
                    [](const boost::property_tree::ptree::value_type &p) {
                      return is_parameter_node(p.second) ||
                             is_alias_node(p.second);
                    });
      if (parameters_exist_here)
        {
          out << "\\begin{itemize}" << '\n';

          // print entries one by one
          for (const auto &p : current_section)
            if (is_parameter_node(p.second) == true)
              {
                const std::string value = p.second.get<std::string>("value");

                // print name
                out << "\\item {\\it Parameter name:} {\\tt "
                    << escape(demangle(p.first)) << "}\n"
                    << "\\phantomsection";
                {
                  // create label: labels are not to be escaped but
                  // mangled
                  std::string label = "parameters:";
                  for (const auto &path : target_subsection_path)
                    {
                      label.append(mangle(path));
                      label.append("/");
                    }
                  label.append(p.first);
                  // Backwards-compatibility. Output the label with and
                  // without escaping whitespace:
                  if (label.find("_20") != std::string::npos)
                    out << "\\label{"
                        << Utilities::replace_in_string(label, "_20", " ")
                        << "}\n";
                  out << "\\label{" << label << "}\n";
                }
                out << "\n\n";

                out << "\\index[prmindex]{" << escape(demangle(p.first))
                    << "}\n";
                out << "\\index[prmindexfull]{";
                for (const auto &path : target_subsection_path)
                  out << escape(path) << "!";
                out << escape(demangle(p.first)) << "}\n";

                // finally print value and default
                out << "{\\it Value:} " << escape(value) << "\n\n"
                    << '\n'
                    << "{\\it Default:} "
                    << escape(p.second.get<std::string>("default_value"))
                    << "\n\n"
                    << '\n';

                // if there is a documenting string, print it as well but
                // don't escape to allow formatting/formulas
                if (!is_short &&
                    !p.second.get<std::string>("documentation").empty())
                  out << "{\\it Description:} "
                      << p.second.get<std::string>("documentation")
                      << ((p.second.get_optional<std::string>(
                             "deprecation_status") &&
                           p.second.get<std::string>("deprecation_status") ==
                             "true") ?
                            " This parameter is deprecated.\n\n" :
                            "\n\n")
                      << '\n';
                if (!is_short)
                  {
                    // also output possible values, do not escape because the
                    // description internally will use LaTeX formatting
                    const unsigned int pattern_index =
                      p.second.get<unsigned int>("pattern");
                    const std::string desc_str =
                      patterns[pattern_index]->description(
                        Patterns::PatternBase::LaTeX);
                    out << "{\\it Possible values:} " << desc_str << '\n';
                  }
              }
            else if (is_alias_node(p.second) == true)
              {
                const std::string alias = p.second.get<std::string>("alias");

                // print name
                out << "\\item {\\it Parameter name:} {\\tt "
                    << escape(demangle(p.first)) << "}\n"
                    << "\\phantomsection";
                {
                  // create label: labels are not to be escaped but
                  // mangled
                  std::string label = "parameters:";
                  for (const auto &path : target_subsection_path)
                    {
                      label.append(mangle(path));
                      label.append("/");
                    }
                  label.append(p.first);
                  // Backwards-compatibility. Output the label with and
                  // without escaping whitespace:
                  if (label.find("_20") != std::string::npos)
                    out << "\\label{"
                        << Utilities::replace_in_string(label, "_20", " ")
                        << "}\n";
                  out << "\\label{" << label << "}\n";
                }
                out << "\n\n";

                out << "\\index[prmindex]{" << escape(demangle(p.first))
                    << "}\n";
                out << "\\index[prmindexfull]{";
                for (const auto &path : target_subsection_path)
                  out << escape(path) << "!";
                out << escape(demangle(p.first)) << "}\n";

                // finally print alias and indicate if it is deprecated
                out
                  << "This parameter is an alias for the parameter ``\\texttt{"
                  << escape(alias) << "}''."
                  << (p.second.get<std::string>("deprecation_status") ==
                          "true" ?
                        " Its use is deprecated." :
                        "")
                  << "\n\n"
                  << '\n';
              }
          out << "\\end{itemize}" << '\n';
        }
    }
  else if ((style & Description) != 0)
    {
      // first find out the longest entry name to be able to align the
      // equal signs
      std::size_t longest_name = 0;
      for (const auto &p : current_section)
        if (is_parameter_node(p.second) == true)
          longest_name = std::max(longest_name, demangle(p.first).size());

      // print entries one by one
      for (const auto &p : current_section)
        if (is_parameter_node(p.second) == true)
          {
            // print name and value
            out << std::setw(overall_indent_level * 2) << ""
                << "set " << demangle(p.first)
                << std::setw(longest_name - demangle(p.first).size() + 1) << " "
                << " = ";

            // print possible values:
            const unsigned int pattern_index =
              p.second.get<unsigned int>("pattern");
            const std::string full_desc_str =
              patterns[pattern_index]->description(Patterns::PatternBase::Text);
            const std::vector<std::string> description_str =
              Utilities::break_text_into_lines(
                full_desc_str, 78 - overall_indent_level * 2 - 2, '|');
            if (description_str.size() > 1)
              {
                out << '\n';
                for (const auto &description : description_str)
                  out << std::setw(overall_indent_level * 2 + 6) << ""
                      << description << '\n';
              }
            else if (description_str.empty() == false)
              out << "  " << description_str[0] << '\n';
            else
              out << '\n';

            // if there is a documenting string, print it as well
            if (!is_short &&
                p.second.get<std::string>("documentation").size() != 0)
              out << std::setw(overall_indent_level * 2 + longest_name + 10)
                  << ""
                  << "(" << p.second.get<std::string>("documentation") << ")"
                  << '\n';
          }
    }
  else
    {
      DEAL_II_NOT_IMPLEMENTED();
    }


  // if there was text before and there are sections to come, put two
  // newlines between the last entry and the first subsection
  {
    unsigned int n_parameters = 0;
    unsigned int n_sections   = 0;
    for (const auto &p : current_section)
      if (is_parameter_node(p.second) == true)
        ++n_parameters;
      else if (is_alias_node(p.second) == false)
        ++n_sections;

    if (((style & Description) == 0) && (!(((style & PRM) != 0) && is_short)) &&
        (n_parameters != 0) && (n_sections != 0))
      out << "\n\n";
  }

  // now transverse subsections tree
  for (const auto &p : current_section)
    if ((is_parameter_node(p.second) == false) &&
        (is_alias_node(p.second) == false))
      {
        // first print the subsection header
        if (((style & PRM) != 0) || ((style & Description) != 0))
          {
            out << std::setw(overall_indent_level * 2) << ""
                << "subsection " << demangle(p.first) << '\n';
          }
        else if ((style & LaTeX) != 0)
          {
            auto escape = [](const std::string &input) {
              return Patterns::internal::escape(input,
                                                Patterns::PatternBase::LaTeX);
            };

            out << '\n' << "\\subsection{Parameters in section \\tt ";

            // find the path to the current section so that we can
            // print it in the \subsection{...} heading
            for (const auto &path : target_subsection_path)
              out << escape(path) << "/";
            out << escape(demangle(p.first));

            out << "}" << '\n';
            out << "\\label{parameters:";
            for (const auto &path : target_subsection_path)
              out << mangle(path) << "/";
            out << p.first << "}";
            out << '\n';

            out << '\n';
          }
        else
          {
            DEAL_II_NOT_IMPLEMENTED();
          }

        // then the contents of the subsection
        const std::string        subsection     = demangle(p.first);
        std::vector<std::string> directory_path = target_subsection_path;
        directory_path.emplace_back(subsection);

        recursively_print_parameters(
          tree, directory_path, style, overall_indent_level + 1, out);

        if (is_short && ((style & PRM) != 0))
          {
            // write end of subsection.
            out << std::setw(overall_indent_level * 2) << ""
                << "end" << '\n';
          }
        else if ((style & PRM) != 0)
          {
            // write end of subsection. one blank line after each
            // subsection
            out << std::setw(overall_indent_level * 2) << ""
                << "end" << '\n'
                << '\n';

            // if this is a toplevel subsection, then have two
            // newlines
            if (overall_indent_level == 0)
              out << '\n';
          }
        else if ((style & Description) != 0)
          {
            // nothing to do
          }
        else if ((style & LaTeX) != 0)
          {
            // nothing to do
          }
        else
          {
            DEAL_II_NOT_IMPLEMENTED();
          }
      }
}



void
ParameterHandler::log_parameters(LogStream &out, const OutputStyle style)
{
  out.push("parameters");
  // dive recursively into the subsections
  log_parameters_section(out, style);

  out.pop();
}


void
ParameterHandler::log_parameters_section(LogStream        &out,
                                         const OutputStyle style)
{
  // Create entries copy and sort it, if needed.
  // In this way we ensure that the class state is never
  // modified by this function.
  boost::property_tree::ptree  sorted_entries;
  boost::property_tree::ptree *current_entries = entries.get();

  // Sort parameters alphabetically, if needed.
  if ((style & KeepDeclarationOrder) == 0)
    {
      sorted_entries  = *entries;
      current_entries = &sorted_entries;

      // Dive recursively into the subsections,
      // starting from the current level.
      recursively_sort_parameters(path_separator,
                                  subsection_path,
                                  sorted_entries);
    }

  const boost::property_tree::ptree &current_section =
    current_entries->get_child(get_current_path());

  // print entries one by one
  for (const auto &p : current_section)
    if (is_parameter_node(p.second) == true)
      out << demangle(p.first) << ": " << p.second.get<std::string>("value")
          << std::endl;

  // now transverse subsections tree
  for (const auto &p : current_section)
    if (is_parameter_node(p.second) == false)
      {
        out.push(demangle(p.first));
        enter_subsection(demangle(p.first));
        log_parameters_section(out, style);
        leave_subsection();
        out.pop();
      }
}



void
ParameterHandler::scan_line(std::string        line,
                            const std::string &input_filename,
                            const unsigned int current_line_n,
                            const bool         skip_undefined)
{
  // save a copy for some error messages
  const std::string original_line = line;

  // if there is a comment, delete it
  if (line.find('#') != std::string::npos)
    line.erase(line.find('#'), std::string::npos);

  // replace \t by space:
  while (line.find('\t') != std::string::npos)
    line.replace(line.find('\t'), 1, " ");

  // trim start and end:
  line = Utilities::trim(line);

  // if line is now empty: leave
  if (line.empty())
    {
      return;
    }
  // enter subsection
  else if (Utilities::match_at_string_start(line, "SUBSECTION ") ||
           Utilities::match_at_string_start(line, "subsection "))
    {
      // delete this prefix
      line.erase(0, std::string("subsection").size() + 1);

      const std::string subsection = Utilities::trim(line);

      // check whether subsection exists
      AssertThrow(skip_undefined || entries->get_child_optional(
                                      get_current_full_path(subsection)),
                  ExcNoSubsection(current_line_n,
                                  input_filename,
                                  demangle(get_current_full_path(subsection))));

      // subsection exists
      subsection_path.push_back(subsection);
    }
  // exit subsection
  else if (Utilities::match_at_string_start(line, "END") ||
           Utilities::match_at_string_start(line, "end"))
    {
      line.erase(0, 3);
      while ((line.size() > 0) && ((std::isspace(line[0])) != 0))
        line.erase(0, 1);

      AssertThrow(
        line.empty(),
        ExcCannotParseLine(current_line_n,
                           input_filename,
                           "Invalid content after 'end' or 'END' statement."));
      AssertThrow(subsection_path.size() != 0,
                  ExcCannotParseLine(current_line_n,
                                     input_filename,
                                     "There is no subsection to leave here."));
      leave_subsection();
    }
  // regular entry
  else if (Utilities::match_at_string_start(line, "SET ") ||
           Utilities::match_at_string_start(line, "set "))
    {
      // erase "set" statement
      line.erase(0, 4);

      std::string::size_type pos = line.find('=');
      AssertThrow(
        pos != std::string::npos,
        ExcCannotParseLine(current_line_n,
                           input_filename,
                           "Invalid format of 'set' or 'SET' statement."));

      // extract entry name and value and trim
      std::string entry_name = Utilities::trim(std::string(line, 0, pos));
      std::string entry_value =
        Utilities::trim(std::string(line, pos + 1, std::string::npos));

      // resolve aliases before we look up the entry. if necessary, print
      // a warning that the alias is deprecated
      std::string path = get_current_full_path(entry_name);
      if (entries->get_optional<std::string>(path + path_separator + "alias"))
        {
          if (entries->get<std::string>(path + path_separator +
                                        "deprecation_status") == "true")
            {
              std::cerr << "Warning in line <" << current_line_n
                        << "> of file <" << input_filename
                        << ">: You are using the deprecated spelling <"
                        << entry_name << "> of the parameter <"
                        << entries->get<std::string>(path + path_separator +
                                                     "alias")
                        << ">." << std::endl;
            }
          path = get_current_full_path(
            entries->get<std::string>(path + path_separator + "alias"));
        }

      // get the node for the entry. if it doesn't exist, then we end up
      // in the else-branch below, which asserts that the entry is indeed
      // declared
      if (entries->get_optional<std::string>(path + path_separator + "value"))
        {
          // if entry was declared: does it match the regex? if not, don't enter
          // it into the database exception: if it contains characters which
          // specify it as a multiple loop entry, then ignore content
          if (entry_value.find('{') == std::string::npos)
            {
              // verify that the new value satisfies the provided pattern
              const unsigned int pattern_index =
                entries->get<unsigned int>(path + path_separator + "pattern");
              AssertThrow(patterns[pattern_index]->match(entry_value),
                          ExcInvalidEntryForPattern(
                            current_line_n,
                            input_filename,
                            entry_value,
                            entry_name,
                            patterns[pattern_index]->description()));

              // then also execute the actions associated with this
              // parameter (if any have been provided)
              const boost::optional<std::string> action_indices_as_string =
                entries->get_optional<std::string>(path + path_separator +
                                                   "actions");
              if (action_indices_as_string)
                {
                  std::vector<int> action_indices =
                    Utilities::string_to_int(Utilities::split_string_list(
                      action_indices_as_string.get()));
                  for (const unsigned int index : action_indices)
                    if (actions.size() >= index + 1)
                      actions[index](entry_value);
                }
            }

          // finally write the new value into the database
          entries->put(path + path_separator + "value", entry_value);

          // record that the entry has been set manually
          auto map_iter = entries_set_status.find(path);
          if (map_iter != entries_set_status.end())
            map_iter->second =
              std::pair<bool, bool>(map_iter->second.first, true);
          else
            AssertThrow(false,
                        ExcMessage("Could not find parameter " + path +
                                   " in map entries_set_status."));

          // Check if the entry (or the resolved alias) is deprecated,
          // throw an exception if it is
          if (entries->get_optional<std::string>(path + path_separator +
                                                 "deprecation_status") &&
              entries->get<std::string>(path + path_separator +
                                        "deprecation_status") == "true")
            {
              AssertThrow(false,
                          ExcEntryIsDeprecated(current_line_n,
                                               input_filename,
                                               entry_name));
            }
        }
      else
        {
          AssertThrow(
            skip_undefined,
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
      line.erase(0, 7);
      while ((line.size() > 0) && (line[0] == ' '))
        line.erase(0, 1);

      // the remainder must then be a filename
      AssertThrow(line.size() != 0,
                  ExcCannotParseLine(current_line_n,
                                     input_filename,
                                     "The current line is an 'include' or "
                                     "'INCLUDE' statement, but it does not "
                                     "name a file for inclusion."));

      std::ifstream input(line);
      AssertThrow(input,
                  ExcCannotOpenIncludeStatementFile(current_line_n,
                                                    input_filename,
                                                    line));
      parse_input(input, line, "", skip_undefined);
    }
  else
    {
      AssertThrow(
        false,
        ExcCannotParseLine(current_line_n,
                           input_filename,
                           "The line\n\n"
                           "        <" +
                             original_line +
                             ">\n\n"
                             "could not be parsed: please check to "
                             "make sure that the file is not missing a "
                             "'set', 'include', 'subsection', or 'end' "
                             "statement."));
    }
}



std::size_t
ParameterHandler::memory_consumption() const
{
  // TODO: add to this an estimate of the memory in the property_tree
  return (MemoryConsumption::memory_consumption(subsection_path));
}



bool
ParameterHandler::operator==(const ParameterHandler &prm2) const
{
  if (patterns.size() != prm2.patterns.size())
    return false;

  for (unsigned int j = 0; j < patterns.size(); ++j)
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
  write_json(o1, *entries);
  write_json(o2, *prm2.entries);
  return (o1.str() == o2.str());
}



std::set<std::string>
ParameterHandler::get_entries_wrongly_not_set() const
{
  std::set<std::string> entries_wrongly_not_set;

  for (const auto &it : entries_set_status)
    if (it.second.first == true && it.second.second == false)
      entries_wrongly_not_set.insert(it.first);

  return entries_wrongly_not_set;
}



void
ParameterHandler::assert_that_entries_have_been_set() const
{
  const std::set<std::string> entries_wrongly_not_set =
    this->get_entries_wrongly_not_set();

  if (entries_wrongly_not_set.size() > 0)
    {
      std::string list_of_missing_parameters = "\n\n";
      for (const auto &it : entries_wrongly_not_set)
        list_of_missing_parameters += "  " + it + "\n";
      list_of_missing_parameters += "\n";

      AssertThrow(
        entries_wrongly_not_set.empty(),
        ExcMessage(
          "Not all entries of the parameter handler that were declared with "
          "`has_to_be_set = true` have been set. The following parameters " +
          list_of_missing_parameters +
          " have not been set. "
          "A possible reason might be that you did not add these parameter to "
          "the input file or that their spelling is not correct."));
    }
}



MultipleParameterLoop::MultipleParameterLoop()
  : n_branches(0)
{}



void
MultipleParameterLoop::parse_input(std::istream      &input,
                                   const std::string &filename,
                                   const std::string &last_line,
                                   const bool         skip_undefined)
{
  AssertThrow(input.fail() == false, ExcIO());

  // Note that (to avoid infinite recursion) we have to explicitly call the
  // base class version of parse_input and *not* a wrapper (which may be
  // virtual and lead us back here)
  ParameterHandler::parse_input(input, filename, last_line, skip_undefined);
  init_branches();
}



void
MultipleParameterLoop::loop(MultipleParameterLoop::UserClass &uc)
{
  for (unsigned int run_no = 0; run_no < n_branches; ++run_no)
    {
      // give create_new one-based numbers
      uc.create_new(run_no + 1);
      fill_entry_values(run_no);
      uc.run(*this);
    }
}



void
MultipleParameterLoop::init_branches()
{
  multiple_choices.clear();
  init_branches_current_section();

  // split up different values
  for (auto &multiple_choice : multiple_choices)
    multiple_choice.split_different_values();

  // finally calculate number of branches
  n_branches = 1;
  for (const auto &multiple_choice : multiple_choices)
    if (multiple_choice.type == Entry::variant)
      n_branches *= multiple_choice.different_values.size();

  // check whether array entries have the correct
  // number of entries
  for (const auto &multiple_choice : multiple_choices)
    if (multiple_choice.type == Entry::array)
      if (multiple_choice.different_values.size() != n_branches)
        std::cerr << "    The entry value" << std::endl
                  << "        " << multiple_choice.entry_value << std::endl
                  << "    for the entry named" << std::endl
                  << "        " << multiple_choice.entry_name << std::endl
                  << "    does not have the right number of entries for the "
                  << std::endl
                  << "        " << n_branches
                  << " variant runs that will be performed." << std::endl;


  // do a first run on filling the values to
  // check for the conformance with the regexp
  // (later on, this will be lost in the whole
  // other output)
  for (unsigned int i = 0; i < n_branches; ++i)
    fill_entry_values(i);
}



void
MultipleParameterLoop::init_branches_current_section()
{
  const boost::property_tree::ptree &current_section =
    entries->get_child(get_current_path());

  // check all entries in the present
  // subsection whether they are
  // multiple entries
  for (const auto &p : current_section)
    if (is_parameter_node(p.second) == true)
      {
        const std::string value = p.second.get<std::string>("value");
        if (value.find('{') != std::string::npos)
          multiple_choices.emplace_back(subsection_path,
                                        demangle(p.first),
                                        value);
      }

  // then loop over all subsections
  for (const auto &p : current_section)
    if (is_parameter_node(p.second) == false)
      {
        enter_subsection(demangle(p.first));
        init_branches_current_section();
        leave_subsection();
      }
}



void
MultipleParameterLoop::fill_entry_values(const unsigned int run_no)
{
  unsigned int possibilities = 1;

  std::vector<Entry>::iterator choice;
  for (choice = multiple_choices.begin(); choice != multiple_choices.end();
       ++choice)
    {
      const unsigned int selection =
        (run_no / possibilities) % choice->different_values.size();
      std::string entry_value;
      if (choice->type == Entry::variant)
        entry_value = choice->different_values[selection];
      else
        {
          if (run_no >= choice->different_values.size())
            {
              std::cerr
                << "The given array for entry <" << choice->entry_name
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
      subsection_path.swap(choice->subsection_path);
      set(choice->entry_name, entry_value);
      subsection_path.swap(choice->subsection_path);

      // move ahead if it was a variant entry
      if (choice->type == Entry::variant)
        possibilities *= choice->different_values.size();
    }
}



std::size_t
MultipleParameterLoop::memory_consumption() const
{
  std::size_t mem = ParameterHandler::memory_consumption();
  for (const auto &multiple_choice : multiple_choices)
    mem += multiple_choice.memory_consumption();

  return mem;
}



MultipleParameterLoop::Entry::Entry(const std::vector<std::string> &ssp,
                                    const std::string              &Name,
                                    const std::string              &Value)
  : subsection_path(ssp)
  , entry_name(Name)
  , entry_value(Value)
  , type(Entry::array)
{}



void
MultipleParameterLoop::Entry::split_different_values()
{
  // split string into three parts:
  // part before the opening "{",
  // the selection itself, final
  // part after "}"
  std::string prefix(entry_value, 0, entry_value.find('{'));
  std::string multiple(entry_value,
                       entry_value.find('{') + 1,
                       entry_value.rfind('}') - entry_value.find('{') - 1);
  std::string postfix(entry_value,
                      entry_value.rfind('}') + 1,
                      std::string::npos);
  // if array entry {{..}}: delete inner
  // pair of braces
  if (multiple[0] == '{')
    multiple.erase(0, 1);
  if (multiple.back() == '}')
    multiple.erase(multiple.size() - 1, 1);
  // erase leading and trailing spaces
  // in multiple
  while (std::isspace(multiple[0]) != 0)
    multiple.erase(0, 1);
  while (std::isspace(multiple.back()) != 0)
    multiple.erase(multiple.size() - 1, 1);

  // delete spaces around '|'
  while (multiple.find(" |") != std::string::npos)
    multiple.replace(multiple.find(" |"), 2, "|");
  while (multiple.find("| ") != std::string::npos)
    multiple.replace(multiple.find("| "), 2, "|");

  while (multiple.find('|') != std::string::npos)
    {
      different_values.push_back(
        prefix + std::string(multiple, 0, multiple.find('|')) + postfix);
      multiple.erase(0, multiple.find('|') + 1);
    }
  // make up the last selection ("while" broke
  // because there was no '|' any more
  different_values.push_back(prefix + multiple + postfix);
  // finally check whether this was a variant
  // entry ({...}) or an array ({{...}})
  if ((entry_value.find("{{") != std::string::npos) &&
      (entry_value.find("}}") != std::string::npos))
    type = Entry::array;
  else
    type = Entry::variant;
}


std::size_t
MultipleParameterLoop::Entry::memory_consumption() const
{
  return (MemoryConsumption::memory_consumption(subsection_path) +
          MemoryConsumption::memory_consumption(entry_name) +
          MemoryConsumption::memory_consumption(entry_value) +
          MemoryConsumption::memory_consumption(different_values) +
          sizeof(type));
}

DEAL_II_NAMESPACE_CLOSE
