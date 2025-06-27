// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/exceptions.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/path_search.h>
#include <deal.II/base/utilities.h>

#include <algorithm>
#include <cstdio>
#include <iostream>

DEAL_II_NAMESPACE_OPEN


std::map<std::string, std::vector<std::string>> PathSearch::path_lists;
std::map<std::string, std::vector<std::string>> PathSearch::suffix_lists;
std::string                                     PathSearch::empty("");

void
PathSearch::initialize_classes()
{
  std::vector<std::string> v;
  v.emplace_back();
  path_lists.insert(map_type("PARAMETER", v));

  /*
   * TODO: re-enable some sensible default paths. Maier, 2012
   */
  path_lists.insert(map_type("MESH", v));

  v.clear();
  v.emplace_back();
  v.emplace_back(".prm");
  suffix_lists.insert(map_type("PARAMETER", v));

  /*
   * TODO: "Would require linking with the deal.II libraries"? This .cc
   * file gets compiled into the library... maier, 2012
   */
  // We cannot use the GridIn class
  // to query the formats, since this
  // would require linking with the
  // deal.II libraries.
  v.clear();
  v.emplace_back();
  v.emplace_back(".inp");
  v.emplace_back(".xda");
  v.emplace_back(".dbmesh");
  v.emplace_back(".dat");
  v.emplace_back(".plt");
  v.emplace_back(".nc");
  v.emplace_back(".msh");
  suffix_lists.insert(map_type("MESH", v));
}

std::vector<std::string> &
PathSearch::get_path_list(const std::string &cls)
{
  if (path_lists.empty())
    initialize_classes();

  // Modified by Luca Heltai. If a class is not there, add it
  if (path_lists.count(cls) == 0)
    add_class(cls);

  // Assert(path_lists.count(cls) != 0, ExcNoClass(cls));
  Assert(path_lists.count(cls) != 0, ExcInternalError());

  return path_lists.find(cls)->second;
}


std::vector<std::string> &
PathSearch::get_suffix_list(const std::string &cls)
{
  // This is redundant. The constructor should have already called the
  // add_path function with the path_list bit...

  // Modified by Luca Heltai. If a class is not there, add it
  if (suffix_lists.count(cls) == 0)
    add_class(cls);

  // Assert(suffix_lists.count(cls) != 0, ExcNoClass(cls));
  Assert(suffix_lists.count(cls) != 0, ExcInternalError());

  return suffix_lists.find(cls)->second;
}


PathSearch::PathSearch(const std::string &cls, const unsigned int debug)
  : cls(cls)
  , my_path_list(get_path_list(cls))
  , my_suffix_list(get_suffix_list(cls))
  , debug(debug)
{}


std::string
PathSearch::find(const std::string &filename,
                 const std::string &suffix,
                 const char        *open_mode)
{
  std::vector<std::string>::const_iterator       path;
  const std::vector<std::string>::const_iterator endp = my_path_list.end();

  std::string real_name;

  if (debug > 2)
    deallog << "PathSearch[" << cls << "] " << my_path_list.size()
            << " directories " << std::endl;

  // Try to open file in the various directories we have
  for (path = my_path_list.begin(); path != endp; ++path)
    {
      // see if the file exists as given, i.e., with
      // the whole filename specified, including (possibly)
      // the suffix
      {
        real_name = *path + filename;
        if (debug > 1)
          deallog << "PathSearch[" << cls << "] trying " << real_name
                  << std::endl;
        std::FILE *fp = std::fopen(real_name.c_str(), open_mode);
        if (fp != nullptr)
          {
            if (debug > 0)
              deallog << "PathSearch[" << cls << "] opened " << real_name
                      << std::endl;
            std::fclose(fp);
            return real_name;
          }
      }

      // try again with the suffix appended, unless there is
      // no suffix
      if (!suffix.empty())
        {
          real_name = *path + filename + suffix;
          if (debug > 1)
            deallog << "PathSearch[" << cls << "] trying " << real_name
                    << std::endl;
          std::FILE *fp = std::fopen(real_name.c_str(), open_mode);
          if (fp != nullptr)
            {
              if (debug > 0)
                deallog << "PathSearch[" << cls << "] opened " << real_name
                        << std::endl;
              std::fclose(fp);
              return real_name;
            }
        }
    }
  AssertThrow(false, ExcFileNotFound(filename, cls));
  return "";
}

std::string
PathSearch::find(const std::string &filename, const char *open_mode)
{
  std::vector<std::string>::const_iterator       suffix;
  const std::vector<std::string>::const_iterator ends = my_suffix_list.end();

  if (debug > 2)
    deallog << "PathSearch[" << cls << "] " << my_path_list.size()
            << " directories " << my_suffix_list.size() << " suffixes"
            << std::endl;

  for (suffix = my_suffix_list.begin(); suffix != ends; ++suffix)
    {
      try
        {
          return find(filename, *suffix, open_mode);
        }
      catch (ExcFileNotFound &)
        {
          continue;
        }
    }
  AssertThrow(false, ExcFileNotFound(filename, cls));
  return "";
}


void
PathSearch::add_class(const std::string &cls)
{
  // Make sure standard classes are
  // initialized first
  if (path_lists.empty())
    initialize_classes();
  // Add empty path and empty suffix
  // for new class
  std::vector<std::string> v;
  v.emplace_back();
  path_lists.insert(map_type(cls, v));
  suffix_lists.insert(map_type(cls, v));
}


void
PathSearch::add_path(const std::string &path, Position pos)
{
  if (pos == back)
    my_path_list.push_back(path);
  else if (pos == front)
    my_path_list.insert(my_path_list.begin(), path);
  else if (pos == after_none)
    {
      std::vector<std::string>::iterator i =
        std::find(my_path_list.begin(), my_path_list.end(), empty);
      if (i != my_path_list.end())
        ++i;
      my_path_list.insert(i, path);
    }
}


void
PathSearch::add_suffix(const std::string &suffix, Position pos)
{
  if (pos == back)
    my_suffix_list.push_back(suffix);
  else if (pos == front)
    my_suffix_list.insert(my_suffix_list.begin(), suffix);
  else if (pos == after_none)
    {
      std::vector<std::string>::iterator i =
        std::find(my_suffix_list.begin(), my_suffix_list.end(), empty);
      if (i != my_suffix_list.end())
        ++i;
      my_suffix_list.insert(i, suffix);
    }
}



DEAL_II_NAMESPACE_CLOSE
