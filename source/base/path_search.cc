// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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


#include <deal.II/base/path_search.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>

#include <iostream>
#include <cstdio>
#include <algorithm>

DEAL_II_NAMESPACE_OPEN


std::map<std::string, std::vector<std::string> > PathSearch::path_lists;
std::map<std::string, std::vector<std::string> > PathSearch::suffix_lists;
std::string PathSearch::empty("");

void
PathSearch::initialize_classes()
{
  std::vector<std::string> v;
  v.push_back(empty);
  path_lists.insert(map_type(std::string("PARAMETER"), v));

  /*
   * TODO: reenable some sensible default paths. Maier, 2012
   */
  path_lists.insert(map_type(std::string("MESH"), v));

  v.clear();
  v.push_back(empty);
  v.push_back(std::string(".prm"));
  suffix_lists.insert(map_type(std::string("PARAMETER"), v));

  /*
   * TODO: "Would require linking with the deal.II libraries"? This .cc
   * file gets compiled into the library... maier, 2012
   */
  // We cannot use the GridIn class
  // to query the formats, since this
  // would require linking with the
  // deal.II libraries.
  v.clear();
  v.push_back(empty);
  v.push_back(std::string(".inp"));
  v.push_back(std::string(".xda"));
  v.push_back(std::string(".dbmesh"));
  v.push_back(std::string(".dat"));
  v.push_back(std::string(".plt"));
  v.push_back(std::string(".nc"));
  v.push_back(std::string(".msh"));
  suffix_lists.insert(map_type(std::string("MESH"), v));
}

std::vector<std::string> &
PathSearch::get_path_list(const std::string &cls)
{
  if (path_lists.empty())
    initialize_classes();

  // Modified by Luca Heltai. If a class is not there, add it
  if (path_lists.count(cls) == 0) add_class(cls);

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
  if (suffix_lists.count(cls) == 0) add_class(cls);

  // Assert(suffix_lists.count(cls) != 0, ExcNoClass(cls));
  Assert(suffix_lists.count(cls) != 0, ExcInternalError());

  return suffix_lists.find(cls)->second;
}


PathSearch::PathSearch(const std::string &cls,
                       const unsigned int debug)
  :
  cls(cls),
  my_path_list(get_path_list(cls)),
  my_suffix_list(get_suffix_list(cls)),
  debug(debug)
{}


std::string
PathSearch::find (const std::string &filename,
                  const std::string &suffix,
                  const char *open_mode)
{
  std::vector<std::string>::const_iterator path;
  const std::vector<std::string>::const_iterator endp = my_path_list.end();

  std::string real_name;

  if (debug > 2)
    deallog << "PathSearch[" << cls << "] "
            << my_path_list.size() << " directories "
            << std::endl;

  // Try to open file in the various directories we have
  for (path = my_path_list.begin(); path != endp; ++path)
    {
      // see if the file exists as given, i.e., with
      // the whole filename specified, including (possibly)
      // the suffix
      {
        real_name = *path + filename;
        if (debug > 1)
          deallog << "PathSearch[" << cls << "] trying "
                  << real_name << std::endl;
        FILE *fp = fopen(real_name.c_str(), open_mode);
        if (fp != 0)
          {
            if (debug > 0)
              deallog << "PathSearch[" << cls << "] opened "
                      << real_name << std::endl;
            fclose(fp);
            return real_name;
          }
      }

      // try again with the suffix appended, unless there is
      // no suffix
      if (suffix != "")
        {
          real_name = *path + filename + suffix;
          if (debug > 1)
            deallog << "PathSearch[" << cls << "] trying "
                    << real_name << std::endl;
          FILE *fp = fopen(real_name.c_str(), open_mode);
          if (fp != 0)
            {
              if (debug > 0)
                deallog << "PathSearch[" << cls << "] opened "
                        << real_name << std::endl;
              fclose(fp);
              return real_name;
            }
        }
    }
  AssertThrow(false, ExcFileNotFound(filename, cls));
  return std::string("");
}

std::string
PathSearch::find (const std::string &filename,
                  const char *open_mode)
{
  std::vector<std::string>::const_iterator suffix;
  const std::vector<std::string>::const_iterator ends = my_suffix_list.end();

  if (debug > 2)
    deallog << "PathSearch[" << cls << "] "
            << my_path_list.size() << " directories "
            << my_suffix_list.size() << " suffixes"
            << std::endl;

  for (suffix = my_suffix_list.begin(); suffix != ends; ++suffix)
    {
      try
        {
          return find(filename, *suffix, open_mode);
        }
      catch (ExcFileNotFound)
        {
          continue;
        }

    }
  AssertThrow(false, ExcFileNotFound(filename, cls));
  return std::string("");
}


void
PathSearch::add_class (const std::string &cls)
{
  // Make sure standard classes are
  // initialized first
  if (path_lists.empty())
    initialize_classes();
  // Add empty path and empty suffix
  // for new class
  std::vector<std::string> v;
  v.push_back(empty);
  path_lists.insert(map_type(cls, v));
  suffix_lists.insert(map_type(cls, v));
}


void
PathSearch::add_path (const std::string &path,
                      Position pos)
{
  if (pos == back)
    my_path_list.push_back(path);
  else if (pos == front)
    my_path_list.insert(my_path_list.begin(), path);
  else if (pos == after_none)
    {
      std::vector<std::string>::iterator
      i = std::find(my_path_list.begin(), my_path_list.end(), empty);
      if (i != my_path_list.end())
        ++i;
      my_path_list.insert(i, path);
    }
}


void
PathSearch::add_suffix (const std::string &suffix,
                        Position pos)
{
  if (pos == back)
    my_suffix_list.push_back(suffix);
  else if (pos == front)
    my_suffix_list.insert(my_suffix_list.begin(), suffix);
  else if (pos == after_none)
    {
      std::vector<std::string>::iterator
      i = std::find(my_suffix_list.begin(), my_suffix_list.end(), empty);
      if (i != my_suffix_list.end())
        ++i;
      my_suffix_list.insert(i, suffix);
    }
}



DEAL_II_NAMESPACE_CLOSE
