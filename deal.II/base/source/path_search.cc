//---------------------------------------------------------------------------
//      $Id$   
//    Version: $Name$
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <base/path_search.h>
#include <base/logstream.h>

#include <iostream>
#include <algorithm>

//TODO:[GK] Clean up open functions, reuse code!

std::map<std::string, std::vector<std::string> > PathSearch::path_lists;
std::map<std::string, std::vector<std::string> > PathSearch::suffix_lists;
std::string PathSearch::empty("");

void
PathSearch::initialize_classes()
{
  std::vector<std::string> v;
  v.push_back(empty);
  path_lists.insert(map_type(std::string("PARAMETER"), v));

  v.push_back(std::string(DEAL_II_PATH "/lib/meshes/"));
  path_lists.insert(map_type(std::string("MESH"), v));

  v.clear();
  v.push_back(empty);
  v.push_back(std::string(".prm"));
  suffix_lists.insert(map_type(std::string("PARAMETER"), v));
  
  v.clear();
  v.push_back(empty);
  v.push_back(std::string(".inp"));
  v.push_back(std::string(".xda"));
  v.push_back(std::string(".dbmesh"));
  suffix_lists.insert(map_type(std::string("MESH"), v));
}

std::vector<std::string>&
PathSearch::get_path_list(const std::string& cls)
{
  if (path_lists.empty())
    initialize_classes();

  Assert(path_lists.count(cls) != 0, ExcNoClass(cls));
  
  return path_lists.find(cls)->second;
}


std::vector<std::string>&
PathSearch::get_suffix_list(const std::string& cls)
{
  Assert(suffix_lists.count(cls) != 0, ExcNoClass(cls));
  
  return suffix_lists.find(cls)->second;
}


PathSearch::PathSearch(const std::string& cls,
		       const unsigned int debug)
		:
		cls(cls),
		my_path_list(get_path_list(cls)),
		my_suffix_list(get_suffix_list(cls)),
		debug(debug)
{}


std::istream&
PathSearch::open (const std::string& filename,
		  const std::string& suffix)
{
  std::string name;
  std::vector<std::string>::const_iterator path;
  const std::vector<std::string>::const_iterator endp = my_path_list.end();

  if (debug > 2)
    deallog << "PathSearch[" << cls << "] "
	    << my_path_list.size() << " directories "
	    << std::endl;
  
				   // Try without suffix first
  for (path = my_path_list.begin(); path != endp; ++path)
    {
      name = *path + filename;
      if (debug > 1)
	deallog << "PathSearch[" << cls << "] trying "
		<< name << std::endl;
      stream.reset(new std::ifstream(name.c_str()));
      if (stream->is_open())
	{
	  if (debug > 0)
	    deallog << "PathSearch[" << cls << "] opened "
		    << name << std::endl;
	  return *stream;
	}
    }

				   // Now try with given suffix
  for (path = my_path_list.begin(); path != endp; ++path)
    {
      name = *path + filename + suffix;
      if (debug > 1)
	deallog << "PathSearch[" << cls << "] trying "
		<< name << std::endl;
      stream.reset(new std::ifstream(name.c_str()));
      if (stream->is_open())
	{
	  if (debug > 0)
	    deallog << "PathSearch[" << cls << "] opened "
		    << name << std::endl;
	  return *stream;
	}
    }
  AssertThrow(false, ExcFileNotFound(filename, cls));
  return *stream;
}


std::istream&
PathSearch::open (const std::string& filename)
{
  std::string name;
  std::vector<std::string>::const_iterator suffix;
  std::vector<std::string>::const_iterator path;
  const std::vector<std::string>::const_iterator ends = my_suffix_list.end();
  const std::vector<std::string>::const_iterator endp = my_path_list.end();

  if (debug > 2)
    deallog << "PathSearch[" << cls << "] "
	    << my_path_list.size() << " directories "
	    << my_suffix_list.size() << " suffixes"
	    << std::endl;
  
  for (suffix = my_suffix_list.begin(); suffix != ends; ++suffix)
    {
      for (path = my_path_list.begin(); path != endp; ++path)
	{
	  name = *path + filename + *suffix;
	  if (debug > 1)
	    deallog << "PathSearch[" << cls << "] trying "
		    << name << std::endl;
	  stream.reset(new std::ifstream(name.c_str()));
	  if (stream->is_open())
	    {
	      if (debug > 0)
		deallog << "PathSearch[" << cls << "] opened "
			<< name << std::endl;
	      return *stream;
	    }
	}
    }
  AssertThrow(false, ExcFileNotFound(filename, cls));
  return *stream;
}


void
PathSearch::add_class (const std::string& cls)
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
PathSearch::add_path (const std::string& path,
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
PathSearch::add_suffix (const std::string& suffix,
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


