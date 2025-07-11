// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/path_search.h>
#include <deal.II/base/utilities.h>

#include <boost/core/demangle.hpp>

#include <fstream>
#include <set>


DEAL_II_NAMESPACE_OPEN


namespace internal
{
  struct ParameterAcceptorCompare
  {
    bool
    operator()(const ParameterAcceptor *p1, const ParameterAcceptor *p2) const
    {
      return p1->get_acceptor_id() < p2->get_acceptor_id();
    }
  };
} // namespace internal


// Static mutex for the class list
std::mutex ParameterAcceptor::class_list_mutex;

// Static empty class set
std::set<ParameterAcceptor *, internal::ParameterAcceptorCompare>
  ParameterAcceptor::class_list;

// Static parameter handler
ParameterHandler ParameterAcceptor::prm;

ParameterAcceptor::ParameterAcceptor(const std::string &name)
  : acceptor_id(ParameterAcceptor::get_next_free_id())
  , section_name(name)
{
  std::lock_guard<std::mutex> l(class_list_mutex);
  class_list.insert(this);
}



ParameterAcceptor::~ParameterAcceptor()
{
  std::lock_guard<std::mutex> l(class_list_mutex);
  // Notice that it is possible that the class is no longer in the static list.
  // This happens when the clear() method has been called. erase() does the
  // righy thing anyway by only removing this class if it's still in the list.
  class_list.erase(this);
}



std::string
ParameterAcceptor::get_section_name() const
{
  return (!section_name.empty() ? section_name :
                                  boost::core::demangle(typeid(*this).name()));
}



void
ParameterAcceptor::initialize(
  const std::string                  &filename,
  const std::string                  &output_filename,
  const ParameterHandler::OutputStyle output_style_for_output_filename,
  ParameterHandler                   &prm,
  const ParameterHandler::OutputStyle output_style_for_filename)
{
  declare_all_parameters(prm);
  if (!filename.empty())
    {
      try
        {
          prm.parse_input(filename);
        }
      catch (const dealii::ExcFileNotOpen &)
        {
          prm.print_parameters(filename, output_style_for_filename);
          AssertThrow(false,
                      ExcMessage("You specified <" + filename + "> as input " +
                                 "parameter file, but it does not exist. " +
                                 "We created it for you."));
        }
    }

  if (!output_filename.empty())
    prm.print_parameters(output_filename, output_style_for_output_filename);

  // Finally do the parsing.
  parse_all_parameters(prm);
}



void
ParameterAcceptor::initialize(std::istream &input_stream, ParameterHandler &prm)

{
  AssertThrow(input_stream.fail() == false, ExcIO());
  declare_all_parameters(prm);
  prm.parse_input(input_stream);
  parse_all_parameters(prm);
}



void
ParameterAcceptor::clear()
{
  std::lock_guard<std::mutex> l(class_list_mutex);
  class_list.clear();
  prm.clear();
}



void
ParameterAcceptor::declare_parameters(ParameterHandler &)
{}



void
ParameterAcceptor::parse_parameters(ParameterHandler &)
{}



void
ParameterAcceptor::parse_all_parameters(ParameterHandler &prm)
{
  for (const auto &instance : class_list)
    {
      instance->enter_my_subsection(prm);
      instance->parse_parameters(prm);
      instance->parse_parameters_call_back();
      instance->leave_my_subsection(prm);
    }
}



void
ParameterAcceptor::declare_all_parameters(ParameterHandler &prm)
{
  for (const auto &instance : class_list)
    {
      instance->enter_my_subsection(prm);
      instance->declare_parameters(prm);
      instance->declare_parameters_call_back();
      instance->leave_my_subsection(prm);
    }
}



std::vector<std::string>
ParameterAcceptor::get_section_path() const
{
  const auto my_section_name = get_section_name();
  const bool is_absolute     = (my_section_name.front() == sep);

  std::vector<std::string> sections =
    Utilities::split_string_list(my_section_name, sep);

  // Split string list removes trailing empty strings, but not
  // preceding ones. Make sure that if we had an absolute path,
  // we don't store as first section the empty string.
  if (is_absolute)
    sections.erase(sections.begin());
  else
    {
      // If we have a relative path, we prepend the path of the previous class
      // to ours. This is tricky. If the previous class has a path with a
      // trailing /, then the full path is used, else only the path except the
      // last one
      for (auto acceptor_it = class_list.rbegin();
           acceptor_it != class_list.rend();
           ++acceptor_it)
        {
          const auto *const acceptor = *acceptor_it;
          if (acceptor->get_acceptor_id() >= get_acceptor_id())
            continue;
          bool has_trailing  = acceptor->get_section_name().back() == sep;
          auto previous_path = acceptor->get_section_path();

          // See if we need to remove last piece of the path
          if ((previous_path.size() > 0) && has_trailing == false)
            previous_path.resize(previous_path.size() - 1);

          sections.insert(sections.begin(),
                          previous_path.begin(),
                          previous_path.end());
          // Exit the for cycle
          break;
        }
    }
  // Finally, insert the remaining subsections
  sections.insert(sections.end(), subsections.begin(), subsections.end());
  return sections;
}



void
ParameterAcceptor::enter_subsection(const std::string &subsection)
{
  AssertThrow(subsection.find(sep) == std::string::npos,
              ExcMessage(
                "A subsection name cannot contain the special character '/'"));

  AssertThrow(subsection != "",
              ExcMessage("Cannot create an empty subsection."));

  subsections.push_back(subsection);
}



void
ParameterAcceptor::leave_subsection()
{
  AssertThrow(subsections.size() > 0,
              ExcMessage("There is no subsection to leave here."));
  subsections.pop_back();
}



void
ParameterAcceptor::enter_my_subsection(
  ParameterHandler &prm = ParameterAcceptor::prm)
{
  const auto sections = get_section_path();
  for (const auto &sec : sections)
    {
      prm.enter_subsection(sec);
    }
}



void
ParameterAcceptor::leave_my_subsection(
  ParameterHandler &prm = ParameterAcceptor::prm)
{
  const auto sections = get_section_path();
  for (unsigned int i = 0; i < sections.size(); ++i)
    {
      prm.leave_subsection();
    }
}



inline unsigned int
ParameterAcceptor::get_acceptor_id() const
{
  return acceptor_id;
}



unsigned int
ParameterAcceptor::get_next_free_id()
{
  static std::mutex           id_mutex;
  std::lock_guard<std::mutex> lock(id_mutex);
  static int                  current_id = 0;
  return current_id++;
}

DEAL_II_NAMESPACE_CLOSE
