//-----------------------------------------------------------
//
//    Copyright (C) 2017 - 2020 by the deal.II authors
//
//    This file is part of the deal.II library.
//
//    The deal.II library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE.md at
//    the top level directory of deal.II.
//
//-----------------------------------------------------------

#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/path_search.h>
#include <deal.II/base/utilities.h>

#include <boost/core/demangle.hpp>

#include <fstream>

DEAL_II_NAMESPACE_OPEN


// Static empty class list
std::vector<SmartPointer<ParameterAcceptor>> ParameterAcceptor::class_list;
// Static parameter handler
ParameterHandler ParameterAcceptor::prm;

ParameterAcceptor::ParameterAcceptor(const std::string &name)
  : acceptor_id(class_list.size())
  , section_name(name)
{
  SmartPointer<ParameterAcceptor> pt(
    this, boost::core::demangle(typeid(*this).name()).c_str());
  class_list.push_back(pt);
}


ParameterAcceptor::~ParameterAcceptor()
{
  class_list[acceptor_id] = nullptr;
}

std::string
ParameterAcceptor::get_section_name() const
{
  return (!section_name.empty() ? section_name :
                                  boost::core::demangle(typeid(*this).name()));
}


void
ParameterAcceptor::initialize(
  const std::string &                 filename,
  const std::string &                 output_filename,
  const ParameterHandler::OutputStyle output_style_for_output_filename,
  ParameterHandler &                  prm,
  const ParameterHandler::OutputStyle output_style_for_filename)
{
  declare_all_parameters(prm);
  if (!filename.empty())
    {
      try
        {
          prm.parse_input(filename);
        }
      catch (const dealii::PathSearch::ExcFileNotFound &)
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
  AssertThrow(input_stream, ExcIO());
  declare_all_parameters(prm);
  prm.parse_input(input_stream);
  parse_all_parameters(prm);
}

void
ParameterAcceptor::clear()
{
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
    if (instance != nullptr)
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
    if (instance != nullptr)
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
  Assert(acceptor_id < class_list.size(), ExcInternalError());
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
      for (int i = acceptor_id - 1; i >= 0; --i)
        if (class_list[i] != nullptr)
          {
            bool has_trailing = class_list[i]->get_section_name().back() == sep;
            auto previous_path = class_list[i]->get_section_path();

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
  return sections;
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



DEAL_II_NAMESPACE_CLOSE
