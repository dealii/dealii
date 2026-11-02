// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// XML round trips execute the destination's actions, including when its action
// registration order differs from that recorded in the XML document.

#include <deal.II/base/parameter_handler.h>

#include <sstream>
#include <string>
#include <vector>

#include "../tests.h"

void
check(const ParameterHandler::OutputStyle style)
{
  ParameterHandler source;
  source.declare_entry("left", "11", Patterns::Integer());
  source.declare_entry("right", "22", Patterns::Integer());
  source.enter_subsection("nested");
  source.declare_entry("value", "33", Patterns::Integer());
  source.leave_subsection();

  unsigned int source_calls = 0;
  source.add_action(
    "left", [&](const std::string &) { ++source_calls; }, false);
  source.add_action(
    "right", [&](const std::string &) { ++source_calls; }, false);
  source.enter_subsection("nested");
  source.add_action(
    "value", [&](const std::string &) { ++source_calls; }, false);
  source.leave_subsection();

  std::stringstream xml;
  source.print_parameters(xml, style);
  AssertThrow(source_calls == 0, ExcInternalError());

  ParameterHandler destination;
  destination.declare_entry("right", "0", Patterns::Integer());
  destination.declare_entry("left", "0", Patterns::Integer());
  destination.enter_subsection("nested");
  destination.declare_entry("value", "0", Patterns::Integer());
  destination.leave_subsection();

  std::vector<std::string> left_calls;
  std::vector<std::string> right_calls;
  std::vector<std::string> nested_calls;
  destination.add_action(
    "right",
    [&](const std::string &value) { right_calls.push_back(value); },
    false);
  destination.add_action(
    "left",
    [&](const std::string &value) { left_calls.push_back("first:" + value); },
    false);
  destination.add_action(
    "left",
    [&](const std::string &value) { left_calls.push_back("second:" + value); },
    false);
  destination.enter_subsection("nested");
  destination.add_action(
    "value",
    [&](const std::string &value) { nested_calls.push_back(value); },
    false);
  destination.leave_subsection();

  destination.parse_input_from_xml(xml);

  AssertThrow(source_calls == 0, ExcInternalError());
  AssertThrow((left_calls == std::vector<std::string>{"first:11", "second:11"}),
              ExcInternalError());
  AssertThrow((right_calls == std::vector<std::string>{"22"}),
              ExcInternalError());
  AssertThrow((nested_calls == std::vector<std::string>{"33"}),
              ExcInternalError());
  AssertThrow(destination.get_integer("left") == 11, ExcInternalError());
  AssertThrow(destination.get_integer("right") == 22, ExcInternalError());
  destination.enter_subsection("nested");
  AssertThrow(destination.get_integer("value") == 33, ExcInternalError());
  destination.leave_subsection();
  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();
  check(ParameterHandler::XML);
  check(ParameterHandler::ShortXML);
}
