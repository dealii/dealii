// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check use of path in ParameterHandler::Get

#include <deal.II/base/parameter_handler.h>

#include <sstream>

#include "../tests.h"

void
log_entry(const ParameterHandler        &prm,
          const std::vector<std::string> path,
          const std::string             &entry)
{
  for (unsigned int i = 0; i < path.size(); ++i)
    {
      deallog << path[i] << '.';
    }

  deallog << entry << " = ";
  deallog << prm.get(path, entry) << std::endl;
}

void
log_exception(const StandardExceptions::ExcMessage &exc)
{
  std::ostringstream string_stream;
  exc.print_info(string_stream);
  std::string exc_info = string_stream.str();

  const std::string whitespace_chars = " \t\r\n\v\f";
  exc_info.erase(0, exc_info.find_first_not_of(whitespace_chars));
  exc_info.erase(exc_info.find_last_not_of(whitespace_chars) + 1);

  deallog << exc_info << std::endl;
}

void
test_basic()
{
  ParameterHandler prm;
  prm.declare_entry("a", "1");
  prm.enter_subsection("x");
  prm.declare_entry("a", "2");
  prm.enter_subsection("y");
  prm.declare_entry("a", "3");
  prm.declare_entry("b", "4");
  prm.leave_subsection();
  prm.enter_subsection("v");
  prm.enter_subsection("v");
  prm.enter_subsection("v");
  prm.declare_entry("a", "5");
  prm.declare_entry("b", "6");
  prm.leave_subsection();
  prm.declare_entry("a", "7");
  prm.leave_subsection();
  prm.leave_subsection();
  prm.enter_subsection("z");
  prm.declare_entry("a", "8");
  prm.leave_subsection();
  prm.leave_subsection();
  prm.declare_entry("b", "9");

  /*
    prm should now look like this

    # Listing of Parameters
    # ---------------------
    set a = 1
    set b = 9

    subsection x
      set a = 2


      subsection v
        subsection v
          set a = 7


          subsection v
            set a = 5
            set b = 6
          end

        end

      end

      subsection y
        set a = 3
        set b = 4
      end

      subsection z
        set a = 8
      end

    end

  */

  // reading from top level
  std::vector<std::string> path = {};
  log_entry(prm, path, "a");
  path.push_back("x");
  log_entry(prm, path, "a");
  path.push_back("y");
  log_entry(prm, path, "a");
  log_entry(prm, path, "b");
  path.pop_back();
  path.push_back("v");
  path.push_back("v");
  path.push_back("v");
  log_entry(prm, path, "a");
  log_entry(prm, path, "b");
  path.pop_back();
  log_entry(prm, path, "a");
  path.pop_back();
  path.pop_back();
  path.push_back("z");
  log_entry(prm, path, "a");
  path.pop_back();
  path.pop_back();
  log_entry(prm, path, "b");

  // reading from other subsections
  prm.enter_subsection("x");
  deallog << "in x" << std::endl;
  log_entry(prm, path, "a");
  path.push_back("y");
  log_entry(prm, path, "a");
  log_entry(prm, path, "b");
  path.pop_back();
  path.push_back("v");
  path.push_back("v");
  path.push_back("v");
  log_entry(prm, path, "a");
  log_entry(prm, path, "b");
  path.pop_back();
  log_entry(prm, path, "a");
  path.pop_back();
  path.pop_back();
  path.push_back("z");
  log_entry(prm, path, "a");
  path.pop_back();

  prm.enter_subsection("y");
  deallog << "in y" << std::endl;
  log_entry(prm, path, "a");
  log_entry(prm, path, "b");

  prm.leave_subsection();
  prm.enter_subsection("v");
  deallog << "in v (first)" << std::endl;
  path.push_back("v");
  path.push_back("v");
  log_entry(prm, path, "a");
  log_entry(prm, path, "b");
  path.pop_back();
  log_entry(prm, path, "a");
  path.pop_back();

  prm.enter_subsection("v");
  deallog << "in v (second)" << std::endl;
  path.push_back("v");
  log_entry(prm, path, "a");
  log_entry(prm, path, "b");
  path.pop_back();
  log_entry(prm, path, "a");
}

void
test_getting_types()
{
  ParameterHandler prm;
  prm.enter_subsection("foo");
  prm.enter_subsection("bar");
  prm.declare_entry("jik", "1");
  prm.declare_entry("baz", "77.3");
  prm.declare_entry("yox", "true");
  prm.declare_entry("spam", "eggs");
  prm.leave_subsection();
  prm.leave_subsection();

  deallog << "in test_getting_types top level" << std::endl;

  if (prm.get_integer({"foo", "bar"}, "jik") != 1)
    deallog << "unexpected failure to read \"1\" as an integer" << std::endl;

  if (prm.get_double({"foo", "bar"}, "jik") != 1)
    deallog << "unexpected failure to read \"1\" as a double" << std::endl;

  if (prm.get_double({"foo", "bar"}, "baz") != 77.3)
    deallog << "unexpected failure to read \"77.3\" as a double" << std::endl;

  if (prm.get_bool({"foo", "bar"}, "yox") != true)
    deallog << "unexpected failure to read \"true\" as a boolean" << std::endl;

  try
    {
      prm.get_integer({"foo", "bar"}, "baz");
    }
  catch (const StandardExceptions::ExcMessage &e)
    {
      log_exception(e);
    }

  try
    {
      prm.get_integer({"foo", "bar"}, "yox");
    }
  catch (const StandardExceptions::ExcMessage &e)
    {
      log_exception(e);
    }

  try
    {
      prm.get_double({"foo", "bar"}, "yox");
    }
  catch (const StandardExceptions::ExcMessage &e)
    {
      log_exception(e);
    }

  try
    {
      prm.get_bool({"foo", "bar"}, "jik");
    }
  catch (const StandardExceptions::ExcMessage &e)
    {
      log_exception(e);
    }

  try
    {
      prm.get_bool({"foo", "bar"}, "spam");
    }
  catch (const StandardExceptions::ExcMessage &e)
    {
      log_exception(e);
    }
}

void
test_weird_strings()
{
  const std::vector<std::string> weird_strings = {
    ".,/<[\"';:=-_)*&~\t`/.\\",
    "$.7.%...  =   . -",
    "`/=/a/!",
    "<<>>>]]]\t\n\n   ",
    "****&//&%.^$!.$@$%@^*&(*)_/*-`~",
    ".,/<[\"';:=-_)*&~\t`/.\\",
    "value",
    " set value = 5 \n\n."};

  ParameterHandler prm;

  for (int i = 0; i < 6; ++i)
    {
      prm.enter_subsection(weird_strings[i]);
    }

  prm.declare_entry(weird_strings[6], weird_strings[7]);

  for (int i = 0; i < 6; ++i)
    {
      prm.leave_subsection();
    }

  std::vector<std::string> path = {};
  for (int i = 0; i < 3; ++i)
    {
      prm.enter_subsection(weird_strings[i]);
      path.push_back(weird_strings[3 + i]);
    }

  deallog << "in <weird_strings[3]>" << std::endl;
  log_entry(prm, path, weird_strings[6]);
}

int
main()
{
  initlog();

  try
    {
      test_basic();
      test_getting_types();
      test_weird_strings();
    }
  catch (const std::exception &exc)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;

      return 1;
    }
  catch (...)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };

  return 0;
}
