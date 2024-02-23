// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
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

#include "../tests.h"

class FirstClass : public ParameterAcceptor
{
public:
  FirstClass(const std::string &name = "First Class")
    : ParameterAcceptor(name)
  {
    add_parameter("First int", f_i);
    add_parameter("First double", f_d);
    add_parameter("First bool", f_b);
    add_parameter("First string", f_s);
  };

private:
  int         f_i = 3;
  double      f_d = 7.7;
  bool        f_b = true;
  std::string f_s = "hello";
};

class SecondClass : public ParameterAcceptor
{
public:
  SecondClass(const std::string &name = "Second Class")
    : ParameterAcceptor(name)
  {
    add_parameter("Second int", s_i);
    add_parameter("Second double", s_d);
    add_parameter("Second bool", s_b);
    add_parameter("Second string", s_s);
  };

private:
  int         s_i = 5;
  double      s_d = 9.9;
  bool        s_b = false;
  std::string s_s = "bye bye";
};

int
main()
{
  initlog();

  FirstClass  f;
  SecondClass s;
  std::string output_name = "used_parameter_acceptor_06.tex";
  ParameterAcceptor::initialize(
    SOURCE_DIR "/parameter_acceptor_parameters/parameter_acceptor_05.prm",
    output_name);
  ParameterAcceptor::prm.log_parameters(deallog);
  std::ifstream file(output_name);

  std::string str;
  deallog << "reading " << output_name << std::endl;
  while (std::getline(file, str))
    deallog << str << std::endl;
}
