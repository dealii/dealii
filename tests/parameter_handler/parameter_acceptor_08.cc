// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
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

#include "../tests.h"

// Test subsectioning

class Test : public ParameterAcceptor
{
public:
  Test(const std::string &sec_name  = "First Class",
       const std::string &par_name  = "Parameter name",
       const std::string &par_value = "Parameter value")
    : ParameterAcceptor(sec_name)
    , par_name(par_name)
    , par_value(par_value)
  {
    add_parameter(par_name, this->par_value);

    deallog << "Section Name    : " << sec_name << std::endl;
    deallog << "N sections      : " << get_section_path().size() << std::endl;
    deallog << "Sections        : ";
    for (auto s : get_section_path())
      deallog << "\"" << s << "\"    ";
    deallog << std::endl << std::endl;
  };

private:
  std::string par_name;
  std::string par_value;
};

int
main()
{
  initlog();
  auto &prm = ParameterAcceptor::prm;

  // Relative
  Test a("Class A", "Parameter A", "a");

  // Absolute (== relative if no other path is added)
  Test c("/Class B", "Parameter C", "c");


  // Absolute (== relative if no other path is added)
  Test c2("/Class C/", "Parameter C", "c");

  // Absolute and no name
  Test a0("/", "Absolute parameter", "a");

  // Relative to absolute = absolute
  Test d("Class C", "Parameter C", "d");

  // Explicit, with absolute path
  Test d1("/Class D/Class D1", "Parameter D1", "d");

  // Relative to previous, with trailing
  Test d2("Class D2/", "Parameter D2", "d");

  // Relative to previous
  Test d3("Class D2D3", "Parameter D2D3", "e");

  // Another relative to reset
  Test d4("/Class E", "Parameter E", "e");

  // Same level of E1
  Test d5("Class E2/", "Parameter E", "e");

  // One level below E2
  Test d6("Class E2E3", "Parameter E", "e");

  ParameterAcceptor::declare_all_parameters();
  prm.log_parameters(deallog);
}
