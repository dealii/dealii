//-----------------------------------------------------------
//
//    Copyright (C) 2017 by the deal.II authors
//
//    This file is part of the deal.II library.
//
//    The deal.II library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE at
//    the top level of the deal.II distribution.
//
//-----------------------------------------------------------



#include "../tests.h"
#include <deal.II/base/path_search.h>
#include <deal.II/base/parameter_acceptor.h>

// Test subsectioning

class Test : public ParameterAcceptor
{
public:
  Test(const std::string &sec_name  = "First Class",
       const std::string &par_name  = "Parameter name",
       const std::string &par_value = "Parameter value"):
    ParameterAcceptor(sec_name),
    par_name(par_name),
    par_value(par_value)
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
main ()
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
