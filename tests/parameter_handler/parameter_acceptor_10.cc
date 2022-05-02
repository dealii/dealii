//-----------------------------------------------------------
//
//    Copyright (C) 2017 - 2021 by the deal.II authors
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


// Test subsectioning within parameter acceptor itself

#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/path_search.h>

#include "../tests.h"

class Test : public ParameterAcceptor
{
public:
  Test(const std::string &sec_name        = "First Class",
       const std::string &par_name        = "Parameter name",
       const std::string &par_value       = "Parameter value",
       const std::string &subsection_name = "")
    : ParameterAcceptor(sec_name)
    , par_name(par_name)
    , par_value(par_value)
  {
    if (subsection_name != "")
      enter_subsection(subsection_name);
    add_parameter(par_name, this->par_value);
    if (subsection_name != "")
      leave_subsection();

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
  std::string subsection_name;
};

int
main()
{
  initlog();
  auto &prm = ParameterAcceptor::prm;

  // Relative, no subsections
  Test a("Class A", "Parameter A", "a");

  // Relative, with subsection
  Test c("Class B", "Parameter C", "c", "Subsection B");

  // Explicit, with absolute path
  Test d1("/Class D/Class D1", "Parameter D1", "d");

  // Relative to previous, with trailing
  Test d2("Class D2/", "Parameter D2", "d");

  // Relative to previous
  Test d3("Class D2D3", "Parameter D2D3", "e");

  // Relative to previous, with subsection
  Test d4("Class D2D3", "Parameter D2D3", "e", "Subsection D2D3");

  ParameterAcceptor::declare_all_parameters();
  prm.log_parameters(deallog);
}
