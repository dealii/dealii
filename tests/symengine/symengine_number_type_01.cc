// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check that the wrapper for symengine expressions works as expected

#include <deal.II/base/signaling_nan.h>

#include <deal.II/differentiation/sd.h>

#include <symengine/complex_double.h>
#include <symengine/expression.h>
#include <symengine/integer.h>
#include <symengine/logic.h>
#include <symengine/parser.h>
#include <symengine/rational.h>
#include <symengine/real_double.h>
#include <symengine/symbol.h>

#include <fstream>
#include <iomanip>

#include "../tests.h"

#include "../serialization/serialization.h"

namespace SD = Differentiation::SD;
namespace SE = ::SymEngine;


// Overload the comparison function in serialization.h
template <>
bool
compare(const SD::Expression &t1, const SD::Expression &t2)
{
  // Could also do the following if it ever becomes necessary:
  // (t1.get_value().__cmp__(t2.get_value()) == 0)
  return t1.get_expression() == t2.get_expression();
}


int
main()
{
  initlog();

  using SD_number_t = SD::Expression;

  deallog << "Constructors" << std::endl;
  {
    std::cout << "Constructor: Boolean" << std::endl;
    {
      const SD_number_t se_bool_true(true);
      const SD_number_t se_bool_false(false);
      Assert(SE::eq(se_bool_true.get_value(), *SE::boolTrue),
             ExcMessage("Problem with constructor"));
      Assert(SE::eq(se_bool_false.get_value(), *SE::boolFalse),
             ExcMessage("Problem with constructor"));
    }

    std::cout << "Constructor: Integer" << std::endl;
    {
      const int         x = -1;
      const SD_number_t se_number(x);
      Assert(SE::eq(se_number.get_value(), *SE::integer(x)),
             ExcMessage("Problem with constructor"));
    }

    std::cout << "Constructor: Unsigned integer" << std::endl;
    {
      const unsigned int x = 2;
      const SD_number_t  se_number(x);
      Assert(SE::eq(se_number.get_value(), *SE::integer(x)),
             ExcMessage("Problem with constructor"));
    }

    std::cout << "Constructor: Rational" << std::endl;
    {
      const unsigned int x = 2;
      const unsigned int y = 3;
      const SD_number_t  se_number(x, y);
      Assert(SE::eq(se_number.get_value(), *SE::rational(x, y)),
             ExcMessage("Problem with constructor"));
    }

    std::cout << "Constructor: Double" << std::endl;
    {
      const double      x = 2.5;
      const SD_number_t se_number(x);
      Assert(SE::eq(se_number.get_value(), *SE::real_double(x)),
             ExcMessage("Problem with constructor"));
    }

    std::cout << "Constructor: Complex double" << std::endl;
    {
      const std::complex<double> x(2.5, 1.5);
      const SD_number_t          se_number(x);
      Assert(SE::eq(se_number.get_value(), *SE::complex_double(x)),
             ExcMessage("Problem with constructor"));
    }

    std::cout << "Constructor: Float" << std::endl;
    {
      const float       x = 2.5;
      const SD_number_t se_number(x);
      Assert(SE::eq(se_number.get_value(), *SE::real_double(x)),
             ExcMessage("Problem with constructor"));
    }

    std::cout << "Constructor: Complex float" << std::endl;
    {
      const std::complex<float> x(2.5, 1.5);
      const SD_number_t         se_number(x);
      Assert(SE::eq(se_number.get_value(), *SE::complex_double(x)),
             ExcMessage("Problem with constructor"));
    }

    std::cout << "Constructor: Symbol" << std::endl;
    {
      const std::string f                   = "f";
      const bool        parse_as_expression = false;
      const SD_number_t se_number(f, parse_as_expression);
      Assert(SE::eq(se_number.get_value(), *SE::symbol(f)),
             ExcMessage("Problem with constructor"));
    }

    std::cout << "Constructor: Symbol" << std::endl;
    {
      const std::string f = "func";
      const SD_number_t se_number(f.c_str());
      Assert(SE::eq(se_number.get_value(), *SE::symbol(f)),
             ExcMessage("Problem with constructor"));
    }

    std::cout << "Constructor: Symbol" << std::endl;
    {
      const std::string f                   = "func";
      const bool        parse_as_expression = false;
      const SD_number_t se_number(f, parse_as_expression);
      Assert(SE::eq(se_number.get_value(), *SE::symbol(f)),
             ExcMessage("Problem with constructor"));
    }

    std::cout << "Constructor: Symbolic expression" << std::endl;
    {
      const std::string f                   = "x + y";
      const bool        parse_as_expression = true;
      const SD_number_t se_number(f, parse_as_expression);
      Assert(SE::eq(se_number.get_value(), *SE::parse(f)),
             ExcMessage("Problem with constructor"));
    }

    std::cout << "Copy constructor" << std::endl;
    {
      const double      x = 2.5;
      const SD_number_t se_number_1(x);
      const SD_number_t se_number_2(se_number_1);
      Assert(SE::eq(se_number_2.get_value(), *SE::real_double(x)),
             ExcMessage("Problem with constructor"));
    }
  }

  deallog << "Utility functions" << std::endl;
  {
    std::cout << "Parse" << std::endl;
    {
      const double      x = 2.5;
      SD_number_t       se_number(x);
      const std::string f = "x + y";
      se_number.parse(f);
      Assert(SE::eq(se_number.get_value(), *SE::parse(f)),
             ExcMessage("Problem with utility function: Parse"));
    }

    std::cout << "Save/load" << std::endl;
    {
      // Output
      const SD_number_t se_number_x_out("x");
      const SD_number_t se_number_y_out("y");
      const SD_number_t se_number_f_out("x + y + 2.5",
                                        true /*parse_as_expression*/);
      {
        std::ofstream o_str("save_load.txt");
        Assert(o_str.is_open(), ExcInternalError());
        se_number_x_out.save(o_str);
        se_number_y_out.save(o_str);
        se_number_f_out.save(o_str);
        o_str.close();
      }

      // Input
      SD_number_t se_number_x_in;
      SD_number_t se_number_y_in;
      SD_number_t se_number_f_in;
      {
        std::ifstream i_str("save_load.txt");
        Assert(i_str.is_open(), ExcInternalError());
        se_number_x_in.load(i_str);
        se_number_y_in.load(i_str);
        se_number_f_in.load(i_str);
      }

      Assert(static_cast<bool>(se_number_x_in == se_number_x_out),
             ExcMessage("Problem with utility function: Save/load"));
      Assert(static_cast<bool>(se_number_y_in == se_number_y_out),
             ExcMessage("Problem with utility function: Save/load"));
      Assert(static_cast<bool>(se_number_f_in == se_number_f_out),
             ExcMessage("Problem with utility function: Save/load"));
    }

    deallog.push("Serialization");
    std::cout << "Serialization" << std::endl;
    {
      const SD_number_t se_number_x_out("x");
      const SD_number_t se_number_y_out("y");
      const SD_number_t se_number_f_out("x + y + 2.5",
                                        true /*parse_as_expression*/);

      // From serialization.h
      SD_number_t se_number_x_in;
      SD_number_t se_number_y_in;
      SD_number_t se_number_f_in;
      verify(se_number_x_out, se_number_x_in);
      verify(se_number_y_out, se_number_y_in);
      verify(se_number_f_out, se_number_f_in);

      Assert(static_cast<bool>(se_number_x_in == se_number_x_out),
             ExcMessage("Problem with utility function: Serialization"));
      Assert(static_cast<bool>(se_number_y_in == se_number_y_out),
             ExcMessage("Problem with utility function: Serialization"));
      Assert(static_cast<bool>(se_number_f_in == se_number_f_out),
             ExcMessage("Problem with utility function: Serialization"));
    }
    deallog.pop();
  }

  deallog << "Comparison operators" << std::endl;
  {
    std::cout << "Equality (true)" << std::endl;
    {
      const int         x = 2;
      const int         y = 2;
      const SD_number_t se_number_1(x);
      const SD_number_t se_number_2(y);
      Assert(static_cast<bool>(se_number_1 == se_number_2) == true,
             ExcMessage("Problem with comparison operator"));
    }

    std::cout << "Equality (false)" << std::endl;
    {
      const int         x = 2;
      const int         y = 3;
      const SD_number_t se_number_1(x);
      const SD_number_t se_number_2(y);
      Assert(static_cast<bool>(se_number_1 == se_number_2) == false,
             ExcMessage("Problem with comparison operator"));
    }

    std::cout << "Non-equality (true)" << std::endl;
    {
      const int         x = 2;
      const int         y = 3;
      const SD_number_t se_number_1(x);
      const SD_number_t se_number_2(y);
      Assert(static_cast<bool>(se_number_1 != se_number_2) == true,
             ExcMessage("Problem with comparison operator"));
    }

    std::cout << "Non-equality (false)" << std::endl;
    {
      const int         x = 2;
      const int         y = 2;
      const SD_number_t se_number_1(x);
      const SD_number_t se_number_2(y);
      Assert(static_cast<bool>(se_number_1 != se_number_2) == false,
             ExcMessage("Problem with comparison operator"));
    }

    std::cout << "Less-than operator (true)" << std::endl;
    {
      const int         x = 2;
      const int         y = 3;
      const SD_number_t se_number_1(x);
      const SD_number_t se_number_2(y);
      Assert(static_cast<bool>(se_number_1 < se_number_2) == true,
             ExcMessage("Problem with comparison operator"));
    }

    std::cout << "Less-than operator (false)" << std::endl;
    {
      const int         x = 3;
      const int         y = 2;
      const int         z = 3;
      const SD_number_t se_number_1(x);
      const SD_number_t se_number_2(y);
      const SD_number_t se_number_3(z);
      Assert(static_cast<bool>(se_number_1 < se_number_2) == false,
             ExcMessage("Problem with comparison operator"));
      Assert(static_cast<bool>(se_number_1 < se_number_3) == false,
             ExcMessage("Problem with comparison operator"));
    }

    std::cout << "Less-than-or-equals operator (true)" << std::endl;
    {
      const int         x = 2;
      const int         y = 3;
      const int         z = 2;
      const SD_number_t se_number_1(x);
      const SD_number_t se_number_2(y);
      const SD_number_t se_number_3(z);
      Assert(static_cast<bool>(se_number_1 <= se_number_2) == true,
             ExcMessage("Problem with comparison operator"));
      Assert(static_cast<bool>(se_number_1 <= se_number_3) == true,
             ExcMessage("Problem with comparison operator"));
    }

    std::cout << "Less-than-or-equals operator (false)" << std::endl;
    {
      const int         x = 3;
      const int         y = 2;
      const SD_number_t se_number_1(x);
      const SD_number_t se_number_2(y);
      Assert(static_cast<bool>(se_number_1 <= se_number_2) == false,
             ExcMessage("Problem with comparison operator"));
    }

    std::cout << "Greater-than operator (true)" << std::endl;
    {
      const int         x = 3;
      const int         y = 2;
      const SD_number_t se_number_1(x);
      const SD_number_t se_number_2(y);
      Assert(static_cast<bool>(se_number_1 > se_number_2) == true,
             ExcMessage("Problem with comparison operator"));
    }

    std::cout << "Greater-than operator (false)" << std::endl;
    {
      const int         x = 2;
      const int         y = 3;
      const int         z = 2;
      const SD_number_t se_number_1(x);
      const SD_number_t se_number_2(y);
      const SD_number_t se_number_3(z);
      Assert(static_cast<bool>(se_number_1 > se_number_2) == false,
             ExcMessage("Problem with comparison operator"));
      Assert(static_cast<bool>(se_number_1 > se_number_3) == false,
             ExcMessage("Problem with comparison operator"));
    }

    std::cout << "Greater-than-or-equals operator (true)" << std::endl;
    {
      const int         x = 3;
      const int         y = 2;
      const int         z = 3;
      const SD_number_t se_number_1(x);
      const SD_number_t se_number_2(y);
      const SD_number_t se_number_3(z);
      Assert(static_cast<bool>(se_number_1 >= se_number_2) == true,
             ExcMessage("Problem with comparison operator"));
      Assert(static_cast<bool>(se_number_1 >= se_number_3) == true,
             ExcMessage("Problem with comparison operator"));
    }

    std::cout << "Greater-than-or-equals operator (false)" << std::endl;
    {
      const int         x = 2;
      const int         y = 3;
      const SD_number_t se_number_1(x);
      const SD_number_t se_number_2(y);
      Assert(static_cast<bool>(se_number_1 >= se_number_2) == false,
             ExcMessage("Problem with comparison operator"));
    }
  }

  deallog << "Logical operators" << std::endl;
  {
    std::cout << "Logical not" << std::endl;
    {
      const bool        x = true;
      const SD_number_t se_number(x);
      Assert(static_cast<bool>(!se_number) == (!x),
             ExcMessage("Problem with logical operator"));
      Assert(static_cast<bool>(!!se_number) == (!!x),
             ExcMessage("Problem with logical operator"));
    }

    std::cout << "Logical and" << std::endl;
    {
      auto test_logical_op = [](const bool x, const bool y) {
        const SD_number_t se_number_1(x);
        const SD_number_t se_number_2(y);
        Assert(static_cast<bool>(se_number_1 & se_number_2) == (x & y),
               ExcMessage("Problem with logical operator"));
      };

      test_logical_op(true, true);
      test_logical_op(true, false);
      test_logical_op(false, true);
      test_logical_op(false, false);
    }

    std::cout << "Logical inclusive or" << std::endl;
    {
      auto test_logical_op = [](const bool x, const bool y) {
        const SD_number_t se_number_1(x);
        const SD_number_t se_number_2(y);
        Assert(static_cast<bool>(se_number_1 | se_number_2) == (x | y),
               ExcMessage("Problem with logical operator"));
      };

      test_logical_op(true, true);
      test_logical_op(true, false);
      test_logical_op(false, true);
      test_logical_op(false, false);
    }

    std::cout << "Logical xor" << std::endl;
    {
      auto test_logical_op = [](const bool x, const bool y) {
        const SD_number_t se_number_1(x);
        const SD_number_t se_number_2(y);
        Assert(static_cast<bool>(se_number_1 ^ se_number_2) == (x ^ y),
               ExcMessage("Problem with logical operator"));
      };

      test_logical_op(true, true);
      test_logical_op(true, false);
      test_logical_op(false, true);
      test_logical_op(false, false);
    }

    std::cout << "And" << std::endl;
    {
      auto test_logical_op = [](const bool x, const bool y) {
        const SD_number_t se_number_1(x);
        const SD_number_t se_number_2(y);
        Assert(static_cast<bool>(se_number_1 && se_number_2) == (x && y),
               ExcMessage("Problem with logical operator"));
      };

      test_logical_op(true, true);
      test_logical_op(true, false);
      test_logical_op(false, true);
      test_logical_op(false, false);
    }

    std::cout << "Inclusive or" << std::endl;
    {
      auto test_logical_op = [](const bool x, const bool y) {
        const SD_number_t se_number_1(x);
        const SD_number_t se_number_2(y);
        Assert(static_cast<bool>(se_number_1 || se_number_2) == (x || y),
               ExcMessage("Problem with logical operator"));
      };

      test_logical_op(true, true);
      test_logical_op(true, false);
      test_logical_op(false, true);
      test_logical_op(false, false);
    }
  }

  deallog << "Assignment operators" << std::endl;
  {
    std::cout << "Assignment operator (SD_number_t)" << std::endl;
    {
      const int         x = 1;
      const int         y = 2;
      SD_number_t       se_number_1(x);
      const SD_number_t se_number_2(y);
      se_number_1 = se_number_2;
      Assert(SE::eq(se_number_1.get_value(), *SE::integer(y)),
             ExcMessage("Problem with assignment operator"));
    }

    std::cout << "Assignment operator (arithmetic type)" << std::endl;
    {
      const int    x = 1;
      const double y = 2;
      SD_number_t  se_number(x);
      se_number = y;
      Assert(SE::eq(se_number.get_value(), *SE::real_double(y)),
             ExcMessage("Problem with assignment operator"));
    }

    std::cout << "Addition assignment (SD_number_t)" << std::endl;
    {
      const int         x = 1;
      const int         y = 2;
      SD_number_t       se_number_1(x);
      const SD_number_t se_number_2(y);
      se_number_1 += se_number_2;
      Assert(SE::eq(se_number_1.get_value(), *SE::integer(x + y)),
             ExcMessage("Problem with assignment operator"));
    }

    std::cout << "Addition assignment (arithmetic type; same type)"
              << std::endl;
    {
      const int   x = 1;
      const int   y = 2;
      SD_number_t se_number(x);
      se_number += y;
      Assert(SE::eq(se_number.get_value(), *SE::integer(x + y)),
             ExcMessage("Problem with assignment operator"));
    }

    std::cout << "Addition assignment (arithmetic type; different type)"
              << std::endl;
    {
      const int    x = 1;
      const double y = 2;
      SD_number_t  se_number(x);
      se_number += y;
      Assert(SE::eq(se_number.get_value(), *SE::real_double(x + y)),
             ExcMessage("Problem with assignment operator"));
    }

    std::cout << "Subtraction assignment (SD_number_t)" << std::endl;
    {
      const int         x = 1;
      const int         y = 2;
      SD_number_t       se_number_1(x);
      const SD_number_t se_number_2(y);
      se_number_1 -= se_number_2;
      Assert(SE::eq(se_number_1.get_value(), *SE::integer(x - y)),
             ExcMessage("Problem with assignment operator"));
    }

    std::cout << "Subtraction assignment (arithmetic type; same type)"
              << std::endl;
    {
      const int   x = 1;
      const int   y = 2;
      SD_number_t se_number(x);
      se_number -= y;
      Assert(SE::eq(se_number.get_value(), *SE::integer(x - y)),
             ExcMessage("Problem with assignment operator"));
    }

    std::cout << "Subtraction assignment (arithmetic type; different type)"
              << std::endl;
    {
      const int    x = 1;
      const double y = 2;
      SD_number_t  se_number(x);
      se_number -= y;
      Assert(SE::eq(se_number.get_value(), *SE::real_double(x - y)),
             ExcMessage("Problem with assignment operator"));
    }

    std::cout << "Multiplication assignment (SD_number_t)" << std::endl;
    {
      const int         x = 2;
      const int         y = 3;
      SD_number_t       se_number_1(x);
      const SD_number_t se_number_2(y);
      se_number_1 *= se_number_2;
      Assert(SE::eq(se_number_1.get_value(), *SE::integer(x * y)),
             ExcMessage("Problem with assignment operator"));
    }

    std::cout << "Multiplication assignment (arithmetic type; same type)"
              << std::endl;
    {
      const int   x = 2;
      const int   y = 3;
      SD_number_t se_number(x);
      se_number *= y;
      Assert(SE::eq(se_number.get_value(), *SE::integer(x * y)),
             ExcMessage("Problem with assignment operator"));
    }

    std::cout << "Multiplication assignment (arithmetic type; different type)"
              << std::endl;
    {
      const int    x = 2;
      const double y = 3;
      SD_number_t  se_number(x);
      se_number *= y;
      Assert(SE::eq(se_number.get_value(), *SE::real_double(x * y)),
             ExcMessage("Problem with assignment operator"));
    }

    std::cout << "Division assignment (SD_number_t)" << std::endl;
    {
      const int         x = 8;
      const int         y = 2;
      SD_number_t       se_number_1(x);
      const SD_number_t se_number_2(y);
      se_number_1 /= se_number_2;
      Assert(SE::eq(se_number_1.get_value(), *SE::integer(x / y)),
             ExcMessage("Problem with assignment operator"));
    }

    std::cout << "Division assignment (arithmetic type; same type)"
              << std::endl;
    {
      const int   x = 8;
      const int   y = 2;
      SD_number_t se_number(x);
      se_number /= y;
      Assert(SE::eq(se_number.get_value(), *SE::integer(x / y)),
             ExcMessage("Problem with assignment operator"));
    }

    std::cout << "Division assignment (arithmetic type; different type)"
              << std::endl;
    {
      const int    x = 8;
      const double y = 2;
      SD_number_t  se_number(x);
      se_number /= y;
      Assert(SE::eq(se_number.get_value(), *SE::real_double(x / y)),
             ExcMessage("Problem with assignment operator"));
    }
  }

  deallog << "Math operators" << std::endl;
  {
    std::cout << "Addition" << std::endl;
    {
      const int         x = 1;
      const int         y = 2;
      const SD_number_t se_number_1(x);
      const SD_number_t se_number_2(y);
      const SD_number_t se_number_3 = se_number_1 + se_number_2;
      Assert(SE::eq(se_number_3.get_value(), *SE::integer(x + y)),
             ExcMessage("Problem with math operator"));
    }

    std::cout << "Subtraction" << std::endl;
    {
      const int         x = 1;
      const int         y = 2;
      const SD_number_t se_number_1(x);
      const SD_number_t se_number_2(y);
      const SD_number_t se_number_3 = se_number_1 - se_number_2;
      Assert(SE::eq(se_number_3.get_value(), *SE::integer(x - y)),
             ExcMessage("Problem with math operator"));
    }

    std::cout << "Multiplication" << std::endl;
    {
      const int         x = 2;
      const int         y = 3;
      const SD_number_t se_number_1(x);
      const SD_number_t se_number_2(y);
      const SD_number_t se_number_3 = se_number_1 * se_number_2;
      Assert(SE::eq(se_number_3.get_value(), *SE::integer(x * y)),
             ExcMessage("Problem with math operator"));
    }

    std::cout << "Division" << std::endl;
    {
      const int         x = 8;
      const int         y = 2;
      const SD_number_t se_number_1(x);
      const SD_number_t se_number_2(y);
      const SD_number_t se_number_3 = se_number_1 / se_number_2;
      Assert(SE::eq(se_number_3.get_value(), *SE::integer(x / y)),
             ExcMessage("Problem with math operator"));
    }
  }

  deallog << "Differentiation" << std::endl;
  {
    std::cout << "Basic differentiation (SD_number_t)" << std::endl;
    {
      const SD_number_t x("x");
      SD_number_t       f(1.0);
      f *= x;
      f *= x; // f = x^2
      const SD_number_t df_dx = f.differentiate(x);
      Assert(static_cast<bool>(df_dx == (2.0 * x)),
             ExcMessage("Problem with differentiation"));
    }

    std::cout << "Basic differentiation (Symbol)" << std::endl;
    {
      const auto        symb_x = SE::symbol("x");
      const SD_number_t x(symb_x);
      SD_number_t       f(1.0);
      f *= x;
      f *= x; // f = x^2
      const SD_number_t df_dx = f.differentiate(symb_x);
      Assert(static_cast<bool>(df_dx == SD_number_t(2.0 * x)),
             ExcMessage("Problem with differentiation"));
    }

    std::cout << "Basic differentiation (Basic representing a symbolic type)"
              << std::endl;
    {
      const auto symb_x_plus_y = SE::add(SE::symbol("x"), SE::symbol("y"));
      const SD_number_t x_plus_y(symb_x_plus_y);
      SD_number_t       f(1.0);
      f *= x_plus_y;
      f *= x_plus_y; // f = (x+y)^2
      const SD_number_t df_dxpy = f.differentiate(symb_x_plus_y);
      // For some reason, the next check is considered to be "symbolic" and
      // cannot be evaluated.
      //                          df_dxpy = 2*(x + y)*1.0
      // 2.0 * SD_number_t(symb_x_plus_y) = 2.0*(x + y)
      // Assert(static_cast<bool>(df_dxpy == (2.0 *
      // SD_number_t(symb_x_plus_y))),
      //        ExcMessage("Problem with differentiation"));
    }
  }

  deallog << "Substitution" << std::endl;
  {
    std::cout << "Basic substitution" << std::endl;
    {
      const SD_number_t x("x");
      const SD_number_t y("y");
      const SD_number_t x_plus_y = x + y;
      const SD_number_t f        = x_plus_y * x_plus_y;
      Assert(static_cast<bool>(f == SD_number_t("(x + y)**2", true)),
             ExcMessage("Problem with substitution"));

      const double      x_val   = 1.0;
      const SD_number_t f_sub_x = f.substitute(x, SD_number_t(x_val));
      Assert(static_cast<bool>(f_sub_x == SD_number_t("(1.0 + y)**2", true)),
             ExcMessage("Problem with substitution"));

      const double      y_val    = 2.0;
      const SD_number_t f_sub_xy = f_sub_x.substitute(y, SD_number_t(y_val));
      Assert(static_cast<double>(f_sub_xy) ==
               (x_val * x_val + 2.0 * x_val * y_val + y_val * y_val),
             ExcMessage("Problem with substitution"));
    }

    std::cout << "Templated-deduced substitution" << std::endl;
    {
      const SD_number_t x("x");
      const SD_number_t y("y");
      const SD_number_t x_plus_y = x + y;
      const SD_number_t f        = x_plus_y * x_plus_y;
      Assert(static_cast<bool>(f == SD_number_t("(x + y)**2", true)),
             ExcMessage("Problem with substitution"));

      const double      x_val   = 1.0;
      const SD_number_t f_sub_x = f.substitute(x, x_val);
      Assert(static_cast<bool>(f_sub_x == SD_number_t("(1.0 + y)**2", true)),
             ExcMessage("Problem with substitution"));

      const double      y_val    = 2.0;
      const SD_number_t f_sub_xy = f_sub_x.substitute(y, y_val);
      Assert(static_cast<double>(f_sub_xy) ==
               (x_val * x_val + 2.0 * x_val * y_val + y_val * y_val),
             ExcMessage("Problem with substitution"));
    }

    std::cout << "Map substitution" << std::endl;
    {
      const SD_number_t x("x");
      const SD_number_t y("y");
      const SD_number_t x_plus_y = x + y;
      const SD_number_t f        = x_plus_y * x_plus_y;

      const double                x_val = 1.0;
      const double                y_val = 2.0;
      SD::types::substitution_map sub_map;
      sub_map[x]                 = SD_number_t(x_val);
      sub_map[y]                 = SD_number_t(y_val);
      const SD_number_t f_sub_xy = f.substitute(sub_map);
      Assert(static_cast<double>(f_sub_xy) ==
               (x_val * x_val + 2.0 * x_val * y_val + y_val * y_val),
             ExcMessage("Problem with substitution"));
    }
  }

  deallog << "Evaluation" << std::endl;
  {
    std::cout << "Substitution with evaluation (integer)" << std::endl;
    {
      const SD_number_t x("x");
      const SD_number_t y("y");
      const SD_number_t x_plus_y = x + y;
      const SD_number_t f        = x_plus_y * x_plus_y;

      const int                   x_val = 1;
      const int                   y_val = 2;
      SD::types::substitution_map sub_map;
      sub_map[x]      = SD_number_t(x_val);
      sub_map[y]      = SD_number_t(y_val);
      const int f_val = f.substitute_and_evaluate<int>(sub_map);
      Assert(f_val == (x_val * x_val + 2 * x_val * y_val + y_val * y_val),
             ExcMessage("Problem with evaluation"));
    }

    std::cout << "Substitution with evaluation (double)" << std::endl;
    {
      const SD_number_t x("x");
      const SD_number_t y("y");
      const SD_number_t x_plus_y = x + y;
      const SD_number_t f        = x_plus_y * x_plus_y;

      const double                x_val = 1.0;
      const double                y_val = 2.0;
      SD::types::substitution_map sub_map;
      sub_map[x]         = SD_number_t(x_val);
      sub_map[y]         = SD_number_t(y_val);
      const double f_val = f.substitute_and_evaluate<double>(sub_map);
      Assert(f_val == (x_val * x_val + 2.0 * x_val * y_val + y_val * y_val),
             ExcMessage("Problem with evaluation"));
    }

    std::cout << "Substitution with evaluation (float)" << std::endl;
    {
      const SD_number_t x("x");
      const SD_number_t y("y");
      const SD_number_t x_plus_y = x + y;
      const SD_number_t f        = x_plus_y * x_plus_y;

      const float                 x_val = 1.0;
      const float                 y_val = 2.0;
      SD::types::substitution_map sub_map;
      sub_map[x]        = SD_number_t(x_val);
      sub_map[y]        = SD_number_t(y_val);
      const float f_val = f.substitute_and_evaluate<float>(sub_map);
      Assert(f_val == (x_val * x_val + 2.0f * x_val * y_val + y_val * y_val),
             ExcMessage("Problem with evaluation"));
    }

    std::cout << "Substitution with evaluation (complex double)" << std::endl;
    {
      const SD_number_t x("x");
      const SD_number_t y("y");
      const SD_number_t x_plus_y = x + y;
      const SD_number_t f        = x_plus_y * x_plus_y;

      const std::complex<double>  x_val(1.0, 3.0);
      const std::complex<double>  y_val(2.0, -2.0);
      SD::types::substitution_map sub_map;
      sub_map[x] = SD_number_t(x_val);
      sub_map[y] = SD_number_t(y_val);
      const std::complex<double> f_val =
        f.substitute_and_evaluate<std::complex<double>>(sub_map);
      Assert(std::abs(f_val - (x_val * x_val + 2.0 * x_val * y_val +
                               y_val * y_val)) < 1e-12,
             ExcMessage("Problem with evaluation"));
    }

    std::cout << "Substitution with evaluation (complex float)" << std::endl;
    {
      const SD_number_t x("x");
      const SD_number_t y("y");
      const SD_number_t x_plus_y = x + y;
      const SD_number_t f        = x_plus_y * x_plus_y;

      const std::complex<float>   x_val(1.0, 3.0);
      const std::complex<float>   y_val(2.0, -2.0);
      SD::types::substitution_map sub_map;
      sub_map[x] = SD_number_t(x_val);
      sub_map[y] = SD_number_t(y_val);
      const std::complex<float> f_val =
        f.substitute_and_evaluate<std::complex<float>>(sub_map);
      Assert(std::abs(f_val - (x_val * x_val + 2.0f * x_val * y_val +
                               y_val * y_val)) < 1e-12,
             ExcMessage("Problem with evaluation"));
    }

    std::cout << "Substitution with evaluation (mixed arithmetic types)"
              << std::endl;
    {
      const SD_number_t x("x");
      const SD_number_t y("y");
      const SD_number_t x_plus_y = x + y;
      const SD_number_t f        = x_plus_y * x_plus_y;

      const int                   x_val = 1;
      const int                   y_val = 2;
      SD::types::substitution_map sub_map;
      sub_map[x]         = SD_number_t(x_val);
      sub_map[y]         = SD_number_t(y_val);
      const double f_val = f.substitute_and_evaluate<double>(sub_map);
      Assert(f_val == static_cast<double>(x_val * x_val + 2.0f * x_val * y_val +
                                          y_val * y_val),
             ExcMessage("Problem with evaluation"));
    }

    std::cout << "Substitution with evaluation (mixed complex types)"
              << std::endl;
    {
      const SD_number_t x("x");
      const SD_number_t y("y");
      const SD_number_t x_plus_y = x + y;
      const SD_number_t f        = x_plus_y * x_plus_y;

      const std::complex<float>   x_val(1.0, 3.0);
      const std::complex<float>   y_val(2.0, -2.0);
      SD::types::substitution_map sub_map;
      sub_map[x] = SD_number_t(x_val);
      sub_map[y] = SD_number_t(y_val);
      const std::complex<double> f_val =
        f.substitute_and_evaluate<std::complex<double>>(sub_map);
      Assert(std::abs(f_val - static_cast<std::complex<double>>(
                                x_val * x_val + 2.0f * x_val * y_val +
                                y_val * y_val)) < 1e-12,
             ExcMessage("Problem with evaluation"));
    }

    std::cout << "Substitution with evaluation (single piecewise function)"
              << std::endl;
    {
      const SD_number_t x("x");
      const SD_number_t y("y");
      const SD_number_t x_plus_y  = x + y;
      const SD_number_t x_minus_y = x - y;
      const SD_number_t f_plus    = x_plus_y * x_plus_y;
      const SD_number_t f_minus   = x_minus_y * x_minus_y;
      const SD_number_t f((x > SD_number_t(0.0)), f_plus, f_minus);

      // Condition is met
      {
        const double                x_val = 1.0;
        const double                y_val = 2.0;
        SD::types::substitution_map sub_map;
        sub_map[x]               = SD_number_t(x_val);
        sub_map[y]               = SD_number_t(y_val);
        const SD_number_t f_subs = f.substitute(sub_map);
        const double      f_val  = static_cast<double>(f_subs);
        Assert(f_val == (x_val * x_val + 2.0 * x_val * y_val + y_val * y_val),
               ExcMessage(
                 "Problem with evaluation of single piecewise function"));
      }

      // Condition is not met
      {
        const double                x_val = -1.0;
        const double                y_val = 2.0;
        SD::types::substitution_map sub_map;
        sub_map[x]         = SD_number_t(x_val);
        sub_map[y]         = SD_number_t(y_val);
        const double f_val = f.substitute_and_evaluate<double>(sub_map);
        Assert(f_val == (x_val * x_val - 2.0 * x_val * y_val + y_val * y_val),
               ExcMessage(
                 "Problem with evaluation of single piecewise function"));
      }
    }

    std::cout << "Substitution with evaluation (double piecewise function)"
              << std::endl;
    {
      const SD_number_t x("x");
      const SD_number_t y("y");
      const SD_number_t x_plus_y    = x + y;
      const SD_number_t x_minus_y   = x - y;
      const SD_number_t f_plus      = x_plus_y * x_plus_y;
      const SD_number_t f_plus_plus = 2.0 * x_plus_y * x_plus_y;
      const SD_number_t f_minus     = x_minus_y * x_minus_y;
      const SD_number_t f((x > SD_number_t(0.0)),
                          SD_number_t(x > SD_number_t(2.0),
                                      f_plus_plus,
                                      f_plus),
                          f_minus);

      // Outer condition is met, inner condition is met
      {
        const double                x_val = 3.0;
        const double                y_val = 2.0;
        SD::types::substitution_map sub_map;
        sub_map[x]               = SD_number_t(x_val);
        sub_map[y]               = SD_number_t(y_val);
        const SD_number_t f_subs = f.substitute(sub_map);
        const double      f_val  = static_cast<double>(f_subs);
        Assert(f_val == (2.0 * x_val * x_val + 4.0 * x_val * y_val +
                         2.0 * y_val * y_val),
               ExcMessage(
                 "Problem with evaluation of double piecewise function"));
      }

      // Outer condition is met, inner condition is not met
      {
        const double                x_val = 1.0;
        const double                y_val = 2.0;
        SD::types::substitution_map sub_map;
        sub_map[x]               = SD_number_t(x_val);
        sub_map[y]               = SD_number_t(y_val);
        const SD_number_t f_subs = f.substitute(sub_map);
        const double      f_val  = static_cast<double>(f_subs);
        Assert(f_val == (x_val * x_val + 2.0 * x_val * y_val + y_val * y_val),
               ExcMessage(
                 "Problem with evaluation of double piecewise function"));
      }

      // Condition is not met
      {
        const double                x_val = -1.0;
        const double                y_val = 2.0;
        SD::types::substitution_map sub_map;
        sub_map[x]         = SD_number_t(x_val);
        sub_map[y]         = SD_number_t(y_val);
        const double f_val = f.substitute_and_evaluate<double>(sub_map);
        Assert(f_val == (x_val * x_val - 2.0 * x_val * y_val + y_val * y_val),
               ExcMessage(
                 "Problem with evaluation of double piecewise function"));
      }
    }

    std::cout << "Substitution with evaluation (if-else_if-else type function)"
              << std::endl;
    {
      const SD_number_t x("x");
      const SD_number_t f_1 = x;
      const SD_number_t f_2 = 2 * x;
      const SD_number_t f_3 = 3 * x;
      const SD_number_t f_4 = 4 * x;
      const SD_number_t f({{(x > SD_number_t(4.0)), f_1},
                           {(x > SD_number_t(2.0)), f_2},
                           {(x > SD_number_t(0.0)), f_3}},
                          f_4);

      // Condition "if" is met
      {
        const double                x_val = 5.0;
        SD::types::substitution_map sub_map;
        sub_map[x]               = SD_number_t(x_val);
        const SD_number_t f_subs = f.substitute(sub_map);
        const double      f_val  = static_cast<double>(f_subs);
        Assert(
          f_val == x_val,
          ExcMessage(
            "Problem with evaluation of 'if' branch of terminated piecewise function"));
      }

      // Condition "else if (1)" is met
      {
        const double                x_val = 3.0;
        SD::types::substitution_map sub_map;
        sub_map[x]               = SD_number_t(x_val);
        const SD_number_t f_subs = f.substitute(sub_map);
        const double      f_val  = static_cast<double>(f_subs);
        Assert(
          f_val == (2 * x_val),
          ExcMessage(
            "Problem with evaluation of 'else if (1)' branch of terminated piecewise function"));
      }

      // Condition "else if (2)" is met
      {
        const double                x_val = 1.0;
        SD::types::substitution_map sub_map;
        sub_map[x]               = SD_number_t(x_val);
        const SD_number_t f_subs = f.substitute(sub_map);
        const double      f_val  = static_cast<double>(f_subs);
        Assert(
          f_val == (3 * x_val),
          ExcMessage(
            "Problem with evaluation of 'else if (2)' branch of terminated piecewise function"));
      }

      // Condition "else" is met
      {
        const double                x_val = -1.0;
        SD::types::substitution_map sub_map;
        sub_map[x]               = SD_number_t(x_val);
        const SD_number_t f_subs = f.substitute(sub_map);
        const double      f_val  = static_cast<double>(f_subs);
        Assert(
          f_val == (4 * x_val),
          ExcMessage(
            "Problem with evaluation of 'else' branch of terminated piecewise function"));
      }
    }

    std::cout << "Substitution with evaluation (if-else_if type function)"
              << std::endl;
    {
      const SD_number_t x("x");
      const SD_number_t f_1 = x;
      const SD_number_t f_2 = 2 * x;
      const SD_number_t f_3 = 3 * x;
      const SD_number_t f({{(x > SD_number_t(4.0)), f_1},
                           {(x > SD_number_t(2.0)), f_2},
                           {(x > SD_number_t(0.0)), f_3}});

      // Condition "if" is met
      {
        const double                x_val = 5.0;
        SD::types::substitution_map sub_map;
        sub_map[x]               = SD_number_t(x_val);
        const SD_number_t f_subs = f.substitute(sub_map);
        const double      f_val  = static_cast<double>(f_subs);
        Assert(
          f_val == x_val,
          ExcMessage(
            "Problem with evaluation of 'if' branch of non-terminated piecewise function"));
      }

      // Condition "else if (1)" is met
      {
        const double                x_val = 3.0;
        SD::types::substitution_map sub_map;
        sub_map[x]               = SD_number_t(x_val);
        const SD_number_t f_subs = f.substitute(sub_map);
        const double      f_val  = static_cast<double>(f_subs);
        Assert(
          f_val == (2 * x_val),
          ExcMessage(
            "Problem with evaluation of 'else if (1)' branch of non-terminated piecewise function"));
      }

      // Condition "else if (2)" is met
      {
        const double                x_val = 1.0;
        SD::types::substitution_map sub_map;
        sub_map[x]               = SD_number_t(x_val);
        const SD_number_t f_subs = f.substitute(sub_map);
        const double      f_val  = static_cast<double>(f_subs);
        Assert(
          f_val == (3 * x_val),
          ExcMessage(
            "Problem with evaluation of 'else if (2)' branch of non-terminated piecewise function"));
      }

      // No condition is met
      {
#if defined(DEBUG) && defined(DEAL_II_HAVE_FP_EXCEPTIONS)
        fedisableexcept(FE_INVALID);
#endif

        const double                x_val = -1.0;
        SD::types::substitution_map sub_map;
        sub_map[x]               = SD_number_t(x_val);
        const SD_number_t f_subs = f.substitute(sub_map);
        const double      f_val  = static_cast<double>(f_subs);
        Assert(
          numbers::is_nan(f_val) == true,
          ExcMessage(
            "Problem with evaluation of 'else' branch of non-terminated piecewise function"));
      }
    }
  }

  deallog << "OK" << std::endl;
}
