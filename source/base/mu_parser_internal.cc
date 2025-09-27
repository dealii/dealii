// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/mu_parser_internal.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/utilities.h>

#include <cmath>
#include <ctime>
#include <limits>
#include <map>
#include <mutex>
#include <random>
#include <vector>

#ifdef DEAL_II_WITH_MUPARSER
#  include <muParser.h>
#endif

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace FunctionParser
  {
    int
    mu_round(const double val)
    {
      return static_cast<int>(val + ((val >= 0.0) ? 0.5 : -0.5));
    }



    double
    mu_if(const double condition,
          const double thenvalue,
          const double elsevalue)
    {
      if (mu_round(condition) != 0)
        return thenvalue;
      else
        return elsevalue;
    }



    double
    mu_or(const double left, const double right)
    {
      return static_cast<double>((mu_round(left) != 0) ||
                                 (mu_round(right) != 0));
    }



    double
    mu_and(const double left, const double right)
    {
      return static_cast<double>((mu_round(left) != 0) &&
                                 (mu_round(right) != 0));
    }



    double
    mu_int(const double value)
    {
      return static_cast<double>(mu_round(value));
    }



    double
    mu_ceil(const double value)
    {
      return std::ceil(value);
    }



    double
    mu_floor(const double value)
    {
      return std::floor(value);
    }



    double
    mu_cot(const double value)
    {
      return 1.0 / std::tan(value);
    }



    double
    mu_csc(const double value)
    {
      return 1.0 / std::sin(value);
    }



    double
    mu_sec(const double value)
    {
      return 1.0 / std::cos(value);
    }



    double
    mu_log(const double value)
    {
      return std::log(value);
    }



    double
    mu_pow(const double a, const double b)
    {
      return std::pow(a, b);
    }



    double
    mu_erf(const double value)
    {
      return std::erf(value);
    }



    double
    mu_erfc(const double value)
    {
      return std::erfc(value);
    }



    // Returns a random value in the range [0,1], after initializing the
    // generator with the given seed
    double
    mu_rand_seed(const double seed)
    {
      static std::mutex           rand_mutex;
      std::lock_guard<std::mutex> lock(rand_mutex);

      std::uniform_real_distribution<> uniform_distribution(0., 1.);

      // for each seed a unique random number generator is created,
      // which is initialized with the seed itself
      static std::map<double, std::mt19937> rng_map;

      return uniform_distribution(
        rng_map.try_emplace(seed, std::mt19937(static_cast<unsigned int>(seed)))
          .first->second);
    }


    // Returns a random value in the range [0,1]
    double
    mu_rand()
    {
      static std::mutex                rand_mutex;
      std::lock_guard<std::mutex>      lock(rand_mutex);
      std::uniform_real_distribution<> uniform_distribution(0., 1.);
      const unsigned int  seed = static_cast<unsigned long>(std::time(nullptr));
      static std::mt19937 rng(seed);
      return uniform_distribution(rng);
    }



    std::vector<std::string>
    get_function_names()
    {
      return {// functions predefined by muparser
              "sin",
              "cos",
              "tan",
              "asin",
              "acos",
              "atan",
              "sinh",
              "cosh",
              "tanh",
              "asinh",
              "acosh",
              "atanh",
              "atan2",
              "log2",
              "log10",
              "log",
              "ln",
              "exp",
              "sqrt",
              "sign",
              "rint",
              "abs",
              "min",
              "max",
              "sum",
              "avg",
              // functions we define ourselves above
              "if",
              "int",
              "ceil",
              "cot",
              "csc",
              "floor",
              "sec",
              "pow",
              "erf",
              "erfc",
              "rand",
              "rand_seed"};
    }

#ifdef DEAL_II_WITH_MUPARSER
    /**
     * PIMPL for mu::Parser.
     */
    class Parser : public muParserBase
    {
    public:
      operator mu::Parser &()
      {
        return parser;
      }

      operator const mu::Parser &() const
      {
        return parser;
      }

    protected:
      mu::Parser parser;
    };
#endif



    template <int dim, typename Number>
    ParserImplementation<dim, Number>::ParserImplementation()
      : initialized(false)
      , n_vars(0)
    {}



    template <int dim, typename Number>
    void
    ParserImplementation<dim, Number>::initialize(
      const std::string                   &variables,
      const std::vector<std::string>      &expressions,
      const std::map<std::string, double> &constants,
      const bool                           time_dependent)
    {
      this->parser_data.clear(); // this will reset all thread-local objects

      this->constants   = constants;
      this->var_names   = Utilities::split_string_list(variables, ',');
      this->expressions = expressions;
      AssertThrow(((time_dependent) ? dim + 1 : dim) == this->var_names.size(),
                  ExcMessage("Wrong number of variables"));

      // Now we define how many variables we expect to read in. We distinguish
      // between two cases: Time dependent problems, and not time dependent
      // problems. In the first case the number of variables is given by the
      // dimension plus one. In the other case, the number of variables is equal
      // to the dimension. Once we parsed the variables string, if none of this
      // is the case, then an exception is thrown.
      if (time_dependent)
        this->n_vars = dim + 1;
      else
        this->n_vars = dim;

      // create a parser object for the current thread we can then query in
      // value() and vector_value(). this is not strictly necessary because a
      // user may never call these functions on the current thread, but it gets
      // us error messages about wrong formulas right away
      this->init_muparser();
      this->initialized = true;
    }



    template <int dim, typename Number>
    void
    ParserImplementation<dim, Number>::init_muparser() const
    {
#ifdef DEAL_II_WITH_MUPARSER
      // check that we have not already initialized the parser on the
      // current thread, i.e., that the current function is only called
      // once per thread
      ParserData &data = this->parser_data.get();
      Assert(data.parsers.empty() && data.vars.empty(), ExcInternalError());
      const unsigned int n_components = expressions.size();

      // initialize the objects for the current thread
      data.parsers.reserve(n_components);
      data.vars.resize(this->var_names.size());
      for (unsigned int component = 0; component < n_components; ++component)
        {
          data.parsers.emplace_back(std::make_unique<Parser>());
          mu::Parser &parser = dynamic_cast<Parser &>(*data.parsers.back());

          for (const auto &constant : this->constants)
            parser.DefineConst(constant.first, constant.second);

          for (unsigned int iv = 0; iv < this->var_names.size(); ++iv)
            parser.DefineVar(this->var_names[iv], &data.vars[iv]);

          // define some compatibility functions:
          parser.DefineFun("if", mu_if, true);
          parser.DefineOprt("|", mu_or, 1);
          parser.DefineOprt("&", mu_and, 2);
          parser.DefineFun("int", mu_int, true);
          parser.DefineFun("ceil", mu_ceil, true);
          parser.DefineFun("cot", mu_cot, true);
          parser.DefineFun("csc", mu_csc, true);
          parser.DefineFun("floor", mu_floor, true);
          parser.DefineFun("sec", mu_sec, true);
          parser.DefineFun("log", mu_log, true);
          parser.DefineFun("pow", mu_pow, true);
          parser.DefineFun("erfc", mu_erfc, true);
          // Disable optimizations (by passing false) that assume the functions
          // will always return the same value:
          parser.DefineFun("rand_seed", mu_rand_seed, false);
          parser.DefineFun("rand", mu_rand, false);

          try
            {
              // muparser expects that functions have no
              // space between the name of the function and the opening
              // parenthesis. this is awkward because it is not backward
              // compatible to the library we used to use before muparser
              // (the fparser library) but also makes no real sense.
              // consequently, in the expressions we set, remove any space
              // we may find after function names
              std::string transformed_expression = this->expressions[component];

              for (const auto &current_function_name : get_function_names())
                {
                  const unsigned int function_name_length =
                    current_function_name.size();

                  std::string::size_type pos = 0;
                  while (true)
                    {
                      // try to find any occurrences of the function name
                      pos =
                        transformed_expression.find(current_function_name, pos);
                      if (pos == std::string::npos)
                        break;

                      // replace whitespace until there no longer is any
                      while (
                        (pos + function_name_length <
                         transformed_expression.size()) &&
                        ((transformed_expression[pos + function_name_length] ==
                          ' ') ||
                         (transformed_expression[pos + function_name_length] ==
                          '\t')))
                        transformed_expression.erase(
                          transformed_expression.begin() + pos +
                          function_name_length);

                      // move the current search position by the size of the
                      // actual function name
                      pos += function_name_length;
                    }
                }

              // now use the transformed expression
              parser.SetExpr(transformed_expression);
            }
          catch (mu::ParserError &e)
            {
              std::cerr << "Message:  <" << e.GetMsg() << ">\n";
              std::cerr << "Formula:  <" << e.GetExpr() << ">\n";
              std::cerr << "Token:    <" << e.GetToken() << ">\n";
              std::cerr << "Position: <" << e.GetPos() << ">\n";
              std::cerr << "Errc:     <" << e.GetCode() << ">" << std::endl;
              AssertThrow(false, ExcParseError(e.GetCode(), e.GetMsg()));
            }
        }
#else
      AssertThrow(false, ExcNeedsFunctionparser());
#endif
    }

    template <int dim, typename Number>
    Number
    ParserImplementation<dim, Number>::do_value(const Point<dim> &p,
                                                const double      time,
                                                unsigned int component) const
    {
#ifdef DEAL_II_WITH_MUPARSER
      Assert(this->initialized == true, ExcNotInitialized());

      // initialize the parser if that hasn't happened yet on the current
      // thread
      internal::FunctionParser::ParserData &data = this->parser_data.get();
      if (data.vars.empty())
        init_muparser();

      for (unsigned int i = 0; i < dim; ++i)
        data.vars[i] = p[i];
      if (dim != this->n_vars)
        data.vars[dim] = time;

      try
        {
          Assert(dynamic_cast<Parser *>(data.parsers[component].get()),
                 ExcInternalError());
          // NOLINTNEXTLINE don't warn about using static_cast once we check
          mu::Parser &parser = static_cast<Parser &>(*data.parsers[component]);
          return parser.Eval();
        } // try
      catch (mu::ParserError &e)
        {
          std::cerr << "Message:  <" << e.GetMsg() << ">\n";
          std::cerr << "Formula:  <" << e.GetExpr() << ">\n";
          std::cerr << "Token:    <" << e.GetToken() << ">\n";
          std::cerr << "Position: <" << e.GetPos() << ">\n";
          std::cerr << "Errc:     <" << e.GetCode() << ">" << std::endl;
          AssertThrow(false, ExcParseError(e.GetCode(), e.GetMsg()));
        } // catch

#else
      (void)p;
      (void)time;
      (void)component;
      AssertThrow(false, ExcNeedsFunctionparser());
#endif
      return std::numeric_limits<double>::signaling_NaN();
    }

    template <int dim, typename Number>
    void
    ParserImplementation<dim, Number>::do_all_values(
      const Point<dim>  &p,
      const double       time,
      ArrayView<Number> &values) const
    {
#ifdef DEAL_II_WITH_MUPARSER
      Assert(this->initialized == true, ExcNotInitialized());

      // initialize the parser if that hasn't happened yet on the current
      // thread
      internal::FunctionParser::ParserData &data = this->parser_data.get();
      if (data.vars.empty())
        init_muparser();

      for (unsigned int i = 0; i < dim; ++i)
        data.vars[i] = p[i];
      if (dim != this->n_vars)
        data.vars[dim] = time;

      AssertDimension(values.size(), data.parsers.size());
      try
        {
          for (unsigned int component = 0; component < data.parsers.size();
               ++component)
            {
              Assert(dynamic_cast<Parser *>(data.parsers[component].get()),
                     ExcInternalError());
              mu::Parser &parser =
                // We just checked that the pointer is valid so suppress the
                // clang-tidy check
                static_cast<Parser &>(*data.parsers[component]); // NOLINT
              values[component] = parser.Eval();
            }
        } // try
      catch (mu::ParserError &e)
        {
          std::cerr << "Message:  <" << e.GetMsg() << ">\n";
          std::cerr << "Formula:  <" << e.GetExpr() << ">\n";
          std::cerr << "Token:    <" << e.GetToken() << ">\n";
          std::cerr << "Position: <" << e.GetPos() << ">\n";
          std::cerr << "Errc:     <" << e.GetCode() << ">" << std::endl;
          AssertThrow(false, ExcParseError(e.GetCode(), e.GetMsg()));
        } // catch
#else
      (void)p;
      (void)time;
      (void)values;
      AssertThrow(false, ExcNeedsFunctionparser());
#endif
    }

// explicit instantiations
#include "base/mu_parser_internal.inst"

  } // namespace FunctionParser
} // namespace internal

DEAL_II_NAMESPACE_CLOSE
