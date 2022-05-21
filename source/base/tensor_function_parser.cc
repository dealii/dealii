// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


#include <deal.II/base/array_view.h>
#include <deal.II/base/mu_parser_internal.h>
#include <deal.II/base/patterns.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/tensor_function_parser.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/utilities.h>

#include <cmath>
#include <map>

#ifdef DEAL_II_WITH_MUPARSER
#  include <muParser.h>
#endif

DEAL_II_NAMESPACE_OPEN


template <int rank, int dim, typename Number>
const std::vector<std::string> &
TensorFunctionParser<rank, dim, Number>::get_expressions() const
{
  return expressions;
}



template <int rank, int dim, typename Number>
TensorFunctionParser<rank, dim, Number>::TensorFunctionParser(
  const double initial_time)
  : TensorFunction<rank, dim, Number>(initial_time)
  , initialized(false)
  , n_vars(0)
  , n_components(Utilities::pow(dim, rank))
{}


template <int rank, int dim, typename Number>
TensorFunctionParser<rank, dim, Number>::TensorFunctionParser(
  const std::string &expression,
  const std::string &constants,
  const std::string &variable_names)
  : TensorFunction<rank, dim, Number>()
  , initialized(false)
  , n_vars(0)
  , n_components(Utilities::pow(dim, rank))
{
  auto constants_map = Patterns::Tools::Convert<ConstMap>::to_value(
    constants,
    Patterns::Map(Patterns::Anything(),
                  Patterns::Double(),
                  0,
                  Patterns::Map::max_int_value,
                  ",",
                  "="));
  initialize(variable_names,
             expression,
             constants_map,
             Utilities::split_string_list(variable_names, ",").size() ==
               dim + 1);
}


#ifdef DEAL_II_WITH_MUPARSER

template <int rank, int dim, typename Number>
void
TensorFunctionParser<rank, dim, Number>::initialize(
  const std::string &                  variables,
  const std::vector<std::string> &     expressions,
  const std::map<std::string, double> &constants,
  const bool                           time_dependent)
{
  this->parser_data.clear(); // this will reset all thread-local objects

  this->constants   = constants;
  this->var_names   = Utilities::split_string_list(variables, ',');
  this->expressions = expressions;
  AssertThrow(((time_dependent) ? dim + 1 : dim) == var_names.size(),
              ExcMessage("Wrong number of variables"));

  // We check that the number of
  // components of this function
  // matches the number of components
  // passed in as a vector of
  // strings.
  AssertThrow(this->n_components == expressions.size(),
              ExcInvalidExpressionSize(this->n_components, expressions.size()));

  // Now we define how many variables
  // we expect to read in.  We
  // distinguish between two cases:
  // Time dependent problems, and not
  // time dependent problems. In the
  // first case the number of
  // variables is given by the
  // dimension plus one. In the other
  // case, the number of variables is
  // equal to the dimension. Once we
  // parsed the variables string, if
  // none of this is the case, then
  // an exception is thrown.
  if (time_dependent)
    n_vars = dim + 1;
  else
    n_vars = dim;

  // create a parser object for the current thread we can then query
  // in value() and vector_value(). this is not strictly necessary
  // because a user may never call these functions on the current
  // thread, but it gets us error messages about wrong formulas right
  // away
  init_muparser();

  // finally set the initialization bit
  initialized = true;
}



namespace internal
{
  namespace TensorFunctionParserImplementation
  {
    namespace
    {
      /**
       * PIMPL for mu::Parser.
       */
      class Parser : public internal::FunctionParser::muParserBase
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
    } // namespace
  }   // namespace TensorFunctionParserImplementation
} // namespace internal



template <int rank, int dim, typename Number>
void
TensorFunctionParser<rank, dim, Number>::init_muparser() const
{
  // check that we have not already initialized the parser on the
  // current thread, i.e., that the current function is only called
  // once per thread
  internal::FunctionParser::ParserData &data = parser_data.get();
  Assert(data.parsers.size() == 0 && data.vars.size() == 0, ExcInternalError());

  // initialize the objects for the current thread
  data.parsers.reserve(this->n_components);
  data.vars.resize(var_names.size());
  for (unsigned int component = 0; component < this->n_components; ++component)
    {
      data.parsers.emplace_back(
        std::make_unique<
          internal::TensorFunctionParserImplementation::Parser>());
      mu::Parser &parser =
        dynamic_cast<internal::TensorFunctionParserImplementation::Parser &>(
          *data.parsers.back());

      for (const auto &constant : constants)
        parser.DefineConst(constant.first, constant.second);

      for (unsigned int iv = 0; iv < var_names.size(); ++iv)
        parser.DefineVar(var_names[iv], &data.vars[iv]);

      // define some compatibility functions:
      parser.DefineFun("if", internal::FunctionParser::mu_if, true);
      parser.DefineOprt("|", internal::FunctionParser::mu_or, 1);
      parser.DefineOprt("&", internal::FunctionParser::mu_and, 2);
      parser.DefineFun("int", internal::FunctionParser::mu_int, true);
      parser.DefineFun("ceil", internal::FunctionParser::mu_ceil, true);
      parser.DefineFun("cot", internal::FunctionParser::mu_cot, true);
      parser.DefineFun("csc", internal::FunctionParser::mu_csc, true);
      parser.DefineFun("floor", internal::FunctionParser::mu_floor, true);
      parser.DefineFun("sec", internal::FunctionParser::mu_sec, true);
      parser.DefineFun("log", internal::FunctionParser::mu_log, true);
      parser.DefineFun("pow", internal::FunctionParser::mu_pow, true);
      parser.DefineFun("erfc", internal::FunctionParser::mu_erfc, true);
      parser.DefineFun("rand_seed",
                       internal::FunctionParser::mu_rand_seed,
                       true);
      parser.DefineFun("rand", internal::FunctionParser::mu_rand, true);

      try
        {
          // muparser expects that functions have no
          // space between the name of the function and the opening
          // parenthesis. this is awkward because it is not backward
          // compatible to the library we used to use before muparser
          // (the fparser library) but also makes no real sense.
          // consequently, in the expressions we set, remove any space
          // we may find after function names
          std::string transformed_expression = expressions[component];

          for (const auto &current_function_name :
               internal::FunctionParser::get_function_names())
            {
              const unsigned int function_name_length =
                current_function_name.size();

              std::string::size_type pos = 0;
              while (true)
                {
                  // try to find any occurrences of the function name
                  pos = transformed_expression.find(current_function_name, pos);
                  if (pos == std::string::npos)
                    break;

                  // replace whitespace until there no longer is any
                  while ((pos + function_name_length <
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
}



template <int rank, int dim, typename Number>
void
TensorFunctionParser<rank, dim, Number>::initialize(
  const std::string &                  vars,
  const std::string &                  expression,
  const std::map<std::string, double> &constants,
  const bool                           time_dependent)
{
  initialize(vars,
             Utilities::split_string_list(expression, ';'),
             constants,
             time_dependent);
}



template <int rank, int dim, typename Number>
Tensor<rank, dim, Number>
TensorFunctionParser<rank, dim, Number>::value(const Point<dim> &p) const
{
  Assert(initialized == true, ExcNotInitialized());

  // initialize the parser if that hasn't happened yet on the current thread
  internal::FunctionParser::ParserData &data = parser_data.get();
  if (data.vars.size() == 0)
    init_muparser();

  for (unsigned int i = 0; i < dim; ++i)
    data.vars[i] = p(i);
  if (dim != n_vars)
    data.vars[dim] = this->get_time();

  std::array<Number, Tensor<rank, dim, Number>::n_independent_components>
    values;

  try
    {
      for (unsigned int component = 0; component < values.size(); ++component)
        {
          Assert(dynamic_cast<
                   internal::TensorFunctionParserImplementation::Parser *>(
                   data.parsers[component].get()),
                 ExcInternalError());
          mu::Parser &parser = static_cast< // NOLINT
            internal::TensorFunctionParserImplementation::Parser &>(
            *data.parsers[component]);
          values[component] = parser.Eval();
        } // for
    }     // try
  catch (mu::ParserError &e)
    {
      std::cerr << "Message:  <" << e.GetMsg() << ">\n";
      std::cerr << "Formula:  <" << e.GetExpr() << ">\n";
      std::cerr << "Token:    <" << e.GetToken() << ">\n";
      std::cerr << "Position: <" << e.GetPos() << ">\n";
      std::cerr << "Errc:     <" << e.GetCode() << ">" << std::endl;
      AssertThrow(false, ExcParseError(e.GetCode(), e.GetMsg()));

      return Tensor<rank, dim, Number>();
    } // catch

  return Tensor<rank, dim, Number>(
    make_array_view(values.begin(), values.end()));
}



template <int rank, int dim, typename Number>
void
TensorFunctionParser<rank, dim, Number>::value_list(
  const std::vector<Point<dim>> &         p,
  std::vector<Tensor<rank, dim, Number>> &values) const
{
  Assert(p.size() == values.size(),
         ExcDimensionMismatch(p.size(), values.size()));

  for (unsigned int i = 0; i < p.size(); ++i)
    {
      values[i] = value(p[i]);
    }
}

#else


template <int rank, int dim, typename Number>
void
TensorFunctionParser<rank, dim, Number>::initialize(
  const std::string &,
  const std::vector<std::string> &,
  const std::map<std::string, double> &,
  const bool)
{
  AssertThrow(false, ExcNeedsFunctionparser());
}

template <int rank, int dim, typename Number>
void
TensorFunctionParser<rank, dim, Number>::initialize(
  const std::string &,
  const std::string &,
  const std::map<std::string, double> &,
  const bool)
{
  AssertThrow(false, ExcNeedsFunctionparser());
}



template <int rank, int dim, typename Number>
Tensor<rank, dim, Number>
TensorFunctionParser<rank, dim, Number>::value(const Point<dim> &) const
{
  AssertThrow(false, ExcNeedsFunctionparser());
  return Tensor<rank, dim, Number>();
}



template <int rank, int dim, typename Number>
void
TensorFunctionParser<rank, dim, Number>::value_list(
  const std::vector<Point<dim>> &,
  std::vector<Tensor<rank, dim, Number>> &) const
{
  AssertThrow(false, ExcNeedsFunctionparser());
}


#endif

// explicit instantiations
#include "tensor_function_parser.inst"

DEAL_II_NAMESPACE_CLOSE
