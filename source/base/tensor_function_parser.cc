// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2020 by the deal.II authors
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
    std::make_unique<Patterns::Map>(Patterns::Anything(),
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


// We deliberately delay the definition of the default destructor
// so that we don't need to include the definition of mu::Parser
// in the header file.
template <int rank, int dim, typename Number>
TensorFunctionParser<rank, dim, Number>::~TensorFunctionParser() = default;


#ifdef DEAL_II_WITH_MUPARSER

template <int rank, int dim, typename Number>
void
TensorFunctionParser<rank, dim, Number>::initialize(
  const std::string &                  variables,
  const std::vector<std::string> &     expressions,
  const std::map<std::string, double> &constants,
  const bool                           time_dependent)
{
  this->tfp.clear(); // this will reset all thread-local objects

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



template <int rank, int dim, typename Number>
void
TensorFunctionParser<rank, dim, Number>::init_muparser() const
{
  // check that we have not already initialized the parser on the
  // current thread, i.e., that the current function is only called
  // once per thread
  Assert(tfp.get().size() == 0, ExcInternalError());

  // initialize the objects for the current thread (tfp.get() and
  // vars.get())
  tfp.get().reserve(this->n_components);
  vars.get().resize(var_names.size());
  for (unsigned int component = 0; component < this->n_components; ++component)
    {
      tfp.get().emplace_back(new mu::Parser());

      for (const auto &constant : constants)
        {
          tfp.get()[component]->DefineConst(constant.first, constant.second);
        }

      for (unsigned int iv = 0; iv < var_names.size(); ++iv)
        tfp.get()[component]->DefineVar(var_names[iv], &vars.get()[iv]);

      // define some compatibility functions:
      tfp.get()[component]->DefineFun("if",
                                      internal::FunctionParser::mu_if,
                                      true);
      tfp.get()[component]->DefineOprt("|", internal::FunctionParser::mu_or, 1);
      tfp.get()[component]->DefineOprt("&",
                                       internal::FunctionParser::mu_and,
                                       2);
      tfp.get()[component]->DefineFun("int",
                                      internal::FunctionParser::mu_int,
                                      true);
      tfp.get()[component]->DefineFun("ceil",
                                      internal::FunctionParser::mu_ceil,
                                      true);
      tfp.get()[component]->DefineFun("cot",
                                      internal::FunctionParser::mu_cot,
                                      true);
      tfp.get()[component]->DefineFun("csc",
                                      internal::FunctionParser::mu_csc,
                                      true);
      tfp.get()[component]->DefineFun("floor",
                                      internal::FunctionParser::mu_floor,
                                      true);
      tfp.get()[component]->DefineFun("sec",
                                      internal::FunctionParser::mu_sec,
                                      true);
      tfp.get()[component]->DefineFun("log",
                                      internal::FunctionParser::mu_log,
                                      true);
      tfp.get()[component]->DefineFun("pow",
                                      internal::FunctionParser::mu_pow,
                                      true);
      tfp.get()[component]->DefineFun("erfc",
                                      internal::FunctionParser::mu_erfc,
                                      true);
      tfp.get()[component]->DefineFun("rand_seed",
                                      internal::FunctionParser::mu_rand_seed,
                                      true);
      tfp.get()[component]->DefineFun("rand",
                                      internal::FunctionParser::mu_rand,
                                      true);

      try
        {
          // muparser expects that functions have no
          // space between the name of the function and the opening
          // parenthesis. this is awkward because it is not backward
          // compatible to the library we used to use before muparser
          // (the tfparser library) but also makes no real sense.
          // consequently, in the expressions we set, remove any space
          // we may find after function names
          std::string transformed_expression = expressions[component];

          for (const auto &current_function_name :
               internal::FunctionParser::function_names)
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
          tfp.get()[component]->SetExpr(transformed_expression);
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
  if (tfp.get().size() == 0)
    init_muparser();

  for (unsigned int i = 0; i < dim; ++i)
    vars.get()[i] = p(i);
  if (dim != n_vars)
    vars.get()[dim] = this->get_time();

  // initialize tensor with zeros
  Tensor<rank, dim, Number> value;

  try
    {
      unsigned int component = 0;
      for (Number *value_ptr = value.begin_raw(); value_ptr != value.end_raw();
           ++value_ptr)
        {
          *value_ptr = tfp.get()[component]->Eval();
          ++component;
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

  return value;
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
