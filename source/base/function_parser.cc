// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#include <deal.II/base/function_parser.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/vector.h>


#ifdef DEAL_II_WITH_MUPARSER
#include <muParser.h>
#else

namespace fparser
{
  class FunctionParser
  {};
}
#endif

DEAL_II_NAMESPACE_OPEN



template <int dim>
FunctionParser<dim>::FunctionParser(const unsigned int n_components,
                                    const double       initial_time)
  :
  Function<dim>(n_components, initial_time)
{}



template <int dim>
FunctionParser<dim>::~FunctionParser()
{}

#ifdef DEAL_II_WITH_MUPARSER
template <int dim>
void FunctionParser<dim>::initialize (const std::string                   &variables,
                                      const std::vector<std::string>      &expressions,
                                      const std::map<std::string, double> &constants,
                                      const bool time_dependent,
                                      const bool use_degrees)
{
  initialize (variables,
              expressions,
              constants,
              std::map< std::string, double >(),
              time_dependent,
              use_degrees);
}


template <int dim>
void FunctionParser<dim>::initialize (const std::string              &vars,
                                      const std::vector<std::string> &expressions,
                                      const std::map<std::string, double> &constants,
                                      const bool time_dependent)
{
  initialize(vars, expressions, constants, time_dependent, false);
}



namespace internal
{
  // convert double into int
  int mu_round(double val)
  {
    return static_cast<int>(val + ((val>=0.0) ? 0.5 : -0.5) );
  }

  double mu_if(double condition, double thenvalue, double elsevalue)
  {
    if (mu_round(condition))
      return thenvalue;
    else
      return elsevalue;
  }

  double mu_or(double left, double right)
  {
    return (mu_round(left)) || (mu_round(right));
  }

  double mu_and(double left, double right)
  {
    return (mu_round(left)) && (mu_round(right));
  }

  double mu_int(double value)
  {
    return static_cast<double>(mu_round(value));
  }

  double mu_ceil(double value)
  {
    return ceil(value);
  }

  double mu_floor(double value)
  {
    return floor(value);
  }

  double mu_cot(double value)
  {
    return 1.0/tan(value);
  }

  double mu_csc(double value)
  {
    return 1.0/sin(value);
  }

  double mu_sec(double value)
  {
    return 1.0/cos(value);
  }

  double mu_log(double value)
  {
    return log(value);
  }
}


template <int dim>
void FunctionParser<dim>:: init_muparser() const
{
  if (fp.get().size()>0)
    return;

  fp.get().resize(this->n_components);
  vars.get().resize(var_names.size());
  for (unsigned int component=0; component<this->n_components; ++component)
    {
      for (std::map< std::string, double >::const_iterator constant = constants.begin();
           constant != constants.end(); ++constant)
        {
          fp.get()[component].DefineConst(constant->first.c_str(), constant->second);
        }

      for (unsigned int iv=0; iv<var_names.size(); ++iv)
        fp.get()[component].DefineVar(var_names[iv].c_str(), &vars.get()[iv]);

      // define some compatibility functions:
      fp.get()[component].DefineFun("if",internal::mu_if, true);
      fp.get()[component].DefineOprt("|", internal::mu_or, 1);
      fp.get()[component].DefineOprt("&", internal::mu_and, 2);
      fp.get()[component].DefineFun("int", internal::mu_int, true);
      fp.get()[component].DefineFun("ceil", internal::mu_ceil, true);
      fp.get()[component].DefineFun("cot", internal::mu_cot, true);
      fp.get()[component].DefineFun("csc", internal::mu_csc, true);
      fp.get()[component].DefineFun("floor", internal::mu_floor, true);
      fp.get()[component].DefineFun("sec", internal::mu_sec, true);
      fp.get()[component].DefineFun("log", internal::mu_log, true);

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

          const char *function_names[] =
          {
            // functions predefined by muparser
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
            "sec"
          };
          for (unsigned int f=0; f<sizeof(function_names)/sizeof(function_names[0]); ++f)
            {
              const std::string  function_name        = function_names[f];
              const unsigned int function_name_length = function_name.size();

              std::string::size_type pos = 0;
              while (true)
                {
                  // try to find any occurrences of the function name
                  pos = transformed_expression.find (function_name, pos);
                  if (pos == std::string::npos)
                    break;

                  // replace whitespace until there no longer is any
                  while ((pos+function_name_length<transformed_expression.size())
                         &&
                         ((transformed_expression[pos+function_name_length] == ' ')
                          ||
                          (transformed_expression[pos+function_name_length] == '\t')))
                    transformed_expression.erase (transformed_expression.begin()+pos+function_name_length);

                  // move the current search position by the size of the
                  // actual function name
                  pos += function_name_length;
                }
            }

          // now use the transformed expression
          fp.get()[component].SetExpr(transformed_expression);
        }
      catch (mu::ParserError &e)
        {
          std::cerr << "Message:  " << e.GetMsg() << "\n";
          std::cerr << "Formula:  " << e.GetExpr() << "\n";
          std::cerr << "Token:    " << e.GetToken() << "\n";
          std::cerr << "Position: " << e.GetPos() << "\n";
          std::cerr << "Errc:     " << e.GetCode() << std::endl;
          AssertThrow(false, ExcParseError(e.GetCode(), e.GetMsg().c_str()));
        }
    }
}


template <int dim>
void FunctionParser<dim>::initialize (const std::string   &variables,
                                      const std::vector<std::string>      &expressions,
                                      const std::map<std::string, double> &constants,
                                      const std::map<std::string, double> &units,
                                      const bool time_dependent,
                                      const bool use_degrees)
{
  this->fp.clear(); // this will reset all thread-local objects

  this->constants = constants;
  this->var_names = Utilities::split_string_list(variables, ',');
  this->expressions = expressions;
  AssertThrow(((time_dependent)?dim+1:dim) == var_names.size(),
              ExcMessage("wrong number of variables"));
  AssertThrow(!use_degrees, ExcNotImplemented());

  // We check that the number of
  // components of this function
  // matches the number of components
  // passed in as a vector of
  // strings.
  AssertThrow(this->n_components == expressions.size(),
              ExcInvalidExpressionSize(this->n_components,
                                       expressions.size()) );

  // we no longer support units:
  AssertThrow(units.size()==0, ExcNotImplemented());

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
    n_vars = dim+1;
  else
    n_vars = dim;

  init_muparser();

  // Now set the initialization bit.
  initialized = true;
}



template <int dim>
void FunctionParser<dim>::initialize (const std::string &vars,
                                      const std::string &expression,
                                      const std::map<std::string, double> &constants,
                                      const bool time_dependent)
{
  initialize(vars, expression, constants, time_dependent, false);
}



template <int dim>
void FunctionParser<dim>::initialize (const std::string &variables,
                                      const std::string &expression,
                                      const std::map<std::string, double> &constants,
                                      const bool time_dependent,
                                      const bool use_degrees)
{
  // initialize with the things
  // we got.
  initialize (variables,
              Utilities::split_string_list (expression, ';'),
              constants,
              time_dependent,
              use_degrees);
}



template <int dim>
void FunctionParser<dim>::initialize (const std::string &variables,
                                      const std::string &expression,
                                      const std::map<std::string, double> &constants,
                                      const std::map<std::string, double> &units,
                                      const bool time_dependent,
                                      const bool use_degrees)
{
  // initialize with the things
  // we got.
  initialize (variables,
              Utilities::split_string_list (expression, ';'),
              constants,
              units,
              time_dependent,
              use_degrees);
}



template <int dim>
double FunctionParser<dim>::value (const Point<dim>  &p,
                                   const unsigned int component) const
{
  Assert (initialized==true, ExcNotInitialized());
  Assert (component < this->n_components,
          ExcIndexRange(component, 0, this->n_components));

  // initialize if not done so on this thread yet:
  init_muparser();

  for (unsigned int i=0; i<dim; ++i)
    vars.get()[i] = p(i);
  if (dim != n_vars)
    vars.get()[dim] = this->get_time();

  try
    {
      return fp.get()[component].Eval();
    }
  catch (mu::ParserError &e)
    {
      std::cerr << "Message:  " << e.GetMsg() << "\n";
      std::cerr << "Formula:  " << e.GetExpr() << "\n";
      std::cerr << "Token:    " << e.GetToken() << "\n";
      std::cerr << "Position: " << e.GetPos() << "\n";
      std::cerr << "Errc:     " << e.GetCode() << std::endl;
      AssertThrow(false, ExcParseError(e.GetCode(), e.GetMsg().c_str()));
      return 0.0;
    }
}



template <int dim>
void FunctionParser<dim>::vector_value (const Point<dim> &p,
                                        Vector<double>   &values) const
{
  Assert (initialized==true, ExcNotInitialized());
  Assert (values.size() == this->n_components,
          ExcDimensionMismatch (values.size(), this->n_components));


  // initialize if not done so on this thread yet:
  init_muparser();

  for (unsigned int i=0; i<dim; ++i)
    vars.get()[i] = p(i);
  if (dim != n_vars)
    vars.get()[dim] = this->get_time();

  for (unsigned int component = 0; component < this->n_components;
       ++component)
    values(component) = fp.get()[component].Eval();
}

#else


template <int dim>
void
FunctionParser<dim>::initialize(const std::string &,
                                const std::vector<std::string> &,
                                const std::map<std::string, double> &,
                                const bool,
                                const bool)
{
  Assert(false, ExcNeedsFunctionparser());
}


template <int dim>
void
FunctionParser<dim>::initialize(const std::string &,
                                const std::vector<std::string> &,
                                const std::map<std::string, double> &,
                                const std::map<std::string, double> &,
                                const bool,
                                const bool)
{
  Assert(false, ExcNeedsFunctionparser());
}


template <int dim>
void
FunctionParser<dim>::initialize(const std::string &,
                                const std::string &,
                                const std::map<std::string, double> &,
                                const bool,
                                const bool)
{
  Assert(false, ExcNeedsFunctionparser());
}


template <int dim>
void
FunctionParser<dim>::initialize(const std::string &,
                                const std::string &,
                                const std::map<std::string, double> &,
                                const std::map<std::string, double> &,
                                const bool,
                                const bool)
{
  Assert(false, ExcNeedsFunctionparser());
}

template <int dim>
void
FunctionParser<dim>::initialize(const std::string &,
                                const std::vector<std::string> &,
                                const std::map<std::string, double> &,
                                const bool)
{
  Assert(false, ExcNeedsFunctionparser());
}

template <int dim>
void
FunctionParser<dim>::initialize(const std::string &,
                                const std::string &,
                                const std::map<std::string, double> &,
                                const bool)
{
  Assert(false, ExcNeedsFunctionparser());
}



template <int dim>
double FunctionParser<dim>::value (
  const Point<dim> &, unsigned int) const
{
  Assert(false, ExcNeedsFunctionparser());
  return 0.;
}


template <int dim>
void FunctionParser<dim>::vector_value (
  const Point<dim> &, Vector<double> &) const
{
  Assert(false, ExcNeedsFunctionparser());
}


#endif

// Explicit Instantiations.

template class FunctionParser<1>;
template class FunctionParser<2>;
template class FunctionParser<3>;

DEAL_II_NAMESPACE_CLOSE
