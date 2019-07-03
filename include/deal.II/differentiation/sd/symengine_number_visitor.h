// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2017 by the deal.II authors
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

#ifndef dealii_differentiation_sd_symengine_number_visitor_h
#define dealii_differentiation_sd_symengine_number_visitor_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SYMENGINE

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
// Low level
#  include <symengine/basic.h>
#  include <symengine/dict.h>
#  include <symengine/symengine_exception.h>
#  include <symengine/symengine_rcp.h>

// Number types
#  include <symengine/complex_double.h>
#  include <symengine/integer.h>
#  include <symengine/number.h>
#  include <symengine/rational.h>
#  include <symengine/real_double.h>
#  include <symengine/symbol.h>

// Visitor
#  include <symengine/visitor.h>

//// Number operations
//#  include <symengine/add.h>
//#  include <symengine/functions.h>
//#  include <symengine/mul.h>
#  include <symengine/pow.h>

//// Evaluation
//#  include <symengine/eval.h>
//#  include <symengine/eval_arb.h>
//#  include <symengine/eval_double.h>

//// Differentiation
//#  include <symengine/derivative.h>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#  include <deal.II/base/exceptions.h>
//#  include <deal.II/base/numbers.h>
//#  include <deal.II/base/symmetric_tensor.h>
//#  include <deal.II/base/table_indices.h>
//#  include <deal.II/base/tensor.h>

//#  include <deal.II/differentiation/sd/symengine_ad_number_types.h>
//#  include <deal.II/differentiation/sd/symengine_la_vector_types.h>
//#  include <deal.II/differentiation/sd/symengine_linear_operator_types.h>
#  include <deal.II/differentiation/sd/symengine_number_types.h>
//#  include <deal.II/differentiation/sd/symengine_number_traits.h>
//#  include <deal.II/differentiation/sd/symengine_utilities.h>

//#  include <deal.II/lac/linear_operator.h>
//
//#  include <algorithm>
//#  include <memory>
//#  include <type_traits>
//#  include <utility>
//#  include <vector>


DEAL_II_NAMESPACE_OPEN


namespace Differentiation
{
  namespace SD
  {
    /* --------------- CommonSubexpressionEliminationVisitor ------------- */


    /**
     * A class that implements common subexpression elimination
     * for dictionary visitor classes.
     *
     * @todo Move to symengine_cse.[h/cc]
     *
     * @author Jean-Paul Pelteret, Isuru Fernando, 2017
     */
    template <typename ReturnType, typename ExpressionType>
    class CSEDictionaryVisitor
    {
      /**
       * Intermediate symbols and their definition
       */
      SymEngine::vec_pair intermediate_symbols_exprs;

      /**
       * Final reduced expressions
       */
      SymEngine::vec_basic reduced_exprs;

    public:
      /**
       * Perform common subexpression elimination
       *
       * Here we build the reduced expressions for the @p dependent_functions
       * as well as a list of intermediate, repeating symbolic expressions
       * that
       * are extracted @p dependent_functions. This operation leads to the
       * elimination of repeated expressions, so they only have to be
       * evaluated once.
       */
      void
      init(const SymEngine::vec_basic &dependent_functions)
      {
        intermediate_symbols_exprs.clear();
        reduced_exprs.clear();

        // After the next call, the data stored in replacements is structured
        // as follows:
        //
        // replacements[i] := [f, f(x)]
        // replacements[i].first  = intermediate function label "f"
        // replacements[i].second = intermediate function definition "f(x)"
        //
        // It is to be evaluated top down (i.e. index 0 to
        // replacements.size()), with the results going back into the
        // substitution map for the next levels. So for each "i", "x" are the
        // superset of the input values and the previously evaluated [f_0(x),
        // f_1(x), ..., f_{i-1}(x)].
        //
        // The final result is a set of reduced expressions
        // that must be computed after the replacement
        // values have been computed.
        SymEngine::cse(intermediate_symbols_exprs,
                       reduced_exprs,
                       dependent_functions);
      }

      void
      call(ReturnType *                output_values,
           const SymEngine::vec_basic &independent_variables,
           const ReturnType *          substitution_values)
      {
        Assert(n_reduced_expressions() > 0, ExcInternalError());

        // First we add the input values into the substitution map...
        SymEngine::map_basic_basic substitution_value_map;
        for (unsigned i = 0; i < independent_variables.size(); ++i)
          substitution_value_map[independent_variables[i]] =
            static_cast<const SymEngine::RCP<const SymEngine::Basic> &>(
              ExpressionType(substitution_values[i]));

        // ... followed by any intermediate evaluations due to the application
        // of CSE. These are fed directly back into the substitution map...
        for (unsigned i = 0; i < intermediate_symbols_exprs.size(); ++i)
          {
            const SymEngine::RCP<const SymEngine::Basic> &cse_symbol =
              intermediate_symbols_exprs[i].first;
            const SymEngine::RCP<const SymEngine::Basic> &cse_expr =
              intermediate_symbols_exprs[i].second;
            Assert(substitution_value_map.find(cse_symbol) ==
                     substitution_value_map.end(),
                   ExcMessage(
                     "Reduced symbol already appears in substitution map. "
                     "Is there a clash between the reduced symbol name and "
                     "the symbol used for an independent variable?"));
            substitution_value_map[cse_symbol] =
              static_cast<const SymEngine::RCP<const SymEngine::Basic> &>(
                ExpressionType(ExpressionType(cse_expr)
                                 .template substitute_and_evaluate<ReturnType>(
                                   substitution_value_map)));
          }

        // ... followed by the final reduced expressions
        for (unsigned i = 0; i < reduced_exprs.size(); ++i)
          output_values[i] = ExpressionType(reduced_exprs[i])
                               .template substitute_and_evaluate<ReturnType>(
                                 substitution_value_map);
      }

      template <typename StreamType>
      void
      print(StreamType &stream,
            const bool /*print_independent_symbols*/ = false,
            const bool /*print_dependent_functions*/ = false,
            const bool /*print_cse_reductions*/      = false) const
      {
        stream << "Common subexpression elimination: \n";
        stream << "  Intermediate reduced expressions: \n";
        for (unsigned i = 0; i < intermediate_symbols_exprs.size(); ++i)
          {
            const SymEngine::RCP<const SymEngine::Basic> &cse_symbol =
              intermediate_symbols_exprs[i].first;
            const SymEngine::RCP<const SymEngine::Basic> &cse_expr =
              intermediate_symbols_exprs[i].second;
            stream << "  " << i << ": " << cse_symbol << " = " << cse_expr
                   << "\n";
          }

        stream << "  Final reduced expressions for dependent variables: \n";
        for (unsigned i = 0; i < reduced_exprs.size(); ++i)
          stream << "  " << i << ": " << reduced_exprs[i] << "\n";

        stream << std::flush;
      }

      /**
       * Return a flag stating whether we've performed CSE or not.
       */
      bool
      executed() const
      {
        // For dictionary substitution, the CSE algorithm moves
        // ownership of the dependent function expression definition
        // to the entries in reduced_exprs. So its size thus determines
        // whether CSE has been executed or not.
        return (n_reduced_expressions() > 0) ||
               (n_intermediate_expressions() > 0);
      }

      /**
       * The number of intermediate expressions that must
       * be evaluated as part of the collection of common
       * subexpressions.
       */
      unsigned int
      n_intermediate_expressions() const
      {
        return intermediate_symbols_exprs.size();
      }

      /**
       * The size of the final set of reduced expressions.
       */
      unsigned int
      n_reduced_expressions() const
      {
        return reduced_exprs.size();
      }
    };


    /* ------------------- DictionarySubstitutionVisitor ----------------- */

    //#define DEBUG_EXTRA_CSE

    /**
     * A class to perform dictionary-based substitution as
     * if it were an optimizer of the "lambda" or "LLVM"
     * variety. This class is only really useful to assist
     * in the easy switching between different optimizers,
     * and more importantly for integrating CSE into a
     * dictionary substitution scheme.
     */
    template <typename ReturnType, typename ExpressionType>
    class DictionarySubstitutionVisitor
      : public SymEngine::BaseVisitor<
          DictionarySubstitutionVisitor<ReturnType, ExpressionType>>
    {
    protected:
      SymEngine::vec_basic                             independent_variables;
      SymEngine::vec_basic                             dependent_variables;
      CSEDictionaryVisitor<ReturnType, ExpressionType> cse;

    public:
      void
      init(const SymEngine::vec_basic &x,
           const SymEngine::Basic &    b,
           const bool                  use_cse = false)
      {
        init(x, {b.rcp_from_this()}, use_cse);
      }

      void
      init(const SymEngine::vec_basic &inputs,
           const SymEngine::vec_basic &outputs,
           const bool                  use_cse = false)
      {
        independent_variables.clear();
        dependent_variables.clear();

        independent_variables = inputs;

        // Perform common subexpression elimination if requested
        // Note: After this is done, the results produced by
        // dependent_variables and cse.reduced_exprs should be
        // the same. We could keep the former so that we can print
        // out the original expressions if we wish to do so.
        if (use_cse == false)
          dependent_variables = outputs;
        else
          {
            cse.init(outputs);

#  if defined(DEBUG) && defined(DEBUG_EXTRA_CSE)
            // Check the result of CSE by directly evaluating
            // the dependent variables on an individual basis
            dependent_variables = outputs;
#  endif
          }
      }

      ReturnType
      call(const std::vector<ReturnType> &substitution_values)
      {
        Assert(
          dependent_variables.size() == 1,
          ExcMessage(
            "Cannot use this call function when more than one symbolic expression is to be evaluated."));
        Assert(
          substitution_values.size() == independent_variables.size(),
          ExcMessage(
            "Input substitution vector does not match size of symbol vector."));

        ReturnType out = dealii::internal::NumberType<ReturnType>::value(0.0);
        call(&out, substitution_values.data());
        return out;
      }

      void
      call(ReturnType *output_values, const ReturnType *substitution_values)
      {
        // Check to see if CSE has been performed
        if (cse.executed())
          {
            //            Assert(cse.n_reduced_expressions() ==
            //            dependent_variables.size(),
            //                ExcDimensionMismatch(n_reduced_expressions.size(),
            //                dependent_variables.size()));
            cse.call(output_values, independent_variables, substitution_values);

#  if defined(DEBUG) && defined(DEBUG_EXTRA_CSE)
            // Directly check the result of the CSE calculations
            Assert(dependent_variables.size() > 0, ExcInternalError());
            SymEngine::map_basic_basic substitution_value_map;
            for (unsigned i = 0; i < independent_variables.size(); ++i)
              substitution_value_map[independent_variables[i]] =
                static_cast<const SymEngine::RCP<const SymEngine::Basic> &>(
                  ExpressionType(substitution_values[i]));
            const ReturnType tol = ReturnType(1e-9);
            for (unsigned i = 0; i < dependent_variables.size(); ++i)
              {
                const ReturnType result =
                  ExpressionType(dependent_variables[i])
                    .template substitute_and_evaluate<ReturnType>(
                      substitution_value_map);
                if (std::abs(output_values[i] - result) / std::abs(result) >=
                    tol)
                  {
                    std::ostringstream stream;
                    stream
                      << "A CSE based evaluation is incorrect. The result should have been "
                      << result << " but was rather calculated to be "
                      << output_values[i] << ".\n Substitution map: \n";
                    for (typename SymEngine::map_basic_basic::const_iterator
                           it = substitution_value_map.begin();
                         it != substitution_value_map.end();
                         ++it)
                      {
                        stream << *(it->first) << " = " << *(it->second)
                               << "\n";
                      }
                    stream << std::endl;
                    const std::string msg = stream.str();
                    Assert(std::abs(output_values[i] - result) /
                               std::abs(result) <
                             tol,
                           ExcMessage(msg));
                  }
              }
#  endif
          }
        else
          {
            // Build a substitution map.
            SymEngine::map_basic_basic substitution_value_map;
            for (unsigned i = 0; i < independent_variables.size(); ++i)
              substitution_value_map[independent_variables[i]] =
                static_cast<const SymEngine::RCP<const SymEngine::Basic> &>(
                  ExpressionType(substitution_values[i]));

            // Since we don't know how to definitively evaluate the
            // input number type, we create a generic Expression
            // with the given symbolic expression and ask it to perform
            // substitution and evaluation for us.
            Assert(dependent_variables.size() > 0, ExcInternalError());
            for (unsigned i = 0; i < dependent_variables.size(); ++i)
              output_values[i] =
                ExpressionType(dependent_variables[i])
                  .template substitute_and_evaluate<ReturnType>(
                    substitution_value_map);
          }
      }

      template <typename StreamType>
      void
      print(StreamType &stream,
            const bool  print_independent_symbols = false,
            const bool  print_dependent_functions = false,
            const bool  print_cse_reductions      = false) const
      {
        if (print_independent_symbols)
          {
            stream << "Independent variables: \n";
            for (unsigned i = 0; i < independent_variables.size(); ++i)
              stream << "  " << i << ": " << independent_variables[i] << "\n";

            stream << std::flush;
          }

        // Check to see if CSE has been performed
        if (print_cse_reductions && cse.executed())
          {
            cse.print(stream);
          }
        else
          {
            Assert(dependent_variables.size() > 0, ExcInternalError());

            if (print_dependent_functions)
              {
                stream << "Dependent variables: \n";
                for (unsigned i = 0; i < dependent_variables.size(); ++i)
                  stream << "  " << i << dependent_variables[i] << "\n";

                stream << std::flush;
              }
          }
      }

      void
      bvisit(const SymEngine::Basic &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::Symbol &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::Constant &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::Integer &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::Rational &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::RealDouble &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::Complex &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::ComplexDouble &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::Add &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::Mul &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::Pow &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::Log &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::Sin &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::Cos &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::Tan &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::Csc &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::Sec &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::Cot &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::ASin &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::ACos &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::ATan &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::ATan2 &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::ACsc &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::ASec &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::ACot &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::Sinh &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::Cosh &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::Tanh &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::Csch &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::Sech &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::Coth &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::ASinh &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::ACosh &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::ATanh &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::ACsch &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::ACoth &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::ASech &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::Abs &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::Gamma &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::LogGamma &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::Erf &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::Erfc &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::Max &)
      {
        AssertThrow(false, ExcNotImplemented());
      }

      void
      bvisit(const SymEngine::Min &)
      {
        AssertThrow(false, ExcNotImplemented());
      }
    };

  } // namespace SD
} // namespace Differentiation


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SYMENGINE

#endif // dealii_differentiation_sd_symengine_number_visitor_h
