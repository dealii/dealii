// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SYMENGINE

#  include <deal.II/differentiation/sd/symengine_optimizer.h>
#  include <deal.II/differentiation/sd/symengine_utilities.h>

#  include <boost/archive/text_iarchive.hpp>
#  include <boost/archive/text_oarchive.hpp>

#  include <utility>

DEAL_II_NAMESPACE_OPEN


namespace Differentiation
{
  namespace SD
  {
    template <typename ReturnType>
    BatchOptimizer<ReturnType>::BatchOptimizer()
      : method(OptimizerType::dictionary)
      , flags(OptimizationFlags::optimize_default)
      , ready_for_value_extraction(false)
      , has_been_serialized(false)
    {}



    template <typename ReturnType>
    BatchOptimizer<ReturnType>::BatchOptimizer(
      const enum OptimizerType &    optimization_method,
      const enum OptimizationFlags &optimization_flags)
      : BatchOptimizer()
    {
      set_optimization_method(optimization_method, optimization_flags);
    }



    template <typename ReturnType>
    BatchOptimizer<ReturnType>::BatchOptimizer(
      const BatchOptimizer<ReturnType> &other)
      : method(other.method)
      , flags(other.flags)
      , independent_variables_symbols(other.independent_variables_symbols)
      , dependent_variables_functions(other.dependent_variables_functions)
      , dependent_variables_output(0)
      , map_dep_expr_vec_entry(other.map_dep_expr_vec_entry)
      , ready_for_value_extraction(false)
    {}



    template <typename ReturnType>
    void
    BatchOptimizer<ReturnType>::set_optimization_method(
      const enum OptimizerType &    optimization_method,
      const enum OptimizationFlags &optimization_flags)
    {
      Assert(
        optimized() == false,
        ExcMessage(
          "Cannot call set_optimization_method() once the optimizer is finalized."));

#  ifndef HAVE_SYMENGINE_LLVM
      if (optimization_method == OptimizerType::llvm)
        {
          AssertThrow(false, ExcSymEngineLLVMNotAvailable());
        }
#  endif
      method = optimization_method;
      flags  = optimization_flags;
    }



    template <typename ReturnType>
    enum OptimizerType
    BatchOptimizer<ReturnType>::optimization_method() const
    {
      return method;
    }



    template <typename ReturnType>
    enum OptimizationFlags
    BatchOptimizer<ReturnType>::optimization_flags() const
    {
      return flags;
    }



    template <typename ReturnType>
    bool
    BatchOptimizer<ReturnType>::use_symbolic_CSE() const
    {
      return internal::use_symbolic_CSE(flags);
    }



    template <typename ReturnType>
    bool
    BatchOptimizer<ReturnType>::optimized() const
    {
      if (dependent_variables_output.size() > 0)
        {
          Assert(dependent_variables_output.size() ==
                   dependent_variables_functions.size(),
                 ExcInternalError());
          return true;
        }

      return false;
    }



    template <typename ReturnType>
    bool
    BatchOptimizer<ReturnType>::values_substituted() const
    {
      return ready_for_value_extraction;
    }



    template <typename ReturnType>
    void
    BatchOptimizer<ReturnType>::register_symbols(
      const SD::types::substitution_map &substitution_map)
    {
      Assert(optimized() == false,
             ExcMessage(
               "Cannot register symbols once the optimizer is finalized."));

#  ifdef DEBUG
      // Ensure that all of the keys in the map are actually symbolic
      // in nature
      for (const auto &entry : substitution_map)
        {
          const SD::Expression &symbol = entry.first;
          Assert(SymEngine::is_a<SymEngine::Symbol>(*(symbol.get_RCP())),
                 ExcMessage("Key entry in map is not a symbol."));
        }
#  endif
      // Merge the two maps, in the process ensuring that there is no
      // duplication of symbols
      independent_variables_symbols.insert(substitution_map.begin(),
                                           substitution_map.end());
    }



    template <typename ReturnType>
    void
    BatchOptimizer<ReturnType>::register_symbols(
      const SymEngine::map_basic_basic &substitution_map)
    {
      register_symbols(
        SD::Utilities::convert_basic_map_to_expression_map(substitution_map));
    }



    template <typename ReturnType>
    void
    BatchOptimizer<ReturnType>::register_symbols(
      const SD::types::symbol_vector &symbols)
    {
      Assert(optimized() == false,
             ExcMessage(
               "Cannot register symbols once the optimizer is finalized."));

      for (const auto &symbol : symbols)
        {
          Assert(independent_variables_symbols.find(symbol) ==
                   independent_variables_symbols.end(),
                 ExcMessage("Symbol is already in the map."));
          independent_variables_symbols.insert(
            std::make_pair(symbol, SD::Expression(0.0)));
        }
    }



    template <typename ReturnType>
    void
    BatchOptimizer<ReturnType>::register_symbols(
      const SymEngine::vec_basic &symbols)
    {
      register_symbols(
        SD::Utilities::convert_basic_vector_to_expression_vector(symbols));
    }



    template <typename ReturnType>
    SD::types::symbol_vector
    BatchOptimizer<ReturnType>::get_independent_symbols(void) const
    {
      return Utilities::extract_symbols(independent_variables_symbols);
    }



    template <typename ReturnType>
    std::size_t
    BatchOptimizer<ReturnType>::n_independent_variables(void) const
    {
      return independent_variables_symbols.size();
    }



    template <typename ReturnType>
    void
    BatchOptimizer<ReturnType>::register_function(const Expression &function)
    {
      Assert(optimized() == false,
             ExcMessage(
               "Cannot register functions once the optimizer is finalized."));

      register_scalar_function(function);
    }



    template <typename ReturnType>
    void
    BatchOptimizer<ReturnType>::register_functions(
      const SD::types::symbol_vector &functions)
    {
      Assert(optimized() == false,
             ExcMessage(
               "Cannot register functions once the optimizer is finalized."));

      register_vector_functions(functions);
    }



    template <typename ReturnType>
    void
    BatchOptimizer<ReturnType>::register_functions(
      const SymEngine::vec_basic &functions)
    {
      register_functions(
        Utilities::convert_basic_vector_to_expression_vector(functions));
    }



    template <typename ReturnType>
    const SD::types::symbol_vector &
    BatchOptimizer<ReturnType>::get_dependent_functions(void) const
    {
      return dependent_variables_functions;
    }



    template <typename ReturnType>
    std::size_t
    BatchOptimizer<ReturnType>::n_dependent_variables(void) const
    {
      if (has_been_serialized == false)
        {
          // If we've had to augment our map after serialization, then
          // this check, unfortunately, cannot be performed.
          Assert(map_dep_expr_vec_entry.size() ==
                   dependent_variables_functions.size(),
                 ExcInternalError());
        }
      return dependent_variables_functions.size();
    }



    template <typename ReturnType>
    void
    BatchOptimizer<ReturnType>::optimize()
    {
      Assert(optimized() == false,
             ExcMessage("Cannot call optimize() more than once."));

      // Create and configure the optimizer
      create_optimizer(optimizer);
      Assert(optimizer, ExcNotInitialized());

      const SD::types::symbol_vector symbol_vec =
        Utilities::extract_symbols(independent_variables_symbols);
      if (typename internal::DictionaryOptimizer<ReturnType>::OptimizerType
            *opt = dynamic_cast<typename internal::DictionaryOptimizer<
              ReturnType>::OptimizerType *>(optimizer.get()))
        {
          Assert(optimization_method() == OptimizerType::dictionary,
                 ExcInternalError());
          internal::OptimizerHelper<ReturnType,
                                    internal::DictionaryOptimizer<ReturnType>>::
            initialize(opt,
                       Utilities::convert_expression_vector_to_basic_vector(
                         symbol_vec),
                       Utilities::convert_expression_vector_to_basic_vector(
                         dependent_variables_functions),
                       optimization_flags());
        }
      else if (typename internal::LambdaOptimizer<ReturnType>::OptimizerType
                 *opt = dynamic_cast<typename internal::LambdaOptimizer<
                   ReturnType>::OptimizerType *>(optimizer.get()))
        {
          Assert(optimization_method() == OptimizerType::lambda,
                 ExcInternalError());
          internal::OptimizerHelper<ReturnType,
                                    internal::LambdaOptimizer<ReturnType>>::
            initialize(opt,
                       Utilities::convert_expression_vector_to_basic_vector(
                         symbol_vec),
                       Utilities::convert_expression_vector_to_basic_vector(
                         dependent_variables_functions),
                       optimization_flags());
        }
#  ifdef HAVE_SYMENGINE_LLVM
      else if (typename internal::LLVMOptimizer<ReturnType>::OptimizerType
                 *opt = dynamic_cast<typename internal::LLVMOptimizer<
                   ReturnType>::OptimizerType *>(optimizer.get()))
        {
          Assert(optimization_method() == OptimizerType::llvm,
                 ExcInternalError());
          internal::OptimizerHelper<ReturnType,
                                    internal::LLVMOptimizer<ReturnType>>::
            initialize(opt,
                       Utilities::convert_expression_vector_to_basic_vector(
                         symbol_vec),
                       Utilities::convert_expression_vector_to_basic_vector(
                         dependent_variables_functions),
                       optimization_flags());
        }
#  endif
      else
        {
          AssertThrow(false, ExcMessage("Unknown optimizer type."));
        }

      // The size of the outputs is now fixed, as is the number and
      // order of the symbols to be substituted.
      // Note: When no optimisation is actually used (i.e. optimization_method()
      // == off and use_symbolic_CSE() == false), we could conceptually go
      // without this data structure. However, since the user expects to perform
      // substitution of all dependent variables in one go, we still require it
      // for intermediate storage of results.
      dependent_variables_output.resize(n_dependent_variables());
    }



    template <typename ReturnType>
    void
    BatchOptimizer<ReturnType>::substitute(
      const SD::types::substitution_map &substitution_map) const
    {
      Assert(
        optimized() == true,
        ExcMessage(
          "The optimizer is not configured to perform substitution. "
          "This action can only performed after optimize() has been called."));
      Assert(optimizer, ExcNotInitialized());

      // Check that the registered symbol map and the input map are compatible
      // with one another
#  ifdef DEBUG
      const SD::types::symbol_vector symbol_sub_vec =
        Utilities::extract_symbols(substitution_map);
      const SD::types::symbol_vector symbol_vec =
        Utilities::extract_symbols(independent_variables_symbols);
      Assert(symbol_sub_vec.size() == symbol_vec.size(),
             ExcDimensionMismatch(symbol_sub_vec.size(), symbol_vec.size()));
      for (unsigned int i = 0; i < symbol_sub_vec.size(); ++i)
        {
          Assert(numbers::values_are_equal(symbol_sub_vec[i], symbol_vec[i]),
                 ExcMessage(
                   "The input substitution map is either incomplete, or does "
                   "not match that used in the register_symbols() call."));
        }
#  endif

      // Extract the values from the substitution map, and use the other
      // function
      const std::vector<ReturnType> values =
        Utilities::extract_values<ReturnType>(substitution_map);
      substitute(values);
    }



    template <typename ReturnType>
    void
    BatchOptimizer<ReturnType>::substitute(
      const SymEngine::map_basic_basic &substitution_map) const
    {
      substitute(
        SD::Utilities::convert_basic_map_to_expression_map(substitution_map));
    }



    template <typename ReturnType>
    void
    BatchOptimizer<ReturnType>::substitute(
      const SD::types::symbol_vector &symbols,
      const std::vector<ReturnType> & values) const
    {
      // Zip the two vectors and use the other function call
      // This ensures the ordering of the input vectors matches that of the
      // stored map.
      substitute(make_substitution_map(symbols, values));
    }



    template <typename ReturnType>
    void
    BatchOptimizer<ReturnType>::substitute(
      const SymEngine::vec_basic &   symbols,
      const std::vector<ReturnType> &values) const
    {
      substitute(SD::Utilities::convert_basic_vector_to_expression_vector(
                   symbols),
                 values);
    }



    template <typename ReturnType>
    void
    BatchOptimizer<ReturnType>::substitute(
      const std::vector<ReturnType> &substitution_values) const
    {
      Assert(
        optimized() == true,
        ExcMessage(
          "The optimizer is not configured to perform substitution. "
          "This action can only performed after optimize() has been called."));
      Assert(optimizer, ExcNotInitialized());
      Assert(substitution_values.size() == independent_variables_symbols.size(),
             ExcDimensionMismatch(substitution_values.size(),
                                  independent_variables_symbols.size()));

      if (typename internal::DictionaryOptimizer<ReturnType>::OptimizerType
            *opt = dynamic_cast<typename internal::DictionaryOptimizer<
              ReturnType>::OptimizerType *>(optimizer.get()))
        {
          Assert(optimization_method() == OptimizerType::dictionary,
                 ExcInternalError());
          internal::OptimizerHelper<ReturnType,
                                    internal::DictionaryOptimizer<ReturnType>>::
            substitute(opt, dependent_variables_output, substitution_values);
        }
      else if (typename internal::LambdaOptimizer<ReturnType>::OptimizerType
                 *opt = dynamic_cast<typename internal::LambdaOptimizer<
                   ReturnType>::OptimizerType *>(optimizer.get()))
        {
          Assert(optimization_method() == OptimizerType::lambda,
                 ExcInternalError());
          internal::OptimizerHelper<ReturnType,
                                    internal::LambdaOptimizer<ReturnType>>::
            substitute(opt, dependent_variables_output, substitution_values);
        }
#  ifdef HAVE_SYMENGINE_LLVM
      else if (typename internal::LLVMOptimizer<ReturnType>::OptimizerType
                 *opt = dynamic_cast<typename internal::LLVMOptimizer<
                   ReturnType>::OptimizerType *>(optimizer.get()))
        {
          Assert(optimization_method() == OptimizerType::llvm,
                 ExcInternalError());
          internal::OptimizerHelper<ReturnType,
                                    internal::LLVMOptimizer<ReturnType>>::
            substitute(opt, dependent_variables_output, substitution_values);
        }
#  endif
      else
        {
          AssertThrow(false, ExcNotImplemented());
        }

      ready_for_value_extraction = true;
    }



    template <typename ReturnType>
    const std::vector<ReturnType> &
    BatchOptimizer<ReturnType>::evaluate() const
    {
      Assert(
        values_substituted() == true,
        ExcMessage(
          "The optimizer is not configured to perform evaluation. "
          "This action can only performed after substitute() has been called."));

      return dependent_variables_output;
    }



    template <typename ReturnType>
    ReturnType
    BatchOptimizer<ReturnType>::evaluate(const Expression &func) const
    {
      Assert(
        values_substituted() == true,
        ExcMessage(
          "The optimizer is not configured to perform evaluation. "
          "This action can only performed after substitute() has been called."));

      // TODO[JPP]: Find a way to fix this bug that crops up in serialization
      // cases, e.g. symengine/batch_optimizer_05. Even though the entry is
      // in the map, it can only be found by an exhaustive search and string
      // comparison. Why? Because the leading zero coefficient may seemingly
      // be dropped (or added) at any time.
      //
      // Just this should theoretically work:
      const typename map_dependent_expression_to_vector_entry_t::const_iterator
        it = map_dep_expr_vec_entry.find(func);

      // But instead we are forced to live with this abomination, and its
      // knock-on effects:
      if (has_been_serialized && it == map_dep_expr_vec_entry.end())
        {
          // Some SymEngine operations might return results with a zero leading
          // coefficient. Upon serialization, this might be dropped, meaning
          // that when we reload the expressions they now look somewhat
          // different to as before. If all data that the user uses is
          // guaranteed to either have been serialized or never serialized, then
          // there would be no problem. However, users might rebuild their
          // dependent expression and just reload the optimizer. This is
          // completely legitimate. But in this scenario we might be out of sync
          // with the expressions. This is not great. So we take the nuclear
          // approach, and run everything through a serialization operation to
          // see if we can homogenize all of the expressions such that they look
          // the same in string form.
          auto serialize_and_deserialize_expression =
            [](const Expression &old_expr) {
              std::ostringstream oss;
              {
                boost::archive::text_oarchive oa(oss,
                                                 boost::archive::no_header);
                oa << old_expr;
              }

              Expression new_expr;
              {
                std::istringstream            iss(oss.str());
                boost::archive::text_iarchive ia(iss,
                                                 boost::archive::no_header);

                ia >> new_expr;
              }

              return new_expr;
            };

          const Expression new_func =
            serialize_and_deserialize_expression(func);

          // Find this in the map, while also making sure to compactify all map
          // entries. If we find the entry that we're looking for, then we
          // (re-)add the input expression into the map, and do the proper
          // search again. We should only need to do this once per invalid
          // entry, as the corrected entry is then cached in the map.
          for (const auto &e : map_dep_expr_vec_entry)
            {
              const Expression new_map_expr =
                serialize_and_deserialize_expression(e.first);

              // Add a new map entry and re-search. This is guaranteed to
              // return a valid entry. Note that we must do a string comparison,
              // because the data structures that form the expressions might
              // still be different.
              if (new_func.get_value().__str__() ==
                  new_map_expr.get_value().__str__())
                {
                  map_dep_expr_vec_entry[func] = e.second;
                  return evaluate(func);
                }
            }

          AssertThrow(
            false,
            ExcMessage(
              "Still cannot find map entry, and there's no hope to recover from this situation."));
        }

      Assert(it != map_dep_expr_vec_entry.end(),
             ExcMessage("Function has not been registered."));
      Assert(it->second < n_dependent_variables(), ExcInternalError());

      return dependent_variables_output[it->second];
    }



    template <typename ReturnType>
    std::vector<ReturnType>
    BatchOptimizer<ReturnType>::evaluate(
      const std::vector<Expression> &funcs) const
    {
      std::vector<ReturnType> out;
      out.reserve(funcs.size());

      for (const auto &func : funcs)
        out.emplace_back(evaluate(func));

      return out;
    }



    template <typename ReturnType>
    bool
    BatchOptimizer<ReturnType>::is_valid_nonunique_dependent_variable(
      const SD::Expression &func) const
    {
      return is_valid_nonunique_dependent_variable(func.get_RCP());
    }



    template <typename ReturnType>
    bool
    BatchOptimizer<ReturnType>::is_valid_nonunique_dependent_variable(
      const SymEngine::RCP<const SymEngine::Basic> &func) const
    {
      // SymEngine's internal constants are the valid
      // reusable return types for various derivative operations
      // See
      // https://github.com/symengine/symengine/blob/master/symengine/constants.h
      if (SymEngine::is_a<SymEngine::Constant>(*func))
        return true;
      if (&*func == &*SymEngine::zero)
        return true;
      if (&*func == &*SymEngine::one)
        return true;
      if (&*func == &*SymEngine::minus_one)
        return true;
      if (&*func == &*SymEngine::I)
        return true;
      if (&*func == &*SymEngine::Inf)
        return true;
      if (&*func == &*SymEngine::NegInf)
        return true;
      if (&*func == &*SymEngine::ComplexInf)
        return true;
      if (&*func == &*SymEngine::Nan)
        return true;

      return false;
    }



    template <typename ReturnType>
    void
    BatchOptimizer<ReturnType>::register_scalar_function(
      const SD::Expression &func)
    {
      Assert(
        dependent_variables_output.size() == 0,
        ExcMessage(
          "Cannot register function as the optimizer has already been finalized."));
      dependent_variables_output.reserve(n_dependent_variables() + 1);
      const bool entry_registered =
        (map_dep_expr_vec_entry.find(func) != map_dep_expr_vec_entry.end());
#  ifdef DEBUG
      if (entry_registered == true &&
          is_valid_nonunique_dependent_variable(func) == false)
        Assert(entry_registered,
               ExcMessage("Function has already been registered."));
#  endif
      if (entry_registered == false)
        {
          dependent_variables_functions.push_back(func);
          map_dep_expr_vec_entry[func] =
            dependent_variables_functions.size() - 1;
        }
    }



    template <typename ReturnType>
    void
    BatchOptimizer<ReturnType>::register_vector_functions(
      const SD::types::symbol_vector &funcs)
    {
      Assert(
        dependent_variables_output.size() == 0,
        ExcMessage(
          "Cannot register function as the optimizer has already been finalized."));
      const std::size_t n_dependents_old = n_dependent_variables();
      dependent_variables_output.reserve(n_dependents_old + funcs.size());
      dependent_variables_functions.reserve(n_dependents_old + funcs.size());

      for (const auto &func : funcs)
        {
          const bool entry_registered =
            (map_dep_expr_vec_entry.find(func) != map_dep_expr_vec_entry.end());
#  ifdef DEBUG
          if (entry_registered == true &&
              is_valid_nonunique_dependent_variable(func) == false)
            Assert(entry_registered,
                   ExcMessage("Function has already been registered."));
#  endif
          if (entry_registered == false)
            {
              dependent_variables_functions.push_back(func);
              map_dep_expr_vec_entry[func] =
                dependent_variables_functions.size() - 1;
            }
        }
    }



    template <typename ReturnType>
    void
    BatchOptimizer<ReturnType>::create_optimizer(
      std::unique_ptr<SymEngine::Visitor> &optimizer)
    {
      Assert(!optimizer, ExcMessage("Optimizer has already been created."));

      if (optimization_method() == OptimizerType::dictionary ||
          optimization_method() == OptimizerType::dictionary)
        {
          using Optimizer_t =
            typename internal::DictionaryOptimizer<ReturnType>::OptimizerType;
          optimizer.reset(new Optimizer_t());
        }
      else if (optimization_method() == OptimizerType::lambda)
        {
          using Optimizer_t =
            typename internal::LambdaOptimizer<ReturnType>::OptimizerType;
          optimizer.reset(new Optimizer_t());
        }
      else if (optimization_method() == OptimizerType::llvm)
        {
#  ifdef HAVE_SYMENGINE_LLVM
          if (internal::LLVMOptimizer<ReturnType>::supported_by_LLVM)
            {
              using Optimizer_t =
                typename internal::LLVMOptimizer<ReturnType>::OptimizerType;
              optimizer.reset(new Optimizer_t());
            }
          else
            {
              AssertThrow(false, ExcSymEngineLLVMReturnTypeNotSupported());
            }
#  else
          AssertThrow(false, ExcSymEngineLLVMNotAvailable());
#  endif
        }
      else
        {
          AssertThrow(false, ExcMessage("Unknown optimizer selected."));
        }
    }

  } // namespace SD
} // namespace Differentiation


/* --- Explicit instantiations --- */
#  include "symengine_optimizer.inst"


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SYMENGINE
