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

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SYMENGINE

#  include <deal.II/differentiation/sd/symengine_optimizer.h>
#  include <deal.II/differentiation/sd/symengine_utilities.h>

#  include <utility>

DEAL_II_NAMESPACE_OPEN


namespace Differentiation
{
  namespace SD
  {
    template <typename ReturnType>
    BatchOptimizer<ReturnType>::BatchOptimizer()
      : method(OptimizerType::off)
      , flags(OptimizationFlags::optimize_default)
      , ready_for_value_extraction(false)
    {}


    template <typename ReturnType>
    BatchOptimizer<ReturnType>::BatchOptimizer(
      const enum OptimizerType &    method,
      const enum OptimizationFlags &flags)
      : BatchOptimizer()
    {
      set_optimization_method(method, flags);
    }


    template <typename ReturnType>
    BatchOptimizer<ReturnType>::BatchOptimizer(const BatchOptimizer<ReturnType> &other/*,
                                               const bool                        copy_initialized*/)
      : method (other.method),
        flags (other.flags),
        independent_variables_symbols (other.independent_variables_symbols),
        dependent_variables_functions (other.dependent_variables_functions),
        dependent_variables_output (0),
        map_dep_ptr_basic_vec_entry (other.map_dep_ptr_basic_vec_entry),
        ready_for_value_extraction(false)
    {
      //      if (copy_initialized && other.optimized() == true)
      //      {
      //        dependent_variables_output = other.dependent_variables_output;
      //        ready_for_value_extraction = other.ready_for_value_extraction;
      //        other.clone_optimizer(optimizer);
      //      }
    }


    template <typename ReturnType>
    void
    BatchOptimizer<ReturnType>::set_optimization_method(
      const enum OptimizerType &    optimization_method,
      const enum OptimizationFlags &optimization_flags)
    {
      Assert(
        optimized() == false,
        ExcMessage(
          "Cannot call set_optimization_method() once the optimizer is finalised."));
      //      Assert (this->optimization_method() == OptimizerType::off,
      //              ExcMessage("Cannot call set_optimization_method() more
      //              than once."));

      method = optimization_method;
#  ifndef HAVE_SYMENGINE_LLVM
      if (this->optimization_method() == OptimizerType::llvm)
        {
          // Fall-back if the LLVM JIT compiler is not available
          deallog
            << "Warning: The LLVM is not available, so the batch optimizer "
            << "is using a lambda optimizer instead." << std::endl;
          method = OptimizerType::lambda;
        }
#  endif
      flags = optimization_flags;
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
      //      return (static_cast<int>(method != OptimizerType::llvm) &&
      //      static_cast<int>(flags & OptimizationFlags::optimize_cse));
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
      const SD::types::substitution_map &symbol_values)
    {
      Assert(optimized() == false,
             ExcMessage(
               "Cannot register symbols once the optimizer is finalised."));

#  ifdef DEBUG
      // Ensure that all of the keys in the map are actually symbolic
      // in nature
      for (SD::types::substitution_map::const_iterator it =
             symbol_values.begin();
           it != symbol_values.end();
           ++it)
        {
          Assert(SymEngine::is_a<SymEngine::Symbol>(*(it->first.get_RCP())),
                 ExcMessage("Key entry in map is not a symbol"));
        }
#  endif
      // Merge the two maps, in the process ensuring that there is no
      // duplication of symbols
      independent_variables_symbols.insert(symbol_values.begin(),
                                           symbol_values.end());
    }


    template <typename ReturnType>
    void
    BatchOptimizer<ReturnType>::register_symbols(
      const SymEngine::map_basic_basic &symbol_values)
    {
      register_symbols(SD::Utilities::convert_basic_map_to_expression_map(symbol_values));
    }


    template <typename ReturnType>
    void
    BatchOptimizer<ReturnType>::register_symbols(
      const SD::types::symbol_vector &symbols)
    {
      Assert(optimized() == false,
             ExcMessage(
               "Cannot register symbols once the optimizer is finalised."));

      for (typename SD::types::symbol_vector::const_iterator it =
             symbols.begin();
           it != symbols.end();
           ++it)
        {
          Assert(independent_variables_symbols.find(*it) ==
                   independent_variables_symbols.end(),
                 ExcMessage("Symbol is already in the map."));
          independent_variables_symbols.insert(
            std::make_pair(*it, SD::Expression(0.0)));
        }
    }


    template <typename ReturnType>
    void
    BatchOptimizer<ReturnType>::register_symbols(
      const SymEngine::vec_basic &symbols)
    {
      register_symbols(SD::Utilities::convert_basic_vector_to_expression_vector(symbols));
    }


    template <typename ReturnType>
    SD::types::symbol_vector
    BatchOptimizer<ReturnType>::get_independent_symbols(void) const
    {
      //      SD::types::symbol_vector independent_symbols;
      //      independent_symbols.reserve(independent_variables_symbols.size());
      //
      //      for (SD::types::substitution_map::const_iterator
      //          it = independent_variables_symbols.begin();
      //          it != independent_variables_symbols.end(); ++it)
      //        {
      //          independent_symbols.emplace_back(it->first);
      //        }
      //
      //      return independent_symbols;
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
    BatchOptimizer<ReturnType>::register_function(const Expression &func)
    {
      Assert(optimized() == false,
             ExcMessage(
               "Cannot register functions once the optimizer is finalised."));

      register_scalar_function(func.get_RCP());
    }


    template <typename ReturnType>
    void
    BatchOptimizer<ReturnType>::register_functions(
      const SD::types::symbol_vector &functions)
    {
      Assert(optimized() == false,
             ExcMessage(
               "Cannot register functions once the optimizer is finalised."));

      register_vector_functions(functions);
    }


    template <typename ReturnType>
    void
    BatchOptimizer<ReturnType>::register_functions(
      const SymEngine::vec_basic &functions)
    {
      register_functions(Utilities::convert_basic_vector_to_expression_vector(functions));
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
      Assert(map_dep_ptr_basic_vec_entry.size() ==
               dependent_variables_functions.size(),
             ExcInternalError());
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
          Assert(optimization_method() == OptimizerType::off,
                 ExcInternalError());
          //          Assert(use_symbolic_CSE() == false, ExcInternalError());
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
          //          Assert(use_symbolic_CSE() == false, ExcInternalError());
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
      //      // Perform CSE if necessary
      //      // Note: The LLVM JIT compiler does its own expression
      //      optimisation,
      //      //       so we don't want to do anything to get in the way there
      //      else if (typename internal::CSEOptimizer<ReturnType,
      //      internal::DictionaryOptimizer>::OptimizerType*
      //               opt = dynamic_cast<typename
      //               internal::CSEOptimizer<ReturnType,
      //               internal::DictionaryOptimizer>::OptimizerType *>
      //               (optimizer.get()))
      //        {
      //          Assert(optimization_method() == OptimizerType::off,
      //          ExcInternalError()); Assert(use_symbolic_CSE() == true,
      //          ExcInternalError()); internal::OptimizerHelper<ReturnType,
      //          typename internal::CSEOptimizer<ReturnType,
      //          internal::DictionaryOptimizer> >::initialize(opt, symbol_vec,
      //          dependent_variables_functions);
      //        }
      //      else if (typename internal::CSEOptimizer<ReturnType,
      //      internal::LambdaOptimizer>::OptimizerType*
      //               opt = dynamic_cast<typename
      //               internal::CSEOptimizer<ReturnType,
      //               internal::LambdaOptimizer>::OptimizerType *>
      //               (optimizer.get()))
      //        {
      //          Assert(optimization_method() == OptimizerType::lambda,
      //          ExcInternalError()); Assert(use_symbolic_CSE() == true,
      //          ExcInternalError()); internal::OptimizerHelper<ReturnType,
      //          typename internal::CSEOptimizer<ReturnType,
      //          internal::LambdaOptimizer> >::initialize(opt, symbol_vec,
      //          dependent_variables_functions);
      //        }
      else
        {
          AssertThrow(false, ExcMessage("Unknown optimizer type."));
        }

      // The size of the outputs is now fixed, as is the number and
      // order of the symbols to be substituted.
      // Note: When no optimisation is actually used (i.e. optimization_method()
      // == off), we could conceptually go without this data structure. However,
      // since the user expects to perform substitution of all dependent
      // variables in one go, we still require it for intermediate
      // storage of results.
      dependent_variables_output.resize(n_dependent_variables());
    }


    template <typename ReturnType>
    void
    BatchOptimizer<ReturnType>::substitute(
      const SD::types::substitution_map &substitution_values) const
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
        Utilities::extract_symbols(substitution_values);
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
        Utilities::extract_values<ReturnType>(substitution_values);
      substitute(values);
    }


    template <typename ReturnType>
    void
    BatchOptimizer<ReturnType>::substitute(
      const SymEngine::map_basic_basic &substitution_values) const
    {
      substitute(SD::Utilities::convert_basic_map_to_expression_map(substitution_values));
    }


    template <typename ReturnType>
    void
    BatchOptimizer<ReturnType>::substitute(
      const SD::types::symbol_vector &symbols,
      const std::vector<ReturnType> &          substitution_values) const
    {
      // Zip the two vectors and use the other function call
      // This ensures the ordering of the input vectors matches that of the
      // stored map.
      substitute(make_substitution_map(symbols, substitution_values));
    }


    template <typename ReturnType>
    void
    BatchOptimizer<ReturnType>::substitute(
      const SymEngine::vec_basic &symbols,
      const std::vector<ReturnType> &          substitution_values) const
    {
      substitute(SD::Utilities::convert_basic_vector_to_expression_vector(symbols), substitution_values);
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
          Assert(optimization_method() == OptimizerType::off,
                 ExcInternalError());
          //          Assert(use_symbolic_CSE() == false, ExcInternalError());
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
          //          Assert(use_symbolic_CSE() == false, ExcInternalError());
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
      //      else if (typename internal::CSEOptimizer<ReturnType,
      //      internal::DictionaryOptimizer>::OptimizerType*
      //               opt = dynamic_cast<typename
      //               internal::CSEOptimizer<ReturnType,
      //               internal::DictionaryOptimizer>::OptimizerType *>
      //               (optimizer.get()))
      //        {
      //          Assert(optimization_method() == OptimizerType::off,
      //          ExcInternalError()); Assert(use_symbolic_CSE() == true,
      //          ExcInternalError()); internal::OptimizerHelper<ReturnType,
      //          typename internal::CSEOptimizer<ReturnType,
      //          internal::DictionaryOptimizer> >::substitute(opt,
      //          dependent_variables_output, substitution_values);
      //        }
      //      else if (typename internal::CSEOptimizer<ReturnType,
      //      internal::LambdaOptimizer>::OptimizerType*
      //               opt = dynamic_cast<typename
      //               internal::CSEOptimizer<ReturnType,
      //               internal::LambdaOptimizer>::OptimizerType *>
      //               (optimizer.get()))
      //        {
      //          Assert(optimization_method() == OptimizerType::lambda,
      //          ExcInternalError()); Assert(use_symbolic_CSE() == true,
      //          ExcInternalError()); internal::OptimizerHelper<ReturnType,
      //          typename internal::CSEOptimizer<ReturnType,
      //          internal::LambdaOptimizer> >::substitute(opt,
      //          dependent_variables_output, substitution_values);
      //        }
      else
        {
          AssertThrow(false, ExcNotImplemented());
        }

      ready_for_value_extraction = true;
    }


    template <typename ReturnType>
    const std::vector<ReturnType> &
    BatchOptimizer<ReturnType>::get_evaluated_functions() const
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

      const typename map_dependent_expression_to_vector_entry_t::const_iterator
        it = map_dep_ptr_basic_vec_entry.find(func.get_RCP());
      Assert(it != map_dep_ptr_basic_vec_entry.end(),
             ExcMessage("Function has not been registered."));
      Assert(it->second < n_dependent_variables(), ExcInternalError());

      return dependent_variables_output[it->second];
    }


#  ifdef DEBUG
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
#  endif


    template <typename ReturnType>
    void
    BatchOptimizer<ReturnType>::register_scalar_function(
      const SymEngine::RCP<const SymEngine::Basic> &func)
    {
      Assert(
        dependent_variables_output.size() == 0,
        ExcMessage(
          "Cannot register function as the optimizer has already been finalized."));
      dependent_variables_output.reserve(n_dependent_variables() + 1);
      const bool entry_registered = (map_dep_ptr_basic_vec_entry.find(func) !=
                                     map_dep_ptr_basic_vec_entry.end());
#  ifdef DEBUG
      if (entry_registered == true &&
          is_valid_nonunique_dependent_variable(func) == false)
        Assert(entry_registered,
               ExcMessage("Function has already been registered."));
#  endif
      if (entry_registered == false)
        {
          dependent_variables_functions.push_back(func);
          map_dep_ptr_basic_vec_entry[func] =
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

      for (typename SD::types::symbol_vector::const_iterator it = funcs.begin();
           it != funcs.end();
           ++it)
        {
          const SymEngine::RCP<const SymEngine::Basic> &func = *it;
          const bool                                    entry_registered =
            (map_dep_ptr_basic_vec_entry.find(func) !=
             map_dep_ptr_basic_vec_entry.end());
#  ifdef DEBUG
          if (entry_registered == true &&
              is_valid_nonunique_dependent_variable(func) == false)
            Assert(entry_registered,
                   ExcMessage("Function has already been registered."));
#  endif
          if (entry_registered == false)
            {
              dependent_variables_functions.push_back(func);
              map_dep_ptr_basic_vec_entry[func] =
                dependent_variables_functions.size() - 1;
            }
        }
    }

    //    namespace internal
    //    {
    //      namespace
    //      {
    //        template<typename NumberType>
    //        struct CSEValueConversion
    //        {
    //          static inline NumberType
    //          basic_to_value (const SymEngine::RCP<const SymEngine::Basic> &x)
    //          {
    //            return internal::FloatingPointNumber<NumberType>::evaluate(x);
    //          };
    //
    //          static inline
    //          SymEngine::RCP<const SymEngine::Basic> value_to_basic(const
    //          NumberType &x)
    //          {
    //            return internal::make_symengine_rcp(x);
    //          };
    //        };
    //      }
    //    }

    template <typename ReturnType>
    void
    BatchOptimizer<ReturnType>::create_optimizer(
      std::unique_ptr<SymEngine::Visitor> &optimizer)
    {
      Assert(!optimizer, ExcMessage("Optimizer has already been created."));

      //      if (use_symbolic_CSE() == false)
      //        {
      if (optimization_method() == OptimizerType::off ||
          optimization_method() == OptimizerType::dictionary)
        {
          optimizer.reset(new typename internal::DictionaryOptimizer<
                          ReturnType>::OptimizerType());
        }
      else if (optimization_method() == OptimizerType::lambda)
        {
          optimizer.reset(
            new
            typename internal::LambdaOptimizer<ReturnType>::OptimizerType());
        }
      else if (optimization_method() == OptimizerType::llvm)
        {
#  ifdef HAVE_SYMENGINE_LLVM
          optimizer.reset(
            new typename internal::LLVMOptimizer<ReturnType>::OptimizerType());
#  else
          AssertThrow(false, ExcMessage("The LLVM compiler is not available."));
#  endif
        }
      else
        {
          AssertThrow(false, ExcMessage("Unknown optimizer selected."));
          //              AssertThrow(false, ExcMessage("Unknown non-CSE
          //              optimizer selected."));
        }
      //        }
      //      else
      //        {
      //          AssertThrow(false, ExcMessage("Unknown CSE optimizer
      //          selected."));
      //        }
      //      else
      //        {
      //          if (optimization_method() == OptimizerType::off)
      //            {
      //              typedef typename internal::CSEOptimizer<ReturnType,
      //              internal::DictionaryOptimizer>::OptimizerType
      //              CSEOptimizerType; typedef typename
      //              internal::CSEOptimizer<ReturnType,
      //              internal::DictionaryOptimizer>::ReturnType CSEReturnType;
      //
      //              optimizer.reset(new
      //              CSEOptimizerType(internal::CSEValueConversion<CSEReturnType>::basic_to_value,
      //                                                   internal::CSEValueConversion<CSEReturnType>::value_to_basic));
      //            }
      //          else if (optimization_method() == OptimizerType::lambda)
      //            {
      //              typedef typename internal::CSEOptimizer<ReturnType,
      //              internal::LambdaOptimizer>::OptimizerType
      //              CSEOptimizerType; typedef typename
      //              internal::CSEOptimizer<ReturnType,
      //              internal::LambdaOptimizer>::ReturnType    CSEReturnType;
      //
      //              optimizer.reset(new
      //              CSEOptimizerType(internal::CSEValueConversion<CSEReturnType>::basic_to_value,
      //                                                   internal::CSEValueConversion<CSEReturnType>::value_to_basic));
      //            }
      //          else
      //            {
      //              AssertThrow(false, ExcMessage("Unknown CSE optimizer
      //              selected."));
      //            }
      //        }
    }

    template <typename ReturnType>
    void
    BatchOptimizer<ReturnType>::clone_optimizer(
      std::unique_ptr<SymEngine::Visitor> & /*new_optimizer*/) const
    {
      // Shared memory pool for AD numbers cannot be copied...
      AssertThrow(false, ExcNotImplemented());

      //      Assert(optimizer, ExcNotInitialized());
      //
      ////      if (use_symbolic_CSE() == false)
      ////        {
      //      if (optimization_method() == OptimizerType::off ||
      //          optimization_method() == OptimizerType::dictionary)
      //        {
      //          typedef typename
      //          internal::DictionaryOptimizer<ReturnType>::OptimizerType
      //          OptType; new_optimizer.reset(new OptType(dynamic_cast<const
      //          OptType &>(*optimizer)));
      //        }
      //      else if (optimization_method() == OptimizerType::lambda)
      //        {
      //          typedef typename
      //          internal::LambdaOptimizer<ReturnType>::OptimizerType OptType;
      //          new_optimizer.reset(new OptType(dynamic_cast<const OptType
      //          &>(*optimizer)));
      //        }
      //      else if (optimization_method() == OptimizerType::llvm)
      //        {
      //#ifdef HAVE_SYMENGINE_LLVM
      //          typedef typename
      //          internal::LLVMOptimizer<ReturnType>::OptimizerType OptType;
      //          new_optimizer.reset(new OptType(dynamic_cast<const OptType
      //          &>(*optimizer)));
      //#else
      //          AssertThrow(false, ExcMessage("The LLVM compiler is not
      //          available."));
      //#endif
      //        }
      //      else
      //        {
      //          AssertThrow(false, ExcMessage("Unknown optimizer selected."));
      //        }
    }

  } // namespace SD
} // namespace Differentiation


/* --- Explicit instantiations --- */
#  include "symengine_optimizer.inst"


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SYMENGINE
