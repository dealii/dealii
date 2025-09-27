// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_differentiation_sd_symengine_optimizer_h
#define dealii_differentiation_sd_symengine_optimizer_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SYMENGINE

// Low level
#  include <symengine/basic.h>
#  include <symengine/dict.h>
#  include <symengine/symengine_exception.h>
#  include <symengine/symengine_rcp.h>

// Optimization
#  include <symengine/lambda_double.h>
#  include <symengine/visitor.h>
#  ifdef HAVE_SYMENGINE_LLVM
#    include <symengine/llvm_double.h>
#  endif

#  include <deal.II/base/exceptions.h>
#  include <deal.II/base/utilities.h>

#  include <deal.II/differentiation/sd/symengine_number_types.h>
#  include <deal.II/differentiation/sd/symengine_number_visitor_internal.h>
#  include <deal.II/differentiation/sd/symengine_scalar_operations.h>
#  include <deal.II/differentiation/sd/symengine_tensor_operations.h>
#  include <deal.II/differentiation/sd/symengine_utilities.h>

#  include <boost/serialization/split_member.hpp>
#  include <boost/type_traits.hpp>

#  include <algorithm>
#  include <map>
#  include <memory>
#  include <type_traits>
#  include <utility>
#  include <vector>


DEAL_II_NAMESPACE_OPEN


namespace Differentiation
{
  namespace SD
  {
    /**
     * @addtogroup Exceptions
     * @{
     */

    /**
     * An exception to indicate that the SymEngine library has not been built
     * against the LLVM compiler.
     */
    DeclExceptionMsg(ExcSymEngineLLVMNotAvailable,
                     "SymEngine has not been built with LLVM support.");

    /**
     * An exception to indicate that SymEngine's LLVM optimizer doesn't work
     * with the desired return type.
     */
    DeclExceptionMsg(ExcSymEngineLLVMReturnTypeNotSupported,
                     "The SymEngine LLVM optimizer does not (yet) support the "
                     "selected return type.");

    /** @} */


    // Forward declarations
    template <typename ReturnType>
    class BatchOptimizer;


    /**
     * An enumeration to distinguish between different optimization methods
     * that can be used by SymEngine to more rapidly evaluate complex
     * symbolic expressions.
     */
    enum class OptimizerType
    {
      /**
       * Use dictionary substitution. This is SymEngine's default method.
       */
      dictionary,
      /**
       * Convert the symbolic expression into a collection of
       * `std::function`s.
       */
      lambda,
      /**
       * Use the LLVM JIT compiler to compile the expression into an
       * aggressively optimized, stand-alone function.
       */
      llvm
    };


    /**
     * Output operator that outputs the selected optimizer type.
     */
    template <typename StreamType>
    inline StreamType &
    operator<<(StreamType &s, OptimizerType o)
    {
      if (o == OptimizerType::dictionary)
        s << "dictionary";
      else if (o == OptimizerType::lambda)
        s << "lambda";
      else if (o == OptimizerType::llvm)
        s << "llvm";
      else
        {
          Assert(false, ExcMessage("Unknown optimization method."));
        }

      return s;
    }


    /**
     * An enumeration to specify which special techniques, over and above
     * those used with the chosen OptimizerType, to be applied to the
     * set of expressions that are to be optimized.
     */
    enum class OptimizationFlags : unsigned char
    {
      /**
       * No additional optimization.
       */
      optimize_default = 0,
      /**
       * Apply common subexpresson elimination.
       */
      optimize_cse = 0x0001,
      /**
       * Employ aggressive optimizations when compiling with the LLVM JIT
       * compiler.
       */
      optimize_aggressive = 0x0002,
      /**
       * Apply all possible optimizations.
       */
      optimize_all = optimize_cse | optimize_aggressive
    };


    /**
     * A global operator that returns an object in which all bits are
     * individually set in the following way:
     * If the corresponding bit in either the first or second argument are set,
     * then the output bit it set. Otherwise the output bit remains unset.
     * This `or` type operation is performed for each bit composing the input
     * arguments (and output) in an individual manner.
     */
    // This operator exists since if it did not then the result of the bit-or
    // <tt>operator |</tt> would be an integer which would in turn trigger a
    // compiler warning when we tried to assign it to an object of type
    // OptimizationFlags.
    inline OptimizationFlags
    operator|(const OptimizationFlags f1, const OptimizationFlags f2)
    {
      return static_cast<OptimizationFlags>(static_cast<unsigned int>(f1) |
                                            static_cast<unsigned int>(f2));
    }


    /**
     * Global operator that sets the bits from the second argument also in the
     * first one.
     */
    inline OptimizationFlags &
    operator|=(OptimizationFlags &f1, const OptimizationFlags f2)
    {
      f1 = f1 | f2;
      return f1;
    }


    /**
     * A global operator that returns an object in which all bits are
     * individually set in the following way:
     * If the corresponding bit in both the first or second argument are set,
     * then the output bit it set. Otherwise the output bit remains unset.
     * This `and` type operation is performed for each bit composing the input
     * arguments (and output) in an individual manner.
     */
    // This operator exists since if it did not then the result of the bit-or
    // <tt>operator |</tt> would be an integer which would in turn trigger a
    // compiler warning when we tried to assign it to an object of type
    // OptimizationFlags.
    inline OptimizationFlags
    operator&(const OptimizationFlags f1, const OptimizationFlags f2)
    {
      return static_cast<OptimizationFlags>(static_cast<unsigned int>(f1) &
                                            static_cast<unsigned int>(f2));
    }


    /**
     * Global operator which clears all the bits in the first argument if they
     * are not also set in the second argument.
     */
    inline OptimizationFlags &
    operator&=(OptimizationFlags &f1, const OptimizationFlags f2)
    {
      f1 = f1 & f2;
      return f1;
    }


    namespace internal
    {
      /**
       * A utility function that checks whether or not CSE
       * has been selected as an optimization flag.
       */
      inline bool
      use_symbolic_CSE(const enum OptimizationFlags &flags)
      {
        return static_cast<int>(flags & OptimizationFlags::optimize_cse);
      }

      /**
       * A utility function that returns the optimization level
       * that is to be employed when the LLVM optimizer is invoked.
       */
      inline int
      get_LLVM_optimization_level(const enum OptimizationFlags &flags)
      {
        // With the LLVM compiler there exists the opportunity to tune
        // the level of optimizations performed during compilation.
        // By default SymEngine sets this at "opt_level=2", which one
        // presumes targets -O2. Here we are a bit more specific about
        // want we want it to do:
        // - Normal compilation: -02 (default settings)
        // - Aggressive mode: -03 (the whole lot!)
        // In theory we could also target
        // - Debug mode: -O0 (no optimizations)
        // but this doesn't make much sense since SymEngine is a
        // tested external library.
        const bool use_agg_opt =
          static_cast<int>(flags & OptimizationFlags::optimize_aggressive);
        const int opt_level = (use_agg_opt ? 3 : 2);
        return opt_level;
      }
    } // namespace internal


    /**
     * Output operator that outputs optimization flags as a set of or'd
     * text values.
     */
    template <typename StreamType>
    inline StreamType &
    operator<<(StreamType &s, OptimizationFlags o)
    {
      s << " OptimizationFlags|";
      if (static_cast<unsigned int>(o & OptimizationFlags::optimize_cse))
        s << "cse|";

      // LLVM optimization level
      s << "-O" + std::to_string(internal::get_LLVM_optimization_level(o)) +
             "|";

      return s;
    }


    namespace internal
    {
      /**
       * A wrapper for dictionary based optimization.
       *
       * @tparam ReturnType The number type that is returned as a result
       *         of operations performed by the optimizer.
       *         Floating point and complex numbers are currently supported.
       * @tparam T An arbitrary type resulting from the application of
       *         the SFINAE idiom to selectively specialize this class.
       */
      template <typename ReturnType, typename T = void>
      struct DictionaryOptimizer;


      /**
       * A wrapper for SymEngine's "lambda" optimizer.
       *
       * @tparam ReturnType The number type that is returned as a result
       *         of operations performed by the optimizer.
       *         Floating point and complex numbers are currently supported.
       * @tparam T An arbitrary type resulting from the application of
       *         the SFINAE idiom to selectively specialize this class.
       */
      template <typename ReturnType, typename T = void>
      struct LambdaOptimizer;


#  ifdef HAVE_SYMENGINE_LLVM
      /**
       * A wrapper for SymEngine's LLVM JIT optimizer.
       *
       * @tparam ReturnType The number type that is returned as a result
       *         of operations performed by the optimizer.
       *         Floating point and complex numbers are currently supported.
       * @tparam T An arbitrary type resulting from the application of
       *         the SFINAE idiom to selectively specialize this class.
       */
      template <typename ReturnType, typename T = void>
      struct LLVMOptimizer;
#  endif // HAVE_SYMENGINE_LLVM


      /**
       * A wrapper class for all supported Optimizer types and
       * @p ReturnTypes. It aims to deal with the case when the
       * @p ReturnType and native return type of the @p Optimizer
       * are not the same.
       *
       * @tparam ReturnType The number type that is returned as a result
       *         of operations performed by the optimizer.
       *         Floating point and complex numbers are currently supported.
       * @tparam Optimizer An internal class that implements a wrapper to a
       *         SymEngine optimizer. Currently, the target classes are the
       *         DictionaryOptimizer, the LambdaOptimizer and the LLVMOptimizer.
       * @tparam T An arbitrary type resulting from the application of
       *         the SFINAE idiom to selectively specialize this class.
       */
      template <typename ReturnType, typename Optimizer, typename T = void>
      struct OptimizerHelper;


#  ifndef DOXYGEN


      /* ----------- Specializations for the Optimizers ----------- */


      // A helper struct to type trait detection for the optimizers that
      // will be defined next.
      template <typename ReturnType_, typename T = void>
      struct SupportedOptimizerTypeTraits
      {
        static const bool is_supported = false;

        using ReturnType = void;
      };



      // Specialization for arithmetic types
      template <typename ReturnType_>
      struct SupportedOptimizerTypeTraits<
        ReturnType_,
        std::enable_if_t<std::is_arithmetic_v<ReturnType_>>>
      {
        static const bool is_supported = true;

        using ReturnType =
          std::conditional_t<std::is_same_v<ReturnType_, float>, float, double>;
      };



      // Specialization for complex arithmetic types
      template <typename ReturnType_>
      struct SupportedOptimizerTypeTraits<
        ReturnType_,
        std::enable_if_t<
          boost::is_complex<ReturnType_>::value &&
          std::is_arithmetic_v<typename ReturnType_::value_type>>>
      {
        static const bool is_supported = true;

        using ReturnType =
          std::conditional_t<std::is_same_v<ReturnType_, std::complex<float>>,
                             std::complex<float>,
                             std::complex<double>>;
      };



      template <typename ReturnType_>
      struct DictionaryOptimizer<ReturnType_,
                                 std::enable_if_t<SupportedOptimizerTypeTraits<
                                   ReturnType_>::is_supported>>
      {
        using ReturnType =
          typename SupportedOptimizerTypeTraits<ReturnType_>::ReturnType;
        using OptimizerType =
          internal::DictionarySubstitutionVisitor<ReturnType, SD::Expression>;


        /**
         * Initialize an instance of an optimizer.
         *
         * @param optimizer The optimizer to be initialized.
         * @param independent_symbols A vector of symbols that represent independent variables.
         * @param dependent_functions A vector of expressions that represent dependent variables.
         * @param optimization_flags A set of flags that indicate the types of optimization to be performed.
         */
        static void
        initialize(OptimizerType                &optimizer,
                   const SymEngine::vec_basic   &independent_symbols,
                   const SymEngine::vec_basic   &dependent_functions,
                   const enum OptimizationFlags &optimization_flags)
        {
          const bool use_symbolic_cse = use_symbolic_CSE(optimization_flags);
          optimizer.init(independent_symbols,
                         dependent_functions,
                         use_symbolic_cse);
        }



        /**
         * Write the data of the @p optimizer to a stream for the purpose
         * of serialization.
         */
        template <class Archive>
        static void
        save(Archive           &archive,
             const unsigned int version,
             OptimizerType     &optimizer)
        {
          optimizer.save(archive, version);
        }



        /**
         * Read the data for the @p optimizer from a stream for the purpose
         * of serialization.
         */
        template <class Archive>
        static void
        load(Archive           &archive,
             const unsigned int version,
             OptimizerType     &optimizer,
             const SymEngine::vec_basic & /*independent_symbols*/,
             const SymEngine::vec_basic & /*dependent_functions*/,
             const enum OptimizationFlags & /*optimization_flags*/)
        {
          optimizer.load(archive, version);
        }



        /**
         * Print some information on state of the internal data
         * structures stored in the @p optimizer.
         *
         * @tparam Stream The type for the output stream.
         * @param stream The output stream to print to.
         * @param optimizer The instance of the optimizer from which to retrieve
         * information to print to the stream.
         * @param print_independent_symbols A flag to indicate if the independent
         * variables should be outputted to the @p stream.
         * @param print_dependent_functions A flag to indicate if the dependent
         * expressions should be outputted to the @p stream.
         * @param print_cse_reductions A flag to indicate whether or not all
         * common subexpressions should be printed to the @p stream.
         */
        template <typename Stream>
        static void
        print(Stream              &stream,
              const OptimizerType &optimizer,
              const bool           print_independent_symbols = false,
              const bool           print_dependent_functions = false,
              const bool           print_cse_reductions      = true)
        {
          optimizer.print(stream,
                          print_independent_symbols,
                          print_dependent_functions,
                          print_cse_reductions);
        }
      };



      template <typename ReturnType_>
      struct LambdaOptimizer<ReturnType_,
                             std::enable_if_t<SupportedOptimizerTypeTraits<
                               ReturnType_>::is_supported>>
      {
        using ReturnType =
          std::conditional_t<!boost::is_complex<ReturnType_>::value,
                             double,
                             std::complex<double>>;
        using OptimizerType =
          std::conditional_t<!boost::is_complex<ReturnType_>::value,
                             SymEngine::LambdaRealDoubleVisitor,
                             SymEngine::LambdaComplexDoubleVisitor>;


        /**
         * Initialize an instance of an optimizer.
         *
         * @param optimizer The optimizer to be initialized.
         * @param independent_symbols A vector of symbols that represent independent variables.
         * @param dependent_functions A vector of expressions that represent dependent variables.
         * @param optimization_flags A set of flags that indicate the types of optimization to be performed.
         */
        static void
        initialize(OptimizerType                &optimizer,
                   const SymEngine::vec_basic   &independent_symbols,
                   const SymEngine::vec_basic   &dependent_functions,
                   const enum OptimizationFlags &optimization_flags)
        {
          const bool use_symbolic_cse = use_symbolic_CSE(optimization_flags);
          optimizer.init(independent_symbols,
                         dependent_functions,
                         use_symbolic_cse);
        }



        /**
         * Write the data of the @p optimizer to a stream for the purpose
         * of serialization.
         */
        template <class Archive>
        static void
        save(Archive & /*archive*/,
             const unsigned int /*version*/,
             OptimizerType & /*optimizer*/)
        {}


        /**
         * Read the data for the @p optimizer from a stream for the purpose
         * of serialization.
         */
        template <class Archive>
        static void
        load(Archive & /*archive*/,
             const unsigned int /*version*/,
             OptimizerType                &optimizer,
             const SymEngine::vec_basic   &independent_symbols,
             const SymEngine::vec_basic   &dependent_functions,
             const enum OptimizationFlags &optimization_flags)
        {
          initialize(optimizer,
                     independent_symbols,
                     dependent_functions,
                     optimization_flags);
        }



        /**
         * Print some information on state of the internal data
         * structures stored in the @p optimizer.
         *
         * @tparam Stream The type for the output stream.
         * @param stream The output stream to print to.
         * @param optimizer The instance of the optimizer from which to retrieve
         * information to print to the stream.
         * @param print_independent_symbols A flag to indicate if the independent
         * variables should be outputted to the @p stream.
         * @param print_dependent_functions A flag to indicate if the dependent
         * expressions should be outputted to the @p stream.
         * @param print_cse_reductions A flag to indicate whether or not all
         * common subexpressions should be printed to the @p stream.
         */
        template <typename StreamType>
        static void
        print(StreamType & /*stream*/,
              const OptimizerType & /*optimizer*/,
              const bool /*print_independent_symbols*/ = false,
              const bool /*print_dependent_functions*/ = false,
              const bool /*print_cse_reductions*/      = true)
        {
          // No built-in print function
        }
      };



#    ifdef HAVE_SYMENGINE_LLVM
      template <typename ReturnType_>
      struct LLVMOptimizer<ReturnType_,
                           std::enable_if_t<std::is_arithmetic_v<ReturnType_>>>
      {
        using ReturnType =
          std::conditional_t<std::is_same_v<ReturnType_, float>, float, double>;
        using OptimizerType =
          std::conditional_t<std::is_same_v<ReturnType_, float>,
                             SymEngine::LLVMFloatVisitor,
                             SymEngine::LLVMDoubleVisitor>;

        /**
         * A flag to indicate if the ReturnType is supported by a
         * SymEngine LLVM wrapper.
         */
        static const bool supported_by_LLVM = true;


        /**
         * Initialize an instance of an optimizer.
         *
         * @param optimizer The optimizer to be initialized.
         * @param independent_symbols A vector of symbols that represent independent variables.
         * @param dependent_functions A vector of expressions that represent dependent variables.
         * @param optimization_flags A set of flags that indicate the types of optimization to be performed.
         */
        static void
        initialize(OptimizerType                &optimizer,
                   const SymEngine::vec_basic   &independent_symbols,
                   const SymEngine::vec_basic   &dependent_functions,
                   const enum OptimizationFlags &optimization_flags)
        {
          const int opt_level = get_LLVM_optimization_level(optimization_flags);
          const bool use_symbolic_cse = use_symbolic_CSE(optimization_flags);
          optimizer.init(independent_symbols,
                         dependent_functions,
                         use_symbolic_cse,
                         opt_level);
        }



        /**
         * Write the data of the @p optimizer to a stream for the purpose
         * of serialization.
         */
        template <class Archive>
        static void
        save(Archive &archive,
             const unsigned int /*version*/,
             OptimizerType &optimizer)
        {
          const std::string llvm_compiled_function = optimizer.dumps();
          archive          &llvm_compiled_function;
        }



        /**
         * Read the data for the @p optimizer from a stream for the purpose
         * of serialization.
         */
        template <class Archive>
        static void
        load(Archive &archive,
             const unsigned int /*version*/,
             OptimizerType &optimizer,
             const SymEngine::vec_basic & /*independent_symbols*/,
             const SymEngine::vec_basic & /*dependent_functions*/,
             const enum OptimizationFlags & /*optimization_flags*/)
        {
          std::string llvm_compiled_function;
          archive    &llvm_compiled_function;
          optimizer.loads(llvm_compiled_function);
        }



        /**
         * Print some information on state of the internal data
         * structures stored in the @p optimizer.
         *
         * @tparam Stream The type for the output stream.
         * @param stream The output stream to print to.
         * @param optimizer The instance of the optimizer from which to retrieve
         * information to print to the stream.
         * @param print_independent_symbols A flag to indicate if the independent
         * variables should be outputted to the @p stream.
         * @param print_dependent_functions A flag to indicate if the dependent
         * expressions should be outputted to the @p stream.
         * @param print_cse_reductions A flag to indicate whether or not all
         * common subexpressions should be printed to the @p stream.
         */
        template <typename StreamType>
        static void
        print(StreamType & /*stream*/,
              const OptimizerType & /*optimizer*/,
              const bool /*print_independent_symbols*/ = false,
              const bool /*print_dependent_functions*/ = false,
              const bool /*print_cse_reductions*/      = true)
        {
          // No built-in print function
        }
      };


      // There is no LLVM optimizer built with complex number support.
      // So we fall back to the LambdaDouble case as a type (required
      // at compile time), but offer no implementation. We expect that
      // the calling class does not create this type: This can be done by
      // checking the `supported_by_LLVM` flag.
      template <typename ReturnType_>
      struct LLVMOptimizer<
        ReturnType_,
        std::enable_if_t<
          boost::is_complex<ReturnType_>::value &&
          std::is_arithmetic_v<typename ReturnType_::value_type>>>
      {
        // Since there is no working implementation, these are dummy types
        // that help with templating in the calling function.
        using ReturnType = typename LambdaOptimizer<ReturnType_>::ReturnType;
        using OptimizerType =
          typename LambdaOptimizer<ReturnType_>::OptimizerType;

        /**
         * A flag to indicate if the ReturnType is supported by a
         * SymEngine LLVM wrapper.
         */
        static const bool supported_by_LLVM = false;


        /**
         * Initialize an instance of an optimizer.
         *
         * @param optimizer The optimizer to be initialized.
         * @param independent_symbols A vector of symbols that represent independent variables.
         * @param dependent_functions A vector of expressions that represent dependent variables.
         * @param optimization_flags A set of flags that indicate the types of optimization to be performed.
         */
        static void
        initialize(OptimizerType & /*optimizer*/,
                   const SymEngine::vec_basic & /*independent_symbols*/,
                   const SymEngine::vec_basic & /*dependent_functions*/,
                   const enum OptimizationFlags & /*optimization_flags*/)
        {
          AssertThrow(false, ExcNotImplemented());
        }



        /**
         * Write the data of the @p optimizer to a stream for the purpose
         * of serialization.
         */
        template <class Archive>
        static void
        save(Archive & /*archive*/,
             const unsigned int /*version*/,
             OptimizerType & /*optimizer*/)
        {
          AssertThrow(false, ExcNotImplemented());
        }



        /**
         * Read the data for the @p optimizer from a stream for the purpose
         * of serialization.
         */
        template <class Archive>
        static void
        load(Archive & /*archive*/,
             const unsigned int /*version*/,
             OptimizerType & /*optimizer*/,
             const SymEngine::vec_basic & /*independent_symbols*/,
             const SymEngine::vec_basic & /*dependent_functions*/,
             const enum OptimizationFlags & /*optimization_flags*/)
        {
          AssertThrow(false, ExcNotImplemented());
        }



        /**
         * Print some information on state of the internal data
         * structures stored in the @p optimizer.
         *
         * @tparam Stream The type for the output stream.
         * @param stream The output stream to print to.
         * @param optimizer The instance of the optimizer from which to retrieve
         * information to print to the stream.
         * @param print_independent_symbols A flag to indicate if the independent
         * variables should be outputted to the @p stream.
         * @param print_dependent_functions A flag to indicate if the dependent
         * expressions should be outputted to the @p stream.
         * @param print_cse_reductions A flag to indicate whether or not all
         * common subexpressions should be printed to the @p stream.
         */
        template <typename StreamType>
        static void
        print(StreamType & /*stream*/,
              const OptimizerType & /*optimizer*/,
              const bool /*print_independent_symbols*/ = false,
              const bool /*print_dependent_functions*/ = false,
              const bool /*print_cse_reductions*/      = true)
        {
          AssertThrow(false, ExcNotImplemented());
        }
      };
#    endif // HAVE_SYMENGINE_LLVM


      /* ----------- Specializations for OptimizerHelper ----------- */


      template <typename ReturnType, typename Optimizer>
      struct OptimizerHelper<
        ReturnType,
        Optimizer,
        std::enable_if_t<
          std::is_same_v<ReturnType, typename Optimizer::ReturnType>>>
      {
        /**
         * Initialize an instance of an optimizer.
         *
         * @param optimizer The optimizer to be initialized.
         * @param independent_symbols A vector of symbols that represent independent variables.
         * @param dependent_functions A vector of expressions that represent dependent variables.
         * @param optimization_flags A set of flags that indicate the types of optimization to be performed.
         */
        static void
        initialize(typename Optimizer::OptimizerType *optimizer,
                   const SymEngine::vec_basic        &independent_symbols,
                   const SymEngine::vec_basic        &dependent_functions,
                   const enum OptimizationFlags      &optimization_flags)
        {
          Assert(optimizer, ExcNotInitialized());

          // Some optimizers don't have the same interface for
          // initialization, we filter them out through the specializations
          // of the Optimizer class
          Optimizer::initialize(*optimizer,
                                independent_symbols,
                                dependent_functions,
                                optimization_flags);
        }



        /**
         * Perform value substitution, evaluating the pre-registered dependent
         * functions with some values associated with all independent symbols.
         *
         * @param optimizer The optimizer on which to perform value substitution
         * for all independent symbols. The values are substituted into the
         * optimized form of the dependent variables that were registered with
         * this class instance.
         * @param output_values The evaluated numerical outcome of the
         * substitution.
         * @param substitution_values The values with which to associate each
         * individual independent symbol.
         */
        static void
        substitute(typename Optimizer::OptimizerType *optimizer,
                   std::vector<ReturnType>           &output_values,
                   const std::vector<ReturnType>     &substitution_values)
        {
          Assert(optimizer, ExcNotInitialized());
          optimizer->call(output_values.data(), substitution_values.data());
        }



        /**
         * Write the data of the @p optimizer to a stream for the purpose
         * of serialization.
         */
        template <class Archive>
        static void
        save(Archive                           &archive,
             const unsigned int                 version,
             typename Optimizer::OptimizerType *optimizer)
        {
          Assert(optimizer, ExcNotInitialized());

          // Some optimizers don't have the same interface for
          // serialization, we filter them out through the specializations
          // of the Optimizer class
          Optimizer::save(archive, version, *optimizer);
        }



        /**
         * Read the data for the @p optimizer from a stream for the purpose
         * of serialization.
         */
        template <class Archive>
        static void
        load(Archive                           &archive,
             const unsigned int                 version,
             typename Optimizer::OptimizerType *optimizer,
             const SymEngine::vec_basic        &independent_symbols,
             const SymEngine::vec_basic        &dependent_functions,
             const enum OptimizationFlags      &optimization_flags)
        {
          Assert(optimizer, ExcNotInitialized());

          // Some optimizers don't have the same interface for
          // serialization, we filter them out through the specializations
          // of the Optimizer class
          Optimizer::load(archive,
                          version,
                          *optimizer,
                          independent_symbols,
                          dependent_functions,
                          optimization_flags);
        }



        /**
         * Print some information on state of the internal data
         * structures stored in the @p optimizer.
         *
         * @tparam Stream The type for the output stream.
         * @param stream The output stream to print to.
         * @param optimizer A pointer to the instance of the optimizer from
         * which to retrieve information to print to the stream.
         * @param print_independent_symbols A flag to indicate if the independent
         * variables should be outputted to the @p stream.
         * @param print_dependent_functions A flag to indicate if the dependent
         * expressions should be outputted to the @p stream.
         * @param print_cse_reductions A flag to indicate whether or not all
         * common subexpressions should be printed to the @p stream.
         */
        template <typename Stream>
        static void
        print(Stream                            &stream,
              typename Optimizer::OptimizerType *optimizer,
              const bool print_independent_symbols = false,
              const bool print_dependent_functions = false,
              const bool print_cse_reductions      = true)
        {
          Assert(optimizer, ExcNotInitialized());

          // Some optimizers don't have a print function, so
          // we filter them out through the specializations of
          // the Optimizer class
          Optimizer::print(stream,
                           *optimizer,
                           print_independent_symbols,
                           print_dependent_functions,
                           print_cse_reductions);
        }
      };

      template <typename ReturnType, typename Optimizer>
      struct OptimizerHelper<
        ReturnType,
        Optimizer,
        std::enable_if_t<
          !std::is_same_v<ReturnType, typename Optimizer::ReturnType>>>
      {
        /**
         * Initialize an instance of an optimizer.
         *
         * @param optimizer The optimizer to be initialized.
         * @param independent_symbols A vector of symbols that represent independent variables.
         * @param dependent_functions A vector of expressions that represent dependent variables.
         * @param optimization_flags A set of flags that indicate the types of optimization to be performed.
         */
        static void
        initialize(typename Optimizer::OptimizerType *optimizer,
                   const SymEngine::vec_basic        &independent_symbols,
                   const SymEngine::vec_basic        &dependent_functions,
                   const enum OptimizationFlags      &optimization_flags)
        {
          Assert(optimizer, ExcNotInitialized());

          const bool use_symbolic_cse = use_symbolic_CSE(optimization_flags);
          optimizer->init(independent_symbols,
                          dependent_functions,
                          use_symbolic_cse);
        }



        /**
         * Perform value substitution, evaluating the pre-registered dependent
         * functions with some values associated with all independent symbols.
         *
         * @param optimizer The optimizer on which to perform value substitution
         * for all independent symbols. The values are substituted into the
         * optimized form of the dependent variables that were registered with
         * this class instance.
         * @param output_values The evaluated numerical outcome of the
         * substitution.
         * @param substitution_values The values with which to associate each
         * individual independent symbol.
         */
        static void
        substitute(typename Optimizer::OptimizerType *optimizer,
                   std::vector<ReturnType>           &output_values,
                   const std::vector<ReturnType>     &substitution_values)
        {
          Assert(optimizer, ExcNotInitialized());

          // Intermediate values to accommodate the difference in
          // value types.
          std::vector<typename Optimizer::ReturnType> int_outputs(
            output_values.size());
          std::vector<typename Optimizer::ReturnType> int_inputs(
            substitution_values.size());

          std::copy(substitution_values.begin(),
                    substitution_values.end(),
                    int_inputs.begin());
          optimizer->call(int_outputs.data(), int_inputs.data());
          std::copy(int_outputs.begin(),
                    int_outputs.end(),
                    output_values.begin());
        }



        /**
         * Write the data of the @p optimizer to a stream for the purpose
         * of serialization.
         */
        template <class Archive>
        static void
        save(Archive                           &archive,
             const unsigned int                 version,
             typename Optimizer::OptimizerType *optimizer)
        {
          Assert(optimizer, ExcNotInitialized());
          Optimizer::save(archive, version, *optimizer);
        }



        /**
         * Read the data for the @p optimizer from a stream for the purpose
         * of serialization.
         */
        template <class Archive>
        static void
        load(Archive                           &archive,
             const unsigned int                 version,
             typename Optimizer::OptimizerType *optimizer,
             const SymEngine::vec_basic        &independent_symbols,
             const SymEngine::vec_basic        &dependent_functions,
             const enum OptimizationFlags      &optimization_flags)
        {
          Assert(optimizer, ExcNotInitialized());

          // Some optimizers don't have the same interface for
          // serialization, we filter them out through the specializations
          // of the Optimizer class
          Optimizer::load(archive,
                          version,
                          *optimizer,
                          independent_symbols,
                          dependent_functions,
                          optimization_flags);
        }



        /**
         * Print some information on state of the internal data
         * structures stored in the @p optimizer.
         *
         * @tparam Stream The type for the output stream.
         * @param stream The output stream to print to.
         * @param optimizer A pointer to the instance of the optimizer from
         * which to retrieve information to print to the stream.
         * @param print_independent_symbols A flag to indicate if the independent
         * variables should be outputted to the @p stream.
         * @param print_dependent_functions A flag to indicate if the dependent
         * expressions should be outputted to the @p stream.
         * @param print_cse_reductions A flag to indicate whether or not all
         * common subexpressions should be printed to the @p stream.
         */
        template <typename Stream>
        static void
        print(Stream                            &stream,
              typename Optimizer::OptimizerType *optimizer,
              const bool                         print_cse_reductions = true,
              const bool print_independent_symbols                    = false,
              const bool print_dependent_functions                    = false)
        {
          Assert(optimizer, ExcNotInitialized());

          optimizer->print(stream,
                           print_independent_symbols,
                           print_dependent_functions,
                           print_cse_reductions);
        }
      };

#  endif // DOXYGEN


      /* -------------------- Utility functions ---------------------- */


      /**
       * A convenience function that returns the numeric equivalent of
       * an input @p symbol_tensor, computed through the @p optimizer.
       *
       * @tparam NumberType The number type that is returned as a result
       *         of operations performed by the optimizer.
       * @tparam rank The rank of the output tensor.
       * @tparam dim The dimension of the output tensor.
       * @tparam TensorType The type of tensor to be evaluated and returned
       *         (i.e. Tensor or SymmetricTensor).
       * @param[in] symbol_tensor The symbolic tensor that is to be evaluated.
       * @param[in] cached_evaluation A vector that stores the numerical
       *            values of all dependent variables referenced by the
       *            @p optimizer. This vector is, most typically, first attained
       *            by a call to the BatchOptimizer::evaluate() variant that
       *            takes no arguments.
       * @param[in] optimizer The optimizer that can evaluate the input
       *       @p symbol_tensor.
       * @return TensorType<rank, dim, NumberType> The numeric result that the
       *         input @p symbol_tensor evaluates to.
       */
      template <typename NumberType,
                int rank,
                int dim,
                template <int, int, typename>
                class TensorType>
      TensorType<rank, dim, NumberType>
      tensor_evaluate_optimized(
        const TensorType<rank, dim, Expression> &symbol_tensor,
        const std::vector<NumberType>           &cached_evaluation,
        const BatchOptimizer<NumberType>        &optimizer)
      {
        TensorType<rank, dim, NumberType> out;
        for (unsigned int i = 0; i < out.n_independent_components; ++i)
          {
            const TableIndices<rank> indices(
              out.unrolled_to_component_indices(i));
            out[indices] =
              optimizer.extract(symbol_tensor[indices], cached_evaluation);
          }
        return out;
      }


      /**
       * A convenience function that returns the numeric equivalent of
       * an input @p symbol_tensor, computed through the @p optimizer.
       * This is a specialization for rank-4 symmetric tensors.
       *
       * @tparam NumberType The number type that is returned as a result
       *         of operations performed by the optimizer.
       * @tparam rank The rank of the output tensor.
       * @tparam dim The dimension of the output tensor.
       * @tparam TensorType The type of tensor to be evaluated and returned
       *         (i.e. Tensor or SymmetricTensor).
       * @param[in] symbol_tensor The symbolic tensor that is to be evaluated.
       * @param[in] cached_evaluation A vector that stores the numerical
       *            values of all dependent variables referenced by the
       *            @p optimizer. This vector is, most typically, first attained
       *            by a call to the BatchOptimizer::evaluate() variant that
       *            takes no arguments.
       * @param[in] optimizer The optimizer that can evaluate the input
       *       @p symbol_tensor.
       * @return TensorType<rank, dim, NumberType> The numeric result that the
       *         input @p symbol_tensor evaluates to.
       */
      template <typename NumberType, int dim>
      SymmetricTensor<4, dim, NumberType>
      tensor_evaluate_optimized(
        const SymmetricTensor<4, dim, Expression> &symbol_tensor,
        const std::vector<NumberType>             &cached_evaluation,
        const BatchOptimizer<NumberType>          &optimizer)
      {
        SymmetricTensor<4, dim, NumberType> out;
        for (unsigned int i = 0;
             i < SymmetricTensor<2, dim>::n_independent_components;
             ++i)
          for (unsigned int j = 0;
               j < SymmetricTensor<2, dim>::n_independent_components;
               ++j)
            {
              const TableIndices<4> indices =
                make_rank_4_tensor_indices<dim>(i, j);
              out[indices] =
                optimizer.extract(symbol_tensor[indices], cached_evaluation);
            }
        return out;
      }


      /**
       * A helper function to register a single @p function with the
       * @p optimizer.
       *
       * @tparam NumberType The number type that is returned as a result
       *         of operations performed by the optimizer.
       * @tparam T A compatible type that may be used to represent a single
       *         dependent variable. This includes scalar Expressions, Tensors
       *         of Expressions and SymmetricTensors of Expressions.
       * @param optimizer The instance of the BatchOptimizer to register the
       *        @p function with.
       * @param function A symbolic expression (scalar or tensor) that
       *        represents a dependent variable.
       *
       * @note This is the end-point for all recursive template functions
       *       with the same name.
       */
      template <typename NumberType, typename T>
      void
      register_functions(BatchOptimizer<NumberType> &optimizer,
                         const T                    &function)
      {
        optimizer.register_function(function);
      }


      /**
       * A helper function to register a vector of @p functions with the
       * @p optimizer.
       *
       * @tparam NumberType The number type that is returned as a result
       *         of operations performed by the optimizer.
       * @tparam T A compatible type that may be used to represent a single
       *         dependent variable. This includes scalar Expressions, Tensors
       *         of Expressions and SymmetricTensors of Expressions.
       * @param optimizer The instance of the BatchOptimizer to register the
       *        @p function with.
       * @param functions A vector of symbolic expressions (scalar or tensor)
       *        that each represent a dependent variable.
       *
       * @note This is the end-point for all recursive template functions
       *       with the same name.
       */
      template <typename NumberType, typename T>
      void
      register_functions(BatchOptimizer<NumberType> &optimizer,
                         const std::vector<T>       &functions)
      {
        for (const auto &function : functions)
          register_functions(optimizer, function);
      }


      /**
       * A helper function to register the symbolic dependent variables
       * collectively given by @p function and @p other_functions with the
       * @p optimizer.
       *
       * @tparam NumberType The number type that is returned as a result
       *         of operations performed by the optimizer.
       * @tparam T A compatible type that may be used to represent a single
       *         dependent variable. This includes scalar Expressions, Tensors
       *         of Expressions and SymmetricTensors of Expressions.
       * @tparam Args The parameter pack that collects all other types of
       *         dependent variables to be registered.
       * @param optimizer The instance of the BatchOptimizer to register the
       *        @p function with.
       * @param function A valid symbolic expression (or collection of symbolic
       *        expression) that represents one (or more) dependent variable.
       * @param other_functions One or more other valid symbolic expression(s)
       *        that represent dependent variable(s).
       */
      template <typename NumberType, typename T, typename... Args>
      void
      register_functions(BatchOptimizer<NumberType> &optimizer,
                         const T                    &function,
                         const Args &...other_functions)
      {
        register_functions(optimizer, function);
        register_functions(optimizer, other_functions...);
      }


      /**
       * A utility function that unrolls the input @p symbol_tensor into
       * a vector of Expressions.
       *
       * @tparam rank The rank of the input tensor.
       * @tparam dim The dimension of the input tensor.
       * @tparam TensorType The type of tensor to be evaluated and returned
       *         (i.e. Tensor or SymmetricTensor).
       * @param symbol_tensor
       * @return A vector of Expressions, with a consistent ordering.
       */
      template <int rank,
                int dim,
                template <int, int, typename>
                class TensorType>
      types::symbol_vector
      unroll_to_expression_vector(
        const TensorType<rank, dim, Expression> &symbol_tensor)
      {
        SD::types::symbol_vector out;
        out.reserve(symbol_tensor.n_independent_components);
        for (unsigned int i = 0; i < symbol_tensor.n_independent_components;
             ++i)
          {
            const TableIndices<rank> indices(
              symbol_tensor.unrolled_to_component_indices(i));
            out.push_back(symbol_tensor[indices].get_RCP());
          }
        return out;
      }


      /**
       * A utility function that unrolls the input @p symbol_tensor into
       * a vector of Expressions.
       * This is a specialization for rank-4 symmetric tensors.
       *
       * @tparam dim The dimension of the input tensor.
       * @param symbol_tensor
       * @return A vector of Expressions, with a consistent ordering.
       */
      template <int dim>
      types::symbol_vector
      unroll_to_expression_vector(
        const SymmetricTensor<4, dim, Expression> &symbol_tensor)
      {
        SD::types::symbol_vector out;
        out.reserve(symbol_tensor.n_independent_components);
        for (unsigned int i = 0;
             i < SymmetricTensor<2, dim>::n_independent_components;
             ++i)
          for (unsigned int j = 0;
               j < SymmetricTensor<2, dim>::n_independent_components;
               ++j)
            {
              const TableIndices<4> indices =
                make_rank_4_tensor_indices<dim>(i, j);
              out.push_back(symbol_tensor[indices].get_RCP());
            }
        return out;
      }

    } // namespace internal



    /**
     * A class that facilitates the optimization of symbol expressions.
     *
     * This expression will be optimized by this class; that is to say that
     * the code path taken to substitute the set of (independent) symbols
     * into a collection of (dependent) symbolic functions will be optimized
     * using a chosen approach.
     *
     * This snippet of pseudo-code describes the general usage of this class:
     * @code
     *
     * // Define some independent variables
     * const Expression x("x");
     * const Expression y("y");
     * ...
     *
     * // Compute some symbolic expressions that are dependent on the
     * // independent variables. These could be, for example, scalar
     * // expressions or tensors of expressions.
     * const auto f = calculate_f(x, y, ...);
     * const auto g = calculate_g(x, y, ...);
     * ...
     *
     * // Now create a optimizer to evaluate the dependent functions.
     * // The numerical result will be of type double, and a "lambda" optimizer,
     * // which employs common subexpression elimination, will be used.
     * using ReturnType = double;
     * BatchOptimizer<ReturnType> optimizer (OptimizerType::lambda,
     *                                       OptimizationFlags::optimize_cse);
     *
     * // Register symbols that represent independent variables...
     * optimizer.register_symbols(x, y, ...);
     * // ... and symbolic expressions that represent dependent functions.
     * optimizer.register_functions(f, g, ...);
     *
     * // Now we determine an equivalent code path that will evaluate
     * // all of the dependent functions at once, but with less computational
     * // cost than when evaluating the symbolic expression directly.
     * optimizer.optimize(); // Note: This is an expensive call.
     *
     * // Next we pass the optimizer the numeric values that we wish the
     * // independent variables to represent.
     * const auto substitution_map
     *   = make_substitution_map({x, ...}, {y, ...}, ...);
     * // When making this next call, the call path used to (numerically)
     * // evaluate the dependent functions is quicker than dictionary
     * // substitution.
     * optimizer.substitute(substitution_map);
     *
     * // Finally, we can get the numeric equivalent of the dependent functions
     * // from the optimizer.
     * const auto result_f = optimizer.evaluate(f);
     * const auto result_g = optimizer.evaluate(g);
     * @endcode
     *
     * Since the call to optimize() may be quite costly, there are a few "best
     * practices" that can be adopted in order to mitigate this cost as much
     * as possible:
     * 1. Reuse a single instance of the class as much as possible.
     *    The most obvious way that this can be achieved would be to place an
     *    instance of this class in a centralized location where it can
     *    potentially be used by multiple calling functions and objects, if
     *    contextually possible.
     * 2. Another form of reuse would entail generalizing the dependent
     *    functions/expressions to be evaluated by the optimizer as much as
     *    possible. For example, material coefficients need not necessarily be
     *    hard-coded, and one generalized statement of a constitutive law could
     *    then be broadly used in other material subdomains governed by the
     *    same class of constitutive law, but with different constitutive
     *    parameters. The same principle applies if using symbolic expressions
     *    to describe boundary conditions, systems of linear equations, etc.
     * 3. When possible, consider using serialization to save and load the state
     *    of an optimizer that has already "optimized", i.e. it has been placed
     *    in a state where it is ready to evaluate expressions.
     *    With the exception of "lambda" optimization, all other forms of
     *    optimization permit checkpointing, meaning that the optimization could
     *    be done up front before executing the main body of code.
     *    It could also be used to duplicate an optimizer in an efficient
     *    manner, should multiple instances of the same optimizer be required.
     *
     * @tparam ReturnType The number type that is to be returned after
     * value substitution and evaluation. Floating point and complex numbers
     * are currently supported.
     *
     * @warning This class is not thread-safe.
     *
     * @warning The LLVM optimizer does not yet support complex numbers. If this
     * incompatible combination of @p ReturnType and optimization method are
     * selected, then an error will be thrown at run time.
     */
    template <typename ReturnType>
    class BatchOptimizer
    {
    public:
      /**
       * Default constructor.
       *
       * By default, dictionary substitution will be selected when this
       * constructor is called. In order to select a specific optimization
       * approach, a call to set_optimization_method() is necessary.
       */
      BatchOptimizer();

      /**
       * Constructor.
       *
       * @param[in] optimization_method The optimization method that is to be
       * employed.
       * @param[in] optimization_flags  The optimization flags that indicate
       * which expression manipulation mechanisms are to be employed.
       *
       * @note As the optimization method is fixed, a further call to
       * set_optimization_method() is not necessary and will result in an
       * error being thrown.
       *
       * @note In the case that the @p optimization_method is not implemented for the
       * required @p ReturnType, or the desired feature is not active,
       * then an error will be thrown. Currently the LLVM optimization method
       * is not compatible with complex numbers.
       */
      BatchOptimizer(const enum OptimizerType     &optimization_method,
                     const enum OptimizationFlags &optimization_flags =
                       OptimizationFlags::optimize_all);

      /**
       * Copy constructor.
       *
       * @note The optimized data and results from previous substitutions
       * executed by the @p other optimizer instance are not copied over.
       * It is therefore necessary to re-optimize the data stored in
       * this class, and it is possible to do so with a different optimization
       * scheme.
       */
      BatchOptimizer(const BatchOptimizer &other);

      /**
       * Move constructor.
       */
      BatchOptimizer(BatchOptimizer &&) noexcept = default;

      /**
       * Destructor.
       */
      ~BatchOptimizer() = default;

      /**
       * Duplicate the data stored in an @p other BatchOptimizer instance.
       *
       * @note The optimized data and results from previous substitutions
       * executed by the @p other optimizer instance are not copied over.
       * It is therefore necessary to call optimize() before it is possible to
       * substitute() values and evaluate() data. One may, however, still
       * extract() values using @p this optimizer instance if those results are
       * stored elsewhere.
       */
      void
      copy_from(const BatchOptimizer &other);

      /**
       * Print some information on state of the internal data
       * structures stored in the class.
       *
       * @tparam Stream The type for the output stream.
       * @param stream The output stream to print to.
       * @param print_cse A flag to indicate whether or not all common
       *                  subexpressions should be printed to the @p stream.
       */
      template <typename Stream>
      void
      print(Stream &stream, const bool print_cse = false) const;

      /**
       * Write the data of this object from a stream for the purpose
       * of serialization using the [BOOST serialization
       * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
       *
       * This effectively saves the value stored into the @p archive with the
       * given @p version number into this object.
       */
      template <class Archive>
      void
      save(Archive &archive, const unsigned int version) const;

      /**
       * Read the data of this object from a stream for the purpose
       * of serialization using the [BOOST serialization
       * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
       *
       * This effectively loads the value stored out of the
       * @p archive with the given @p version number into this object.
       * In doing so, the previous contents of this object are thrown away.
       *
       * @note When deserializing a symbolic expression, it is imperative that
       * you first create or deserialize all of the symbolic variables used in
       * the serialized expression.
       */
      template <class Archive>
      void
      load(Archive &archive, const unsigned int version);

#  ifdef DOXYGEN
      /**
       * Write and read the data of this object from a stream for the purpose
       * of serialization using the [BOOST serialization
       * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
       *
       * This effectively saves or loads the value stored into/out of the
       * @p archive with the given @p version number into this object.
       * If deserializing data, then the previous contents of this object
       * are thrown away.
       *
       * @note When deserializing a batch optimizer, it is imperative that
       * you first create or deserialize all of the symbolic variables and
       * symbolic functions used in the optimizer.
       *
       * @note Complete serialization is not possible when the "lambda"
       * optimization method is invoked. Although the registered symbols
       * and dependent function expressions are stored, the optimization
       * is itself not stored. It might, therefore, take some time for the
       * deserialization when "lambda" optimization is used as the optimization
       * step will be (automatically) performed once more.
       */
      template <class Archive>
      void
      serialize(Archive &archive, const unsigned int version);
#  else
      // This macro defines the serialize() method that is compatible with
      // the templated save() and load() method that have been implemented.
      BOOST_SERIALIZATION_SPLIT_MEMBER()
#  endif

      /**
       * @name Independent variables
       */
      /** @{ */

      /**
       * Register a collection of symbols that represents an independent
       * variable. These symbols are stored as the <tt>key</tt> to
       * the @p substitution_map.
       */
      void
      register_symbols(const types::substitution_map &substitution_map);

      /**
       * Register a collection of symbols that represents an independent
       * variable. These symbols are stored as the <tt>key</tt> to
       * the @p substitution_map.
       */
      void
      register_symbols(const SymEngine::map_basic_basic &substitution_map);

      /**
       * Register a collection of symbols that represents independent variables.
       *
       * @warning When using this function is no mechanism to check that the ordering
       * of the later used @p substitution_values vector or map matches the internal
       * ordering of the registered symbols. This function is therefore
       * typically used in conjunction with the substitute() function that takes
       * in a vector of values. With this pair of functions to the class
       * interface, the management of symbol ordering is maintained by the user.
       */
      void
      register_symbols(const types::symbol_vector &symbols);

      /**
       * Register a collection of symbols that represents independent variables.
       *
       * @warning When using this function is no mechanism to check that the ordering
       * of the later used @p substitution_values vector or map matches the internal
       * ordering of the registered symbols. This function is therefore
       * typically used in conjunction with the substitute() function that takes
       * in a vector of values. With this pair of functions to the class
       * interface, the management of symbol ordering is maintained by the user.
       */
      void
      register_symbols(const SymEngine::vec_basic &symbols);

      /**
       * Return a vector of symbols that have been registered as independent
       * variables.
       */
      types::symbol_vector
      get_independent_symbols() const;

      /**
       * The number of independent variables that this optimizer will recognize.
       * This is equal to the number of unique symbols passed to this class
       * instance through the register_symbols() function.
       */
      std::size_t
      n_independent_variables() const;

      /** @} */

      /**
       * @name Dependent variables
       */
      /** @{ */

      /**
       * Register a scalar symbolic expression that represents a dependent
       * variable.
       */
      void
      register_function(const Expression &function);

      /**
       * Register a tensor of symbolic expressions that represents a dependent
       * variable.
       */
      template <int rank, int dim>
      void
      register_function(const Tensor<rank, dim, Expression> &function_tensor);

      /**
       * Register a symmetric tensor of symbolic expressions that represents a
       * dependent variable.
       */
      template <int rank, int dim>
      void
      register_function(
        const SymmetricTensor<rank, dim, Expression> &function_tensor);

      /**
       * Register a collection of symbolic expressions that represent dependent
       * variables.
       */
      void
      register_functions(const types::symbol_vector &functions);

      /**
       * Register a collection of symbolic expressions that represent multiple
       * dependent variables.
       */
      void
      register_functions(const SymEngine::vec_basic &functions);

      /**
       * Register a collection of symbolic expressions that represent multiple
       * dependent variables.
       *
       * @tparam T A compatible type that may be used to represent a single dependent
       *         variable. This includes scalar Expressions, Tensors of
       *         Expressions and SymmetricTensors of Expressions.
       * @param functions A vector of symbolic dependent variables.
       */
      template <typename T>
      void
      register_functions(const std::vector<T> &functions);

      /**
       * Register a collection of symbolic expressions that represent dependent
       * variables.
       *
       * @tparam T A compatible type that may be used to represent a single dependent
       *         variable. This includes scalar Expressions, Tensors of
       * Expressions, SymmetricTensors of Expressions, and `std::vector`s
       * of Expressions.
       * @tparam Args A variadic template that represents a collection of any compatible
       *         symbolic dependent variable types.
       * @param functions One or more symbolic dependent variables.
       * @param other_functions An arbitrary collection of symbolic dependent variables.
       */
      template <typename T, typename... Args>
      void
      register_functions(const T &functions, const Args &...other_functions);

      /**
       * Return a vector of expressions that have been registered as dependent
       * variables.
       */
      const types::symbol_vector &
      get_dependent_functions() const;

      /**
       * The number of dependent symbolic expressions that this optimizer
       * will optimize. This is equal to the number of unique symbolic
       * functions / expressions passed to this class instance through the
       * register_functions() method.
       */
      std::size_t
      n_dependent_variables() const;

      /** @} */

      /**
       * @name Optimization
       */
      /** @{ */

      /**
       * Select the @p optimization_method for the batch optimizer to
       * employ, in conjunction with the @p optimization_flags.
       *
       * It is required that the this class instance is not yet optimized, i.e.
       * that the optimize() method has not yet been called.
       *
       * @note In the case that the @p method is not implemented for the
       * required @p ReturnType, or the desired feature is not active,
       * then a safe default will be selected.
       */
      void
      set_optimization_method(const enum OptimizerType     &optimization_method,
                              const enum OptimizationFlags &optimization_flags =
                                OptimizationFlags::optimize_all);

      /**
       * Return the optimization method that has been
       * selected for use.
       */
      enum OptimizerType
      optimization_method() const;

      /**
       * Return the optimization flags that have been
       * selected for use.
       */
      enum OptimizationFlags
      optimization_flags() const;

      /**
       * State whether the internal selection of optimization
       * methods and flags will render an optimizer that uses
       * common subexpression elimination (CSE).
       */
      bool
      use_symbolic_CSE() const;

      /**
       * Perform the optimization of all registered dependent functions using
       * the registered symbols.
       *
       * @note This function, which should only be called once per instance of this
       * class, finalizes the set of accepted independent symbols and dependent
       * functions that are recognized and used by the optimizer.
       *
       * @note This may be a time-consuming process, but if the class instance is
       * retained throughout the course of a simulation (and both the
       * independent and dependent variables that are associated with the class
       * instance remain unchanged) then it need only be performed once.
       * Serialization also offers the opportunity to reuse the already computed
       * optimized evaluation call path.
       */
      void
      optimize();

      /**
       * Returns a flag which indicates whether the optimize()
       * function has been called and the class is finalized.
       */
      bool
      optimized() const;

      /** @} */

      /**
       * @name Symbol substitution
       */
      /** @{ */

      /**
       * Perform batch substitution of all of the registered symbols
       * into the registered functions. The result is cached and can
       * be extracted by calls to evaluate().
       *
       * @note Calling substitute() again with a new
       * @p substitution_map overwrites any previously computed
       * results.
       */
      void
      substitute(const types::substitution_map &substitution_map) const;

      /**
       * Perform batch substitution of all of the registered symbols
       * into the registered functions. The result is cached and can
       * be extracted by calls to evaluate().
       *
       * @note Calling substitute() again with a new
       * @p substitution_map overwrites any previously computed
       * results.
       */
      void
      substitute(const SymEngine::map_basic_basic &substitution_map) const;

      /**
       * Perform batch substitution of all of the registered symbols
       * into the registered functions. The result is cached and can
       * be extracted by calls to evaluate().
       * It is expected that there is a 1-1 correspondence between each
       * of the @p symbols and @p values.
       *
       * @note Calling substitute() again with a new set of
       * @p values overwrites any previously computed results.
       */
      void
      substitute(const types::symbol_vector    &symbols,
                 const std::vector<ReturnType> &values) const;

      /**
       * Perform batch substitution of all of the registered symbols
       * into the registered functions. The result is cached and can
       * be extracted by calls to evaluate().
       * It is expected that there is a 1-1 correspondence between each
       * of the @p symbols and @p values.
       *
       * @note Calling substitute() again with a new set of
       * @p values overwrites any previously computed results.
       */
      void
      substitute(const SymEngine::vec_basic    &symbols,
                 const std::vector<ReturnType> &values) const;

      /**
       * Returns a flag to indicate whether the substitute()
       * function has been called and if there are meaningful
       * values that will be returned upon evaluation.
       */
      bool
      values_substituted() const;

      /** @} */

      /**
       * @name Evaluation / data extraction
       */
      /** @{ */

      /**
       * Returns the result of a value substitution into the optimized
       * counterpart of all dependent functions. This function fetches all of
       * those cached values. If these results are stored elsewhere, and you
       * want to use this optimizer to extract components of the cached results
       * then you can make use of the extract() functions listed below.
       *
       * These values were computed by substituting a @p substitution_values map
       * during substitute() call.
       *
       * @note In contrast to the other variants of this function, the order of
       * the returned entries is an internal implementation detail, and cannot
       * be guaranteed under all conditions. The entries in the returned vector
       * are, in general, identical to the order in which the dependent
       * expressions are originally registered. However, when registering
       * tensors and symmetric tensors of expressions, these are "unrolled"
       * and their components are individually registered. If it is necessary to
       * control the order in which results appear in the returned vector, then
       * the register_function() method that takes a Tensor or a SymmetricTensor
       * as an argument should be avoided. Instead, the individual entries of
       * these types of data should be registered one by one.
       */
      const std::vector<ReturnType> &
      evaluate() const;

      /**
       * Returns the result of a value substitution into the optimized
       * counterpart of @p func. This function fetches that one cached value.
       *
       * This value was computed by substituting a @p substitution_values map
       * during substitute() call.
       */
      ReturnType
      evaluate(const Expression &func) const;

      /**
       * Returns the result of a value substitution into the optimized
       * counterpart of @p funcs. This function fetches that subset of cached
       * values.
       *
       * This value was computed by substituting a @p substitution_values map
       * during substitute() call.
       */
      std::vector<ReturnType>
      evaluate(const std::vector<Expression> &funcs) const;

      /**
       * Returns the result of a tensor value substitution into the optimized
       * counterpart of @p funcs. This function fetches those cached tensor
       * components.
       *
       * This value was computed by substituting a @p substitution_values map
       * during substitute() call.
       */
      template <int rank, int dim>
      Tensor<rank, dim, ReturnType>
      evaluate(const Tensor<rank, dim, Expression> &funcs) const;


      /**
       * Returns the result of a tensor value substitution into the optimized
       * counterpart of @p funcs. This function fetches those cached symmetric
       * tensor components.
       *
       * This value was computed by substituting a @p substitution_values map
       * during substitute() call.
       */
      template <int rank, int dim>
      SymmetricTensor<rank, dim, ReturnType>
      evaluate(const SymmetricTensor<rank, dim, Expression> &funcs) const;


      /**
       * Returns the result of a value substitution into the optimized
       * counterpart of @p func. This function fetches that one value
       * stored in the @p cached_evaluation vector. This vector is, most
       * typically, first attained by a call to the evaluate() variant that
       * takes no arguments.
       */
      ReturnType
      extract(const Expression              &func,
              const std::vector<ReturnType> &cached_evaluation) const;


      /**
       * Returns the result of a value substitution into the optimized
       * counterpart of @p funcs. This function fetches that subset of
       * values stored in the @p cached_evaluation vector. This vector is, most
       * typically, first attained by a call to the evaluate() variant that
       * takes no arguments.
       */
      std::vector<ReturnType>
      extract(const std::vector<Expression> &funcs,
              const std::vector<ReturnType> &cached_evaluation) const;


      /**
       * Returns the result of a tensor value substitution into the optimized
       * counterpart of @p funcs. This function fetches those tensor components
       * stored in the @p cached_evaluation vector. This vector is, most
       * typically, first attained by a call to the evaluate() variant that
       * takes no arguments.
       */
      template <int rank, int dim>
      Tensor<rank, dim, ReturnType>
      extract(const Tensor<rank, dim, Expression> &funcs,
              const std::vector<ReturnType>       &cached_evaluation) const;


      /**
       * Returns the result of a tensor value substitution into the optimized
       * counterpart of @p funcs. This function fetches those symmetric
       * tensor components stored in the @p cached_evaluation vector. This vector
       * is, most typically, first attained by a call to the evaluate() variant
       * that takes no arguments.
       */
      template <int rank, int dim>
      SymmetricTensor<rank, dim, ReturnType>
      extract(const SymmetricTensor<rank, dim, Expression> &funcs,
              const std::vector<ReturnType> &cached_evaluation) const;

      /** @} */

    private:
      /**
       * The optimization methods that is to be employed.
       */
      enum OptimizerType method;

      /**
       * The optimization flags that indicate which expression manipulation
       * mechanisms are to be employed.
       */
      enum OptimizationFlags flags;

      /**
       * A map that represents the symbols that form the set of independent
       * variables upon which optimized symbolic expressions are to be based.
       *
       * @note As the ordering of the input symbols is fixed at the time at
       * which optimization is performed, we store all of the entries in
       * a map to ensure that both we and the user never mistakenly
       * swap the order of two or more symbols during evaluation.
       */
      types::substitution_map independent_variables_symbols;

      /**
       * A set of symbolic expressions to be optimized. It is required that
       * the symbols on which these dependent functions are based are
       * registered in the @p independent_variables_symbols map.
       */
      types::symbol_vector dependent_variables_functions;

      /**
       * A check to see if a function is exactly equal to one of the logical
       * results of a differentiation operation.
       */
      bool
      is_valid_nonunique_dependent_variable(
        const SD::Expression &function) const;

      /**
       * A check to see if a function is exactly equal to one of the logical
       * results of a differentiation operation.
       */
      bool
      is_valid_nonunique_dependent_variable(
        const SymEngine::RCP<const SymEngine::Basic> &function) const;

      /**
       * The output of substituting symbolic values with floating point
       * values through the use of the @p optimizer.
       *
       * @p It is necessary to use this intermediate storage mechanism
       * to store the result of a substitution sweep as some optimizers
       * work on all symbolic expressions in a single batch. In this way
       * they can employ methods such as common subexpression elimination
       * to minimize the number of terms evaluated across all symbolic
       * functions.
       *
       * This variable is marked as mutable. This facilitates the substitution
       * functionality being used in logically constant `get_*` functions.
       */
      mutable std::vector<ReturnType> dependent_variables_output;

      /**
       * A map type used to indicate which dependent variable is associated
       * with which entry in the output vector.
       *
       * @note We use a custom comparator here because otherwise we can't use
       * std::map::find; it is sensitive to the some data from the
       * underlying SymEngine::Basic other than the value that it represents.
       */
      using map_dependent_expression_to_vector_entry_t =
        std::map<SD::Expression,
                 std::size_t,
                 SD::types::internal::ExpressionKeyLess>;

      /**
       * A map indicating which dependent variable is associated with which
       * entry in the output vector.
       */
      mutable map_dependent_expression_to_vector_entry_t map_dep_expr_vec_entry;

      /**
       * A pointer to an instance of an optimizer that will be used to
       * reformulate the substitution of symbolic expressions in a
       * manner that is more efficient than plain dictionary-based
       * approach.
       */
      mutable std::unique_ptr<SymEngine::Visitor> optimizer;

      /**
       * A flag to record whether or not substitution has taken place and
       * values can now be extracted.
       *
       * @note This variable is marked as mutable. This facilitates the
       * substitution functionality being used in logically constant
       * `get_*` functions.
       */
      mutable bool ready_for_value_extraction;

      /**
       * A flag to record whether or not this class instance has been serialized
       * in the past.
       */
      mutable bool has_been_serialized;

      /**
       * Register a single symbol that represents a dependent variable.
       */
      void
      register_scalar_function(const SD::Expression &function);

      /**
       * Register a collection of symbols that represent dependent
       * variables.
       */
      void
      register_vector_functions(const types::symbol_vector &functions);

      /**
       * Create an instance of the selected optimizer.
       */
      void
      create_optimizer(std::unique_ptr<SymEngine::Visitor> &optimizer);

      /**
       * Perform batch substitution of all of the registered symbols
       * into the registered functions. The result is cached and can
       * be extracted by calls to evaluate().
       *
       * @note Calling substitute() again with a new set of
       * @p substitution_values overwrites any previously computed
       * results.
       *
       * @warning When using this function there is no mechanism to check that
       * the ordering of the @p substitution_values vector matches the internal
       * ordering of the registered symbols. This function is therefore
       * typically used in conjunction with the register_symbols() function that
       * takes in a vector of symbols. With this pair of functions to the class
       * interface, the management of symbol ordering is maintained by the user.
       */
      void
      substitute(const std::vector<ReturnType> &substitution_values) const;
    };



    /* -------------------- inline and template functions ------------------ */


#  ifndef DOXYGEN


    template <typename ReturnType>
    template <typename Stream>
    void
    BatchOptimizer<ReturnType>::print(Stream &stream,
                                      const bool /*print_cse*/) const
    {
      // Settings
      stream << "Method? " << optimization_method() << '\n';
      stream << "Flags: " << optimization_flags() << '\n';
      stream << "Optimized? " << (optimized() ? "Yes" : "No") << '\n';
      stream << "Values substituted? " << values_substituted() << "\n\n";

      // Independent variables
      stream << "Symbols (" << n_independent_variables()
             << " independent variables):" << '\n';
      int cntr = 0;
      for (SD::types::substitution_map::const_iterator it =
             independent_variables_symbols.begin();
           it != independent_variables_symbols.end();
           ++it, ++cntr)
        {
          stream << cntr << ": " << it->first << '\n';
        }
      stream << '\n' << std::flush;

      // Dependent functions
      stream << "Functions (" << n_dependent_variables()
             << " dependent variables):" << '\n';
      cntr = 0;
      for (typename SD::types::symbol_vector::const_iterator it =
             dependent_variables_functions.begin();
           it != dependent_variables_functions.end();
           ++it, ++cntr)
        {
          stream << cntr << ": " << (*it) << '\n';
        }
      stream << '\n' << std::flush;

      // Common subexpression
      if (optimized() == true && use_symbolic_CSE() == true)
        {
          Assert(optimizer, ExcNotInitialized());
          const bool print_cse_reductions      = true;
          const bool print_independent_symbols = false;
          const bool print_dependent_functions = false;

          if (optimization_method() == OptimizerType::dictionary)
            {
              Assert(dynamic_cast<typename internal::DictionaryOptimizer<
                       ReturnType>::OptimizerType *>(optimizer.get()),
                     ExcMessage("Cannot cast optimizer to Dictionary type."));

              internal::OptimizerHelper<
                ReturnType,
                internal::DictionaryOptimizer<ReturnType>>::
                print(stream,
                      dynamic_cast<typename internal::DictionaryOptimizer<
                        ReturnType>::OptimizerType *>(optimizer.get()),
                      print_independent_symbols,
                      print_dependent_functions,
                      print_cse_reductions);

              stream << '\n' << std::flush;
            }
          else if (optimization_method() == OptimizerType::lambda)
            {
              Assert(dynamic_cast<typename internal::LambdaOptimizer<
                       ReturnType>::OptimizerType *>(optimizer.get()),
                     ExcMessage("Cannot cast optimizer to Lambda type."));

              internal::OptimizerHelper<ReturnType,
                                        internal::LambdaOptimizer<ReturnType>>::
                print(stream,
                      dynamic_cast<typename internal::LambdaOptimizer<
                        ReturnType>::OptimizerType *>(optimizer.get()),
                      print_independent_symbols,
                      print_dependent_functions,
                      print_cse_reductions);
            }
#    ifdef HAVE_SYMENGINE_LLVM
          else if (optimization_method() == OptimizerType::llvm)
            {
              Assert(dynamic_cast<typename internal::LLVMOptimizer<
                       ReturnType>::OptimizerType *>(optimizer.get()),
                     ExcMessage("Cannot cast optimizer to LLVM type."));

              internal::OptimizerHelper<ReturnType,
                                        internal::LLVMOptimizer<ReturnType>>::
                print(stream,
                      dynamic_cast<typename internal::LLVMOptimizer<
                        ReturnType>::OptimizerType *>(optimizer.get()),
                      print_independent_symbols,
                      print_dependent_functions,
                      print_cse_reductions);
            }
#    endif // HAVE_SYMENGINE_LLVM
          else
            {
              AssertThrow(false, ExcMessage("Unknown optimizer type."));
            }
        }

      if (values_substituted())
        {
          stream << "Evaluated functions:" << '\n';
          stream << std::flush;
          cntr = 0;
          for (typename std::vector<ReturnType>::const_iterator it =
                 dependent_variables_output.begin();
               it != dependent_variables_output.end();
               ++it, ++cntr)
            {
              stream << cntr << ": " << (*it) << '\n';
            }
          stream << '\n' << std::flush;
        }
    }



    template <typename ReturnType>
    template <class Archive>
    void
    BatchOptimizer<ReturnType>::save(Archive           &ar,
                                     const unsigned int version) const
    {
      // Serialize enum classes...
      {
        const auto m =
          static_cast<std::underlying_type_t<OptimizerType>>(method);
        ar &m;
      }
      {
        const auto f =
          static_cast<std::underlying_type_t<OptimizationFlags>>(flags);
        ar &f;
      }

      // Important: Independent variables must always be
      // serialized before the dependent variables.
      ar &independent_variables_symbols;
      ar &dependent_variables_functions;

      ar &dependent_variables_output;
      ar &map_dep_expr_vec_entry;
      ar &ready_for_value_extraction;

      // Mark that we've saved this class at some point.
      has_been_serialized = true;
      ar &has_been_serialized;

      // When we serialize the optimizer itself, we have to (unfortunately)
      // provide it with sufficient information to rebuild itself from scratch.
      // This is because only two of the three optimization classes support
      // real serialization (i.e. have save/load capability).
      const SD::types::symbol_vector symbol_vec =
        Utilities::extract_symbols(independent_variables_symbols);
      if (typename internal::DictionaryOptimizer<ReturnType>::OptimizerType
            *opt = dynamic_cast<typename internal::DictionaryOptimizer<
              ReturnType>::OptimizerType *>(optimizer.get()))
        {
          Assert(optimization_method() == OptimizerType::dictionary,
                 ExcInternalError());
          internal::OptimizerHelper<
            ReturnType,
            internal::DictionaryOptimizer<ReturnType>>::save(ar, version, opt);
        }
      else if (typename internal::LambdaOptimizer<ReturnType>::OptimizerType
                 *opt = dynamic_cast<typename internal::LambdaOptimizer<
                   ReturnType>::OptimizerType *>(optimizer.get()))
        {
          Assert(optimization_method() == OptimizerType::lambda,
                 ExcInternalError());
          internal::OptimizerHelper<
            ReturnType,
            internal::LambdaOptimizer<ReturnType>>::save(ar, version, opt);
        }
#    ifdef HAVE_SYMENGINE_LLVM
      else if (typename internal::LLVMOptimizer<ReturnType>::OptimizerType
                 *opt = dynamic_cast<typename internal::LLVMOptimizer<
                   ReturnType>::OptimizerType *>(optimizer.get()))
        {
          Assert(optimization_method() == OptimizerType::llvm,
                 ExcInternalError());
          internal::OptimizerHelper<
            ReturnType,
            internal::LLVMOptimizer<ReturnType>>::save(ar, version, opt);
        }
#    endif
      else
        {
          AssertThrow(false, ExcMessage("Unknown optimizer type."));
        }
    }



    template <typename ReturnType>
    template <class Archive>
    void
    BatchOptimizer<ReturnType>::load(Archive &ar, const unsigned int version)
    {
      Assert(independent_variables_symbols.empty(), ExcInternalError());
      Assert(dependent_variables_functions.empty(), ExcInternalError());
      Assert(dependent_variables_output.empty(), ExcInternalError());
      Assert(map_dep_expr_vec_entry.empty(), ExcInternalError());
      Assert(ready_for_value_extraction == false, ExcInternalError());

      // Deserialize enum classes...
      {
        std::underlying_type_t<OptimizerType> m;
        ar                                   &m;
        method = static_cast<OptimizerType>(m);
      }
      {
        std::underlying_type_t<OptimizationFlags> f;
        ar                                       &f;
        flags = static_cast<OptimizationFlags>(f);
      }

      // Important: Independent variables must always be
      // deserialized before the dependent variables.
      ar &independent_variables_symbols;
      ar &dependent_variables_functions;

      ar &dependent_variables_output;
      ar &map_dep_expr_vec_entry;
      ar &ready_for_value_extraction;

      ar &has_been_serialized;

      // If we're reading in data, then create the optimizer
      // and then deserialize it.
      Assert(!optimizer, ExcInternalError());

      // Create and configure the optimizer
      create_optimizer(optimizer);
      Assert(optimizer, ExcNotInitialized());

      // When we deserialize the optimizer itself, we have to (unfortunately)
      // provide it with sufficient information to rebuild itself from scratch.
      // This is because only two of the three optimization classes support
      // real serialization (i.e. have save/load capability).
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
            load(ar,
                 version,
                 opt,
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
            load(ar,
                 version,
                 opt,
                 Utilities::convert_expression_vector_to_basic_vector(
                   symbol_vec),
                 Utilities::convert_expression_vector_to_basic_vector(
                   dependent_variables_functions),
                 optimization_flags());
        }
#    ifdef HAVE_SYMENGINE_LLVM
      else if (typename internal::LLVMOptimizer<ReturnType>::OptimizerType
                 *opt = dynamic_cast<typename internal::LLVMOptimizer<
                   ReturnType>::OptimizerType *>(optimizer.get()))
        {
          Assert(optimization_method() == OptimizerType::llvm,
                 ExcInternalError());
          internal::OptimizerHelper<ReturnType,
                                    internal::LLVMOptimizer<ReturnType>>::
            load(ar,
                 version,
                 opt,
                 Utilities::convert_expression_vector_to_basic_vector(
                   symbol_vec),
                 Utilities::convert_expression_vector_to_basic_vector(
                   dependent_variables_functions),
                 optimization_flags());
        }
#    endif
      else
        {
          AssertThrow(false, ExcMessage("Unknown optimizer type."));
        }
    }



    template <typename ReturnType>
    template <int rank, int dim>
    void
    BatchOptimizer<ReturnType>::register_function(
      const Tensor<rank, dim, Expression> &function_tensor)
    {
      Assert(optimized() == false,
             ExcMessage(
               "Cannot register functions once the optimizer is finalised."));

      register_vector_functions(
        internal::unroll_to_expression_vector(function_tensor));
    }



    template <typename ReturnType>
    template <int rank, int dim>
    void
    BatchOptimizer<ReturnType>::register_function(
      const SymmetricTensor<rank, dim, Expression> &function_tensor)
    {
      Assert(optimized() == false,
             ExcMessage(
               "Cannot register functions once the optimizer is finalised."));

      register_vector_functions(
        internal::unroll_to_expression_vector(function_tensor));
    }



    template <typename ReturnType>
    template <typename T, typename... Args>
    void
    BatchOptimizer<ReturnType>::register_functions(
      const T &functions,
      const Args &...other_functions)
    {
      internal::register_functions(*this, functions);
      internal::register_functions(*this, other_functions...);
    }



    template <typename ReturnType>
    template <typename T>
    void
    BatchOptimizer<ReturnType>::register_functions(
      const std::vector<T> &functions)
    {
      internal::register_functions(*this, functions);
    }



    template <typename ReturnType>
    template <int rank, int dim>
    Tensor<rank, dim, ReturnType>
    BatchOptimizer<ReturnType>::extract(
      const Tensor<rank, dim, Expression> &funcs,
      const std::vector<ReturnType>       &cached_evaluation) const
    {
      return internal::tensor_evaluate_optimized(funcs,
                                                 cached_evaluation,
                                                 *this);
    }



    template <typename ReturnType>
    template <int rank, int dim>
    Tensor<rank, dim, ReturnType>
    BatchOptimizer<ReturnType>::evaluate(
      const Tensor<rank, dim, Expression> &funcs) const
    {
      Assert(
        values_substituted() == true,
        ExcMessage(
          "The optimizer is not configured to perform evaluation. "
          "This action can only performed after substitute() has been called."));

      return extract(funcs, dependent_variables_output);
    }



    template <typename ReturnType>
    template <int rank, int dim>
    SymmetricTensor<rank, dim, ReturnType>
    BatchOptimizer<ReturnType>::extract(
      const SymmetricTensor<rank, dim, Expression> &funcs,
      const std::vector<ReturnType>                &cached_evaluation) const
    {
      return internal::tensor_evaluate_optimized(funcs,
                                                 cached_evaluation,
                                                 *this);
    }



    template <typename ReturnType>
    template <int rank, int dim>
    SymmetricTensor<rank, dim, ReturnType>
    BatchOptimizer<ReturnType>::evaluate(
      const SymmetricTensor<rank, dim, Expression> &funcs) const
    {
      Assert(
        values_substituted() == true,
        ExcMessage(
          "The optimizer is not configured to perform evaluation. "
          "This action can only performed after substitute() has been called."));

      return extract(funcs, dependent_variables_output);
    }

#  endif // DOXYGEN

  } // namespace SD
} // namespace Differentiation


DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SYMENGINE

#endif
