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

#ifndef dealii_differentiation_sd_symengine_optimizer_h
#define dealii_differentiation_sd_symengine_optimizer_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SYMENGINE

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
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
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#  include <deal.II/base/logstream.h>
#  include <deal.II/base/utilities.h>

#  include <deal.II/differentiation/sd/symengine_number_types.h>
#  include <deal.II/differentiation/sd/symengine_number_visitor.h>
#  include <deal.II/differentiation/sd/symengine_scalar_operations.h>
#  include <deal.II/differentiation/sd/symengine_tensor_operations.h>
#  include <deal.II/differentiation/sd/symengine_utilities.h>

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
    // Forward declarations
    template <typename ReturnType>
    class BatchOptimizer;


    /*
     * An enumeration to distinguish between different optimization methods
     * that can be used by SymEngine to more rapidly evaluate complex
     * symbolic expressions.
     */
    enum class OptimizerType
    {
      /*
       * Use no specific optimization method.
       */
      off,
      /*
       * Use dictionary substitution. This is SymEngine's default method.
       */
      dictionary = off,
      /*
       * Convert the symbolic expression into a collection of
       * <tt>std::function<\tt>s.
       */
      lambda,
      /*
       * Use the LLVM JIT compiler to compile the expression into an
       * aggressively optimized, stand-alone function.
       */
      llvm
    };

    /**
     * Output operator which outputs update flags as a set of or'd text values.
     *
     * @ref UpdateFlags
     */
    template <class StreamType>
    inline StreamType &
    operator<<(StreamType &s, OptimizerType o)
    {
      if (o == OptimizerType::off || o == OptimizerType::dictionary)
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
       * No additional optimization
       */
      optimize_default = 0,
      /**
       * Apply common subexpresson elimination
       */
      optimize_cse = 0x0001,
      //    /**
      //     * Apply expression simplification
      //     */
      //    optimize_simplify = 0x0002,
      /**
       * Employ aggressive optimizations when compiling with the LLVM JIT
       * compiler
       */
      optimize_aggressive = 0x0002,
      /**
       * Apply all possible optimizations
       */
      optimize_all = optimize_cse | optimize_aggressive
    };


    /**
     * Global operator which returns an object in which all bits are set which
     * are either set in the first or the second argument. This operator exists
     * since if it did not then the result of the bit-or <tt>operator |</tt>
     * would be an integer which would in turn trigger a compiler warning when
     * we tried to assign it to an object of type OptimizationFlags.
     *
     * @ref UpdateFlags
     */
    inline OptimizationFlags
    operator|(OptimizationFlags f1, OptimizationFlags f2)
    {
      return static_cast<OptimizationFlags>(static_cast<unsigned int>(f1) |
                                            static_cast<unsigned int>(f2));
    }


    /**
     * Global operator which sets the bits from the second argument also in the
     * first one.
     *
     * @ref UpdateFlags
     */
    inline OptimizationFlags &
    operator|=(OptimizationFlags &f1, OptimizationFlags f2)
    {
      f1 = f1 | f2;
      return f1;
    }


    /**
     * Global operator which returns an object in which all bits are set which
     * are set in the first as well as the second argument. This operator exists
     * since if it did not then the result of the bit-and <tt>operator &</tt>
     * would be an integer which would in turn trigger a compiler warning when
     * we tried to assign it to an object of type OptimizationFlags.
     *
     * @ref UpdateFlags
     */
    inline OptimizationFlags operator&(OptimizationFlags f1,
                                       OptimizationFlags f2)
    {
      return static_cast<OptimizationFlags>(static_cast<unsigned int>(f1) &
                                            static_cast<unsigned int>(f2));
    }


    /**
     * Global operator which clears all the bits in the first argument if they
     * are not also set in the second argument.
     *
     * @ref UpdateFlags
     */
    inline OptimizationFlags &
    operator&=(OptimizationFlags &f1, OptimizationFlags f2)
    {
      f1 = f1 & f2;
      return f1;
    }

    namespace internal
    {
      inline bool
      use_symbolic_CSE(const enum OptimizationFlags &flags)
      {
        return static_cast<int>(flags & OptimizationFlags::optimize_cse);
      }

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
      inline int
      get_LLVM_optimization_level(const enum OptimizationFlags &flags)
      {
        const bool use_agg_opt =
          static_cast<int>(flags & OptimizationFlags::optimize_aggressive);
        const int opt_level = (use_agg_opt ? 3 : 2);
        return opt_level;
      }
    } // namespace internal


    /**
     * Output operator which outputs update flags as a set of or'd text values.
     *
     * @ref UpdateFlags
     */
    template <class StreamType>
    inline StreamType &
    operator<<(StreamType &s, OptimizationFlags o)
    {
      s << " OptimizationFlags|";
      if (static_cast<unsigned int>(o & OptimizationFlags::optimize_cse))
        s << "cse|";

      // LLVM optimization level
      s << "-O" +
             dealii::Utilities::to_string(
               internal::get_LLVM_optimization_level(o)) +
             "|";

      return s;
    }


    namespace internal
    {
      namespace
      {
        // Dummy "optimizer" for dictionary-based
        // substitution
        template <typename ReturnType, typename = void>
        struct DictionaryOptimizer;

        // Define, create and initialise the lambda function
        // optimiser for each number type
        template <typename ReturnType, typename = void>
        struct LambdaOptimizer;

#  ifdef HAVE_SYMENGINE_LLVM
        template <typename ReturnType, typename = void>
        struct LLVMOptimizer;
#  endif

        // Floating point numbers
        template <typename ReturnType_>
        struct DictionaryOptimizer<
          ReturnType_,
          typename std::enable_if<std::is_arithmetic<ReturnType_>::value>::type>
        {
          typedef double ReturnType;
          typedef DictionarySubstitutionVisitor<ReturnType, SD::Expression>
            OptimizerType;

          static void
          initialize(OptimizerType &               optimizer,
                     const SymEngine::vec_basic &  symbol_vec,
                     const SymEngine::vec_basic &  dependent_variables,
                     const enum OptimizationFlags &optimization_flags)
          {
            const bool use_symbolic_cse = use_symbolic_CSE(optimization_flags);
            optimizer.init(symbol_vec, dependent_variables, use_symbolic_cse);
          }

          template <typename Stream>
          static void
          print(Stream &             stream,
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
        struct DictionaryOptimizer<
          ReturnType_,
          typename std::enable_if<
            boost::is_complex<ReturnType_>::value &&
            std::is_arithmetic<typename ReturnType_::value_type>::value>::type>
        {
          typedef std::complex<double> ReturnType;
          typedef DictionarySubstitutionVisitor<ReturnType, SD::Expression>
            OptimizerType;

          static void
          initialize(OptimizerType &               optimizer,
                     const SymEngine::vec_basic &  symbol_vec,
                     const SymEngine::vec_basic &  dependent_variables,
                     const enum OptimizationFlags &optimization_flags)
          {
            const bool use_symbolic_cse = use_symbolic_CSE(optimization_flags);
            optimizer.init(symbol_vec, dependent_variables, use_symbolic_cse);
          }

          template <typename Stream>
          static void
          print(Stream &             stream,
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
        struct LambdaOptimizer<
          ReturnType_,
          typename std::enable_if<std::is_arithmetic<ReturnType_>::value>::type>
        {
          typedef double                             ReturnType;
          typedef SymEngine::LambdaRealDoubleVisitor OptimizerType;

          static void
          initialize(OptimizerType &               optimizer,
                     const SymEngine::vec_basic &  symbol_vec,
                     const SymEngine::vec_basic &  dependent_variables,
                     const enum OptimizationFlags &optimization_flags)
          {
            const bool use_symbolic_cse = use_symbolic_CSE(optimization_flags);
            optimizer.init(symbol_vec, dependent_variables, use_symbolic_cse);
          }

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

        template <typename ReturnType_>
        struct LambdaOptimizer<
          ReturnType_,
          typename std::enable_if<
            boost::is_complex<ReturnType_>::value &&
            std::is_arithmetic<typename ReturnType_::value_type>::value>::type>
        {
          typedef std::complex<double>                  ReturnType;
          typedef SymEngine::LambdaComplexDoubleVisitor OptimizerType;

          static void
          initialize(OptimizerType &               optimizer,
                     const SymEngine::vec_basic &  symbol_vec,
                     const SymEngine::vec_basic &  dependent_variables,
                     const enum OptimizationFlags &optimization_flags)
          {
            const bool use_symbolic_cse = use_symbolic_CSE(optimization_flags);
            optimizer.init(symbol_vec, dependent_variables, use_symbolic_cse);
          }

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

#  ifdef HAVE_SYMENGINE_LLVM
        template <typename ReturnType_>
        struct LLVMOptimizer<
          ReturnType_,
          typename std::enable_if<std::is_arithmetic<ReturnType_>::value>::type>
        {
          typedef double                       ReturnType;
          typedef SymEngine::LLVMDoubleVisitor OptimizerType;

          static void
          initialize(OptimizerType &               optimizer,
                     const SymEngine::vec_basic &  symbol_vec,
                     const SymEngine::vec_basic &  dependent_variables,
                     const enum OptimizationFlags &optimization_flags)
          {
            const int opt_level =
              get_LLVM_optimization_level(optimization_flags);
            const bool use_symbolic_cse = use_symbolic_CSE(optimization_flags);
            optimizer.init(symbol_vec,
                           dependent_variables,
                           use_symbolic_cse,
                           opt_level);
          }

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

        // TODO:
        // There is no LLVM optimiser built with complex number support.
        // So we fall back to the LambdaDouble case instead.
        template <typename ReturnType_>
        struct LLVMOptimizer<
          ReturnType_,
          typename std::enable_if<
            boost::is_complex<ReturnType_>::value &&
            std::is_arithmetic<typename ReturnType_::value_type>::value>::type>
        {
          typedef typename LambdaOptimizer<ReturnType_>::ReturnType ReturnType;
          typedef
            typename LambdaOptimizer<ReturnType_>::OptimizerType OptimizerType;

          static void
          initialize(OptimizerType &               optimizer,
                     const SymEngine::vec_basic &  symbol_vec,
                     const SymEngine::vec_basic &  dependent_variables,
                     const enum OptimizationFlags &optimization_flags)
          {
            deallog
              << "WARNING: No complex number support exists for the LLVM optimizer. "
              << "A fallback option of the LambdaOptimizer has been selected instead."
              << std::endl;
            LambdaOptimizer<ReturnType_>::initialize(optimizer,
                                                     symbol_vec,
                                                     dependent_variables,
                                                     optimization_flags);
          }

          template <typename StreamType>
          static void
          print(StreamType &         stream,
                const OptimizerType &optimizer,
                const bool           print_independent_symbols = false,
                const bool           print_dependent_functions = false,
                const bool           print_cse_reductions      = true)
          {
            LambdaOptimizer<ReturnType_>::print(stream,
                                                optimizer,
                                                print_independent_symbols,
                                                print_dependent_functions,
                                                print_cse_reductions);
          }
        };
#  endif

        template <typename NumberType,
                  int rank,
                  int dim,
                  template <int, int, typename> class TensorType>
        TensorType<rank, dim, NumberType>
        tensor_evaluate_optimized(
          const TensorType<rank, dim, Expression> &symbol_tensor,
          const BatchOptimizer<NumberType> &       optimizer)
        {
          TensorType<rank, dim, NumberType> out;
          for (unsigned int i = 0; i < out.n_independent_components; ++i)
            {
              const TableIndices<rank> indices(
                out.unrolled_to_component_indices(i));
              out[indices] = optimizer.evaluate(symbol_tensor[indices]);
            }
          return out;
        }

        template <typename NumberType, int dim>
        SymmetricTensor<4, dim, NumberType>
        tensor_evaluate_optimized(
          const SymmetricTensor<4, dim, Expression> &symbol_tensor,
          const BatchOptimizer<NumberType> &         optimizer)
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
                out[indices] = optimizer.evaluate(symbol_tensor[indices]);
              }
          return out;
        }

        /**
         * End-point for recursive template function
         */
        template <typename NumberType, typename T>
        void
        register_functions(BatchOptimizer<NumberType> &optimizer,
                           const T &                   function)
        {
          optimizer.register_function(function);
        }

        template <typename NumberType, typename T, typename... Args>
        void
        register_functions(BatchOptimizer<NumberType> &optimizer,
                           const T &                   function,
                           const Args &... other_functions)
        {
          register_functions(optimizer, function);
          register_functions(optimizer, other_functions...);
        }

        template <typename ReturnType, typename Optimizer, typename = void>
        struct OptimizerHelper;

        template <typename ReturnType, typename Optimizer>
        struct OptimizerHelper<ReturnType,
                               Optimizer,
                               typename std::enable_if<std::is_same<
                                 ReturnType,
                                 typename Optimizer::ReturnType>::value>::type>
        {
          static void
          initialize(typename Optimizer::OptimizerType *optimizer,
                     const SymEngine::vec_basic &       symbol_vec,
                     const SymEngine::vec_basic &       dependent_variables,
                     const enum OptimizationFlags &     optimization_flags)
          {
            //            const bool use_symbolic_cse =
            //            use_symbolic_CSE(optimization_flags);
            //            optimizer->init(symbol_vec, dependent_variables,
            //            use_symbolic_cse);

            // Some optimizers don't have the same interface for
            // initialization, we filter them out through the specializations
            // of the Optimizer class
            Optimizer::initialize(*optimizer,
                                  symbol_vec,
                                  dependent_variables,
                                  optimization_flags);
          }

          static void
          substitute(typename Optimizer::OptimizerType *optimizer,
                     std::vector<ReturnType> &          outputs,
                     const std::vector<ReturnType> &    inputs)
          {
            optimizer->call(outputs.data(), inputs.data());
          }

          template <typename Stream>
          static void
          print(Stream &                           stream,
                typename Optimizer::OptimizerType *optimizer,
                const bool print_independent_symbols = false,
                const bool print_dependent_functions = false,
                const bool print_cse_reductions      = true)
          {
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
        struct OptimizerHelper<ReturnType,
                               Optimizer,
                               typename std::enable_if<!std::is_same<
                                 ReturnType,
                                 typename Optimizer::ReturnType>::value>::type>
        {
          static void
          initialize(typename Optimizer::OptimizerType *optimizer,
                     const SymEngine::vec_basic &       symbol_vec,
                     const SymEngine::vec_basic &       dependent_variables,
                     const enum OptimizationFlags &     optimization_flags)
          {
            const bool use_symbolic_cse = use_symbolic_CSE(optimization_flags);
            optimizer->init(symbol_vec, dependent_variables, use_symbolic_cse);
          }

          static void
          substitute(typename Optimizer::OptimizerType *optimizer,
                     std::vector<ReturnType> &          outputs,
                     const std::vector<ReturnType> &    inputs)
          {
            // Intermediate values to accommodate the difference in
            // value types.
            std::vector<typename Optimizer::ReturnType> int_outputs(
              outputs.size());
            std::vector<typename Optimizer::ReturnType> int_inputs(
              inputs.size());

            std::copy(inputs.begin(), inputs.end(), int_inputs.begin());
            optimizer->call(int_outputs.data(), int_inputs.data());
            std::copy(int_outputs.begin(), int_outputs.end(), outputs.begin());
          }

          template <typename Stream>
          static void
          print(Stream &                           stream,
                typename Optimizer::OptimizerType *optimizer,
                const bool                         print_cse_reductions = true,
                const bool print_independent_symbols                    = false,
                const bool print_dependent_functions                    = false)
          {
            optimizer->print(stream,
                             print_independent_symbols,
                             print_dependent_functions,
                             print_cse_reductions);
          }
        };

      } // namespace
    }   // namespace internal



    /**
     * A class that facilitates the optimization of symbol expressions.
     *
     * This expression will be optimized by this class; that is to say that
     * the code path taken to substitute the set of (independent) symbols
     * into a collection of (dependent) symbolic functions will be optimized
     * using a chosen approach.
     *
     * General usage:
     * @code
     * typedef NumberType double;
     * SymEngineWrappers::BatchOptimizer<NumberType> optimizer
     * (SymEngineWrappers::lambda, SymEngineWrappers::optimize_cse);
     *
     * optimizer.register_symbols();
     * optimizer.register_functions();
     *
     * optimizer.optimize(); // Expensive calls
     *
     * optimizer.substitute(); // Now quicker than dictionary substitution
     * optimizer.evaluate();
     * @endcode
     *
     * @tparam ReturnType The number type that is to be returned after
     * value substitution and evaluation.
     *
     * @warning This class is not thread-safe.
     *
     * @author Jean-Paul Pelteret, Isuru Fernando, 2017
     */
    template <typename ReturnType>
    class BatchOptimizer
    {
    public:
      /**
       * Default constructor
       *
       * By default, dictionary substitution will be selected when this
       * constructor is called. In order to select a specific optimization
       * approach, a call set_optimization_method() is necessary.
       */
      BatchOptimizer();

      /**
       * Constructor
       *
       * @param[in] method The optimization method that is to be employed.
       * @param[in] flags  The optimization flags that indicate which
       * expression manipulation mechanisms are to be employed.
       *
       * @note As the optimization method is fixed, a further call to
       * set_optimization_method() is not necessary and will result in an
       * error being thrown.
       *
       * @note In the case that the @p method is not implemented for the
       * required @p ReturnType, or the desired feature is not active,
       * then a safe default will be selected.
       */
      BatchOptimizer(
        const enum OptimizerType &    method,
        const enum OptimizationFlags &flags = OptimizationFlags::optimize_all);

      /**
       * Copy constructor
       *
       * The @p copy_initialized flag, which is set to <code>true</code> by default,
       * determines whether or not all of the optimized data is copied over from
       * the @p other optimizer instance. Only with the flag set to <code>false</code>
       * is it possible to re-optimize the data stored in this class with a
       * different optimization scheme.
       */
      BatchOptimizer(const BatchOptimizer &other/*,
                     const bool            copy_initialized = true*/);

      /**
       * Move constructor
       */
      BatchOptimizer(BatchOptimizer &&) = default;

      /**
       * Destructor
       */
      ~BatchOptimizer() = default;

      /**
       * @note In the case that the @p method is not implemented for the
       * required @p ReturnType, or the desired feature is not active,
       * then a safe default will be selected.
       */
      void
      set_optimization_method(const enum OptimizerType &    optimization_method,
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
       * CSE.
       */
      bool
      use_symbolic_CSE() const;

      /**
       * Returns a flag which indicates whether the optimize()
       * function has been called and the class is finalized.
       */
      bool
      optimized() const;

      /**
       * Returns a flag to indicate whether the substitute()
       * function has been called and if there are meaningful
       * values that will be returned upon evaluation().
       */
      bool
      values_substituted() const;

      /**
       * Print some information on state of the internal data
       * structures stored in the class.
       */
      template <typename Stream>
      void
      print(Stream &stream, const bool print_cse = true) const;


      /**
       * @name Independent variables
       */
      //@{

      /**
       * Register a collection of symbols that represents an independent
       * variable.
       * These symbols are stored as the <tt>key</tt> to the @p symbol_values map.
       */
      void
      register_symbols(const SD::types::substitution_map &symbol_values);

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
      register_symbols(const typename SD::types::symbol_vector &symbols);

      /**
       * Return the list of symbols that have been registered as independent
       * variables.
       */
      SD::types::symbol_vector
      get_independent_symbols(void) const;

      /**
       * The number of independent variables that this optimizer will recognize.
       * This is equal to the number of unique symbols passed to this class
       * instance through the register_symbols() function.
       */
      std::size_t
      n_independent_variables(void) const;

      //@}

      /**
       * @name Dependent variables
       */
      //@{

      /**
       * Register a symbolic that represents a dependent variable.
       */
      void
      register_function(const Expression &func);

      template <int rank, int dim>
      void
      register_function(const Tensor<rank, dim, Expression> &function_tensor);

      template <int rank, int dim>
      void
      register_function(
        const SymmetricTensor<rank, dim, Expression> &function_tensor);

      /**
       * Register a collection of symbols that represents dependent variables.
       */
      void
      register_functions(const typename SD::types::symbol_vector &functions);

      template <typename T, typename... Args>
      void
      register_functions(const T &function, const Args &... other_functions);

      const SD::types::symbol_vector &
      get_dependent_functions(void) const;

      /**
       * The number of dependent symbolic expressions that this optimizer
       * will optimize. This is equal to the number of unique symbolic functions
       * passed to this class instance through the register_functions()
       * function.
       */
      std::size_t
      n_dependent_variables(void) const;

      //@}

      /**
       * @name Optimization
       */
      //@{

      /**
       * Perform the optimization of all registered dependent functions using
       * the registered symbols.
       *
       * @note This may be a time-consuming process, but if the class instance is
       * retained throughout the coarse of a simulation then it need only be
       * performed once.
       *
       * @note This function ,which should only be called once per instance of this
       * class, finalizes the set of accepted independent symbols and dependent
       * functions that are recognized and used by the optimizer.
       */
      void
      optimize();

      //@}

      /**
       * @name Substitution / evaluation and data extraction
       */
      //@{

      /**
       * Perform batch substitution of all of the registered symbols
       * into the registered functions. The result is cached and can
       * be extracted by calls to evaluate().
       *
       * @note Calling substitute() again with a new set of
       * @p substitution_values overwrites any previously computed
       * results.
       */
      void
      substitute(const SD::types::substitution_map &substitution_values) const;

      /**
       * Perform batch substitution of all of the registered symbols
       * into the registered functions. The result is cached and can
       * be extracted by calls to evaluate().
       *
       * @note Calling substitute() again with a new set of
       * @p substitution_values overwrites any previously computed
       * results.
       */
      void
      substitute(const typename SD::types::symbol_vector &symbols,
                 const std::vector<ReturnType> &substitution_values) const;

      /**
       * Returns the result of a value substitution into the optimized
       * counterpart of all dependent functions. This function fetches those
       * cached values.
       *
       * These values were computed by substituting a @p substitution_values map
       * during substitute() call.
       */
      const std::vector<ReturnType> &
      get_evaluated_functions() const;

      /**
       * Returns the result of a value substitution into the optimized
       * counterpart
       * of @p func. This function fetches that cached value.
       *
       * This value was computed by substituting a @p substitution_values map
       * during substitute() call.
       */
      ReturnType
      evaluate(const Expression &func) const;

      /**
       * Returns the result of a tensor value substitution into the optimized
       * counterpart of @p funcs. This function fetches that cached tensor
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
       * counterpart of @p funcs. This function fetches that cached symmetric
       * tensor components.
       *
       * This value was computed by substituting a @p substitution_values map
       * during substitute() call.
       */
      template <int rank, int dim>
      SymmetricTensor<rank, dim, ReturnType>
      evaluate(const SymmetricTensor<rank, dim, Expression> &funcs) const;

      //@}

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
      SD::types::substitution_map independent_variables_symbols;

      /**
       * A set of symbolic expressions to be optimized. It is required that
       * the symbols on which these dependent functions be based are are
       * registered in the @p independent_variables_symbols map.
       */
      SD::types::symbol_vector dependent_variables_functions;

#  ifdef DEBUG
      /**
       * Check to see if a function is exactly equal to one of the logical
       * results of a differentiation operation.
       */
      bool
      is_valid_nonunique_dependent_variable(
        const SymEngine::RCP<const SymEngine::Basic> &func) const;
#  endif

      /**
       * The output of substituting symbolic values with floating point
       * values through the use of the @p optimizer.
       *
       * @p It is necessary to use this intermediate storage mechanism
       * to store the result of a substitution sweep as some optimizers
       * work on all symbolic expressions in a single batch. In this way
       * they can employ methods such as common subexpression elimination
       * to minimise the number of terms evaluated across all symbolic
       * functions.
       *
       * This variable is marked as mutable. This facilitates the substitution
       * functionality being used in logically constant "get_*" functions.
       */
      mutable std::vector<ReturnType> dependent_variables_output;

      /**
       * A map type used to indicate which dependent variable is associated
       * with which entry in the output vector.
       *
       * We use a custom comparator here because otherwise we can't use
       * std::map::find; it is sensitive to the some data from the
       * underlying SymEngine::Basic other than the value that it represents.
       */
      using map_dependent_expression_to_vector_entry_t =
        std::map<SymEngine::RCP<const SymEngine::Basic>,
                 std::size_t,
                 SymEngine::RCPBasicKeyLess>;

      /**
       * A map indicating which dependent variable is associated with which
       * entry in the output vector.
       */
      map_dependent_expression_to_vector_entry_t map_dep_ptr_basic_vec_entry;

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
       * This variable is marked as mutable. This facilitates the substitution
       * functionality being used in logically constant "get_*" functions.
       */
      mutable bool ready_for_value_extraction;

      /**
       * Register a single symbol that represents an dependent variable.
       */
      void
      register_scalar_function(
        const SymEngine::RCP<const SymEngine::Basic> &func);

      /**
       * Register a collection of symbols that represents an dependent
       * variables.
       */
      void
      register_vector_functions(const typename SD::types::symbol_vector &funcs);

      /**
       * Create an instance of the selected optimizer.
       */
      void
      create_optimizer(std::unique_ptr<SymEngine::Visitor> &optimizer);

      /**
       * Clone an instance of the selected optimizer.
       */
      void
      clone_optimizer(std::unique_ptr<SymEngine::Visitor> &new_optimizer) const;

      /**
       * Perform batch substitution of all of the registered symbols
       * into the registered functions. The result is cached and can
       * be extracted by calls to evaluate().
       *
       * @note Calling substitute() again with a new set of
       * @p substitution_values overwrites any previously computed
       * results.
       *
       * @warning When using this function is no mechanism to check that the ordering
       * of the @p substitution_values vector matches the internal ordering of
       * the registered symbols. This function is therefore typically used in
       * conjunction with the register_symbols() function that takes in a vector
       * of symbols. With this pair of functions to the class interface, the
       * management of symbol ordering is maintained by the user.
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
      stream << "Method? " << optimization_method() << "\n";
      stream << "Flags: " << optimization_flags() << "\n";
      stream << "Optimized? " << (optimized() ? "Yes" : "No") << "\n";
      stream << "Values substituted? " << values_substituted() << "\n\n";

      // Independent variables
      stream << "Symbols (" << n_independent_variables()
             << " independent variables):"
             << "\n";
      int cntr = 0;
      for (SD::types::substitution_map::const_iterator it =
             independent_variables_symbols.begin();
           it != independent_variables_symbols.end();
           ++it, ++cntr)
        {
          stream << cntr << ": " << it->first << "\n";
        }
      stream << "\n" << std::flush;

      // Dependent functions
      stream << "Functions (" << n_dependent_variables()
             << " dependent variables):"
             << "\n";
      cntr = 0;
      for (typename SD::types::symbol_vector::const_iterator it =
             dependent_variables_functions.begin();
           it != dependent_variables_functions.end();
           ++it, ++cntr)
        {
          stream << cntr << ": " << (*it) << "\n";
        }
      stream << "\n" << std::flush;

      // Common subexpression
      if (optimized() == true && use_symbolic_CSE() == true)
        {
          Assert(optimizer, ExcNotInitialized());
          const bool print_cse_reductions      = true;
          const bool print_independent_symbols = false;
          const bool print_dependent_functions = false;

          // TODO: Why does this not work with dynamic_cast, but rather requires
          // a static_cast?
          //       Dynamic_casts work perfectly fine in the other functions.
          //       Maybe this is related:
          //       https://stackoverflow.com/questions/590371/dynamic-cast-fails
          //          if (typename
          //          internal::DictionaryOptimizer<ReturnType>::OptimizerType*
          //              opt = dynamic_cast<typename
          //              internal::DictionaryOptimizer<ReturnType>::OptimizerType
          //              *> (optimizer.get()))
          if (optimization_method() == OptimizerType::dictionary)
            {
              internal::OptimizerHelper<
                ReturnType,
                internal::DictionaryOptimizer<ReturnType>>::
                print(stream,
                      static_cast<typename internal::DictionaryOptimizer<
                        ReturnType>::OptimizerType *>(optimizer.get()),
                      print_independent_symbols,
                      print_dependent_functions,
                      print_cse_reductions);
              stream << "\n" << std::flush;
            }
          //          else if (typename
          //          internal::LambdaOptimizer<ReturnType>::OptimizerType*
          //                   opt = dynamic_cast<typename
          //                   internal::LambdaOptimizer<ReturnType>::OptimizerType
          //                   *> (optimizer.get()))
          else if (optimization_method() == OptimizerType::lambda)
            {
              //              opt->print(stream, print_independent_symbols,
              //              print_dependent_functions, print_cse_reductions);
              internal::OptimizerHelper<ReturnType,
                                        internal::LambdaOptimizer<ReturnType>>::
                print(stream,
                      static_cast<typename internal::LambdaOptimizer<
                        ReturnType>::OptimizerType *>(optimizer.get()),
                      print_independent_symbols,
                      print_dependent_functions,
                      print_cse_reductions);
            }
#    ifdef HAVE_SYMENGINE_LLVM
          //          else if (typename
          //          internal::LLVMOptimizer<ReturnType>::OptimizerType*
          //                   opt = dynamic_cast<typename
          //                   internal::LLVMOptimizer<ReturnType>::OptimizerType
          //                   *> (optimizer.get()))
          else if (optimization_method() == OptimizerType::llvm)
            {
              internal::OptimizerHelper<ReturnType,
                                        internal::LLVMOptimizer<ReturnType>>::
                print(stream,
                      static_cast<typename internal::LLVMOptimizer<
                        ReturnType>::OptimizerType *>(optimizer.get()),
                      print_independent_symbols,
                      print_dependent_functions,
                      print_cse_reductions);
            }
#    endif
          else
            {
              AssertThrow(false, ExcMessage("Unknown optimizer type."));
            }
        }

      if (values_substituted())
        {
          stream << "Evaluated functions:"
                 << "\n";
          stream << std::flush;
          cntr = 0;
          for (typename std::vector<ReturnType>::const_iterator it =
                 dependent_variables_output.begin();
               it != dependent_variables_output.end();
               ++it, ++cntr)
            {
              stream << cntr << ": " << (*it) << "\n";
            }
          stream << "\n" << std::flush;
        }
    }


    namespace internal
    {
      template <int rank,
                int dim,
                template <int, int, typename> class TensorType>
      SD::types::symbol_vector
      unroll_to_RCP_vector(
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


      template <int dim>
      SD::types::symbol_vector
      unroll_to_RCP_vector(
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
        internal::unroll_to_RCP_vector(function_tensor));
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
        internal::unroll_to_RCP_vector(function_tensor));
    }


    template <typename ReturnType>
    template <typename T, typename... Args>
    void
    BatchOptimizer<ReturnType>::register_functions(
      const T &function,
      const Args &... other_functions)
    {
      internal::register_functions(*this, function);
      internal::register_functions(*this, other_functions...);
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

      return internal::tensor_evaluate_optimized(funcs, *this);
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

      return internal::tensor_evaluate_optimized(funcs, *this);
    }

#  endif // DOXYGEN

  } // namespace SD
} // namespace Differentiation


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SYMENGINE

#endif
