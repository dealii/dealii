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

#ifndef dealii_differentiation_sd_symengine_number_visitor_internal_h
#define dealii_differentiation_sd_symengine_number_visitor_internal_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SYMENGINE

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
// Low level
#  include <symengine/basic.h>
#  include <symengine/dict.h>
#  include <symengine/symengine_exception.h>
#  include <symengine/symengine_rcp.h>

// Visitor
#  include <symengine/visitor.h>

DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#  include <deal.II/base/exceptions.h>
#  include <deal.II/base/numbers.h>

#  include <deal.II/differentiation/sd/symengine_number_types.h>
#  include <deal.II/differentiation/sd/symengine_utilities.h>

#  include <boost/serialization/split_member.hpp>


DEAL_II_NAMESPACE_OPEN


namespace Differentiation
{
  namespace SD
  {
    namespace internal
    {
      /**
       * A class that implements common subexpression elimination
       * for dictionary visitor classes.
       *
       * It is intended that this class only be used in conjunction
       * with the DictionarySubstitutionVisitor.
       *
       * @author Jean-Paul Pelteret, Isuru Fernando, 2017, 2020
       */
      template <typename ReturnType, typename ExpressionType>
      class CSEDictionaryVisitor
      {
        using symbol_vector_pair =
          std::vector<std::pair<SD::Expression, SD::Expression>>;

      public:
        /*
         * Constructor.
         */
        CSEDictionaryVisitor() = default;

        /*
         * Destructor.
         */
        virtual ~CSEDictionaryVisitor() = default;

        /**
         * Initialize and perform common subexpression elimination.
         *
         * Here we build the reduced expressions for the @p dependent_functions
         * as well as a list of intermediate, repeating symbolic expressions
         * that are extracted @p dependent_functions. This operation leads to
         * the elimination of repeated expressions, so they only have to be
         * evaluated once.
         *
         * @param dependent_functions A vector of expressions that represent
         * the dependent functions that CSE will be performed on.
         */
        void
        init(const types::symbol_vector &dependent_functions);

        /**
         * Evaluate the (reduced) dependent functions using the common
         * subexpressions generated during the call to init().
         *
         * The @p output_values are the numerical result of substituting
         * each of the @p substitution_values for their corresponding entry
         * of the @p independent_symbols. Specifically, we first calculate
         * the numerical values of each reduced subexpression, and feed those
         * into a reduced equivalent of the dependent variables. The result of
         * this is then written into the @p output_values array.
         *
         * @param[out] output_values A pointer to the first element in an array
         * of type @p ReturnType. After this call, the underlying array will
         * hold the numerical result of substituting the @p substitution_values
         * into the pre-registered dependent functions.
         * @param[in]  independent_symbols A vector of symbols that represent
         * the independent variables that are arguments to the previously
         * defined dependent functions.
         * @param[in]  substitution_values A pointer to the first element in an
         * array of type @p ReturnType. Each entry in this array stores the
         * numerical value that an independent variable is to take for the
         * purpose of value substitution.
         *
         * @note It is expected that both the @p output_values and
         * @p substitution_values arrays be correctly dimensioned, as there is
         * no range checking performed on these data structures.
         * The @p output_values array should have the same number of elements as
         * the dependent functions first passed to this class in the call to
         * init(). Similarly, the @p substitution_values array should have the
         * same number of elements as the @p independent_symbols vector has.
         *
         * @note It is expected that there be a 1-1 correspondence between the
         * entries in @p independent_symbols and @p substitution_values. This
         * is not checked within this function, so it is up to the user to
         * ensure that the relationships between the independent variables and
         * their numerical value be correctly set up and maintained.
         */
        void
        call(ReturnType *                output_values,
             const types::symbol_vector &independent_symbols,
             const ReturnType *          substitution_values);

        /**
         * Write the data of this object to a stream for the purpose
         * of serialization.
         */
        template <class Archive>
        void
        save(Archive &archive, const unsigned int version) const;

        /**
         * Read the data for this object from a stream for the purpose
         * of serialization.
         */
        template <class Archive>
        void
        load(Archive &archive, const unsigned int version);

#  ifdef DOXYGEN
        /**
         * Write and read the data of this object from a stream for the purpose
         * of serialization.
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
         * Print all of the intermediate reduced expressions, as well as
         * the final reduced expressions for dependent variables to the
         * @p stream.
         */
        template <typename StreamType>
        void
        print(StreamType &stream) const;

        /**
         * Return a flag stating whether we've performed CSE or not.
         */
        bool
        executed() const;

        /**
         * The number of intermediate expressions that must be evaluated as
         * part of the collection of common subexpressions.
         */
        unsigned int
        n_intermediate_expressions() const;

        /**
         * The size of the final set of reduced expressions.
         */
        unsigned int
        n_reduced_expressions() const;

      protected:
        /**
         * Initialize and perform common subexpression elimination.
         *
         * This function performs the same action as the other init() function,
         * except that it works with native SymEngine data types;
         * the single argument is  a vector of
         * `SymEngine::RCP<const SymEngine::Basic>`.
         *
         * @param dependent_functions A vector of expressions that represent
         * the dependent functions that CSE will be performed on.
         */
        void
        init(const SymEngine::vec_basic &dependent_functions);

        /**
         * Evaluate the (reduced) dependent functions using the common
         * subexpressions generated during the call to init().
         *
         * This function performs the same action as the other init() function,
         * except that it works with native SymEngine data types.
         * The @p independent_symbols is represented by a vector of
         * `SymEngine::RCP<const SymEngine::Basic>`.
         *
         * @param[out] output_values A pointer to the first element in an array
         * of type @p ReturnType. After this call, the underlying array will
         * hold the numerical result of substituting the @p substitution_values
         * into the pre-registered dependent functions.
         * @param[in]  independent_symbols A vector of symbols that represent
         * the independent variables that are arguments to the previously
         * defined dependent functions.
         * @param[in]  substitution_values A pointer to the first element in an
         * array of type @p ReturnType. Each entry in this array stores the
         * numerical value that an independent variable is to take for the
         * purpose of value substitution.
         *
         * @note The caveats described in the documentation of the other call()
         * function apply here as well.
         */
        void
        call(ReturnType *                output_values,
             const SymEngine::vec_basic &independent_symbols,
             const ReturnType *          substitution_values);

      private:
        // Note: It would be more efficient to store this data in native
        // SymEngine types, as it would prevent some copying of the data
        // structures. However, this makes serialization more difficult,
        // so we use our own serializable types instead, and lose a bit
        // of efficiency.

        /**
         * Intermediate symbols and their definition.
         */
        symbol_vector_pair intermediate_symbols_exprs;

        /**
         * Final reduced expressions.
         */
        types::symbol_vector reduced_exprs;
      };



      /**
       * A class to perform dictionary-based substitution as if it were an
       * optimizer of the "lambda" or "LLVM" variety.
       *
       * This class is only really useful to assist in the easy switching
       * between different optimizers and, more importantly, for integrating
       * CSE into a dictionary substitution scheme. It is therefore only
       * intended to be created and used by a BatchOptimizer.
       *
       * @author Jean-Paul Pelteret, Isuru Fernando, 2017, 2020
       */
      template <typename ReturnType, typename ExpressionType>
      class DictionarySubstitutionVisitor
        : public SymEngine::BaseVisitor<
            DictionarySubstitutionVisitor<ReturnType, ExpressionType>>
      {
      public:
        /*
         * Constructor.
         */
        DictionarySubstitutionVisitor() = default;

        /*
         * Destructor.
         */
        virtual ~DictionarySubstitutionVisitor() override = default;

        /**
         * Initialization, and registration of the independent and dependent
         * variables.
         *
         * This variation of the initialization function registers a single
         * dependent expression. If the @p use_cse is set to `true`, then common
         * subexpression elimination is also performed at the time of
         * initialization.
         *
         * @param independent_symbols A vector of symbols that represent
         * the independent variables that are arguments to the
         * @p dependent_function.
         * @param dependent_function A single symbolic expression that
         * represents a dependent variable.
         * @param use_cse A flag to indicate whether or not to use common
         * subexpression elimination.
         *
         * @note After this function is called, no further registration of
         * dependent functions can be performed. If it is desired that multiple
         * dependent expressions be registered, then the other variant of this
         * function that takes in a vector of dependent expressions should be
         * used.
         */
        void
        init(const types::symbol_vector &independent_symbols,
             const Expression &          dependent_function,
             const bool                  use_cse = false);

        /**
         * Initialization, and registration of the independent and dependent
         * variables.
         *
         * This function performs the same action as the other init() function,
         * described above, except that it works with native SymEngine data
         * types. The @p independent_symbols are represented by a vector of
         * `SymEngine::RCP<const SymEngine::Basic>`, and the
         * @p dependent_function is of type
         * `SymEngine::RCP<const SymEngine::Basic>`.
         *
         * @param independent_symbols A vector of symbols that represent
         * the independent variables that are arguments to the
         * @p dependent_function.
         * @param dependent_function A single symbolic expression that
         * represents a dependent variable.
         * @param use_cse A flag to indicate whether or not to use common
         * subexpression elimination.
         *
         * @note The caveats described in the documentation of the other init()
         * function apply here as well.
         */
        // The following definition is required due to base class CRTP.
        void
        init(const SymEngine::vec_basic &independent_symbols,
             const SymEngine::Basic &    dependent_function,
             const bool                  use_cse = false);

        /**
         * Initialization, and registration of the independent and dependent
         * variables.
         *
         * This variation of the initialization function registers a vector of
         * dependent expressions. If the @p use_cse is set to `true`, then
         * common subexpression elimination is also performed at the time of
         * initialization.
         *
         * @param independent_symbols A vector of symbols that represent
         * the independent variables that are arguments to the
         * @p dependent_functions.
         * @param dependent_functions A vector of symbolic expressions that
         * represent the dependent variables.
         * @param use_cse A flag to indicate whether or not to use common
         * subexpression elimination.
         *
         * @note After this function is called, no further registration of
         * dependent functions can be performed.
         */
        void
        init(const types::symbol_vector &independent_symbols,
             const types::symbol_vector &dependent_functions,
             const bool                  use_cse = false);


        /**
         * Initialization, and registration of the independent and dependent
         * variables.
         *
         * This function performs the same action as the other init() function,
         * described above, except that it works with native SymEngine data
         * types. Both the @p independent_symbols and @p dependent_functions
         * are represented by a vector of
         * `SymEngine::RCP<const SymEngine::Basic>`.
         *
         * @param independent_symbols A vector of symbols that represent
         * the independent variables that are arguments to the
         * @p dependent_functions.
         * @param dependent_functions A vector of symbolic expressions that
         * represent the dependent variables.
         * @param use_cse A flag to indicate whether or not to use common
         * subexpression elimination.
         *
         * @note The caveats described in the documentation of the other init()
         * function apply here as well.
         */
        // The following definition is required due to base class CRTP.
        void
        init(const SymEngine::vec_basic &independent_symbols,
             const SymEngine::vec_basic &dependent_functions,
             const bool                  use_cse = false);

        /**
         * Evaluate the dependent functions that were registered at
         * initializtion time.
         *
         * The @p output_values are the numerical result of substituting
         * each of the @p substitution_values for their corresponding entry
         * of the pre-registered independent variables.
         *
         * @param[out] output_values A pointer to the first element in an array
         * of type @p ReturnType. After this call, the underlying array will
         * hold the numerical result of substituting the @p substitution_values
         * into the pre-registered dependent functions.
         * @param[in]  substitution_values A pointer to the first element in an
         * array of type @p ReturnType. Each entry in this array stores the
         * numerical value that an independent variable is to take for the
         * purpose of value substitution.
         *
         * @note It is expected that both the @p output_values and
         * @p substitution_values arrays be correctly dimensioned, as there is
         * no range checking performed on these data structures.
         * The @p output_values array should have the same number of elements as
         * the dependent functions first passed to this class in the call to
         * init(). Similarly, the @p substitution_values array should have the
         * same number of elements as the there were independent variables
         * that were registered at the time of initialization.
         *
         * @note It is expected that there be a 1-1 correspondence between the
         * entries in independent variables and @p substitution_values. This
         * is not checked within this function, so it is up to the user to
         * ensure that the relationships between the independent variables and
         * their numerical value be correctly set up and maintained.
         */
        void
        call(ReturnType *output_values, const ReturnType *substitution_values);

        /**
         * Evaluate the dependent function that were registered at
         * initializtion time.
         *
         * The purpose of this function is the same as the other call()
         * functions, but
         *
         * @param substitution_values A vector that stores the numerical values
         * that each independent variable that was registered with the class
         * instance is to take for the purpose of value substitution.
         * @return ReturnType The numerical value associated with the single
         * dependent function that is registered with this class.
         *
         * @note It is expected that there be a 1-1 correspondence between the
         * entries in independent variables and @p substitution_values.
         */
        // The following definition is required due to base class CRTP.
        ReturnType
        call(const std::vector<ReturnType> &substitution_values);

        /**
         * Write the data of this object to a stream for the purpose
         * of serialization.
         */
        template <class Archive>
        void
        save(Archive &archive, const unsigned int version) const;

        /**
         * Read the data for this object from a stream for the purpose
         * of serialization.
         */
        template <class Archive>
        void
        load(Archive &archive, const unsigned int version);

#  ifdef DOXYGEN
        /**
         * Write and read the data of this object from a stream for the purpose
         * of serialization.
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
         * Print some information on state of the internal data
         * structures stored in the class.
         *
         * @tparam Stream The type for the output stream.
         * @param stream The output stream to print to.
         * @param print_independent_symbols A flag to indicate if the independent
         * variables should be outputted to the @p stream.
         * @param print_dependent_functions A flag to indicate if the dependent
         * expressions should be outputted to the @p stream.
         * @param print_cse_reductions A flag to indicate whether or not all
         * common subexpressions should be printed to the @p stream.
         */
        template <typename StreamType>
        void
        print(StreamType &stream,
              const bool  print_independent_symbols = false,
              const bool  print_dependent_functions = false,
              const bool  print_cse_reductions      = false) const;

#  ifndef DOXYGEN
        // The following definitions are required due to base class CRTP.
        // Since these are not used, and therefore not important to
        // understand, we'll define them in the most concise manner possible.
        // We also won't bother to document their existence, since they cannot
        // be used.
#    define IMPLEMENT_DSV_BVISIT(Argument)       \
      void bvisit(const Argument &)              \
      {                                          \
        AssertThrow(false, ExcNotImplemented()); \
      }

        IMPLEMENT_DSV_BVISIT(SymEngine::Basic)
        IMPLEMENT_DSV_BVISIT(SymEngine::Symbol)
        IMPLEMENT_DSV_BVISIT(SymEngine::Constant)
        IMPLEMENT_DSV_BVISIT(SymEngine::Integer)
        IMPLEMENT_DSV_BVISIT(SymEngine::Rational)
        IMPLEMENT_DSV_BVISIT(SymEngine::RealDouble)
        IMPLEMENT_DSV_BVISIT(SymEngine::ComplexDouble)
        IMPLEMENT_DSV_BVISIT(SymEngine::Add)
        IMPLEMENT_DSV_BVISIT(SymEngine::Mul)
        IMPLEMENT_DSV_BVISIT(SymEngine::Pow)
        IMPLEMENT_DSV_BVISIT(SymEngine::Log)
        IMPLEMENT_DSV_BVISIT(SymEngine::Sin)
        IMPLEMENT_DSV_BVISIT(SymEngine::Cos)
        IMPLEMENT_DSV_BVISIT(SymEngine::Tan)
        IMPLEMENT_DSV_BVISIT(SymEngine::Csc)
        IMPLEMENT_DSV_BVISIT(SymEngine::Sec)
        IMPLEMENT_DSV_BVISIT(SymEngine::Cot)
        IMPLEMENT_DSV_BVISIT(SymEngine::ASin)
        IMPLEMENT_DSV_BVISIT(SymEngine::ACos)
        IMPLEMENT_DSV_BVISIT(SymEngine::ATan)
        IMPLEMENT_DSV_BVISIT(SymEngine::ATan2)
        IMPLEMENT_DSV_BVISIT(SymEngine::ACsc)
        IMPLEMENT_DSV_BVISIT(SymEngine::ASec)
        IMPLEMENT_DSV_BVISIT(SymEngine::ACot)
        IMPLEMENT_DSV_BVISIT(SymEngine::Sinh)
        IMPLEMENT_DSV_BVISIT(SymEngine::Cosh)
        IMPLEMENT_DSV_BVISIT(SymEngine::Tanh)
        IMPLEMENT_DSV_BVISIT(SymEngine::Csch)
        IMPLEMENT_DSV_BVISIT(SymEngine::Sech)
        IMPLEMENT_DSV_BVISIT(SymEngine::Coth)
        IMPLEMENT_DSV_BVISIT(SymEngine::ASinh)
        IMPLEMENT_DSV_BVISIT(SymEngine::ACosh)
        IMPLEMENT_DSV_BVISIT(SymEngine::ATanh)
        IMPLEMENT_DSV_BVISIT(SymEngine::ACsch)
        IMPLEMENT_DSV_BVISIT(SymEngine::ACoth)
        IMPLEMENT_DSV_BVISIT(SymEngine::ASech)
        IMPLEMENT_DSV_BVISIT(SymEngine::Abs)
        IMPLEMENT_DSV_BVISIT(SymEngine::Gamma)
        IMPLEMENT_DSV_BVISIT(SymEngine::LogGamma)
        IMPLEMENT_DSV_BVISIT(SymEngine::Erf)
        IMPLEMENT_DSV_BVISIT(SymEngine::Erfc)
        IMPLEMENT_DSV_BVISIT(SymEngine::Max)
        IMPLEMENT_DSV_BVISIT(SymEngine::Min)

#    undef IMPLEMENT_DSV_BVISIT
#  endif // DOXYGEN

      private:
        // Note: It would be more efficient to store this data in native
        // SymEngine types, as it would prevent some copying of the data
        // structures. However, this makes serialization more difficult,
        // so we use our own serializable types instead, and lose a bit
        // of efficiency.

        /**
         * A vector of symbols that represent the independent variables.
         */
        SD::types::symbol_vector independent_symbols;

        /**
         * A vector of expressions that represent dependent functions.
         */
        SD::types::symbol_vector dependent_functions;

        /**
         * A data structure that may be used to invoke common subexpression
         * elimination on the dependent functions, with the aim to decrease
         * the time taken to evaluate them.
         */
        CSEDictionaryVisitor<ReturnType, ExpressionType> cse;
      };



      /* ------------------ inline and template functions ------------------ */


#  ifndef DOXYGEN

      /* -------------- CommonSubexpressionEliminationVisitor -------------- */


      template <typename ReturnType, typename ExpressionType>
      void
      CSEDictionaryVisitor<ReturnType, ExpressionType>::init(
        const SD::types::symbol_vector &dependent_functions)
      {
        init(Utilities::convert_expression_vector_to_basic_vector(
          dependent_functions));
      }



      template <typename ReturnType, typename ExpressionType>
      void
      CSEDictionaryVisitor<ReturnType, ExpressionType>::init(
        const SymEngine::vec_basic &dependent_functions)
      {
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
        SymEngine::vec_pair  se_replacements;
        SymEngine::vec_basic se_reduced_exprs;
        SymEngine::cse(se_replacements, se_reduced_exprs, dependent_functions);

        intermediate_symbols_exprs =
          Utilities::convert_basic_pair_vector_to_expression_pair_vector(
            se_replacements);
        reduced_exprs = Utilities::convert_basic_vector_to_expression_vector(
          se_reduced_exprs);
      }



      template <typename ReturnType, typename ExpressionType>
      void
      CSEDictionaryVisitor<ReturnType, ExpressionType>::call(
        ReturnType *                    output_values,
        const SD::types::symbol_vector &independent_symbols,
        const ReturnType *              substitution_values)
      {
        call(output_values,
             Utilities::convert_expression_vector_to_basic_vector(
               independent_symbols),
             substitution_values);
      }



      template <typename ReturnType, typename ExpressionType>
      void
      CSEDictionaryVisitor<ReturnType, ExpressionType>::call(
        ReturnType *                output_values,
        const SymEngine::vec_basic &independent_symbols,
        const ReturnType *          substitution_values)
      {
        Assert(n_reduced_expressions() > 0, ExcInternalError());

        // First we add the input values into the substitution map...
        SymEngine::map_basic_basic substitution_value_map;
        for (unsigned i = 0; i < independent_symbols.size(); ++i)
          substitution_value_map[independent_symbols[i]] =
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



      template <typename ReturnType, typename ExpressionType>
      template <class Archive>
      void
      CSEDictionaryVisitor<ReturnType, ExpressionType>::save(
        Archive &ar,
        const unsigned int /*version*/) const
      {
        // The reduced expressions depend on the intermediate expressions,
        // so we serialize the latter before the former.
        ar &intermediate_symbols_exprs;
        ar &reduced_exprs;
      }



      template <typename ReturnType, typename ExpressionType>
      template <class Archive>
      void
      CSEDictionaryVisitor<ReturnType, ExpressionType>::load(
        Archive &ar,
        const unsigned int /*version*/)
      {
        Assert(intermediate_symbols_exprs.empty(), ExcInternalError());
        Assert(reduced_exprs.empty(), ExcInternalError());

        // The reduced expressions depend on the intermediate expressions,
        // so we deserialize the latter before the former.
        ar &intermediate_symbols_exprs;
        ar &reduced_exprs;
      }



      template <typename ReturnType, typename ExpressionType>
      template <typename StreamType>
      void
      CSEDictionaryVisitor<ReturnType, ExpressionType>::print(
        StreamType &stream) const
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



      template <typename ReturnType, typename ExpressionType>
      bool
      CSEDictionaryVisitor<ReturnType, ExpressionType>::executed() const
      {
        // For dictionary substitution, the CSE algorithm moves
        // ownership of the dependent function expression definition
        // to the entries in reduced_exprs. So its size thus determines
        // whether CSE has been executed or not.
        return (n_reduced_expressions() > 0) ||
               (n_intermediate_expressions() > 0);
      }



      template <typename ReturnType, typename ExpressionType>
      unsigned int
      CSEDictionaryVisitor<ReturnType,
                           ExpressionType>::n_intermediate_expressions() const
      {
        return intermediate_symbols_exprs.size();
      }



      template <typename ReturnType, typename ExpressionType>
      unsigned int
      CSEDictionaryVisitor<ReturnType, ExpressionType>::n_reduced_expressions()
        const
      {
        return reduced_exprs.size();
      }



      /* ------------------ DictionarySubstitutionVisitor ------------------ */


      template <typename ReturnType, typename ExpressionType>
      void
      DictionarySubstitutionVisitor<ReturnType, ExpressionType>::init(
        const types::symbol_vector &inputs,
        const SD::Expression &      output,
        const bool                  use_cse)
      {
        init(inputs, types::symbol_vector{output}, use_cse);
      }



      template <typename ReturnType, typename ExpressionType>
      void
      DictionarySubstitutionVisitor<ReturnType, ExpressionType>::init(
        const SymEngine::vec_basic &inputs,
        const SymEngine::Basic &    output,
        const bool                  use_cse)
      {
        init(Utilities::convert_basic_vector_to_expression_vector(inputs),
             SD::Expression(output.rcp_from_this()),
             use_cse);
      }



      template <typename ReturnType, typename ExpressionType>
      void
      DictionarySubstitutionVisitor<ReturnType, ExpressionType>::init(
        const SymEngine::vec_basic &inputs,
        const SymEngine::vec_basic &outputs,
        const bool                  use_cse)
      {
        init(Utilities::convert_basic_vector_to_expression_vector(inputs),
             Utilities::convert_basic_vector_to_expression_vector(outputs),
             use_cse);
      }



      template <typename ReturnType, typename ExpressionType>
      void
      DictionarySubstitutionVisitor<ReturnType, ExpressionType>::init(
        const types::symbol_vector &inputs,
        const types::symbol_vector &outputs,
        const bool                  use_cse)
      {
        independent_symbols.clear();
        dependent_functions.clear();

        independent_symbols = inputs;

        // Perform common subexpression elimination if requested
        // Note: After this is done, the results produced by
        // dependent_functions and cse.reduced_exprs should be
        // the same. We could keep the former so that we can print
        // out the original expressions if we wish to do so.
        if (use_cse == false)
          dependent_functions = outputs;
        else
          {
            cse.init(outputs);
          }
      }



      template <typename ReturnType, typename ExpressionType>
      ReturnType
      DictionarySubstitutionVisitor<ReturnType, ExpressionType>::call(
        const std::vector<ReturnType> &substitution_values)
      {
        Assert(
          dependent_functions.size() == 1,
          ExcMessage(
            "Cannot use this call function when more than one symbolic expression is to be evaluated."));
        Assert(
          substitution_values.size() == independent_symbols.size(),
          ExcMessage(
            "Input substitution vector does not match size of symbol vector."));

        ReturnType out = dealii::internal::NumberType<ReturnType>::value(0.0);
        call(&out, substitution_values.data());
        return out;
      }



      template <typename ReturnType, typename ExpressionType>
      void
      DictionarySubstitutionVisitor<ReturnType, ExpressionType>::call(
        ReturnType *      output_values,
        const ReturnType *substitution_values)
      {
        // Check to see if CSE has been performed
        if (cse.executed())
          {
            cse.call(output_values, independent_symbols, substitution_values);
          }
        else
          {
            // Build a substitution map.
            SymEngine::map_basic_basic substitution_value_map;
            for (unsigned i = 0; i < independent_symbols.size(); ++i)
              substitution_value_map[independent_symbols[i]] =
                static_cast<const SymEngine::RCP<const SymEngine::Basic> &>(
                  ExpressionType(substitution_values[i]));

            // Since we don't know how to definitively evaluate the
            // input number type, we create a generic Expression
            // with the given symbolic expression and ask it to perform
            // substitution and evaluation for us.
            Assert(dependent_functions.size() > 0, ExcInternalError());
            for (unsigned i = 0; i < dependent_functions.size(); ++i)
              output_values[i] =
                ExpressionType(dependent_functions[i])
                  .template substitute_and_evaluate<ReturnType>(
                    substitution_value_map);
          }
      }



      template <typename ReturnType, typename ExpressionType>
      template <class Archive>
      void
      DictionarySubstitutionVisitor<ReturnType, ExpressionType>::save(
        Archive &          ar,
        const unsigned int version) const
      {
        // Add some dynamic information to determine if CSE has been used,
        // without relying on the CSE class when deserializing.
        // const bool used_cse = cse.executed();
        // ar &used_cse;

        // CSE and dependent variables both require the independent
        // symbols, so we serialize them first. The dependent variables
        // might depend on the outcome of CSE, so we have to serialize
        // them last.
        ar &independent_symbols;
        cse.save(ar, version);
        ar &dependent_functions;
      }



      template <typename ReturnType, typename ExpressionType>
      template <class Archive>
      void
      DictionarySubstitutionVisitor<ReturnType, ExpressionType>::load(
        Archive &          ar,
        const unsigned int version)
      {
        Assert(cse.executed() == false, ExcInternalError());
        Assert(cse.n_intermediate_expressions() == 0, ExcInternalError());
        Assert(cse.n_reduced_expressions() == 0, ExcInternalError());

        // CSE and dependent variables both require the independent
        // symbols, so we deserialize them first. The dependent variables
        // might depend on the outcome of CSE, so we have to deserialize
        // them last.
        ar &independent_symbols;
        cse.load(ar, version);
        ar &dependent_functions;
      }



      template <typename ReturnType, typename ExpressionType>
      template <typename StreamType>
      void
      DictionarySubstitutionVisitor<ReturnType, ExpressionType>::print(
        StreamType &stream,
        const bool  print_independent_symbols,
        const bool  print_dependent_functions,
        const bool  print_cse_reductions) const
      {
        if (print_independent_symbols)
          {
            stream << "Independent variables: \n";
            for (unsigned i = 0; i < independent_symbols.size(); ++i)
              stream << "  " << i << ": " << independent_symbols[i] << "\n";

            stream << std::flush;
          }

        // Check to see if CSE has been performed
        if (print_cse_reductions && cse.executed())
          {
            cse.print(stream);
          }
        else
          {
            Assert(dependent_functions.size() > 0, ExcInternalError());

            if (print_dependent_functions)
              {
                stream << "Dependent variables: \n";
                for (unsigned i = 0; i < dependent_functions.size(); ++i)
                  stream << "  " << i << dependent_functions[i] << "\n";

                stream << std::flush;
              }
          }
      }

#  endif // DOXYGEN

    } // namespace internal
  }   // namespace SD
} // namespace Differentiation


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SYMENGINE

#endif // dealii_differentiation_sd_symengine_number_visitor_internal_h
