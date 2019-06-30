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

#ifndef dealii_differentiation_sd_symengine_helpers_h
#define dealii_differentiation_sd_symengine_helpers_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SYMENGINE

#  include <deal.II/base/numbers.h>
#  include <deal.II/base/symmetric_tensor.h>
#  include <deal.II/base/tensor.h>

#  include <deal.II/differentiation/sd/symengine_number_types.h>
#  include <deal.II/differentiation/sd/symengine_optimizer.h>
#  include <deal.II/differentiation/sd/symengine_product_types.h>
#  include <deal.II/differentiation/sd/symengine_scalar_operations.h>
#  include <deal.II/differentiation/sd/symengine_tensor_operations.h>

#  include <deal.II/fe/fe_values.h>
#  include <deal.II/fe/fe_values_extractors.h>

#  include <deal.II/lac/full_matrix.h>
#  include <deal.II/lac/vector.h>

#  include <algorithm>
#  include <iostream>
#  include <iterator>
#  include <numeric>
#  include <set>

DEAL_II_NAMESPACE_OPEN

namespace Differentiation
{
  namespace SD
  {
    namespace internal
    {
      namespace
      {
        template <int spacedim>
        struct CurlType;

        template <>
        struct CurlType<1>
        {
          typedef Tensor<1, 1, Expression> type;

          static type
          make_symbolic_expression(const std::string &symbol)
          {
            return make_tensor_of_symbols<1, 1>(symbol);
          }
        };

        template <>
        struct CurlType<2>
        {
          typedef Tensor<1, 1, Expression> type;

          static type
          make_symbolic_expression(const std::string &symbol)
          {
            return make_tensor_of_symbols<1, 1>(symbol);
          }
        };

        template <>
        struct CurlType<3>
        {
          typedef Tensor<1, 3, Expression> type;

          static type
          make_symbolic_expression(const std::string &symbol)
          {
            return make_tensor_of_symbols<1, 3>(symbol);
          }
        };

        /**
         * A small struct to help create symbolic tensor
         * names based on an input FEValuesViews.
         */
        template <int spacedim,
                  typename SDNumberType,
                  typename FEValuesViewsType>
        struct SymbolCreator;

        template <int spacedim, typename SDNumberType>
        struct SymbolCreator<spacedim,
                             SDNumberType,
                             FEValuesViews::Scalar<spacedim>>
        {
          typedef SDNumberType                      value_type;
          typedef Tensor<1, spacedim, SDNumberType> gradient_type;
          typedef Tensor<2, spacedim, SDNumberType> hessian_type;
          typedef Tensor<3, spacedim, SDNumberType> third_derivative_type;

          static value_type
          value(const std::string &symbol)
          {
            return make_symbol(symbol);
          }

          static gradient_type
          gradient(const std::string &symbol)
          {
            return make_tensor_of_symbols<1, spacedim>(symbol);
          }

          static hessian_type
          hessian(const std::string &symbol)
          {
            return make_tensor_of_symbols<2, spacedim>(symbol);
          }

          static third_derivative_type
          third_derivative(const std::string &symbol)
          {
            return make_tensor_of_symbols<3, spacedim>(symbol);
          }
        };

        template <int spacedim, typename SDNumberType>
        struct SymbolCreator<spacedim,
                             SDNumberType,
                             FEValuesViews::Vector<spacedim>>
        {
          typedef Tensor<1, spacedim, SDNumberType> value_type;
          typedef Tensor<2, spacedim, SDNumberType> gradient_type;
          typedef SymmetricTensor<2, spacedim, SDNumberType>
                                                    symmetric_gradient_type;
          typedef SDNumberType                      divergence_type;
          typedef typename CurlType<spacedim>::type curl_type;
          typedef Tensor<3, spacedim, SDNumberType> hessian_type;
          typedef Tensor<4, spacedim, SDNumberType> third_derivative_type;

          static value_type
          value(const std::string &symbol)
          {
            return make_tensor_of_symbols<1, spacedim>(symbol);
          }

          static gradient_type
          gradient(const std::string &symbol)
          {
            return make_tensor_of_symbols<2, spacedim>(symbol);
          }

          static symmetric_gradient_type
          symmetric_gradient(const std::string &symbol)
          {
            return make_symmetric_tensor_of_symbols<2, spacedim>(symbol);
          }

          static divergence_type
          divergence(const std::string &symbol)
          {
            return make_symbol(symbol);
          }

          static curl_type
          curl(const std::string &symbol)
          {
            return CurlType<spacedim>::make_symbolic_expression(symbol);
          }

          static hessian_type
          hessian(const std::string &symbol)
          {
            return make_tensor_of_symbols<3, spacedim>(symbol);
          }

          static third_derivative_type
          third_derivative(const std::string &symbol)
          {
            return make_tensor_of_symbols<4, spacedim>(symbol);
          }
        };

        //      // Instantiation
        //      template struct SymbolCreator<2, FEValuesViews::Scalar<2> >;
        //      template struct SymbolCreator<3, FEValuesViews::Scalar<3> >;
        //      template struct SymbolCreator<2, FEValuesViews::Vector<2> >;
        //      template struct SymbolCreator<3, FEValuesViews::Vector<3> >;

        // TODO: SymbolCreator<Tensor< 2, dim, spacedim > >
        // TODO: SymbolCreator<SymmetricTensor< 2, dim, spacedim > >
      } // namespace
    }   // namespace internal


    /**
     * The base class for symbolic assembly helpers that are used to
     * symbolically define, and then evaluate, the residual and tangent
     * contributions derived from a single quadrature point.
     *
     * This data structure will keep all of the symbolic data to be
     * fed into the associated optimizer. It will also provide a mechanism
     * to define all of the integration and finite element information
     * symbolically.
     *
     * @author Jean-Paul Pelteret, 2017
     */
    template <int dim,
              typename NumberType,
              typename ExpressionType = SD::Expression>
    class CellLevelBase
    {
    public:
      /**
       * The scalar type of the computed residual and linearization.
       */
      typedef NumberType scalar_type;

      /**
       * Type definition for the symbolic differentiation number type that is
       * used in all computations.
       */
      typedef ExpressionType sd_type;

      // Forward declarations
      struct FESymbolicNames;
      struct AdditionalData;

      /**
       * Default constructor
       */
      CellLevelBase(const unsigned int    n_independent_variables,
                    const unsigned int    n_dependent_variables,
                    const bool            symmetric_system,
                    const AdditionalData &additional_data = AdditionalData());

      /**
       * Copy constructor
       *
       * This constructor is deleted because the optimizer cannot
       * be copied.
       */
      CellLevelBase(const CellLevelBase &) = delete;

      /**
       * Move constructor
       *
       * The optimizer is not copyable so instead when necesasry
       * we invoke a move constructor to transfer ownership of
       * the pointer to the optimizer to another instance of this
       * object.
       */
      CellLevelBase(CellLevelBase &&rhs);

      /**
       * @name Interrogation of internal information
       */
      //@{

      /**
       * Returns the number of independent variables that this object expects to
       * work with.
       */
      std::size_t
      n_independent_variables() const;

      /**
       * Returns the number of dependent variables that this object expects to
       * operate on.
       */
      std::size_t
      n_dependent_variables() const;

      /**
       * Print the symbolic expressions for the
       * residual and tangent contributions from
       * a quadrature point
       */
      template <typename Stream>
      void
      print(Stream &   stream,
            const bool print_residual  = true,
            const bool print_tangent   = true,
            const bool print_optimizer = false) const;

      //@}

      /**
       * @name Independent variables
       */
      //@{

      /**
       * Degree-of-freedom values
       */
      std::vector<sd_type>
      dof_values() const;

      /**
       * Integration factor
       */
      sd_type
      JxW() const;

      //@}

      /**
       * @name Independent variables: Values
       */
      //@{

      /**
       * Return the shape functions values associated with a specific
       * degree-of-freedom.
       *
       * @note When using this function it is typically assumed that
       * the space from which the test function and trial solution
       * derive are the same. If it is necessary to make a distinction
       * between these, then this is provided by
       * shape_function_value_test() and shape_function_value_trial().
       */
      template <typename FEValuesType, typename FEExtractorType>
      auto
      shape_function_value(const FEValuesType &   fe_values,
                           const FEExtractorType &extractor,
                           const std::string &    field_symbol,
                           const unsigned int &   function_no) const ->
        typename internal::SymbolCreator<
          dim,
          sd_type,
          typename std::decay<decltype(fe_values[extractor])>::type>::
          value_type;

      /**
       * Return the shape function values associated with a specific
       * component of the test function.
       */
      template <typename FEValuesType, typename FEExtractorType>
      auto
      shape_function_value_test(const FEValuesType &   fe_values,
                                const FEExtractorType &extractor,
                                const std::string &    field_symbol,
                                const unsigned int &   function_no) const ->
        typename internal::SymbolCreator<
          dim,
          sd_type,
          typename std::decay<decltype(fe_values[extractor])>::type>::
          value_type;

      /**
       * Return the shape function values associated with a specific
       * trial solution (corresponding to a degree-of-freedom).
       */
      template <typename FEValuesType, typename FEExtractorType>
      auto
      shape_function_value_trial(const FEValuesType &   fe_values,
                                 const FEExtractorType &extractor,
                                 const std::string &    field_symbol,
                                 const unsigned int &   function_no) const ->
        typename internal::SymbolCreator<
          dim,
          sd_type,
          typename std::decay<decltype(fe_values[extractor])>::type>::
          value_type;

      /**
       * Return a vector of shape function values associated with each
       * degree-of-freedom.
       *
       * @note When using this function it is typically assumed that
       * the space from which the test function and trial solution
       * derive are the same. If it is necessary to make a distinction
       * between these, then this is provided by
       * shape_functions_value_test() and shape_functions_value_trial().
       */
      template <typename FEValuesType, typename FEExtractorType>
      auto
      shape_function_values(const FEValuesType &   fe_values,
                            const FEExtractorType &extractor,
                            const std::string &    field_symbol) const
        -> std::vector<typename internal::SymbolCreator<
          dim,
          sd_type,
          typename std::decay<decltype(fe_values[extractor])>::type>::
                         value_type>;


      /**
       * Return a vector of shape function values associated with the
       * test function.
       */
      template <typename FEValuesType, typename FEExtractorType>
      auto
      shape_function_values_test(const FEValuesType &   fe_values,
                                 const FEExtractorType &extractor,
                                 const std::string &    field_symbol) const
        -> std::vector<typename internal::SymbolCreator<
          dim,
          sd_type,
          typename std::decay<decltype(fe_values[extractor])>::type>::
                         value_type>;


      /**
       * Return a vector of shape function values associated with the
       * trial solution.
       */
      template <typename FEValuesType, typename FEExtractorType>
      auto
      shape_function_values_trial(const FEValuesType &   fe_values,
                                  const FEExtractorType &extractor,
                                  const std::string &    field_symbol) const
        -> std::vector<typename internal::SymbolCreator<
          dim,
          sd_type,
          typename std::decay<decltype(fe_values[extractor])>::type>::
                         value_type>;

      /**
       * Solution value expressed in terms of the values
       * of all degrees-of-freedom and the shape function values.
       *
       * These values are expressed a product of the degree-of-freedom
       * values and the shape function values derived from the trial
       * solution space.
       */
      template <typename FEValuesType, typename FEExtractorType>
      auto
      solution_value(const FEValuesType &   fe_values,
                     const FEExtractorType &extractor,
                     const std::string &    field_symbol) const ->
        typename internal::SymbolCreator<
          dim,
          sd_type,
          typename std::decay<decltype(fe_values[extractor])>::type>::
          value_type;

      //@}

      /**
       * @name Independent variables: Gradients
       */
      //@{

      /**
       * Return a shape function gradients associated with a specific
       * degree-of-freedom.
       *
       * @note When using this function it is typically assumed that
       * the space from which the test function and trial solution
       * derive are the same. If it is necessary to make a distinction
       * between these, then this is provided by
       * shape_functions_value_test() and shape_functions_value_trial().
       */
      template <typename FEValuesType, typename FEExtractorType>
      auto
      shape_function_gradient(const FEValuesType &   fe_values,
                              const FEExtractorType &extractor,
                              const std::string &    field_symbol,
                              const unsigned int &   function_no) const ->
        typename internal::SymbolCreator<
          dim,
          sd_type,
          typename std::decay<decltype(fe_values[extractor])>::type>::
          gradient_type;

      /**
       * Return the shape function values associated with a specific
       * component of the test function.
       */
      template <typename FEValuesType, typename FEExtractorType>
      auto
      shape_function_gradient_test(const FEValuesType &   fe_values,
                                   const FEExtractorType &extractor,
                                   const std::string &    field_symbol,
                                   const unsigned int &   function_no) const ->
        typename internal::SymbolCreator<
          dim,
          sd_type,
          typename std::decay<decltype(fe_values[extractor])>::type>::
          gradient_type;

      /**
       * Return the shape function values associated with a specific
       * trial solution (corresponding to a degree-of-freedom).
       */
      template <typename FEValuesType, typename FEExtractorType>
      auto
      shape_function_gradient_trial(const FEValuesType &   fe_values,
                                    const FEExtractorType &extractor,
                                    const std::string &    field_symbol,
                                    const unsigned int &   function_no) const ->
        typename internal::SymbolCreator<
          dim,
          sd_type,
          typename std::decay<decltype(fe_values[extractor])>::type>::
          gradient_type;

      /**
       * Return a vector of shape function gradients associated with
       * each degree-of-freedom.
       *
       * @note When using this function it is typically assumed that
       * the space from which the test function and trial solution
       * derive are the same. If it is necessary to make a distinction
       * between these, then this is provided by
       * shape_functions_value_test() and shape_functions_value_trial().
       */
      template <typename FEValuesType, typename FEExtractorType>
      auto
      shape_function_gradients(const FEValuesType &   fe_values,
                               const FEExtractorType &extractor,
                               const std::string &    field_symbol) const
        -> std::vector<typename internal::SymbolCreator<
          dim,
          sd_type,
          typename std::decay<decltype(fe_values[extractor])>::type>::
                         gradient_type>;

      /**
       * Return a vector of shape function values associated with the
       * test function.
       */
      template <typename FEValuesType, typename FEExtractorType>
      auto
      shape_function_gradients_test(const FEValuesType &   fe_values,
                                    const FEExtractorType &extractor,
                                    const std::string &    field_symbol) const
        -> std::vector<typename internal::SymbolCreator<
          dim,
          sd_type,
          typename std::decay<decltype(fe_values[extractor])>::type>::
                         gradient_type>;


      /**
       * Return a vector of shape function values associated with the
       * trial solution.
       */
      template <typename FEValuesType, typename FEExtractorType>
      auto
      shape_function_gradients_trial(const FEValuesType &   fe_values,
                                     const FEExtractorType &extractor,
                                     const std::string &    field_symbol) const
        -> std::vector<typename internal::SymbolCreator<
          dim,
          sd_type,
          typename std::decay<decltype(fe_values[extractor])>::type>::
                         gradient_type>;

      /**
       * Solution value expressed in terms of the values
       * of all degrees-of-freedom and the shape function values.
       */
      template <typename FEValuesType, typename FEExtractorType>
      auto
      solution_gradient(const FEValuesType &   fe_values,
                        const FEExtractorType &extractor,
                        const std::string &    field_symbol) const ->
        typename internal::SymbolCreator<
          dim,
          sd_type,
          typename std::decay<decltype(fe_values[extractor])>::type>::
          gradient_type;

      //@}

      /**
       * @name Independent variables: Curl
       */
      //@{

      // TODO: Curl

      //@}

      /**
       * @name Independent variables: Hessian
       */
      //@{

      // TODO: Hessian

      //@}

      /**
       * @name Independent variables: Divergence
       */
      //@{

      // TODO: Divergence

      //@}

      /**
       * A data structure that stores some settings used
       * by this class.
       */
      const AdditionalData additional_data;

      /**
       * A flag to indicate whether or not the tangent matrix
       * is assumed to be symmetric or not.
       */
      const bool symmetric_system;

      /**
       * @name Symbolic expression optimization and substitution
       */
      //@{

      /**
       * Execute optimisation of all dependent symbolic functions
       */
      void
      optimize(const SD::types::substitution_map &sub_vals_optim);

      /**
       * Collect and substitute all of the substitution values, including
       * - Degree-of-freedom values
       * - Shape function values and those of their gradients etc.
       * - Constitutive parameters and internal variables
       */
      void
      substitute(const SD::types::substitution_map &sub_vals);

      //@}

      /**
       * @name Dependent variables
       */
      //@{

      /**
       * Extract the residual vector from the substituted
       * (and evaluated) optimizer.
       */
      void
      compute_residual(Vector<scalar_type> &residual) const;

      /**
       * Extract the tangent matrix from the substituted
       * (and evaluated) optimizer.
       */
      void
      compute_linearization(FullMatrix<scalar_type> &linearization) const;

      //@}


    protected:
      /**
       * The string to be used as a divider between
       * strings used to compose a symbol.
       */
      std::string
      divider() const;

      /**
       * Invalidate a symbol by giving it a numeric
       * value
       */
      void
      invalidate_symbol(sd_type &symbol) const;

      /**
       * Test if a symbol has been invalidated or not
       */
      bool
      is_valid_symbol(const SymEngine::Basic &symb) const;

      /**
       * Test if a symbol has been invalidated or not
       */
      bool
      is_valid_symbol(const sd_type &symb) const;

      /**
       * @name Optimiser
       */
      //@{

      /**
       * An collection of optimizers that will speed up evaluation of the
       * residual and its linearization.
       */
      std::vector<std::unique_ptr<BatchOptimizer<scalar_type>>> optimizers;

      /**
       * A map indicating which optimizer has been allocated ownership
       * of which dependent symbolic expression. This decomposition is
       * dictated by the AdditionalData::n_batches
       */
      std::map<const sd_type *const, BatchOptimizer<scalar_type> *>
        map_dependent_variable_to_optimizer;

      BatchOptimizer<scalar_type> &
      get_optimizer(const sd_type &dependent_variable) const;

      void
      distribute_dependent_symbols_to_optimizers();

      //@}

      /**
       * @name Dependent variables
       */
      //@{

      /**
       * Residual contribution derived from a single quadrature point
       */
      std::vector<sd_type> residual;

      /**
       * Linearisation of the residual contribution derived from a
       * single quadrature point
       */
      std::vector<std::vector<sd_type>> linearization;

      /**
       * Resize the vector used to store the dependent variables.
       * This includes:
       * - the residual contribution
       * - the tangent contribution
       */
      void
      initialize_dependent_variables(const unsigned int n_independent_variables,
                                     const unsigned int n_dependent_variables);

      /**
       * Initialize the optimizer based on the settings set in the
       * @p additional_data object.
       */
      void
      initialize_optimizer();

      /**
       * Returns the index past the last element stored in the
       * given @p row of the linearization matrix.
       */
      unsigned int
      linearization_row_end(const unsigned int &row) const;

      //@}
    }; // class CellLevelBase


    /**
     * A data structure that defines the labels to be used
     * to contruct symbolic variables identifiers.
     *
     * @note It is critical to ensure that the labels are
     * unique. If not then there is the possibility that one
     * can generate conflicting symbolic expressions that
     * will not be detected during their use.
     */
    template <int dim, typename NumberType, typename ExpressionType>
    struct CellLevelBase<dim, NumberType, ExpressionType>::FESymbolicNames
    {
      /**
       * Default constructor
       */
      FESymbolicNames(const std::string dof_value          = "U",
                      const std::string test_function      = "d",
                      const std::string trial_solution     = "D",
                      const std::string shape_function     = "Nx",
                      const std::string JxW                = "JxW",
                      const std::string gradient           = "Grad",
                      const std::string symmetric_gradient = "symm_Grad",
                      const std::string divergence         = "Div",
                      const std::string curl               = "Curl",
                      const std::string hessian            = "Hessian",
                      const std::string third_derivative   = "3rd_Derivative");

      /**
       * Symbol for a degree-of-freedom value
       */
      const std::string dof_value;

      /**
       * Symbol for the test function
       */
      const std::string test_function;

      /**
       * Symbol for the trial solution
       */
      const std::string trial_solution;

      /**
       * Symbol for a shape function
       */
      const std::string shape_function;

      /**
       * Symbol for the integration constant
       */
      const std::string JxW;

      /**
       * Symbol for the gradient operator
       */
      const std::string gradient;

      /**
       * Symbol for the symmetric gradient operator
       */
      const std::string symmetric_gradient;

      /**
       * Symbol for the divergence operator
       */
      const std::string divergence;

      /**
       * Symbol for the curl operator
       */
      const std::string curl;

      /**
       * Symbol for the hessian
       */
      const std::string hessian;

      /**
       * Symbol for third derivative
       */
      const std::string third_derivative;
    }; // struct CellLevelBase::FESymbolicNames



    /**
     * A struct that defines how the load of doing symbolic
     * expression optimization is to be distributed between
     * the stipulated number of batches.
     */
    enum class LoadBalancing
    {
      /**
       * Perform no load balancing. This is the same as using
       * a single batch to do the work. This option may be highly
       * memory intensive, particularly when using the LLVM JIT
       * compiler.
       */
      none,
      /**
       * Ensure optimal performance by keeping each
       * row of the linearization on a single optimizer
       * instance. This is recommended in 2d, as the symbolic
       * expressions for the linearization are not necessarily
       * lengthy.
       */
      max_performance,
      /**
       * To counter possible memory issues when the symbolic
       * expressions are very long, this method will distribute
       * the elements of the linearization equally amongst all
       * optimizer instances. The residual will, however, still
       * be evaluated by a single optimizer. This is the recommended
       * method for 3d when a finite element has a large number of
       * degrees of freedom (e.g. for a vector-valued problems with
       * higher order elements)
       */
      equal_elements
    };


    /**
     * A data structure that defines the settings to be used
     * for the underlying optimizer.
     */
    template <int dim, typename NumberType, typename ExpressionType>
    struct CellLevelBase<dim, NumberType, ExpressionType>::AdditionalData
    {
      /**
       * Default constructor
       */
      AdditionalData(
        const unsigned int &     n_batches     = 1,
        const LoadBalancing &    LoadBalancing = LoadBalancing::equal_elements,
        const OptimizerType &    optim_method  = OptimizerType::llvm,
        const OptimizationFlags &optim_flags = OptimizationFlags::optimize_all,
        const FESymbolicNames &  symbols     = FESymbolicNames());

      /**
       * The number of batches into which sizable optimization problem will
       * be decomposed. This effectively defines the number of optimizers
       * that will be employed to evaluate a subset of the residual and
       * linearization expressions. In the case that the LLVM optimizer is
       * used, this will decrease the memory consumption during optimization
       * (i.e. conversion from symbolic notation to LLVM bytecode) at the
       * expense of overall performance (due to the few optimizations being
       * performed across the entire range of symbolic expressions).
       *
       * @note When using LLVM optimization, it it may be necessary to tune
       * batch size for each finite element such that the system memory is
       * not exceeded.
       * It is recommended that the number of batches be chosen based on the
       * following criterion, all of which contribute to the number and
       * complexity of the resultant expressions:
       * - Whether the problem is in 1d, 2d or 3d
       * - The number of components in the system
       * - The nature of the system being linearized (does it utilize shape
       * function values, gradients etc.)
       * - The polynomial order of the finite element discretisation
       * - The number of optimizations (note optimizations, not substitutions)
       * simultaneously being performed in parallel .
       *
       * @note  There is no performance / memory benefit to be had when lambda
       * optimization is to be performed, so in this case the a batch size
       * of 1 is recommended to take the best advantage of common
       * subexpression elimination and other symbolic manipulations.
       */
      unsigned int n_batches;

      /**
       * Load balancing method to be applied when distributing work to
       * the optimizers.
       */
      enum LoadBalancing load_balancing;

      /**
       * Optimization method
       */
      enum OptimizerType optim_method;

      /**
       * Symbolic manipulation flags
       */
      enum OptimizationFlags optim_flags;

      /**
       * The textual names given to several fundamental finite-element level
       * operations and values that are to be considered. The option to change
       * these exists in the event of a name clash, but it is recommended that
       * the defaults be used and that a different set of symbols be used
       * by the user when constructing their own symbolic expressions.
       */
      FESymbolicNames symbols;
    }; // struct CellLevelBase::AdditionalData



    /**
     * @brief An assembly helper for the case where the cell contributions
     * to both the residual vector and the tangent stiffness matrix are
     * computed using symbolic differentiation.
     *
     * A symbolic expression for the cell-level potential energy contribution
     * arising from a single quadrature point located on/in a finite element
     * is provided. Both the residual and tangent are then computed
     * symbolically.
     *
     * @author Jean-Paul Pelteret, 2017
     */
    template <int dim,
              typename NumberType,
              typename ExpressionType = SD::Expression>
    class EnergyFunctional
      : public CellLevelBase<dim, NumberType, ExpressionType>
    {
    public:
      using
        typename CellLevelBase<dim, NumberType, ExpressionType>::AdditionalData;
      using typename CellLevelBase<dim, NumberType, ExpressionType>::sd_type;

      /**
       * Constructor
       */
      EnergyFunctional(
        const unsigned int &  n_independent_variables,
        const AdditionalData &additional_data = AdditionalData());

      /**
       * @name Dependent variables
       */
      //@{

      /**
       * Declare the set of dependent variables derived from
       * the quadrature point energy contribution. This approach
       * can be applied if the problem is variational, and is
       * the same as deriving the residual and tangent contributions
       * from a total potential energy.
       *
       * @note In this function, a conversion between the
       * local description of kinematic variables and those
       * defined strictly in terms of values associated with
       * a finite element is performed.
       *
       * @param qp_energy Symbolic expression for the total
       * energy at the quadrature point
       * @param sub_vals_fe Substitution map that converts
       * local field variables to one with a cell-level
       * equivalent (expressed in terms of cell
       * degree-of-freedom values). This ensures that the
       * entire expression is dependent on local finite
       * element cell geometry and global solution.
       */
      void
      register_energy_functional(
        const sd_type &                    qp_energy,
        const SD::types::substitution_map &sub_vals_local_to_fe =
          SD::types::substitution_map());

      //@}
    }; // class EnergyFunctional


    /**
     * @brief An assembly helper for the case where the cell contribution
     * to the tangent stiffness matrix from an individual quadrature point
     * is computed using symbolic differentiation.
     *
     * A symbolic expression for the residual contribution arising from a single
     * quadrature point located on/in a finite element is provided. The tangent
     * is then computed symbolically.
     *
     * @author Jean-Paul Pelteret, 2017
     */
    template <int dim,
              typename NumberType,
              typename ExpressionType = SD::Expression>
    class ResidualLinearization
      : public CellLevelBase<dim, NumberType, ExpressionType>
    {
    public:
      using
        typename CellLevelBase<dim, NumberType, ExpressionType>::AdditionalData;
      using typename CellLevelBase<dim, NumberType, ExpressionType>::sd_type;

      /**
       * Constructor
       */
      ResidualLinearization(
        const unsigned int    n_independent_variables,
        const unsigned int    n_dependent_variables,
        const bool            symmetric_system,
        const AdditionalData &additional_data = AdditionalData());

      /**
       * @name Dependent variables
       */
      //@{

      /**
       * Declare the set of dependent variables derived from
       * the quadrature point residual contribution. This is
       * the most general approach .
       *
       * @note In this function, a conversion between the
       * local description of kinematic variables and those
       * defined strictly in terms of values associated with
       * a finite element is performed.
       *
       * @warning It is currently assumed that the tangent
       * stiffness matrix is symmetric. This holds only
       * for a limited category of problems (e.g. variational
       * problems).
       * @param qp_residual Symbolic expression for the total
       * residual contribution at the quadrature point
       * @param sub_vals_fe Substitution map that converts
       * local field variables to one with a cell-level
       * equivalent (expressed in terms of cell
       * degree-of-freedom values). This ensures that the
       * entire expression is dependent on local finite
       * element cell geometry and global solution.
       */
      void
      register_residual_vector(
        const std::vector<sd_type> &       qp_residual,
        const SD::types::substitution_map &sub_vals_local_to_fe =
          SD::types::substitution_map());

      //@}
    }; // class ResidualLinearization


    /* -------------------- inline and template functions ------------------ */


#  ifndef DOXYGEN


    /* ------------------------- CellLevelBase --------------------- */


    template <int dim, typename NumberType, typename ExpressionType>
    template <typename Stream>
    void
    CellLevelBase<dim, NumberType, ExpressionType>::print(
      Stream &   stream,
      const bool print_residual,
      const bool print_tangent,
      const bool print_optimizer) const
    {
      if (print_residual)
        {
          stream << "Residual contribution" << std::endl;
          for (unsigned int I = 0; I < n_dependent_variables(); ++I)
            {
              Assert(I < residual.size(), ExcIndexRange(I, 0, residual.size()));
              stream << "R[" << I << "]: " << residual[I] << std::endl;
            }
        }
      if (print_tangent)
        {
          stream << "Tangent matrix contribution" << std::endl;
          for (unsigned int I = 0; I < n_dependent_variables(); ++I)
            {
              Assert(I < linearization.size(),
                     ExcIndexRange(I, 0, linearization.size()));
              const unsigned int J_end = linearization_row_end(I);
              for (unsigned int J = 0; J < J_end; ++J)
                {
                  Assert(J < linearization[I].size(),
                         ExcIndexRange(J, 0, linearization[I].size()));
                  stream << "K[" << I << "][" << J
                         << "]: " << linearization[I][J] << std::endl;
                }
            }
        }
      if (print_optimizer)
        for (auto &optimizer : optimizers)
          optimizer->print(stream);
    }

    template <int dim, typename NumberType, typename ExpressionType>
    template <typename FEValuesType, typename FEExtractorType>
    auto
    CellLevelBase<dim, NumberType, ExpressionType>::shape_function_value(
      const FEValuesType &   fe_values,
      const FEExtractorType &extractor,
      const std::string &    field_symbol,
      const unsigned int &   function_no) const ->
      typename internal::SymbolCreator<
        dim,
        sd_type,
        typename std::decay<decltype(fe_values[extractor])>::type>::value_type
    {
      typedef typename std::decay<decltype(fe_values[extractor])>::type
        FEValuesViewType;
#    ifdef DEBUG
      const unsigned int max_f_no =
        std::max(n_independent_variables(), n_dependent_variables());
      Assert(function_no < max_f_no, ExcIndexRange(function_no, 0, max_f_no));
#    endif

      return internal::SymbolCreator<dim, sd_type, FEValuesViewType>::value(
        additional_data.symbols.shape_function + divider() + field_symbol +
        divider() + dealii::Utilities::int_to_string(function_no));
    }


    template <int dim, typename NumberType, typename ExpressionType>
    template <typename FEValuesType, typename FEExtractorType>
    auto
    CellLevelBase<dim, NumberType, ExpressionType>::shape_function_value_test(
      const FEValuesType &   fe_values,
      const FEExtractorType &extractor,
      const std::string &    field_symbol,
      const unsigned int &   function_no) const ->
      typename internal::SymbolCreator<
        dim,
        sd_type,
        typename std::decay<decltype(fe_values[extractor])>::type>::value_type
    {
      Assert(function_no < n_dependent_variables(),
             ExcIndexRange(function_no, 0, n_dependent_variables()));

      // Piggy-back onto the generic function defined above
      return shape_function_value(fe_values,
                                  extractor,
                                  additional_data.symbols.test_function +
                                    field_symbol,
                                  function_no);
    }


    template <int dim, typename NumberType, typename ExpressionType>
    template <typename FEValuesType, typename FEExtractorType>
    auto
    CellLevelBase<dim, NumberType, ExpressionType>::shape_function_value_trial(
      const FEValuesType &   fe_values,
      const FEExtractorType &extractor,
      const std::string &    field_symbol,
      const unsigned int &   function_no) const ->
      typename internal::SymbolCreator<
        dim,
        sd_type,
        typename std::decay<decltype(fe_values[extractor])>::type>::value_type
    {
      Assert(function_no < n_independent_variables(),
             ExcIndexRange(function_no, 0, n_independent_variables()));

      // Piggy-back onto the generic function defined above.
      // To keep things simple and  (i.e. ensure consistency when using
      // solution_value()), we simply assume the default expression for
      // this function.
      return shape_function_value(fe_values,
                                  extractor,
                                  field_symbol,
                                  function_no);
    }


    template <int dim, typename NumberType, typename ExpressionType>
    template <typename FEValuesType, typename FEExtractorType>
    auto
    CellLevelBase<dim, NumberType, ExpressionType>::shape_function_values(
      const FEValuesType &   fe_values,
      const FEExtractorType &extractor,
      const std::string &    field_symbol) const
      -> std::vector<typename internal::SymbolCreator<
        dim,
        sd_type,
        typename std::decay<decltype(fe_values[extractor])>::type>::value_type>
    {
      typedef typename std::decay<decltype(fe_values[extractor])>::type
        FEValuesViewType;
      typedef typename internal::SymbolCreator<dim, sd_type, FEValuesViewType>::
        value_type ValueType;

      // Test function and trial solution must come from same space
      Assert(this->n_dependent_variables() == this->n_independent_variables(),
             ExcInternalError());

      const unsigned int n_independent_variables =
        this->n_independent_variables();
      std::vector<ValueType> values(n_independent_variables);
      for (unsigned int K = 0; K < n_independent_variables; ++K)
        {
          values[K] =
            shape_function_value(fe_values, extractor, field_symbol, K);
        }
      return values;
    }


    template <int dim, typename NumberType, typename ExpressionType>
    template <typename FEValuesType, typename FEExtractorType>
    auto
    CellLevelBase<dim, NumberType, ExpressionType>::shape_function_values_test(
      const FEValuesType &   fe_values,
      const FEExtractorType &extractor,
      const std::string &    field_symbol) const
      -> std::vector<typename internal::SymbolCreator<
        dim,
        sd_type,
        typename std::decay<decltype(fe_values[extractor])>::type>::value_type>
    {
      typedef typename std::decay<decltype(fe_values[extractor])>::type
        FEValuesViewType;
      typedef typename internal::SymbolCreator<dim, sd_type, FEValuesViewType>::
        value_type ValueType;

      const unsigned int n_dependent_variables = this->n_dependent_variables();
      std::vector<ValueType> values(n_dependent_variables);
      for (unsigned int K = 0; K < n_dependent_variables; ++K)
        {
          values[K] =
            shape_function_value_test(fe_values, extractor, field_symbol, K);
        }
      return values;
    }


    template <int dim, typename NumberType, typename ExpressionType>
    template <typename FEValuesType, typename FEExtractorType>
    auto
    CellLevelBase<dim, NumberType, ExpressionType>::shape_function_values_trial(
      const FEValuesType &   fe_values,
      const FEExtractorType &extractor,
      const std::string &    field_symbol) const
      -> std::vector<typename internal::SymbolCreator<
        dim,
        sd_type,
        typename std::decay<decltype(fe_values[extractor])>::type>::value_type>
    {
      typedef typename std::decay<decltype(fe_values[extractor])>::type
        FEValuesViewType;
      typedef typename internal::SymbolCreator<dim, sd_type, FEValuesViewType>::
        value_type ValueType;

      const unsigned int n_independent_variables =
        this->n_independent_variables();
      std::vector<ValueType> values(n_independent_variables);
      for (unsigned int K = 0; K < n_independent_variables; ++K)
        {
          values[K] =
            shape_function_value_trial(fe_values, extractor, field_symbol, K);
        }
      return values;
    }


    template <int dim, typename NumberType, typename ExpressionType>
    template <typename FEValuesType, typename FEExtractorType>
    auto
    CellLevelBase<dim, NumberType, ExpressionType>::solution_value(
      const FEValuesType &   fe_values,
      const FEExtractorType &extractor,
      const std::string &    field_symbol) const ->
      typename internal::SymbolCreator<
        dim,
        sd_type,
        typename std::decay<decltype(fe_values[extractor])>::type>::value_type
    {
      typedef typename std::decay<decltype(fe_values[extractor])>::type
        FEValuesViewType;
      typedef typename internal::SymbolCreator<dim, sd_type, FEValuesViewType>::
        value_type ValueType;

      const unsigned int n_independent_variables =
        this->n_independent_variables();
      const std::vector<sd_type>   dof_values = this->dof_values();
      const std::vector<ValueType> sf_values =
        this->shape_function_values_trial(fe_values, extractor, field_symbol);
      Assert(dof_values.size() == n_independent_variables, ExcInternalError());
      Assert(sf_values.size() == n_independent_variables, ExcInternalError());

      ValueType value(0.0);
      for (unsigned int K = 0; K < n_independent_variables; ++K)
        {
          value += dof_values[K] * sf_values[K];
        }
      return value;
    }


    template <int dim, typename NumberType, typename ExpressionType>
    template <typename FEValuesType, typename FEExtractorType>
    auto
    CellLevelBase<dim, NumberType, ExpressionType>::shape_function_gradient(
      const FEValuesType &   fe_values,
      const FEExtractorType &extractor,
      const std::string &    field_symbol,
      const unsigned int &   function_no) const ->
      typename internal::SymbolCreator<
        dim,
        sd_type,
        typename std::decay<decltype(fe_values[extractor])>::type>::
        gradient_type
    {
      typedef typename std::decay<decltype(fe_values[extractor])>::type
        FEValuesViewType;
      typedef typename internal::SymbolCreator<dim, sd_type, FEValuesViewType>::
        gradient_type GradientType;

      GradientType gradient =
        internal::SymbolCreator<dim, sd_type, FEValuesViewType>::gradient(
          additional_data.symbols.gradient + divider() +
          additional_data.symbols.shape_function + divider() + field_symbol +
          divider() + dealii::Utilities::int_to_string(function_no));

      // As an optimisation, we nullify unnecessary (i.e. guaranteed zero)
      // entries in the shape function gradient tensors.
      // This invalidates these symbols and will be picked up
      // later by the is_valid_symbol() function.
      // Luckily, since we assume a tensor product of the shape functions
      // the gradient will have the same non-zero's at all quadrature points
      // so they can be computed a priori for each shape function index.
      const unsigned int q_point = 0;
      const auto         qp_Grad_Nx =
        fe_values[extractor].gradient(function_no, q_point);
      Assert(qp_Grad_Nx.n_independent_components ==
               gradient.n_independent_components,
             ExcDimensionMismatch(qp_Grad_Nx.n_independent_components,
                                  gradient.n_independent_components));
      for (unsigned int idx = 0; idx < qp_Grad_Nx.n_independent_components;
           ++idx)
        {
          const TableIndices<GradientType::rank> indices =
            qp_Grad_Nx.unrolled_to_component_indices(idx);
          if (qp_Grad_Nx[indices] == 0.0)
            {
              invalidate_symbol(gradient[indices]);
              Assert(is_valid_symbol(gradient[indices]) == false,
                     ExcInternalError());
            }
        }

      return gradient;
    }


    template <int dim, typename NumberType, typename ExpressionType>
    template <typename FEValuesType, typename FEExtractorType>
    auto
    CellLevelBase<dim, NumberType, ExpressionType>::
      shape_function_gradient_test(const FEValuesType &   fe_values,
                                   const FEExtractorType &extractor,
                                   const std::string &    field_symbol,
                                   const unsigned int &   function_no) const ->
      typename internal::SymbolCreator<
        dim,
        sd_type,
        typename std::decay<decltype(fe_values[extractor])>::type>::
        gradient_type
    {
      Assert(function_no < n_dependent_variables(),
             ExcIndexRange(function_no, 0, n_dependent_variables()));

      // Piggy-back onto the generic function defined above
      return shape_function_gradient(fe_values,
                                     extractor,
                                     additional_data.symbols.test_function +
                                       field_symbol,
                                     function_no);
    }


    template <int dim, typename NumberType, typename ExpressionType>
    template <typename FEValuesType, typename FEExtractorType>
    auto
    CellLevelBase<dim, NumberType, ExpressionType>::
      shape_function_gradient_trial(const FEValuesType &   fe_values,
                                    const FEExtractorType &extractor,
                                    const std::string &    field_symbol,
                                    const unsigned int &   function_no) const ->
      typename internal::SymbolCreator<
        dim,
        sd_type,
        typename std::decay<decltype(fe_values[extractor])>::type>::
        gradient_type
    {
      Assert(function_no < n_independent_variables(),
             ExcIndexRange(function_no, 0, n_independent_variables()));

      // Piggy-back onto the generic function defined above.
      // To keep things simple and  (i.e. ensure consistency when using
      // solution_gradient()), we simply assume the default expression for
      // this function.
      return shape_function_gradient(fe_values,
                                     extractor,
                                     field_symbol,
                                     function_no);
    }


    template <int dim, typename NumberType, typename ExpressionType>
    template <typename FEValuesType, typename FEExtractorType>
    auto
    CellLevelBase<dim, NumberType, ExpressionType>::shape_function_gradients(
      const FEValuesType &   fe_values,
      const FEExtractorType &extractor,
      const std::string &    field_symbol) const
      -> std::vector<typename internal::SymbolCreator<
        dim,
        sd_type,
        typename std::decay<decltype(fe_values[extractor])>::type>::
                       gradient_type>
    {
      typedef typename std::decay<decltype(fe_values[extractor])>::type
        FEValuesViewType;
      typedef typename internal::SymbolCreator<dim, sd_type, FEValuesViewType>::
        gradient_type GradientType;

      // Test function and trial solution must come from same space
      Assert(this->n_dependent_variables() == this->n_independent_variables(),
             ExcInternalError());

      const unsigned int n_independent_variables =
        this->n_independent_variables();
      std::vector<GradientType> gradients(n_independent_variables);
      for (unsigned int K = 0; K < n_independent_variables; ++K)
        {
          gradients[K] =
            shape_function_gradient(fe_values, extractor, field_symbol, K);
        }
      return gradients;
    }


    template <int dim, typename NumberType, typename ExpressionType>
    template <typename FEValuesType, typename FEExtractorType>
    auto
    CellLevelBase<dim, NumberType, ExpressionType>::
      shape_function_gradients_test(const FEValuesType &   fe_values,
                                    const FEExtractorType &extractor,
                                    const std::string &    field_symbol) const
      -> std::vector<typename internal::SymbolCreator<
        dim,
        sd_type,
        typename std::decay<decltype(fe_values[extractor])>::type>::
                       gradient_type>
    {
      typedef typename std::decay<decltype(fe_values[extractor])>::type
        FEValuesViewType;
      typedef typename internal::SymbolCreator<dim, sd_type, FEValuesViewType>::
        gradient_type GradientType;

      const unsigned int n_dependent_variables = this->n_dependent_variables();
      std::vector<GradientType> gradients(n_dependent_variables);
      for (unsigned int K = 0; K < n_dependent_variables; ++K)
        {
          gradients[K] =
            shape_function_gradient_test(fe_values, extractor, field_symbol, K);
        }
      return gradients;
    }


    template <int dim, typename NumberType, typename ExpressionType>
    template <typename FEValuesType, typename FEExtractorType>
    auto
    CellLevelBase<dim, NumberType, ExpressionType>::
      shape_function_gradients_trial(const FEValuesType &   fe_values,
                                     const FEExtractorType &extractor,
                                     const std::string &    field_symbol) const
      -> std::vector<typename internal::SymbolCreator<
        dim,
        sd_type,
        typename std::decay<decltype(fe_values[extractor])>::type>::
                       gradient_type>
    {
      typedef typename std::decay<decltype(fe_values[extractor])>::type
        FEValuesViewType;
      typedef typename internal::SymbolCreator<dim, sd_type, FEValuesViewType>::
        gradient_type GradientType;

      const unsigned int n_independent_variables =
        this->n_independent_variables();
      std::vector<GradientType> gradients(n_independent_variables);
      for (unsigned int K = 0; K < n_independent_variables; ++K)
        {
          gradients[K] = shape_function_gradient_trial(fe_values,
                                                       extractor,
                                                       field_symbol,
                                                       K);
        }
      return gradients;
    }


    template <int dim, typename NumberType, typename ExpressionType>
    template <typename FEValuesType, typename FEExtractorType>
    auto
    CellLevelBase<dim, NumberType, ExpressionType>::solution_gradient(
      const FEValuesType &   fe_values,
      const FEExtractorType &extractor,
      const std::string &    field_symbol) const ->
      typename internal::SymbolCreator<
        dim,
        sd_type,
        typename std::decay<decltype(fe_values[extractor])>::type>::
        gradient_type
    {
      typedef typename std::decay<decltype(fe_values[extractor])>::type
        FEValuesViewType;
      typedef typename internal::SymbolCreator<dim, sd_type, FEValuesViewType>::
        gradient_type GradientType;

      const unsigned int n_independent_variables =
        this->n_independent_variables();
      const std::vector<sd_type>      dof_values = this->dof_values();
      const std::vector<GradientType> sf_grads =
        this->shape_function_gradients(fe_values, extractor, field_symbol);
      Assert(dof_values.size() == n_independent_variables, ExcInternalError());
      Assert(sf_grads.size() == n_independent_variables, ExcInternalError());

      GradientType gradient;
      for (unsigned int K = 0; K < n_independent_variables; ++K)
        {
          gradient += dof_values[K] * sf_grads[K];
        }
      return gradient;
    }

#  endif


  } // namespace SD
} // namespace Differentiation

DEAL_II_NAMESPACE_CLOSE

#endif // defined(DEAL_II_WITH_SYMENGINE)

#endif
