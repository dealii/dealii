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

#  include <deal.II/differentiation/sd/symengine_helpers.h>
#  include <deal.II/differentiation/sd/symengine_number_types.h>

#  include <algorithm>
#  include <iterator>
#  include <queue>
#  include <vector>


DEAL_II_NAMESPACE_OPEN


namespace Differentiation
{
  namespace SD
  {
    namespace SE = ::SymEngine;

    // -------------------------- CellLevelBase ----------------------


    template <int dim, typename NumberType, typename ExpressionType>
    CellLevelBase<dim, NumberType, ExpressionType>::FESymbolicNames::
      FESymbolicNames(const std::string dof_value,
                      const std::string test_function,
                      const std::string trial_solution,
                      const std::string shape_function,
                      const std::string JxW,
                      const std::string gradient,
                      const std::string symmetric_gradient,
                      const std::string divergence,
                      const std::string curl,
                      const std::string hessian,
                      const std::string third_derivative)
      : dof_value(dof_value)
      , test_function(test_function)
      , trial_solution(trial_solution)
      , shape_function(shape_function)
      , JxW(JxW)
      , gradient(gradient)
      , symmetric_gradient(symmetric_gradient)
      , divergence(divergence)
      , curl(curl)
      , hessian(hessian)
      , third_derivative(third_derivative)
    {}



    template <int dim, typename NumberType, typename ExpressionType>
    CellLevelBase<dim, NumberType, ExpressionType>::AdditionalData::
      AdditionalData(const unsigned int &      n_batches,
                     const enum LoadBalancing &load_balancing,
                     const OptimizerType &     optim_method,
                     const OptimizationFlags & optim_flags,
                     const FESymbolicNames &   symbols)
      : n_batches(n_batches)
      , load_balancing(load_balancing)
      , optim_method(optim_method)
      , optim_flags(optim_flags)
      , symbols(symbols)
    {}



    template <int dim, typename NumberType, typename ExpressionType>
    CellLevelBase<dim, NumberType, ExpressionType>::CellLevelBase(
      const unsigned int    n_independent_variables,
      const unsigned int    n_dependent_variables,
      const bool            symmetric_system,
      const AdditionalData &additional_data)
      : additional_data(additional_data)
      , symmetric_system(symmetric_system)
    {
      if (symmetric_system)
        AssertThrow(
          n_dependent_variables == n_independent_variables,
          ExcMessage(
            "The number of dependent and independent variables should be equal "
            "when the system is flagged as being symmetric."));

      initialize_dependent_variables(n_independent_variables,
                                     n_dependent_variables);
      initialize_optimizer();
    }


    template <int dim, typename NumberType, typename ExpressionType>
    CellLevelBase<dim, NumberType, ExpressionType>::CellLevelBase(
      CellLevelBase &&rhs)
      : additional_data(std::move(rhs.additional_data))
      , symmetric_system(std::move(rhs.symmetric_system))
      , optimizers(std::move(rhs.optimizers))
      , residual(std::move(rhs.residual))
      , linearization(std::move(rhs.linearization))
    {}


    template <int dim, typename NumberType, typename ExpressionType>
    std::size_t
    CellLevelBase<dim, NumberType, ExpressionType>::n_independent_variables()
      const
    {
      // Regardless of whether the tangent matrix is marked as symmetric or
      // not, the last row of the stored data will have the right size (since
      // we store, at the very least, the lower triangle + diagonal).
      return linearization.back().size();
    }


    template <int dim, typename NumberType, typename ExpressionType>
    std::size_t
    CellLevelBase<dim, NumberType, ExpressionType>::n_dependent_variables()
      const
    {
      return residual.size();
    }


    template <int dim, typename NumberType, typename ExpressionType>
    std::vector<
      typename CellLevelBase<dim, NumberType, ExpressionType>::sd_type>
    CellLevelBase<dim, NumberType, ExpressionType>::dof_values() const
    {
      const unsigned int n_independent_variables =
        this->n_independent_variables();
      std::vector<sd_type> values(n_independent_variables);
      for (unsigned int K = 0; K < n_independent_variables; ++K)
        {
          values[K] =
            make_symbol(additional_data.symbols.dof_value + divider() +
                        dealii::Utilities::int_to_string(K));
        }
      return values;
    }


    template <int dim, typename NumberType, typename ExpressionType>
    typename CellLevelBase<dim, NumberType, ExpressionType>::sd_type
    CellLevelBase<dim, NumberType, ExpressionType>::JxW() const
    {
      return additional_data.symbols.JxW;
    }


    template <int dim, typename NumberType, typename ExpressionType>
    void
    CellLevelBase<dim, NumberType, ExpressionType>::optimize(
      const SD::types::substitution_map &sub_vals_optim)
    {
      // Register independent symbols
      for (auto &optimizer : optimizers)
        optimizer->register_symbols(sub_vals_optim);

      // Decide on the decomposition of symbolic expressions amongst the
      // collection of optimizers
      distribute_dependent_symbols_to_optimizers();

      // Dependent symbolic functions
      const unsigned int n_dependent_variables = this->n_dependent_variables();
      for (unsigned int I = 0; I < n_dependent_variables; ++I)
        {
          Assert(I < residual.size(), ExcIndexRange(I, 0, residual.size()));
          Assert(I < linearization.size(),
                 ExcIndexRange(I, 0, linearization.size()));

          BatchOptimizer<NumberType> &optimizer_res_I =
            get_optimizer(residual[I]);
          optimizer_res_I.register_function(residual[I]);

          const unsigned int J_end = this->linearization_row_end(I);
          for (unsigned int J = 0; J < J_end; ++J)
            {
              Assert(J < linearization[I].size(),
                     ExcIndexRange(J, 0, linearization[I].size()));
              BatchOptimizer<NumberType> &optimizer_lin_IJ =
                get_optimizer(linearization[I][J]);
              optimizer_lin_IJ.register_function(linearization[I][J]);
            }
        }

      // Perform symbolic optimisation
      //    unsigned int o = 0;
      for (auto &optimizer : optimizers)
        {
          //      std::cout << "Optimising " << o++ << " of " <<
          //      optimizers.size() << std::endl;
          optimizer->optimize();
        }
    }


    template <int dim, typename NumberType, typename ExpressionType>
    void
    CellLevelBase<dim, NumberType, ExpressionType>::substitute(
      const SD::types::substitution_map &sub_vals)
    {
      for (auto &optimizer : optimizers)
        optimizer->substitute(sub_vals);
    }


    template <int dim, typename NumberType, typename ExpressionType>
    void
    CellLevelBase<dim, NumberType, ExpressionType>::compute_residual(
      Vector<NumberType> &res_out) const
    {
      const unsigned int n_dependent_variables = this->n_dependent_variables();

      // We can neglect correctly initializing the entries as
      // we'll be overwriting them immediately.
      if (res_out.size() != n_dependent_variables)
        res_out.reinit(n_dependent_variables, true /*omit_zeroing_entries*/);

      for (unsigned int I = 0; I < n_dependent_variables; ++I)
        {
          Assert(I < residual.size(), ExcIndexRange(I, 0, residual.size()));
          Assert(I < res_out.size(), ExcIndexRange(I, 0, res_out.size()));
          BatchOptimizer<NumberType> &optimizer_res_I =
            get_optimizer(residual[I]);
          res_out[I] = optimizer_res_I.evaluate(residual[I]);
        }
    }


    template <int dim, typename NumberType, typename ExpressionType>
    void
    CellLevelBase<dim, NumberType, ExpressionType>::compute_linearization(
      FullMatrix<NumberType> &lin_out) const
    {
      const unsigned int n_dependent_variables = this->n_dependent_variables();
      const unsigned int n_independent_variables =
        this->n_independent_variables();

      // We can neglect correctly initializing the entries as
      // we'll be overwriting them immediately.
      if (lin_out.m() != n_dependent_variables &&
          lin_out.n() != n_independent_variables)
        lin_out.reinit({n_dependent_variables, n_independent_variables},
                       true /*omit_default_initialization*/);

      for (unsigned int I = 0; I < n_dependent_variables; ++I)
        {
          Assert(I < linearization.size(),
                 ExcIndexRange(I, 0, linearization.size()));
          Assert(I < lin_out.m(), ExcIndexRange(I, 0, lin_out.m()));

          const unsigned int J_end = linearization_row_end(I);
          for (unsigned int J = 0; J < J_end; ++J)
            {
              Assert(J < linearization[I].size(),
                     ExcIndexRange(J, 0, linearization[I].size()));
              Assert(J < lin_out.n(), ExcIndexRange(J, 0, lin_out.n()));

              BatchOptimizer<NumberType> &optimizer_lin_IJ =
                get_optimizer(linearization[I][J]);
              lin_out[I][J] = optimizer_lin_IJ.evaluate(linearization[I][J]);
              if (symmetric_system == true && I != J)
                lin_out[J][I] = lin_out[I][J]; // Symmetry
            }
        }
    }


    template <int dim, typename NumberType, typename ExpressionType>
    std::string
    CellLevelBase<dim, NumberType, ExpressionType>::divider() const
    {
      return "_";
    }


    template <int dim, typename NumberType, typename ExpressionType>
    void
    CellLevelBase<dim, NumberType, ExpressionType>::invalidate_symbol(
      sd_type &symbol) const
    {
      symbol = 0.0;
    }


    template <int dim, typename NumberType, typename ExpressionType>
    bool
    CellLevelBase<dim, NumberType, ExpressionType>::is_valid_symbol(
      const SE::Basic &symb) const
    {
      return SE::is_a<SE::Symbol>(symb);
    }


    template <int dim, typename NumberType, typename ExpressionType>
    bool
    CellLevelBase<dim, NumberType, ExpressionType>::is_valid_symbol(
      const sd_type &symb) const
    {
      return is_valid_symbol(symb.get_value());
    }


    template <int dim, typename NumberType, typename ExpressionType>
    BatchOptimizer<NumberType> &
    CellLevelBase<dim, NumberType, ExpressionType>::get_optimizer(
      const sd_type &dependent_variable) const
    {
      const auto it =
        map_dependent_variable_to_optimizer.find(&dependent_variable);
      Assert(it != map_dependent_variable_to_optimizer.end(),
             ExcMessage("This variable is not allocated to any optimizer."));
      return *(it->second);
    }


    template <int dim, typename NumberType, typename ExpressionType>
    void
    CellLevelBase<dim, NumberType, ExpressionType>::
      distribute_dependent_symbols_to_optimizers()
    {
      map_dependent_variable_to_optimizer.clear();
      const unsigned int &n_batches = additional_data.n_batches;
      Assert(
        optimizers.size() == n_batches,
        ExcMessage(
          "Optimizers are not initialized correctly. This should have been done in the call to initialize_optimizer()."));
      const unsigned int n_dependent_variables = this->n_dependent_variables();

      if (n_batches == 1 ||
          additional_data.load_balancing == LoadBalancing::none)
        {
          // Assignment is trivial
          for (unsigned int I = 0; I < n_dependent_variables; ++I)
            {
              Assert(I < residual.size(), ExcIndexRange(I, 0, residual.size()));
              map_dependent_variable_to_optimizer[&(residual[I])] =
                optimizers[0].get();
            }
          // For the same reason we also wish to group together all of the
          // entries in a row of the linearisation.
          for (unsigned int I = 0; I < n_dependent_variables; ++I)
            {
              Assert(I < linearization.size(),
                     ExcIndexRange(I, 0, linearization.size()));
              const unsigned int J_end = linearization_row_end(I);
              for (unsigned int J = 0; J < J_end; ++J)
                {
                  Assert(J < linearization[I].size(),
                         ExcIndexRange(J, 0, linearization[I].size()));
                  map_dependent_variable_to_optimizer[&(linearization[I][J])] =
                    optimizers[0].get();
                }
            }
        }
      else
        {
          // Get a collection of pointers to the symbolic variables to
          // be batched together.
          std::vector<std::vector<sd_type *>> dependent_variable_bin;
          dependent_variable_bin.reserve(1 + n_dependent_variables);
          // We wish to group together all of the residual entries in a
          // single optimizer to take best advantage of CSE.
          dependent_variable_bin.push_back(
            std::vector<sd_type *>(n_dependent_variables));
          unsigned int n_dependents_extracted = 0;
          for (unsigned int I = 0; I < n_dependent_variables; ++I)
            {
              Assert(I < residual.size(), ExcIndexRange(I, 0, residual.size()));
              dependent_variable_bin.back()[I] = &(residual[I]);
              ++n_dependents_extracted;
            }
          // For the same reason we also wish to group together all of the
          // entries in a row of the linearisation.
          // We add from the last row to the first we'll be combining multiple
          // rows together to balance out the workload. This can be leveraged in
          // the symmetric case to add the smaller entries last to get better
          // balancing.
          for (int I = static_cast<int>(n_dependent_variables) - 1; I >= 0; --I)
            {
              Assert(I < static_cast<int>(linearization.size()),
                     ExcIndexRange(I, 0, linearization.size()));
              const unsigned int J_end = linearization_row_end(I);

              dependent_variable_bin.push_back(std::vector<sd_type *>(J_end));
              for (unsigned int J = 0; J < J_end; ++J)
                {
                  Assert(J < static_cast<unsigned int>(linearization[I].size()),
                         ExcIndexRange(J, 0, linearization[I].size()));
                  dependent_variable_bin.back()[J] = &(linearization[I][J]);
                  ++n_dependents_extracted;
                }
            }

          // Perform a quick sanity check
          const unsigned int n_dependents_total =
            this->n_dependent_variables() // residual.size()
            +
            (symmetric_system == false ?
               (this->n_dependent_variables() *
                this->n_independent_variables()) :
               (this->n_dependent_variables() * this->n_dependent_variables() +
                this->n_dependent_variables()) /
                 2);
          Assert(n_dependents_extracted == n_dependents_total,
                 ExcMessage("Not all dependent variables were extracted."));
          (void)n_dependents_total;

          // Perform the assignment to the worker optimizers
          std::vector<std::vector<sd_type *>> assignment_bin(n_batches);
          if (additional_data.load_balancing == LoadBalancing::max_performance)
            {
              auto f_vector_less_than =
                [](const std::vector<sd_type *> &lhs,
                   const std::vector<sd_type *> &rhs) -> bool {
                return lhs.size() < rhs.size();
              };

              // Put an entire row of dependents (that come from either the
              // residual or an entire row of its linearization) on a single
              // optimizer.
              for (unsigned int e = 0; e < dependent_variable_bin.size(); ++e)
                {
                  // Get the entry in the assignment bin with the least amount
                  // of dependent variables assigned to it
                  auto it_assignment = std::min_element(assignment_bin.begin(),
                                                        assignment_bin.end(),
                                                        f_vector_less_than);
                  std::copy(dependent_variable_bin[e].begin(),
                            dependent_variable_bin[e].end(),
                            std::back_inserter(*it_assignment));
                }
            }
          else if (additional_data.load_balancing ==
                   LoadBalancing::equal_elements)
            {
              Assert(n_batches >= 2, ExcInternalError());
              Assert(n_dependents_total > this->n_dependent_variables(),
                     ExcInternalError());
              const unsigned int n_dependents_linearization =
                n_dependents_total - this->n_dependent_variables();
              (void)n_dependents_linearization;

              // We will be performing a lot of data moving from this object
              // (from front to back), so instead of using a vector which would
              // reallocate memory during a move we simply use a FIFO queue
              std::queue<sd_type *> unassigned_dep_lin;
              for (unsigned int e = 0; e < dependent_variable_bin.size(); ++e)
                {
                  if (e == 0)
                    // The whole residual is to be assigned to a single
                    // optimizer
                    std::copy(dependent_variable_bin[e].begin(),
                              dependent_variable_bin[e].end(),
                              std::back_inserter(assignment_bin[0]));
                  else
                    // Take the contents of the linearization and make a single
                    // long list of assignments to be performed.
                    for (auto it = dependent_variable_bin[e].begin();
                         it != dependent_variable_bin[e].end();
                         ++it)
                      unassigned_dep_lin.push(*it);
                }
              Assert(unassigned_dep_lin.size() == n_dependents_linearization,
                     ExcInternalError());

              // Move the data from the queue of unassigned dependents to a bin
              // To do this in a balanced way, we work out the equal
              // distribution of remaining work after each assignment is done.
              // This is easier than trying to do the allocation a-priori and
              // having to work around a non-integer allocation to the bins
              for (unsigned int b = 1; b < assignment_bin.size(); ++b)
                {
                  const unsigned int n_batches_remaining = n_batches - b;
                  const unsigned int n_dep_lin_remaining =
                    unassigned_dep_lin.size();
                  const unsigned int est_n_deps_lin_per_optimizer =
                    static_cast<unsigned int>(
                      static_cast<float>(n_dep_lin_remaining) /
                      n_batches_remaining);

                  assignment_bin[b].reserve(est_n_deps_lin_per_optimizer);
                  for (unsigned int e = 0; e < est_n_deps_lin_per_optimizer;
                       ++e)
                    {
                      assignment_bin[b].push_back(
                        std::move(unassigned_dep_lin.front()));
                      unassigned_dep_lin.pop();
                    }
                }
              Assert(unassigned_dep_lin.empty(),
                     ExcMessage(
                       "Expected all elements to be removed from queue."));
            }
          else
            {
              AssertThrow(false,
                          ExcMessage(
                            "Unknown load balancing algorithm selected."));
            }

          //        for (unsigned int b=0; b<assignment_bin.size(); ++b)
          //          std::cout << "Bin: " << b << "  size: " <<
          //          assignment_bin[b].size() << std::endl;

          // Fill the pointer map with the computed assignment
          for (unsigned int b = 0; b < assignment_bin.size(); ++b)
            for (unsigned int e = 0; e < assignment_bin[b].size(); ++e)
              map_dependent_variable_to_optimizer[assignment_bin[b][e]] =
                optimizers[b].get();

          Assert(map_dependent_variable_to_optimizer.size() ==
                   n_dependents_total,
                 ExcMessage("Not all dependent variables were assigned."));
        }
    }


    template <int dim, typename NumberType, typename ExpressionType>
    void
    CellLevelBase<dim, NumberType, ExpressionType>::
      initialize_dependent_variables(const unsigned int n_independent_variables,
                                     const unsigned int n_dependent_variables)
    {
      // Setup residual and tangent contribution
      // Note that we store only the lower half
      // and diagonal of the tangent contribution
      residual.resize(n_dependent_variables);
      linearization.resize(n_dependent_variables);
      for (unsigned int K = 0; K < n_dependent_variables; ++K)
        if (symmetric_system == true)
          // Store only the bottom half of the symmetric
          // tensor plus the diagonal
          linearization[K].resize(K + 1);
        else
          // Store the entire row
          linearization[K].resize(n_independent_variables);
    }


    template <int dim, typename NumberType, typename ExpressionType>
    void
    CellLevelBase<dim, NumberType, ExpressionType>::initialize_optimizer()
    {
      optimizers.resize(additional_data.n_batches);

      for (unsigned int o = 0; o < additional_data.n_batches; ++o)
        optimizers[o].reset(
          new BatchOptimizer<NumberType>(additional_data.optim_method,
                                         additional_data.optim_flags));
    }


    template <int dim, typename NumberType, typename ExpressionType>
    unsigned int
    CellLevelBase<dim, NumberType, ExpressionType>::linearization_row_end(
      const unsigned int &row) const
    {
      // In the symmetric case, we need to ensure that J's range
      // includes the diagonal
      return (symmetric_system == true ? row + 1 : n_independent_variables());
    }


    /* --------------------- EnergyFunctional --------------------- */


    template <int dim, typename NumberType, typename ExpressionType>
    EnergyFunctional<dim, NumberType, ExpressionType>::EnergyFunctional(
      const unsigned int &  n_independent_variables,
      const AdditionalData &additional_data)
      : CellLevelBase<dim, NumberType, ExpressionType>(
          n_independent_variables,
          n_independent_variables,
          true /*symmetric_system*/,
          additional_data)
    {}


    template <int dim, typename NumberType, typename ExpressionType>
    void
    EnergyFunctional<dim, NumberType, ExpressionType>::
      register_energy_functional(
        const sd_type &                    qp_energy,
        const SD::types::substitution_map &sub_vals_local_to_fe)
    {
      Assert(this->symmetric_system == true,
             ExcMessage(
               "Variational problems always result in a symmetric system."));

      // Convert the local field variable descriptions into those
      // that are dependent on local finite element cell geometry and
      // global solution.
      sd_type qp_energy_fe(0.0);
      if (!sub_vals_local_to_fe.empty())
        {
          qp_energy_fe = substitute(qp_energy, sub_vals_local_to_fe);
          Assert(numbers::values_are_not_equal(qp_energy_fe, sd_type(0.0)),
                 ExcMessage("Bad substitution."));
          Assert(numbers::values_are_not_equal(qp_energy_fe, qp_energy),
                 ExcMessage("No substitution has taken place!"));
        }
      else
        qp_energy_fe = qp_energy;

      const std::vector<sd_type> dof_values = this->dof_values();

      // Now we perform the differentiation operations and compute
      // both the residual and its linearization
      const unsigned int n_dependent_variables = this->n_dependent_variables();
      for (unsigned int I = 0; I < n_dependent_variables; ++I)
        {
          Assert(I < this->residual.size(),
                 ExcIndexRange(I, 0, this->residual.size()));
          Assert(I < this->linearization.size(),
                 ExcIndexRange(I, 0, this->linearization.size()));

          this->residual[I] = qp_energy_fe.differentiate(dof_values[I]);
          for (unsigned int J = 0; J <= I; ++J) // Symmetry
            {
              Assert(J < this->linearization[I].size(),
                     ExcIndexRange(J, 0, this->linearization[I].size()));
              this->linearization[I][J] =
                this->residual[I].differentiate(dof_values[J]);
            }
        }
    }


    /* -------------------- ResidualLinearization ------------------ */


    template <int dim, typename NumberType, typename ExpressionType>
    ResidualLinearization<dim, NumberType, ExpressionType>::
      ResidualLinearization(const unsigned int    n_independent_variables,
                            const unsigned int    n_dependent_variables,
                            const bool            symmetric_system,
                            const AdditionalData &additional_data)
      : CellLevelBase<dim, NumberType, ExpressionType>(n_independent_variables,
                                                       n_dependent_variables,
                                                       symmetric_system,
                                                       additional_data)
    {}


    template <int dim, typename NumberType, typename ExpressionType>
    void
    ResidualLinearization<dim, NumberType, ExpressionType>::
      register_residual_vector(
        const std::vector<sd_type> &       qp_residual,
        const SD::types::substitution_map &sub_vals_local_to_fe)
    {
      const std::vector<sd_type> dof_values    = this->dof_values();
      const unsigned int n_dependent_variables = this->n_dependent_variables();
      const unsigned int n_independent_variables =
        this->n_independent_variables();
      (void)n_independent_variables;

      Assert(qp_residual.size() == n_dependent_variables,
             ExcDimensionMismatch(qp_residual.size(), n_dependent_variables));
      Assert(dof_values.size() == n_independent_variables,
             ExcDimensionMismatch(dof_values.size(), n_independent_variables));

      // Now we perform the differentiation operations and compute
      // both the residual and its linearization
      for (unsigned int I = 0; I < n_dependent_variables; ++I)
        {
          Assert(I < this->residual.size(),
                 ExcIndexRange(I, 0, this->residual.size()));
          Assert(I < this->linearization.size(),
                 ExcIndexRange(I, 0, this->linearization.size()));

          // Convert the local field variable descriptions into those
          // that are dependent on local finite element cell geometry and
          // global solution.
          if (!sub_vals_local_to_fe.empty())
            {
              this->residual[I] =
                substitute(qp_residual[I], sub_vals_local_to_fe);
              Assert(numbers::values_are_not_equal(this->residual[I],
                                                   sd_type(0.0)),
                     ExcMessage("Bad substitution."));
              //            Assert(SymEngine::is_a<SymEngine::Symbol>(this->residual[I].get_value()),
              //            ExcMessage("Bad substitution."));
              //            Assert(this->residual[I] != qp_residual[I],
              //            ExcMessage("No substitution has taken place!"));
            }
          else
            this->residual[I] = qp_residual[I];

          const unsigned int J_end = this->linearization_row_end(I);
          for (unsigned int J = 0; J < J_end; ++J)
            {
              Assert(J < this->linearization[I].size(),
                     ExcIndexRange(J, 0, this->linearization[I].size()));
              this->linearization[I][J] =
                this->residual[I].differentiate(dof_values[J]);
            }
        }
    }

  } // namespace SD
} // namespace Differentiation

/* --- Explicit instantiations --- */

#  include "symengine_helpers.inst"


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SYMENGINE
