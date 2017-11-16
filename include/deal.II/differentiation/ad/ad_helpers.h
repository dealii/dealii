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

#ifndef dealii_differentiation_ad_ad_helpers_h
#define dealii_differentiation_ad_ad_helpers_h

#include <deal.II/base/config.h>

#if defined(DEAL_II_WITH_ADOLC) || defined(DEAL_II_WITH_TRILINOS)

#include <deal.II/base/numbers.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/symmetric_tensor.h>

#include <deal.II/differentiation/ad/ad_number_traits.h>
#include <deal.II/differentiation/ad/adolc_number_types.h>
#include <deal.II/differentiation/ad/adolc_product_types.h>
#include <deal.II/differentiation/ad/sacado_number_types.h>
#include <deal.II/differentiation/ad/sacado_product_types.h>

#include <deal.II/fe/fe_values_extractors.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#ifdef DEAL_II_WITH_ADOLC

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <adolc/drivers/drivers.h>
#include <adolc/taping.h>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#endif // DEAL_II_WITH_ADOLC

#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <set>

DEAL_II_NAMESPACE_OPEN

namespace Differentiation
{
  namespace AD
  {

    namespace internal
    {
      namespace
      {
        // --- Define a whole bunch of helper structs to assist with getting data
        //     related to FEExtractors ---
        template <int dim, typename ExtractorType>
        struct Extractor;

        template <int dim>
        struct Extractor<dim,FEValuesExtractors::Scalar>
        {
          static const unsigned int rank = 0;
          static const unsigned int n_components = 1;
          static bool symmetric_component (const unsigned int idx)
          {
            return false;
          }

          template<typename NumberType>
          using tensor_type = Tensor<rank,dim,NumberType>; // NumberType;

          // Note: FEValuesViews::Scalar::tensor_type is double, so we can't
          // use it (FEValuesViews) in this context.
          // In fact, sadly, all FEValuesViews objects expect doubles as value
          // types
          template<typename NumberType>
          using value_type = NumberType;

          template<typename NumberType>
          using gradient_type = Tensor<rank+1,dim,NumberType>; // NumberType;

          template<typename NumberType, typename IndexType = unsigned int, int rank_in>
          static IndexType
          local_component(const TableIndices<rank_in> &table_indices,
                          const unsigned int           column_offset)
          {
            return 0;
          }
        };

        template <int dim>
        struct Extractor<dim,FEValuesExtractors::Vector>
        {
          static const unsigned int rank = 1; // Tensor<1,dim>::rank;
          static const unsigned int n_components = dim; // Tensor<1,dim>::n_independent_components
          static bool symmetric_component (const unsigned int idx)
          {
            return false;
          }

          template<typename NumberType>
          using tensor_type = Tensor<rank,dim,NumberType>;

          template<typename NumberType>
          using value_type = tensor_type<NumberType>;

          template<typename NumberType>
          using gradient_type = Tensor<rank+1,dim,NumberType>;

          template<int rank_in>
          static TableIndices<rank>
          table_index_view (const TableIndices<rank_in> &table_indices,
                            const unsigned int           column_offset)
          {
            Assert(0+column_offset < rank_in, ExcInternalError());
            return TableIndices<rank>(table_indices[0+column_offset]);
          }

          template<typename NumberType, typename IndexType = unsigned int, int rank_in>
          static IndexType
          local_component(const TableIndices<rank_in> &table_indices,
                          const unsigned int           column_offset)
          {
            typedef tensor_type<NumberType> TensorType;
            return TensorType::component_to_unrolled_index(table_index_view(table_indices,
                                                           column_offset));
          }
        };

        template <int dim>
        struct Extractor< dim,FEValuesExtractors::Tensor<1> >
        {
          static const unsigned int rank = Tensor<1,dim>::rank;
          static const unsigned int n_components = Tensor<1,dim>::n_independent_components;
          static bool symmetric_component (const unsigned int idx)
          {
            return false;
          }

          template<typename NumberType>
          using tensor_type = Tensor<rank,dim,NumberType>;

          template<typename NumberType>
          using value_type = tensor_type<NumberType>;

          template<typename NumberType>
          using gradient_type = Tensor<rank+1,dim,NumberType>;

          template<int rank_in>
          static TableIndices<rank>
          table_index_view (const TableIndices<rank_in> &table_indices,
                            const unsigned int           column_offset)
          {
            Assert(column_offset < rank_in, ExcInternalError());
            return TableIndices<rank>(table_indices[column_offset]);
          }

          template<typename NumberType, typename IndexType = unsigned int, int rank_in>
          static IndexType
          local_component(const TableIndices<rank_in> &table_indices,
                          const unsigned int           column_offset)
          {
            typedef tensor_type<NumberType> TensorType;
            return TensorType::component_to_unrolled_index(table_index_view(table_indices,
                                                           column_offset));
          }
        };

        template <int dim>
        struct Extractor< dim,FEValuesExtractors::Tensor<2> >
        {
          static const unsigned int rank = Tensor<2,dim>::rank;
          static const unsigned int n_components = Tensor<2,dim>::n_independent_components;
          static bool symmetric_component (const unsigned int idx)
          {
            return false;
          }

          template<typename NumberType>
          using tensor_type = Tensor<rank,dim,NumberType>;

          template<typename NumberType>
          using value_type = tensor_type<NumberType>;

          template<typename NumberType>
          using gradient_type = Tensor<rank+1,dim,NumberType>;

          template<int rank_in>
          static TableIndices<rank>
          table_index_view (const TableIndices<rank_in> &table_indices,
                            const unsigned int           column_offset)
          {
            Assert(column_offset < rank_in, ExcInternalError());
            Assert(1+column_offset < rank_in, ExcInternalError());
            return TableIndices<rank>(table_indices[column_offset],
                                      table_indices[1+column_offset]);
          }

          template<typename NumberType, typename IndexType = unsigned int, int rank_in>
          static IndexType
          local_component(const TableIndices<rank_in> &table_indices,
                          const unsigned int           column_offset)
          {
            typedef tensor_type<NumberType> TensorType;
            return TensorType::component_to_unrolled_index(table_index_view(table_indices,
                                                           column_offset));
          }
        };

        template <int dim>
        struct Extractor< dim,FEValuesExtractors::SymmetricTensor<2> >
        {
          static const unsigned int rank = SymmetricTensor<2,dim>::rank;
          static const unsigned int n_components = SymmetricTensor<2,dim>::n_independent_components;
          static bool symmetric_component (const unsigned int idx)
          {
            static const  SymmetricTensor<2,dim> t;
            const TableIndices<2> table_indices = t.unrolled_to_component_indices(idx);
            return table_indices[0] != table_indices[1];
          }

          template<typename NumberType>
          using tensor_type = SymmetricTensor<rank,dim,NumberType>;

          template<typename NumberType>
          using value_type = tensor_type<NumberType>;

          template<typename NumberType>
          using gradient_type = Tensor<rank+1,dim,NumberType>;

          template<int rank_in>
          static TableIndices<rank>
          table_index_view (const TableIndices<rank_in> &table_indices,
                            const unsigned int           column_offset)
          {
            Assert(column_offset < rank_in, ExcInternalError());
            Assert(1+column_offset < rank_in, ExcInternalError());
            return TableIndices<rank>(table_indices[column_offset],
                                      table_indices[1+column_offset]);
          }

          template<typename NumberType, typename IndexType = unsigned int, int rank_in>
          static IndexType
          local_component(const TableIndices<rank_in> &table_indices,
                          const unsigned int           column_offset)
          {
            typedef tensor_type<NumberType> TensorType;
            return TensorType::component_to_unrolled_index(table_index_view(table_indices,
                                                           column_offset));
          }
        };

        // --- Define the return types for Gradient calculations ---
        template <int dim, typename NumberType, typename ExtractorType>
        struct Gradient
        {
          typedef typename Extractor<dim,ExtractorType>::template tensor_type<NumberType> type;
        };

        // --- Define the return types for Hessian calculations ---

        // TODO: Can this be made more elegant?
        template<typename ExtractorType_Row,typename ExtractorType_Col>
        struct HessianType
        {
          // Note: We set the return type for
          // HessianType<FEExtractor::Vector,FEExtractor::Vector>
          // as a normal Tensor. This is because if one has two vector components,
          // the coupling tensor (i.e. hessian component<FE::V_1,FE::V_2>) is in
          // general not symmetric.
          template<int rank, int dim, typename NumberType>
          using type = Tensor<rank,dim,NumberType>;
        };
        template<>
        struct HessianType<FEValuesExtractors::SymmetricTensor<2>, FEValuesExtractors::Scalar>
        {
          template<int /*rank*/, int dim, typename NumberType>
          using type = SymmetricTensor<2 /*rank*/,dim,NumberType>;
        };
        template<>
        struct HessianType< FEValuesExtractors::Scalar, FEValuesExtractors::SymmetricTensor<2> >
        {
          template<int /*rank*/, int dim, typename NumberType>
          using type = SymmetricTensor<2 /*rank*/,dim,NumberType>;
        };
        template<>
        struct HessianType< FEValuesExtractors::SymmetricTensor<2>, FEValuesExtractors::SymmetricTensor<2> >
        {
          template<int /*rank*/, int dim, typename NumberType>
          using type = SymmetricTensor<4 /*rank*/,dim,NumberType>;
        };

        template <int dim,typename NumberType,
                  typename ExtractorType_Row,typename ExtractorType_Col>
        struct Hessian
        {
          static const int rank
            = Extractor<dim,ExtractorType_Row>::rank
              + Extractor<dim,ExtractorType_Col>::rank;

          typedef typename HessianType<ExtractorType_Row,ExtractorType_Col>
          ::template type<rank,dim,NumberType> type;

          // typedef typename HessianType<ExtractorType_Row,ExtractorType_Col>::template
          //          type< Extractor<dim,ExtractorType_Row>::rank
          //               +Extractor<dim,ExtractorType_Col>::rank,
          //                dim,NumberType> type;
        };

        // --- Generic functions (Many of these could be integrated into the
        //     extractors helper structs) --

        static inline
        unsigned int
        first_component (const FEValuesExtractors::Scalar &extractor)
        {
          return extractor.component;
        }
        static inline
        unsigned int
        first_component (const FEValuesExtractors::Vector &extractor)
        {
          return extractor.first_vector_component;
        }
        static inline
        unsigned int
        first_component (const FEValuesExtractors::Tensor<1> &extractor)
        {
          return extractor.first_tensor_component;
        }
        static inline
        unsigned int
        first_component (const FEValuesExtractors::Tensor<2> &extractor)
        {
          return extractor.first_tensor_component;
        }
        static inline
        unsigned int
        first_component (const FEValuesExtractors::SymmetricTensor<2> &extractor)
        {
          return extractor.first_tensor_component;
        }

        template<int dim, typename IndexType = unsigned int, typename ExtractorType>
        std::vector<IndexType>
        extract_index_set (const ExtractorType &extractor,
                           const bool           /*ignore_symmetries*/ = true)
        {
          const IndexType n_components = internal::Extractor<dim,ExtractorType>::n_components;
          const IndexType comp_first = first_component(extractor);
          // const unsigned int comp_last = comp_begin + n_components;
          std::vector<IndexType> indices (n_components);
          std::iota(indices.begin(), indices.end(), comp_first);
          return indices;
        }

        template<int dim, typename IndexType = unsigned int>
        std::vector<IndexType>
        extract_index_set (const FEValuesExtractors::SymmetricTensor<dim> &extractor_symm_tensor,
                           const bool                                      ignore_symmetries = true)
        {
          const IndexType n_components = internal::Extractor<dim,FEValuesExtractors::SymmetricTensor<dim> >::n_components;
          if (ignore_symmetries == true)
            {
              const IndexType comp_first = first_component(extractor_symm_tensor);
              std::vector<IndexType> indices (n_components);
              std::iota(indices.begin(), indices.end(), comp_first);
              return indices;
            }
          else
            {
              // First get all of the indices of the non-symmetric tensor
              const FEValuesExtractors::Tensor<dim> extractor_tensor (extractor_symm_tensor.first_tensor_component);
              std::vector<IndexType> indices = extract_index_set<dim>(extractor_tensor, true);

              // Then we overwrite any illegal entries with the equivalent indices
              // from the symmetric tensor
              for (unsigned int i=0; i<indices.size(); ++i)
                {
                  if (indices[i] >= n_components)
                    {
                      const TableIndices<2> ti_tensor = Tensor<2,dim>::unrolled_to_component_indices(indices[i]);
                      const IndexType sti_new_index = SymmetricTensor<2,dim>::component_to_unrolled_index(ti_tensor);
                      indices[i] = sti_new_index;
                    }
                }

              return indices;
            }
        }

        template<int rank, int dim, typename NumberType>
        TableIndices<rank>
        get_table_indices (const Tensor<rank,dim,NumberType> &t,
                           const unsigned int                &idx)
        {
          return t.unrolled_to_component_indices(idx);
        }

        template<int rank, int dim, typename NumberType>
        TableIndices<rank>
        get_table_indices (const SymmetricTensor<rank,dim,NumberType> &t,
                           const unsigned int                         &idx)
        {
          return t.unrolled_to_component_indices(idx);
        }

        template<typename TensorType, typename NumberType>
        void
        set_tensor_entry(TensorType         &t,
                         const unsigned int &idx,
                         const NumberType   &value)
        {
          // Where possible, set values using TableIndices
          Assert(idx < t.n_independent_components,
                 ExcMessage("Index out of range"));
          t[get_table_indices(t,idx)] = value;
        }

        template<int dim, typename NumberType>
        void
        set_tensor_entry(SymmetricTensor<4,dim,NumberType> &t,
                         const unsigned int                &idx_row,
                         const unsigned int                &idx_col,
                         const NumberType                  &value)
        {
          // Fourth order symmetric tensors require a specialized interface
          // to extract values.
          // NOTE: Don't templatize this on NumberType
          // We have to get rid of all instances of adtl::adouble
          // before calling adtl::setNumDir.
          static const SymmetricTensor<2,dim,double> dummy = SymmetricTensor<2,dim,double>();
          Assert(idx_row < dummy.n_independent_components,
                 ExcMessage("Index out of range"));
          Assert(idx_col < dummy.n_independent_components,
                 ExcMessage("Index out of range"));
          const TableIndices<2> indices_row = get_table_indices(dummy,idx_row);
          const TableIndices<2> indices_col = get_table_indices(dummy,idx_col);
          t[indices_row[0]][indices_row[1]][indices_col[0]][indices_col[1]] = value;
        }

        // Specialisations
        template<int dim, typename NumberType>
        void
        set_tensor_entry(Tensor<0,dim,NumberType> &t,
                         const unsigned int       &idx,
                         const NumberType         &value)
        {
          Assert(idx==0,
                 ExcMessage("Index for rank-0 tensor out of range"));
          (void) idx;
          t = value;
        }

        template<typename NumberType>
        void
        set_tensor_entry(NumberType         &t,
                         const unsigned int &idx,
                         const NumberType   &value)
        {
          Assert(idx==0,
                 ExcMessage("Index for scalar value out of range"));
          (void) idx;
          t = value;
        }

        template<int rank, int dim, typename NumberType,
                 template<int, int, typename> class TensorType>
        NumberType
        get_tensor_entry(const TensorType<rank,dim,NumberType> &t,
                         const unsigned int                    &idx)
        {
          // Where possible, get values using TableIndices
          Assert(idx < t.n_independent_components,
                 ExcMessage("Index out of range"));
          return t[get_table_indices(t,idx)];
        }

        template<int dim, typename NumberType,
                 template<int, int, typename> class TensorType>
        NumberType
        get_tensor_entry(const TensorType<0,dim,NumberType> &t,
                         const unsigned int                 &idx)
        {
          Assert(idx==0,
                 ExcMessage("Index for rank-0 tensor out of range"));
          return t;
        }

        template<typename NumberType>
        const NumberType &
        get_tensor_entry(const NumberType   &t,
                         const unsigned int &idx)
        {
          Assert(idx==0,
                 ExcMessage("Index for scalar value out of range"));
          return t;
        }

        template<int rank, int dim, typename NumberType,
                 template<int, int, typename> class TensorType>
        NumberType &
        get_tensor_entry(TensorType<rank,dim,NumberType> &t,
                         const unsigned int              &idx)
        {
          // Where possible, get values using TableIndices
          Assert(idx < t.n_independent_components,
                 ExcMessage("Index out of range"));
          return t[get_table_indices(t,idx)];
        }

        template<int dim, typename NumberType,
                 template<int, int, typename> class TensorType>
        NumberType &
        get_tensor_entry(TensorType<0,dim,NumberType> &t,
                         const unsigned int           &idx)
        {
          Assert(idx==0,
                 ExcMessage("Index for rank-0 tensor out of range"));
          (void) idx;
          return t;
        }

        template<typename NumberType>
        NumberType &
        get_tensor_entry(NumberType         &t,
                         const unsigned int &idx)
        {
          Assert(idx==0,
                 ExcMessage("Index for scalar value out of range"));
          (void) idx;
          return t;
        }

      }
    } // namespace internal



    /**
     * A base helper class that facilitates the evaluation of the derivitive
     * of a number of user-defined dependent variables $\mathbf{f}(\mathbf{X})$
     * with respect to a set of independent variables $\mathbf{X}$.
     *
     * @warning Adol-C does not support the standard threading models used by
     * deal.II, so this class should \b not be embedded within a multithreaded
     * function. It is, however, suitable for use in both serial and MPI routines.
     *
     * In addition to the dimension @p dim, this class is templated on the
     * floating point type @p scalar_type of the number that we'd like to
     * differentiate, as well as an enumeration indicating the @p ADNumberTypeCode .
     *
     * @author Jean-Paul Pelteret, 2016, 2017
     */
    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType = double>
    class ADHelperBase
    {
    public:

      /**
      * Type definition for the floating point number type that are used in, and
      * result from, all computations.
      */
      typedef typename AD::NumberTraits<ScalarType,ADNumberTypeCode>::scalar_type scalar_type;

      /**
      * Type definition for the auto-differentiation number type that is used
      * in all computations.
      */
      typedef typename AD::NumberTraits<ScalarType,ADNumberTypeCode>::ad_type     ad_type;

      /**
       * @name Constructor / destructor
       */
//@{

      /**
      * The constructor for the class.
      *
      * @param[in] n_independent_variables The number of independent variables
      * that will be used in the definition of the functions with which one
      * desires to compute the sensitivities of.
      * @param[in] n_dependent_variables Number of scalar functions defined
      * with a sensitivity to the given independent variables.
      */
      ADHelperBase(const unsigned int n_independent_variables,
                   const unsigned int n_dependent_variables);

      /**
      * Destructor
      */
      virtual
      ~ADHelperBase();

//@}

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
      * Returns the tape number which is currently activated for recording or
      * reading.
      */
      int
      active_tape() const;

      /**
      * Prints the values currently assigned to the independent variables.
      *
      * @param[in] stream The output stream to which the values are to be written.
      */
      void
      print_values(std::ostream &stream) const;

      /**
      * Prints the statistics regarding the usage of the tapes.
      *
      * @param[in] stream The output stream to which the values are to be written.
      */
      void
      print_tape_stats(std::ostream &stream,
                       const unsigned int &tape_index) const;

//@}

      /**
       * @name Recording tapes
       */
//@{

      /**
      * Reset the state of the helper class.
      *
      * When an instance of an ADHelperBase is stored as a class member object
      * with the intention to reuse its instance, it is necessary to reset the
      * state of the object before use. This is because, internally, there is
      * error checking performed to ensure that the correct auto-differentiable
      * data is being tracked and used only when appropriate. This function
      * clears all member data and, therefore, allows the state of all internal
      * flags to be safely reset to their initial state.
      *
      * In the rare case that the number of independent or dependent variables
      * has changed, this can also reconfigured by passing in the appropriate
      * arguments to the function.
      *
      * @note This also resets the active tape number to an invalid number, and
      * deactivates the recording mode for taped variables.
      */
      virtual void
      reset (const unsigned int n_independent_variables = 0,
             const unsigned int n_dependent_variables   = 0);

      /**
      * Specify the number of @p independent_variables to be used in tapeless
      * mode.
      *
      * Although this function is called internally in the ADHelperBase
      * constructor, there may be occasions when adtl::adoubles are created
      * before an instance of this class is created. This function therefore
      * allows one to declare at the earliest possible instance how many
      * directional derivatives will be considered in tapeless mode.
      *
      * @warning Calling this function leaves the set number of directional
      * derivatives in a persistent state, so it will not be possible to
      * further modify during course of the program's execution.
      */
      static void
      configure_tapeless_mode (const unsigned int &n_independent_variables);

      /**
      * Select a tape to record to or read from.
      *
      * @param[in] tape_index The index of the tape to be written
      *
      * @note The chosen tape index must be greater than zero and less than
      * max_tape_index.
      */
      void
      activate_tape(const unsigned int &tape_index);

      /**
      * Enable recording mode for a given tape.
      *
      * The operations that take place between this function call and that
      * of disable_record_sensitivities are recorded to the tape and can
      * be replayed and reevaluated as necessary.
      *
      * The typical set of operations to be performed during this phase are:
      *   - Definition of some dummy independent variables via
      *     register_independent_variable(). These define the branch of operations
      *     tracked by the tape
      *   - Extraction of a set of independent variables of auto-differentiable
      *     type using get_sensitive_variables(). These are then tracked during
      *     later computations.
      *   - Defining the dependent variables via register_dependent_variable().
      *     These are the
      *
      * @param[in] tape_index The index of the tape to be written
      * @param[in] overwrite_tape Express whether tapes are allowed to be reused
      * @param[in] keep Determines whether the numerical values of all independent
      *                 variables are recorded in the tape buffer. If true, then
      *                 tape can be immediately used to perform computations after
      *                 recording is complete.
      *
      * @note During the recording phase, no value(), gradient() or hessian()
      *       operations can be performed.
      */
      bool
      enable_record_sensitivities(const unsigned int &tape_index,
                                  const bool          overwrite_tape = false,
                                  const bool          keep = true);

      /**
      * Disable recording mode for a given tape.
      *
      * @note After this functiom call, the tape is considered ready for use and
      * operations such as value(), gradient() or hessian() can be executed.
      *
      * @note This operation is only valid in recording mode.
      */
      void
      disable_record_sensitivities(const bool write_tapes_to_file = false);

    protected:

      /**
       * @name Taping
       */
//@{

      /**
      * A tape index that is, in general, reserved by Adol-C and unusable to be
      * safely used under all conditions.
      *
      * Adol-C doesn't allow us to record to this tape (i.e. can't write it to
      * file), so we can safely use it as an invalidation case. In general, we
      * want the user to be able to record to a tape if they'd like.
      */
      static const int invalid_tape_index = 0;

      /**
      * The maximum number of tapes that can be written on one process
      */
      static const int max_tape_index = 32; // TODO: Double check whether this is the maximum number of tapes allowed

      /*-
      * Index of the tape that is currently in use. It is this tape that will be
      * recorded to or read from when performing computations.
      */
      int              active_tape_index;

      /*-
      * A collection of tapes that have been recorded to on this process.
      *
      * It is important to keep track of this so that one doesn't accidentally
      * record over a tape (unless specifically instructed to) and that one
      * doesn't try to use a tape that doesn't exist.
      */
      std::set<int>    registered_tapes;

      /**
      * Mark whether we're we're going to tell Adol-C to keep the values stored
      * on the tape so that they can be evaluated again at a later stage.
      */
      bool             keep_values;

      /**
      * Mark whether we're currently recording a tape. Dependent on the state of
      * this flag, only a restricted set of operations are allowable.
      */
      bool             is_recording;

//@}

      /**
       * @name Independent variables
       */
//@{

      /**
      * A set of independent variables $\mathbf{X}$ that we will be performing
      * symbolic differentiation with respect to.
      *
      * The gradients and hessians of dependent variables will be computed
      * at these finite values.
      */
      mutable std::vector<scalar_type> independent_variable_values;

      /**
      * A set of sensitive independent variables $\mathbf{X}$ that we will be
      * performing symbolic differentiation with respect to.
      *
      * The gradients and hessians of dependent variables will be computed
      * using these configured AD numbers. Note that only reverse-mode AD
      * requires that the sensitive independent variables be stored.
      */
      mutable std::vector<ad_type>     independent_variables;

      /**
      * Registered independent variables that have been manipulated for a given
      * set of operations.
      */
      std::vector<bool>                touched_independent_variables;

      /**
      * Registered independent variables that have been extracted and their
      * sensitivities marked.
      */
      mutable std::vector<bool>        touched_marked_independent_variables;

      /**
      * Reset the boolean vector that indicates which independent varaibles
      * we've been manipulating for the current set of operations
      */
      void
      reset_touched_independent_variables ();

      /**
      * Set the actual value of the independent variable $X_{i}$.
      *
      * @param[in] value The value to set the index'd independent variable to.
      * @param[in] index The index in the vector of independent variables.
      */
      void
      set_sensitivity_value (const scalar_type    &value,
                             const unsigned int  index);

      /**
      * Initialise an independent variable $X_{i}$ and record it to the tape.
      *
      * @note Care must be taken to mark each independent variable only once.
      *
      * @note The order in which the independent variables are marked defined the
      * order of all future internal operations. They must be manipulated in the
      * same order as that in which they are first marked. If not then Adol-C
      * won't throw an error, but rather it might complain nonsensically during
      * later computations or produce garbage results.
      *
      * @param[out] out An auto-differentiable number that is ready for use in
      * computations. The operations that are performed with it are recorded on
      * the tape and will be considered in the computation of dependent variable
      * values.
      * @param[in] index The index in the vector of independent variables.
      */
      void
      mark_independent_variable (ad_type            &out,
                                 const unsigned int  index) const;

      /**
       * Finalize the state of the independent variables before use.
       *
       * This step and the storage of the independent variables is done
       * separately because some derived classes may offer the capability
       * to add independent varialbes in a staggered manner. This function
       * is to be triggered when these values are considered finalized
       * and we can safely initialize the sensitive equivalents of those
       * values.
       */
      void
      finalize_sensitive_independent_variables() const;

      /**
      * Initialise an independent variable $X_{i}$.
      *
      * @param[out] out An auto-differentiable number that is ready for use in
      * standard computations. The operations that are performed with it are not
      * recorded on the tape, and so should only be used when not in recording
      * mode.
      * @param[in] index The index in the vector of independent variables.
      */
      void
      get_independent_variable (ad_type            &out,
                                const unsigned int  index) const;

      /**
      * The number of independent variables that have been manipulated within a
      * set of operations.
      */
      unsigned int
      n_touched_independent_variables () const;
//@}

      /**
       * @name Dependent variables
       */
//@{

      /**
      * A set of dependent variables $\mathbolsymbol{\Psi}(\mathbf{X})$ for which
      * we wish to compute the derivatives of.
      *
      * The gradients and hessians of these dependent variables will be computed
      * at the values $\mathbf{X}$ set with the set_sensitivity_value function.
      *
      * @note These are stored as an @p ad_type so that we can use them to
      * compute function values and directional derivatives in the case that
      * tapeless numbers are used
      */
      std::vector<ad_type> dependent_variables;

      /**
      * A list of dependent variables.
      */
      std::vector<bool>    touched_dependent_variables;

      /**
      * The number of dependent variables that have been registered.
      */
      unsigned int
      n_touched_dependent_variables () const;

      /**
      * Register the definition of the index'th dependent variable
      * $\Psi(\mathbf{X})$.
      *
      * @param[in] func The recorded function that defines a dependent variable.
      * @param[in] index The index of the entry in the global list of dependent
      * variables that this function belongs to.
      *
      * @note Each dependent variable must only be registered once
      */
      void
      register_dependent_variable (const ad_type      &func,
                                   const unsigned int  index);
//@}

    }; // class ADHelperBase



    /**
     * A base helper class that facilitates the evaluation of point-wise defined
     * functions. This class is not complete for use since no derivative
     * computations have been implemeneted. What can be computed and how the
     * computed values are interpreted depends on the number of dependent variables
     * to be considered.
     *
     * @author Jean-Paul Pelteret, 2016
     */
    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType = double>
    class ADHelperPointLevelFunctionsBase
      : public ADHelperBase<dim,ADNumberTypeCode,ScalarType>
    {
    public:

      /**
      * Type definition for the floating point number type that are used in, and
      * result from, all computations.
      */
      using scalar_type = typename ADHelperBase<dim,ADNumberTypeCode,ScalarType>::scalar_type;

      /**
      * Type definition for the auto-differentiation number type that is used
      * in all computations.
      */
      using ad_type = typename ADHelperBase<dim,ADNumberTypeCode,ScalarType>::ad_type;

      /**
       * @name Constructor / destructor
       */
      //@{

      /**
      * The constructor for the class.
      *
      * @param[in] n_independent_variables The number of independent variables
      * that the scalar function is sensitive to.
      */
      ADHelperPointLevelFunctionsBase(const unsigned int n_independent_variables,
                                      const unsigned int n_dependent_variables);

      /**
      * Destructor
      */
      virtual
      ~ADHelperPointLevelFunctionsBase();

      //@}

      /**
       * @name Recording tapes
       */
      //@{

      /**
      * Reset the state of the helper class.
      *
      * When an instance of an ADHelperBase is stored as a class member object
      * with the intention to reuse its instance, it is necessary to reset the
      * state of the object before use. This is because, internally, there is
      * error checking performed to ensure that the correct auto-differentiable
      * data is being tracked and used only when appropriate. This function
      * clears all member data and, therefore, allows the state of all internal
      * flags to be safely reset to their initial state.
      *
      * In the rare case that the number of independent or dependent variables
      * has changed, this can also reconfigured by passing in the appropriate
      * arguments to the function.
      *
      * @note This also resets the active tape number to an invalid number, and
      * deactivates the recording mode for taped variables.
      */
      virtual void
      reset (const unsigned int n_independent_variables = 0,
             const unsigned int n_dependent_variables   = 0);

      /**
      * Register the complete set of independent variables $\mathbf{X}$.
      *
      * @param[in] values A field that defines the values of all independent
      * variables. To avoid potential issues with branch switching, it may be a
      * good idea to choose the field's values close to those that will be later
      * evaluated.
      *
      * @note The input value type must correspond to this class's number type.
      * So if this class is templated on type double then the number type
      * associated with value must also be of type double.
      *
      * @note This operation is only valid in recording mode.
      */
      void
      register_independent_variables (const std::vector<scalar_type> &values);

      /**
      * Register the subset of independent variables
      * $\mathbf{A} \in \mathbf{X}$.
      *
      * @param[in] value A field that defines a number of independent variables.
      * To avoid potential issues with branch switching, it may be a good idea
      * to choose the field's values close to those that will be later
      * evaluated.
      * @param[in] extractor An extractor associated with the input field
      * variables. This effectively defines which components of the global set
      * of independent variables this field is associated with.
      *
      * @note The input value type must correspond to this class's number type.
      * So if this class is templated on type double then the number type
      * associated with value must also be of type double.
      *
      * @note The input extractor must correspond to the input value type.
      * So, for example, if a value is a rank-1 tensor, then the extractor must
      * be an FEValuesExtractors::Vector or FEValuesExtractors::Tensor<1>.
      *
      * @note This operation is only valid in recording mode.
      */
      template<typename ValueType, typename ExtractorType>
      void
      register_independent_variable (const ValueType     &value,
                                     const ExtractorType &extractor);

      /**
      * The complete set of independent variables of auto-differentiable type.
      * It is indicated to Adol-C that operations performed with these numbers
      * are to be tracked, so they are considered "sensitive" variables.
      * The values of the components of the returned object are initialised to
      * the values set with register_independent_variable().
      *
      * @return An object of auto-differentiable type numbers.
      *
      * @note This operation is only valid in recording mode.
      */
      const std::vector<ad_type> &
      get_sensitive_variables ();

      /*
      * Extract of a set of independent variables of auto-differentiable type.
      * It is indicated to Adol-C that operations performed with these numbers
      * are to be tracked, so they are considered "sensitive" variables.
      * The values of the components of the returned object are initialised to
      * the values set with register_independent_variable().
      *
      * @param[in] extractor An extractor associated with the input field
      * variables. This effectively defines which components of the global set
      * of independent variables this field is associated with.
      * @return An object of auto-differentiable type numbers. The return type is
      * based on the input extractor, and will be either a scalar, Tensor<1,dim>,
      * Tensor<2,dim>, or SymmetricTensor<2,dim>.
      *
      * @note This operation is only valid in recording mode.
      */
      template<typename ExtractorType>
      typename internal::Extractor<dim,ExtractorType>::template tensor_type<ad_type>
      get_sensitive_variables (const ExtractorType &extractor);

      //@}

      /**
       * @name Post-processing tapes
       */
      //@{

      /*
      * Extract of a set of independent variables of auto-differentiable type.
      * Operations performed with these numbers are not tracked by Adol-C,
      * so they are considered "non-sensitive" variables.
      * The values of the components of the returned object are initialised to
      * the values set with register_independent_variable().
      *
      * @param[in] extractor An extractor associated with the input field
      * variables. This effectively defines which components of the global set
      * of independent variables this field is associated with.
      * @return An object of auto-differentiable type numbers. The return type is
      * based on the input extractor, and will be either a scalar, Tensor<1,dim>,
      * Tensor<2,dim>, or SymmetricTensor<2,dim>.
      *
      * @note This function is not typically used within the context of automatic
      * differentation computations, but can make performing substitutions in
      * symbolic computations easier.
      *
      * @note This operation is only valid outside recording mode.
      */
      template<typename ExtractorType>
      typename internal::Extractor<dim,ExtractorType>::template tensor_type<ad_type>
      get_non_sensitive_variables (const ExtractorType &extractor);

      //@}

      /**
       * @name Computations using the recorded tape at the point defined by the
       * set independent variable values
       */
      //@{

      /**
      * Set the values for the independent variables $\mathbf{X}$.
      *
      * @param[in] values A vector field that defines the values of all
      * independent variables.
      *
      * @note The input value type must correspond to this class's number type.
      * So if this class is templated on type double then the number type
      * associated with value must also be of type double.
      *
      * @note If the keep flag has been set when enable_record_sensitivities() is
      * called, the tape is immediately used after creation, and the values of
      * the independent variables set by register_independent_variable() are those
      * at which the function is to be evaluated, then a call to this function
      * is not strictly necessary.
      */
      void
      set_independent_variables (const std::vector<scalar_type> &values);

      /**
      * Set the values for a subset of independent variables
      * $\mathbf{A} \in \mathbf{X}$.
      *
      * @param[in] value A field that defines the values of a number of
      * independent variables.
      * @param[in] extractor An extractor associated with the input field
      * variables. This effectively defines which components of the global set
      * of independent variables this field is associated with.
      *
      * @note The input value type must correspond to this class's number type.
      * So if this class is templated on type double then the number type
      * associated with value must also be of type double.
      *
      * @note The input extractor must correspond to the input value type.
      * So, for example, if a value is a rank-1 tensor, then the extractor must
      * be an FEValuesExtractors::Vector or FEValuesExtractors::Tensor<1>.
      *
      * @note If the keep flag has been set when enable_record_sensitivities() is
      * called, the tape is immediately used after creation, and the values of
      * the independent variables set by register_independent_variable() are those
      * at which the function is to be evaluated, then a call to this function
      * is not strictly necessary.
      */
      template<typename ValueType, typename ExtractorType>
      void
      set_independent_variable (const ValueType     &value,
                                const ExtractorType &extractor);

      //@}

    protected:

      /**
       * @name Independent variables
       */
      //@{

      /**
      * Set the actual value of the independent variable $X_{i}$.
      *
      * @param[in] value The value to set the index'd independent variable to.
      * @param[in] index The index in the vector of independent variables.
      * @param[in] symmetric_dof mark whether this index relates to a component
      * of a symmetric field.
      */
      void
      set_sensitivity_value (const scalar_type  &value,
                             const unsigned int  index,
                             const bool          symmetric_dof);

      /**
      * The independent variables for which we must take into account symmetry
      * when extracting their gradient or hessian values.
      */
      bool
      is_symmetric_independent_variable (const unsigned int index) const;

      /**
      * The number of independent variables that have been marked as being
      * components of a symmetric field
      */
      unsigned int
      n_symmetric_independent_variables () const;


      //@}

    private:

      /**
       * @name Independent variables
       */
      //@{

      /**
      * The independent variables for which we must take into account symmetry
      * when extracting their gradient or hessian values.
      */
      std::vector<bool> symmetric_independent_variables;

      //@}

    }; // class ADHelperPointLevelFunctionsBase



    /**
     * A helper class that facilitates the evaluation of a scalar function,
     * its first and second derivatives (gradient and hessian). This class
     * would typically be used to compute the first and second derivatives
     * of a <b>stored energy function<\b> defined at a quadrature point.
     * It can also be used to compute derivatives of any other scalar field
     * so long as all its dependence on the independent variables are explicit
     * (that is to say, no independent variables may have some implicit
     * dependence on one another).
     *
     * An example of its usage in the case of a multi-field constitutive laws
     * might be as follows:
     * @code
     * // Define a tape; this pre-records the steps necessary to perform AD
     * // computations. Reusing it may be be a time-saver in some scenarios.
     * // Instead of having to perform these operations each time a AD helper is
     * // used, one could run them once up front and then recall a set of taped
     * // operations as necessary.
     * const int tape_no = 1;
     * const bool is_recording = ad_helper.enable_record_sensitivities(tape_no, //  material_id
     *                                                                 true,    // overwrite_tape
     *                                                                 true);   // keep
     *
     * // Define some extractors that will help us set independent variables and
     * // and later get the computed values related to the dependent variables
     * const FEValuesExtractors::SymmetricTensor<2> C_dofs (0);
     * const FEValuesExtractors::Vector             H_dofs (dealii::SymmetricTensor<2,dim>::n_independent_components);
     * const unsigned int n_independent_variables = SymmetricTensor<2,dim>::n_independent_components
     *                                            + Tensor<1,dim>::n_independent_components;
     *
     * // Define the helper that we will use in the AD computations for our scalar
     * // energy function. Note that we expect it to return values of type double.
     * ADHelperPointLevelFunctionsBase<dim,double> ad_helper (n_independent_variables);
     *
     * // Fields that provide the independent values.
     * // When the tape is being played, these should be set to Something
     * // meaningful.
     * const Tensor<1,dim> H = ...;
     * const SymmetricTensor<2,dim> C = ...;
     *
     * if (is_recording == true)
     * {
     *   // These could happily be set to anything, unless the function will
     *   // be evaluated along a branch not otherwise traversed during later
     *   // use. For this reason, in this example instead of using some dummy
     *   // values, we'll actually map out the function at the same point around
     *   // which we'll later linearise it.
     *   ad_helper.register_independent_variable(H, H_dofs);
     *   ad_helper.register_independent_variable(C, C_dofs);
     *
     *   // NOTE: We have to extract the sensitivities in the order we wish to
     *   //       introduce them. So this means we have to do it by logical order
     *   //       of the extractors that we've created.
     *   const SymmetricTensor<2,dim,ADNumberType> C_AD = ad_helper.get_sensitive_variables(C_dofs);
     *   const Tensor<1,dim,ADNumberType>          H_AD = ad_helper.get_sensitive_variables(H_dofs);
     *
     *   // Here we define the material stored energy function.
     *   // Insert your own fancy function here...
     *   ADNumberType psi = 0.5*(1.0+(H_AD*H_AD)/100.0)*(trace(C_AD) - dim);
     *
     *   ad_helper.register_dependent_variable(psi_CH);
     *   ad_helper.disable_record_sensitivities(false); // write_tapes_to_file
     * }
     * else
     * {
     *   Assert(is_recording==true, ExcInternalError());
     * }
     *
     * // Now we do some work...
     * // First we activate a tape and set the values of the
     * // independent variables. Note that since we've set the keep
     * // flag and their values have not changed, it is not strictly
     * // necessary to reset the values of the independent variables.
     * ad_helper.activate_tape(tape_no);
     * ad_helper.set_independent_variable(C, C_dofs);
     * ad_helper.set_independent_variable(H, H_dofs);
     *
     * // Play the tape and record the output function value, its gradient and
     * // linearisation. These are expensive to compute, so we'll do this once
     * // and extract the desired values from tme
     * const double psi               = ad_helper.compute_value();
     * const Vector<double> Dpsi      = ad_helper.compute_gradient();
     * const FullMatrix<double> D2psi = ad_helper.compute_hessian();
     *
     * // Extract the desired components of the gradient vector and Hessian matrix
     * const SymmetricTensor<2,dim> S = 2.0*ad_helper.extract_gradient_component(Dpsi,C_dofs);
     * const SymmetricTensor<4,dim> H = 4.0*ad_helper.extract_hessian_component(D2psi,C_dofs,C_dofs);
     * @endcode
     *
     * @warning Adol-C does not support the standard threading models used by
     * deal.II, so this class should \b not be embedded within a multithreaded
     * function. It is, however, suitable for use in both serial and MPI routines.
     *
     * @author Jean-Paul Pelteret, 2016
     */
    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType = double>
    class ADHelperScalarFunction
      : public ADHelperPointLevelFunctionsBase<dim,ADNumberTypeCode,ScalarType>
    {
    public:

      /**
      * Type definition for the floating point number type that are used in, and
      * result from, all computations.
      */
      using scalar_type = typename ADHelperBase<dim,ADNumberTypeCode,ScalarType>::scalar_type;

      /**
      * Type definition for the auto-differentiation number type that is used
      * in all computations.
      */
      using ad_type = typename ADHelperBase<dim,ADNumberTypeCode,ScalarType>::ad_type;

      /**
       * @name Constructor / destructor
       */
//@{

      /**
      * The constructor for the class.
      *
      * @param[in] n_independent_variables The number of independent variables
      * that the scalar function is sensitive to.
      */
      ADHelperScalarFunction(const unsigned int n_independent_variables);

      /**
      * Destructor
      */
      virtual
      ~ADHelperScalarFunction();

//@}

      /**
       * @name Dependent variables
       */
//@{

      /**
      * Register the definition of the scalar field $\Psi(\mathbf{X})$.
      *
      * @param[in] func The recorded function that defines a dependent variable.
      *
      * @note For this class that expects only one dependent variable, this
      * function must only be called once per tape.
      *
      * @note This operation is only valid in recording mode.
      */
      void
      register_dependent_variable (const ad_type &func);

//@}

      /**
       * @name Computations using the recorded tape at the point defined by the
       * set independent variable values
       */
//@{

      /**
      * Computes the value of the scalar field $\Psi(\mathbf{X})$ using the tape
      * as opposed to executing the source code.
      *
      * @return A scalar object with the value for the scalar field evaluated
      * at the point defined by the independent variable values.
      */
      scalar_type
      compute_value() const;

      /**
      * Computes the gradient (first derivative) of the scalar field with respect
      * to all independent variables, i.e.
      * @f[
      *   \frac{\partial\Psi(\mathbf{X})}{\partial\mathbf{X}}
      * @f]
      *
      * @return A Vector with the values for the scalar field gradient evaluated
      * at the point defined by the independent variable values.
      */
      Vector<scalar_type>
      compute_gradient() const;

      /**
      * Computes the hessian (second derivative)  of the scalar field with respect
      * to all independent variables, i.e.
      * @f[
      *   \frac{\partial^{2}\Psi(\mathbf{X})}{\partial\mathbf{X} \otimes \partial\mathbf{X}}
      * @f]
      *
      * @return A FullMatrix with the values for the scalar field gradient evaluated
      * at the point defined by the independent variable values.
      */
      FullMatrix<scalar_type>
      compute_hessian() const;

      /**
      * Extract the function gradient for a subset of independent variables
      * $\mathbf{A} \in \mathbf{X}$, i.e.
      * @f[
      *   \frac{\partial\Psi(\mathbf{X})}{\partial\mathbf{A}}
      * @f]
      *
      * @param[in] gradient The gradient of the scalar function with respect to
      * all independent variables, i.e. that returned by compute_gradient().
      * @param[in] extractor_row An extractor associated with the input field
      * variables. This effectively defines which components of the global set
      * of independent variables this field is associated with.
      */
      template<typename ExtractorType_Row>
      typename internal::Gradient<dim,scalar_type,ExtractorType_Row>::type
      extract_gradient_component(const Vector<scalar_type> &gradient,
                                 const ExtractorType_Row   &extractor_row) const;

      /**
      * Extract the function hessian for a subset of independent variables
      * $\mathbf{A},\mathbf{B} \in \mathbf{X}$, i.e.
      * @f[
      *   \frac{}{\partial\mathbf{B}} \left[ \frac{\partial\Psi(\mathbf{X})}{\partial\mathbf{A}} \right]
      *   = \frac{\partial^{2}\Psi(\mathbf{X})}{\partial\mathbf{B} \otimes \partial\mathbf{A}}
      * @f]
      *
      * @param[in] hessian The hessian of the scalar function with respect to
      * all independent variables, i.e. that returned by compute_hessian().
      * @param[in] extractor_row An extractor associated with the input field
      * variables for which the first index of the hessian is extracted.
      * @param[in] extractor_col An extractor associated with the input field
      * variables for which the second index of the hessian is extracted.
      */
      template<typename ExtractorType_Row, typename ExtractorType_Col>
      typename internal::Hessian<dim,scalar_type,ExtractorType_Row,ExtractorType_Col>::type
      extract_hessian_component(const FullMatrix<scalar_type> &hessian,
                                const ExtractorType_Row       &extractor_row,
                                const ExtractorType_Col       &extractor_col) const;

      /**
      * Extract the function hessian for a subset of independent variables
      * $\mathbf{A},\mathbf{B} \in \mathbf{X}$, i.e.
      * @f[
      *   \frac{}{\partial\mathbf{B}} \left[ \frac{\partial\Psi(\mathbf{X})}{\partial\mathbf{A}} \right]
      * @f]
      *
      * This function is a specialisation of the above for rank-0 tensors (scalars)
      */
      Tensor<0,dim,scalar_type>
      extract_hessian_component(const FullMatrix<scalar_type>    &hessian,
                                const FEValuesExtractors::Scalar &extractor_row,
                                const FEValuesExtractors::Scalar &extractor_col) const;

      /**
      * Extract the function hessian for a subset of independent variables
      * $\mathbf{A},\mathbf{B} \in \mathbf{X}$, i.e.
      * @f[
      *   \frac{}{\partial\mathbf{B}} \left[ \frac{\partial\Psi(\mathbf{X})}{\partial\mathbf{A}} \right]
      * @f]
      *
      * This function is a specialisation of the above for rank-4 symmetric tensors
      */
      SymmetricTensor<4,dim,scalar_type>
      extract_hessian_component(const FullMatrix<scalar_type>                &hessian,
                                const FEValuesExtractors::SymmetricTensor<2> &extractor_row,
                                const FEValuesExtractors::SymmetricTensor<2> &extractor_col) const;

//@}

    }; // class ADHelperScalarFunction



    /**
     * A helper class that facilitates the evaluation of a vector of functions,
     * typically one that represents a collection of coupled, multi-dimensional
     * fields.
     * This class would typically be used to compute the linearisation a set of
     * kinetic field variables defined at the quadrature point level.
     *
     * An example of its usage in the case of a quadrature point data linearisation
     * might be as follows:
     * @code
     * @endcode
     *
     * @warning Adol-C does not support the standard threading models used by
     * deal.II, so this class should \b not be embedded within a multithreaded
     * function. It is, however, suitable for use in both serial and MPI routines.
     *
     * @author Jean-Paul Pelteret, 2016
     */
    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType = double>
    class ADHelperVectorFunction
      : public ADHelperPointLevelFunctionsBase<dim,ADNumberTypeCode,ScalarType>
    {
    public:

      /**
      * Type definition for the floating point number type that are used in, and
      * result from, all computations.
      */
      using scalar_type = typename ADHelperBase<dim,ADNumberTypeCode,ScalarType>::scalar_type;

      /**
      * Type definition for the auto-differentiation number type that is used
      * in all computations.
      */
      using ad_type = typename ADHelperBase<dim,ADNumberTypeCode,ScalarType>::ad_type;

      /**
       * @name Constructor / destructor
       */
//@{

      /**
      * The constructor for the class.
      *
      * @param[in] n_independent_variables The number of independent variables
      * that the scalar function is sensitive to.
      */
      ADHelperVectorFunction(const unsigned int n_independent_variables,
                             const unsigned int n_dependent_variables);

      /**
      * Destructor
      */
      virtual
      ~ADHelperVectorFunction();

//@}

      /**
       * @name Dependent variables
       */
//@{

      /**
      * Register the definition of the vector field $\boldsymbol{\Psi}(\mathbf{X})$.
      *
      * @param[in] func A vector of recorded functions that defines the dependent
      * variables.
      *
      * @note For this class that expects only vector field of dependent
      * variables, this function must only be called once per tape.
      *
      * @note This operation is only valid in recording mode.
      */
      void
      register_dependent_variables (const std::vector<ad_type> &funcs);

      /**
      * Register the definition of the scalar field $\Psi(\mathbf{X})$.
      *
      * @param[in] funcs The recorded functions that define a set of dependent variable.
      * @param[in] extractor An extractor associated with the input field
      * variables. This effectively defines which components of the global set
      * of dependent variables this field is associated with.
      *
      * @note The input extractor must correspond to the input value type.
      * So, for example, if a value is a rank-1 tensor, then the extractor must
      * be an FEValuesExtractors::Vector or FEValuesExtractors::Tensor<1>.
      *
      * @note This operation is only valid in recording mode.
      */
      template<typename ValueType, typename ExtractorType>
      void
      register_dependent_variable (const ValueType     &funcs,
                                   const ExtractorType &extractor);

//@}

      /**
       * @name Computations using the recorded tape at the point defined by the
       * set independent variable values
       */
//@{

      /**
      * Computes the value of the vector field $\mathbf{f} = \boldsymbol{\Psi}(\mathbf{X})$
      * using the tape as opposed to executing the source code.
      *
      * @return A Vector object with the value for each component of the vector
      * field evaluated at the point defined by the independent variable values.
      */
      Vector<scalar_type>
      compute_values() const;

      /**
      * Computes the jacobian (first derivative) of the vector field with respect
      * to all independent variables, i.e.
      * @f[
      *   \mathbf{J}(\boldsymbol{\Psi})
      *      = \frac{\partial\boldsymbol{\Psi}(\mathbf{X})}{\partial\mathbf{X}}
      * @f]
      *
      * @return A FullMatrix with the gradient of each component of the vector
      * field evaluated at the point defined by the independent variable values.
      */
      FullMatrix<scalar_type>
      compute_jacobian() const;


      /**
      * Extract the set of function' values for a subset of independent variables
      * $\mathbf{A} \in \mathbf{X}$, i.e.
      * @f[
      *   \mathbf{f}(\boldsymbol{\Psi})
      *      = \frac{\partial\boldsymbol{\Psi}(\mathbf{X})}{\partial\mathbf{A}}
      * @f]
      *
      * @param[in] gradient The gradient of the scalar function with respect to
      * all independent variables, i.e. that returned by compute_gradient().
      * @param[in] extractor_row An extractor associated with the input field
      * variables. This effectively defines which components of the global set
      * of independent variables this field is associated with.
      */
      template<typename ExtractorType_Row>
      typename internal::Gradient<dim,scalar_type,ExtractorType_Row>::type
      extract_value_component(const Vector<scalar_type> &values,
                              const ExtractorType_Row   &extractor_row) const;

      /**
      * Extract the set of functions' jacobian for a subset of independent
      * variables $\mathbf{A} \in \mathbf{X}$, i.e.
      * @f[
      *   \mathbf{J}(\boldsymbol{\Psi})
      *      = \frac{\partial\boldsymbol{\Psi}(\mathbf{X})}{\partial\mathbf{A}}
      * @f]
      *
      * @param[in] jacobian The jacobian of the vector function with respect to
      * all independent variables, i.e. that returned by compute_jacobian().
      * @param[in] extractor_row An extractor associated with the input field
      * variables for which the first index of the jacobian is extracted.
      * @param[in] extractor_col An extractor associated with the input field
      * variables for which the second index of the jacobian is extracted.
      */
      template<typename ExtractorType_Row, typename ExtractorType_Col>
      typename internal::Hessian<dim,scalar_type,ExtractorType_Row,ExtractorType_Col>::type
      extract_jacobian_component(const FullMatrix<scalar_type> &jacobian,
                                 const ExtractorType_Row       &extractor_row,
                                 const ExtractorType_Col       &extractor_col) const;

      /**
      * Extract the set of functions' jacobian for a subset of independent
      * variables $\mathbf{A},\mathbf{B} \in \mathbf{X}$, i.e.
      * @f[
      *   \mathbf{J}(\boldsymbol{\Psi})
      *      = \frac{}{\partial\mathbf{B}} \left[ \frac{\partial\boldsymbol{\Psi}(\mathbf{X})}{\partial\mathbf{A}} \right]
      * @f]
      *
      * This function is a specialisation of the above for rank-0 tensors (scalars)
      */
      Tensor<0,dim,scalar_type>
      extract_jacobian_component(const FullMatrix<scalar_type>    &jacobian,
                                 const FEValuesExtractors::Scalar &extractor_row,
                                 const FEValuesExtractors::Scalar &extractor_col) const;

      /**
      * Extract the set of functions' jacobian for a subset of independent
      * variables $\mathbf{A},\mathbf{B} \in \mathbf{X}$, i.e.
      * @f[
      *   \mathbf{J}(\boldsymbol{\Psi})
      *      = \frac{}{\partial\mathbf{B}} \left[ \frac{\partial\boldsymbol{\Psi}(\mathbf{X})}{\partial\mathbf{A}} \right]
      * @f]
      *
      * This function is a specialisation of the above for rank-4 symmetric tensors
      */
      SymmetricTensor<4,dim,scalar_type>
      extract_jacobian_component(const FullMatrix<scalar_type>                &jacobian,
                                 const FEValuesExtractors::SymmetricTensor<2> &extractor_row,
                                 const FEValuesExtractors::SymmetricTensor<2> &extractor_col) const;

//@}

    }; // class ADHelperVectorFunction



    /**
     * A general helper class that facilitates the evaluation of a vector of
     * functions, as well as its first derivatives (their Jacobian).
     * This class would typically be used to compute the linearisation of a
     * set of local nonlinear equations, but can also be used as the basis of
     * the linearisation of the residual vector defined on the level of a finite
     * element.
     *
     * @note When using the cell-level AD methods in 3d and/or with higher
     * order elements, it is incredibly easy to exceed the tape buffer size.
     * These buffer variables dictate the amount of memory allocated to a tape
     * before it is written to file (at a significant performance loss).
     * Therefore, as stated by the manual, it may be desirable to create a file
     * ".adolcrc" in the program run directory and set the buffer size therein.
     * For example, the following settings increase the default buffer size by
     * 128 times:
     * <code>
       "OBUFSIZE" "67108864"
       "LBUFSIZE" "67108864"
       "VBUFSIZE" "67108864"
    "     TBUFSIZE" "67108864"
     * </code>
     * Note that the quotation marks are mandatory.
     *
     * @warning Adol-C does not support the standard threading models used by
     * deal.II, so this class should \b not be embedded within a multithreaded
     * function. It is, however, suitable for use in both serial and MPI routines.
     *
     * @author Jean-Paul Pelteret, 2016
     */
    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType = double>
    class ADHelperCellLevelBase
      : public ADHelperBase<dim,ADNumberTypeCode,ScalarType>
    {
    public:

      /**
      * Type definition for the floating point number type that are used in, and
      * result from, all computations.
      */
      using scalar_type = typename ADHelperBase<dim,ADNumberTypeCode,ScalarType>::scalar_type;

      /**
      * Type definition for the auto-differentiation number type that is used
      * in all computations.
      */
      using ad_type = typename ADHelperBase<dim,ADNumberTypeCode,ScalarType>::ad_type;

      /**
       * @name Constructor / destructor
       */
//@{

      /**
      * The constructor for the class.
      *
      * @param[in] n_independent_variables The number of independent variables
      * that the vector of functions is sensitive to.
      * @param[in] n_dependent_variables The number of dependent variables, i.e.
      * the number of functions in the vector that one wishes to compute the
      * sensitivities of.
      */
      ADHelperCellLevelBase(const unsigned int n_independent_variables,
                            const unsigned int n_dependent_variables);

      /**
      * Destructor
      */
      virtual
      ~ADHelperCellLevelBase();

//@}

      /**
       * @name Recording tapes
       */
//@{

      /**
      * Register the complete set of independent variables $\mathbf{X}$ that
      * represent the local degree-of-freedom values.
      *
      * @param[in] values A field that defines the values of all degrees-of-freedom.
      * To avoid potential issues with branch switching, it may be a good idea to
      * choose these values close to those that will be later evaluated and linearized
      * around.
      *
      * @note The input value type must correspond to this class's number type.
      * So if this class is templated on type double then the number type
      * associated with value must also be of type double.
      *
      * @note This operation is only valid in recording mode.
      */
      void
      register_dof_values (const std::vector<scalar_type> &dof_values);

      /**
      * Register the complete set of independent variables $\mathbf{X}$ that
      * represent the local degree-of-freedom values.
      *
      * @param[in] values A global field from which the values of all independent
      * variables will be extracted. This typically will be the solution vector
      * around at which point a residual vector is to be computed and around which
      * linearisation is to occur. To avoid potential issues with branch
      * switching, it may be a good idea to choose the field's values close to
      * those that will be later evaluated.
      * @param[in] local_dof_indices A vector of degree of freedom indices from
      * which to extract the local degree of freedom values.
      *
      * @note This operation is only valid in recording mode.
      */
      template<typename VectorType>
      void
      register_dof_values (const VectorType                           &values,
                           const std::vector<types::global_dof_index> &local_dof_indices);

      /**
      * The complete set of DoF values of auto-differentiable type.
      * It is indicated to Adol-C that operations performed with these numbers
      * are to be tracked, so they are considered "sensitive" variables.
      * The values of the components of the returned object are initialised to
      * the values set with register_independent_variable().
      *
      * @return An object of auto-differentiable type numbers.
      *
      * @note This operation is only valid in recording mode.
      */
      const std::vector<ad_type> &
      get_sensitive_dof_values ();

//@}

      /**
       * @name Post-processing tapes
       */
//@{

      /*
      * The complete set of DoF values of auto-differentiable type.
      * Operations performed with these numbers are not tracked by Adol-C,
      * so they are considered "non-sensitive" variables.
      * The values of the components of the returned object are initialised to
      * the values set with register_independent_variable().
      *
      * @return An object of auto-differentiable type numbers.
      *
      *
      * @note This function is not typically used within the context of automatic
      * differentation computations, but can make performing substutitions in
      * symbolic computations easier.
      *
      * @note This operation is only valid outside recording mode.
      */
      std::vector<ad_type>
      get_non_sensitive_dof_values ();

//@}

      /**
       * @name Computations using the recorded tape at the point defined by the
       * set independent variable values
       */
//@{

      /**
      * Set the values for the independent variables $\mathbf{X}$.
      *
      * @param[in] values A vector field that defines the values of all
      * independent variables.
      *
      * @note The input value type must correspond to this class's number type.
      * So if this class is templated on type double then the number type
      * associated with value must also be of type double.
      *
      * @note If the keep flag has been set when enable_record_sensitivities() is
      * called, the tape is immediately used after creation, and the values of
      * the independent variables set by register_independent_variable() are those
      * at which the function is to be evaluated, then a call to this function
      * is not strictly necessary.
      */
      void
      set_dof_values (const std::vector<scalar_type> &values);

      /**
      * Set the values for the independent variables $\mathbf{X}$.
      *
      * @param[in] values A vector field from which the values of all
      * independent variables is to be extracted.
      * @param[in] local_dof_indices A vector of degree of freedom indices from
      * which to extract the local degree of freedom values.
      *
      * @note If the keep flag has been set when enable_record_sensitivities() is
      * called, the tape is immediately used after creation, and the values of
      * the independent variables set by register_independent_variable() are those
      * at which the function is to be evaluated, then a call to this function
      * is not strictly necessary.
      */
      template<typename VectorType>
      void
      set_dof_values (const VectorType                           &values,
                      const std::vector<types::global_dof_index> &local_dof_indices);

      /**
      * Computes the value of the vector field $\boldsymbol{\Psi}(\mathbf{X})$
      * using the tape as opposed to executing the source code.
      *
      * @return A Vector object with the value for each component of the vector
      * field evaluated at the point defined by the independent variable values.
      */
      virtual Vector<scalar_type>
      compute_residual() const = 0;

      /**
      * Computes the gradient (first derivative) of the vector field with respect
      * to all independent variables, i.e.
      * @f[
      *   \frac{\partial\boldsymbol{\Psi}(\mathbf{X})}{\partial\mathbf{X}}
      * @f]
      *
      * @return A FullMatrix with the gradient of each component of the vector
      * field evaluated at the point defined by the independent variable values.
      */
      virtual FullMatrix<scalar_type>
      compute_linearization() const = 0;

      //@}

    protected:

      /**
       * @name Dependent variables
       */
//@{

      /**
      * Register the definition of the vector field $\boldsymbol{\Psi}(\mathbf{X})$.
      *
      * @param[in] func A vector of recorded functions that defines the dependent
      * variables.
      *
      * @note For this class that expects only vector field of dependent
      * variables, this function must only be called once per tape.
      *
      * @note This operation is only valid in recording mode.
      */
      void
      register_dependent_variables (const std::vector<ad_type> &funcs);

//@}

    }; // class ADHelperCellLevelBase



    /**
       * A helper class that facilitates the implementation of a generic (incremental)
       * variational formulation from which the computation of the residual vector, as
       * well as its linearisation, is automated. This class would typically be used to
       * derive the residual vector and tangent matrix, defined on the level of a finite
       * element, or a linearized system of equations, starting from a scalar energy
       * functional.
       *
       * An example of its usage in the case of a residual and tangent computations
       * might be as follows:
       * @code
       * @endcode
       *
       * @warning Adol-C does not support the standard threading models used by
       * deal.II, so this class should \b not be embedded within a multithreaded
       * function. It is, however, suitable for use in both serial and MPI routines.
       *
       * @author Jean-Paul Pelteret, 2016
       */
    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType = double>
    class ADHelperVariationalFormulation
      : public ADHelperCellLevelBase<dim,ADNumberTypeCode,ScalarType>
    {
    public:

      /**
      * Type definition for the floating point number type that are used in, and
      * result from, all computations.
      */
      using scalar_type = typename ADHelperBase<dim,ADNumberTypeCode,ScalarType>::scalar_type;

      /**
      * Type definition for the auto-differentiation number type that is used
      * in all computations.
      */
      using ad_type = typename ADHelperBase<dim,ADNumberTypeCode,ScalarType>::ad_type;

      /**
       * @name Constructor / destructor
       */
      //@{

      /**
      * The constructor for the class.
      *
      * @param[in] n_independent_variables The number of independent variables
      * that the vector of functions is sensitive to.
      *
      * @note These is only one dependent variable associated with the total
      * energy attributed to the local finite element.
      */
      ADHelperVariationalFormulation(const unsigned int n_independent_variables);

      /**
      * Destructor
      */
      virtual
      ~ADHelperVariationalFormulation();

      //@}

      /**
       * @name Dependent variables
       */
      //@{

      /**
      * Register the definition of the total potential $\boldsymbol{\Psi}(\mathbf{X})$.
      *
      * @param[in] func A vector of recorded functions that defines the residual.
      * The components of this vector represents the dependent variables.
      *
      * @note For this class that expects only a single scalar dependent
      * variable, this function must only be called once per tape.
      *
      * @note This operation is only valid in recording mode.
      */
      void
      register_energy_functional (const ad_type &energy);

      //@}

      /**
       * @name Computations using the recorded tape at the point defined by the
       * set independent variable values
       */
      //@{

      /**
      * Evaluation of the residual for a chosen set of degree-of-freedom values.
      * Underlying this is the computation of the gradient (first derivative) of
      * the scalar energy field with respect to all independent variables, i.e.
      * @f[
      *   \mathbf{r}(\mathbf{X}) = \frac{\partial\boldsymbol{\Psi}(\mathbf{X})}{\partial\mathbf{X}}
      * @f]
      *
      * @return A Vector object, for which the value for each entry represents the
      * residual value for the corresponding local degree-of freedom.
      */
      virtual Vector<scalar_type>
      compute_residual() const;

      /**
      * Computes the linearization the residual vector around a chosen set of
      * degree-of-freedom values.
      * Underlying this is the computation of the hessian (second derivative) of
      * the vector field with respect to all independent variables, i.e.
      * @f[
      *   \frac{\partial\mathbf{r}(\mathbf{X})}{\partial\mathbf{X}}
      *     = \frac{\partial^{2}\boldsymbol{\Psi}(\mathbf{X})}{\partial\mathbf{X} \otimes \partial\mathbf{X}}
      * @f]
      *
      * @return A FullMatrix representing the linearization of the residual vector.
      */
      virtual FullMatrix<scalar_type>
      compute_linearization() const;

      //@}

    }; // class ADHelperVariationalFormulation



    /**
     * A helper class that facilitates the evaluation of a vector of functions,
     * that represent the local degree-of-freedom values corresponding to the
     * residual vector, as well as its linearisation. This class
     * would typically be used to compute the linearisation of a residual vector
     * defined on the level of a finite element, or for local nonlinear equations.
     *
     * An example of its usage in the case of a residual linearisation
     * might be as follows:
     * @code
     * @endcode
     *
     * @warning Adol-C does not support the standard threading models used by
     * deal.II, so this class should \b not be embedded within a multithreaded
     * function. It is, however, suitable for use in both serial and MPI routines.
     *
     * @author Jean-Paul Pelteret, 2016
     */
    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType = double>
    class ADHelperResidualLinearisation
      : public ADHelperCellLevelBase<dim,ADNumberTypeCode,ScalarType>
    {
    public:

      /**
      * Type definition for the floating point number type that are used in, and
      * result from, all computations.
      */
      using scalar_type = typename ADHelperBase<dim,ADNumberTypeCode,ScalarType>::scalar_type;

      /**
      * Type definition for the auto-differentiation number type that is used
      * in all computations.
      */
      using ad_type = typename ADHelperBase<dim,ADNumberTypeCode,ScalarType>::ad_type;

      /**
       * @name Constructor / destructor
       */
//@{

      /**
      * The constructor for the class.
      *
      * @param[in] n_independent_variables The number of independent variables
      * that the vector of functions is sensitive to.
      * @param[in] n_dependent_variables The number of dependent variables, i.e.
      * the number of functions in the vector that one wishes to compute the
      * sensitivities of.
      */
      ADHelperResidualLinearisation(const unsigned int n_independent_variables,
                                    const unsigned int n_dependent_variables);

      /**
      * Destructor
      */
      virtual
      ~ADHelperResidualLinearisation();

//@}

      /**
       * @name Dependent variables
       */
//@{

      /**
      * Register the definition of the residual vector field $\mathbf{r}(\mathbf{X})$.
      *
      * @param[in] func A vector of recorded functions that defines the residual.
      * The components of this vector represents the dependent variables.
      *
      * @note For this class that expects only vector field of dependent
      * variables, this function must only be called once per tape.
      *
      * @note This operation is only valid in recording mode.
      */
      void
      register_residual_vector (const std::vector<ad_type> &residual);

//@}

      /**
       * @name Computations using the recorded tape at the point defined by the
       * set independent variable values
       */
//@{

      /**
      * Evaluation of the residual for a chosen set of degree-of-freedom values.
      * This corresponds to the computation the vector-valued residual field
      * $\mathbf{r}(\mathbf{X})$ using the tape as opposed to executing the source code.
      *
      * @return A Vector object, for which the value for each entry represents the
      * residual value for the corresponding local degree-of freedom.
      */
      virtual Vector<scalar_type>
      compute_residual() const;

      /**
      * Computes the linearization the residual vector around a chosen set of
      * degree-of-freedom values.
      * Underlying this is the computation of the gradient (first derivative) of
      * the vector field with respect to all independent variables, i.e.
      * @f[
      *   \frac{\partial\mathbf{r}(\mathbf{X})}{\partial\mathbf{X}}
      * @f]
      *
      * @return A FullMatrix representing the linearization of the residual vector.
      */
      virtual FullMatrix<scalar_type>
      compute_linearization() const;

//@}

    }; // class ADHelperResidualLinearisation


  } // namespace AD
} // namespace Differentiation


/* --------------------------- inline and template functions ------------------------- */


#ifndef DOXYGEN

namespace Differentiation
{
  namespace AD
  {

    namespace internal
    {
      namespace
      {

        template<typename ADNumberType>
        static typename std::enable_if<
        ADNumberTraits<ADNumberType>::is_taped>
        ::type
        configure_adtl (const unsigned int)
        {

        }

#ifdef DEAL_II_WITH_ADOLC

        template<typename ADNumberType>
        static typename std::enable_if<
        ADNumberTraits<ADNumberType>::is_tapeless &&
        !is_adolc_number<ADNumberType>::value>
        ::type
        configure_adtl (const unsigned int)
        {

        }

        template<typename ADNumberType>
        static typename std::enable_if<
        ADNumberTraits<ADNumberType>::is_tapeless &&
        is_adolc_number<ADNumberType>::value>
        ::type
        configure_adtl (const unsigned int n_directional_derivatives)
        {
          // Enable vector mode
          // See Adol-C manual section 7.1
          // NOTE: It is critical that this is done for tapeless mode BEFORE
          // any adtl::adouble are created. If this is not done, then we see
          // this scary warning:
          //
          // ADOL-C Warning: Tapeless: Setting numDir could change memory
          // allocation of derivatives in existing adoubles and may lead to
          // erronious results or memory corruption
          //
          // So we use this dummy function to configure this setting before
          // we create and initialize our class data
          const std::size_t n_live_variables = adtl::refcounter::getNumLiveVar();
          if (n_live_variables == 0)
            {
              adtl::setNumDir(n_directional_derivatives);
            }
          else
            {
              // So there are some live active variables floating around. Here we
              // check if we ask to increase the number of number of computable
              // directional derivatives. If this really is necessary then its
              // absolutely vital that there exist no live variables before doing
              // so.
              const std::size_t n_set_directional_derivatives = adtl::getNumDir();
              if (n_directional_derivatives > n_set_directional_derivatives)
                AssertThrow(n_live_variables == 0,
                            ExcMessage("There are currently " +
                                       Utilities::to_string(n_live_variables) + " "
                                       "live adtl::adouble variables in existence. They currently "
                                       "assume " +
                                       Utilities::to_string(n_set_directional_derivatives) + " directional derivatives "
                                       "but you wish to increase this to " +
                                       Utilities::to_string(n_directional_derivatives) + ". \n"
                                       "To safely change (or more specifically in this case, "
                                       "increase) the number of directional derivatives, there "
                                       "must be no tapeless doubles in local/global scope."));
            }
        }

#else

        template<typename ad_type>
        static typename std::enable_if<
        ADNumberTraits<ad_type>::is_tapeless>
        ::type
        configure_adtl (const unsigned int)
        {

        }

#endif

        template<typename ADNumberType>
        static int
        configure_adtl_and_return_tape_index (const unsigned int n_independent_variables,
                                              const int invalid_tape_index)
        {
          configure_adtl<ADNumberType>(n_independent_variables);
          return invalid_tape_index;
        }


        /**
         * Define the active dependent variable when using reverse-mode AD.
         *
         * If there are multiple dependent variables then it is necessary to
         * inform the independent variables, from which the adjoints are computed,
         * which dependent variable they are computing the gradients with respect
         * to. This function broadcasts this information.
         */
        template<typename ADNumberType>
        static typename std::enable_if<
        is_sacado_rad_number<ADNumberType>::value
        >
        ::type
        reverse_mode_dependent_variable_activation (ADNumberType &dependent_variable)
        {
          // Compute all gradients (adjoints) for this
          // reverse-mode Sacado dependent variable.
          // For reverse-mode Sacado numbers it is necessary to broadcast to
          // all independent variables that it is time to compute gradients.
          // For one dependent variable one would just need to all
          // ad_type::Gradcomp(), but since we have a more
          // generic implementation for vectors of dependent variables
          // (vector mode) we default to the complex case.
          ADNumberType::Outvar_Gradcomp(dependent_variable);
        }


        template<typename ADNumberType>
        static typename std::enable_if<
        !is_sacado_rad_number<ADNumberType>::value>
        ::type
        reverse_mode_dependent_variable_activation (ADNumberType &)
        {

        }

      }
    } // namespace internal



// -------------------------- ADHelperBase ----------------------

// -------------------------- ADHelperPointLevelFunctionsBase ----------------------



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    template<typename ValueType, typename ExtractorType>
    void
    ADHelperPointLevelFunctionsBase<dim,ADNumberTypeCode,ScalarType>::register_independent_variable (
      const ValueType     &value,
      const ExtractorType &extractor)
    {
      // This is actually the same thing the set_independent_variable function,
      // in the sense that we simply populate our array of independent values
      // with a meaningful number. However, in this case we need to double check
      // that we're not registering these variables twice
#ifdef DEBUG
      const std::vector<unsigned int> index_set (internal::extract_index_set<dim>(extractor));
      for (unsigned int i=0; i<index_set.size(); ++i)
        {
          Assert(this->touched_independent_variables[index_set[i]] == false,
                 ExcMessage("Overlapping indices for independent variables."));
        }
#endif
      set_independent_variable(value,extractor);
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    template<typename ValueType, typename ExtractorType>
    void
    ADHelperPointLevelFunctionsBase<dim,ADNumberTypeCode,ScalarType>::set_independent_variable (
      const ValueType     &value,
      const ExtractorType &extractor)
    {
      if (ADNumberTraits<ad_type>::is_tapeless == true)
        {
          Assert(this->is_recording==true,
                 ExcMessage("Cannot change the value of an independent variable "
                            "of the tapeless variety."));
        }
      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          Assert(this->active_tape()!=this->invalid_tape_index,
                 ExcMessage("Invalid tape index"));
        }

      const std::vector<unsigned int> index_set (internal::extract_index_set<dim>(extractor));
      for (unsigned int i=0; i<index_set.size(); ++i)
        {
          set_sensitivity_value(internal::get_tensor_entry(value,i),
                                index_set[i],
                                internal::Extractor<dim,ExtractorType>::symmetric_component(i));
        }
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    template<typename ExtractorType>
    typename internal::Extractor<dim,ExtractorType>::template tensor_type<typename ADHelperBase<dim,ADNumberTypeCode,ScalarType>::ad_type>
    ADHelperPointLevelFunctionsBase<dim,ADNumberTypeCode,ScalarType>::get_sensitive_variables (const ExtractorType &extractor)
    {
      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          Assert(this->active_tape()!=this->invalid_tape_index,
                 ExcMessage("Invalid tape index"));
        }

      // If necessary, finalize the internally stored vector of
      // AD numbers that represents the independent variables
      this->finalize_sensitive_independent_variables();
      Assert(this->independent_variables.size()==this->n_independent_variables(),
             ExcInternalError());

      const std::vector<unsigned int> index_set (internal::extract_index_set<dim>(extractor));
      typename internal::Extractor<dim,ExtractorType>::template tensor_type<ad_type> out;

      for (unsigned int i=0; i<index_set.size(); ++i)
        {
          const unsigned int index = index_set[i];
          Assert(index < this->n_independent_variables(), ExcInternalError());
          Assert(this->touched_independent_variables[index] == true, ExcInternalError());
          internal::get_tensor_entry(out,i) = this->independent_variables[index];
        }

      return out;
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    template<typename ExtractorType>
    typename internal::Extractor<dim,ExtractorType>::template tensor_type<typename ADHelperBase<dim,ADNumberTypeCode,ScalarType>::ad_type>
    ADHelperPointLevelFunctionsBase<dim,ADNumberTypeCode,ScalarType>::get_non_sensitive_variables (const ExtractorType &extractor)
    {
      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          Assert(this->active_tape()!=this->invalid_tape_index,
                 ExcMessage("Invalid tape index"));
        }

      const std::vector<unsigned int> index_set (internal::extract_index_set<dim>(extractor));
      typename internal::Extractor<dim,ExtractorType>::template tensor_type<ad_type> out;

      for (unsigned int i=0; i<index_set.size(); ++i)
        this->get_independent_variable(internal::get_tensor_entry(out,i),
                                       index_set[i]);

      return out;
    }



// -------------------------- ADHelperScalarFunction ----------------------



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    template<typename ExtractorType_Row>
    typename internal::Gradient<dim,typename ADHelperScalarFunction<dim,ADNumberTypeCode,ScalarType>::scalar_type,ExtractorType_Row>::type
    ADHelperScalarFunction<dim,ADNumberTypeCode,ScalarType>::extract_gradient_component(
      const Vector<scalar_type>  &gradient,
      const ExtractorType_Row  &extractor_row) const
    {
      // NOTE: The order of components must be consistently defined throughout this class.
      typename internal::Gradient<dim,scalar_type,ExtractorType_Row>::type out;

      // Get indexsets for the subblock from which we wish to extract the gradient values
      const std::vector<unsigned int> row_index_set (internal::extract_index_set<dim>(extractor_row));
      Assert(out.n_independent_components == row_index_set.size(),
             ExcMessage("Not all tensor components have been extracted!"));
      for (unsigned int r=0; r<row_index_set.size(); ++r)
        internal::set_tensor_entry(out, r,
                                   gradient[row_index_set[r]]);

      return out;
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    template<typename ExtractorType_Row, typename ExtractorType_Col>
    typename internal::Hessian<dim,typename ADHelperScalarFunction<dim,ADNumberTypeCode,ScalarType>::scalar_type,ExtractorType_Row,ExtractorType_Col>::type
    ADHelperScalarFunction<dim,ADNumberTypeCode,ScalarType>::extract_hessian_component(
      const FullMatrix<scalar_type>  &hessian,
      const ExtractorType_Row      &extractor_row,
      const ExtractorType_Col      &extractor_col) const
    {
      typedef internal::Hessian<dim,scalar_type,ExtractorType_Row,ExtractorType_Col> InternalHessian;
      typedef internal::Extractor<dim,ExtractorType_Row> InternalExtractorRow;
      typedef internal::Extractor<dim,ExtractorType_Col> InternalExtractorCol;
      typedef typename InternalHessian::type HessianType;

      // NOTE: The order of components must be consistently defined throughout this class.
      HessianType out;

      // Get indexsets for the subblocks from which we wish to extract the hessian values
      // NOTE: Here we have to do some clever accounting when the one extractor is a symmetric Tensor
      // and the other is not, e.g. <SymmTensor,Vector>. In this scenario the return type is a
      // non-symmetric Tensor<3,dim> but we have to fetch information from a SymmTensor row/column
      // that has too few entries to fill the output tensor. So we must duplicate the relevant
      // entries in the row/column indexset to fetch off-diagonal components that are Otherwise
      // non-existent in a SymmTensor.
      const std::vector<unsigned int> row_index_set (internal::extract_index_set<dim>(extractor_row,false /*ignore_symmetries*/));
      const std::vector<unsigned int> col_index_set (internal::extract_index_set<dim>(extractor_col,false /*ignore_symmetries*/));

      for (unsigned int idx=0; idx<HessianType::n_independent_components; ++idx)
        {
          const TableIndices<HessianType::rank> ti_out = HessianType::unrolled_to_component_indices(idx);
          const unsigned int r = InternalExtractorRow::template local_component<scalar_type>(ti_out, 0);
          const unsigned int c = InternalExtractorCol::template local_component<scalar_type>(ti_out, InternalExtractorRow::rank);

          internal::set_tensor_entry(out, idx,
                                     hessian[row_index_set[r]][col_index_set[c]]);
        }

      return out;
    }



// -------------------------- ADHelperVectorFunction ----------------------



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    template<typename ValueType, typename ExtractorType>
    void
    ADHelperVectorFunction<dim,ADNumberTypeCode,ScalarType>::register_dependent_variable (
      const ValueType     &funcs,
      const ExtractorType &extractor)
    {
      const std::vector<unsigned int> index_set (internal::extract_index_set<dim>(extractor));
      for (unsigned int i=0; i<index_set.size(); ++i)
        {
          Assert(this->touched_dependent_variables[index_set[i]] == false,
                 ExcMessage("Overlapping indices for independent variables."));
          ADHelperBase<dim,ADNumberTypeCode,ScalarType>::register_dependent_variable(
            internal::get_tensor_entry(funcs,i),
            index_set[i]);
        }
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    template<typename ExtractorType_Row>
    typename internal::Gradient<dim,typename ADHelperVectorFunction<dim,ADNumberTypeCode,ScalarType>::scalar_type,ExtractorType_Row>::type
    ADHelperVectorFunction<dim,ADNumberTypeCode,ScalarType>::extract_value_component(
      const Vector<scalar_type>  &values,
      const ExtractorType_Row  &extractor_row) const
    {
      // NOTE: The order of components must be consistently defined throughout this class.
      typename internal::Gradient<dim,scalar_type,ExtractorType_Row>::type out;

      // Get indexsets for the subblock from which we wish to extract the gradient values
      const std::vector<unsigned int> row_index_set (internal::extract_index_set<dim>(extractor_row));
      Assert(out.n_independent_components == row_index_set.size(),
             ExcMessage("Not all tensor components have been extracted!"));
      for (unsigned int r=0; r<row_index_set.size(); ++r)
        internal::set_tensor_entry(out, r,
                                   values[row_index_set[r]]);

      return out;
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    template<typename ExtractorType_Row, typename ExtractorType_Col>
    typename internal::Hessian<dim,typename ADHelperVectorFunction<dim,ADNumberTypeCode,ScalarType>::scalar_type,ExtractorType_Row,ExtractorType_Col>::type
    ADHelperVectorFunction<dim,ADNumberTypeCode,ScalarType>::extract_jacobian_component(
      const FullMatrix<scalar_type>  &jacobian,
      const ExtractorType_Row      &extractor_row,
      const ExtractorType_Col      &extractor_col) const
    {
      typedef internal::Hessian<dim,scalar_type,ExtractorType_Row,ExtractorType_Col> InternalHessian;
      typedef internal::Extractor<dim,ExtractorType_Row> InternalExtractorRow;
      typedef internal::Extractor<dim,ExtractorType_Col> InternalExtractorCol;
      typedef typename InternalHessian::type HessianType;

      // NOTE: The order of components must be consistently defined throughout this class.
      HessianType out;

      // Get indexsets for the subblocks from which we wish to extract the hessian values
      // NOTE: Here we have to do some clever accounting when the one extractor is a symmetric Tensor
      // and the other is not, e.g. <SymmTensor,Vector>. In this scenario the return type is a
      // non-symmetric Tensor<3,dim> but we have to fetch information from a SymmTensor row/column
      // that has too few entries to fill the output tensor. So we must duplicate the relevant
      // entries in the row/column indexset to fetch off-diagonal components that are Otherwise
      // non-existent in a SymmTensor.
      const std::vector<unsigned int> row_index_set (internal::extract_index_set<dim>(extractor_row,false /*ignore_symmetries*/));
      const std::vector<unsigned int> col_index_set (internal::extract_index_set<dim>(extractor_col,false /*ignore_symmetries*/));

      for (unsigned int idx=0; idx<HessianType::n_independent_components; ++idx)
        {
          const TableIndices<HessianType::rank> ti_out = HessianType::unrolled_to_component_indices(idx);
          const unsigned int r = InternalExtractorRow::template local_component<scalar_type>(ti_out, 0);
          const unsigned int c = InternalExtractorCol::template local_component<scalar_type>(ti_out, InternalExtractorRow::rank);

          internal::set_tensor_entry(out, idx,
                                     jacobian[row_index_set[r]][col_index_set[c]]);
        }

      return out;
    }



// -------------------------- ADHelperCellLevelBase ----------------------



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    template<typename VectorType>
    void
    ADHelperCellLevelBase<dim,ADNumberTypeCode,ScalarType>::register_dof_values (
      const VectorType                           &values,
      const std::vector<types::global_dof_index> &local_dof_indices)
    {
      // This is actually the same thing the set_dof_values function,
      // in the sense that we simply populate our array of independent values
      // with a meaningful number. However, in this case we need to double check
      // that we're not registering these variables twice
      Assert(local_dof_indices.size() == this->n_independent_variables(),
             ExcMessage("DoF index vector size does not match number of independent variables"));
#ifdef DEBUG
      for (unsigned int i=0; i<this->n_independent_variables(); ++i)
        {
          Assert(this->touched_independent_variables[i] == false,
                 ExcMessage("Independent variables already registered."));
        }
#endif
      set_dof_values(values, local_dof_indices);
    }



    template<int dim, enum AD::NumberTypes ADNumberTypeCode, typename ScalarType>
    template<typename VectorType>
    void
    ADHelperCellLevelBase<dim,ADNumberTypeCode,ScalarType>::set_dof_values (
      const VectorType                           &values,
      const std::vector<types::global_dof_index> &local_dof_indices)
    {
      if (ADNumberTraits<ad_type>::is_taped == true)
        {
          Assert(this->active_tape()!=this->invalid_tape_index,
                 ExcMessage("Invalid tape index"));
        }
      Assert(local_dof_indices.size() == this->n_independent_variables(),
             ExcMessage("Vector size does not match number of independent variables"));
      for (unsigned int i=0; i<this->n_independent_variables(); ++i)
        ADHelperBase<dim,ADNumberTypeCode,ScalarType>::set_sensitivity_value(values[local_dof_indices[i]], i);
    }



// -------------------------- ADHelperVariationalFormulation ----------------------



// -------------------------- ADHelperResidualLinearisation ----------------------



  } // namespace AD
} // namespace Differentiation


#endif // DOXYGEN


DEAL_II_NAMESPACE_CLOSE

#endif // defined(DEAL_II_WITH_ADOLC) || defined(DEAL_II_WITH_TRILINOS)

#endif // dealii__adolc_helpers_h
