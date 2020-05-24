// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_data_out_dof_data_templates_h
#define dealii_data_out_dof_data_templates_h


#include <deal.II/base/config.h>

#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_out_dof_data.h>

#include <memory>
#include <string>
#include <utility>
#include <vector>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace DataOutImplementation
  {
    template <int dim, int spacedim>
    ParallelDataBase<dim, spacedim>::ParallelDataBase(
      const unsigned int               n_datasets,
      const unsigned int               n_subdivisions,
      const std::vector<unsigned int> &n_postprocessor_outputs,
      const Mapping<dim, spacedim> &   mapping,
      const std::vector<
        std::shared_ptr<dealii::hp::FECollection<dim, spacedim>>>
        &               finite_elements,
      const UpdateFlags update_flags,
      const bool        use_face_values)
      : n_datasets(n_datasets)
      , n_subdivisions(n_subdivisions)
      , postprocessed_values(n_postprocessor_outputs.size())
      , mapping_collection(mapping)
      , finite_elements(finite_elements)
      , update_flags(update_flags)
    {
      unsigned int n_q_points = 0;
      if (use_face_values == false)
        {
          dealii::hp::QCollection<dim> quadrature(
            QIterated<dim>(QTrapez<1>(), n_subdivisions));
          n_q_points = quadrature[0].size();
          x_fe_values.resize(this->finite_elements.size());
          for (unsigned int i = 0; i < this->finite_elements.size(); ++i)
            {
              // check if there is a finite element that is equal to the present
              // one, then we can re-use the FEValues object
              for (unsigned int j = 0; j < i; ++j)
                if (this->finite_elements[i].get() ==
                    this->finite_elements[j].get())
                  {
                    x_fe_values[i] = x_fe_values[j];
                    break;
                  }
              if (x_fe_values[i].get() == nullptr)
                x_fe_values[i] =
                  std::make_shared<dealii::hp::FEValues<dim, spacedim>>(
                    this->mapping_collection,
                    *this->finite_elements[i],
                    quadrature,
                    this->update_flags);
            }
        }
      else
        {
          dealii::hp::QCollection<dim - 1> quadrature(
            QIterated<dim - 1>(QTrapez<1>(), n_subdivisions));
          n_q_points = quadrature[0].size();
          x_fe_face_values.resize(this->finite_elements.size());
          for (unsigned int i = 0; i < this->finite_elements.size(); ++i)
            {
              // check if there is a finite element that is equal to the present
              // one, then we can re-use the FEValues object
              for (unsigned int j = 0; j < i; ++j)
                if (this->finite_elements[i].get() ==
                    this->finite_elements[j].get())
                  {
                    x_fe_face_values[i] = x_fe_face_values[j];
                    break;
                  }
              if (x_fe_face_values[i].get() == nullptr)
                x_fe_face_values[i] =
                  std::make_shared<dealii::hp::FEFaceValues<dim, spacedim>>(
                    this->mapping_collection,
                    *this->finite_elements[i],
                    quadrature,
                    this->update_flags);
            }
        }

      patch_values_scalar.solution_values.resize(n_q_points);
      patch_values_scalar.solution_gradients.resize(n_q_points);
      patch_values_scalar.solution_hessians.resize(n_q_points);
      patch_values_system.solution_values.resize(n_q_points);
      patch_values_system.solution_gradients.resize(n_q_points);
      patch_values_system.solution_hessians.resize(n_q_points);

      for (unsigned int dataset = 0; dataset < n_postprocessor_outputs.size();
           ++dataset)
        if (n_postprocessor_outputs[dataset] != 0)
          postprocessed_values[dataset].resize(
            n_q_points,
            dealii::Vector<double>(n_postprocessor_outputs[dataset]));
    }



    // implement copy constructor to create a thread's own version of
    // x_fe_values
    template <int dim, int spacedim>
    ParallelDataBase<dim, spacedim>::ParallelDataBase(
      const ParallelDataBase<dim, spacedim> &data)
      : n_datasets(data.n_datasets)
      , n_subdivisions(data.n_subdivisions)
      , patch_values_scalar(data.patch_values_scalar)
      , patch_values_system(data.patch_values_system)
      , postprocessed_values(data.postprocessed_values)
      , mapping_collection(data.mapping_collection)
      , finite_elements(data.finite_elements)
      , update_flags(data.update_flags)
    {
      if (data.x_fe_values.empty() == false)
        {
          Assert(data.x_fe_face_values.empty() == true, ExcInternalError());
          dealii::hp::QCollection<dim> quadrature(
            QIterated<dim>(QTrapez<1>(), n_subdivisions));
          x_fe_values.resize(this->finite_elements.size());
          for (unsigned int i = 0; i < this->finite_elements.size(); ++i)
            {
              // check if there is a finite element that is equal to the present
              // one, then we can re-use the FEValues object
              for (unsigned int j = 0; j < i; ++j)
                if (this->finite_elements[i].get() ==
                    this->finite_elements[j].get())
                  {
                    x_fe_values[i] = x_fe_values[j];
                    break;
                  }
              if (x_fe_values[i].get() == nullptr)
                x_fe_values[i] =
                  std::make_shared<dealii::hp::FEValues<dim, spacedim>>(
                    this->mapping_collection,
                    *this->finite_elements[i],
                    quadrature,
                    this->update_flags);
            }
        }
      else
        {
          dealii::hp::QCollection<dim - 1> quadrature(
            QIterated<dim - 1>(QTrapez<1>(), n_subdivisions));
          x_fe_face_values.resize(this->finite_elements.size());
          for (unsigned int i = 0; i < this->finite_elements.size(); ++i)
            {
              // check if there is a finite element that is equal to the present
              // one, then we can re-use the FEValues object
              for (unsigned int j = 0; j < i; ++j)
                if (this->finite_elements[i].get() ==
                    this->finite_elements[j].get())
                  {
                    x_fe_face_values[i] = x_fe_face_values[j];
                    break;
                  }
              if (x_fe_face_values[i].get() == nullptr)
                x_fe_face_values[i] =
                  std::make_shared<dealii::hp::FEFaceValues<dim, spacedim>>(
                    this->mapping_collection,
                    *this->finite_elements[i],
                    quadrature,
                    this->update_flags);
            }
        }
    }



    template <int dim, int spacedim>
    template <typename DoFHandlerType>
    void
    ParallelDataBase<dim, spacedim>::reinit_all_fe_values(
      std::vector<std::shared_ptr<DataEntryBase<DoFHandlerType>>> &dof_data,
      const typename dealii::Triangulation<dim, spacedim>::cell_iterator &cell,
      const unsigned int                                                  face)
    {
      for (unsigned int dataset = 0; dataset < dof_data.size(); ++dataset)
        {
          const bool is_duplicate = std::any_of(
            finite_elements.cbegin(),
            finite_elements.cbegin() + dataset,
            [&](const std::shared_ptr<dealii::hp::FECollection<dim, spacedim>>
                  &fe) { return finite_elements[dataset].get() == fe.get(); });
          if (is_duplicate == false)
            {
              if (cell->active())
                {
                  typename DoFHandlerType::active_cell_iterator dh_cell(
                    &cell->get_triangulation(),
                    cell->level(),
                    cell->index(),
                    dof_data[dataset]->dof_handler);
                  if (x_fe_values.empty())
                    {
                      AssertIndexRange(face, GeometryInfo<dim>::faces_per_cell);
                      x_fe_face_values[dataset]->reinit(dh_cell, face);
                    }
                  else
                    x_fe_values[dataset]->reinit(dh_cell);
                }
              else
                x_fe_values[dataset]->reinit(cell);
            }
        }
      if (dof_data.empty())
        {
          if (x_fe_values.empty())
            {
              AssertIndexRange(face, GeometryInfo<dim>::faces_per_cell);
              x_fe_face_values[0]->reinit(cell, face);
            }
          else
            x_fe_values[0]->reinit(cell);
        }
    }



    template <int dim, int spacedim>
    const FEValuesBase<dim, spacedim> &
    ParallelDataBase<dim, spacedim>::get_present_fe_values(
      const unsigned int dataset) const
    {
      AssertIndexRange(dataset, finite_elements.size());
      if (x_fe_values.empty())
        return x_fe_face_values[dataset]->get_present_fe_values();
      else
        return x_fe_values[dataset]->get_present_fe_values();
    }



    template <int dim, int spacedim>
    void
    ParallelDataBase<dim, spacedim>::resize_system_vectors(
      const unsigned int n_components)
    {
      Assert(patch_values_system.solution_values.size() > 0,
             ExcInternalError());
      AssertDimension(patch_values_system.solution_values.size(),
                      patch_values_system.solution_gradients.size());
      AssertDimension(patch_values_system.solution_values.size(),
                      patch_values_system.solution_hessians.size());
      if (patch_values_system.solution_values[0].size() == n_components)
        return;
      for (unsigned int k = 0; k < patch_values_system.solution_values.size();
           ++k)
        {
          patch_values_system.solution_values[k].reinit(n_components);
          patch_values_system.solution_gradients[k].resize(n_components);
          patch_values_system.solution_hessians[k].resize(n_components);
        }
    }



    /**
     * In a WorkStream context, use this function to append the patch computed
     * by the parallel stage to the array of patches.
     */
    template <int dim, int spacedim>
    void
    append_patch_to_list(
      const DataOutBase::Patch<dim, spacedim> &       patch,
      std::vector<DataOutBase::Patch<dim, spacedim>> &patches)
    {
      patches.push_back(patch);
      patches.back().patch_index = patches.size() - 1;
    }
  } // namespace DataOutImplementation
} // namespace internal

namespace internal
{
  namespace DataOutImplementation
  {
    /**
     * Extract the specified component of a number. This template is used when
     * the given value is assumed to be a real scalar, so asking for the real
     * part is the only valid choice for the second argument.
     */
    template <typename NumberType>
    double
    get_component(const NumberType         value,
                  const ComponentExtractor extract_component)
    {
      (void)extract_component;
      static_assert(
        numbers::NumberTraits<NumberType>::is_complex == false,
        "This function must not be called for complex-valued data types.");
      Assert(extract_component == ComponentExtractor::real_part,
             ExcMessage("You cannot extract anything other than the real "
                        "part from a real number."));
      return value;
    }



    /**
     * Extract the specified component of a number. This template is used when
     * the given value is a complex number
     */
    template <typename NumberType>
    double
    get_component(const std::complex<NumberType> &value,
                  const ComponentExtractor        extract_component)
    {
      switch (extract_component)
        {
          case ComponentExtractor::real_part:
            return value.real();

          case ComponentExtractor::imaginary_part:
            return value.imag();

          default:
            Assert(false, ExcInternalError());
        }

      return numbers::signaling_nan<double>();
    }



    template <int rank, int dim, typename NumberType>
    Tensor<rank, dim>
    get_component(const Tensor<rank, dim, NumberType> &value,
                  const ComponentExtractor             extract_component)
    {
      Assert(extract_component == ComponentExtractor::real_part,
             ExcMessage("You cannot extract anything other than the real "
                        "part from a real number."));

      Tensor<rank, dim, double> t;
      for (unsigned int d = 0; d < dim; ++d)
        t[d] = get_component(value[d], extract_component);

      return t;
    }


    /**
     * Helper class templated on vector type to allow different implementations
     * to extract information from a vector.
     */
    template <typename VectorType>
    struct VectorHelper
    {
      /**
       * extract the @p indices from @p vector and put them into @p values.
       */
      static void
      extract(const VectorType &                          vector,
              const std::vector<types::global_dof_index> &indices,
              const ComponentExtractor                    extract_component,
              std::vector<double> &                       values);
    };



    template <typename VectorType>
    void
    VectorHelper<VectorType>::extract(
      const VectorType &                          vector,
      const std::vector<types::global_dof_index> &indices,
      const ComponentExtractor                    extract_component,
      std::vector<double> &                       values)
    {
      for (unsigned int i = 0; i < values.size(); ++i)
        values[i] = get_component(vector[indices[i]], extract_component);
    }



#ifdef DEAL_II_WITH_TRILINOS
    template <>
    inline void
    VectorHelper<LinearAlgebra::EpetraWrappers::Vector>::extract(
      const LinearAlgebra::EpetraWrappers::Vector & /*vector*/,
      const std::vector<types::global_dof_index> & /*indices*/,
      const ComponentExtractor /*extract_component*/,
      std::vector<double> & /*values*/)
    {
      // TODO: we don't have element access
      Assert(false, ExcNotImplemented());
    }
#endif



#if defined(DEAL_II_TRILINOS_WITH_TPETRA) && defined(DEAL_II_WITH_MPI)
    template <>
    inline void
    VectorHelper<LinearAlgebra::TpetraWrappers::Vector<double>>::extract(
      const LinearAlgebra::TpetraWrappers::Vector<double> & /*vector*/,
      const std::vector<types::global_dof_index> & /*indices*/,
      const ComponentExtractor /*extract_component*/,
      std::vector<double> & /*values*/)
    {
      // TODO: we don't have element access
      Assert(false, ExcNotImplemented());
    }

    template <>
    inline void
    VectorHelper<LinearAlgebra::TpetraWrappers::Vector<float>>::extract(
      const LinearAlgebra::TpetraWrappers::Vector<float> & /*vector*/,
      const std::vector<types::global_dof_index> & /*indices*/,
      const ComponentExtractor /*extract_component*/,
      std::vector<double> & /*values*/)
    {
      // TODO: we don't have element access
      Assert(false, ExcNotImplemented());
    }
#endif



    template <typename DoFHandlerType>
    DataEntryBase<DoFHandlerType>::DataEntryBase(
      const DoFHandlerType *          dofs,
      const std::vector<std::string> &names_in,
      const std::vector<
        DataComponentInterpretation::DataComponentInterpretation>
        &data_component_interpretation)
      : dof_handler(dofs,
                    typeid(
                      dealii::DataOut_DoFData<DoFHandlerType,
                                              DoFHandlerType::dimension,
                                              DoFHandlerType::space_dimension>)
                      .name())
      , names(names_in)
      , data_component_interpretation(data_component_interpretation)
      , postprocessor(nullptr, typeid(*this).name())
      , n_output_variables(names.size())
    {
      Assert(names.size() == data_component_interpretation.size(),
             ExcDimensionMismatch(data_component_interpretation.size(),
                                  names.size()));

      // check that the names use only allowed characters
      for (const auto &name : names)
        {
          (void)name;
          Assert(name.find_first_not_of("abcdefghijklmnopqrstuvwxyz"
                                        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                        "0123456789_<>()") == std::string::npos,
                 Exceptions::DataOutImplementation::ExcInvalidCharacter(
                   name,
                   name.find_first_not_of("abcdefghijklmnopqrstuvwxyz"
                                          "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                          "0123456789_<>()")));
        }
    }



    template <typename DoFHandlerType>
    DataEntryBase<DoFHandlerType>::DataEntryBase(
      const DoFHandlerType *dofs,
      const DataPostprocessor<DoFHandlerType::space_dimension>
        *data_postprocessor)
      : dof_handler(dofs,
                    typeid(
                      dealii::DataOut_DoFData<DoFHandlerType,
                                              DoFHandlerType::dimension,
                                              DoFHandlerType::space_dimension>)
                      .name())
      , names(data_postprocessor->get_names())
      , data_component_interpretation(
          data_postprocessor->get_data_component_interpretation())
      , postprocessor(data_postprocessor, typeid(*this).name())
      , n_output_variables(names.size())
    {
      Assert(data_postprocessor->get_names().size() ==
               data_postprocessor->get_data_component_interpretation().size(),
             ExcDimensionMismatch(
               data_postprocessor->get_names().size(),
               data_postprocessor->get_data_component_interpretation().size()));

      // check that the names use only allowed characters
      for (const auto &name : names)
        {
          (void)name;
          Assert(name.find_first_not_of("abcdefghijklmnopqrstuvwxyz"
                                        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                        "0123456789_<>()") == std::string::npos,
                 Exceptions::DataOutImplementation::ExcInvalidCharacter(
                   name,
                   name.find_first_not_of("abcdefghijklmnopqrstuvwxyz"
                                          "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                          "0123456789_<>()")));
        }
    }



    /**
     * Class that stores a pointer to a vector of type equal to the template
     * argument, and provides the functions to extract data from it.
     *
     * @author Wolfgang Bangerth, 2004
     */
    template <typename DoFHandlerType, typename VectorType>
    class DataEntry : public DataEntryBase<DoFHandlerType>
    {
    public:
      /**
       * Constructor. Give a list of names for the individual components of
       * the vector and their interpretation as scalar or vector data. This
       * constructor assumes that no postprocessor is going to be used.
       */
      DataEntry(const DoFHandlerType *          dofs,
                const VectorType *              data,
                const std::vector<std::string> &names,
                const std::vector<
                  DataComponentInterpretation::DataComponentInterpretation>
                  &data_component_interpretation);

      /**
       * Constructor when a data postprocessor is going to be used. In that
       * case, the names and vector declarations are going to be acquired from
       * the postprocessor.
       */
      DataEntry(const DoFHandlerType *dofs,
                const VectorType *    data,
                const DataPostprocessor<DoFHandlerType::space_dimension>
                  *data_postprocessor);

      /**
       * Assuming that the stored vector is a cell vector, extract the given
       * element from it.
       */
      virtual double
      get_cell_data_value(
        const unsigned int       cell_number,
        const ComponentExtractor extract_component) const override;

      /**
       * Given a FEValuesBase object, extract the values on the present cell
       * from the vector we actually store.
       */
      virtual void
      get_function_values(
        const FEValuesBase<DoFHandlerType::dimension,
                           DoFHandlerType::space_dimension> &fe_patch_values,
        const ComponentExtractor                             extract_component,
        std::vector<double> &patch_values) const override;

      /**
       * Given a FEValuesBase object, extract the values on the present cell
       * from the vector we actually store. This function does the same as the
       * one above but for vector-valued finite elements.
       */
      virtual void
      get_function_values(
        const FEValuesBase<DoFHandlerType::dimension,
                           DoFHandlerType::space_dimension> &fe_patch_values,
        const ComponentExtractor                             extract_component,
        std::vector<dealii::Vector<double>> &patch_values_system)
        const override;

      /**
       * Given a FEValuesBase object, extract the gradients on the present
       * cell from the vector we actually store.
       */
      virtual void
      get_function_gradients(
        const FEValuesBase<DoFHandlerType::dimension,
                           DoFHandlerType::space_dimension> &fe_patch_values,
        const ComponentExtractor                             extract_component,
        std::vector<Tensor<1, DoFHandlerType::space_dimension>>
          &patch_gradients) const override;

      /**
       * Given a FEValuesBase object, extract the gradients on the present
       * cell from the vector we actually store. This function does the same
       * as the one above but for vector-valued finite elements.
       */
      virtual void
      get_function_gradients(
        const FEValuesBase<DoFHandlerType::dimension,
                           DoFHandlerType::space_dimension> &fe_patch_values,
        const ComponentExtractor                             extract_component,
        std::vector<std::vector<Tensor<1, DoFHandlerType::space_dimension>>>
          &patch_gradients_system) const override;

      /**
       * Given a FEValuesBase object, extract the second derivatives on the
       * present cell from the vector we actually store.
       */
      virtual void
      get_function_hessians(
        const FEValuesBase<DoFHandlerType::dimension,
                           DoFHandlerType::space_dimension> &fe_patch_values,
        const ComponentExtractor                             extract_component,
        std::vector<Tensor<2, DoFHandlerType::space_dimension>> &patch_hessians)
        const override;

      /**
       * Given a FEValuesBase object, extract the second derivatives on the
       * present cell from the vector we actually store. This function does
       * the same as the one above but for vector-valued finite elements.
       */
      virtual void
      get_function_hessians(
        const FEValuesBase<DoFHandlerType::dimension,
                           DoFHandlerType::space_dimension> &fe_patch_values,
        const ComponentExtractor                             extract_component,
        std::vector<std::vector<Tensor<2, DoFHandlerType::space_dimension>>>
          &patch_hessians_system) const override;

      /**
       * Return whether the data represented by (a derived class of) this object
       * represents a complex-valued (as opposed to real-valued) information.
       */
      virtual bool
      is_complex_valued() const override;

      /**
       * Clear all references to the vectors.
       */
      virtual void
      clear() override;

      /**
       * Determine an estimate for the memory consumption (in bytes) of this
       * object.
       */
      virtual std::size_t
      memory_consumption() const override;

    private:
      /**
       * Pointer to the data vector. Note that ownership of the vector pointed
       * to remains with the caller of this class.
       */
      const VectorType *vector;
    };



    template <typename DoFHandlerType, typename VectorType>
    DataEntry<DoFHandlerType, VectorType>::DataEntry(
      const DoFHandlerType *          dofs,
      const VectorType *              data,
      const std::vector<std::string> &names,
      const std::vector<
        DataComponentInterpretation::DataComponentInterpretation>
        &data_component_interpretation)
      : DataEntryBase<DoFHandlerType>(dofs,
                                      names,
                                      data_component_interpretation)
      , vector(data)
    {}



    template <typename DoFHandlerType, typename VectorType>
    DataEntry<DoFHandlerType, VectorType>::DataEntry(
      const DoFHandlerType *dofs,
      const VectorType *    data,
      const DataPostprocessor<DoFHandlerType::space_dimension>
        *data_postprocessor)
      : DataEntryBase<DoFHandlerType>(dofs, data_postprocessor)
      , vector(data)
    {}



    template <typename DoFHandlerType, typename VectorType>
    double
    DataEntry<DoFHandlerType, VectorType>::get_cell_data_value(
      const unsigned int       cell_number,
      const ComponentExtractor extract_component) const
    {
      return get_component(
        internal::ElementAccess<VectorType>::get(*vector, cell_number),
        extract_component);
    }



    template <typename DoFHandlerType, typename VectorType>
    void
    DataEntry<DoFHandlerType, VectorType>::get_function_values(
      const FEValuesBase<DoFHandlerType::dimension,
                         DoFHandlerType::space_dimension> &fe_patch_values,
      const ComponentExtractor                             extract_component,
      std::vector<dealii::Vector<double>> &patch_values_system) const
    {
      if (typeid(typename VectorType::value_type) == typeid(double))
        {
          Assert(extract_component == ComponentExtractor::real_part,
                 ExcMessage("You cannot extract anything other than the real "
                            "part from a real number."));

          fe_patch_values.get_function_values(
            *vector,
            // reinterpret output argument type; because of the 'if' statement
            // above, this is the identity cast whenever the code is executed,
            // but the cast is necessary to allow compilation even if we don't
            // get here
            reinterpret_cast<
              std::vector<dealii::Vector<typename VectorType::value_type>> &>(
              patch_values_system));
        }
      else
        {
          // OK, so we know that the data type stored by the user-provided
          // vector is not simply a double. In that case, we need to ask the
          // FEValuesBase object to extract function values for us from the
          // evaluation points in the provided data type, and then copy them by
          // hand into the output location.
          const unsigned int n_components =
            fe_patch_values.get_fe().n_components();
          const unsigned int n_eval_points =
            fe_patch_values.n_quadrature_points;

          std::vector<dealii::Vector<typename VectorType::value_type>> tmp(
            n_eval_points);
          for (unsigned int i = 0; i < n_eval_points; i++)
            tmp[i].reinit(n_components);

          fe_patch_values.get_function_values(*vector, tmp);

          AssertDimension(patch_values_system.size(), n_eval_points);
          for (unsigned int i = 0; i < n_eval_points; i++)
            {
              AssertDimension(patch_values_system[i].size(), n_components);

              for (unsigned int j = 0; j < n_components; ++j)
                patch_values_system[i](j) =
                  get_component(tmp[i](j), extract_component);
            }
        }
    }



    template <typename DoFHandlerType, typename VectorType>
    void
    DataEntry<DoFHandlerType, VectorType>::get_function_values(
      const FEValuesBase<DoFHandlerType::dimension,
                         DoFHandlerType::space_dimension> &fe_patch_values,
      const ComponentExtractor                             extract_component,
      std::vector<double> &                                patch_values) const
    {
      if (typeid(typename VectorType::value_type) == typeid(double))
        {
          Assert(extract_component == ComponentExtractor::real_part,
                 ExcMessage("You cannot extract anything other than the real "
                            "part from a real number."));

          fe_patch_values.get_function_values(
            *vector,
            // reinterpret output argument type; because of the 'if' statement
            // above, this is the identity cast whenever the code is executed,
            // but the cast is necessary to allow compilation even if we don't
            // get here
            reinterpret_cast<std::vector<typename VectorType::value_type> &>(
              patch_values));
        }
      else
        {
          std::vector<typename VectorType::value_type> tmp(patch_values.size());

          fe_patch_values.get_function_values(*vector, tmp);

          for (unsigned int i = 0; i < tmp.size(); i++)
            patch_values[i] = get_component(tmp[i], extract_component);
        }
    }



    template <typename DoFHandlerType, typename VectorType>
    void
    DataEntry<DoFHandlerType, VectorType>::get_function_gradients(
      const FEValuesBase<DoFHandlerType::dimension,
                         DoFHandlerType::space_dimension> &fe_patch_values,
      const ComponentExtractor                             extract_component,
      std::vector<std::vector<Tensor<1, DoFHandlerType::space_dimension>>>
        &patch_gradients_system) const
    {
      if (typeid(typename VectorType::value_type) == typeid(double))
        {
          Assert(extract_component == ComponentExtractor::real_part,
                 ExcMessage("You cannot extract anything other than the real "
                            "part from a real number."));

          fe_patch_values.get_function_gradients(
            *vector,
            // reinterpret output argument type; because of the 'if' statement
            // above, this is the identity cast whenever the code is executed,
            // but the cast is necessary to allow compilation even if we don't
            // get here
            reinterpret_cast<std::vector<
              std::vector<Tensor<1,
                                 DoFHandlerType::space_dimension,
                                 typename VectorType::value_type>>> &>(
              patch_gradients_system));
        }
      else
        {
          // OK, so we know that the data type stored by the user-provided
          // vector is not simply a double. In that case, we need to ask the
          // FEValuesBase object to extract function values for us from the
          // evaluation points in the provided data type, and then copy them by
          // hand into the output location.
          const unsigned int n_components =
            fe_patch_values.get_fe().n_components();
          const unsigned int n_eval_points =
            fe_patch_values.n_quadrature_points;

          std::vector<std::vector<Tensor<1,
                                         DoFHandlerType::space_dimension,
                                         typename VectorType::value_type>>>
            tmp(n_eval_points);
          for (unsigned int i = 0; i < n_eval_points; i++)
            tmp[i].resize(n_components);

          fe_patch_values.get_function_gradients(*vector, tmp);

          AssertDimension(patch_gradients_system.size(), n_eval_points);
          for (unsigned int i = 0; i < n_eval_points; i++)
            {
              AssertDimension(patch_gradients_system[i].size(), n_components);

              for (unsigned int j = 0; j < n_components; j++)
                patch_gradients_system[i][j] =
                  get_component(tmp[i][j], extract_component);
            }
        }
    }



    template <typename DoFHandlerType, typename VectorType>
    void
    DataEntry<DoFHandlerType, VectorType>::get_function_gradients(
      const FEValuesBase<DoFHandlerType::dimension,
                         DoFHandlerType::space_dimension> &fe_patch_values,
      const ComponentExtractor                             extract_component,
      std::vector<Tensor<1, DoFHandlerType::space_dimension>> &patch_gradients)
      const
    {
      if (typeid(typename VectorType::value_type) == typeid(double))
        {
          Assert(extract_component == ComponentExtractor::real_part,
                 ExcMessage("You cannot extract anything other than the real "
                            "part from a real number."));

          fe_patch_values.get_function_gradients(
            *vector,
            // reinterpret output argument type; because of the 'if' statement
            // above, this is the identity cast whenever the code is executed,
            // but the cast is necessary to allow compilation even if we don't
            // get here
            reinterpret_cast<
              std::vector<Tensor<1,
                                 DoFHandlerType::space_dimension,
                                 typename VectorType::value_type>> &>(
              patch_gradients));
        }
      else
        {
          std::vector<Tensor<1,
                             DoFHandlerType::space_dimension,
                             typename VectorType::value_type>>
            tmp;
          tmp.resize(patch_gradients.size());

          fe_patch_values.get_function_gradients(*vector, tmp);

          for (unsigned int i = 0; i < tmp.size(); i++)
            patch_gradients[i] = get_component(tmp[i], extract_component);
        }
    }



    template <typename DoFHandlerType, typename VectorType>
    void
    DataEntry<DoFHandlerType, VectorType>::get_function_hessians(
      const FEValuesBase<DoFHandlerType::dimension,
                         DoFHandlerType::space_dimension> &fe_patch_values,
      const ComponentExtractor                             extract_component,
      std::vector<std::vector<Tensor<2, DoFHandlerType::space_dimension>>>
        &patch_hessians_system) const
    {
      if (typeid(typename VectorType::value_type) == typeid(double))
        {
          Assert(extract_component == ComponentExtractor::real_part,
                 ExcMessage("You cannot extract anything other than the real "
                            "part from a real number."));

          fe_patch_values.get_function_hessians(
            *vector,
            // reinterpret output argument type; because of the 'if' statement
            // above, this is the identity cast whenever the code is executed,
            // but the cast is necessary to allow compilation even if we don't
            // get here
            reinterpret_cast<std::vector<
              std::vector<Tensor<2,
                                 DoFHandlerType::space_dimension,
                                 typename VectorType::value_type>>> &>(
              patch_hessians_system));
        }
      else
        {
          // OK, so we know that the data type stored by the user-provided
          // vector is not simply a double. In that case, we need to ask the
          // FEValuesBase object to extract function values for us from the
          // evaluation points in the provided data type, and then copy them by
          // hand into the output location.
          const unsigned int n_components =
            fe_patch_values.get_fe().n_components();
          const unsigned int n_eval_points =
            fe_patch_values.n_quadrature_points;

          std::vector<std::vector<Tensor<2,
                                         DoFHandlerType::space_dimension,
                                         typename VectorType::value_type>>>
            tmp(n_eval_points);
          for (unsigned int i = 0; i < n_eval_points; i++)
            tmp[i].resize(n_components);

          fe_patch_values.get_function_hessians(*vector, tmp);

          AssertDimension(patch_hessians_system.size(), n_eval_points);
          for (unsigned int i = 0; i < n_eval_points; i++)
            {
              AssertDimension(patch_hessians_system[i].size(), n_components);

              for (unsigned int j = 0; j < n_components; j++)
                patch_hessians_system[i][j] =
                  get_component(tmp[i][j], extract_component);
            }
        }
    }



    template <typename DoFHandlerType, typename VectorType>
    void
    DataEntry<DoFHandlerType, VectorType>::get_function_hessians(
      const FEValuesBase<DoFHandlerType::dimension,
                         DoFHandlerType::space_dimension> &fe_patch_values,
      const ComponentExtractor                             extract_component,
      std::vector<Tensor<2, DoFHandlerType::space_dimension>> &patch_hessians)
      const
    {
      if (typeid(typename VectorType::value_type) == typeid(double))
        {
          Assert(extract_component == ComponentExtractor::real_part,
                 ExcMessage("You cannot extract anything other than the real "
                            "part from a real number."));

          fe_patch_values.get_function_hessians(
            *vector,
            // reinterpret output argument type; because of the 'if' statement
            // above, this is the identity cast whenever the code is executed,
            // but the cast is necessary to allow compilation even if we don't
            // get here
            reinterpret_cast<
              std::vector<Tensor<2,
                                 DoFHandlerType::space_dimension,
                                 typename VectorType::value_type>> &>(
              patch_hessians));
        }
      else
        {
          std::vector<Tensor<2,
                             DoFHandlerType::space_dimension,
                             typename VectorType::value_type>>
            tmp(patch_hessians.size());

          fe_patch_values.get_function_hessians(*vector, tmp);

          for (unsigned int i = 0; i < tmp.size(); i++)
            patch_hessians[i] = get_component(tmp[i], extract_component);
        }
    }



    template <typename DoFHandlerType, typename VectorType>
    bool
    DataEntry<DoFHandlerType, VectorType>::is_complex_valued() const
    {
      return numbers::NumberTraits<typename VectorType::value_type>::is_complex;
    }



    template <typename DoFHandlerType, typename VectorType>
    std::size_t
    DataEntry<DoFHandlerType, VectorType>::memory_consumption() const
    {
      return (sizeof(vector) +
              MemoryConsumption::memory_consumption(this->names));
    }



    template <typename DoFHandlerType, typename VectorType>
    void
    DataEntry<DoFHandlerType, VectorType>::clear()
    {
      vector            = nullptr;
      this->dof_handler = nullptr;
    }



    /**
     * Like DataEntry, but used to look up data from multigrid computations.
     * Data will use level-DoF indices to look up in a
     * MGLevelObject<VectorType> given on the specific level instead of
     * interpolating data to coarser cells.
     */
    template <typename DoFHandlerType, typename VectorType>
    class MGDataEntry : public DataEntryBase<DoFHandlerType>
    {
    public:
      MGDataEntry(const DoFHandlerType *           dofs,
                  const MGLevelObject<VectorType> *vectors,
                  const std::vector<std::string> & names,
                  const std::vector<
                    DataComponentInterpretation::DataComponentInterpretation>
                    &data_component_interpretation)
        : DataEntryBase<DoFHandlerType>(dofs,
                                        names,
                                        data_component_interpretation)
        , vectors(vectors)
      {}

      virtual double
      get_cell_data_value(
        const unsigned int       cell_number,
        const ComponentExtractor extract_component) const override;

      virtual void
      get_function_values(
        const FEValuesBase<DoFHandlerType::dimension,
                           DoFHandlerType::space_dimension> &fe_patch_values,
        const ComponentExtractor                             extract_component,
        std::vector<double> &patch_values) const override;

      /**
       * Given a FEValuesBase object, extract the values on the present cell
       * from the vector we actually store. This function does the same as the
       * one above but for vector-valued finite elements.
       */
      virtual void
      get_function_values(
        const FEValuesBase<DoFHandlerType::dimension,
                           DoFHandlerType::space_dimension> &fe_patch_values,
        const ComponentExtractor                             extract_component,
        std::vector<dealii::Vector<double>> &patch_values_system)
        const override;

      /**
       * Given a FEValuesBase object, extract the gradients on the present
       * cell from the vector we actually store.
       */
      virtual void
      get_function_gradients(
        const FEValuesBase<DoFHandlerType::dimension,
                           DoFHandlerType::space_dimension>
          & /*fe_patch_values*/,
        const ComponentExtractor /*extract_component*/,
        std::vector<Tensor<1, DoFHandlerType::space_dimension>>
          & /*patch_gradients*/) const override
      {
        Assert(false, ExcNotImplemented());
      }

      /**
       * Given a FEValuesBase object, extract the gradients on the present
       * cell from the vector we actually store. This function does the same
       * as the one above but for vector-valued finite elements.
       */
      virtual void
      get_function_gradients(
        const FEValuesBase<DoFHandlerType::dimension,
                           DoFHandlerType::space_dimension>
          & /*fe_patch_values*/,
        const ComponentExtractor /*extract_component*/,
        std::vector<std::vector<Tensor<1, DoFHandlerType::space_dimension>>>
          & /*patch_gradients_system*/) const override
      {
        Assert(false, ExcNotImplemented());
      }


      /**
       * Given a FEValuesBase object, extract the second derivatives on the
       * present cell from the vector we actually store.
       */
      virtual void
      get_function_hessians(
        const FEValuesBase<DoFHandlerType::dimension,
                           DoFHandlerType::space_dimension>
          & /*fe_patch_values*/,
        const ComponentExtractor /*extract_component*/,
        std::vector<Tensor<2, DoFHandlerType::space_dimension>>
          & /*patch_hessians*/) const override
      {
        Assert(false, ExcNotImplemented());
      }

      /**
       * Given a FEValuesBase object, extract the second derivatives on the
       * present cell from the vector we actually store. This function does
       * the same as the one above but for vector-valued finite elements.
       */
      virtual void
      get_function_hessians(
        const FEValuesBase<DoFHandlerType::dimension,
                           DoFHandlerType::space_dimension>
          & /*fe_patch_values*/,
        const ComponentExtractor /*extract_component*/,
        std::vector<std::vector<Tensor<2, DoFHandlerType::space_dimension>>>
          & /*patch_hessians_system*/) const override
      {
        Assert(false, ExcNotImplemented());
      }

      /**
       * Return whether the data represented by (a derived class of) this object
       * represents a complex-valued (as opposed to real-valued) information.
       */
      virtual bool
      is_complex_valued() const override
      {
        Assert(
          numbers::NumberTraits<typename VectorType::value_type>::is_complex ==
            false,
          ExcNotImplemented());
        return numbers::NumberTraits<
          typename VectorType::value_type>::is_complex;
      }

      /**
       * Clear all references to the vectors.
       */
      virtual void
      clear() override
      {
        vectors = nullptr;
      }

      /**
       * Determine an estimate for the memory consumption (in bytes) of this
       * object.
       */
      virtual std::size_t
      memory_consumption() const override
      {
        return sizeof(vectors);
      }

    private:
      const MGLevelObject<VectorType> *vectors;
    };



    template <typename DoFHandlerType, typename VectorType>
    double
    MGDataEntry<DoFHandlerType, VectorType>::get_cell_data_value(
      const unsigned int       cell_number,
      const ComponentExtractor extract_component) const
    {
      Assert(false, ExcNotImplemented());

      (void)cell_number;
      (void)extract_component;
      return 0.0;
    }



    template <typename DoFHandlerType, typename VectorType>
    void
    MGDataEntry<DoFHandlerType, VectorType>::get_function_values(
      const FEValuesBase<DoFHandlerType::dimension,
                         DoFHandlerType::space_dimension> &fe_patch_values,
      const ComponentExtractor                             extract_component,
      std::vector<double> &                                patch_values) const
    {
      Assert(extract_component == ComponentExtractor::real_part,
             ExcNotImplemented());

      const typename DoFHandlerType::level_cell_iterator dof_cell(
        &fe_patch_values.get_cell()->get_triangulation(),
        fe_patch_values.get_cell()->level(),
        fe_patch_values.get_cell()->index(),
        this->dof_handler);

      const VectorType *vector = &((*vectors)[dof_cell->level()]);

      const unsigned int dofs_per_cell =
        this->dof_handler->get_fe(0).dofs_per_cell;

      std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
      dof_cell->get_mg_dof_indices(dof_indices);
      std::vector<double> values(dofs_per_cell);
      VectorHelper<VectorType>::extract(*vector,
                                        dof_indices,
                                        extract_component,
                                        values);
      std::fill(patch_values.begin(), patch_values.end(), 0.0);

      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        for (unsigned int q_point = 0; q_point < patch_values.size(); ++q_point)
          patch_values[q_point] +=
            values[i] * fe_patch_values.shape_value(i, q_point);
    }



    template <typename DoFHandlerType, typename VectorType>
    void
    MGDataEntry<DoFHandlerType, VectorType>::get_function_values(
      const FEValuesBase<DoFHandlerType::dimension,
                         DoFHandlerType::space_dimension> &fe_patch_values,
      const ComponentExtractor                             extract_component,
      std::vector<dealii::Vector<double>> &patch_values_system) const
    {
      Assert(extract_component == ComponentExtractor::real_part,
             ExcNotImplemented());

      typename DoFHandlerType::level_cell_iterator dof_cell(
        &fe_patch_values.get_cell()->get_triangulation(),
        fe_patch_values.get_cell()->level(),
        fe_patch_values.get_cell()->index(),
        this->dof_handler);

      const VectorType *vector = &((*vectors)[dof_cell->level()]);

      const unsigned int dofs_per_cell =
        this->dof_handler->get_fe(0).dofs_per_cell;

      std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
      dof_cell->get_mg_dof_indices(dof_indices);
      std::vector<double> values(dofs_per_cell);
      VectorHelper<VectorType>::extract(*vector,
                                        dof_indices,
                                        extract_component,
                                        values);

      const unsigned int n_components = fe_patch_values.get_fe().n_components();
      const unsigned int n_eval_points = fe_patch_values.n_quadrature_points;

      AssertDimension(patch_values_system.size(), n_eval_points);
      for (unsigned int q = 0; q < n_eval_points; q++)
        {
          AssertDimension(patch_values_system[q].size(), n_components);
          patch_values_system[q] = 0.0;

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            for (unsigned int c = 0; c < n_components; ++c)
              patch_values_system[q](c) +=
                values[i] * fe_patch_values.shape_value_component(i, q, c);
        }
    }

  } // namespace DataOutImplementation
} // namespace internal



template <typename DoFHandlerType, int patch_dim, int patch_space_dim>
DataOut_DoFData<DoFHandlerType, patch_dim, patch_space_dim>::DataOut_DoFData()
  : triangulation(nullptr, typeid(*this).name())
  , dofs(nullptr, typeid(*this).name())
{}



template <typename DoFHandlerType, int patch_dim, int patch_space_dim>
DataOut_DoFData<DoFHandlerType, patch_dim, patch_space_dim>::~DataOut_DoFData()
{
  // virtual functions called in constructors and destructors never use the
  // override in a derived class for clarity be explicit on which function is
  // called
  DataOut_DoFData<DoFHandlerType, patch_dim, patch_space_dim>::clear();
}



template <typename DoFHandlerType, int patch_dim, int patch_space_dim>
void
DataOut_DoFData<DoFHandlerType, patch_dim, patch_space_dim>::attach_dof_handler(
  const DoFHandlerType &d)
{
  Assert(dof_data.size() == 0,
         Exceptions::DataOutImplementation::ExcOldDataStillPresent());
  Assert(cell_data.size() == 0,
         Exceptions::DataOutImplementation::ExcOldDataStillPresent());

  triangulation =
    SmartPointer<const Triangulation<DoFHandlerType::dimension,
                                     DoFHandlerType::space_dimension>>(
      &d.get_triangulation(), typeid(*this).name());
  dofs = SmartPointer<const DoFHandlerType>(&d, typeid(*this).name());
}



template <typename DoFHandlerType, int patch_dim, int patch_space_dim>
void
DataOut_DoFData<DoFHandlerType, patch_dim, patch_space_dim>::
  attach_triangulation(
    const Triangulation<DoFHandlerType::dimension,
                        DoFHandlerType::space_dimension> &tria)
{
  Assert(dof_data.size() == 0,
         Exceptions::DataOutImplementation::ExcOldDataStillPresent());
  Assert(cell_data.size() == 0,
         Exceptions::DataOutImplementation::ExcOldDataStillPresent());

  triangulation =
    SmartPointer<const Triangulation<DoFHandlerType::dimension,
                                     DoFHandlerType::space_dimension>>(
      &tria, typeid(*this).name());
}



template <typename DoFHandlerType, int patch_dim, int patch_space_dim>
template <typename VectorType>
void
DataOut_DoFData<DoFHandlerType, patch_dim, patch_space_dim>::add_data_vector(
  const DoFHandlerType &                                    dof_handler,
  const VectorType &                                        vec,
  const DataPostprocessor<DoFHandlerType::space_dimension> &data_postprocessor)
{
  // this is a specialized version of the other function where we have a
  // postprocessor. if we do, we know that we have type_dof_data, which makes
  // things a bit simpler, we also don't need to deal with some of the other
  // stuff and use a different constructor of DataEntry
  if (triangulation != nullptr)
    {
      Assert(&dof_handler.get_triangulation() == triangulation,
             ExcMessage("The triangulation attached to the DoFHandler does not "
                        "match with the one set previously"));
    }
  else
    {
      triangulation =
        SmartPointer<const Triangulation<DoFHandlerType::dimension,
                                         DoFHandlerType::space_dimension>>(
          &dof_handler.get_triangulation(), typeid(*this).name());
    }

  Assert(vec.size() == dof_handler.n_dofs(),
         Exceptions::DataOutImplementation::ExcInvalidVectorSize(
           vec.size(),
           dof_handler.n_dofs(),
           dof_handler.get_triangulation().n_active_cells()));


  auto new_entry = std::make_unique<
    internal::DataOutImplementation::DataEntry<DoFHandlerType, VectorType>>(
    &dof_handler, &vec, &data_postprocessor);
  dof_data.emplace_back(std::move(new_entry));
}



template <typename DoFHandlerType, int patch_dim, int patch_space_dim>
template <typename VectorType>
void
DataOut_DoFData<DoFHandlerType, patch_dim, patch_space_dim>::
  add_data_vector_internal(
    const DoFHandlerType *          dof_handler,
    const VectorType &              data_vector,
    const std::vector<std::string> &names,
    const DataVectorType            type,
    const std::vector<DataComponentInterpretation::DataComponentInterpretation>
      &        data_component_interpretation_,
    const bool deduce_output_names)
{
  // Check available mesh information:
  if (triangulation == nullptr)
    {
      Assert(dof_handler != nullptr, ExcInternalError());
      triangulation =
        SmartPointer<const Triangulation<DoFHandlerType::dimension,
                                         DoFHandlerType::space_dimension>>(
          &dof_handler->get_triangulation(), typeid(*this).name());
    }

  if (dof_handler != nullptr)
    {
      Assert(&dof_handler->get_triangulation() == triangulation,
             ExcMessage("The triangulation attached to the DoFHandler does not "
                        "match with the one set previously"));
    }

  // Figure out the data type:
  DataVectorType actual_type = type;
  if (type == type_automatic)
    {
      Assert(
        (dof_handler == nullptr) ||
          (triangulation->n_active_cells() != dof_handler->n_dofs()),
        ExcMessage(
          "Unable to determine the type of vector automatically because the number of DoFs "
          "is equal to the number of cells. Please specify DataVectorType."));

      if (data_vector.size() == triangulation->n_active_cells())
        actual_type = type_cell_data;
      else
        actual_type = type_dof_data;
    }
  Assert(actual_type == type_cell_data || actual_type == type_dof_data,
         ExcInternalError());

  // If necessary, append '_1', '_2', etc. to component names:
  std::vector<std::string> deduced_names;
  if (deduce_output_names && actual_type == type_dof_data)
    {
      Assert(names.size() == 1, ExcInternalError());
      Assert(dof_handler != nullptr, ExcInternalError());
      Assert(dof_handler->n_dofs() > 0,
             ExcMessage("The DoF handler attached to the current output vector "
                        "does not have any degrees of freedom, so it is not "
                        "possible to output DoF data in this context."));
      const std::string  name         = names[0];
      const unsigned int n_components = dof_handler->get_fe(0).n_components();
      deduced_names.resize(n_components);
      if (n_components > 1)
        {
          for (unsigned int i = 0; i < n_components; ++i)
            {
              deduced_names[i] = name + '_' + std::to_string(i);
            }
        }
      else
        {
          deduced_names[0] = name;
        }
    }
  else
    {
      deduced_names = names;
    }

  // Check that things have the right sizes for the data type:
  switch (actual_type)
    {
      case type_cell_data:
        Assert(data_vector.size() == triangulation->n_active_cells(),
               ExcDimensionMismatch(data_vector.size(),
                                    triangulation->n_active_cells()));
        Assert(deduced_names.size() == 1,
               Exceptions::DataOutImplementation::ExcInvalidNumberOfNames(
                 deduced_names.size(), 1));
        break;
      case type_dof_data:
        Assert(dof_handler != nullptr,
               Exceptions::DataOutImplementation::ExcNoDoFHandlerSelected());
        Assert(data_vector.size() == dof_handler->n_dofs(),
               Exceptions::DataOutImplementation::ExcInvalidVectorSize(
                 data_vector.size(),
                 dof_handler->n_dofs(),
                 triangulation->n_active_cells()));
        Assert(deduced_names.size() == dof_handler->get_fe(0).n_components(),
               Exceptions::DataOutImplementation::ExcInvalidNumberOfNames(
                 deduced_names.size(), dof_handler->get_fe(0).n_components()));
        break;
      default:
        Assert(false, ExcInternalError());
    }

  const auto &data_component_interpretation =
    (data_component_interpretation_.size() != 0 ?
       data_component_interpretation_ :
       std::vector<DataComponentInterpretation::DataComponentInterpretation>(
         deduced_names.size(),
         DataComponentInterpretation::component_is_scalar));

  // finally, add the data vector:
  auto new_entry = std::make_unique<
    internal::DataOutImplementation::DataEntry<DoFHandlerType, VectorType>>(
    dof_handler, &data_vector, deduced_names, data_component_interpretation);

  if (actual_type == type_dof_data)
    dof_data.emplace_back(std::move(new_entry));
  else
    cell_data.emplace_back(std::move(new_entry));
}



template <typename DoFHandlerType, int patch_dim, int patch_space_dim>
template <class VectorType>
void
DataOut_DoFData<DoFHandlerType, patch_dim, patch_space_dim>::add_mg_data_vector(
  const DoFHandlerType &           dof_handler,
  const MGLevelObject<VectorType> &data,
  const std::string &              name)
{
  // forward the call to the vector version:
  std::vector<std::string> names(1, name);
  add_mg_data_vector(dof_handler, data, names);
}



template <typename DoFHandlerType, int patch_dim, int patch_space_dim>
template <class VectorType>
void
DataOut_DoFData<DoFHandlerType, patch_dim, patch_space_dim>::add_mg_data_vector(
  const DoFHandlerType &           dof_handler,
  const MGLevelObject<VectorType> &data,
  const std::vector<std::string> & names,
  const std::vector<DataComponentInterpretation::DataComponentInterpretation>
    &data_component_interpretation_)
{
  if (triangulation == nullptr)
    triangulation =
      SmartPointer<const Triangulation<DoFHandlerType::dimension,
                                       DoFHandlerType::space_dimension>>(
        &dof_handler.get_triangulation(), typeid(*this).name());

  Assert(&dof_handler.get_triangulation() == triangulation,
         ExcMessage("The triangulation attached to the DoFHandler does not "
                    "match with the one set previously"));

  const unsigned int       n_components  = dof_handler.get_fe(0).n_components();
  std::vector<std::string> deduced_names = names;

  if (names.size() == 1 && n_components > 1)
    {
      deduced_names.resize(n_components);
      for (unsigned int i = 0; i < n_components; ++i)
        {
          deduced_names[i] = names[0] + '_' + std::to_string(i);
        }
    }

  Assert(deduced_names.size() == n_components,
         ExcMessage("Invalid number of names given."));

  const std::vector<DataComponentInterpretation::DataComponentInterpretation>
    &data_component_interpretation =
      (data_component_interpretation_.size() != 0 ?
         data_component_interpretation_ :
         std::vector<DataComponentInterpretation::DataComponentInterpretation>(
           n_components, DataComponentInterpretation::component_is_scalar));

  Assert(data_component_interpretation.size() == n_components,
         ExcMessage(
           "Invalid number of entries in data_component_interpretation."));

  auto new_entry = std::make_unique<
    internal::DataOutImplementation::MGDataEntry<DoFHandlerType, VectorType>>(
    &dof_handler, &data, deduced_names, data_component_interpretation);
  dof_data.emplace_back(std::move(new_entry));
}



template <typename DoFHandlerType, int patch_dim, int patch_space_dim>
void
DataOut_DoFData<DoFHandlerType, patch_dim, patch_space_dim>::
  clear_data_vectors()
{
  dof_data.erase(dof_data.begin(), dof_data.end());
  cell_data.erase(cell_data.begin(), cell_data.end());

  // delete patches
  std::vector<Patch> dummy;
  patches.swap(dummy);
}



template <typename DoFHandlerType, int patch_dim, int patch_space_dim>
void
DataOut_DoFData<DoFHandlerType, patch_dim, patch_space_dim>::
  clear_input_data_references()
{
  for (unsigned int i = 0; i < dof_data.size(); ++i)
    dof_data[i]->clear();

  for (unsigned int i = 0; i < cell_data.size(); ++i)
    cell_data[i]->clear();

  if (dofs != nullptr)
    dofs = nullptr;
}



template <typename DoFHandlerType, int patch_dim, int patch_space_dim>
void
DataOut_DoFData<DoFHandlerType, patch_dim, patch_space_dim>::clear()
{
  dof_data.erase(dof_data.begin(), dof_data.end());
  cell_data.erase(cell_data.begin(), cell_data.end());

  if (dofs != nullptr)
    dofs = nullptr;

  // delete patches
  std::vector<Patch> dummy;
  patches.swap(dummy);
}



template <typename DoFHandlerType, int patch_dim, int patch_space_dim>
std::vector<std::string>
DataOut_DoFData<DoFHandlerType, patch_dim, patch_space_dim>::get_dataset_names()
  const
{
  std::vector<std::string> names;

  // Loop over all DoF-data datasets and push the names. If the
  // vector underlying a data set is complex-valued, then
  // expand it into its real and imaginary part. Note, however,
  // that what comes back from a postprocessor is *always*
  // real-valued, regardless of what goes in, so we don't
  // have this to do this name expansion for data sets that
  // have a postprocessor.
  //
  // The situation is made complicated when considering vector- and
  // tensor-valued component sets. This is because if, for example, we have a
  // complex-valued vector, we don't want to output Re(u_x), then Im(u_x), then
  // Re(u_y), etc. That's because if we did this, then visualization programs
  // will not easily be able to patch together the 1st, 3rd, 5th components into
  // the vector representing the real part of a vector field, and similarly for
  // the 2nd, 4th, 6th component for the imaginary part of the vector field.
  // Rather, we need to put all real components of the same vector field into
  // consecutive components.
  //
  // This sort of logic is also explained in some detail in
  //   DataOut::build_one_patch().
  for (const auto &input_data : dof_data)
    if (input_data->is_complex_valued() == false ||
        (input_data->postprocessor != nullptr))
      {
        for (const auto &name : input_data->names)
          names.push_back(name);
      }
    else
      {
        // OK, so we have a complex-valued vector. We then need to go through
        // all components and order them appropriately
        for (unsigned int i = 0; i < input_data->names.size();
             /* increment of i happens below */)
          {
            switch (input_data->data_component_interpretation[i])
              {
                case DataComponentInterpretation::component_is_scalar:
                  {
                    // It's a scalar. Just output real and imaginary parts one
                    // after the other:
                    names.push_back(input_data->names[i] + "_re");
                    names.push_back(input_data->names[i] + "_im");

                    // Move forward by one component
                    ++i;

                    break;
                  }

                case DataComponentInterpretation::component_is_part_of_vector:
                  {
                    // It's a vector. First output all real parts, then all
                    // imaginary parts:
                    const unsigned int size = patch_space_dim;
                    for (unsigned int vec_comp = 0; vec_comp < size; ++vec_comp)
                      names.push_back(input_data->names[i + vec_comp] + "_re");
                    for (unsigned int vec_comp = 0; vec_comp < size; ++vec_comp)
                      names.push_back(input_data->names[i + vec_comp] + "_im");

                    // Move forward by dim components
                    i += size;

                    break;
                  }

                case DataComponentInterpretation::component_is_part_of_tensor:
                  {
                    // It's a tensor. First output all real parts, then all
                    // imaginary parts:
                    const unsigned int size = patch_space_dim * patch_space_dim;
                    for (unsigned int tensor_comp = 0; tensor_comp < size;
                         ++tensor_comp)
                      names.push_back(input_data->names[i + tensor_comp] +
                                      "_re");
                    for (unsigned int tensor_comp = 0; tensor_comp < size;
                         ++tensor_comp)
                      names.push_back(input_data->names[i + tensor_comp] +
                                      "_im");

                    // Move forward by dim*dim components
                    i += size;

                    break;
                  }

                default:
                  Assert(false, ExcInternalError());
              }
          }
      }

  // Do the same as above for cell-type data. This is simpler because it
  // is always scalar, and so we don't have to worry about whether some
  // components together form vectors or tensors.
  for (const auto &input_data : cell_data)
    {
      Assert(input_data->names.size() == 1, ExcInternalError());
      if ((input_data->is_complex_valued() == false) ||
          (input_data->postprocessor != nullptr))
        names.push_back(input_data->names[0]);
      else
        {
          names.push_back(input_data->names[0] + "_re");
          names.push_back(input_data->names[0] + "_im");
        }
    }

  return names;
}



template <typename DoFHandlerType, int patch_dim, int patch_space_dim>
std::vector<
  std::tuple<unsigned int,
             unsigned int,
             std::string,
             DataComponentInterpretation::DataComponentInterpretation>>
DataOut_DoFData<DoFHandlerType, patch_dim, patch_space_dim>::
  get_nonscalar_data_ranges() const
{
  std::vector<
    std::tuple<unsigned int,
               unsigned int,
               std::string,
               DataComponentInterpretation::DataComponentInterpretation>>
    ranges;

  // collect the ranges of dof and cell data
  unsigned int output_component = 0;
  for (const auto &input_data : dof_data)
    for (unsigned int i = 0; i < input_data->n_output_variables;
         /* i is updated below */)
      // see what kind of data we have here. note that for the purpose of the
      // current function all we care about is vector data
      switch (input_data->data_component_interpretation[i])
        {
          case DataComponentInterpretation::component_is_scalar:
            {
              // Just move component forward by one (or two if the
              // component happens to be complex-valued and we don't use a
              // postprocessor)
              // -- postprocessors always return real-valued things)
              ++i;
              output_component += (input_data->is_complex_valued() &&
                                       (input_data->postprocessor == nullptr) ?
                                     2 :
                                     1);

              break;
            }

          case DataComponentInterpretation::component_is_part_of_vector:
            {
              // ensure that there is a continuous number of next space_dim
              // components that all deal with vectors
              Assert(
                i + patch_space_dim <= input_data->n_output_variables,
                Exceptions::DataOutImplementation::ExcInvalidVectorDeclaration(
                  i, input_data->names[i]));
              for (unsigned int dd = 1; dd < patch_space_dim; ++dd)
                Assert(
                  input_data->data_component_interpretation[i + dd] ==
                    DataComponentInterpretation::component_is_part_of_vector,
                  Exceptions::DataOutImplementation::
                    ExcInvalidVectorDeclaration(i, input_data->names[i]));

              // all seems right, so figure out whether there is a common
              // name to these components. if not, leave the name empty and
              // let the output format writer decide what to do here
              std::string name = input_data->names[i];
              for (unsigned int dd = 1; dd < patch_space_dim; ++dd)
                if (name != input_data->names[i + dd])
                  {
                    name = "";
                    break;
                  }

              // Finally add a corresponding range. If this is a real-valued
              // vector, then we only need to do this once. But if it is a
              // complex-valued vector and it is not postprocessed, then we need
              // to do it twice -- once for the real parts and once for the
              // imaginary parts
              //
              // This sort of logic is also explained in some detail in
              //   DataOut::build_one_patch().
              if (input_data->is_complex_valued() == false ||
                  (input_data->postprocessor != nullptr))
                {
                  ranges.emplace_back(std::forward_as_tuple(
                    output_component,
                    output_component + patch_space_dim - 1,
                    name,
                    DataComponentInterpretation::component_is_part_of_vector));

                  // increase the 'component' counter by the appropriate amount,
                  // same for 'i', since we have already dealt with all these
                  // components
                  output_component += patch_space_dim;
                  i += patch_space_dim;
                }
              else
                {
                  ranges.emplace_back(std::forward_as_tuple(
                    output_component,
                    output_component + patch_space_dim - 1,
                    name + "_re",
                    DataComponentInterpretation::component_is_part_of_vector));
                  output_component += patch_space_dim;

                  ranges.emplace_back(std::forward_as_tuple(
                    output_component,
                    output_component + patch_space_dim - 1,
                    name + "_im",
                    DataComponentInterpretation::component_is_part_of_vector));
                  output_component += patch_space_dim;

                  i += patch_space_dim;
                }


              break;
            }

          case DataComponentInterpretation::component_is_part_of_tensor:
            {
              const unsigned int size = patch_space_dim * patch_space_dim;
              // ensure that there is a continuous number of next
              // space_dim*space_dim components that all deal with tensors
              Assert(
                i + size <= input_data->n_output_variables,
                Exceptions::DataOutImplementation::ExcInvalidTensorDeclaration(
                  i, input_data->names[i]));
              for (unsigned int dd = 1; dd < size; ++dd)
                Assert(
                  input_data->data_component_interpretation[i + dd] ==
                    DataComponentInterpretation::component_is_part_of_tensor,
                  Exceptions::DataOutImplementation::
                    ExcInvalidTensorDeclaration(i, input_data->names[i]));

              // all seems right, so figure out whether there is a common
              // name to these components. if not, leave the name empty and
              // let the output format writer decide what to do here
              std::string name = input_data->names[i];
              for (unsigned int dd = 1; dd < size; ++dd)
                if (name != input_data->names[i + dd])
                  {
                    name = "";
                    break;
                  }

              // Finally add a corresponding range. If this is a real-valued
              // tensor, then we only need to do this once. But if it is a
              // complex-valued tensor and it is not postprocessed, then we need
              // to do it twice -- once for the real parts and once for the
              // imaginary parts
              //
              // This sort of logic is also explained in some detail in
              //   DataOut::build_one_patch().
              if (input_data->is_complex_valued() == false ||
                  (input_data->postprocessor != nullptr))
                {
                  ranges.emplace_back(std::forward_as_tuple(
                    output_component,
                    output_component + size - 1,
                    name,
                    DataComponentInterpretation::component_is_part_of_tensor));

                  // increase the 'component' counter by the appropriate amount,
                  // same for 'i', since we have already dealt with all these
                  // components
                  output_component += size;
                  i += size;
                }
              else
                {
                  ranges.emplace_back(std::forward_as_tuple(
                    output_component,
                    output_component + size - 1,
                    name + "_re",
                    DataComponentInterpretation::component_is_part_of_tensor));
                  output_component += size;

                  ranges.emplace_back(std::forward_as_tuple(
                    output_component,
                    output_component + size - 1,
                    name + "_im",
                    DataComponentInterpretation::component_is_part_of_tensor));
                  output_component += size;

                  i += size;
                }
              break;
            }

          default:
            Assert(false, ExcNotImplemented());
        }

  // note that we do not have to traverse the list of cell data here because
  // cell data is one value per (logical) cell and therefore cannot be a
  // vector

  return ranges;
}



template <typename DoFHandlerType, int patch_dim, int patch_space_dim>
const std::vector<dealii::DataOutBase::Patch<patch_dim, patch_space_dim>> &
DataOut_DoFData<DoFHandlerType, patch_dim, patch_space_dim>::get_patches() const
{
  return patches;
}



template <typename DoFHandlerType, int patch_dim, int patch_space_dim>
std::vector<
  std::shared_ptr<dealii::hp::FECollection<DoFHandlerType::dimension,
                                           DoFHandlerType::space_dimension>>>
DataOut_DoFData<DoFHandlerType, patch_dim, patch_space_dim>::get_fes() const
{
  const unsigned int dhdim      = DoFHandlerType::dimension;
  const unsigned int dhspacedim = DoFHandlerType::space_dimension;
  std::vector<std::shared_ptr<dealii::hp::FECollection<dhdim, dhspacedim>>>
    finite_elements(this->dof_data.size());
  for (unsigned int i = 0; i < this->dof_data.size(); ++i)
    {
      Assert(dof_data[i]->dof_handler != nullptr,
             Exceptions::DataOutImplementation::ExcNoDoFHandlerSelected());

      // avoid creating too many finite elements and doing a lot of work on
      // initializing FEValues downstream: if two DoFHandlers are the same
      // (checked by pointer comparison), we can re-use the shared_ptr object
      // for the second one. We cannot check for finite element equalities
      // because we need different FEValues objects for different dof
      // handlers.
      bool duplicate = false;
      for (unsigned int j = 0; j < i; ++j)
        if (dof_data[i]->dof_handler == dof_data[j]->dof_handler)
          {
            finite_elements[i] = finite_elements[j];
            duplicate          = true;
          }
      if (duplicate == false)
        finite_elements[i] =
          std::make_shared<dealii::hp::FECollection<dhdim, dhspacedim>>(
            this->dof_data[i]->dof_handler->get_fe_collection());
    }
  if (this->dof_data.empty())
    {
      finite_elements.resize(1);
      finite_elements[0] =
        std::make_shared<dealii::hp::FECollection<dhdim, dhspacedim>>(
          FE_DGQ<dhdim, dhspacedim>(0));
    }
  return finite_elements;
}



template <typename DoFHandlerType, int patch_dim, int patch_space_dim>
std::size_t
DataOut_DoFData<DoFHandlerType, patch_dim, patch_space_dim>::
  memory_consumption() const
{
  return (DataOutInterface<patch_dim, patch_space_dim>::memory_consumption() +
          MemoryConsumption::memory_consumption(dofs) +
          MemoryConsumption::memory_consumption(patches));
}

DEAL_II_NAMESPACE_CLOSE

#endif
