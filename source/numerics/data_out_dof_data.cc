// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2013 by the deal.II authors
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

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/parallel_vector.h>
#include <deal.II/lac/parallel_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_out_dof_data.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <sstream>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace DataOut
  {
    template <int dim, int spacedim>
    ParallelDataBase<dim,spacedim>::
    ParallelDataBase (const unsigned int n_datasets,
                      const unsigned int n_subdivisions,
                      const std::vector<unsigned int> &n_postprocessor_outputs,
                      const Mapping<dim,spacedim> &mapping,
                      const std::vector<std_cxx11::shared_ptr<dealii::hp::FECollection<dim,spacedim> > > &finite_elements,
                      const UpdateFlags update_flags,
                      const bool        use_face_values)
      :
      n_datasets (n_datasets),
      n_subdivisions (n_subdivisions),
      postprocessed_values (n_postprocessor_outputs.size()),
      mapping_collection (mapping),
      finite_elements (finite_elements),
      update_flags (update_flags)
    {
      unsigned int n_q_points = 0;
      if (use_face_values == false)
        {
          dealii::hp::QCollection<dim>
          quadrature(QIterated<dim>(QTrapez<1>(), n_subdivisions));
          n_q_points = quadrature[0].size();
          x_fe_values.resize(this->finite_elements.size());
          for (unsigned int i=0; i<this->finite_elements.size(); ++i)
            {
              // check if there is a finite element that is equal to the
              // present one, then we can re-use the FEValues object
              for (unsigned int j=0; j<i; ++j)
                if (this->finite_elements[i].get() ==
                    this->finite_elements[j].get())
                  {
                    x_fe_values[i] = x_fe_values[j];
                    break;
                  }
              if (x_fe_values[i].get() == 0)
                x_fe_values[i].reset(new dealii::hp::FEValues<dim,spacedim>
                                     (this->mapping_collection,
                                      *this->finite_elements[i],
                                      quadrature,
                                      this->update_flags));
            }
        }
      else
        {
          dealii::hp::QCollection<dim-1>
          quadrature(QIterated<dim-1>(QTrapez<1>(), n_subdivisions));
          n_q_points = quadrature[0].size();
          x_fe_face_values.resize(this->finite_elements.size());
          for (unsigned int i=0; i<this->finite_elements.size(); ++i)
            {
              // check if there is a finite element that is equal to the
              // present one, then we can re-use the FEValues object
              for (unsigned int j=0; j<i; ++j)
                if (this->finite_elements[i].get() ==
                    this->finite_elements[j].get())
                  {
                    x_fe_face_values[i] = x_fe_face_values[j];
                    break;
                  }
              if (x_fe_face_values[i].get() == 0)
                x_fe_face_values[i].reset(new dealii::hp::FEFaceValues<dim,spacedim>
                                          (this->mapping_collection,
                                           *this->finite_elements[i],
                                           quadrature,
                                           this->update_flags));
            }
        }

      patch_values.resize (n_q_points);
      patch_values_system.resize (n_q_points);
      patch_gradients.resize (n_q_points);
      patch_gradients_system.resize (n_q_points);
      patch_hessians.resize (n_q_points);
      patch_hessians_system.resize (n_q_points);

      for (unsigned int dataset=0; dataset<n_postprocessor_outputs.size(); ++dataset)
        if (n_postprocessor_outputs[dataset] != 0)
          postprocessed_values[dataset]
          .resize(n_q_points,
                  dealii::Vector<double>(n_postprocessor_outputs[dataset]));
    }





    // implement copy constructor to create a thread's own version of
    // x_fe_values
    template <int dim, int spacedim>
    ParallelDataBase<dim,spacedim>::
    ParallelDataBase (const ParallelDataBase<dim,spacedim> &data)
      :
      n_datasets (data.n_datasets),
      n_subdivisions (data.n_subdivisions),
      patch_values (data.patch_values),
      patch_values_system (data.patch_values_system),
      patch_gradients (data.patch_gradients),
      patch_gradients_system (data.patch_gradients_system),
      patch_hessians (data.patch_hessians),
      patch_hessians_system (data.patch_hessians_system),
      postprocessed_values (data.postprocessed_values),
      mapping_collection (data.mapping_collection),
      finite_elements (data.finite_elements),
      update_flags (data.update_flags)
    {
      if (data.x_fe_values.empty() == false)
        {
          Assert(data.x_fe_face_values.empty() == true, ExcInternalError());
          dealii::hp::QCollection<dim>
          quadrature(QIterated<dim>(QTrapez<1>(), n_subdivisions));
          x_fe_values.resize(this->finite_elements.size());
          for (unsigned int i=0; i<this->finite_elements.size(); ++i)
            {
              // check if there is a finite element that is equal to the
              // present one, then we can re-use the FEValues object
              for (unsigned int j=0; j<i; ++j)
                if (this->finite_elements[i].get() ==
                    this->finite_elements[j].get())
                  {
                    x_fe_values[i] = x_fe_values[j];
                    break;
                  }
              if (x_fe_values[i].get() == 0)
                x_fe_values[i].reset(new dealii::hp::FEValues<dim,spacedim>
                                     (this->mapping_collection,
                                      *this->finite_elements[i],
                                      quadrature,
                                      this->update_flags));
            }
        }
      else
        {
          dealii::hp::QCollection<dim-1>
          quadrature(QIterated<dim-1>(QTrapez<1>(), n_subdivisions));
          x_fe_face_values.resize(this->finite_elements.size());
          for (unsigned int i=0; i<this->finite_elements.size(); ++i)
            {
              // check if there is a finite element that is equal to the
              // present one, then we can re-use the FEValues object
              for (unsigned int j=0; j<i; ++j)
                if (this->finite_elements[i].get() ==
                    this->finite_elements[j].get())
                  {
                    x_fe_face_values[i] = x_fe_face_values[j];
                    break;
                  }
              if (x_fe_face_values[i].get() == 0)
                x_fe_face_values[i].reset(new dealii::hp::FEFaceValues<dim,spacedim>
                                          (this->mapping_collection,
                                           *this->finite_elements[i],
                                           quadrature,
                                           this->update_flags));
            }
        }
    }



    template <int dim, int spacedim>
    template <typename DH>
    void
    ParallelDataBase<dim,spacedim>::
    reinit_all_fe_values(std::vector<std_cxx11::shared_ptr<DataEntryBase<DH> > > &dof_data,
                         const typename dealii::Triangulation<dim,spacedim>::cell_iterator &cell,
                         const unsigned int face)
    {
      for (unsigned int dataset=0; dataset<dof_data.size(); ++dataset)
        {
          bool duplicate = false;
          for (unsigned int j=0; j<dataset; ++j)
            if (finite_elements[dataset].get() == finite_elements[j].get())
              duplicate = true;
          if (duplicate == false)
            {
              typename DH::active_cell_iterator dh_cell(&cell->get_triangulation(),
                                                        cell->level(),
                                                        cell->index(),
                                                        dof_data[dataset]->dof_handler);
              if (x_fe_values.empty())
                {
                  AssertIndexRange(face,
                                   GeometryInfo<dim>::faces_per_cell);
                  x_fe_face_values[dataset]->reinit(dh_cell, face);
                }
              else
                x_fe_values[dataset]->reinit (dh_cell);
            }
        }
      if (dof_data.empty())
        {
          if (x_fe_values.empty())
            {
              AssertIndexRange(face,
                               GeometryInfo<dim>::faces_per_cell);
              x_fe_face_values[0]->reinit(cell, face);
            }
          else
            x_fe_values[0]->reinit (cell);
        }
    }



    template <int dim, int spacedim>
    const FEValuesBase<dim,spacedim> &
    ParallelDataBase<dim,spacedim>::
    get_present_fe_values(const unsigned int dataset) const
    {
      AssertIndexRange(dataset, finite_elements.size());
      if (x_fe_values.empty())
        return x_fe_face_values[dataset]->get_present_fe_values();
      else
        return x_fe_values[dataset]->get_present_fe_values();
    }



    template <int dim, int spacedim>
    void
    ParallelDataBase<dim,spacedim>::
    resize_system_vectors(const unsigned int n_components)
    {
      Assert(patch_values_system.size() > 0, ExcInternalError());
      AssertDimension(patch_values_system.size(),
                      patch_gradients_system.size());
      AssertDimension(patch_values_system.size(),
                      patch_hessians_system.size());
      if (patch_values_system[0].size() == n_components)
        return;
      for (unsigned int k=0; k<patch_values_system.size(); ++k)
        {
          patch_values_system[k].reinit(n_components);
          patch_gradients_system[k].resize(n_components);
          patch_hessians_system[k].resize(n_components);
        }
    }




    /**
     * In a WorkStream context, use this function to append the patch computed
     * by the parallel stage to the array of patches.
     */
    template <int dim, int spacedim>
    void
    append_patch_to_list (const DataOutBase::Patch<dim,spacedim> &patch,
                          std::vector<DataOutBase::Patch<dim,spacedim> > &patches)
    {
      patches.push_back (patch);
      patches.back().patch_index = patches.size()-1;
    }
  }
}

namespace internal
{
  namespace DataOut
  {
    template <class DH>
    DataEntryBase<DH>::DataEntryBase (const DH                       *dofs,
                                      const std::vector<std::string> &names_in,
                                      const std::vector<DataComponentInterpretation::DataComponentInterpretation> &data_component_interpretation)
      :
      dof_handler (dofs, typeid(dealii::DataOut_DoFData<DH,DH::dimension,DH::space_dimension>).name()),
      names(names_in),
      data_component_interpretation (data_component_interpretation),
      postprocessor(0, typeid(*this).name()),
      n_output_variables(names.size())
    {
      Assert (names.size() == data_component_interpretation.size(),
              ExcDimensionMismatch(data_component_interpretation.size(),
                                   names.size()));

      // check that the names use only allowed characters
      for (unsigned int i=0; i<names.size(); ++i)
        Assert (names[i].find_first_not_of("abcdefghijklmnopqrstuvwxyz"
                                           "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                           "0123456789_<>()") == std::string::npos,
                typename dealii::DataOut<DH::dimension>::
                ExcInvalidCharacter (names[i],
                                     names[i].find_first_not_of("abcdefghijklmnopqrstuvwxyz"
                                                                "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                                                "0123456789_<>()")));
    }



    template <class DH>
    DataEntryBase<DH>::DataEntryBase (const DH *dofs,
                                      const DataPostprocessor<DH::space_dimension> *data_postprocessor)
      :
      dof_handler (dofs, typeid(dealii::DataOut_DoFData<DH,DH::dimension,DH::space_dimension>).name()),
      names(data_postprocessor->get_names()),
      data_component_interpretation (data_postprocessor->get_data_component_interpretation()),
      postprocessor(data_postprocessor, typeid(*this).name()),
      n_output_variables(names.size())
    {
      Assert (data_postprocessor->get_names().size()
              ==
              data_postprocessor->get_data_component_interpretation().size(),
              ExcDimensionMismatch (data_postprocessor->get_names().size(),
                                    data_postprocessor->get_data_component_interpretation().size()));

      // check that the names use only allowed characters
      for (unsigned int i=0; i<names.size(); ++i)
        Assert (names[i].find_first_not_of("abcdefghijklmnopqrstuvwxyz"
                                           "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                           "0123456789_<>()") == std::string::npos,
                typename dealii::DataOut<DH::dimension>::
                ExcInvalidCharacter (names[i],
                                     names[i].find_first_not_of("abcdefghijklmnopqrstuvwxyz"
                                                                "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                                                "0123456789_<>()")));
    }



    template <class DH>
    DataEntryBase<DH>::~DataEntryBase ()
    {}



    /**
     * Class that stores a pointer to a vector of type equal to the template
     * argument, and provides the functions to extract data from it.
     *
     * @author Wolfgang Bangerth, 2004
     */
    template <class DH, typename VectorType>
    class DataEntry : public DataEntryBase<DH>
    {
    public:
      /**
       * Constructor. Give a list of names for the individual components of
       * the vector and their interpretation as scalar or vector data. This
       * constructor assumes that no postprocessor is going to be used.
       */
      DataEntry (const DH                       *dofs,
                 const VectorType               *data,
                 const std::vector<std::string> &names,
                 const std::vector<DataComponentInterpretation::DataComponentInterpretation> &data_component_interpretation);

      /**
       * Constructor when a data postprocessor is going to be used. In that
       * case, the names and vector declarations are going to be acquired from
       * the postprocessor.
       */
      DataEntry (const DH                                     *dofs,
                 const VectorType                             *data,
                 const DataPostprocessor<DH::space_dimension> *data_postprocessor);

      /**
       * Assuming that the stored vector is a cell vector, extract the given
       * element from it.
       */
      virtual
      double
      get_cell_data_value (const unsigned int cell_number) const;

      /**
       * Given a FEValuesBase object, extract the values on the present cell
       * from the vector we actually store.
       */
      virtual
      void
      get_function_values (const FEValuesBase<DH::dimension,DH::space_dimension> &fe_patch_values,
                           std::vector<double>             &patch_values) const;

      /**
       * Given a FEValuesBase object, extract the values on the present cell
       * from the vector we actually store. This function does the same as the
       * one above but for vector-valued finite elements.
       */
      virtual
      void
      get_function_values (const FEValuesBase<DH::dimension,DH::space_dimension> &fe_patch_values,
                           std::vector<dealii::Vector<double> >    &patch_values_system) const;

      /**
       * Given a FEValuesBase object, extract the gradients on the present
       * cell from the vector we actually store.
       */
      virtual
      void
      get_function_gradients (const FEValuesBase<DH::dimension,DH::space_dimension> &fe_patch_values,
                              std::vector<Tensor<1,DH::space_dimension> >       &patch_gradients) const;

      /**
       * Given a FEValuesBase object, extract the gradients on the present
       * cell from the vector we actually store. This function does the same
       * as the one above but for vector-valued finite elements.
       */
      virtual
      void
      get_function_gradients (const FEValuesBase<DH::dimension,DH::space_dimension> &fe_patch_values,
                              std::vector<std::vector<Tensor<1,DH::space_dimension> > > &patch_gradients_system) const;

      /**
       * Given a FEValuesBase object, extract the second derivatives on the
       * present cell from the vector we actually store.
       */
      virtual
      void
      get_function_hessians (const FEValuesBase<DH::dimension,DH::space_dimension> &fe_patch_values,
                             std::vector<Tensor<2,DH::space_dimension> >       &patch_hessians) const;

      /**
       * Given a FEValuesBase object, extract the second derivatives on the
       * present cell from the vector we actually store. This function does
       * the same as the one above but for vector-valued finite elements.
       */
      virtual
      void
      get_function_hessians (const FEValuesBase<DH::dimension,DH::space_dimension> &fe_patch_values,
                             std::vector<std::vector< Tensor<2,DH::space_dimension> > > &patch_hessians_system) const;

      /**
       * Clear all references to the vectors.
       */
      virtual void clear ();

      /**
       * Determine an estimate for the memory consumption (in bytes) of this
       * object.
       */
      virtual std::size_t memory_consumption () const;

    private:
      /**
       * Pointer to the data vector. Note that ownership of the vector pointed
       * to remains with the caller of this class.
       */
      const VectorType *vector;
    };



    template <class DH, class VectorType>
    DataEntry<DH,VectorType>::
    DataEntry (const DH                               *dofs,
               const VectorType                       *data,
               const std::vector<std::string>         &names,
               const std::vector<DataComponentInterpretation::DataComponentInterpretation> &data_component_interpretation)
      :
      DataEntryBase<DH> (dofs, names, data_component_interpretation),
      vector (data)
    {}



    template <class DH, class VectorType>
    DataEntry<DH,VectorType>::
    DataEntry (const DH                                     *dofs,
               const VectorType                             *data,
               const DataPostprocessor<DH::space_dimension> *data_postprocessor)
      :
      DataEntryBase<DH> (dofs, data_postprocessor),
      vector (data)
    {}


    namespace
    {
      template <class VectorType>
      double
      get_vector_element (const VectorType &vector,
                          const unsigned int cell_number)
      {
        return vector[cell_number];
      }


      double
      get_vector_element (const IndexSet &is,
                          const unsigned int cell_number)
      {
        return (is.is_element(cell_number) ? 1 : 0);
      }
    }



    template <class DH, class VectorType>
    double
    DataEntry<DH,VectorType>::
    get_cell_data_value (const unsigned int cell_number) const
    {
      return get_vector_element(*vector, cell_number);
    }



    template <class DH, class VectorType>
    void
    DataEntry<DH,VectorType>::
    get_function_values (const FEValuesBase<DH::dimension,DH::space_dimension> &fe_patch_values,
                         std::vector<dealii::Vector<double> >    &patch_values_system) const
    {
      fe_patch_values.get_function_values (*vector, patch_values_system);
    }



    template <class DH, typename VectorType>
    void
    DataEntry<DH,VectorType>::
    get_function_values (const FEValuesBase<DH::dimension,DH::space_dimension> &fe_patch_values,
                         std::vector<double>             &patch_values) const
    {
      fe_patch_values.get_function_values (*vector, patch_values);
    }



    template <class DH, class VectorType>
    void
    DataEntry<DH,VectorType>::
    get_function_gradients (const FEValuesBase<DH::dimension,DH::space_dimension> &fe_patch_values,
                            std::vector<std::vector<Tensor<1,DH::space_dimension> > >   &patch_gradients_system) const
    {
      fe_patch_values.get_function_gradients (*vector, patch_gradients_system);
    }



    template <class DH, typename VectorType>
    void
    DataEntry<DH,VectorType>::
    get_function_gradients (const FEValuesBase<DH::dimension,DH::space_dimension> &fe_patch_values,
                            std::vector<Tensor<1,DH::space_dimension> >       &patch_gradients) const
    {
      fe_patch_values.get_function_gradients (*vector, patch_gradients);
    }



    template <class DH, class VectorType>
    void
    DataEntry<DH,VectorType>::
    get_function_hessians (const FEValuesBase<DH::dimension,DH::space_dimension> &fe_patch_values,
                           std::vector<std::vector<Tensor<2,DH::space_dimension> > >   &patch_hessians_system) const
    {
      fe_patch_values.get_function_hessians (*vector, patch_hessians_system);
    }



    template <class DH, typename VectorType>
    void
    DataEntry<DH,VectorType>::
    get_function_hessians (const FEValuesBase<DH::dimension,DH::space_dimension> &fe_patch_values,
                           std::vector<Tensor<2,DH::space_dimension> >       &patch_hessians) const
    {
      fe_patch_values.get_function_hessians (*vector, patch_hessians);
    }



    template <class DH, typename VectorType>
    std::size_t
    DataEntry<DH,VectorType>::memory_consumption () const
    {
      return (sizeof (vector) +
              MemoryConsumption::memory_consumption (this->names));
    }



    template <class DH, class VectorType>
    void
    DataEntry<DH,VectorType>::clear ()
    {
      vector = 0;
      this->dof_handler = 0;
    }
  }
}



template <class DH,
          int patch_dim, int patch_space_dim>
DataOut_DoFData<DH,patch_dim,patch_space_dim>::DataOut_DoFData ()
  :
  triangulation(0,typeid(*this).name()),
  dofs(0,typeid(*this).name())
{}



template <class DH, int patch_dim, int patch_space_dim>
DataOut_DoFData<DH,patch_dim,patch_space_dim>::~DataOut_DoFData ()
{
  clear ();
}



template <class DH, int patch_dim, int patch_space_dim>
void
DataOut_DoFData<DH,patch_dim,patch_space_dim>::
attach_dof_handler (const DH &d)
{
  Assert (dof_data.size() == 0, ExcOldDataStillPresent());
  Assert (cell_data.size() == 0, ExcOldDataStillPresent());

  triangulation = SmartPointer<const Triangulation<DH::dimension,DH::space_dimension> >(&d.get_tria(), typeid(*this).name());
  dofs = SmartPointer<const DH>(&d, typeid(*this).name());
}



template <class DH, int patch_dim, int patch_space_dim>
void
DataOut_DoFData<DH,patch_dim,patch_space_dim>::
attach_triangulation (const Triangulation<DH::dimension,DH::space_dimension> &tria)
{
  Assert (dof_data.size() == 0, ExcOldDataStillPresent());
  Assert (cell_data.size() == 0, ExcOldDataStillPresent());

  triangulation = SmartPointer<const Triangulation<DH::dimension,DH::space_dimension> >(&tria, typeid(*this).name());
}




template <class DH,
          int patch_dim, int patch_space_dim>
template <class VECTOR>
void
DataOut_DoFData<DH,patch_dim,patch_space_dim>::
add_data_vector (const VECTOR                             &vec,
                 const std::string                        &name,
                 const DataVectorType                      type,
                 const std::vector<DataComponentInterpretation::DataComponentInterpretation> &data_component_interpretation)
{
  Assert (triangulation != 0, ExcNoTriangulationSelected ());
  const unsigned int n_components =
    dofs != 0 ? dofs->get_fe().n_components () : 1;

  std::vector<std::string> names;
  // if only one component or vector is cell vector: we only need one name
  if ((n_components == 1) ||
      (vec.size() == triangulation->n_active_cells()))
    {
      names.resize (1, name);
    }
  else
    // otherwise append _i to the given name
    {
      names.resize (n_components);
      for (unsigned int i=0; i<n_components; ++i)
        {
          std::ostringstream namebuf;
          namebuf << '_' << i;
          names[i] = name + namebuf.str();
        }
    }

  add_data_vector (vec, names, type, data_component_interpretation);
}



template <class DH,
          int patch_dim, int patch_space_dim>
template <class VECTOR>
void
DataOut_DoFData<DH,patch_dim,patch_space_dim>::
add_data_vector (const VECTOR                             &vec,
                 const std::vector<std::string>           &names,
                 const DataVectorType                      type,
                 const std::vector<DataComponentInterpretation::DataComponentInterpretation> &data_component_interpretation_)
{
  Assert (triangulation != 0, ExcNoTriangulationSelected ());

  const std::vector<DataComponentInterpretation::DataComponentInterpretation> &
  data_component_interpretation
    = (data_component_interpretation_.size() != 0
       ?
       data_component_interpretation_
       :
       std::vector<DataComponentInterpretation::DataComponentInterpretation>
       (names.size(), DataComponentInterpretation::component_is_scalar));

  // either cell data and one name,
  // or dof data and n_components names
  DataVectorType actual_type = type;
  if (type == type_automatic)
    {
      // in the rare case that someone has a DGP(0) attached, we can not decide what she wants here:
      Assert((dofs == 0) || (triangulation->n_active_cells() != dofs->n_dofs()),
             ExcMessage("Unable to determine the type of vector automatically because the number of DoFs "
                        "is equal to the number of cells. Please specify DataVectorType."));

      if (vec.size() == triangulation->n_active_cells())
        actual_type = type_cell_data;
      else
        actual_type = type_dof_data;
    }

  switch (actual_type)
    {
    case type_cell_data:
      Assert (vec.size() == triangulation->n_active_cells(),
              ExcDimensionMismatch (vec.size(),
                                    triangulation->n_active_cells()));
      Assert (names.size() == 1,
              ExcInvalidNumberOfNames (names.size(), 1));
      break;

    case type_dof_data:
      Assert (dofs != 0, ExcNoDoFHandlerSelected ());
      Assert (vec.size() == dofs->n_dofs(),
              ExcInvalidVectorSize (vec.size(),
                                    dofs->n_dofs(),
                                    triangulation->n_active_cells()));
      Assert (names.size() == dofs->get_fe().n_components(),
              ExcInvalidNumberOfNames (names.size(), dofs->get_fe().n_components()));
      break;

    case type_automatic:
      // this case should have been handled above...
      Assert (false, ExcInternalError());
    }

  internal::DataOut::DataEntryBase<DH> *new_entry
    = new internal::DataOut::DataEntry<DH,VECTOR>(dofs, &vec, names,
                                                  data_component_interpretation);
  if (actual_type == type_dof_data)
    dof_data.push_back (std_cxx11::shared_ptr<internal::DataOut::DataEntryBase<DH> >(new_entry));
  else
    cell_data.push_back (std_cxx11::shared_ptr<internal::DataOut::DataEntryBase<DH> >(new_entry));
}



template <class DH,
          int patch_dim, int patch_space_dim>
template <class VECTOR>
void
DataOut_DoFData<DH,patch_dim,patch_space_dim>::
add_data_vector (const VECTOR                           &vec,
                 const DataPostprocessor<DH::space_dimension> &data_postprocessor)
{
  // this is a specialized version of the other function where we have a
  // postprocessor. if we do, we know that we have type_dof_data, which makes
  // things a bit simpler, we also don't need to deal with some of the other
  // stuff and use a different constructor of DataEntry

  Assert (dofs != 0, ExcNoDoFHandlerSelected ());

  Assert (vec.size() == dofs->n_dofs(),
          ExcInvalidVectorSize (vec.size(),
                                dofs->n_dofs(),
                                dofs->get_tria().n_active_cells()));

  internal::DataOut::DataEntryBase<DH> *new_entry
    = new internal::DataOut::DataEntry<DH,VECTOR>(dofs, &vec, &data_postprocessor);
  dof_data.push_back (std_cxx11::shared_ptr<internal::DataOut::DataEntryBase<DH> >(new_entry));
}



template <class DH,
          int patch_dim, int patch_space_dim>
template <class VECTOR>
void
DataOut_DoFData<DH,patch_dim,patch_space_dim>::
add_data_vector (const DH                               &dof_handler,
                 const VECTOR                           &vec,
                 const DataPostprocessor<DH::space_dimension> &data_postprocessor)
{
  // this is a specialized version of the other function where we have a
  // postprocessor. if we do, we know that we have type_dof_data, which makes
  // things a bit simpler, we also don't need to deal with some of the other
  // stuff and use a different constructor of DataEntry

  AssertDimension (vec.size(), dof_handler.n_dofs());

  internal::DataOut::DataEntryBase<DH> *new_entry
    = new internal::DataOut::DataEntry<DH,VECTOR>(&dof_handler, &vec, &data_postprocessor);
  dof_data.push_back (std_cxx11::shared_ptr<internal::DataOut::DataEntryBase<DH> >(new_entry));
}



template <class DH,
          int patch_dim, int patch_space_dim>
template <class VECTOR>
void
DataOut_DoFData<DH,patch_dim,patch_space_dim>::
add_data_vector (const DH                       &dof_handler,
                 const VECTOR                   &data,
                 const std::string              &name,
                 const std::vector<DataComponentInterpretation::DataComponentInterpretation> &data_component_interpretation)
{
  const unsigned int n_components = dof_handler.get_fe().n_components ();

  std::vector<std::string> names;
  // if only one component: we only need one name
  if (n_components == 1)
    names.resize (1, name);
  else
    // otherwise append _i to the given name
    {
      names.resize (n_components);
      for (unsigned int i=0; i<n_components; ++i)
        {
          std::ostringstream namebuf;
          namebuf << '_' << i;
          names[i] = name + namebuf.str();
        }
    }

  add_data_vector (dof_handler, data, names, data_component_interpretation);
}



template <class DH,
          int patch_dim, int patch_space_dim>
template <class VECTOR>
void
DataOut_DoFData<DH,patch_dim,patch_space_dim>::
add_data_vector (const DH                       &dof_handler,
                 const VECTOR                   &data,
                 const std::vector<std::string> &names,
                 const std::vector<DataComponentInterpretation::DataComponentInterpretation> &data_component_interpretation_)
{
  // this is an extended version of the other functions where we pass a vector
  // together with its DoFHandler. if we do, we know that we have
  // type_dof_data, which makes things a bit simpler
  if (triangulation == 0)
    triangulation = SmartPointer<const Triangulation<DH::dimension,DH::space_dimension> >(&dof_handler.get_tria(), typeid(*this).name());

  Assert (&dof_handler.get_tria() == triangulation,
          ExcMessage("The triangulation attached to the DoFHandler does not "
                     "match with the one set previously"));

  Assert (data.size() == dof_handler.n_dofs(),
          ExcDimensionMismatch (data.size(), dof_handler.n_dofs()));

  const std::vector<DataComponentInterpretation::DataComponentInterpretation> &
  data_component_interpretation
    = (data_component_interpretation_.size() != 0
       ?
       data_component_interpretation_
       :
       std::vector<DataComponentInterpretation::DataComponentInterpretation>
       (names.size(), DataComponentInterpretation::component_is_scalar));

  internal::DataOut::DataEntryBase<DH> *new_entry
    = new internal::DataOut::DataEntry<DH,VECTOR>(&dof_handler, &data, names,
                                                  data_component_interpretation);
  dof_data.push_back (std_cxx11::shared_ptr<internal::DataOut::DataEntryBase<DH> >(new_entry));
}



template <class DH,
          int patch_dim, int patch_space_dim>
void DataOut_DoFData<DH,patch_dim,patch_space_dim>::clear_data_vectors ()
{
  dof_data.erase (dof_data.begin(), dof_data.end());
  cell_data.erase (cell_data.begin(), cell_data.end());

  // delete patches
  std::vector<Patch> dummy;
  patches.swap (dummy);
}



template <class DH,
          int patch_dim, int patch_space_dim>
void
DataOut_DoFData<DH,patch_dim,patch_space_dim>::
clear_input_data_references ()
{
  for (unsigned int i=0; i<dof_data.size(); ++i)
    dof_data[i]->clear ();

  for (unsigned int i=0; i<cell_data.size(); ++i)
    cell_data[i]->clear ();

  if (dofs != 0)
    dofs = 0;
}



template <class DH,
          int patch_dim, int patch_space_dim>
void
DataOut_DoFData<DH,patch_dim,patch_space_dim>::clear ()
{
  dof_data.erase (dof_data.begin(), dof_data.end());
  cell_data.erase (cell_data.begin(), cell_data.end());

  if (dofs != 0)
    dofs = 0;

  // delete patches
  std::vector<Patch> dummy;
  patches.swap (dummy);
}



template <class DH,
          int patch_dim, int patch_space_dim>
std::vector<std::string>
DataOut_DoFData<DH,patch_dim,patch_space_dim>::
get_dataset_names () const
{
  std::vector<std::string> names;
  // collect the names of dof
  // and cell data
  typedef
  typename std::vector<std_cxx11::shared_ptr<internal::DataOut::DataEntryBase<DH> > >::const_iterator
  data_iterator;

  for (data_iterator  d=dof_data.begin();
       d!=dof_data.end(); ++d)
    for (unsigned int i=0; i<(*d)->names.size(); ++i)
      names.push_back ((*d)->names[i]);
  for (data_iterator d=cell_data.begin(); d!=cell_data.end(); ++d)
    {
      Assert ((*d)->names.size() == 1, ExcInternalError());
      names.push_back ((*d)->names[0]);
    }

  return names;
}



template <class DH,
          int patch_dim, int patch_space_dim>
std::vector<std_cxx11::tuple<unsigned int, unsigned int, std::string> >
DataOut_DoFData<DH,patch_dim,patch_space_dim>::get_vector_data_ranges () const
{
  std::vector<std_cxx11::tuple<unsigned int, unsigned int, std::string> >
  ranges;

  // collect the ranges of dof
  // and cell data
  typedef
  typename std::vector<std_cxx11::shared_ptr<internal::DataOut::DataEntryBase<DH> > >::const_iterator
  data_iterator;

  unsigned int output_component = 0;
  for (data_iterator  d=dof_data.begin();
       d!=dof_data.end(); ++d)
    for (unsigned int i=0; i<(*d)->n_output_variables;
         ++i, ++output_component)
      // see what kind of data we have
      // here. note that for the purpose of
      // the current function all we care
      // about is vector data
      if ((*d)->data_component_interpretation[i] ==
          DataComponentInterpretation::component_is_part_of_vector)
        {
          // ensure that there is a
          // continuous number of next
          // space_dim components that all
          // deal with vectors
          Assert (i+patch_space_dim <=
                  (*d)->n_output_variables,
                  ExcInvalidVectorDeclaration (i,
                                               (*d)->names[i]));
          for (unsigned int dd=1; dd<patch_space_dim; ++dd)
            Assert ((*d)->data_component_interpretation[i+dd]
                    ==
                    DataComponentInterpretation::component_is_part_of_vector,
                    ExcInvalidVectorDeclaration (i,
                                                 (*d)->names[i]));

          // all seems alright, so figure out
          // whether there is a common name
          // to these components. if not,
          // leave the name empty and let the
          // output format writer decide what
          // to do here
          std::string name = (*d)->names[i];
          for (unsigned int dd=1; dd<patch_space_dim; ++dd)
            if (name != (*d)->names[i+dd])
              {
                name = "";
                break;
              }

          // finally add a corresponding
          // range
          std_cxx11::tuple<unsigned int, unsigned int, std::string>
          range (output_component,
                 output_component+patch_space_dim-1,
                 name);

          ranges.push_back (range);

          // increase the 'component' counter
          // by the appropriate amount, same
          // for 'i', since we have already
          // dealt with all these components
          output_component += patch_space_dim-1;
          i += patch_space_dim-1;
        }

  // note that we do not have to traverse the
  // list of cell data here because cell data
  // is one value per (logical) cell and
  // therefore cannot be a vector

  // as a final check, the 'component'
  // counter should be at the total number of
  // components added up now
#ifdef DEBUG
  unsigned int n_output_components = 0;
  for (data_iterator  d=dof_data.begin();
       d!=dof_data.end(); ++d)
    n_output_components += (*d)->n_output_variables;
  Assert (output_component == n_output_components,
          ExcInternalError());
#endif

  return ranges;
}



template <class DH,
          int patch_dim, int patch_space_dim>
const std::vector< dealii::DataOutBase::Patch<patch_dim, patch_space_dim> > &
DataOut_DoFData<DH,patch_dim,patch_space_dim>::get_patches () const
{
  return patches;
}



template <class DH,
          int patch_dim, int patch_space_dim>
std::vector<std_cxx11::shared_ptr<dealii::hp::FECollection<DH::dimension,DH::space_dimension> > >
DataOut_DoFData<DH,patch_dim,patch_space_dim>::get_finite_elements() const
{
  const unsigned int dhdim = DH::dimension;
  const unsigned int dhspacedim = DH::space_dimension;
  std::vector<std_cxx11::shared_ptr<dealii::hp::FECollection<dhdim,dhspacedim> > >
  finite_elements(this->dof_data.size());
  for (unsigned int i=0; i<this->dof_data.size(); ++i)
    {
      Assert (dof_data[i]->dof_handler != 0, ExcNoDoFHandlerSelected ());

      // avoid creating too many finite elements and doing a lot of work on
      // initializing FEValues downstream: if two DoFHandlers are the same
      // (checked by pointer comparison), we can re-use the shared_ptr object
      // for the second one. We cannot check for finite element equalities
      // because we need different FEValues objects for different dof
      // handlers.
      bool duplicate = false;
      for (unsigned int j=0; j<i; ++j)
        if (dof_data[i]->dof_handler == dof_data[j]->dof_handler)
          {
            finite_elements[i] = finite_elements[j];
            duplicate = true;
          }
      if (duplicate == false)
        finite_elements[i].reset(new dealii::hp::FECollection<dhdim,dhspacedim>
                                 (this->dof_data[i]->dof_handler->get_fe()));
    }
  if (this->dof_data.empty())
    {
      finite_elements.resize(1);
      finite_elements[0].reset(new dealii::hp::FECollection<dhdim,dhspacedim>
                               (FE_DGQ<dhdim,dhspacedim>(0)));
    }
  return finite_elements;
}



template <class DH,
          int patch_dim, int patch_space_dim>
std::size_t
DataOut_DoFData<DH,patch_dim,patch_space_dim>::memory_consumption () const
{
  return (DataOutInterface<patch_dim,patch_space_dim>::memory_consumption () +
          MemoryConsumption::memory_consumption (dofs) +
          MemoryConsumption::memory_consumption (patches));
}



// explicit instantiations
#include "data_out_dof_data.inst"

DEAL_II_NAMESPACE_CLOSE
