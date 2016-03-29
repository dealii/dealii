// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2015 by the deal.II authors
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

#ifdef DEAL_II_WITH_P4EST

#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/parallel_vector.h>
#include <deal.II/lac/parallel_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>

#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/base/std_cxx11/bind.h>

DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  namespace distributed
  {

    template<int dim, typename VectorType, typename DoFHandlerType>
    SolutionTransfer<dim, VectorType, DoFHandlerType>::SolutionTransfer (const DoFHandlerType &dof)
      :
      dof_handler(&dof, typeid(*this).name())
    {
      Assert (dynamic_cast<const parallel::distributed::Triangulation<dim>*>
              (&dof_handler->get_triangulation()) != 0,
              ExcMessage("parallel::distributed::SolutionTransfer requires a parallel::distributed::Triangulation object."));
    }



    template<int dim, typename VectorType, typename DoFHandlerType>
    SolutionTransfer<dim, VectorType, DoFHandlerType>::~SolutionTransfer ()
    {}



    template<int dim, typename VectorType, typename DoFHandlerType>
    void
    SolutionTransfer<dim, VectorType, DoFHandlerType>::
    prepare_for_coarsening_and_refinement (const std::vector<const VectorType *> &all_in)
    {
      input_vectors = all_in;
      register_data_attach( get_data_size() * input_vectors.size() );
    }



    template<int dim, typename VectorType, typename DoFHandlerType>
    void
    SolutionTransfer<dim, VectorType, DoFHandlerType>::register_data_attach (const std::size_t size)
    {
      Assert(size > 0, ExcMessage("Please transfer at least one vector!"));

//TODO: casting away constness is bad
      parallel::distributed::Triangulation<dim,DoFHandlerType::space_dimension> *tria
        = (dynamic_cast<parallel::distributed::Triangulation<dim,DoFHandlerType::space_dimension>*>
           (const_cast<dealii::Triangulation<dim,DoFHandlerType::space_dimension>*>
            (&dof_handler->get_triangulation())));
      Assert (tria != 0, ExcInternalError());

      offset
        = tria->register_data_attach(size,
                                     std_cxx11::bind(&SolutionTransfer<dim, VectorType,
                                                     DoFHandlerType>::pack_callback,
                                                     this,
                                                     std_cxx11::_1,
                                                     std_cxx11::_2,
                                                     std_cxx11::_3));

    }



    template<int dim, typename VectorType, typename DoFHandlerType>
    void
    SolutionTransfer<dim, VectorType, DoFHandlerType>::
    prepare_for_coarsening_and_refinement (const VectorType &in)
    {
      std::vector<const VectorType *> all_in(1, &in);
      prepare_for_coarsening_and_refinement(all_in);
    }



    template<int dim, typename VectorType, typename DoFHandlerType>
    void
    SolutionTransfer<dim, VectorType, DoFHandlerType>::prepare_serialization (const VectorType &in)
    {
      std::vector<const VectorType *> all_in(1, &in);
      prepare_serialization(all_in);
    }



    template<int dim, typename VectorType, typename DoFHandlerType>
    void
    SolutionTransfer<dim, VectorType, DoFHandlerType>::prepare_serialization
    (const std::vector<const VectorType *> &all_in)
    {
      prepare_for_coarsening_and_refinement (all_in);
    }



    template<int dim, typename VectorType, typename DoFHandlerType>
    void
    SolutionTransfer<dim, VectorType, DoFHandlerType>::deserialize (VectorType &in)
    {
      std::vector<VectorType *> all_in(1, &in);
      deserialize(all_in);
    }



    template<int dim, typename VectorType, typename DoFHandlerType>
    void
    SolutionTransfer<dim, VectorType, DoFHandlerType>::deserialize (std::vector<VectorType *> &all_in)
    {
      register_data_attach( get_data_size() * all_in.size() );

      // this makes interpolate() happy
      input_vectors.resize(all_in.size());

      interpolate(all_in);
    }


    template<int dim, typename VectorType, typename DoFHandlerType>
    void
    SolutionTransfer<dim, VectorType, DoFHandlerType>::interpolate (std::vector<VectorType *> &all_out)
    {
      Assert(input_vectors.size()==all_out.size(),
             ExcDimensionMismatch(input_vectors.size(), all_out.size()) );

//TODO: casting away constness is bad
      parallel::distributed::Triangulation<dim,DoFHandlerType::space_dimension> *tria
        = (dynamic_cast<parallel::distributed::Triangulation<dim,DoFHandlerType::space_dimension>*>
           (const_cast<dealii::Triangulation<dim,DoFHandlerType::space_dimension>*>
            (&dof_handler->get_triangulation())));
      Assert (tria != 0, ExcInternalError());

      tria->notify_ready_to_unpack(offset,
                                   std_cxx11::bind(&SolutionTransfer<dim, VectorType,
                                                   DoFHandlerType>::unpack_callback,
                                                   this,
                                                   std_cxx11::_1,
                                                   std_cxx11::_2,
                                                   std_cxx11::_3,
                                                   std_cxx11::ref(all_out)));


      for (typename std::vector<VectorType *>::iterator it=all_out.begin();
           it !=all_out.end();
           ++it)
        (*it)->compress(::dealii::VectorOperation::insert);

      input_vectors.clear();
    }



    template<int dim, typename VectorType, typename DoFHandlerType>
    void
    SolutionTransfer<dim, VectorType, DoFHandlerType>::interpolate (VectorType &out)
    {
      std::vector<VectorType *> all_out(1, &out);
      interpolate(all_out);
    }



    template<int dim, typename VectorType, typename DoFHandlerType>
    unsigned int
    SolutionTransfer<dim, VectorType, DoFHandlerType>::get_data_size() const
    {
      return sizeof(typename VectorType::value_type)* DoFTools::max_dofs_per_cell(*dof_handler);
    }


    template<int dim, typename VectorType, typename DoFHandlerType>
    void
    SolutionTransfer<dim, VectorType, DoFHandlerType>::
    pack_callback(const typename Triangulation<dim,DoFHandlerType::space_dimension>::cell_iterator &cell_,
                  const typename Triangulation<dim,DoFHandlerType::space_dimension>::CellStatus /*status*/,
                  void *data)
    {
      typename VectorType::value_type *data_store = reinterpret_cast<typename VectorType::value_type *>(data);

      typename DoFHandlerType::cell_iterator cell(*cell_, dof_handler);

      const unsigned int dofs_per_cell=cell->get_fe().dofs_per_cell;
      ::dealii::Vector<typename VectorType::value_type> dofvalues(dofs_per_cell);
      for (typename std::vector<const VectorType *>::iterator it=input_vectors.begin();
           it !=input_vectors.end();
           ++it)
        {
          cell->get_interpolated_dof_values(*(*it), dofvalues);
          std::memcpy(data_store, &dofvalues(0), sizeof(typename VectorType::value_type)*dofs_per_cell);
          data_store += dofs_per_cell;
        }
    }


    template<int dim, typename VectorType, typename DoFHandlerType>
    void
    SolutionTransfer<dim, VectorType, DoFHandlerType>::unpack_callback
    (const typename Triangulation<dim,DoFHandlerType::space_dimension>::cell_iterator &cell_,
     const typename Triangulation<dim,DoFHandlerType::space_dimension>::CellStatus    /*status*/,
     const void                                           *data,
     std::vector<VectorType *>                            &all_out)
    {
      typename DoFHandlerType::cell_iterator
      cell(*cell_, dof_handler);

      const unsigned int dofs_per_cell=cell->get_fe().dofs_per_cell;
      ::dealii::Vector<typename VectorType::value_type> dofvalues(dofs_per_cell);
      const typename VectorType::value_type *data_store = reinterpret_cast<const typename VectorType::value_type *>(data);

      for (typename std::vector<VectorType *>::iterator it = all_out.begin();
           it != all_out.end();
           ++it)
        {
          std::memcpy(&dofvalues(0), data_store, sizeof(typename VectorType::value_type)*dofs_per_cell);
          cell->set_dof_values_by_interpolation(dofvalues, *(*it));
          data_store += dofs_per_cell;
        }
    }


  }
}


// explicit instantiations
#include "solution_transfer.inst"

DEAL_II_NAMESPACE_CLOSE

#endif
