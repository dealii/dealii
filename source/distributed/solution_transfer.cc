// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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

    template<int dim, typename VECTOR, class DH>
    SolutionTransfer<dim, VECTOR, DH>::SolutionTransfer(const DH &dof)
      :
      dof_handler(&dof, typeid(*this).name())
    {}



    template<int dim, typename VECTOR, class DH>
    SolutionTransfer<dim, VECTOR, DH>::~SolutionTransfer()
    {}



    template<int dim, typename VECTOR, class DH>
    void
    SolutionTransfer<dim, VECTOR, DH>::
    prepare_for_coarsening_and_refinement (const std::vector<const VECTOR *> &all_in)
    {
      input_vectors = all_in;
      register_data_attach( get_data_size() * input_vectors.size() );
    }



    template<int dim, typename VECTOR, class DH>
    void
    SolutionTransfer<dim, VECTOR, DH>::
    register_data_attach(const std::size_t size)
    {
      Assert(size > 0, ExcMessage("Please transfer at least one vector!"));

//TODO: casting away constness is bad
      parallel::distributed::Triangulation<dim> *tria
        = (dynamic_cast<parallel::distributed::Triangulation<dim>*>
           (const_cast<dealii::Triangulation<dim>*>
            (&dof_handler->get_tria())));
      Assert (tria != 0, ExcInternalError());

      offset
        = tria->register_data_attach(size,
                                     std_cxx11::bind(&SolutionTransfer<dim, VECTOR, DH>::pack_callback,
                                                     this,
                                                     std_cxx11::_1,
                                                     std_cxx11::_2,
                                                     std_cxx11::_3));

    }



    template<int dim, typename VECTOR, class DH>
    void
    SolutionTransfer<dim, VECTOR, DH>::
    prepare_for_coarsening_and_refinement (const VECTOR &in)
    {
      std::vector<const VECTOR *> all_in(1, &in);
      prepare_for_coarsening_and_refinement(all_in);
    }



    template<int dim, typename VECTOR, class DH>
    void
    SolutionTransfer<dim, VECTOR, DH>::
    prepare_serialization(const VECTOR &in)
    {
      std::vector<const VECTOR *> all_in(1, &in);
      prepare_serialization(all_in);
    }



    template<int dim, typename VECTOR, class DH>
    void
    SolutionTransfer<dim, VECTOR, DH>::
    prepare_serialization(const std::vector<const VECTOR *> &all_in)
    {
      prepare_for_coarsening_and_refinement (all_in);
    }



    template<int dim, typename VECTOR, class DH>
    void
    SolutionTransfer<dim, VECTOR, DH>::
    deserialize(VECTOR &in)
    {
      std::vector<VECTOR *> all_in(1, &in);
      deserialize(all_in);
    }



    template<int dim, typename VECTOR, class DH>
    void
    SolutionTransfer<dim, VECTOR, DH>::
    deserialize(std::vector<VECTOR *> &all_in)
    {
      register_data_attach( get_data_size() * all_in.size() );

      // this makes interpolate() happy
      input_vectors.resize(all_in.size());

      interpolate(all_in);
    }


    template<int dim, typename VECTOR, class DH>
    void
    SolutionTransfer<dim, VECTOR, DH>::
    interpolate (std::vector<VECTOR *> &all_out)
    {
      Assert(input_vectors.size()==all_out.size(),
             ExcDimensionMismatch(input_vectors.size(), all_out.size()) );

//TODO: casting away constness is bad
      parallel::distributed::Triangulation<dim> *tria
        = (dynamic_cast<parallel::distributed::Triangulation<dim>*>
           (const_cast<dealii::Triangulation<dim>*>
            (&dof_handler->get_tria())));
      Assert (tria != 0, ExcInternalError());

      tria->notify_ready_to_unpack(offset,
                                   std_cxx11::bind(&SolutionTransfer<dim, VECTOR, DH>::unpack_callback,
                                                   this,
                                                   std_cxx11::_1,
                                                   std_cxx11::_2,
                                                   std_cxx11::_3,
                                                   std_cxx11::ref(all_out)));


      for (typename std::vector<VECTOR *>::iterator it=all_out.begin();
           it !=all_out.end();
           ++it)
        (*it)->compress(::dealii::VectorOperation::insert);

      input_vectors.clear();
    }



    template<int dim, typename VECTOR, class DH>
    void
    SolutionTransfer<dim, VECTOR, DH>::
    interpolate (VECTOR &out)
    {
      std::vector<VECTOR *> all_out(1, &out);
      interpolate(all_out);
    }



    template<int dim, typename VECTOR, class DH>
    unsigned int
    SolutionTransfer<dim, VECTOR, DH>::
    get_data_size() const
    {
      return sizeof(double)* DoFTools::max_dofs_per_cell(*dof_handler);
    }


    template<int dim, typename VECTOR, class DH>
    void
    SolutionTransfer<dim, VECTOR, DH>::
    pack_callback(const typename Triangulation<dim,dim>::cell_iterator &cell_,
                  const typename Triangulation<dim,dim>::CellStatus /*status*/,
                  void *data)
    {
      double *data_store = reinterpret_cast<double *>(data);

      typename DH::cell_iterator cell(*cell_, dof_handler);

      const unsigned int dofs_per_cell=cell->get_fe().dofs_per_cell;
      ::dealii::Vector<double> dofvalues(dofs_per_cell);
      for (typename std::vector<const VECTOR *>::iterator it=input_vectors.begin();
           it !=input_vectors.end();
           ++it)
        {
          cell->get_interpolated_dof_values(*(*it), dofvalues);
          std::memcpy(data_store, &dofvalues(0), sizeof(double)*dofs_per_cell);
          data_store += dofs_per_cell;
        }
    }


    template<int dim, typename VECTOR, class DH>
    void
    SolutionTransfer<dim, VECTOR, DH>::
    unpack_callback(const typename Triangulation<dim,dim>::cell_iterator &cell_,
                    const typename Triangulation<dim,dim>::CellStatus /*status*/,
                    const void *data,
                    std::vector<VECTOR *> &all_out)
    {
      typename DH::cell_iterator
      cell(*cell_, dof_handler);

      const unsigned int dofs_per_cell=cell->get_fe().dofs_per_cell;
      ::dealii::Vector<double> dofvalues(dofs_per_cell);
      const double *data_store = reinterpret_cast<const double *>(data);

      for (typename std::vector<VECTOR *>::iterator it = all_out.begin();
           it != all_out.end();
           ++it)
        {
          std::memcpy(&dofvalues(0), data_store, sizeof(double)*dofs_per_cell);
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
