// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2017 by the deal.II authors
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


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_P4EST

#  include <deal.II/distributed/solution_transfer.h>
#  include <deal.II/distributed/tria.h>

#  include <deal.II/dofs/dof_accessor.h>
#  include <deal.II/dofs/dof_tools.h>

#  include <deal.II/grid/tria_accessor.h>
#  include <deal.II/grid/tria_iterator.h>

#  include <deal.II/lac/block_vector.h>
#  include <deal.II/lac/la_parallel_block_vector.h>
#  include <deal.II/lac/la_parallel_vector.h>
#  include <deal.II/lac/petsc_parallel_block_vector.h>
#  include <deal.II/lac/petsc_parallel_vector.h>
#  include <deal.II/lac/trilinos_parallel_block_vector.h>
#  include <deal.II/lac/trilinos_vector.h>
#  include <deal.II/lac/vector.h>

#  include <functional>

DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  namespace distributed
  {
    template <int dim, typename VectorType, typename DoFHandlerType>
    SolutionTransfer<dim, VectorType, DoFHandlerType>::SolutionTransfer(
      const DoFHandlerType &dof)
      : dof_handler(&dof, typeid(*this).name())
      , handle(numbers::invalid_unsigned_int)
    {
      Assert(
        (dynamic_cast<const parallel::distributed::
                        Triangulation<dim, DoFHandlerType::space_dimension> *>(
           &dof_handler->get_triangulation()) != nullptr),
        ExcMessage(
          "parallel::distributed::SolutionTransfer requires a parallel::distributed::Triangulation object."));
    }



    template <int dim, typename VectorType, typename DoFHandlerType>
    void
    SolutionTransfer<dim, VectorType, DoFHandlerType>::
      prepare_for_coarsening_and_refinement(
        const std::vector<const VectorType *> &all_in)
    {
      input_vectors = all_in;
      register_data_attach();
    }



    template <int dim, typename VectorType, typename DoFHandlerType>
    void
    SolutionTransfer<dim, VectorType, DoFHandlerType>::register_data_attach()
    {
      // TODO: casting away constness is bad
      parallel::distributed::Triangulation<dim, DoFHandlerType::space_dimension>
        *tria = (dynamic_cast<parallel::distributed::Triangulation<
                   dim,
                   DoFHandlerType::space_dimension> *>(
          const_cast<dealii::Triangulation<dim, DoFHandlerType::space_dimension>
                       *>(&dof_handler->get_triangulation())));
      Assert(tria != nullptr, ExcInternalError());

      handle = tria->register_data_attach(std::bind(
        &SolutionTransfer<dim, VectorType, DoFHandlerType>::pack_callback,
        this,
        std::placeholders::_1,
        std::placeholders::_2));
    }



    template <int dim, typename VectorType, typename DoFHandlerType>
    void
    SolutionTransfer<dim, VectorType, DoFHandlerType>::
      prepare_for_coarsening_and_refinement(const VectorType &in)
    {
      std::vector<const VectorType *> all_in(1, &in);
      prepare_for_coarsening_and_refinement(all_in);
    }



    template <int dim, typename VectorType, typename DoFHandlerType>
    void
    SolutionTransfer<dim, VectorType, DoFHandlerType>::prepare_serialization(
      const VectorType &in)
    {
      std::vector<const VectorType *> all_in(1, &in);
      prepare_serialization(all_in);
    }



    template <int dim, typename VectorType, typename DoFHandlerType>
    void
    SolutionTransfer<dim, VectorType, DoFHandlerType>::prepare_serialization(
      const std::vector<const VectorType *> &all_in)
    {
      prepare_for_coarsening_and_refinement(all_in);
    }



    template <int dim, typename VectorType, typename DoFHandlerType>
    void
    SolutionTransfer<dim, VectorType, DoFHandlerType>::deserialize(
      VectorType &in)
    {
      std::vector<VectorType *> all_in(1, &in);
      deserialize(all_in);
    }



    template <int dim, typename VectorType, typename DoFHandlerType>
    void
    SolutionTransfer<dim, VectorType, DoFHandlerType>::deserialize(
      std::vector<VectorType *> &all_in)
    {
      register_data_attach();

      // this makes interpolate() happy
      input_vectors.resize(all_in.size());

      interpolate(all_in);
    }


    template <int dim, typename VectorType, typename DoFHandlerType>
    void
    SolutionTransfer<dim, VectorType, DoFHandlerType>::interpolate(
      std::vector<VectorType *> &all_out)
    {
      Assert(input_vectors.size() == all_out.size(),
             ExcDimensionMismatch(input_vectors.size(), all_out.size()));

      // TODO: casting away constness is bad
      parallel::distributed::Triangulation<dim, DoFHandlerType::space_dimension>
        *tria = (dynamic_cast<parallel::distributed::Triangulation<
                   dim,
                   DoFHandlerType::space_dimension> *>(
          const_cast<dealii::Triangulation<dim, DoFHandlerType::space_dimension>
                       *>(&dof_handler->get_triangulation())));
      Assert(tria != nullptr, ExcInternalError());

      tria->notify_ready_to_unpack(
        handle,
        std::bind(
          &SolutionTransfer<dim, VectorType, DoFHandlerType>::unpack_callback,
          this,
          std::placeholders::_1,
          std::placeholders::_2,
          std::placeholders::_3,
          std::ref(all_out)));


      for (typename std::vector<VectorType *>::iterator it = all_out.begin();
           it != all_out.end();
           ++it)
        (*it)->compress(::dealii::VectorOperation::insert);

      input_vectors.clear();
    }



    template <int dim, typename VectorType, typename DoFHandlerType>
    void
    SolutionTransfer<dim, VectorType, DoFHandlerType>::interpolate(
      VectorType &out)
    {
      std::vector<VectorType *> all_out(1, &out);
      interpolate(all_out);
    }



    template <int dim, typename VectorType, typename DoFHandlerType>
    std::vector<char>
    SolutionTransfer<dim, VectorType, DoFHandlerType>::pack_callback(
      const typename Triangulation<dim, DoFHandlerType::space_dimension>::
        cell_iterator &cell_,
      const typename Triangulation<dim, DoFHandlerType::space_dimension>::
        CellStatus /*status*/)
    {
      typename DoFHandlerType::cell_iterator cell(*cell_, dof_handler);

      const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;

      // create buffer for each individual object
      std::vector<::dealii::Vector<typename VectorType::value_type>> dofvalues(
        input_vectors.size());

      auto cit_input = input_vectors.cbegin();
      auto it_output = dofvalues.begin();
      for (; cit_input != input_vectors.cend(); ++cit_input, ++it_output)
        {
          it_output->reinit(dofs_per_cell);
          cell->get_interpolated_dof_values(*(*cit_input), *it_output);
        }

      // to get consistent data sizes on each cell for the fixed size transfer,
      // we won't allow compression
      return Utilities::pack(dofvalues, /*allow_compression=*/false);
    }



    template <int dim, typename VectorType, typename DoFHandlerType>
    void
    SolutionTransfer<dim, VectorType, DoFHandlerType>::unpack_callback(
      const typename Triangulation<dim, DoFHandlerType::space_dimension>::
        cell_iterator &cell_,
      const typename Triangulation<dim, DoFHandlerType::space_dimension>::
        CellStatus /*status*/,
      const boost::iterator_range<std::vector<char>::const_iterator>
        &                        data_range,
      std::vector<VectorType *> &all_out)
    {
      typename DoFHandlerType::cell_iterator cell(*cell_, dof_handler);

      const std::vector<::dealii::Vector<typename VectorType::value_type>>
        dofvalues = Utilities::unpack<
          std::vector<::dealii::Vector<typename VectorType::value_type>>>(
          data_range.begin(), data_range.end(), /*allow_compression=*/false);

      // check if sizes match
      Assert(dofvalues.size() == all_out.size(), ExcInternalError());

      // check if we have enough dofs provided by the FE object
      // to interpolate the transferred data correctly
      for (auto it_dofvalues = dofvalues.begin();
           it_dofvalues != dofvalues.end();
           ++it_dofvalues)
        Assert(
          cell->get_fe().dofs_per_cell == it_dofvalues->size(),
          ExcMessage(
            "The transferred data was packed with a different number of dofs than the"
            "currently registered FE object assigned to the DoFHandler has."));

      // distribute data for each registered vector on mesh
      auto it_input  = dofvalues.cbegin();
      auto it_output = all_out.begin();
      for (; it_input != dofvalues.cend(); ++it_input, ++it_output)
        {
          cell->set_dof_values_by_interpolation(*it_input, *(*it_output));
        }
    }


  } // namespace distributed
} // namespace parallel


// explicit instantiations
#  include "solution_transfer.inst"

DEAL_II_NAMESPACE_CLOSE

#endif
