// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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

#  include <deal.II/distributed/error_predictor.h>
#  include <deal.II/distributed/tria.h>

#  include <deal.II/dofs/dof_accessor.h>
#  include <deal.II/dofs/dof_tools.h>

#  include <deal.II/grid/tria_accessor.h>
#  include <deal.II/grid/tria_iterator.h>

#  include <deal.II/hp/dof_handler.h>

#  include <deal.II/lac/block_vector.h>
#  include <deal.II/lac/la_parallel_block_vector.h>
#  include <deal.II/lac/la_parallel_vector.h>
#  include <deal.II/lac/petsc_block_vector.h>
#  include <deal.II/lac/petsc_vector.h>
#  include <deal.II/lac/trilinos_parallel_block_vector.h>
#  include <deal.II/lac/trilinos_vector.h>
#  include <deal.II/lac/vector.h>

#  include <functional>
#  include <numeric>


DEAL_II_NAMESPACE_OPEN


namespace parallel
{
  namespace distributed
  {
    template <int dim, int spacedim>
    ErrorPredictor<dim, spacedim>::ErrorPredictor(
      const hp::DoFHandler<dim, spacedim> &dof)
      : dof_handler(&dof, typeid(*this).name())
      , handle(numbers::invalid_unsigned_int)
      , gamma_p(0.)
      , gamma_h(0.)
      , gamma_n(0.)
    {
      Assert(
        (dynamic_cast<const parallel::distributed::Triangulation<dim, spacedim>
                        *>(&dof_handler->get_triangulation()) != nullptr),
        ExcMessage(
          "parallel::distributed::ErrorPredictor requires a parallel::distributed::Triangulation object."));
    }



    template <int dim, int spacedim>
    void
    ErrorPredictor<dim, spacedim>::prepare_for_coarsening_and_refinement(
      const std::vector<const Vector<float> *> &all_in,
      const double                              gamma_p,
      const double                              gamma_h,
      const double                              gamma_n)
    {
      error_indicators = all_in;
      this->gamma_p    = gamma_p;
      this->gamma_h    = gamma_h;
      this->gamma_n    = gamma_n;
      register_data_attach();
    }



    template <int dim, int spacedim>
    void
    ErrorPredictor<dim, spacedim>::register_data_attach()
    {
      // TODO: casting away constness is bad
      parallel::distributed::Triangulation<dim, spacedim> *tria =
        (dynamic_cast<parallel::distributed::Triangulation<dim, spacedim> *>(
          const_cast<dealii::Triangulation<dim, spacedim> *>(
            &dof_handler->get_triangulation())));
      Assert(tria != nullptr, ExcInternalError());

      handle = tria->register_data_attach(
        [this](const typename Triangulation<dim, spacedim>::cell_iterator &cell,
               const typename Triangulation<dim, spacedim>::CellStatus status) {
          return this->pack_callback(cell, status);
        },
        /*returns_variable_size_data=*/false);
    }



    template <int dim, int spacedim>
    void
    ErrorPredictor<dim, spacedim>::prepare_for_coarsening_and_refinement(
      const Vector<float> &in,
      const double         gamma_p,
      const double         gamma_h,
      const double         gamma_n)
    {
      std::vector<const Vector<float> *> all_in(1, &in);
      prepare_for_coarsening_and_refinement(all_in, gamma_p, gamma_h, gamma_n);
    }



    template <int dim, int spacedim>
    void
    ErrorPredictor<dim, spacedim>::unpack(std::vector<Vector<float> *> &all_out)
    {
      Assert(error_indicators.size() == all_out.size(),
             ExcDimensionMismatch(error_indicators.size(), all_out.size()));

      // TODO: casting away constness is bad
      parallel::distributed::Triangulation<dim, spacedim> *tria =
        (dynamic_cast<parallel::distributed::Triangulation<dim, spacedim> *>(
          const_cast<dealii::Triangulation<dim, spacedim> *>(
            &dof_handler->get_triangulation())));
      Assert(tria != nullptr, ExcInternalError());

      tria->notify_ready_to_unpack(
        handle,
        [this, &all_out](
          const typename Triangulation<dim, spacedim>::cell_iterator &cell,
          const typename Triangulation<dim, spacedim>::CellStatus     status,
          const boost::iterator_range<std::vector<char>::const_iterator>
            &data_range) {
          this->unpack_callback(cell, status, data_range, all_out);
        });

      for (const auto &out : all_out)
        out->compress(::dealii::VectorOperation::insert);

      error_indicators.clear();
    }



    template <int dim, int spacedim>
    void
    ErrorPredictor<dim, spacedim>::unpack(Vector<float> &out)
    {
      std::vector<Vector<float> *> all_out(1, &out);
      unpack(all_out);
    }



    template <int dim, int spacedim>
    std::vector<char>
    ErrorPredictor<dim, spacedim>::pack_callback(
      const typename Triangulation<dim, spacedim>::cell_iterator &cell_,
      const typename Triangulation<dim, spacedim>::CellStatus     status)
    {
      typename hp::DoFHandler<dim, spacedim>::cell_iterator cell(*cell_,
                                                                 dof_handler);

      // create buffer for each individual input vector
      std::vector<float> predicted_errors(error_indicators.size());

      auto predicted_error_it = predicted_errors.begin();
      auto estimated_error_it = error_indicators.cbegin();
      for (; estimated_error_it != error_indicators.cend();
           ++predicted_error_it, ++estimated_error_it)
        switch (status)
          {
            case parallel::distributed::Triangulation<dim,
                                                      spacedim>::CELL_PERSIST:
              {
                if (cell->future_fe_index_set())
                  {
                    const int degree_difference =
                      dof_handler->get_fe_collection()[cell->future_fe_index()]
                        .degree -
                      cell->get_fe().degree;

                    *predicted_error_it =
                      (**estimated_error_it)[cell->active_cell_index()] *
                      std::pow(gamma_p, degree_difference);
                  }
                else
                  {
                    *predicted_error_it =
                      (**estimated_error_it)[cell->active_cell_index()] *
                      gamma_n;
                  }
                break;
              }

            case parallel::distributed::Triangulation<dim,
                                                      spacedim>::CELL_REFINE:
              {
                // Determine the exponent by the finite element degree on the
                // adapted mesh.
                const unsigned int future_fe_degree =
                  dof_handler->get_fe_collection()[cell->future_fe_index()]
                    .degree;

                *predicted_error_it =
                  (**estimated_error_it)[cell->active_cell_index()] *
                  (gamma_h * std::pow(.5, future_fe_degree));

                // If the future fe index differs from the active one, also take
                // into account p-adaptation.
                if (cell->future_fe_index_set())
                  *predicted_error_it *=
                    std::pow(gamma_p,
                             static_cast<int>(future_fe_degree -
                                              cell->get_fe().degree));

                break;
              }

            case parallel::distributed::Triangulation<dim,
                                                      spacedim>::CELL_COARSEN:
              {
                // First figure out which finite element will be assigned to the
                // parent cell after h-adaptation analogously to
                // dealii::internal::hp::DoFHandlerImplementation::Implementation::
                //   collect_fe_indices_on_cells_to_be_refined()
                std::set<unsigned int> fe_indices_children;
                for (unsigned int child_index = 0;
                     child_index < cell->n_children();
                     ++child_index)
                  {
                    const auto child = cell->child(child_index);
                    Assert(child->is_active() && child->coarsen_flag_set(),
                           typename dealii::Triangulation<
                             dim>::ExcInconsistentCoarseningFlags());

                    fe_indices_children.insert(child->future_fe_index());
                  }

                const unsigned int future_fe_index =
                  dof_handler->get_fe_collection().find_dominated_fe_extended(
                    fe_indices_children, /*codim=*/0);

                const unsigned int future_fe_degree =
                  dof_handler->get_fe_collection()[future_fe_index].degree;

                // Then determine the actually contirbution to the predicted
                // error of every single cells that is about to be coarsened.
                float sqrsum_of_predicted_errors = 0.;
                float predicted_error            = 0.;
                int   degree_difference          = 0;
                for (unsigned int child_index = 0;
                     child_index < cell->n_children();
                     ++child_index)
                  {
                    const auto child = cell->child(child_index);

                    predicted_error =
                      (**estimated_error_it)[child->active_cell_index()] /
                      (gamma_h * std::pow(.5, future_fe_degree));

                    degree_difference =
                      future_fe_degree - child->get_fe().degree;

                    if (degree_difference != 0)
                      predicted_error *= std::pow(gamma_p, degree_difference);

                    sqrsum_of_predicted_errors +=
                      predicted_error * predicted_error;
                  }
                *predicted_error_it = std::sqrt(sqrsum_of_predicted_errors);

                break;
              }

            default:
              Assert(false, ExcInternalError());
              break;
          }

      // We don't have to pack the whole container if there is just one entry.
      if (error_indicators.size() == 1)
        return Utilities::pack(predicted_errors[0],
                               /*allow_compression=*/false);
      else
        return Utilities::pack(predicted_errors, /*allow_compression=*/false);
    }



    template <int dim, int spacedim>
    void
    ErrorPredictor<dim, spacedim>::unpack_callback(
      const typename Triangulation<dim, spacedim>::cell_iterator &cell,
      const typename Triangulation<dim, spacedim>::CellStatus     status,
      const boost::iterator_range<std::vector<char>::const_iterator>
        &                           data_range,
      std::vector<Vector<float> *> &all_out)
    {
      std::vector<float> predicted_errors;

      if (all_out.size() == 1)
        predicted_errors.push_back(
          Utilities::unpack<float>(data_range.begin(),
                                   data_range.end(),
                                   /*allow_compression=*/false));
      else
        predicted_errors =
          Utilities::unpack<std::vector<float>>(data_range.begin(),
                                                data_range.end(),
                                                /*allow_compression=*/false);

      Assert(predicted_errors.size() == all_out.size(), ExcInternalError());

      auto it_input  = predicted_errors.cbegin();
      auto it_output = all_out.begin();
      for (; it_input != predicted_errors.cend(); ++it_input, ++it_output)
        switch (status)
          {
            case parallel::distributed::Triangulation<dim,
                                                      spacedim>::CELL_PERSIST:
            case parallel::distributed::Triangulation<dim,
                                                      spacedim>::CELL_COARSEN:
              (**it_output)[cell->active_cell_index()] = *it_input;
              break;


            case parallel::distributed::Triangulation<dim,
                                                      spacedim>::CELL_REFINE:
              for (unsigned int child_index = 0;
                   child_index < cell->n_children();
                   ++child_index)
                (**it_output)[cell->child(child_index)->active_cell_index()] =
                  (*it_input) / std::sqrt(cell->n_children());
              break;

            default:
              Assert(false, ExcInternalError());
              break;
          }
    }
  } // namespace distributed
} // namespace parallel


// explicit instantiations
#  include "error_predictor.inst"

DEAL_II_NAMESPACE_CLOSE

#endif /* DEAL_II_WITH_P4EST */
