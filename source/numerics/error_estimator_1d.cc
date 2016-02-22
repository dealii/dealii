// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2016 by the deal.II authors
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

#include <deal.II/base/thread_management.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/parallel_vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/parallel_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/base/std_cxx11/bind.h>

#include <numeric>
#include <algorithm>
#include <cmath>
#include <vector>

DEAL_II_NAMESPACE_OPEN



template <int spacedim>
template <typename InputVector, typename DoFHandlerType>
void
KellyErrorEstimator<1,spacedim>::
estimate (const Mapping<1,spacedim>                  &mapping,
          const DoFHandlerType                       &dof_handler,
          const Quadrature<0>                        &quadrature,
          const typename FunctionMap<spacedim>::type &neumann_bc,
          const InputVector                          &solution,
          Vector<float>                              &error,
          const ComponentMask                        &component_mask,
          const Function<spacedim>                   *coefficients,
          const unsigned int                          n_threads,
          const types::subdomain_id                   subdomain_id,
          const types::material_id                    material_id)
{
  // just pass on to the other function
  const std::vector<const InputVector *> solutions (1, &solution);
  std::vector<Vector<float>*>              errors (1, &error);
  estimate (mapping, dof_handler, quadrature, neumann_bc, solutions, errors,
            component_mask, coefficients, n_threads, subdomain_id, material_id);
}



template <int spacedim>
template <typename InputVector, typename DoFHandlerType>
void
KellyErrorEstimator<1,spacedim>::
estimate (const DoFHandlerType                       &dof_handler,
          const Quadrature<0>                        &quadrature,
          const typename FunctionMap<spacedim>::type &neumann_bc,
          const InputVector                          &solution,
          Vector<float>                              &error,
          const ComponentMask                        &component_mask,
          const Function<spacedim>                   *coefficients,
          const unsigned int                          n_threads,
          const types::subdomain_id                   subdomain_id,
          const types::material_id                    material_id)
{
  estimate(StaticMappingQ1<1,spacedim>::mapping, dof_handler, quadrature, neumann_bc, solution,
           error, component_mask, coefficients, n_threads, subdomain_id, material_id);
}



template <int spacedim>
template <typename InputVector, typename DoFHandlerType>
void
KellyErrorEstimator<1,spacedim>::
estimate (const DoFHandlerType                       &dof_handler,
          const Quadrature<0>                        &quadrature,
          const typename FunctionMap<spacedim>::type &neumann_bc,
          const std::vector<const InputVector *>     &solutions,
          std::vector<Vector<float>*>                &errors,
          const ComponentMask                        &component_mask,
          const Function<spacedim>                   *coefficients,
          const unsigned int                          n_threads,
          const types::subdomain_id                   subdomain_id,
          const types::material_id                    material_id)
{
  estimate(StaticMappingQ1<1,spacedim>::mapping, dof_handler, quadrature, neumann_bc, solutions,
           errors, component_mask, coefficients, n_threads, subdomain_id, material_id);
}



template <int spacedim>
template <typename InputVector, typename DoFHandlerType>
void
KellyErrorEstimator<1,spacedim>::
estimate (const Mapping<1,spacedim>                  &mapping,
          const DoFHandlerType                       &dof_handler,
          const hp::QCollection<0>                   &quadrature,
          const typename FunctionMap<spacedim>::type &neumann_bc,
          const InputVector                          &solution,
          Vector<float>                              &error,
          const ComponentMask                        &component_mask,
          const Function<spacedim>                   *coefficients,
          const unsigned int                          n_threads,
          const types::subdomain_id                   subdomain_id,
          const types::material_id                    material_id)
{
  // just pass on to the other function
  const std::vector<const InputVector *> solutions (1, &solution);
  std::vector<Vector<float>*>              errors (1, &error);
  estimate (mapping, dof_handler, quadrature, neumann_bc, solutions, errors,
            component_mask, coefficients, n_threads, subdomain_id, material_id);
}


template <int spacedim>
template <typename InputVector, typename DoFHandlerType>
void
KellyErrorEstimator<1,spacedim>::
estimate (const DoFHandlerType                       &dof_handler,
          const hp::QCollection<0>                   &quadrature,
          const typename FunctionMap<spacedim>::type &neumann_bc,
          const InputVector                          &solution,
          Vector<float>                              &error,
          const ComponentMask                        &component_mask,
          const Function<spacedim>                   *coefficients,
          const unsigned int                          n_threads,
          const types::subdomain_id                   subdomain_id,
          const types::material_id                    material_id)
{
  estimate(StaticMappingQ1<1,spacedim>::mapping, dof_handler, quadrature, neumann_bc, solution,
           error, component_mask, coefficients, n_threads, subdomain_id, material_id);
}



template <int spacedim>
template <typename InputVector, typename DoFHandlerType>
void
KellyErrorEstimator<1,spacedim>::
estimate (const DoFHandlerType                       &dof_handler,
          const hp::QCollection<0>                   &quadrature,
          const typename FunctionMap<spacedim>::type &neumann_bc,
          const std::vector<const InputVector *>     &solutions,
          std::vector<Vector<float>*>                &errors,
          const ComponentMask                        &component_mask,
          const Function<spacedim>                   *coefficients,
          const unsigned int                          n_threads,
          const types::subdomain_id                   subdomain_id,
          const types::material_id                    material_id)
{
  estimate(StaticMappingQ1<1,spacedim>::mapping, dof_handler, quadrature, neumann_bc, solutions,
           errors, component_mask, coefficients, n_threads, subdomain_id, material_id);
}




template <int spacedim>
template <typename InputVector, typename DoFHandlerType>
void KellyErrorEstimator<1,spacedim>::
estimate (const Mapping<1,spacedim>                  & /*mapping*/,
          const DoFHandlerType                       & /*dof_handler*/,
          const hp::QCollection<0> &,
          const typename FunctionMap<spacedim>::type & /*neumann_bc*/,
          const std::vector<const InputVector *>     & /*solutions*/,
          std::vector<Vector<float>*>                & /*errors*/,
          const ComponentMask                        & /*component_mask_*/,
          const Function<spacedim>                   * /*coefficient*/,
          const unsigned int,
          const types::subdomain_id                    /*subdomain_id*/,
          const types::material_id                     /*material_id*/)
{
  Assert (false, ExcInternalError());
}



template <int spacedim>
template <typename InputVector, typename DoFHandlerType>
void KellyErrorEstimator<1,spacedim>::
estimate (const Mapping<1,spacedim>                  &mapping,
          const DoFHandlerType                       &dof_handler,
          const Quadrature<0> &,
          const typename FunctionMap<spacedim>::type &neumann_bc,
          const std::vector<const InputVector *>     &solutions,
          std::vector<Vector<float>*>                &errors,
          const ComponentMask                        &component_mask,
          const Function<spacedim>                   *coefficient,
          const unsigned int,
          const types::subdomain_id                   subdomain_id_,
          const types::material_id                    material_id)
{
  typedef typename InputVector::value_type number;
#ifdef DEAL_II_WITH_P4EST
  if (dynamic_cast<const parallel::distributed::Triangulation<1,spacedim>*>
      (&dof_handler.get_triangulation())
      != 0)
    Assert ((subdomain_id_ == numbers::invalid_subdomain_id)
            ||
            (subdomain_id_ ==
             dynamic_cast<const parallel::distributed::Triangulation<1,spacedim>&>
             (dof_handler.get_triangulation()).locally_owned_subdomain()),
            ExcMessage ("For parallel distributed triangulations, the only "
                        "valid subdomain_id that can be passed here is the "
                        "one that corresponds to the locally owned subdomain id."));

  const types::subdomain_id subdomain_id
    = ((dynamic_cast<const parallel::distributed::Triangulation<1,spacedim>*>
        (&dof_handler.get_triangulation())
        != 0)
       ?
       dynamic_cast<const parallel::distributed::Triangulation<1,spacedim>&>
       (dof_handler.get_triangulation()).locally_owned_subdomain()
       :
       subdomain_id_);
#else
  const types::subdomain_id subdomain_id
    = subdomain_id_;
#endif

  const unsigned int n_components       = dof_handler.get_fe().n_components();
  const unsigned int n_solution_vectors = solutions.size();

  // sanity checks
  Assert (neumann_bc.find(numbers::internal_face_boundary_id) == neumann_bc.end(),
          ExcMessage("You are not allowed to list the special boundary "
                     "indicator for internal boundaries in your boundary "
                     "value map."));

  for (typename FunctionMap<spacedim>::type::const_iterator i=neumann_bc.begin();
       i!=neumann_bc.end(); ++i)
    Assert (i->second->n_components == n_components,
            ExcInvalidBoundaryFunction(i->first,
                                       i->second->n_components,
                                       n_components));

  Assert (component_mask.represents_n_components(n_components),
          ExcInvalidComponentMask());
  Assert (component_mask.n_selected_components(n_components) > 0,
          ExcInvalidComponentMask());

  Assert ((coefficient == 0) ||
          (coefficient->n_components == n_components) ||
          (coefficient->n_components == 1),
          ExcInvalidCoefficient());

  Assert (solutions.size() > 0,
          ExcNoSolutions());
  Assert (solutions.size() == errors.size(),
          ExcIncompatibleNumberOfElements(solutions.size(), errors.size()));
  for (unsigned int n=0; n<solutions.size(); ++n)
    Assert (solutions[n]->size() == dof_handler.n_dofs(),
            ExcDimensionMismatch(solutions[n]->size(),
                                 dof_handler.n_dofs()));

  Assert ((coefficient == 0) ||
          (coefficient->n_components == n_components) ||
          (coefficient->n_components == 1),
          ExcInvalidCoefficient());

  for (typename FunctionMap<spacedim>::type::const_iterator i=neumann_bc.begin();
       i!=neumann_bc.end(); ++i)
    Assert (i->second->n_components == n_components,
            ExcInvalidBoundaryFunction(i->first,
                                       i->second->n_components,
                                       n_components));

  // reserve one slot for each cell and set it to zero
  for (unsigned int n=0; n<n_solution_vectors; ++n)
    (*errors[n]).reinit (dof_handler.get_triangulation().n_active_cells());

  // fields to get the gradients on the present and the neighbor cell.
  //
  // for the neighbor gradient, we need several auxiliary fields, depending on
  // the way we get it (see below)
  std::vector<std::vector<std::vector<Tensor<1,spacedim,number> > > >
  gradients_here (n_solution_vectors,
                  std::vector<std::vector<Tensor<1,spacedim,number> > >(2, std::vector<Tensor<1,spacedim,number> >(n_components)));
  std::vector<std::vector<std::vector<Tensor<1,spacedim,number> > > >
  gradients_neighbor (gradients_here);
  std::vector<Vector<number> >
  grad_neighbor (n_solution_vectors, Vector<number>(n_components));

  // reserve some space for coefficient values at one point.  if there is no
  // coefficient, then we fill it by unity once and for all and don't set it
  // any more
  Vector<double> coefficient_values (n_components);
  if (coefficient == 0)
    for (unsigned int c=0; c<n_components; ++c)
      coefficient_values(c) = 1;

  const QTrapez<1> quadrature;
  const hp::QCollection<1> q_collection(quadrature);
  const QGauss<0> face_quadrature(1);
  const hp::QCollection<0> q_face_collection(face_quadrature);

  const hp::FECollection<1,spacedim> fe (dof_handler.get_fe());

  hp::MappingCollection<1,spacedim> mapping_collection;
  mapping_collection.push_back (mapping);

  hp::FEValues<1,spacedim> fe_values (mapping_collection, fe, q_collection,
                                      update_gradients);
  hp::FEFaceValues<1,spacedim> fe_face_values (/*mapping_collection,*/ fe, q_face_collection,
      update_normal_vectors);

  // loop over all cells and do something on the cells which we're told to
  // work on. note that the error indicator is only a sum over the two
  // contributions from the two vertices of each cell.
  for (typename DoFHandlerType::active_cell_iterator cell = dof_handler.begin_active();
       cell != dof_handler.end();
       ++cell)
    if (((subdomain_id == numbers::invalid_subdomain_id)
         ||
         (cell->subdomain_id() == subdomain_id))
        &&
        ((material_id == numbers::invalid_material_id)
         ||
         (cell->material_id() == material_id)))
      {
        for (unsigned int n=0; n<n_solution_vectors; ++n)
          (*errors[n])(cell->active_cell_index()) = 0;

        fe_values.reinit (cell);
        for (unsigned int s=0; s<n_solution_vectors; ++s)
          fe_values.get_present_fe_values()
          .get_function_gradients (*solutions[s], gradients_here[s]);

        // loop over the two points bounding this line. n==0 is left point,
        // n==1 is right point
        for (unsigned int n=0; n<2; ++n)
          {
            // find left or right active neighbor
            typename DoFHandlerType::cell_iterator neighbor = cell->neighbor(n);
            if (neighbor.state() == IteratorState::valid)
              while (neighbor->has_children())
                neighbor = neighbor->child(n==0 ? 1 : 0);

            fe_face_values.reinit (cell, n);
            Tensor<1,spacedim> normal =
              fe_face_values.get_present_fe_values().get_all_normal_vectors()[0];

            if (neighbor.state() == IteratorState::valid)
              {
                fe_values.reinit (neighbor);

                for (unsigned int s=0; s<n_solution_vectors; ++s)
                  fe_values.get_present_fe_values()
                  .get_function_gradients (*solutions[s],
                                           gradients_neighbor[s]);

                fe_face_values.reinit (neighbor, n==0 ? 1 : 0);
                Tensor<1,spacedim> neighbor_normal =
                  fe_face_values.get_present_fe_values().get_all_normal_vectors()[0];

                // extract the gradient in normal direction of all the components.
                for (unsigned int s=0; s<n_solution_vectors; ++s)
                  for (unsigned int c=0; c<n_components; ++c)
                    grad_neighbor[s](c)
                      = - (gradients_neighbor[s][n==0 ? 1 : 0][c]*neighbor_normal);
              }
            else if (neumann_bc.find(n) != neumann_bc.end())
              // if Neumann b.c., then fill the gradients field which will be
              // used later on.
              {
                if (n_components==1)
                  {
                    const double
                    v = neumann_bc.find(n)->second->value(cell->vertex(n));

                    for (unsigned int s=0; s<n_solution_vectors; ++s)
                      grad_neighbor[s](0) = v;
                  }
                else
                  {
                    Vector<double> v(n_components);
                    neumann_bc.find(n)->second->vector_value(cell->vertex(n), v);

                    for (unsigned int s=0; s<n_solution_vectors; ++s)
                      grad_neighbor[s] = v;
                  }
              }
            else
              // fill with zeroes.
              for (unsigned int s=0; s<n_solution_vectors; ++s)
                grad_neighbor[s] = 0;

            // if there is a coefficient, then evaluate it at the present
            // position. if there is none, reuse the preset values.
            if (coefficient != 0)
              {
                if (coefficient->n_components == 1)
                  {
                    const double c_value = coefficient->value (cell->vertex(n));
                    for (unsigned int c=0; c<n_components; ++c)
                      coefficient_values(c) = c_value;
                  }
                else
                  coefficient->vector_value(cell->vertex(n),
                                            coefficient_values);
              }


            for (unsigned int s=0; s<n_solution_vectors; ++s)
              for (unsigned int component=0; component<n_components; ++component)
                if (component_mask[component] == true)
                  {
                    // get gradient here
                    const number grad_here = gradients_here[s][n][component]
                                             * normal;

                    const number jump = ((grad_here - grad_neighbor[s](component)) *
                                         coefficient_values(component));
                    (*errors[s])(cell->active_cell_index()) += numbers::NumberTraits<number>::abs_square(jump) * cell->diameter();
                  }
          }

        for (unsigned int s=0; s<n_solution_vectors; ++s)
          (*errors[s])(cell->active_cell_index()) = std::sqrt((*errors[s])(cell->active_cell_index()));
      }
}


// explicit instantiations
#include "error_estimator_1d.inst"


DEAL_II_NAMESPACE_CLOSE
