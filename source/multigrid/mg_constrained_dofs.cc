// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/dofs/dof_tools.h>

#include <deal.II/grid/manifold.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_tools.h>

#include <set>


DEAL_II_NAMESPACE_OPEN


template <int dim, int spacedim>
void
MGConstrainedDoFs::initialize(
  const DoFHandler<dim, spacedim> &dof,
  const MGLevelObject<IndexSet>   &level_relevant_dofs)
{
  boundary_indices.clear();
  refinement_edge_indices.clear();
  level_constraints.clear();
  user_constraints.clear();

  const unsigned int nlevels   = dof.get_triangulation().n_global_levels();
  const unsigned int min_level = level_relevant_dofs.min_level();
  const unsigned int max_level = (level_relevant_dofs.max_level() == 0) ?
                                   nlevels - 1 :
                                   level_relevant_dofs.max_level();
  const bool         use_provided_level_relevant_dofs =
    (level_relevant_dofs.max_level() > 0);

  // At this point level_constraint and refinement_edge_indices are empty.
  refinement_edge_indices.resize(nlevels);
  level_constraints.resize(nlevels);
  user_constraints.resize(nlevels);
  for (unsigned int l = min_level; l <= max_level; ++l)
    {
      if (use_provided_level_relevant_dofs)
        {
          level_constraints[l].reinit(dof.locally_owned_mg_dofs(l),
                                      level_relevant_dofs[l]);
        }
      else
        {
          const IndexSet relevant_dofs =
            DoFTools::extract_locally_relevant_level_dofs(dof, l);
          level_constraints[l].reinit(dof.locally_owned_mg_dofs(l),
                                      relevant_dofs);
        }
    }

  // TODO: currently we only consider very basic periodic constraints
  const IdentityMatrix transformation(dof.get_fe().n_dofs_per_face());
  const ComponentMask  component_mask;
  const double         periodicity_factor = 1.0;

  for (const auto &[first_cell, second_cell] :
       dof.get_triangulation().get_periodic_face_map())
    {
      // only consider non-artificial cells
      if (first_cell.first->is_artificial_on_level())
        continue;
      if (second_cell.first.first->is_artificial_on_level())
        continue;

      // consider cell pairs with the same level
      if (first_cell.first->level() != second_cell.first.first->level())
        continue;

      DoFTools::internal::set_periodicity_constraints(
        first_cell.first->as_dof_handler_level_iterator(dof)->face(
          first_cell.second),
        second_cell.first.first->as_dof_handler_level_iterator(dof)->face(
          second_cell.first.second),
        transformation,
        level_constraints[first_cell.first->level()],
        component_mask,
        second_cell.second,
        periodicity_factor,
        first_cell.first->level());
    }

  for (unsigned int l = min_level; l <= max_level; ++l)
    {
      level_constraints[l].close();

      // Initialize with empty IndexSet of correct size
      refinement_edge_indices[l] = IndexSet(dof.n_dofs(l));
    }

  MGTools::extract_inner_interface_dofs(dof, refinement_edge_indices);
}



template <int dim, int spacedim>
void
MGConstrainedDoFs::make_zero_boundary_constraints(
  const DoFHandler<dim, spacedim>    &dof,
  const std::set<types::boundary_id> &boundary_ids,
  const ComponentMask                &component_mask)
{
  // allocate an IndexSet for each global level. Contents will be
  // overwritten inside make_boundary_list.
  const unsigned int n_levels = dof.get_triangulation().n_global_levels();
  Assert(boundary_indices.empty() || boundary_indices.size() == n_levels,
         ExcInternalError());
  boundary_indices.resize(n_levels);

  MGTools::make_boundary_list(dof,
                              boundary_ids,
                              boundary_indices,
                              component_mask);
}



template <int dim, int spacedim>
void
MGConstrainedDoFs::add_boundary_indices(const DoFHandler<dim, spacedim> &dof,
                                        const unsigned int               level,
                                        const IndexSet &level_boundary_indices)
{
  const unsigned int n_levels = dof.get_triangulation().n_global_levels();
  if (boundary_indices.empty())
    {
      boundary_indices.resize(n_levels);
      for (unsigned int i = 0; i < n_levels; ++i)
        boundary_indices[i] = IndexSet(dof.n_dofs(i));
    }
  AssertDimension(boundary_indices.size(), n_levels);
  boundary_indices[level].add_indices(level_boundary_indices);
}



template <int dim, int spacedim>
void
MGConstrainedDoFs::make_no_normal_flux_constraints(
  const DoFHandler<dim, spacedim> &dof,
  const types::boundary_id         bid,
  const unsigned int               first_vector_component)
{
  // For a given boundary id, find which vector component is on the boundary
  // and set a zero boundary constraint for those degrees of freedom.
  const unsigned int n_components = dof.get_fe_collection().n_components();
  AssertIndexRange(first_vector_component + dim - 1, n_components);

  ComponentMask comp_mask(n_components, false);


  typename Triangulation<dim>::face_iterator
    face = dof.get_triangulation().begin_face(),
    endf = dof.get_triangulation().end_face();
  for (; face != endf; ++face)
    if (face->at_boundary() && face->boundary_id() == bid)
      for (unsigned int d = 0; d < dim; ++d)
        {
          Tensor<1, dim, double> unit_vec;
          unit_vec[d] = 1.0;

          const Tensor<1, dim> normal_vec =
            face->get_manifold().normal_vector(face, face->center());

          if (std::abs(std::abs(unit_vec * normal_vec) - 1.0) < 1e-10)
            comp_mask.set(d + first_vector_component, true);
          else
            Assert(
              std::abs(unit_vec * normal_vec) < 1e-10,
              ExcMessage(
                "We can currently only support no normal flux conditions "
                "for a specific boundary id if all faces are normal to the "
                "x, y, or z axis."));
        }

  Assert(comp_mask.n_selected_components() == 1,
         ExcMessage(
           "We can currently only support no normal flux conditions "
           "for a specific boundary id if all faces are facing in the "
           "same direction, i.e., a boundary normal to the x-axis must "
           "have a different boundary id than a boundary normal to the "
           "y- or z-axis and so on. If the mesh here was produced using "
           "GridGenerator::..., setting colorize=true during mesh generation "
           "and calling make_no_normal_flux_constraints() for each no normal "
           "flux boundary will fulfill the condition."));

  this->make_zero_boundary_constraints(dof, {bid}, comp_mask);
}



void
MGConstrainedDoFs::add_user_constraints(
  const unsigned int               level,
  const AffineConstraints<double> &constraints_on_level)
{
  AssertIndexRange(level, user_constraints.size());

  // Get the relevant DoFs from level_constraints if
  // the user constraint matrix has not been initialized
  if (user_constraints[level].get_local_lines().size() == 0)
    user_constraints[level].reinit(
      level_constraints[level].get_locally_owned_indices(),
      level_constraints[level].get_local_lines());

  user_constraints[level].merge(
    constraints_on_level,
    AffineConstraints<double>::MergeConflictBehavior::right_object_wins);
  user_constraints[level].close();
}



void
MGConstrainedDoFs::clear_user_constraints()
{
  for (auto &constraint : user_constraints)
    constraint.clear();
}



void
MGConstrainedDoFs::clear()
{
  boundary_indices.clear();
  refinement_edge_indices.clear();
  user_constraints.clear();
}



template <typename Number>
void
MGConstrainedDoFs::merge_constraints(AffineConstraints<Number> &constraints,
                                     const unsigned int         level,
                                     const bool add_boundary_indices,
                                     const bool add_refinement_edge_indices,
                                     const bool add_level_constraints,
                                     const bool add_user_constraints) const
{
  constraints.clear();

  // determine local lines
  IndexSet index_set(this->get_refinement_edge_indices(level).size());

  if (add_boundary_indices && this->have_boundary_indices())
    index_set.add_indices(this->get_boundary_indices(level));

  if (add_refinement_edge_indices)
    index_set.add_indices(this->get_refinement_edge_indices(level));

  if (add_level_constraints)
    index_set.add_indices(this->get_level_constraints(level).get_local_lines());

  if (add_user_constraints)
    index_set.add_indices(
      this->get_user_constraint_matrix(level).get_local_lines());

  constraints.reinit(level_constraints[level].get_locally_owned_indices(),
                     index_set);

  // merge constraints
  if (add_boundary_indices && this->have_boundary_indices())
    for (const auto i : this->get_boundary_indices(level))
      constraints.constrain_dof_to_zero(i);

  if (add_refinement_edge_indices)
    for (const auto i : this->get_refinement_edge_indices(level))
      constraints.constrain_dof_to_zero(i);

  if (add_level_constraints)
    constraints.merge(this->get_level_constraints(level),
                      AffineConstraints<Number>::left_object_wins,
                      true);

  if (add_user_constraints)
    constraints.merge(this->get_user_constraint_matrix(level),
                      AffineConstraints<Number>::left_object_wins,
                      true);

  // finalize setup
  constraints.close();
}

#include "multigrid/mg_constrained_dofs.inst"


DEAL_II_NAMESPACE_CLOSE
