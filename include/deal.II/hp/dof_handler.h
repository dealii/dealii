// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2021 by the deal.II authors
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

#ifndef dealii_hp_dof_handler_h
#define dealii_hp_dof_handler_h

#include <deal.II/dofs/dof_handler.h>

DEAL_II_NAMESPACE_OPEN

namespace hp
{
  /**
   * Manage the distribution and numbering of the degrees of freedom for hp-
   * FEM algorithms. This class satisfies the
   * @ref ConceptMeshType "MeshType concept"
   * requirements.
   *
   * The purpose of this class is to allow for an enumeration of degrees of
   * freedom in the same way as the ::DoFHandler class, but it allows to use a
   * different finite element on every cell. To this end, one assigns an
   * <code>active_fe_index</code> to every cell that indicates which element
   * within a collection of finite elements (represented by an object of type
   * hp::FECollection) is the one that lives on this cell. The class then
   * enumerates the degree of freedom associated with these finite elements on
   * each cell of a triangulation and, if possible, identifies degrees of
   * freedom at the interfaces of cells if they match. If neighboring cells
   * have degrees of freedom along the common interface that do not immediate
   * match (for example, if you have $Q_2$ and $Q_3$ elements meeting at a
   * common face), then one needs to compute constraints to ensure that the
   * resulting finite element space on the mesh remains conforming.
   *
   * The whole process of working with objects of this type is explained in
   * step-27. Many of the algorithms this class implements are described in
   * the
   * @ref hp_paper "hp-paper".
   *
   *
   * <h3>Active FE indices and their behavior under mesh refinement</h3>
   *
   * The typical workflow for using this class is to create a mesh, assign an
   * active FE index to every active cell, calls
   * hp::DoFHandler::distribute_dofs(), and then assemble a linear system and
   * solve a problem on this finite element space. However, one can skip
   * assigning active FE indices upon mesh refinement in certain
   * circumstances. In particular, the following rules apply:
   * - Upon mesh refinement, child cells inherit the active FE index of
   *   the parent.
   * - When coarsening cells, the (now active) parent cell will be assigned
   *   an active FE index that is determined from its (no longer active)
   *   children, following the FiniteElementDomination logic: Out of the set of
   *   elements previously assigned to the former children, we choose the one
   *   dominated by all children for the parent cell. If none was found, we pick
   *   the most dominant element in the whole collection that is dominated by
   *   all former children. See hp::FECollection::find_dominated_fe_extended()
   *   for further information on this topic.
   *
   * @note Finite elements need to be assigned to each cell by either calling
   * set_fe() or distribute_dofs() first to make this functionality available.
   *
   *
   * <h3>Active FE indices and parallel meshes</h3>
   *
   * When this class is used with either a parallel::shared::Triangulation
   * or a parallel::distributed::Triangulation, you can only set active
   * FE indices on cells that are locally owned,
   * using a call such as <code>cell-@>set_active_fe_index(...)</code>.
   * On the other hand, setting the active FE index on ghost
   * or artificial cells is not allowed.
   *
   * Ghost cells do acquire the information what element
   * is active on them, however: whenever
   * you call hp::DoFHandler::distribute_dofs(), all processors that
   * participate in the parallel mesh exchange information in such a way
   * that the active FE index on ghost cells equals the active FE index
   * that was set on that processor that owned that particular ghost cell.
   * Consequently, one can <i>query</i> the @p active_fe_index on ghost
   * cells, just not set it by hand.
   *
   * On artificial cells, no information is available about the
   * @p active_fe_index used there. That's because we don't even know
   * whether these cells exist at all, and even if they did, the
   * current processor does not know anything specific about them.
   * See
   * @ref GlossArtificialCell "the glossary entry on artificial cells"
   * for more information.
   *
   * During refinement and coarsening, information about the @p active_fe_index
   * of each cell will be automatically transferred.
   *
   * However, using a parallel::distributed::Triangulation with an
   * hp::DoFHandler requires additional attention during serialization, since no
   * information on active FE indices will be automatically transferred. This
   * has to be done manually using the
   * prepare_for_serialization_of_active_fe_indices() and
   * deserialize_active_fe_indices() functions. The former has to be called
   * before parallel::distributed::Triangulation::save() is invoked, and the
   * latter needs to be run after parallel::distributed::Triangulation::load().
   * If further data will be attached to the triangulation via the
   * parallel::distributed::CellDataTransfer,
   * parallel::distributed::SolutionTransfer, or Particles::ParticleHandler
   * classes, all corresponding preparation and deserialization function calls
   * need to happen in the same order. Consult the documentation of
   * parallel::distributed::SolutionTransfer for more information.
   *
   *
   * @ingroup dofs
   *
   * @deprecated The basic dealii::DoFHandler is capable of hp-adaptation now.
   */
  template <int dim, int spacedim = dim>
  using DoFHandler DEAL_II_DEPRECATED = dealii::DoFHandler<dim, spacedim>;
} // namespace hp

DEAL_II_NAMESPACE_CLOSE

#endif
