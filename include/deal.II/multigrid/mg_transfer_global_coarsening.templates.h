// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_mg_transfer_global_coarsening_templates_h
#define dealii_mg_transfer_global_coarsening_templates_h

#include <deal.II/base/config.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/repartitioning_policy_tools.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/multigrid/mg_transfer_global_coarsening.h>

#include <limits>

DEAL_II_NAMESPACE_OPEN



namespace MGTransferGlobalCoarseningTools
{
  template <int dim, int spacedim>
  std::vector<std::shared_ptr<const Triangulation<dim, spacedim>>>
  create_geometric_coarsening_sequence(
    const Triangulation<dim, spacedim> &fine_triangulation_in)
  {
    std::vector<std::shared_ptr<const Triangulation<dim, spacedim>>>
      coarse_grid_triangulations(fine_triangulation_in.n_global_levels());

    coarse_grid_triangulations.back().reset(&fine_triangulation_in, [](auto *) {
      // empty deleter, since fine_triangulation_in is an external field
      // and its destructor is called somewhere else
    });

    // for a single level nothing has to be done
    if (fine_triangulation_in.n_global_levels() == 1)
      return coarse_grid_triangulations;

    Assert(
      (dynamic_cast<
         const parallel::fullydistributed::Triangulation<dim, spacedim> *>(
         &fine_triangulation_in) == nullptr),
      ExcMessage(
        "Triangulations of type parallel::fullydistributed::Triangulation are "
        "not supported by this function!"));

    const auto create_new_empty_triangulation =
      [&]() -> std::shared_ptr<Triangulation<dim, spacedim>> {
#ifdef DEAL_II_WITH_P4EST
      if (const auto fine_triangulation = dynamic_cast<
            const parallel::distributed::Triangulation<dim, spacedim> *>(
            &fine_triangulation_in))
        return std::make_shared<
          parallel::distributed::Triangulation<dim, spacedim>>(
          fine_triangulation->get_mpi_communicator());
      else
#endif
#ifdef DEAL_II_WITH_MPI
        if (const auto fine_triangulation = dynamic_cast<
              const parallel::shared::Triangulation<dim, spacedim> *>(
              &fine_triangulation_in))
        return std::make_shared<parallel::shared::Triangulation<dim, spacedim>>(
          fine_triangulation->get_mpi_communicator(),
          Triangulation<dim, spacedim>::none,
          fine_triangulation->with_artificial_cells());
      else
#endif
        return std::make_shared<Triangulation<dim, spacedim>>();
    };

    const unsigned int max_level = fine_triangulation_in.n_global_levels() - 1;

    // clear 'eliminate_unrefined_islands' from MeshSmoothing flags
    // to prevent unintentional refinement during coarsen_global()
    const auto mesh_smoothing =
      static_cast<typename Triangulation<dim, spacedim>::MeshSmoothing>(
        fine_triangulation_in.get_mesh_smoothing() &
        ~(Triangulation<dim, spacedim>::eliminate_unrefined_islands));

    // create coarse meshes
    for (unsigned int l = max_level; l > 0; --l)
      {
        // copy triangulation
        auto new_tria = create_new_empty_triangulation();
        new_tria->copy_triangulation(*coarse_grid_triangulations[l]);
        new_tria->set_mesh_smoothing(mesh_smoothing);

        // coarsen mesh
        new_tria->coarsen_global();

        // save mesh
        coarse_grid_triangulations[l - 1] = new_tria;
      }

    AssertDimension(coarse_grid_triangulations[0]->n_global_levels(), 1);

    return coarse_grid_triangulations;
  }



  template <int dim, int spacedim>
  std::vector<std::shared_ptr<const Triangulation<dim, spacedim>>>
  create_geometric_coarsening_sequence(
    Triangulation<dim, spacedim>                         &fine_triangulation_in,
    const RepartitioningPolicyTools::Base<dim, spacedim> &policy,
    const bool keep_fine_triangulation,
    const bool repartition_fine_triangulation)
  {
    std::vector<std::shared_ptr<const Triangulation<dim, spacedim>>>
      coarse_grid_triangulations(fine_triangulation_in.n_global_levels());

#ifndef DEAL_II_WITH_P4EST
    DEAL_II_NOT_IMPLEMENTED();
    (void)policy;
    (void)keep_fine_triangulation;
    (void)repartition_fine_triangulation;
#else
    const auto fine_triangulation =
      dynamic_cast<parallel::distributed::Triangulation<dim, spacedim> *>(
        &fine_triangulation_in);

    Assert(fine_triangulation, ExcNotImplemented());

    const auto comm = fine_triangulation->get_mpi_communicator();

    if (keep_fine_triangulation == true &&
        repartition_fine_triangulation == false)
      {
        coarse_grid_triangulations.back().reset(&fine_triangulation_in,
                                                [](auto *) {
                                                  // empty deleter, since
                                                  // fine_triangulation_in is an
                                                  // external field and its
                                                  // destructor is called
                                                  // somewhere else
                                                });
      }
    else
      {
        // create triangulation description
        const auto construction_data =
          repartition_fine_triangulation ?
            TriangulationDescription::Utilities::
              create_description_from_triangulation(
                *fine_triangulation, policy.partition(*fine_triangulation)) :
            TriangulationDescription::Utilities::
              create_description_from_triangulation(*fine_triangulation, comm);

        // create new triangulation
        const auto new_fine_triangulation = std::make_shared<
          parallel::fullydistributed::Triangulation<dim, spacedim>>(comm);

        for (const auto i : fine_triangulation->get_manifold_ids())
          if (i != numbers::flat_manifold_id)
            new_fine_triangulation->set_manifold(
              i, fine_triangulation->get_manifold(i));

        new_fine_triangulation->create_triangulation(construction_data);

        // save mesh
        coarse_grid_triangulations.back() = new_fine_triangulation;
      }

    // for a single level nothing has to be done
    if (fine_triangulation_in.n_global_levels() == 1)
      return coarse_grid_triangulations;

    parallel::distributed::Triangulation<dim, spacedim> temp_triangulation(
      comm);

    if (keep_fine_triangulation == true)
      temp_triangulation.copy_triangulation(*fine_triangulation);

    auto *temp_triangulation_ptr =
      keep_fine_triangulation ? &temp_triangulation : fine_triangulation;

    // clear 'eliminate_unrefined_islands' from MeshSmoothing flags
    // to prevent unintentional refinement during coarsen_global()
    const auto mesh_smoothing =
      static_cast<typename Triangulation<dim, spacedim>::MeshSmoothing>(
        temp_triangulation_ptr->get_mesh_smoothing() &
        ~(Triangulation<dim, spacedim>::eliminate_unrefined_islands));
    temp_triangulation_ptr->set_mesh_smoothing(mesh_smoothing);

    const unsigned int max_level = fine_triangulation->n_global_levels() - 1;

    // create coarse meshes
    for (unsigned int l = max_level; l > 0; --l)
      {
        // coarsen mesh
        temp_triangulation_ptr->coarsen_global();

        // create triangulation description
        const auto construction_data = TriangulationDescription::Utilities::
          create_description_from_triangulation(
            *temp_triangulation_ptr, policy.partition(*temp_triangulation_ptr));

        // create new triangulation
        const auto level_triangulation = std::make_shared<
          parallel::fullydistributed::Triangulation<dim, spacedim>>(comm);

        for (const auto i : fine_triangulation->get_manifold_ids())
          if (i != numbers::flat_manifold_id)
            level_triangulation->set_manifold(
              i, fine_triangulation->get_manifold(i));

        level_triangulation->create_triangulation(construction_data);

        // save mesh
        coarse_grid_triangulations[l - 1] = level_triangulation;
      }

    // recover MeshSmoothing flags in case we used the fine_triangulation
    // to build the sequence
    if (keep_fine_triangulation == false)
      fine_triangulation->set_mesh_smoothing(
        coarse_grid_triangulations.back()->get_mesh_smoothing());
#endif

    AssertDimension(coarse_grid_triangulations[0]->n_global_levels(), 1);

    return coarse_grid_triangulations;
  }



  template <int dim, int spacedim>
  std::vector<std::shared_ptr<const Triangulation<dim, spacedim>>>
  create_geometric_coarsening_sequence(
    const Triangulation<dim, spacedim>                   &fine_triangulation_in,
    const RepartitioningPolicyTools::Base<dim, spacedim> &policy,
    const bool repartition_fine_triangulation)
  {
    // remove const and convert it to flag
    return create_geometric_coarsening_sequence(
      const_cast<Triangulation<dim, spacedim> &>(fine_triangulation_in),
      policy,
      true,
      repartition_fine_triangulation);
  }

} // namespace MGTransferGlobalCoarseningTools


DEAL_II_NAMESPACE_CLOSE

#endif
