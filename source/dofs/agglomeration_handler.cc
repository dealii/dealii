// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2001 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



#include <deal.II/dofs/agglomeration_handler.h>

#include <deal.II/lac/sparsity_tools.h>

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim>
AgglomerationHandler<dim, spacedim>::AgglomerationHandler(
  const GridTools::Cache<dim, spacedim> &cache_tria)
  : cached_tria(std::make_unique<GridTools::Cache<dim, spacedim>>(
      cache_tria.get_triangulation(),
      cache_tria.get_mapping()))
  , communicator(cache_tria.get_triangulation().get_communicator())
{
  Assert(dim == spacedim,
         ExcNotImplemented("Not available with codimension greater than 0"));
  Assert(dim == 2 || dim == 3, ExcImpossibleInDim(1));
  Assert((dynamic_cast<const parallel::shared::Triangulation<dim, spacedim> *>(
            &cached_tria->get_triangulation()) == nullptr),
         ExcNotImplemented());
  Assert(cached_tria->get_triangulation().n_active_cells() > 0,
         ExcMessage(
           "The triangulation must not be empty upon calling this function."));

  n_agglomerations = 0;
  hybrid_mesh      = false;
  initialize_agglomeration_data(cached_tria);
}



template <int dim, int spacedim>
typename AgglomerationHandler<dim, spacedim>::agglomeration_iterator
AgglomerationHandler<dim, spacedim>::define_agglomerate(
  const AgglomerationContainer &cells)
{
  Assert(cells.size() > 0, ExcMessage("No cells to be agglomerated."));

  if (cells.size() == 1)
    hybrid_mesh = true; // mesh is made also by classical cells

  // First index drives the selection of the master cell. After that, store the
  // master cell.
  const types::global_cell_index global_master_idx =
    cells[0]->global_active_cell_index();
  const types::global_cell_index master_idx = cells[0]->active_cell_index();
  master_cells_container.push_back(cells[0]);
  master_slave_relationships[global_master_idx] = -1;

  const typename DoFHandler<dim>::active_cell_iterator cell_dh =
    cells[0]->as_dof_handler_iterator(agglo_dh);
  cell_dh->set_active_fe_index(CellAgglomerationType::master);

  // Store slave cells and save the relationship with the parent
  std::vector<typename Triangulation<dim, spacedim>::active_cell_iterator>
    slaves;
  slaves.reserve(cells.size() - 1);
  // exclude first cell since it's the master cell
  for (auto it = ++cells.begin(); it != cells.end(); ++it)
    {
      slaves.push_back(*it);
      master_slave_relationships[(*it)->global_active_cell_index()] =
        global_master_idx; // mark each slave
      master_slave_relationships_iterators[(*it)->active_cell_index()] =
        cells[0];

      const typename DoFHandler<dim>::active_cell_iterator cell =
        (*it)->as_dof_handler_iterator(agglo_dh);
      cell->set_active_fe_index(CellAgglomerationType::slave); // slave cell

      // If we have a p::d::T, check that all cells are in the same subdomain.
      // If serial, just check that the subdomain_id is invalid.
      Assert(((*it)->subdomain_id() == tria->locally_owned_subdomain() ||
              tria->locally_owned_subdomain() == numbers::invalid_subdomain_id),
             ExcInternalError());
    }

  master_slave_relationships_iterators[master_idx] =
    cells[0]; // set iterator to master cell

  // Store the slaves of each master
  master2slaves[master_idx] = slaves;
  // Save to which polygon this agglomerate correspond
  master2polygon[master_idx] = n_agglomerations;

  ++n_agglomerations; // an agglomeration has been performed, record it

  create_bounding_box(cells); // fill the vector of bboxes

  // Finally, return a polygonal iterator to the polytope just constructed.
  return {cells[0], this};
}



template <int dim, int spacedim>
void
AgglomerationHandler<dim, spacedim>::initialize_fe_values(
  const Quadrature<dim>     &cell_quadrature,
  const UpdateFlags         &flags,
  const Quadrature<dim - 1> &face_quadrature,
  const UpdateFlags         &face_flags)
{
  agglomeration_quad       = cell_quadrature;
  agglomeration_flags      = flags;
  agglomeration_face_quad  = face_quadrature;
  agglomeration_face_flags = face_flags | internal_agglomeration_face_flags;


  no_values = std::make_unique<FEValues<dim>>(
    *mapping,
    dummy_fe,
    agglomeration_quad,
    update_quadrature_points | update_JxW_values); // needed only for quadrature
  no_face_values = std::make_unique<FEFaceValues<dim>>(
    *mapping,
    dummy_fe,
    agglomeration_face_quad,
    update_quadrature_points | update_JxW_values |
      update_normal_vectors); // needed only for quadrature (faces)
}



template <int dim, int spacedim>
void
AgglomerationHandler<dim, spacedim>::initialize_agglomeration_data(
  const std::unique_ptr<GridTools::Cache<dim, spacedim>> &cache_tria)
{
  tria    = &cache_tria->get_triangulation();
  mapping = &cache_tria->get_mapping();

  agglo_dh.reinit(*tria);

  if (const auto parallel_tria = dynamic_cast<
        const dealii::parallel::TriangulationBase<dim, spacedim> *>(&*tria))
    {
      const std::weak_ptr<const Utilities::MPI::Partitioner> cells_partitioner =
        parallel_tria->global_active_cell_index_partitioner();
      master_slave_relationships.reinit(
        cells_partitioner.lock()->locally_owned_range(), communicator);
    }
  else
    {
      master_slave_relationships.reinit(tria->n_active_cells(), MPI_COMM_SELF);
    }

  polytope_cache.clear();
  bboxes.clear();

  // First, update the pointer
  cached_tria = std::make_unique<GridTools::Cache<dim, spacedim>>(
    cache_tria->get_triangulation(), cache_tria->get_mapping());


  n_agglomerations = 0;
}



template <int dim, int spacedim>
void
AgglomerationHandler<dim, spacedim>::distribute_agglomerated_dofs(
  const FiniteElement<dim> &fe_space)
{
  if (dynamic_cast<const FE_DGQ<dim> *>(&fe_space))
    fe = std::make_unique<FE_DGQ<dim>>(fe_space.degree);
  else if (dynamic_cast<const FE_SimplexDGP<dim> *>(&fe_space))
    fe = std::make_unique<FE_SimplexDGP<dim>>(fe_space.degree);
  else
    AssertThrow(
      false,
      ExcNotImplemented(
        "Currently, this interface supports only DGQ and DGP bases."));

  box_mapping = std::make_unique<MappingBox<dim>>(
    bboxes,
    master2polygon); // construct bounding box mapping

  if (hybrid_mesh)
    {
      // the mesh is composed by standard and agglomerate cells. initialize
      // classes needed for standard cells in order to treat that finite
      // element space as defined on a standard shape and not on the
      // BoundingBox.
      standard_scratch =
        std::make_unique<ScratchData>(*mapping,
                                      *fe,
                                      QGauss<dim>(2 * fe_space.degree + 2),
                                      internal_agglomeration_flags);
    }


  fe_collection.push_back(*fe); // master
  fe_collection.push_back(
    FE_Nothing<dim, spacedim>(fe->reference_cell())); // slave

  initialize_hp_structure();

  // in case the tria is distributed, communicate ghost information with
  // neighboring ranks
  const bool needs_ghost_info =
    dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(&*tria) !=
    nullptr;
  if (needs_ghost_info)
    setup_ghost_polytopes();

  setup_connectivity_of_agglomeration();

  if (needs_ghost_info)
    exchange_interface_values();
}



template <int dim, int spacedim>
void
AgglomerationHandler<dim, spacedim>::create_bounding_box(
  const AgglomerationContainer &polytope)
{
  Assert(n_agglomerations > 0,
         ExcMessage("No agglomeration has been performed."));
  Assert(dim > 1, ExcNotImplemented());

  std::vector<Point<spacedim>> pts; // store all the vertices
  for (const auto &cell : polytope)
    for (const auto i : cell->vertex_indices())
      pts.push_back(cell->vertex(i));

  bboxes.emplace_back(pts);
}



template <int dim, int spacedim>
void
AgglomerationHandler<dim, spacedim>::setup_connectivity_of_agglomeration()
{
  Assert(master_cells_container.size() > 0,
         ExcMessage("No agglomeration has been performed."));
  Assert(
    agglo_dh.n_dofs() > 0,
    ExcMessage(
      "The DoFHandler associated to the agglomeration has not been initialized."
      "It's likely that you forgot to distribute the DoFs. You may want"
      "to check if a call to `initialize_hp_structure()` has been done."));

  number_of_agglomerated_faces.resize(master2polygon.size(), 0);
  for (const auto &cell : master_cells_container)
    {
      internal::AgglomerationHandlerImplementation<dim, spacedim>::
        setup_master_neighbor_connectivity(cell, *this);
    }

  if (Utilities::MPI::job_supports_mpi())
    {
      // communicate the number of faces
      recv_n_faces = Utilities::MPI::some_to_some(communicator, local_n_faces);

      // send information about boundaries and neighboring polytopes id
      recv_bdary_info =
        Utilities::MPI::some_to_some(communicator, local_bdary_info);

      recv_ghosted_master_id =
        Utilities::MPI::some_to_some(communicator, local_ghosted_master_id);
    }
}



template <int dim, int spacedim>
void
AgglomerationHandler<dim, spacedim>::exchange_interface_values()
{
  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  for (const auto &polytope : polytope_iterators())
    {
      if (polytope->is_locally_owned())
        {
          const unsigned int n_faces = polytope->n_faces();
          for (unsigned int f = 0; f < n_faces; ++f)
            {
              if (!polytope->at_boundary(f))
                {
                  const auto &neigh_polytope = polytope->neighbor(f);
                  if (!neigh_polytope->is_locally_owned())
                    {
                      // Neighboring polytope is ghosted.

                      // Compute shape functions at the interface
                      const auto &current_fe = reinit(polytope, f);

                      std::vector<Point<spacedim>> qpoints_to_send =
                        current_fe.get_quadrature_points();

                      const std::vector<double> &jxws_to_send =
                        current_fe.get_JxW_values();

                      const std::vector<Tensor<1, spacedim>> &normals_to_send =
                        current_fe.get_normal_vectors();


                      const types::subdomain_id neigh_rank =
                        neigh_polytope->subdomain_id();

                      std::pair<CellId, unsigned int> cell_and_face{
                        polytope->id(), f};
                      // Prepare data to send
                      local_qpoints[neigh_rank].emplace(cell_and_face,
                                                        qpoints_to_send);

                      local_jxws[neigh_rank].emplace(cell_and_face,
                                                     jxws_to_send);

                      local_normals[neigh_rank].emplace(cell_and_face,
                                                        normals_to_send);


                      const unsigned int n_qpoints = qpoints_to_send.size();

                      // TODO: check `agglomeration_flags` before computing
                      // values and gradients.
                      std::vector<std::vector<double>> values_per_qpoints(
                        dofs_per_cell);

                      std::vector<std::vector<Tensor<1, spacedim>>>
                        gradients_per_qpoints(dofs_per_cell);

                      for (unsigned int i = 0; i < dofs_per_cell; ++i)
                        {
                          values_per_qpoints[i].resize(n_qpoints);
                          gradients_per_qpoints[i].resize(n_qpoints);
                          for (unsigned int q = 0; q < n_qpoints; ++q)
                            {
                              values_per_qpoints[i][q] =
                                current_fe.shape_value(i, q);
                              gradients_per_qpoints[i][q] =
                                current_fe.shape_grad(i, q);
                            }
                        }

                      local_values[neigh_rank].emplace(cell_and_face,
                                                       values_per_qpoints);
                      local_gradients[neigh_rank].emplace(
                        cell_and_face, gradients_per_qpoints);
                    }
                }
            }
        }
    }

  // Finally, exchange with neighboring ranks
  recv_qpoints   = Utilities::MPI::some_to_some(communicator, local_qpoints);
  recv_jxws      = Utilities::MPI::some_to_some(communicator, local_jxws);
  recv_normals   = Utilities::MPI::some_to_some(communicator, local_normals);
  recv_values    = Utilities::MPI::some_to_some(communicator, local_values);
  recv_gradients = Utilities::MPI::some_to_some(communicator, local_gradients);
}



template <int dim, int spacedim>
Quadrature<dim>
AgglomerationHandler<dim, spacedim>::agglomerated_quadrature(
  const typename AgglomerationHandler<dim, spacedim>::AgglomerationContainer
    &cells,
  const typename Triangulation<dim, spacedim>::active_cell_iterator
    &master_cell) const
{
  Assert(is_master_cell(master_cell),
         ExcMessage("This must be a master cell."));

  std::vector<Point<dim>> vec_pts;
  std::vector<double>     vec_JxWs;
  for (const auto &dummy_cell : cells)
    {
      no_values->reinit(dummy_cell);
      auto        q_points = no_values->get_quadrature_points(); // real qpoints
      const auto &JxWs     = no_values->get_JxW_values();

      std::transform(q_points.begin(),
                     q_points.end(),
                     std::back_inserter(vec_pts),
                     [&](const Point<spacedim> &p) { return p; });
      std::transform(JxWs.begin(),
                     JxWs.end(),
                     std::back_inserter(vec_JxWs),
                     [&](const double w) { return w; });
    }

  // Map back each point in real space by using the map associated to the
  // bounding box.
  std::vector<Point<dim>> unit_points(vec_pts.size());
  const auto             &bbox =
    bboxes[master2polygon.at(master_cell->active_cell_index())];
  unit_points.reserve(vec_pts.size());

  for (unsigned int i = 0; i < vec_pts.size(); i++)
    unit_points[i] = bbox.real_to_unit(vec_pts[i]);

  return Quadrature<dim>(unit_points, vec_JxWs);
}



template <int dim, int spacedim>
void
AgglomerationHandler<dim, spacedim>::initialize_hp_structure()
{
  Assert(agglo_dh.get_triangulation().n_cells() > 0,
         ExcMessage(
           "Triangulation must not be empty upon calling this function."));
  Assert(n_agglomerations > 0,
         ExcMessage("No agglomeration has been performed."));

  agglo_dh.distribute_dofs(fe_collection);
}



template <int dim, int spacedim>
const FEValues<dim, spacedim> &
AgglomerationHandler<dim, spacedim>::reinit(
  const AgglomerationIterator<dim, spacedim> &polytope) const
{
  const auto &deal_cell = polytope->as_dof_handler_iterator(agglo_dh);

  const auto &agglo_cells = polytope->get_agglomerate();

  Quadrature<dim> agglo_quad = agglomerated_quadrature(agglo_cells, deal_cell);

  agglomerated_scratch = std::make_unique<ScratchData>(*box_mapping,
                                                       fe_collection[0],
                                                       agglo_quad,
                                                       agglomeration_flags);
  return agglomerated_scratch->reinit(deal_cell);
}



template <int dim, int spacedim>
const FEValuesBase<dim, spacedim> &
AgglomerationHandler<dim, spacedim>::reinit_master(
  const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
  const unsigned int                                              face_index,
  std::unique_ptr<NonMatching::FEImmersedSurfaceValues<spacedim>>
    &agglo_isv_ptr) const
{
  return internal::AgglomerationHandlerImplementation<dim, spacedim>::
    reinit_master(cell, face_index, agglo_isv_ptr, *this);
}



template <int dim, int spacedim>
const FEValuesBase<dim, spacedim> &
AgglomerationHandler<dim, spacedim>::reinit(
  const AgglomerationIterator<dim, spacedim> &polytope,
  const unsigned int                          face_index) const
{
  const auto &deal_cell = polytope->as_dof_handler_iterator(agglo_dh);
  Assert(is_master_cell(deal_cell), ExcMessage("This should be true."));

  return internal::AgglomerationHandlerImplementation<dim, spacedim>::
    reinit_master(deal_cell, face_index, agglomerated_isv_bdary, *this);
}



template <int dim, int spacedim>
std::pair<const FEValuesBase<dim, spacedim> &,
          const FEValuesBase<dim, spacedim> &>
AgglomerationHandler<dim, spacedim>::reinit_interface(
  const AgglomerationIterator<dim, spacedim> &polytope_in,
  const AgglomerationIterator<dim, spacedim> &neigh_polytope,
  const unsigned int                          local_in,
  const unsigned int                          local_neigh) const
{
  // If current and neighboring polytopes are both locally owned, then compute
  // the jump in the classical way without needing information about ghosted
  // entities.
  if (polytope_in->is_locally_owned() && neigh_polytope->is_locally_owned())
    {
      const auto &cell_in = polytope_in->as_dof_handler_iterator(agglo_dh);
      const auto &neigh_cell =
        neigh_polytope->as_dof_handler_iterator(agglo_dh);

      const auto &fe_in =
        internal::AgglomerationHandlerImplementation<dim, spacedim>::
          reinit_master(cell_in, local_in, agglomerated_isv, *this);
      const auto &fe_out =
        internal::AgglomerationHandlerImplementation<dim, spacedim>::
          reinit_master(neigh_cell, local_neigh, agglomerated_isv_neigh, *this);
      std::pair<const FEValuesBase<dim, spacedim> &,
                const FEValuesBase<dim, spacedim> &>
        my_p(fe_in, fe_out);

      return my_p;
    }
  else
    {
      Assert((polytope_in->is_locally_owned() &&
              !neigh_polytope->is_locally_owned()),
             ExcInternalError());

      const auto &cell = polytope_in->as_dof_handler_iterator(agglo_dh);
      const auto &bbox = bboxes[master2polygon.at(cell->active_cell_index())];
      // const double bbox_measure = bbox.volume();

      const unsigned int neigh_rank = neigh_polytope->subdomain_id();
      const CellId      &neigh_id   = neigh_polytope->id();

      // Retrieve qpoints,JxWs, normals sent previously from the neighboring
      // rank.
      std::vector<Point<spacedim>> &real_qpoints =
        recv_qpoints.at(neigh_rank).at({neigh_id, local_neigh});

      const auto &JxWs = recv_jxws.at(neigh_rank).at({neigh_id, local_neigh});

      std::vector<Tensor<1, spacedim>> &normals =
        recv_normals.at(neigh_rank).at({neigh_id, local_neigh});

      // Apply the necessary scalings due to the fact that functions are defined
      // on a bbox.
      std::vector<Point<spacedim>> final_unit_q_points;
      std::transform(real_qpoints.begin(),
                     real_qpoints.end(),
                     std::back_inserter(final_unit_q_points),
                     [&](const Point<spacedim> &p) {
                       return bbox.real_to_unit(p);
                     });

      // Since we received normal vectors from a neighbor, we have to swap
      // the
      // // sign of the vector in order to have outward normals.
      for (unsigned int q = 0; q < final_unit_q_points.size(); ++q)
        normals[q] *= -1;


      NonMatching::ImmersedSurfaceQuadrature<dim, spacedim> surface_quad(
        final_unit_q_points, JxWs, normals);

      agglomerated_isv =
        std::make_unique<NonMatching::FEImmersedSurfaceValues<spacedim>>(
          *box_mapping, *fe, surface_quad, agglomeration_face_flags);


      agglomerated_isv->reinit(cell);

      std::pair<const FEValuesBase<dim, spacedim> &,
                const FEValuesBase<dim, spacedim> &>
        my_p(*agglomerated_isv, *agglomerated_isv);

      return my_p;
    }
}



template <int dim, int spacedim>
template <typename SparsityPatternType, typename Number>
void
AgglomerationHandler<dim, spacedim>::create_agglomeration_sparsity_pattern(
  SparsityPatternType             &dsp,
  const AffineConstraints<Number> &constraints,
  const bool                       keep_constrained_dofs,
  const types::subdomain_id        subdomain_id)
{
  Assert(n_agglomerations > 0,
         ExcMessage("The agglomeration has not been set up correctly."));
  Assert(dsp.empty(),
         ExcMessage(
           "The Sparsity pattern must be empty upon calling this function."));

  const IndexSet &locally_owned_dofs = agglo_dh.locally_owned_dofs();
  const IndexSet  locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(agglo_dh);

  if constexpr (std::is_same_v<SparsityPatternType, DynamicSparsityPattern>)
    dsp.reinit(locally_owned_dofs.size(),
               locally_owned_dofs.size(),
               locally_relevant_dofs);
  else if constexpr (std::is_same_v<SparsityPatternType,
                                    TrilinosWrappers::SparsityPattern>)
    dsp.reinit(locally_owned_dofs, communicator);
  else
    AssertThrow(false, ExcNotImplemented());

  // Create the sparsity pattern corresponding only to volumetric terms. The
  // fluxes needed by DG methods will be filled later.
  DoFTools::make_sparsity_pattern(
    agglo_dh, dsp, constraints, keep_constrained_dofs, subdomain_id);


  const unsigned int dofs_per_cell = agglo_dh.get_fe(0).n_dofs_per_cell();
  std::vector<types::global_dof_index> current_dof_indices(dofs_per_cell);
  std::vector<types::global_dof_index> neighbor_dof_indices(dofs_per_cell);

  // Loop over all locally owned polytopes, find the neighbor (also ghosted)
  // and add fluxes to the sparsity pattern.
  for (const auto &polytope : polytope_iterators())
    {
      if (polytope->is_locally_owned())
        {
          const unsigned int n_current_faces = polytope->n_faces();
          polytope->get_dof_indices(current_dof_indices);
          for (unsigned int f = 0; f < n_current_faces; ++f)
            {
              const auto &neigh_polytope = polytope->neighbor(f);
              if (neigh_polytope.state() == IteratorState::valid)
                {
                  neigh_polytope->get_dof_indices(neighbor_dof_indices);
                  constraints.add_entries_local_to_global(current_dof_indices,
                                                          neighbor_dof_indices,
                                                          dsp,
                                                          keep_constrained_dofs,
                                                          {});
                }
            }
        }
    }



  if constexpr (std::is_same_v<SparsityPatternType,
                               TrilinosWrappers::SparsityPattern>)
    dsp.compress();
}



template <int dim, int spacedim>
void
AgglomerationHandler<dim, spacedim>::setup_ghost_polytopes()
{
  [[maybe_unused]] const auto parallel_triangulation =
    dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(&*tria);
  Assert(parallel_triangulation != nullptr, ExcInternalError());

  const unsigned int                   n_dofs_per_cell = fe->dofs_per_cell;
  std::vector<types::global_dof_index> global_dof_indices(n_dofs_per_cell);
  for (const auto &polytope : polytope_iterators())
    if (polytope->is_locally_owned())
      {
        const CellId &master_cell_id = polytope->id();

        const auto polytope_dh = polytope->as_dof_handler_iterator(agglo_dh);
        polytope_dh->get_dof_indices(global_dof_indices);


        const auto &agglomerate = polytope->get_agglomerate();

        for (const auto &cell : agglomerate)
          {
            // interior, locally owned, cell
            for (const auto &f : cell->face_indices())
              {
                if (!cell->at_boundary(f))
                  {
                    const auto &neighbor = cell->neighbor(f);
                    if (neighbor->is_ghost())
                      {
                        // key of the map: the rank to which send the data
                        const types::subdomain_id neigh_rank =
                          neighbor->subdomain_id();

                        // inform the "standard" neighbor about the neighboring
                        // id and its master cell
                        local_cell_ids_neigh_cell[neigh_rank].emplace(
                          cell->id(), master_cell_id);

                        // inform the neighboring rank that this master cell
                        // (hence polytope) has the following DoF indices
                        local_ghost_dofs[neigh_rank].emplace(
                          master_cell_id, global_dof_indices);

                        // ...same for bounding boxes
                        const auto &bbox = bboxes[polytope->index()];
                        local_ghosted_bbox[neigh_rank].emplace(master_cell_id,
                                                               bbox);
                      }
                  }
              }
          }
      }

  recv_cell_ids_neigh_cell =
    Utilities::MPI::some_to_some(communicator, local_cell_ids_neigh_cell);

  // Exchange with neighboring ranks the neighboring bounding boxes
  recv_ghosted_bbox =
    Utilities::MPI::some_to_some(communicator, local_ghosted_bbox);

  // Exchange with neighboring ranks the neighboring ghosted DoFs
  recv_ghost_dofs =
    Utilities::MPI::some_to_some(communicator, local_ghost_dofs);
}



namespace internal
{
  template <int dim, int spacedim>
  class AgglomerationHandlerImplementation
  {
  public:
    static const FEValuesBase<dim, spacedim> &
    reinit_master(
      const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
      const unsigned int face_index,
      std::unique_ptr<NonMatching::FEImmersedSurfaceValues<spacedim>>
                                                &agglo_isv_ptr,
      const AgglomerationHandler<dim, spacedim> &handler)
    {
      Assert(handler.is_master_cell(cell),
             ExcMessage("This cell must be a master one."));

      AgglomerationIterator<dim, spacedim> it{cell, &handler};
      const auto &neigh_polytope = it->neighbor(face_index);

      const CellId polytope_in_id = cell->id();

      // Retrieve the bounding box of the agglomeration
      const auto &bbox =
        handler.bboxes[handler.master2polygon.at(cell->active_cell_index())];

      CellId polytope_out_id;
      if (neigh_polytope.state() == IteratorState::valid)
        polytope_out_id = neigh_polytope->id();
      else
        polytope_out_id = polytope_in_id; // on the boundary. Same id

      const auto &common_face =
        handler.polytope_cache.interface.at({polytope_in_id, polytope_out_id});

      std::vector<Point<spacedim>> final_unit_q_points;
      std::vector<double>          final_weights;
      std::vector<Tensor<1, dim>>  final_normals;

      const unsigned int expected_qpoints =
        common_face.size() * handler.agglomeration_face_quad.size();
      final_unit_q_points.reserve(expected_qpoints);
      final_weights.reserve(expected_qpoints);
      final_normals.reserve(expected_qpoints);


      for (const auto &[deal_cell, local_face_idx] : common_face)
        {
          handler.no_face_values->reinit(deal_cell, local_face_idx);

          const auto &q_points =
            handler.no_face_values->get_quadrature_points();
          const auto &JxWs    = handler.no_face_values->get_JxW_values();
          const auto &normals = handler.no_face_values->get_normal_vectors();

          const unsigned int n_qpoints_agglo = q_points.size();

          for (unsigned int q = 0; q < n_qpoints_agglo; ++q)
            {
              final_unit_q_points.push_back(bbox.real_to_unit(q_points[q]));
              final_weights.push_back(JxWs[q]);
              final_normals.push_back(normals[q]);
            }
        }

      NonMatching::ImmersedSurfaceQuadrature<dim, spacedim> surface_quad(
        final_unit_q_points, final_weights, final_normals);

      agglo_isv_ptr =
        std::make_unique<NonMatching::FEImmersedSurfaceValues<spacedim>>(
          *(handler.box_mapping),
          *(handler.fe),
          surface_quad,
          handler.agglomeration_face_flags);

      agglo_isv_ptr->reinit(cell);

      return *agglo_isv_ptr;
    }



    /**
     * Given an agglomeration described by the master cell `master_cell`,
     * this function:
     * - enumerates the faces of the agglomeration
     * - stores who is the neighbor, the local face indices from outside and
     * inside*/
    static void
    setup_master_neighbor_connectivity(
      const typename Triangulation<dim, spacedim>::active_cell_iterator
                                                &master_cell,
      const AgglomerationHandler<dim, spacedim> &handler)
    {
      Assert(
        handler.master_slave_relationships[master_cell
                                             ->global_active_cell_index()] ==
          -1,
        ExcMessage("The present cell with index " +
                   std::to_string(master_cell->global_active_cell_index()) +
                   "is not a master one."));

      const auto &agglomeration = handler.get_agglomerate(master_cell);
      const types::global_cell_index current_polytope_index =
        handler.master2polygon.at(master_cell->active_cell_index());

      CellId current_polytope_id = master_cell->id();


      std::set<types::global_cell_index> visited_polygonal_neighbors;

      std::map<unsigned int, CellId> face_to_neigh_id;

      std::map<unsigned int, bool> is_face_at_boundary;

      // same as above, but with CellId
      std::set<CellId> visited_polygonal_neighbors_id;
      unsigned int     ghost_counter = 0;

      for (const auto &cell : agglomeration)
        {
          const types::global_cell_index cell_index = cell->active_cell_index();

          const CellId cell_id = cell->id();

          for (const auto f : cell->face_indices())
            {
              const auto &neighboring_cell = cell->neighbor(f);

              const bool valid_neighbor =
                neighboring_cell.state() == IteratorState::valid;

              if (valid_neighbor)
                {
                  if (neighboring_cell->is_locally_owned() &&
                      !handler.are_cells_agglomerated(cell, neighboring_cell))
                    {
                      // - cell is not on the boundary,
                      // - it's not agglomerated with the neighbor. If so,
                      // it's a neighbor of the present agglomeration

                      // a new face of the agglomeration has been
                      // discovered.
                      handler.polygon_boundary[master_cell].push_back(
                        cell->face(f));

                      // global index of neighboring deal.II cell
                      const types::global_cell_index neighboring_cell_index =
                        neighboring_cell->active_cell_index();

                      // master cell for the neighboring polytope
                      const auto &master_of_neighbor =
                        handler.master_slave_relationships_iterators.at(
                          neighboring_cell_index);

                      const auto nof = cell->neighbor_of_neighbor(f);

                      if (handler.is_slave_cell(neighboring_cell))
                        {
                          // index of the neighboring polytope
                          const types::global_cell_index
                            neighbor_polytope_index = handler.master2polygon.at(
                              master_of_neighbor->active_cell_index());

                          CellId neighbor_polytope_id =
                            master_of_neighbor->id();

                          if (visited_polygonal_neighbors.find(
                                neighbor_polytope_index) ==
                              std::end(visited_polygonal_neighbors))
                            {
                              // found a neighbor

                              const unsigned int n_face =
                                handler.number_of_agglomerated_faces
                                  [current_polytope_index];

                              handler.polytope_cache.cell_face_at_boundary[{
                                current_polytope_index, n_face}] = {
                                false, master_of_neighbor};

                              is_face_at_boundary[n_face] = true;

                              ++handler.number_of_agglomerated_faces
                                  [current_polytope_index];

                              visited_polygonal_neighbors.insert(
                                neighbor_polytope_index);
                            }


                          if (handler.polytope_cache.visited_cell_and_faces
                                .find({cell_index, f}) ==
                              std::end(
                                handler.polytope_cache.visited_cell_and_faces))
                            {
                                handler.polytope_cache
                                  .interface[{
                                  current_polytope_id, neighbor_polytope_id}]
                                  .emplace_back(cell, f);

                                handler.polytope_cache.visited_cell_and_faces
                                  .insert({cell_index, f});
                            }


                          if (handler.polytope_cache.visited_cell_and_faces
                                .find({neighboring_cell_index, nof}) ==
                              std::end(
                                handler.polytope_cache.visited_cell_and_faces))
                            {
                                handler.polytope_cache
                                  .interface[{
                                  neighbor_polytope_id, current_polytope_id}]
                                  .emplace_back(neighboring_cell, nof);

                                handler.polytope_cache.visited_cell_and_faces
                                  .insert({neighboring_cell_index, nof});
                            }
                        }
                      else
                        {
                          // neighboring cell is a master

                          // save the pair of neighboring cells
                          const types::global_cell_index
                            neighbor_polytope_index =
                              handler.master2polygon.at(neighboring_cell_index);

                          CellId neighbor_polytope_id = neighboring_cell->id();

                          if (visited_polygonal_neighbors.find(
                                neighbor_polytope_index) ==
                              std::end(visited_polygonal_neighbors))
                            {
                                // found a neighbor
                                const unsigned int n_face =
                                  handler.number_of_agglomerated_faces
                                    [current_polytope_index];


                                handler.polytope_cache.cell_face_at_boundary[{
                                  current_polytope_index, n_face}] = {
                                  false, neighboring_cell};

                                is_face_at_boundary[n_face] = true;

                                ++handler.number_of_agglomerated_faces
                                    [current_polytope_index];

                                visited_polygonal_neighbors.insert(
                                  neighbor_polytope_index);
                            }



                          if (handler.polytope_cache.visited_cell_and_faces
                                .find({cell_index, f}) ==
                              std::end(
                                handler.polytope_cache.visited_cell_and_faces))
                            {
                                handler.polytope_cache
                                  .interface[{
                                  current_polytope_id, neighbor_polytope_id}]
                                  .emplace_back(cell, f);

                                handler.polytope_cache.visited_cell_and_faces
                                  .insert({cell_index, f});
                            }

                          if (handler.polytope_cache.visited_cell_and_faces
                                .find({neighboring_cell_index, nof}) ==
                              std::end(
                                handler.polytope_cache.visited_cell_and_faces))
                            {
                                handler.polytope_cache
                                  .interface[{
                                  neighbor_polytope_id, current_polytope_id}]
                                  .emplace_back(neighboring_cell, nof);

                                handler.polytope_cache.visited_cell_and_faces
                                  .insert({neighboring_cell_index, nof});
                            }
                        }
                    }
                  else if (neighboring_cell->is_ghost())
                    {
                      const auto nof = cell->neighbor_of_neighbor(f);

                      // from neighboring rank,receive the association
                      // between standard cell ids and neighboring polytope.
                      // This tells to the current rank that the
                      // neighboring cell has the following CellId as master
                      // cell.
                      const auto &check_neigh_poly_ids =
                        handler.recv_cell_ids_neigh_cell.at(
                          neighboring_cell->subdomain_id());

                      const CellId neighboring_cell_id = neighboring_cell->id();

                      const CellId &check_neigh_polytope_id =
                        check_neigh_poly_ids.at(neighboring_cell_id);

                      // const auto master_index =
                      // master_indices[ghost_counter];

                      if (visited_polygonal_neighbors_id.find(
                            check_neigh_polytope_id) ==
                          std::end(visited_polygonal_neighbors_id))
                        {
                          handler.polytope_cache.cell_face_at_boundary[{
                            current_polytope_index,
                            handler.number_of_agglomerated_faces
                              [current_polytope_index]}] = {false,
                                                            neighboring_cell};


                          // record the cell id of the neighboring polytope
                          handler.polytope_cache.ghosted_master_id[{
                            current_polytope_id,
                            handler.number_of_agglomerated_faces
                              [current_polytope_index]}] =
                            check_neigh_polytope_id;


                          const unsigned int n_face =
                            handler.number_of_agglomerated_faces
                              [current_polytope_index];

                          face_to_neigh_id[n_face] = check_neigh_polytope_id;

                          is_face_at_boundary[n_face] = false;


                          // increment number of faces
                          ++handler.number_of_agglomerated_faces
                              [current_polytope_index];

                          visited_polygonal_neighbors_id.insert(
                            check_neigh_polytope_id);

                          // ghosted polytope has been found, increment
                          // ghost counter
                          ++ghost_counter;
                        }



                      if (handler.polytope_cache.visited_cell_and_faces_id.find(
                            {cell_id, f}) ==
                          std::end(
                            handler.polytope_cache.visited_cell_and_faces_id))
                        {
                            handler.polytope_cache
                              .interface[{
                              current_polytope_id, check_neigh_polytope_id}]
                              .emplace_back(cell, f);

                            handler.polytope_cache.visited_cell_and_faces_id
                              .insert({cell_id, f});
                        }


                      if (handler.polytope_cache.visited_cell_and_faces_id.find(
                            {neighboring_cell_id, nof}) ==
                          std::end(
                            handler.polytope_cache.visited_cell_and_faces_id))
                        {
                            handler.polytope_cache
                              .interface[{
                              check_neigh_polytope_id, current_polytope_id}]
                              .emplace_back(neighboring_cell, nof);

                            handler.polytope_cache.visited_cell_and_faces_id
                              .insert({neighboring_cell_id, nof});
                        }
                    }
                }
              else if (cell->face(f)->at_boundary())
                {
                  // Boundary face of a boundary cell.
                  // Note that the neighboring cell must be invalid.

                  handler.polygon_boundary[master_cell].push_back(
                    cell->face(f));

                  if (visited_polygonal_neighbors.find(
                        std::numeric_limits<unsigned int>::max()) ==
                      std::end(visited_polygonal_neighbors))
                    {
                      // boundary face. Notice that `neighboring_cell` is
                      // invalid here.
                      handler.polytope_cache.cell_face_at_boundary[{
                        current_polytope_index,
                        handler.number_of_agglomerated_faces
                          [current_polytope_index]}] = {true, neighboring_cell};

                      const unsigned int n_face =
                        handler
                          .number_of_agglomerated_faces[current_polytope_index];

                      is_face_at_boundary[n_face] = true;

                      ++handler
                          .number_of_agglomerated_faces[current_polytope_index];

                      visited_polygonal_neighbors.insert(
                        std::numeric_limits<unsigned int>::max());
                    }



                  if (handler.polytope_cache.visited_cell_and_faces.find(
                        {cell_index, f}) ==
                      std::end(handler.polytope_cache.visited_cell_and_faces))
                    {
                        handler.polytope_cache
                          .interface[{
                          current_polytope_id, current_polytope_id}]
                          .emplace_back(cell, f);

                        handler.polytope_cache.visited_cell_and_faces.insert(
                          {cell_index, f});
                    }
                }
            } // loop over faces
        }     // loop over all cells of agglomerate



      if (ghost_counter > 0)
        {
          const auto parallel_triangulation = dynamic_cast<
            const dealii::parallel::TriangulationBase<dim, spacedim> *>(
            &(*handler.tria));

          const unsigned int n_faces_current_poly =
            handler.number_of_agglomerated_faces[current_polytope_index];

          // Communicate to neighboring ranks that current_polytope_id has
          // a number of faces equal to n_faces_current_poly faces:
          // current_polytope_id -> n_faces_current_poly
          for (const unsigned int neigh_rank :
               parallel_triangulation->ghost_owners())
            {
              handler.local_n_faces[neigh_rank].emplace(current_polytope_id,
                                                        n_faces_current_poly);

              handler.local_bdary_info[neigh_rank].emplace(current_polytope_id,
                                                           is_face_at_boundary);

              handler.local_ghosted_master_id[neigh_rank].emplace(
                current_polytope_id, face_to_neigh_id);
            }
        }
    }
  };



} // namespace internal


// explicit instantiations
#include "dofs/agglomeration_handler.inst"

DEAL_II_NAMESPACE_CLOSE
