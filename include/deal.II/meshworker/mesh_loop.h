// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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


#ifndef dealii__mesh_worker_mesh_loop_h
#define dealii__mesh_worker_mesh_loop_h

#include <deal.II/base/config.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/meshworker/loop.h>
#include <deal.II/meshworker/local_integrator.h>
#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/assemble_flags.h>

#include <functional>

DEAL_II_NAMESPACE_OPEN

template <typename> class TriaActiveIterator;

namespace MeshWorker
{
  /**
   * A low level work function of this namespace.
   *
   * @ingroup MeshWorker
   * @author Luca Heltai, 2017
   */
  template <class CellIteratorType,
            class ScratchData, class CopyData>
  void mesh_loop(const CellIteratorType &begin,
                 const typename identity<CellIteratorType>::type &end,

                 const typename identity<std::function<void (const CellIteratorType &, ScratchData &, CopyData &)>>::type &cell_worker,
                 const typename identity<std::function<void (const CopyData &)>>::type &copyer,

                 ScratchData &scratch_data,
                 CopyData &copy_data,

                 const AssembleFlags flags = assemble_own_cells,

                 const typename identity<std::function<void (const CellIteratorType &, const unsigned int&, ScratchData &, CopyData &)>>::type &boundary_worker=
                  std::function<void (const CellIteratorType &, const unsigned int&, ScratchData &, CopyData &)>(),

                 const typename identity<std::function<void (const CellIteratorType &, const unsigned int &, const unsigned int &,
                                           const CellIteratorType &, const unsigned int &, const unsigned int &,
                                           ScratchData &, CopyData &)>>::type &face_worker=
                  std::function<void (const CellIteratorType &, const unsigned int &, const unsigned int &,
                                      const CellIteratorType &, const unsigned int &, const unsigned int &,
                                      ScratchData &, CopyData &)>(),

                 const unsigned int   queue_length = 2*MultithreadInfo::n_threads(),
                 const unsigned int   chunk_size = 8)
  {
    auto cell_action = [&] (const CellIteratorType &cell, ScratchData &scratch, CopyData &copy)
    {
      const bool ignore_subdomain = (cell->get_triangulation().locally_owned_subdomain()
                                     == numbers::invalid_subdomain_id);

      types::subdomain_id current_subdomain_id = /*(cell->is_level_cell())
                                 ? cell->level_subdomain_id()
                                 : */cell->subdomain_id();

      const bool own_cell = ignore_subdomain || (current_subdomain_id == cell->get_triangulation().locally_owned_subdomain());

      if ((!ignore_subdomain) && (current_subdomain_id == numbers::artificial_subdomain_id))
        return;

      if ( (flags & (assemble_cells_first)) &&
           ( ((flags & (assemble_own_cells)) && own_cell)
             || ( (flags & assemble_ghost_cells) && !own_cell) ) )
        cell_worker(cell, scratch, copy);

      if (flags & assemble_own_faces)
        for (unsigned int face_no=0; face_no < GeometryInfo<CellIteratorType::AccessorType::Container::dimension>::faces_per_cell; ++face_no)
          {
            typename CellIteratorType::AccessorType::Container::face_iterator face = cell->face(face_no);
            if (cell->at_boundary(face_no) && !cell->has_periodic_neighbor(face_no))
              {
                // only integrate boundary faces of own cells
                if ( (flags & assemble_boundary_faces) && own_cell)
                  {
                    boundary_worker(cell, face_no, scratch, copy);
                  }
              }
            else if (flags & assemble_own_interior_faces)
              {
                // Interior face
                TriaIterator<typename CellIteratorType::AccessorType> neighbor = cell->neighbor_or_periodic_neighbor(face_no);

                types::subdomain_id neighbor_subdomain_id = numbers::artificial_subdomain_id;
//                if (neighbor->is_level_cell())
//                  neighbid = neighbor->level_subdomain_id();
                //subdomain id is only valid for active cells
                /*else */if (neighbor->active())
                  neighbor_subdomain_id = neighbor->subdomain_id();

                const bool own_neighbor = ignore_subdomain ||
                                          (neighbor_subdomain_id == cell->get_triangulation().locally_owned_subdomain());

                // skip all faces between two ghost cells
                if (!own_cell && !own_neighbor)
                  continue;

                // skip if the user doesn't want faces between own cells
                if (own_cell && own_neighbor && !(flags & assemble_own_interior_faces))
                  continue;

                // skip face to ghost
                if (own_cell != own_neighbor && !(flags & assemble_ghost_faces))
                  continue;

                // Deal with refinement edges from the refined side. Assuming one-irregular
                // meshes, this situation should only occur if both cells are active.
                const bool periodic_neighbor = cell->has_periodic_neighbor(face_no);

                if ((!periodic_neighbor && cell->neighbor_is_coarser(face_no))
                    || (periodic_neighbor && cell->periodic_neighbor_is_coarser(face_no)))
                  {
                    Assert(!cell->has_children(), ExcInternalError());
                    Assert(!neighbor->has_children(), ExcInternalError());

                    // skip if only one processor needs to assemble the face
                    // to a ghost cell and the fine cell is not ours.
                    if (!own_cell && (flags & assemble_ghost_faces_once))
                      continue;

                    const std::pair<unsigned int, unsigned int> neighbor_face_no
                      = periodic_neighbor?
                        cell->periodic_neighbor_of_coarser_periodic_neighbor(face_no):
                        cell->neighbor_of_coarser_neighbor(face_no);
                    const typename CellIteratorType::AccessorType::Container::face_iterator nface
                      = neighbor->face(neighbor_face_no.first);

                    face_worker(cell, face_no, numbers::invalid_unsigned_int,
                                neighbor, neighbor_face_no.first, neighbor_face_no.second,
                                scratch, copy);

                    if (flags & assemble_own_interior_faces_both)
                      {
                        face_worker(neighbor, neighbor_face_no.first, neighbor_face_no.second,
                                    cell, face_no, numbers::invalid_unsigned_int,
                                    scratch, copy);
                      }
                  }
                else
                  {
                    // If iterator is active and neighbor is refined, skip
                    // internal face.
                    if (internal::is_active_iterator(cell) && neighbor->has_children())
                      continue;

                    // Now neighbor is on same level, double-check this:
                    Assert(cell->level()==neighbor->level(), ExcInternalError());

                    // If we own both cells only do faces from one side (unless
                    // AssembleFlags says otherwise). Here, we rely on cell comparison
                    // that will look at cell->index().
                    if (own_cell && own_neighbor
                        && (flags & assemble_own_interior_faces_once)
                        && (neighbor < cell))
                      continue;

                    // independent of loop_control.faces_to_ghost,
                    // we only look at faces to ghost on the same level once
                    // (only where own_cell=true and own_neighbor=false)
                    if (!own_cell)
                      continue;

                    // now only one processor assembles faces_to_ghost. We let the
                    // processor with the smaller (level-)subdomain id assemble the
                    // face.
                    if (own_cell && !own_neighbor
                        && (flags & assemble_ghost_faces_once)
                        && (neighbor_subdomain_id < current_subdomain_id))
                      continue;

                    const unsigned int neighbor_face_no = periodic_neighbor?
                                                          cell->periodic_neighbor_face_no(face_no):
                                                          cell->neighbor_face_no(face_no);
                    Assert (periodic_neighbor || neighbor->face(neighbor_face_no) == face, ExcInternalError());

                    face_worker(cell, face_no, numbers::invalid_unsigned_int,
                                neighbor, neighbor_face_no, numbers::invalid_unsigned_int,
                                scratch, copy);

                  }
              }
          } // faces

      // Execute this, if faces
      // have to be handled first
      if ((flags & assemble_own_cells) && !(flags & assemble_cells_first) &&
          ( ((flags & assemble_own_cells) && own_cell) || ((flags & assemble_ghost_cells) && !own_cell)))
        cell_worker(cell, scratch, copy);
    };


    // Loop over all cells
    WorkStream::run(begin, end,
                    cell_action, copyer,
                    scratch_data, copy_data,
                    queue_length, chunk_size);
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
