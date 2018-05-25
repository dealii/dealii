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


#ifndef dealii_mesh_worker_mesh_loop_h
#define dealii_mesh_worker_mesh_loop_h

#include <deal.II/base/config.h>

#include <deal.II/base/template_constraints.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/meshworker/assemble_flags.h>
#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/local_integrator.h>
#include <deal.II/meshworker/loop.h>

#include <functional>

DEAL_II_NAMESPACE_OPEN

template <typename>
class TriaActiveIterator;

namespace MeshWorker
{
  /**
   * This function extends the WorkStream concept to work on meshes
   * (cells and/or faces) and handles the complicated logic for
   * work on adaptively refined faces
   * and parallel computation (work on faces to ghost neighbors for example).
   * The @p mesh_loop can be used to simplify operations on cells (for example
   * assembly), on boundaries (Neumann type boundary conditions), or on
   * interior faces (for example in discontinuous Galerkin methods).
   *
   * For uniformly refined meshes, it would be relatively easy to use
   * WorkStream::run() with a @p cell_worker that also loops over faces, and
   * takes care of assembling face terms depending on the current and neighbor
   * cell. All user codes that do these loops would then need to insert
   * manually the logic that identifies, for every face of the current cell,
   * the neighboring cell, and the face index on the neighboring cell that
   * corresponds to the current face.
   *
   * This is more complicated if local refinement is enabled and the current or
   * neighbor cells have hanging nodes. In this case it is also necessary to
   * identify the corresponding subface on either the current or the neighbor
   * faces.
   *
   * This method externalises that logic (which is independent from user codes)
   * and separates the assembly of face terms (internal faces, boundary faces,
   * or faces between different subdomain ids on parallel computations) from
   * the assembling on cells, allowing the user to specify two additional
   * workers (a @p cell_worker, a @p boundary_worker, and a @p face_worker) that
   * are called automatically in each @p cell, according to the specific
   * AssembleFlags @p flags that are passed. The @p cell_worker is passed the
   * cell identifier, a ScratchData object, and a CopyData object, following
   * the same principles of WorkStream::run. Internally the function passes to
   * @p boundary_worker, in addition to the above, also a @p face_no parameter
   * that identifies the face on which the integration should be performed. The
   * @p face_worker instead needs to identify the current face unambiguously both on
   * the cell and on the neighboring cell, and it is therefore called with six
   * arguments (three for each cell: the actual cell, the face index, and
   * the subface_index. If no subface integration is needed, then the
   * subface_index is numbers::invalid_unsigned_int) in addition to the usual
   * ScratchData and CopyData objects.
   *
   * If the flag AssembleFlags::assemble_own_cells is passed, then the default
   * behavior is to first loop over faces and do the work there, and then
   * compute the actual work on the cell. It is possible to perform the
   * integration on the cells after working on faces, by adding the flag
   * AssembleFlags::cells_after_faces.
   *
   * If the flag AssembleFlags::assemble_own_interior_faces_once is specified,
   * then each interior face is visited only once, and the @p face_worker is
   * assumed to integrate all face terms at once (and add contributions to both
   * sides of the face in a discontinuous Galerkin setting).
   *
   * This method is equivalent to the WorkStream::run() method when
   * AssembleFlags contains only @p assemble_own_cells, and can be used as a
   * drop-in replacement for that method.
   *
   * The two data types ScratchData and CopyData need to have a working copy
   * constructor. ScratchData is only used in the worker function, while
   * CopyData is the object passed from the worker to the copier.
   *
   * The queue_length argument indicates the number of items that can be live at
   * any given time. Each item consists of chunk_size elements of the input
   * stream that will be worked on by the worker and copier functions one after
   * the other on the same thread.
   *
   * If your data objects are large, or their constructors are expensive, it is
   * helpful to keep in mind that queue_length copies of the ScratchData object
   * and queue_length*chunk_size copies of the CopyData object are generated.
   *
   * @note More information about requirements on template types and meaning
   * of @p queue_length and @p chunk_size can be found in the documentation of the
   * WorkStream namespace and its members.
   *
   * @ingroup MeshWorker
   * @author Luca Heltai and Timo Heister, 2017
   */
  template <class CellIteratorType, class ScratchData, class CopyData>
  void
  mesh_loop(const CellIteratorType &                         begin,
            const typename identity<CellIteratorType>::type &end,

            const typename identity<std::function<
              void(const CellIteratorType &, ScratchData &, CopyData &)>>::type
              &cell_worker,
            const typename identity<std::function<void(const CopyData &)>>::type
              &copier,

            const ScratchData &sample_scratch_data,
            const CopyData &   sample_copy_data,

            const AssembleFlags flags = assemble_own_cells,

            const typename identity<std::function<void(const CellIteratorType &,
                                                       const unsigned int &,
                                                       ScratchData &,
                                                       CopyData &)>>::type
              &boundary_worker = std::function<void(const CellIteratorType &,
                                                    const unsigned int &,
                                                    ScratchData &,
                                                    CopyData &)>(),

            const typename identity<std::function<void(const CellIteratorType &,
                                                       const unsigned int &,
                                                       const unsigned int &,
                                                       const CellIteratorType &,
                                                       const unsigned int &,
                                                       const unsigned int &,
                                                       ScratchData &,
                                                       CopyData &)>>::type
              &face_worker = std::function<void(const CellIteratorType &,
                                                const unsigned int &,
                                                const unsigned int &,
                                                const CellIteratorType &,
                                                const unsigned int &,
                                                const unsigned int &,
                                                ScratchData &,
                                                CopyData &)>(),

            const unsigned int queue_length = 2 * MultithreadInfo::n_threads(),
            const unsigned int chunk_size   = 8)
  {
    Assert(
      (!cell_worker) == !(flags & work_on_cells),
      ExcMessage(
        "If you specify a cell_worker, you need to set assemble_own_cells or assemble_ghost_cells."));

    Assert(
      (flags &
       (assemble_own_interior_faces_once | assemble_own_interior_faces_both)) !=
        (assemble_own_interior_faces_once | assemble_own_interior_faces_both),
      ExcMessage(
        "You can only specify assemble_own_interior_faces_once OR assemble_own_interior_faces_both."));

    Assert(
      (flags & (assemble_ghost_faces_once | assemble_ghost_faces_both)) !=
        (assemble_ghost_faces_once | assemble_ghost_faces_both),
      ExcMessage(
        "You can only specify assemble_ghost_faces_once OR assemble_ghost_faces_both."));

    Assert(
      !(flags & cells_after_faces) ||
        (flags & (assemble_own_cells | assemble_ghost_cells)),
      ExcMessage(
        "The option cells_after_faces only makes sense if you assemble on cells."));

    Assert((!face_worker) == !(flags & work_on_faces),
           ExcMessage(
             "If you specify a face_worker, assemble_face_* needs to be set."));

    Assert(
      (!boundary_worker) == !(flags & assemble_boundary_faces),
      ExcMessage(
        "If you specify a boundary_worker, assemble_boundary_faces needs to be set."));

    auto cell_action = [&](const CellIteratorType &cell,
                           ScratchData &           scratch,
                           CopyData &              copy) {
      // First reset the CopyData class to the empty copy_data given by the
      // user.
      copy = sample_copy_data;

      const bool ignore_subdomain =
        (cell->get_triangulation().locally_owned_subdomain() ==
         numbers::invalid_subdomain_id);

      types::subdomain_id current_subdomain_id =
        (cell->is_level_cell() ? cell->level_subdomain_id() :
                                 cell->subdomain_id());

      const bool own_cell =
        ignore_subdomain ||
        (current_subdomain_id ==
         cell->get_triangulation().locally_owned_subdomain());

      if ((!ignore_subdomain) &&
          (current_subdomain_id == numbers::artificial_subdomain_id))
        return;

      if (!(flags & (cells_after_faces)) &&
          (((flags & (assemble_own_cells)) && own_cell) ||
           ((flags & assemble_ghost_cells) && !own_cell)))
        cell_worker(cell, scratch, copy);

      if (flags & (work_on_faces | work_on_boundary))
        for (unsigned int face_no = 0;
             face_no < GeometryInfo<CellIteratorType::AccessorType::Container::
                                      dimension>::faces_per_cell;
             ++face_no)
          {
            typename CellIteratorType::AccessorType::Container::face_iterator
              face = cell->face(face_no);
            if (cell->at_boundary(face_no) &&
                !cell->has_periodic_neighbor(face_no))
              {
                // only integrate boundary faces of own cells
                if ((flags & assemble_boundary_faces) && own_cell)
                  boundary_worker(cell, face_no, scratch, copy);
              }
            else
              {
                // interior face, potentially assemble
                TriaIterator<typename CellIteratorType::AccessorType> neighbor =
                  cell->neighbor_or_periodic_neighbor(face_no);

                types::subdomain_id neighbor_subdomain_id =
                  numbers::artificial_subdomain_id;
                if (neighbor->is_level_cell())
                  neighbor_subdomain_id = neighbor->level_subdomain_id();
                // subdomain id is only valid for active cells
                else if (neighbor->active())
                  neighbor_subdomain_id = neighbor->subdomain_id();

                const bool own_neighbor =
                  ignore_subdomain ||
                  (neighbor_subdomain_id ==
                   cell->get_triangulation().locally_owned_subdomain());

                // skip all faces between two ghost cells
                if (!own_cell && !own_neighbor)
                  continue;

                // skip if the user doesn't want faces between own cells
                if (own_cell && own_neighbor &&
                    !(flags & (assemble_own_interior_faces_both |
                               assemble_own_interior_faces_once)))
                  continue;

                // skip face to ghost
                if (own_cell != own_neighbor &&
                    !(flags &
                      (assemble_ghost_faces_both | assemble_ghost_faces_once)))
                  continue;

                // Deal with refinement edges from the refined side. Assuming
                // one-irregular meshes, this situation should only occur if
                // both cells are active.
                const bool periodic_neighbor =
                  cell->has_periodic_neighbor(face_no);

                if ((!periodic_neighbor &&
                     cell->neighbor_is_coarser(face_no)) ||
                    (periodic_neighbor &&
                     cell->periodic_neighbor_is_coarser(face_no)))
                  {
                    Assert(!cell->has_children(), ExcInternalError());
                    Assert(!neighbor->has_children(), ExcInternalError());

                    // skip if only one processor needs to assemble the face
                    // to a ghost cell and the fine cell is not ours.
                    if (!own_cell && (flags & assemble_ghost_faces_once))
                      continue;

                    const std::pair<unsigned int, unsigned int>
                      neighbor_face_no =
                        periodic_neighbor ?
                          cell->periodic_neighbor_of_coarser_periodic_neighbor(
                            face_no) :
                          cell->neighbor_of_coarser_neighbor(face_no);

                    face_worker(cell,
                                face_no,
                                numbers::invalid_unsigned_int,
                                neighbor,
                                neighbor_face_no.first,
                                neighbor_face_no.second,
                                scratch,
                                copy);

                    if (flags & assemble_own_interior_faces_both)
                      {
                        // If own faces are to be assembled from both sides,
                        // call the faceworker again with swapped arguments.
                        // This is because we won't be looking at an adaptively
                        // refined edge coming from the other side.
                        face_worker(neighbor,
                                    neighbor_face_no.first,
                                    neighbor_face_no.second,
                                    cell,
                                    face_no,
                                    numbers::invalid_unsigned_int,
                                    scratch,
                                    copy);
                      }
                  }
                else
                  {
                    // If iterator is active and neighbor is refined, skip
                    // internal face.
                    if (internal::is_active_iterator(cell) &&
                        neighbor->has_children())
                      continue;

                    // Now neighbor is on same level, double-check this:
                    Assert(cell->level() == neighbor->level(),
                           ExcInternalError());

                    // If we own both cells only do faces from one side (unless
                    // AssembleFlags says otherwise). Here, we rely on cell
                    // comparison that will look at cell->index().
                    if (own_cell && own_neighbor &&
                        (flags & assemble_own_interior_faces_once) &&
                        (neighbor < cell))
                      continue;

                    // We only look at faces to ghost on the same level once
                    // (only where own_cell=true and own_neighbor=false)
                    if (!own_cell)
                      continue;

                    // now only one processor assembles faces_to_ghost. We let
                    // the processor with the smaller (level-)subdomain id
                    // assemble the face.
                    if (own_cell && !own_neighbor &&
                        (flags & assemble_ghost_faces_once) &&
                        (neighbor_subdomain_id < current_subdomain_id))
                      continue;

                    const unsigned int neighbor_face_no =
                      periodic_neighbor ?
                        cell->periodic_neighbor_face_no(face_no) :
                        cell->neighbor_face_no(face_no);
                    Assert(periodic_neighbor ||
                             neighbor->face(neighbor_face_no) == face,
                           ExcInternalError());

                    face_worker(cell,
                                face_no,
                                numbers::invalid_unsigned_int,
                                neighbor,
                                neighbor_face_no,
                                numbers::invalid_unsigned_int,
                                scratch,
                                copy);
                  }
              }
          } // faces

      // Execute the cell_worker if faces are handled before cells
      if ((flags & cells_after_faces) &&
          (((flags & assemble_own_cells) && own_cell) ||
           ((flags & assemble_ghost_cells) && !own_cell)))
        cell_worker(cell, scratch, copy);
    };

    // Submit to workstream
    WorkStream::run(begin,
                    end,
                    cell_action,
                    copier,
                    sample_scratch_data,
                    sample_copy_data,
                    queue_length,
                    chunk_size);
  }
} // namespace MeshWorker

DEAL_II_NAMESPACE_CLOSE

#endif
