// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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


#ifndef dealii_mesh_worker_mesh_loop_h
#define dealii_mesh_worker_mesh_loop_h

#include <deal.II/base/config.h>

#include <deal.II/base/template_constraints.h>
#include <deal.II/base/types.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/meshworker/assemble_flags.h>
#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/local_integrator.h>
#include <deal.II/meshworker/loop.h>

#include <functional>
#include <type_traits>

DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
template <typename>
class TriaActiveIterator;
#endif

namespace MeshWorker
{
  namespace internal
  {
    /**
     * A helper class to provide a type definition for the underlying cell
     * iterator type.
     */
    template <class CellIteratorType>
    struct CellIteratorBaseType
    {
      /**
       * Type definition for the cell iterator type.
       */
      using type = CellIteratorType;
    };

    /**
     * A helper class to provide a type definition for the underlying cell
     * iterator type.
     *
     * This specialization is for IteratorRange, which may have either a
     * TriaActiveIterator or a FilteredIterator as its base type.
     */
    template <class CellIteratorType>
    struct CellIteratorBaseType<IteratorOverIterators<CellIteratorType>>
    {
      /**
       * Type definition for the cell iterator type.
       */
      // Since we can have filtered iterators and the like as template
      // arguments, we recursivelyremove the template layers to retrieve the
      // underlying iterator type.
      using type = typename CellIteratorBaseType<CellIteratorType>::type;
    };

    /**
     * A helper class to provide a type definition for the underlying cell
     * iterator type.
     *
     * This specialization is for FilteredIterator, which may have either a
     * TriaActiveIterator as its base type, or may be nested with another
     * FilteredIterator as the type to iterate over.
     */
    template <class CellIteratorType>
    struct CellIteratorBaseType<FilteredIterator<CellIteratorType>>
    {
      /**
       * Type definition for the cell iterator type.
       */
      // Since we can have nested filtered iterators, we recursively
      // remove the template layers to retrieve the underlying iterator type.
      using type = typename CellIteratorBaseType<CellIteratorType>::type;
    };
  } // namespace internal

#ifdef DOXYGEN
  /**
   * This alias introduces a friendly and short name for the function type
   * for the cell worker used in mesh_loop().
   */
  using CellWorkerFunctionType = std::function<
    void(const CellIteratorBaseType &, ScratchData &, CopyData &)>;

  /**
   * This alias introduces a friendly and short name for the function type
   * for the cell worker used in mesh_loop().
   */
  using CopierFunctionType = std::function<void(const CopyData &)>;

  /**
   * This alias introduces a friendly and short name for the function type
   * for the boundary worker used in mesh_loop().
   */
  using BoundaryWorkerFunctionType =
    std::function<void(const CellIteratorBaseType &,
                       const unsigned int,
                       ScratchData &,
                       CopyData &)>;

  /**
   * This alias introduces a friendly and short name for the function type
   * for the face worker used in mesh_loop().
   */
  using FaceWorkerFunctionType =
    std::function<void(const CellIteratorBaseType &,
                       const unsigned int,
                       const unsigned int,
                       const CellIteratorBaseType &,
                       const unsigned int,
                       const unsigned int,
                       ScratchData &,
                       CopyData &)>;
#endif

  /**
   * This function extends the WorkStream concept to work on meshes
   * (cells and/or faces) and handles the complicated logic for
   * work on adaptively refined faces
   * and parallel computation (work on faces to ghost neighbors for example).
   * The @p mesh_loop can be used to simplify operations on cells (for example
   * assembly), on boundaries (Neumann type boundary conditions), or on
   * interior faces (for example in discontinuous Galerkin methods). The
   * function is used in a number of tutorials, including step-12, step-16,
   * and step-47, to name just a few.
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
   * This method externalizes that logic (which is independent from user codes)
   * and separates the assembly of face terms (internal faces, boundary faces,
   * or faces between different subdomain ids on parallel computations) from
   * the assembling on cells, allowing the user to specify two additional
   * workers (a @p cell_worker, a @p boundary_worker, and a @p face_worker) that
   * are called automatically in each @p cell, according to the specific
   * AssembleFlags @p flags that are passed. The @p cell_worker is passed the
   * cell identifier, a ScratchData object, and a CopyData object, following
   * the same principles of WorkStream::run(). Internally the function passes to
   * @p boundary_worker, in addition to the above, also a @p face_no parameter
   * that identifies the face on which the integration should be performed. The
   * @p face_worker instead needs to identify the current face unambiguously
   * both on the cell and on the neighboring cell, and it is therefore called
   * with six arguments (three for each cell: the actual cell, the face index,
   * and the subface_index. If no subface integration is needed, then the
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
   * CopyData is the object passed from the worker to the copier. The CopyData
   * object is reset to the value provided to this function every time this
   * function visits a new cell (where it then calls the cell and face
   * workers). In other words, no state carries over between calling the
   * `copier` on one cell and the `cell_worker`/`face_worker`/`boundary_worker`
   * functions on the next cell, and user code needs not reset the copy
   * object either at the beginning of the cell integration or end of the
   * copy operation. Resetting the state of the `copier ` inside of a
   * `face_worker` or `boundary_worker` constitutes a bug, and may lead to
   * some unexpected results. The following example shows what is not
   * permissible, as the copier is potentially shared among numerous faces
   * on a cell:
   * @code
   *
   * using ScratchData      = MeshWorker::ScratchData<dim, spacedim>;
   * using CopyData         = MeshWorker::CopyData<1, 1, 1>;
   * using CellIteratorType = decltype(dof_handler.begin_active());
   *
   * ScratchData            scratch(...);
   * CopyData               copy(...);
   *
   * std::function<void(const CellIteratorType &, ScratchData &, CopyData &)>
   *   empty_cell_worker;
   *
   * auto boundary_worker = [...] (
   *   const CellIteratorType &cell,
   *   const unsigned int      face,
   *   ScratchData            &scratch_data,
   *   CopyData               &copy_data)
   * {
   *  const auto &fe_face_values = scratch_data.reinit(cell, face);
   *  copy_data = CopyData(...); // This is an error, as we lose the
   *                             // accumulation that has been performed on
   *                             // other boundary faces of the same cell.
   *
   *  for (unsigned int q_point = 0;
   *       q_point < fe_face_values.n_quadrature_points;
   *       ++q_point)
   *    {
   *      copy_data.vectors[0][0] += 1.0 * fe_face_values.JxW(q_point);
   *    }
   * };
   *
   * double value = 0;
   * auto copier = [...](const CopyData &copy_data)
   * {
   *   value += copy_data.vectors[0][0]; // Contributions from some faces may
   *                                     // be missing.
   * };
   *
   * MeshWorker::mesh_loop(dof_handler.active_cell_iterators(),
   *                       empty_cell_worker, copier,
   *                       scratch, copy,
   *                       MeshWorker::assemble_boundary_faces,
   *                       boundary_worker);
   * @endcode
   *
   * The queue_length argument indicates the number of items that can be live at
   * any given time. Each item consists of chunk_size elements of the input
   * stream that will be worked on by the worker and copier functions one after
   * the other on the same thread.
   *
   * If your data objects are large, or their constructors are expensive, it is
   * helpful to keep in mind that queue_length copies of the ScratchData object
   * and `queue_length*chunk_size` copies of the CopyData object are generated.
   *
   * @note The types of the function arguments and the default values (empty worker functions)
   * displayed in the Doxygen documentation here are slightly simplified
   * compared to the real types.
   *
   * @note More information about requirements on template types and meaning
   * of @p queue_length and @p chunk_size can be found in the documentation of the
   * WorkStream namespace and its members.
   *
   * @ingroup MeshWorker
   */
  template <class CellIteratorType,
            class ScratchData,
            class CopyData,
            class CellIteratorBaseType =
              typename internal::CellIteratorBaseType<CellIteratorType>::type>
  void
  mesh_loop(
#ifdef DOXYGEN
    const CellIteratorType &begin,
    const CellIteratorType &end,

    const CellWorkerFunctionType &cell_worker,
    const CopierType &            copier,

    const ScratchData &sample_scratch_data,
    const CopyData &   sample_copy_data,

    const AssembleFlags flags = assemble_own_cells,

    const BoundaryWorkerFunctionType &boundary_worker =
      BoundaryWorkerFunctionType(),

    const FaceWorkerFunctionType &face_worker = FaceWorkerFunctionType(),
    const unsigned int queue_length = 2 * MultithreadInfo::n_threads(),
    const unsigned int chunk_size   = 8
#else
    const CellIteratorType &                         begin,
    const typename identity<CellIteratorType>::type &end,

    const typename identity<std::function<
      void(const CellIteratorBaseType &, ScratchData &, CopyData &)>>::type
      &cell_worker,
    const typename identity<std::function<void(const CopyData &)>>::type
      &copier,

    const ScratchData &sample_scratch_data,
    const CopyData &   sample_copy_data,

    const AssembleFlags flags = assemble_own_cells,

    const typename identity<std::function<void(const CellIteratorBaseType &,
                                               const unsigned int,
                                               ScratchData &,
                                               CopyData &)>>::type
      &boundary_worker = std::function<void(const CellIteratorBaseType &,
                                            const unsigned int,
                                            ScratchData &,
                                            CopyData &)>(),

    const typename identity<std::function<void(const CellIteratorBaseType &,
                                               const unsigned int,
                                               const unsigned int,
                                               const CellIteratorBaseType &,
                                               const unsigned int,
                                               const unsigned int,
                                               ScratchData &,
                                               CopyData &)>>::type
      &face_worker = std::function<void(const CellIteratorBaseType &,
                                        const unsigned int,
                                        const unsigned int,
                                        const CellIteratorBaseType &,
                                        const unsigned int,
                                        const unsigned int,
                                        ScratchData &,
                                        CopyData &)>(),

    const unsigned int queue_length = 2 * MultithreadInfo::n_threads(),
    const unsigned int chunk_size   = 8
#endif
  )
  {
    Assert(
      (!cell_worker) == !(flags & work_on_cells),
      ExcMessage(
        "If you provide a cell worker function, you also need to request "
        "that work should be done on cells by setting the 'work_on_cells' flag. "
        "Conversely, if you don't provide a cell worker function, you "
        "cannot set the 'work_on_cells' flag. One of these two "
        "conditions is not satisfied."));

    Assert((flags & (assemble_own_interior_faces_once |
                     assemble_own_interior_faces_both)) !=
             (assemble_own_interior_faces_once |
              assemble_own_interior_faces_both),
           ExcMessage(
             "If you provide a face worker function, you also need to request "
             "that work should be done on interior faces by setting either the "
             "'assemble_own_interior_faces_once' flag or the "
             "'assemble_own_interior_faces_both' flag. "
             "Conversely, if you don't provide a face worker function, you "
             "cannot set either of these two flags. One of these two "
             "conditions is not satisfied."));

    Assert((flags & (assemble_ghost_faces_once | assemble_ghost_faces_both)) !=
             (assemble_ghost_faces_once | assemble_ghost_faces_both),
           ExcMessage(
             "You can only 'specify assemble_ghost_faces_once' "
             "OR 'assemble_ghost_faces_both', but not both of these flags."));

    Assert(
      !(flags & cells_after_faces) ||
        (flags & (assemble_own_cells | assemble_ghost_cells)),
      ExcMessage(
        "The option 'cells_after_faces' only makes sense if you assemble on cells."));

    Assert(
      (!face_worker) == !(flags & work_on_faces),
      ExcMessage(
        "If you provide a face worker function, you also need to request "
        "that work should be done on faces by setting the 'work_on_faces' flag. "
        "Conversely, if you don't provide a face worker function, you "
        "cannot set the 'work_on_faces' flag. One of these two "
        "conditions is not satisfied."));

    Assert(
      (!boundary_worker) == !(flags & assemble_boundary_faces),
      ExcMessage(
        "If you provide a boundary face worker function, you also need to request "
        "that work should be done on boundary faces by setting the 'assemble_boundary_faces' flag. "
        "Conversely, if you don't provide a boundary face worker function, you "
        "cannot set the 'assemble_boundary_faces' flag. One of these two "
        "conditions is not satisfied."));

    auto cell_action = [&](const CellIteratorBaseType &cell,
                           ScratchData &               scratch,
                           CopyData &                  copy) {
      // First reset the CopyData class to the empty copy_data given by the
      // user.
      copy = sample_copy_data;

      // Store the dimension in which we are working for later use
      const auto dim = cell->get_triangulation().dimension;

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
        for (const unsigned int face_no : cell->face_indices())
          {
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
                TriaIterator<typename CellIteratorBaseType::AccessorType>
                  neighbor = cell->neighbor_or_periodic_neighbor(face_no);

                types::subdomain_id neighbor_subdomain_id =
                  numbers::artificial_subdomain_id;
                if (neighbor->is_level_cell())
                  neighbor_subdomain_id = neighbor->level_subdomain_id();
                // subdomain id is only valid for active cells
                else if (neighbor->is_active())
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

                if (dim > 1 && ((!periodic_neighbor &&
                                 cell->neighbor_is_coarser(face_no) &&
                                 neighbor->is_active()) ||
                                (periodic_neighbor &&
                                 cell->periodic_neighbor_is_coarser(face_no) &&
                                 neighbor->is_active())))
                  {
                    Assert(cell->is_active(), ExcInternalError());

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
                else if (dim == 1 && cell->level() > neighbor->level())
                  {
                    // In one dimension, there is no other check to do
                    const unsigned int neighbor_face_no =
                      periodic_neighbor ?
                        cell->periodic_neighbor_face_no(face_no) :
                        cell->neighbor_face_no(face_no);
                    Assert(periodic_neighbor ||
                             neighbor->face(neighbor_face_no) ==
                               cell->face(face_no),
                           ExcInternalError());

                    face_worker(cell,
                                face_no,
                                numbers::invalid_unsigned_int,
                                neighbor,
                                neighbor_face_no,
                                numbers::invalid_unsigned_int,
                                scratch,
                                copy);

                    if (flags & assemble_own_interior_faces_both)
                      {
                        // If own faces are to be assembled from both sides,
                        // call the faceworker again with swapped arguments.
                        face_worker(neighbor,
                                    neighbor_face_no,
                                    numbers::invalid_unsigned_int,
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
                    if (dealii::internal::is_active_iterator(cell) &&
                        neighbor->has_children())
                      continue;

                    // Now neighbor is on the same refinement level.
                    // Double check.
                    Assert(!cell->neighbor_is_coarser(face_no),
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
                             neighbor->face(neighbor_face_no) ==
                               cell->face(face_no),
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

  /**
   * Same as the function above, but for iterator ranges (and, therefore,
   * filtered iterators).
   *
   * An example usage of the function for the serial case is given by
   * @code
   *
   * using ScratchData      = MeshWorker::ScratchData<dim, spacedim>;
   * using CopyData         = MeshWorker::CopyData<1, 1, 1>;
   * using CellIteratorType = decltype(dof_handler.begin_active());
   *
   * ScratchData            scratch(...);
   * CopyData               copy(...);
   *
   * auto cell_worker = [...] (
   *   const CellIteratorType &cell,
   *   ScratchData            &scratch_data,
   *   CopyData               &copy_data)
   * {
   *   ...
   * };
   *
   * auto copier = [...](const CopyData &copy_data)
   * {
   *   ...
   * };
   *
   * MeshWorker::mesh_loop(dof_handler.active_cell_iterators(),
   *                       cell_worker, copier,
   *                       scratch, copy,
   *                       MeshWorker::assemble_own_cells);
   * @endcode
   *
   * and an example usage of the function for the parallel distributed case,
   * where the copier is only to be called on locally owned cells, is given by
   * @code
   *
   * using ScratchData      = MeshWorker::ScratchData<dim, spacedim>;
   * using CopyData         = MeshWorker::CopyData<1, 1, 1>;
   * using CellIteratorType = decltype(dof_handler.begin_active());
   *
   * ScratchData            scratch(...);
   * CopyData               copy(...);
   *
   * auto cell_worker = [...] (
   *   const CellIteratorType &cell,
   *   ScratchData            &scratch_data,
   *   CopyData               &copy_data)
   * {
   *   ...
   * };
   *
   * auto copier = [...](const CopyData &copy_data)
   * {
   *   ...
   * };
   *
   * const auto filtered_iterator_range =
   *   filter_iterators(dof_handler.active_cell_iterators(),
   *                    IteratorFilters::LocallyOwnedCell());
   *
   * MeshWorker::mesh_loop(filtered_iterator_range,
   *                       cell_worker, copier,
   *                       scratch, copy,
   *                       MeshWorker::assemble_own_cells);
   * @endcode
   *
   * @ingroup MeshWorker
   */
  template <class CellIteratorType,
            class ScratchData,
            class CopyData,
            class CellIteratorBaseType =
              typename internal::CellIteratorBaseType<CellIteratorType>::type>
  void
  mesh_loop(
    IteratorRange<CellIteratorType> iterator_range,
    const typename identity<std::function<
      void(const CellIteratorBaseType &, ScratchData &, CopyData &)>>::type
      &cell_worker,
    const typename identity<std::function<void(const CopyData &)>>::type
      &copier,

    const ScratchData &sample_scratch_data,
    const CopyData &   sample_copy_data,

    const AssembleFlags flags = assemble_own_cells,

    const typename identity<std::function<void(const CellIteratorBaseType &,
                                               const unsigned int,
                                               ScratchData &,
                                               CopyData &)>>::type
      &boundary_worker = std::function<void(const CellIteratorBaseType &,
                                            const unsigned int,
                                            ScratchData &,
                                            CopyData &)>(),

    const typename identity<std::function<void(const CellIteratorBaseType &,
                                               const unsigned int,
                                               const unsigned int,
                                               const CellIteratorBaseType &,
                                               const unsigned int,
                                               const unsigned int,
                                               ScratchData &,
                                               CopyData &)>>::type
      &face_worker = std::function<void(const CellIteratorBaseType &,
                                        const unsigned int,
                                        const unsigned int,
                                        const CellIteratorBaseType &,
                                        const unsigned int,
                                        const unsigned int,
                                        ScratchData &,
                                        CopyData &)>(),

    const unsigned int queue_length = 2 * MultithreadInfo::n_threads(),
    const unsigned int chunk_size   = 8)
  {
    // Call the function above
    mesh_loop<typename IteratorRange<CellIteratorType>::IteratorOverIterators,
              ScratchData,
              CopyData,
              CellIteratorBaseType>(iterator_range.begin(),
                                    iterator_range.end(),
                                    cell_worker,
                                    copier,
                                    sample_scratch_data,
                                    sample_copy_data,
                                    flags,
                                    boundary_worker,
                                    face_worker,
                                    queue_length,
                                    chunk_size);
  }

  /**
   * This is a variant of the mesh_loop() function that can be used for worker
   * and copier functions that are member functions of a class.
   *
   * The argument passed as @p end must be convertible to the same type as @p
   * begin, but doesn't have to be of the same type itself. This allows to
   * write code like <code>mesh_loop(dof_handler.begin_active(),
   * dof_handler.end(), ...)</code> where the first is of type
   * DoFHandler::active_cell_iterator whereas the second is of type
   * DoFHandler::raw_cell_iterator.
   *
   * The @p queue_length argument indicates the number of items that can be
   * live at any given time. Each item consists of @p chunk_size elements of
   * the input stream that will be worked on by the worker and copier
   * functions one after the other on the same thread.
   *
   * @note If your data objects are large, or their constructors are
   * expensive, it is helpful to keep in mind that <tt>queue_length</tt>
   * copies of the <tt>ScratchData</tt> object and
   * <tt>queue_length*chunk_size</tt> copies of the <tt>CopyData</tt> object
   * are generated.
   *
   * An example usage of the function is given by
   * @code
   *
   * struct ScratchData;
   * struct CopyData;
   *
   * template <int dim, int spacedim>
   * class MyClass
   * {
   * public:
   *   void
   *   cell_worker(const CellIteratorType &cell, ScratchData &, CopyData &);
   *
   *   void
   *   copier(const CopyData &);
   *
   *   ...
   * };
   *
   * ...
   *
   * MyClass<dim, spacedim> my_class;
   * ScratchData            scratch;
   * CopyData               copy;
   *
   * mesh_loop(tria.begin_active(),
   *           tria.end(),
   *           my_class,
   *           &MyClass<dim, spacedim>::cell_worker,
   *           &MyClass<dim, spacedim>::copier,
   *           scratch,
   *           copy,
   *           assemble_own_cells);
   * @endcode
   *
   * @ingroup MeshWorker
   */
  template <class CellIteratorType,
            class ScratchData,
            class CopyData,
            class MainClass>
  void
  mesh_loop(const CellIteratorType &                         begin,
            const typename identity<CellIteratorType>::type &end,
            MainClass &                                      main_class,
            void (MainClass::*cell_worker)(const CellIteratorType &,
                                           ScratchData &,
                                           CopyData &),
            void (MainClass::*copier)(const CopyData &),
            const ScratchData & sample_scratch_data,
            const CopyData &    sample_copy_data,
            const AssembleFlags flags                      = assemble_own_cells,
            void (MainClass::*boundary_worker)(const CellIteratorType &,
                                               const unsigned int,
                                               ScratchData &,
                                               CopyData &) = nullptr,
            void (MainClass::*face_worker)(const CellIteratorType &,
                                           const unsigned int,
                                           const unsigned int,
                                           const CellIteratorType &,
                                           const unsigned int,
                                           const unsigned int,
                                           ScratchData &,
                                           CopyData &)     = nullptr,
            const unsigned int queue_length = 2 * MultithreadInfo::n_threads(),
            const unsigned int chunk_size   = 8)
  {
    std::function<void(const CellIteratorType &, ScratchData &, CopyData &)>
      f_cell_worker;

    std::function<void(
      const CellIteratorType &, const unsigned int, ScratchData &, CopyData &)>
      f_boundary_worker;

    std::function<void(const CellIteratorType &,
                       const unsigned int,
                       const unsigned int,
                       const CellIteratorType &,
                       const unsigned int,
                       const unsigned int,
                       ScratchData &,
                       CopyData &)>
      f_face_worker;

    if (cell_worker != nullptr)
      f_cell_worker = [&main_class,
                       cell_worker](const CellIteratorType &cell_iterator,
                                    ScratchData &           scratch_data,
                                    CopyData &              copy_data) {
        (main_class.*cell_worker)(cell_iterator, scratch_data, copy_data);
      };

    if (boundary_worker != nullptr)
      f_boundary_worker =
        [&main_class, boundary_worker](const CellIteratorType &cell_iterator,
                                       const unsigned int      face_no,
                                       ScratchData &           scratch_data,
                                       CopyData &              copy_data) {
          (main_class.*
           boundary_worker)(cell_iterator, face_no, scratch_data, copy_data);
        };

    if (face_worker != nullptr)
      f_face_worker = [&main_class,
                       face_worker](const CellIteratorType &cell_iterator_1,
                                    const unsigned int      face_index_1,
                                    const unsigned int      subface_index_1,
                                    const CellIteratorType &cell_iterator_2,
                                    const unsigned int      face_index_2,
                                    const unsigned int      subface_index_2,
                                    ScratchData &           scratch_data,
                                    CopyData &              copy_data) {
        (main_class.*face_worker)(cell_iterator_1,
                                  face_index_1,
                                  subface_index_1,
                                  cell_iterator_2,
                                  face_index_2,
                                  subface_index_2,
                                  scratch_data,
                                  copy_data);
      };

    mesh_loop(begin,
              end,
              f_cell_worker,
              [&main_class, copier](const CopyData &copy_data) {
                (main_class.*copier)(copy_data);
              },
              sample_scratch_data,
              sample_copy_data,
              flags,
              f_boundary_worker,
              f_face_worker,
              queue_length,
              chunk_size);
  }

  /**
   * Same as the function above, but for iterator ranges (and, therefore,
   * filtered iterators).
   *
   * An example usage of the function for the serial case is given by
   * @code
   *
   * struct ScratchData;
   * struct CopyData;
   *
   * template <int dim, int spacedim>
   * class MyClass
   * {
   * public:
   *   void
   *   cell_worker(const CellIteratorType &cell, ScratchData &, CopyData &);
   *
   *   void
   *   copier(const CopyData &);
   *
   *   ...
   * };
   *
   * ...
   *
   * MyClass<dim, spacedim> my_class;
   * ScratchData            scratch;
   * CopyData               copy;
   *
   * mesh_loop(tria.active_cell_iterators(),
   *           my_class,
   *           &MyClass<dim, spacedim>::cell_worker,
   *           &MyClass<dim, spacedim>::copier,
   *           scratch,
   *           copy,
   *           assemble_own_cells);
   * @endcode
   *
   * and an example usage of the function for the parallel distributed case,
   * where the copier is only to be called on locally owned cells, is given by
   * @code
   *
   * struct ScratchData;
   * struct CopyData;
   *
   * template <int dim, int spacedim>
   * class MyClass
   * {
   * public:
   *   void
   *   cell_worker(const CellIteratorType &cell, ScratchData &, CopyData &);
   *
   *   void
   *   copier(const CopyData &);
   *
   *   ...
   * };
   *
   * ...
   *
   * MyClass<dim, spacedim> my_class;
   * ScratchData            scratch;
   * CopyData               copy;
   *
   * const auto filtered_iterator_range =
   *   filter_iterators(distributed_tria.active_cell_iterators(),
   *                    IteratorFilters::LocallyOwnedCell());
   *
   * mesh_loop(filtered_iterator_range,
   *           my_class,
   *           &MyClass<dim, spacedim>::cell_worker,
   *           &MyClass<dim, spacedim>::copier,
   *           scratch,
   *           copy,
   *           assemble_own_cells);
   * @endcode
   *
   * @ingroup MeshWorker
   */
  template <class CellIteratorType,
            class ScratchData,
            class CopyData,
            class MainClass,
            class CellIteratorBaseType =
              typename internal::CellIteratorBaseType<CellIteratorType>::type>
  void
  mesh_loop(IteratorRange<CellIteratorType> iterator_range,
            MainClass &                     main_class,
            void (MainClass::*cell_worker)(const CellIteratorBaseType &,
                                           ScratchData &,
                                           CopyData &),
            void (MainClass::*copier)(const CopyData &),
            const ScratchData & sample_scratch_data,
            const CopyData &    sample_copy_data,
            const AssembleFlags flags                      = assemble_own_cells,
            void (MainClass::*boundary_worker)(const CellIteratorBaseType &,
                                               const unsigned int,
                                               ScratchData &,
                                               CopyData &) = nullptr,
            void (MainClass::*face_worker)(const CellIteratorBaseType &,
                                           const unsigned int,
                                           const unsigned int,
                                           const CellIteratorBaseType &,
                                           const unsigned int,
                                           const unsigned int,
                                           ScratchData &,
                                           CopyData &)     = nullptr,
            const unsigned int queue_length = 2 * MultithreadInfo::n_threads(),
            const unsigned int chunk_size   = 8)
  {
    // Call the function above
    mesh_loop<typename IteratorRange<CellIteratorType>::IteratorOverIterators,
              ScratchData,
              CopyData,
              MainClass,
              CellIteratorBaseType>(iterator_range.begin(),
                                    iterator_range.end(),
                                    main_class,
                                    cell_worker,
                                    copier,
                                    sample_scratch_data,
                                    sample_copy_data,
                                    flags,
                                    boundary_worker,
                                    face_worker,
                                    queue_length,
                                    chunk_size);
  }
} // namespace MeshWorker

DEAL_II_NAMESPACE_CLOSE

#endif
