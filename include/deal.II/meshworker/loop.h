// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2014 by the deal.II authors
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


#ifndef __deal2__mesh_worker_loop_h
#define __deal2__mesh_worker_loop_h

#include <deal.II/base/config.h>
#include <deal.II/base/std_cxx11/function.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/grid/tria.h>
#include <deal.II/meshworker/local_integrator.h>
#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>


#define DEAL_II_MESHWORKER_PARALLEL 1

DEAL_II_NAMESPACE_OPEN

template <typename> class TriaActiveIterator;
template <typename> class FilteredIterator;

namespace internal
{
  /**
   * Find out if an iterator supports inactive cells.
   */
  template <class DI>
  inline bool is_active_iterator(const DI &)
  {
    return false;
  }

  template <class ACCESSOR>
  inline bool is_active_iterator(const TriaActiveIterator<ACCESSOR> &)
  {
    return true;
  }

  template <class ACCESSOR>
  inline bool is_active_iterator(const FilteredIterator<TriaActiveIterator<ACCESSOR> > &)
  {
    return true;
  }

  template<int dim, class DOFINFO, class A>
  void assemble(const MeshWorker::DoFInfoBox<dim, DOFINFO> &dinfo, A *assembler)
  {
    dinfo.assemble(*assembler);
  }
}



namespace MeshWorker
{
  /**
   * Collection of parameters for execution of MeshWorker loops.
   */
  class LoopControl
  {
  public:

    /**
     * Constructor.
     */
    LoopControl()
      : own_cells(true), ghost_cells(false),
        faces_to_ghost(LoopControl::one), own_faces(LoopControl::one),
        cells_first(true)
    {
    }

    /**
     * Loop over cells owned by this process. Defaults to <code>true</code>.
     */
    bool own_cells;
    /**
     * Loop over cells not owned by this process. Defaults to <code>false</code>.
     */
    bool ghost_cells;

    enum FaceOption
    {
      never,
      one,
      both
    };

    /**
     * Loop over faces between a locally owned cell and a ghost cell:
     * - never: do not assembly these faces
     * - one: only one of the processes will assemble these faces (
     * from the finer side or the process with the lower mpi rank)
     * - both: both processes will assemble these faces
     * Note that these faces are never assembled from both sides on a single
     * process.
     * Default is one.
     */
    FaceOption faces_to_ghost;

    /**
     * Loop over faces between two locally owned cells:
     * - never: do not assemble face terms
     * - one: assemble once (always coming from the finer side)
     * - both: assemble each face twice (not implemented for hanging nodes!)
     * Default is one.
     */
    FaceOption own_faces;


    /**
     * Flag to determine if cells integrals should be done before or after
     * face integrals. Default is t
     */
    bool cells_first;
  };



  /**
   * The function called by loop() to perform the required actions on a
   * cell and its faces. The three functions <tt>cell_worker</tt>,
   * <tt>boundary_worker</tt> and <tt>face_worker</tt> are the same ones
   * handed to loop(). While there we only run the loop over all cells,
   * here, we do a single cell and, if necessary, its faces, interior
   * and boundary.
   *
   * Upon return, the DoFInfo objects in the DoFInfoBox are filled with
   * the data computed on the cell and each of the faces. Thus, after
   * the execution of this function, we are ready to call
   * DoFInfoBox::assemble() to distribute the local data into global
   * data.
   *
   * @param cell is the cell we work on
   * @param dof_info is the object into which local results are
   * entered. It is expected to have been set up for the right types of
   * data.
   * @param info is the object containing additional data only needed
   * for internal processing.
   * @param cell_worker defines the local action on each cell.
   * @param boundary_worker defines the local action on boundary faces
   * @param face_worker defines the local action on interior faces.
   * @param loop_control control structure to specify what actions should be performed.
   *
   * @ingroup MeshWorker
   * @author Guido Kanschat
   * @date 2010
   */
  template<class INFOBOX, class DOFINFO, int dim, int spacedim, class ITERATOR>
  void cell_action(
    ITERATOR cell,
    DoFInfoBox<dim, DOFINFO> &dof_info,
    INFOBOX &info,
    const std_cxx11::function<void (DOFINFO &, typename INFOBOX::CellInfo &)> &cell_worker,
    const std_cxx11::function<void (DOFINFO &, typename INFOBOX::CellInfo &)> &boundary_worker,
    const std_cxx11::function<void (DOFINFO &, DOFINFO &,
                                    typename INFOBOX::CellInfo &,
                                    typename INFOBOX::CellInfo &)> &face_worker,
    const LoopControl &loop_control)
  {
    const bool ignore_subdomain = (cell->get_triangulation().locally_owned_subdomain()
                                   == numbers::invalid_subdomain_id);

    types::subdomain_id csid = (cell->is_level_cell())
                               ? cell->level_subdomain_id()
                               : cell->subdomain_id();

    const bool own_cell = ignore_subdomain || (csid == cell->get_triangulation().locally_owned_subdomain());

    dof_info.reset();

    if ((!ignore_subdomain) && (csid == numbers::artificial_subdomain_id))
      return;

    dof_info.cell.reinit(cell);
    dof_info.cell_valid = true;

    const bool integrate_cell          = (cell_worker != 0);
    const bool integrate_boundary      = (boundary_worker != 0);
    const bool integrate_interior_face = (face_worker != 0);

    if (integrate_cell)
      info.cell.reinit(dof_info.cell);
    // Execute this, if cells
    // have to be dealt with
    // before faces
    if (integrate_cell && loop_control.cells_first &&
        ((loop_control.own_cells && own_cell) || (loop_control.ghost_cells && !own_cell)))
      cell_worker(dof_info.cell, info.cell);

    // Call the callback function in
    // the info box to do
    // computations between cell and
    // face action.
    info.post_cell(dof_info);

    if (integrate_interior_face || integrate_boundary)
      for (unsigned int face_no=0; face_no < GeometryInfo<ITERATOR::AccessorType::Container::dimension>::faces_per_cell; ++face_no)
        {
          typename ITERATOR::AccessorType::Container::face_iterator face = cell->face(face_no);
          if (cell->at_boundary(face_no))
            {
              // only integrate boundary faces of own cells
              if (integrate_boundary && own_cell)
                {
                  dof_info.interior_face_available[face_no] = true;
                  dof_info.interior[face_no].reinit(cell, face, face_no);
                  info.boundary.reinit(dof_info.interior[face_no]);
                  boundary_worker(dof_info.interior[face_no], info.boundary);
                }
            }
          else if (integrate_interior_face)
            {
              // Interior face
              TriaIterator<typename ITERATOR::AccessorType> neighbor = cell->neighbor(face_no);

              types::subdomain_id neighbid = numbers::artificial_subdomain_id;
              if (neighbor->is_level_cell())
                neighbid = neighbor->level_subdomain_id();
              //subdomain id is only valid for active cells
              else if (neighbor->active())
                neighbid = neighbor->subdomain_id();

              const bool own_neighbor = ignore_subdomain ||
                                        (neighbid == cell->get_triangulation().locally_owned_subdomain());

              // skip all faces between two ghost cells
              if (!own_cell && !own_neighbor)
                continue;

              // skip if the user doesn't want faces between own cells
              if (own_cell && own_neighbor && loop_control.own_faces==LoopControl::never)
                continue;

              // skip face to ghost
              if (own_cell != own_neighbor && loop_control.faces_to_ghost==LoopControl::never)
                continue;

              // Deal with
              // refinement edges
              // from the refined
              // side. Assuming
              // one-irregular
              // meshes, this
              // situation should
              // only occur if
              // both cells are
              // active.
              if (cell->neighbor_is_coarser(face_no))
                {
                  Assert(!cell->has_children(), ExcInternalError());
                  Assert(!neighbor->has_children(), ExcInternalError());

                  // skip if only one processor needs to assemble the face
                  // to a ghost cell and the fine cell is not ours.
                  if (!own_cell
                      && loop_control.faces_to_ghost == LoopControl::one)
                    continue;

                  const std::pair<unsigned int, unsigned int> neighbor_face_no
                    = cell->neighbor_of_coarser_neighbor(face_no);
                  const typename ITERATOR::AccessorType::Container::face_iterator nface
                    = neighbor->face(neighbor_face_no.first);

                  dof_info.interior_face_available[face_no] = true;
                  dof_info.exterior_face_available[face_no] = true;
                  dof_info.interior[face_no].reinit(cell, face, face_no);
                  info.face.reinit(dof_info.interior[face_no]);
                  dof_info.exterior[face_no].reinit(
                    neighbor, nface, neighbor_face_no.first, neighbor_face_no.second);
                  info.subface.reinit(dof_info.exterior[face_no]);

                  face_worker(dof_info.interior[face_no], dof_info.exterior[face_no],
                              info.face, info.subface);
                }
              else
                {
                  // If iterator is active and neighbor is refined, skip
                  // internal face.
                  if (internal::is_active_iterator(cell) && neighbor->has_children())
                    {
                      Assert(loop_control.own_faces != LoopControl::both, ExcMessage(
                               "Assembling from both sides for own_faces is not "
                               "supported with hanging nodes!"));
                      continue;
                    }

                  // Now neighbor is on same level, double-check this:
                  Assert(cell->level()==neighbor->level(), ExcInternalError());

                  // If we own both cells only do faces from one side (unless
                  // LoopControl says otherwise). Here, we rely on cell comparison
                  // that will look at cell->index().
                  if (own_cell && own_neighbor
                      && loop_control.own_faces == LoopControl::one
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
                      && loop_control.faces_to_ghost == LoopControl::one
                      && (neighbid < csid))
                    continue;

                  const unsigned int neighbor_face_no = cell->neighbor_face_no(face_no);
                  Assert (neighbor->face(neighbor_face_no) == face, ExcInternalError());
                  // Regular interior face
                  dof_info.interior_face_available[face_no] = true;
                  dof_info.exterior_face_available[face_no] = true;
                  dof_info.interior[face_no].reinit(cell, face, face_no);
                  info.face.reinit(dof_info.interior[face_no]);
                  dof_info.exterior[face_no].reinit(
                    neighbor, neighbor->face(neighbor_face_no), neighbor_face_no);
                  info.neighbor.reinit(dof_info.exterior[face_no]);

                  face_worker(dof_info.interior[face_no], dof_info.exterior[face_no],
                              info.face, info.neighbor);
                }
            }
        } // faces
    // Call the callback function in
    // the info box to do
    // computations between face and
    // cell action.
    info.post_faces(dof_info);

    // Execute this, if faces
    // have to be handled first
    if (integrate_cell && !loop_control.cells_first &&
        ((loop_control.own_cells && own_cell) || (loop_control.ghost_cells && !own_cell)))
      cell_worker(dof_info.cell, info.cell);
  }


  /**
   * The main work function of this namespace. It is a loop over all
   * cells in an iterator range, in which cell_action() is called for
   * each cell. Unilaterally refined interior faces are handled
   * automatically by the loop.
   * Most of the work in this loop is done in cell_action(), which also
   * receives most of the parameters of this function. See the
   * documentation there for more details.
   *
   * If you don't want anything to be done on cells, interior or boundary faces
   * to happen, simply pass the Null pointer to one of the function
   * arguments.
   *
   * @ingroup MeshWorker
   * @author Guido Kanschat, 2009
   */
  template<int dim, int spacedim, class DOFINFO, class INFOBOX, class ASSEMBLER, class ITERATOR>
  void loop(ITERATOR begin,
            typename identity<ITERATOR>::type end,
            DOFINFO &dinfo,
            INFOBOX &info,
            const std_cxx11::function<void (DOFINFO &, typename INFOBOX::CellInfo &)> &cell_worker,
            const std_cxx11::function<void (DOFINFO &, typename INFOBOX::CellInfo &)> &boundary_worker,
            const std_cxx11::function<void (DOFINFO &, DOFINFO &,
                                            typename INFOBOX::CellInfo &,
                                            typename INFOBOX::CellInfo &)> &face_worker,
            ASSEMBLER &assembler,
            const LoopControl &lctrl = LoopControl())
  {
    DoFInfoBox<dim, DOFINFO> dof_info(dinfo);

    assembler.initialize_info(dof_info.cell, false);
    for (unsigned int i=0; i<GeometryInfo<dim>::faces_per_cell; ++i)
      {
        assembler.initialize_info(dof_info.interior[i], true);
        assembler.initialize_info(dof_info.exterior[i], true);
      }

    // Loop over all cells
#ifdef DEAL_II_MESHWORKER_PARALLEL
    WorkStream::run(begin, end,
                    std_cxx11::bind(&cell_action<INFOBOX, DOFINFO, dim, spacedim, ITERATOR>,
                                    std_cxx11::_1, std_cxx11::_3, std_cxx11::_2,
                                    cell_worker, boundary_worker, face_worker, lctrl),
                    std_cxx11::bind(&internal::assemble<dim,DOFINFO,ASSEMBLER>, std_cxx11::_1, &assembler),
                    info, dof_info);
#else
    for (ITERATOR cell = begin; cell != end; ++cell)
      {
        cell_action<INFOBOX,DOFINFO,dim,spacedim>(cell, dof_info,
                                                  info, cell_worker,
                                                  boundary_worker, face_worker,
                                                  lctrl);
        dof_info.assemble(assembler);
      }
#endif
  }

  template<int dim, int spacedim, class DOFINFO, class INFOBOX, class ASSEMBLER, class ITERATOR>
  void loop(ITERATOR begin,
            typename identity<ITERATOR>::type end,
            DOFINFO &dinfo,
            INFOBOX &info,
            const std_cxx11::function<void (DOFINFO &, typename INFOBOX::CellInfo &)> &cell_worker,
            const std_cxx11::function<void (DOFINFO &, typename INFOBOX::CellInfo &)> &boundary_worker,
            const std_cxx11::function<void (DOFINFO &, DOFINFO &,
                                            typename INFOBOX::CellInfo &,
                                            typename INFOBOX::CellInfo &)> &face_worker,
            ASSEMBLER &assembler,
            bool cells_first,
            bool unique_faces_only = true) DEAL_II_DEPRECATED;

  template<int dim, int spacedim, class DOFINFO, class INFOBOX, class ASSEMBLER, class ITERATOR>
  void loop(ITERATOR begin,
            typename identity<ITERATOR>::type end,
            DOFINFO &dinfo,
            INFOBOX &info,
            const std_cxx11::function<void (DOFINFO &, typename INFOBOX::CellInfo &)> &cell_worker,
            const std_cxx11::function<void (DOFINFO &, typename INFOBOX::CellInfo &)> &boundary_worker,
            const std_cxx11::function<void (DOFINFO &, DOFINFO &,
                                            typename INFOBOX::CellInfo &,
                                            typename INFOBOX::CellInfo &)> &face_worker,
            ASSEMBLER &assembler,
            bool cells_first,
            bool unique_faces_only)
  {
    LoopControl lctrl;
    lctrl.cells_first = cells_first;
    lctrl.own_faces = (unique_faces_only)
                      ? LoopControl::one
                      : LoopControl::both;

    loop<dim,spacedim>(begin, end, dinfo, info, cell_worker, boundary_worker, face_worker, assembler, lctrl);
  }

  /**
   * @deprecated The simplification in this loop is
   * insignificant. Therefore, it is recommended to use loop() instead.
   *
   * Simplified interface for loop() if specialized for integration.
   *
   * @ingroup MeshWorker
   * @author Guido Kanschat, 2009
   */
  template<int dim, int spacedim, class ITERATOR, class ASSEMBLER>
  void integration_loop(ITERATOR begin,
                        typename identity<ITERATOR>::type end,
                        DoFInfo<dim, spacedim> &dof_info,
                        IntegrationInfoBox<dim, spacedim> &box,
                        const std_cxx11::function<void (DoFInfo<dim>&, IntegrationInfo<dim, spacedim>&)> &cell_worker,
                        const std_cxx11::function<void (DoFInfo<dim>&, IntegrationInfo<dim, spacedim>&)> &boundary_worker,
                        const std_cxx11::function<void (DoFInfo<dim> &, DoFInfo<dim> &,
                                                        IntegrationInfo<dim, spacedim> &,
                                                        IntegrationInfo<dim, spacedim> &)> &face_worker,
                        ASSEMBLER &assembler,
                        bool cells_first) DEAL_II_DEPRECATED;


  template<int dim, int spacedim, class ITERATOR, class ASSEMBLER>
  void integration_loop(ITERATOR begin,
                        typename identity<ITERATOR>::type end,
                        DoFInfo<dim, spacedim> &dof_info,
                        IntegrationInfoBox<dim, spacedim> &box,
                        const std_cxx11::function<void (DoFInfo<dim>&, IntegrationInfo<dim, spacedim>&)> &cell_worker,
                        const std_cxx11::function<void (DoFInfo<dim>&, IntegrationInfo<dim, spacedim>&)> &boundary_worker,
                        const std_cxx11::function<void (DoFInfo<dim> &, DoFInfo<dim> &,
                                                        IntegrationInfo<dim, spacedim> &,
                                                        IntegrationInfo<dim, spacedim> &)> &face_worker,
                        ASSEMBLER &assembler,
                        bool cells_first)
  {
    LoopControl lctrl;

    loop<dim, spacedim>
    (begin, end,
     dof_info,
     box,
     cell_worker,
     boundary_worker,
     face_worker,
     assembler,
     cells_first,
     lctrl);
  }


  /**
   * Simplified interface for loop() if specialized for integration,
   * using the virtual functions in LocalIntegrator.
   *
   * @ingroup MeshWorker
   * @author Guido Kanschat, 2009
   */
  template<int dim, int spacedim, class ITERATOR, class ASSEMBLER>
  void integration_loop(ITERATOR begin,
                        typename identity<ITERATOR>::type end,
                        DoFInfo<dim, spacedim> &dof_info,
                        IntegrationInfoBox<dim, spacedim> &box,
                        const LocalIntegrator<dim, spacedim> &integrator,
                        ASSEMBLER &assembler,
                        bool cells_first)
  DEAL_II_DEPRECATED;

  template<int dim, int spacedim, class ITERATOR, class ASSEMBLER>
  void integration_loop(ITERATOR begin,
                        typename identity<ITERATOR>::type end,
                        DoFInfo<dim, spacedim> &dof_info,
                        IntegrationInfoBox<dim, spacedim> &box,
                        const LocalIntegrator<dim, spacedim> &integrator,
                        ASSEMBLER &assembler,
                        bool cells_first)
  {
    std_cxx11::function<void (DoFInfo<dim>&, IntegrationInfo<dim, spacedim>&)> cell_worker;
    std_cxx11::function<void (DoFInfo<dim>&, IntegrationInfo<dim, spacedim>&)> boundary_worker;
    std_cxx11::function<void (DoFInfo<dim> &, DoFInfo<dim> &,
                              IntegrationInfo<dim, spacedim> &,
                              IntegrationInfo<dim, spacedim> &)> face_worker;
    if (integrator.use_cell)
      cell_worker = std_cxx11::bind(&LocalIntegrator<dim, spacedim>::cell, &integrator, std_cxx11::_1, std_cxx11::_2);
    if (integrator.use_boundary)
      boundary_worker = std_cxx11::bind(&LocalIntegrator<dim, spacedim>::boundary, &integrator, std_cxx11::_1, std_cxx11::_2);
    if (integrator.use_face)
      face_worker = std_cxx11::bind(&LocalIntegrator<dim, spacedim>::face, &integrator, std_cxx11::_1, std_cxx11::_2, std_cxx11::_3, std_cxx11::_4);

    loop<dim, spacedim>
    (begin, end,
     dof_info,
     box,
     cell_worker,
     boundary_worker,
     face_worker,
     assembler,
     cells_first);
  }

  /**
   * As above but using LoopControl
   */
  template<int dim, int spacedim, class ITERATOR, class ASSEMBLER>
  void integration_loop(ITERATOR begin,
                        typename identity<ITERATOR>::type end,
                        DoFInfo<dim, spacedim> &dof_info,
                        IntegrationInfoBox<dim, spacedim> &box,
                        const LocalIntegrator<dim, spacedim> &integrator,
                        ASSEMBLER &assembler,
                        const LoopControl &lctrl = LoopControl())
  {
    std_cxx11::function<void (DoFInfo<dim>&, IntegrationInfo<dim, spacedim>&)> cell_worker;
    std_cxx11::function<void (DoFInfo<dim>&, IntegrationInfo<dim, spacedim>&)> boundary_worker;
    std_cxx11::function<void (DoFInfo<dim> &, DoFInfo<dim> &,
                              IntegrationInfo<dim, spacedim> &,
                              IntegrationInfo<dim, spacedim> &)> face_worker;
    if (integrator.use_cell)
      cell_worker = std_cxx11::bind(&LocalIntegrator<dim, spacedim>::cell, &integrator, std_cxx11::_1, std_cxx11::_2);
    if (integrator.use_boundary)
      boundary_worker = std_cxx11::bind(&LocalIntegrator<dim, spacedim>::boundary, &integrator, std_cxx11::_1, std_cxx11::_2);
    if (integrator.use_face)
      face_worker = std_cxx11::bind(&LocalIntegrator<dim, spacedim>::face, &integrator, std_cxx11::_1, std_cxx11::_2, std_cxx11::_3, std_cxx11::_4);

    loop<dim, spacedim>
    (begin, end,
     dof_info,
     box,
     cell_worker,
     boundary_worker,
     face_worker,
     assembler,
     lctrl);
  }

}

DEAL_II_NAMESPACE_CLOSE

#endif
