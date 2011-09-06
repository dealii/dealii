//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2006, 2007, 2008, 2009, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef __deal2__mesh_worker_loop_h
#define __deal2__mesh_worker_loop_h

#include <deal.II/base/config.h>
#include <deal.II/base/std_cxx1x/function.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/grid/tria.h>
#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>

#define DEAL_II_MESHWORKER_PARALLEL 1

DEAL_II_NAMESPACE_OPEN

template <typename> class TriaActiveIterator;

namespace internal
{
/**
 * Find out if an iterator supports inactive cells.
 */
  template <class DI>
  inline bool is_active_iterator(const DI&)
  {
    return false;
  }

  template <class ACCESSOR>
  inline bool is_active_iterator(const TriaActiveIterator<ACCESSOR>&)
  {
    return true;
  }

  template<int dim, class DOFINFO, class A>
  void assemble(const MeshWorker::DoFInfoBox<dim, DOFINFO>& dinfo, A* assembler)
  {
    dinfo.assemble(*assembler);
  }
}



namespace MeshWorker
{
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
 * @param cells_first determines, whether, on a given cell, face or cell
 *        integrals are to be  dealt with first. Note that independent of the
 *        value of this flag, cell and face integrals of a given cell are
 *        all taken care of before moving to the next cell.
 * @param unique_faces_only determines, that a face between two cells
 * of the same level is processed only from the cell which is less
 * than its neighbor. If this parameter is <tt>false</tt> these faces
 * are processed from both cells.
 *
 * @ingroup MeshWorker
 * @author Guido Kanschat
 * @date 2010
 */
  template<class INFOBOX, class DOFINFO, int dim, int spacedim, class ITERATOR>
  void cell_action(
    ITERATOR cell,
    DoFInfoBox<dim, DOFINFO>& dof_info,
    INFOBOX& info,
    const std_cxx1x::function<void (DOFINFO&, typename INFOBOX::CellInfo&)>& cell_worker,
    const std_cxx1x::function<void (DOFINFO&, typename INFOBOX::CellInfo&)>& boundary_worker,
    const std_cxx1x::function<void (DOFINFO&, DOFINFO&,
				    typename INFOBOX::CellInfo&,
				    typename INFOBOX::CellInfo&)>& face_worker,
    const bool cells_first,
    const bool unique_faces_only)
  {
    const bool integrate_cell          = (cell_worker != 0);
    const bool integrate_boundary      = (boundary_worker != 0);
    const bool integrate_interior_face = (face_worker != 0);

    dof_info.reset();

    dof_info.cell.reinit(cell);
    if (integrate_cell)
      info.cell.reinit(dof_info.cell);
				     // Execute this, if cells
				     // have to be dealt with
				     // before faces
    if (integrate_cell && cells_first)
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
	      if (integrate_boundary)
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
	      typename ITERATOR::AccessorType::Container::cell_iterator
		neighbor = cell->neighbor(face_no);

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
	      if (neighbor->level() < cell->level())
		{
		  Assert(!cell->has_children(), ExcInternalError());
		  Assert(!neighbor->has_children(), ExcInternalError());

		  std::pair<unsigned int, unsigned int> neighbor_face_no
		    = cell->neighbor_of_coarser_neighbor(face_no);
		  typename ITERATOR::AccessorType::Container::face_iterator nface
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
						   // Neighbor is
						   // on same
						   // level, but
						   // only do this
						   // from one side.
		  if (unique_faces_only && (neighbor < cell)) continue;

						   // If iterator
						   // is active
						   // and neighbor
						   // is refined,
						   // skip
						   // internal face.
		  if (internal::is_active_iterator(cell) && neighbor->has_children())
		    continue;

		  unsigned int neighbor_face_no = cell->neighbor_of_neighbor(face_no);
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
    if (integrate_cell && !cells_first)
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
	    DOFINFO& dinfo,
	    INFOBOX& info,
	    const std_cxx1x::function<void (DOFINFO&, typename INFOBOX::CellInfo&)>& cell_worker,
	    const std_cxx1x::function<void (DOFINFO&, typename INFOBOX::CellInfo&)>& boundary_worker,
	    const std_cxx1x::function<void (DOFINFO&, DOFINFO&,
					    typename INFOBOX::CellInfo&,
					    typename INFOBOX::CellInfo&)>& face_worker,
	    ASSEMBLER& assembler,
	    bool cells_first = true)
  {
    DoFInfoBox<dim, DOFINFO> dof_info(dinfo);

    assembler.initialize_info(dof_info.cell, false);
    for (unsigned int i=0;i<GeometryInfo<dim>::faces_per_cell;++i)
      {
	assembler.initialize_info(dof_info.interior[i], true);
	assembler.initialize_info(dof_info.exterior[i], true);
      }

				     // Loop over all cells
#ifdef DEAL_II_MESHWORKER_PARALLEL
    WorkStream::run(begin, end,
 		    std_cxx1x::bind(&cell_action<INFOBOX, DOFINFO, dim, spacedim, ITERATOR>, 
				    std_cxx1x::_1, std_cxx1x::_3, std_cxx1x::_2,
				    cell_worker, boundary_worker, face_worker, cells_first, true),
 		    std_cxx1x::bind(&internal::assemble<dim,DOFINFO,ASSEMBLER>, std_cxx1x::_1, &assembler),
 		    info, dof_info);
#else
    for (ITERATOR cell = begin; cell != end; ++cell)
      {
	cell_action<INFOBOX,DOFINFO,dim,spacedim>(cell, dof_info,
						  info, cell_worker,
						  boundary_worker, face_worker,
						  cells_first,
						  true);
	dof_info.assemble(assembler);
      }
#endif
  }

/**
 * Simplified interface for loop() if specialized for integration.
 *
 * @ingroup MeshWorker
 * @author Guido Kanschat, 2009
 */
  template<int dim, int spacedim, class ITERATOR, class ASSEMBLER>
  void integration_loop(ITERATOR begin,
			typename identity<ITERATOR>::type end,
			DoFInfo<dim, spacedim>& dof_info,
			IntegrationInfoBox<dim, spacedim>& box,
			const std_cxx1x::function<void (DoFInfo<dim>&, IntegrationInfo<dim, spacedim>&)> &cell_worker,
			const std_cxx1x::function<void (DoFInfo<dim>&, IntegrationInfo<dim, spacedim>&)> &boundary_worker,
			const std_cxx1x::function<void (DoFInfo<dim>&, DoFInfo<dim>&,
							IntegrationInfo<dim, spacedim>&,
							IntegrationInfo<dim, spacedim>&)> &face_worker,
			ASSEMBLER &assembler,
			bool cells_first = true)
  {
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
}

DEAL_II_NAMESPACE_CLOSE

#endif
