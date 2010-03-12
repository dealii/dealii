//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2006, 2007, 2008, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef __deal2__mesh_worker_loop_h
#define __deal2__mesh_worker_loop_h

#include <base/config.h>
#include <base/std_cxx1x/function.h>
#include <base/work_stream.h>
#include <base/template_constraints.h>
#include <grid/tria.h>
#include <numerics/mesh_worker_info.h>


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

  template<int dim, int sdim, class A>
  void assemble(const MeshWorker::DoFInfoBox<dim, sdim>& dinfo, A& assembler)
  {
    dinfo.assemble(assembler);
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
 * @param cells_first determines, whether faces or cells are to be
 * dealt with first.
 * @param unique_faces_only determines, that a face between two cells
 * of the same level is processed only from the cell which is less
 * than its neighbor. If this parameter is <tt>false</tt> these faces
 * are processed from both cells.
 *
 * @author Guido Kanschat, 2010
 */
  template<class INFOBOX, int dim, int spacedim, class ITERATOR>
  void cell_action(
    ITERATOR cell,
    DoFInfoBox<dim, spacedim>& dof_info,
    INFOBOX& info,
    const std_cxx1x::function<void (DoFInfo<dim, spacedim>&, typename INFOBOX::CellInfo &)> &cell_worker,
    const std_cxx1x::function<void (DoFInfo<dim, spacedim>&, typename INFOBOX::FaceInfo &)> &boundary_worker,
    const std_cxx1x::function<void (DoFInfo<dim, spacedim>&, DoFInfo<dim, spacedim>&,
				    typename INFOBOX::FaceInfo &,
				    typename INFOBOX::FaceInfo &)> &face_worker,
    const bool cells_first,
    const bool unique_faces_only)
  {
    const bool integrate_cell          = (cell_worker != 0);
    const bool integrate_boundary      = (boundary_worker != 0);
    const bool integrate_interior_face = (face_worker != 0);
    
    dof_info.reset();
				     // Execute this, if cells
				     // have to be dealt with
				     // before faces
    if (integrate_cell && cells_first)
      {
	dof_info.cell.reinit(cell);
	info.cell.reinit(dof_info.cell);
	cell_worker(dof_info.cell, info.cell);
      }
    
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
				     // Execute this, if faces
				     // have to be handled first
    if (integrate_cell && !cells_first)
      {
	dof_info.cell.reinit(cell);
	info.cell.reinit(dof_info.cell);
	cell_worker(dof_info.cell, info.cell);
      }  
  }
  
  
/**
 * The main work function of this namespace. Its action consists of two
 * loops.
 *
 * First, a loop over all cells in the iterator range is performed, in
 * each step updating the CellInfo object, then calling
 * cell_worker() with this object.
 *
 * In the second loop, we walk through all the faces of cells in the
 * iterator range. The functions boundary_worker() and
 * face_worker() are called for each boundary and interior face,
 * respectively. Unilaterally refined interior faces are handled
 * automatically by the loop.
 *
 * If you don't want integration on cells, interior or boundary faces
 * to happen, simply pass the Null pointer to one of the function
 * arguments.
 *
 * @ingroup MeshWorker
 * @author Guido Kanschat, 2009
 */
  template<class INFOBOX, int dim, int spacedim, typename ASSEMBLER, class ITERATOR>
  void loop(ITERATOR begin,
	    typename identity<ITERATOR>::type end,
	    DoFInfo<dim, spacedim> dinfo,
	    INFOBOX& info,
	    const std_cxx1x::function<void (DoFInfo<dim, spacedim>&, typename INFOBOX::CellInfo &)> &cell_worker,
	    const std_cxx1x::function<void (DoFInfo<dim, spacedim>&, typename INFOBOX::FaceInfo &)> &boundary_worker,
	    const std_cxx1x::function<void (DoFInfo<dim, spacedim>&, DoFInfo<dim, spacedim>&,
					    typename INFOBOX::FaceInfo &,
					    typename INFOBOX::FaceInfo &)> &face_worker,
	    ASSEMBLER &assembler,
	    bool cells_first = true)
  {
    DoFInfoBox<dim, spacedim> dof_info(dinfo);
    
    assembler.initialize_info(dof_info.cell, false);
    for (unsigned int i=0;i<GeometryInfo<dim>::faces_per_cell;++i)
      {
	assembler.initialize_info(dof_info.interior[i], true);
	assembler.initialize_info(dof_info.exterior[i], true);
      }
    
				     // Loop over all cells
//     WorkStream::run(begin, end,
// 		    std_cxx1x::bind(cell_action<INFOBOX, dim, spacedim, ITERATOR>, _1, _3, _2,
// 				    cell_worker, boundary_worker, face_worker, cells_first, true),
// 		    std_cxx1x::bind(internal::assemble<dim,spacedim,ASSEMBLER>, _1, assembler),
// 		    info, dof_info);
    
    for (ITERATOR cell = begin; cell != end; ++cell)
      {
	cell_action(cell, dof_info, info, cell_worker, boundary_worker, face_worker, cells_first, true);
	dof_info.assemble(assembler);
      }
  }

/**
 * Simplified interface for loop() if specialized for integration.
 *
 * @ingroup MeshWorker
 * @author Guido Kanschat, 2009
 */
  template<class CELLINFO, class FACEINFO, int dim, class ITERATOR, typename ASSEMBLER>
  void integration_loop(ITERATOR begin,
			typename identity<ITERATOR>::type end,
			DoFInfo<dim>& dof_info,
			IntegrationInfoBox<dim, dim>& box,
			const std_cxx1x::function<void (DoFInfo<dim>&, CELLINFO &)> &cell_worker,
			const std_cxx1x::function<void (DoFInfo<dim>&, FACEINFO &)> &boundary_worker,
			const std_cxx1x::function<void (DoFInfo<dim>&, DoFInfo<dim>&, FACEINFO &, FACEINFO &)> &face_worker,
			ASSEMBLER &assembler,
			bool cells_first = true)
  {
    loop<DoFInfo<dim>,IntegrationInfoBox<dim, dim> >
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
