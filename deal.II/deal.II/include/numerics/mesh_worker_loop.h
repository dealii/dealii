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

}



namespace MeshWorker
{
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
  template<class DOFINFO, class INFOBOX, class ITERATOR, typename ASSEMBLER>
  void loop(ITERATOR begin,
	    typename identity<ITERATOR>::type end,
	    DOFINFO cell_dof_info,
	    INFOBOX& info,
	    const std_cxx1x::function<void (DOFINFO&, typename INFOBOX::CellInfo &)> &cell_worker,
	    const std_cxx1x::function<void (DOFINFO&, typename INFOBOX::FaceInfo &)> &boundary_worker,
	    const std_cxx1x::function<void (DOFINFO&, DOFINFO&,
					    typename INFOBOX::FaceInfo &,
					    typename INFOBOX::FaceInfo &)> &face_worker,
	    ASSEMBLER &assembler,
	    bool cells_first = true)
  {
    const bool integrate_cell          = (cell_worker != 0);
    const bool integrate_boundary      = (boundary_worker != 0);
    const bool integrate_interior_face = (face_worker != 0);

    ;
    DOFINFO face_dof_info = cell_dof_info;
    DOFINFO neighbor_dof_info = cell_dof_info;
    assembler.initialize_info(cell_dof_info, false);
    assembler.initialize_info(face_dof_info, true);
    assembler.initialize_info(neighbor_dof_info, true);
    
				     // Loop over all cells
    for (ITERATOR cell = begin; cell != end; ++cell)
      {
					 // Execute this, if cells
					 // have to be dealt with
					 // before faces
	if (integrate_cell && cells_first)
	  {
	    cell_dof_info.reinit(cell);
	    info.cell.reinit(cell_dof_info);
	    cell_worker(cell_dof_info, info.cell);
	    assembler.assemble(cell_dof_info);
	  }

	if (integrate_interior_face || integrate_boundary)
	  for (unsigned int face_no=0; face_no < GeometryInfo<ITERATOR::AccessorType::Container::dimension>::faces_per_cell; ++face_no)
	    {
	      typename ITERATOR::AccessorType::Container::face_iterator face = cell->face(face_no);
	      if (cell->at_boundary(face_no))
		{
		  if (integrate_boundary)
		    {
		      cell_dof_info.reinit(cell, face, face_no);
		      info.boundary.reinit(cell_dof_info);
		      boundary_worker(cell_dof_info, info.boundary);
		      assembler.assemble(cell_dof_info);
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
		      
		      face_dof_info.reinit(cell, face, face_no);
		      info.face.reinit(face_dof_info);
		      neighbor_dof_info.reinit(neighbor, nface,
					       neighbor_face_no.first, neighbor_face_no.second);
		      info.subface.reinit(neighbor_dof_info);
						       // Neighbor
						       // first to
						       // conform to
						       // old version
		      face_worker(face_dof_info, neighbor_dof_info, info.face, info.subface);
		      assembler.assemble (face_dof_info, neighbor_dof_info);
		    }
		  else
		    {
						       // Neighbor is
						       // on same
						       // level, but
						       // only do this
						       // from one side.
		      if (neighbor < cell) continue;
		      
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
		      face_dof_info.reinit(cell, face, face_no);
		      info.face.reinit(face_dof_info);
		      neighbor_dof_info.reinit(neighbor, neighbor->face(neighbor_face_no),
					       neighbor_face_no);
		      info.neighbor.reinit(neighbor_dof_info);
		      face_worker(face_dof_info, neighbor_dof_info, info.face, info.neighbor);
		      assembler.assemble(face_dof_info, neighbor_dof_info);
		    }
		}
	    } // faces
					 // Execute this, if faces
					 // have to be handled first
	if (integrate_cell && !cells_first)
	  {
	    cell_dof_info.reinit(cell);
	    info.cell.reinit(cell_dof_info);
	    cell_worker(cell_dof_info, info.cell);
	    assembler.assemble (cell_dof_info);
	  }
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
