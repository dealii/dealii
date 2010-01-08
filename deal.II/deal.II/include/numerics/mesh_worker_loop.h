//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2006, 2007, 2008, 2009, 2010 by Guido Kanschat
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
#include <base/template_constraints.h>
#include <grid/tria.h>
#include <numerics/mesh_worker_info.h>


DEAL_II_NAMESPACE_OPEN

template <typename> class TriaActiveIterator;

namespace MeshWorker
{
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


/**
 * The main work function of this namespace. Its action consists of two
 * loops.
 *
 * First, a loop over all cells in the iterator range is performed, in
 * each step updating the CellInfo object, then calling
 * LocalWorker::cell() with this object.
 *
 * In the second loop, we work through all the faces of cells in the
 * iterator range. The functions LocalWorker::bdry() and
 * LocalWorker::face() are called for each boundary and interior face,
 * respectively. Unilaterally refined interior faces are handled
 * automatically by the loop.
 *
 * The extend of the second loop will be determined by the two control
 * variables LocalWorker::boundary_fluxes and
 * LocalWorker::interior_fluxes.
 *
 * @ingroup MeshWorker
 * @author Guido Kanschat, 2009
 */
  template<class ITERATOR, class ENDITERATOR, class CELLINFO, class FACEINFO, class LOCALWORKER>
  void loop(ITERATOR begin, ENDITERATOR end,
	    CELLINFO& cellinfo, FACEINFO& bdryinfo,
	    FACEINFO& faceinfo, FACEINFO& subfaceinfo, FACEINFO& ngbrinfo,
	    LOCALWORKER& localworker,
	    bool cells_first = true)
  {
    const_cast<Triangulation<ITERATOR::AccessorType::Container::dimension>&>(
      begin->get_triangulation()).clear_user_flags();

				     // Loop over all cells
    for (ITERATOR cell = begin; cell != end; ++cell)
      {
					 // Execute this, if cells
					 // have to be dealt with
					 // before faces
	if (cells_first)
	  {
	    cellinfo.reinit(cell);
	    localworker.cell(cellinfo);
	  }

	if (localworker.interior_fluxes || localworker.boundary_fluxes)
	  for (unsigned int face_no=0; face_no < GeometryInfo<ITERATOR::AccessorType::Container::dimension>::faces_per_cell; ++face_no)
	    {
	      typename ITERATOR::AccessorType::Container::face_iterator face = cell->face(face_no);
					       // Treat every face only once
	      if (face->user_flag_set ()) continue;

	      if (cell->at_boundary(face_no))
		{
		  if (localworker.boundary_fluxes)
		    {
		      bdryinfo.reinit(cell, face, face_no);
		      localworker.bdry(bdryinfo);
		    }
		}
	      else if (localworker.interior_fluxes)
		{
		  if (face->user_flag_set ()) continue;
		  face->set_user_flag ();
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
		      faceinfo.reinit(cell, face, face_no);
		      subfaceinfo.reinit(neighbor, nface, neighbor_face_no.first, neighbor_face_no.second);
						       // Neighbor
						       // first to
						       // conform to
						       // old version
		      localworker.face(faceinfo, subfaceinfo);
		    }
		  else
		    {
						       // Neighbor is
						       // on same
						       // level

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
		      faceinfo.reinit(cell, face, face_no);
		      ngbrinfo.reinit(neighbor, neighbor->face(neighbor_face_no),
				      neighbor_face_no);

		      localworker.face(faceinfo, ngbrinfo);
		      neighbor->face(neighbor_face_no)->set_user_flag ();
		    }

		}
	    } // faces
					 // Execute this, if faces
					 // have to be handled first
	if (!cells_first)
	  {
	    cellinfo.reinit(cell);
	    localworker.cell(cellinfo);
	  }
      }
  }

/**
 * Simplified interface for loop() if specialized for integration.
 *
 * @ingroup MeshWorker
 * @author Guido Kanschat, 2009
 */
  template<int dim, class ITERATOR, class LOCALWORKER>
  void integration_loop(ITERATOR begin,
			typename identity<ITERATOR>::type end,
			IntegrationInfoBox<dim, dim>& box,
			LOCALWORKER& localworker,
			bool cells_first = true)
  {
    loop(begin, end, box.cell_info, box.bdry_info,
	 box.face_info, box.subface_info, box.neighbor_info,
	 localworker, cells_first);
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
