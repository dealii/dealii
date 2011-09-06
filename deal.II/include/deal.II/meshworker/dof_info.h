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

#ifndef __deal2__mesh_worker_dof_info_h
#define __deal2__mesh_worker_dof_info_h

#include <deal.II/base/config.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/std_cxx1x/shared_ptr.h>
#include <deal.II/dofs/block_info.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/meshworker/local_results.h>
#include <deal.II/meshworker/vector_selector.h>

DEAL_II_NAMESPACE_OPEN

namespace MeshWorker
{
  template <int dim, class DOFINFO> class DoFInfoBox;
/**
 * A class containing information on geometry and degrees of freedom
 * of a mesh object.
 *
 * The information in these objects is usually used by one of the
 * Assembler classes. It is also the kind of information which is
 * needed in mesh based matrices (often referred to as matrix free
 * methods).
 *
 * In addition to the information on degrees of freedom stored in this
 * class, it also provides the local computation space for the worker
 * object operating on it in LocalResults. This base class will automatically
 * reinitialized on each cell, but initial setup is up to the user and
 * should be done when initialize() for this class is called.
 *
 * This class operates in two different modes, corresponding to the
 * data models discussed in the Assembler namespace documentation.
 *
 * The choice of the local data model is triggered by the vector
 * BlockInfo::local_renumbering, which in turn is usually filled by
 * BlockInfo::initialize_local(). If this function has been used, or
 * the vector has been changed from zero-length, then local dof
 * indices stored in this object will automatically be renumbered to
 * reflect local block structure. This means, the first entries in
 * #indices will refer to the first block of the system, then comes
 * the second block and so on.
 *
 * The BlockInfo object is stored as a pointer. Therefore, if the
 * block structure changes, for instance because of mesh refinement,
 * the DoFInfo class will automatically use the new structures.
 *
 * @ingroup MeshWorker
 * @author Guido Kanschat, 2009
 */
  template<int dim, int spacedim = dim, typename number = double>
  class DoFInfo : public LocalResults<number>
  {
    public:
				       /// The current cell
      typename Triangulation<dim, spacedim>::cell_iterator cell;

				       /// The current face
      typename Triangulation<dim, spacedim>::face_iterator face;

				       /**
					* The number of the current
					* face on the current cell.
					*
					* This number is
					* deal_II_numbers::invalid_unsigned_int
					* if the info object was
					* initialized with a cell.
					*/

      unsigned int face_number;
				       /**
					* The number of the current
					* subface on the current
					* face
					*
					* This number is
					* deal_II_numbers::invalid_unsigned_int
					* if the info object was not
					* initialized with a subface.
					*/
      unsigned int sub_number;

				       /// The DoF indices of the current cell
      std::vector<unsigned int> indices;

				       /**
					* Constructor setting the
					* #block_info pointer.
					*/
      DoFInfo(const BlockInfo& block_info);

				       /**
					* Constructor
					* leaving the #block_info
					* pointer empty, but setting
					* the #aux_local_indices.
					*/
      template <class DH>
      DoFInfo (const DH& dof_handler);

				       /**
					* Set the current cell and
					* fill #indices.
					*/
      template <class DHCellIterator>
      void reinit(const DHCellIterator& c);

				       /**
					* Set the current face and
					* fill #indices if the #cell
					* changed.
					*/
      template <class DHCellIterator, class DHFaceIterator>
      void reinit(const DHCellIterator& c,
		  const DHFaceIterator& f,
		  const unsigned int n);

				       /**
					* Set the current subface
					* and fill #indices if the
					* #cell changed.
					*/
      template <class DHCellIterator, class DHFaceIterator>
      void reinit(const DHCellIterator& c,
		  const DHFaceIterator& f,
		  const unsigned int n,
		  const unsigned int s);

				       /**
					* Switch to a new face of the
					* same cell. Does not change
					* #indices and does not reset
					* data in LocalResults.
					*/
      template <class DHFaceIterator>
      void set_face (const DHFaceIterator& f,
		     const unsigned int n);
				       /**
					* Switch to a new subface of the
					* same cell. Does not change
					* #indices and does not reset
					* data in LocalResults.
					*/
      template <class DHFaceIterator>
      void set_subface (const DHFaceIterator& f,
			const unsigned int n,
			const unsigned int s);

      const BlockIndices& local_indices() const;


				       /// The block structure of the system
      SmartPointer<const BlockInfo,DoFInfo<dim,spacedim> > block_info;

      bool level_cell;
    private:
				       /**
					* Standard constructor, not
					* setting any block
					* indices. Use of this
					* constructor is not
					* recommended, but it is
					* needed for the arrays in
					* DoFInfoBox.
					*/
      DoFInfo ();

      				       /// Fill index vector with active indices
      void get_indices(const typename DoFHandler<dim, spacedim>::cell_iterator& c);

				       /// Fill index vector with level indices
      void get_indices(const typename MGDoFHandler<dim, spacedim>::cell_iterator& c);

				       /// Auxiliary vector
      std::vector<unsigned int> indices_org;

				       /**
					* An auxiliary local
					* BlockIndices object created
					* if #block_info is not set.
					* It contains just a single
					* block of the size of
					* degrees of freedom per cell.
					*/
      BlockIndices aux_local_indices;

      friend class DoFInfoBox<dim, DoFInfo<dim, spacedim, number> >;
  };


  /**
 * A class bundling the MeshWorker::DoFInfo objects used on a cell.
 *
 * @todo Currently, we are storing an object for the cells and two for
 * each face. We could gather all face data pertaining to the cell
 * itself in one object, saving a bit of memory and a few operations,
 * but sacrificing some cleanliness.
 *
 * @ingroup MeshWorker
 * @author Guido Kanschat, 2010
 */
  template <int dim, class DOFINFO>
  class DoFInfoBox
  {
    public:
				       /**
					* Constructor copying the seed
					* into all other objects.
					*/
      DoFInfoBox(const DOFINFO& seed);

				       /**
					* Copy constructor, taking
					* #cell and using it as a seed
					* in the other constructor.
					*/
      DoFInfoBox(const DoFInfoBox<dim, DOFINFO>&);

				       /**
					* Reset all the availability flags.
					*/
      void reset();

				       /**
					* After all info objects have
					* been filled appropriately,
					* use the ASSEMBLER object
					* to assemble them into the
					* global data. See
					* MeshWorker::Assembler for
					* available classes.
					*/
      template <class ASSEMBLER>
      void assemble(ASSEMBLER& ass) const;

				       /**
					* The memory used by this object.
					*/
      std::size_t memory_consumption () const;
      

				       /**
					* The data for the cell.
					*/
      DOFINFO cell;
				       /**
					* The data for the faces from inside.
					*/
      DOFINFO interior[GeometryInfo<dim>::faces_per_cell];
				       /**
					* The data for the faces from outside.
					*/
      DOFINFO exterior[GeometryInfo<dim>::faces_per_cell];

				       /**
					* A set of flags, indicating
					* whether data on an interior
					* face is available.
					*/
      bool interior_face_available[GeometryInfo<dim>::faces_per_cell];
				       /**
					* A set of flags, indicating
					* whether data on an exterior
					* face is available.
					*/
      bool exterior_face_available[GeometryInfo<dim>::faces_per_cell];
  };

//----------------------------------------------------------------------//

  template <int dim, int spacedim, typename number>
  template <class DH>
  DoFInfo<dim,spacedim,number>::DoFInfo(const DH& dof_handler)
  {
    std::vector<unsigned int> aux(1);
    aux[0] = dof_handler.get_fe().dofs_per_cell;
    aux_local_indices.reinit(aux);
  }


  template <int dim, int spacedim, typename number>
  template <class DHCellIterator>
  inline void
  DoFInfo<dim,spacedim,number>::reinit(const DHCellIterator& c)
  {
    get_indices(c);
    cell = static_cast<typename Triangulation<dim,spacedim>::cell_iterator> (c);
    face_number = deal_II_numbers::invalid_unsigned_int;
    sub_number = deal_II_numbers::invalid_unsigned_int;
    if (block_info)
      LocalResults<number>::reinit(block_info->local());
    else
      LocalResults<number>::reinit(aux_local_indices);
  }


  template<int dim, int spacedim, typename number>
  template <class DHFaceIterator>
  inline void
  DoFInfo<dim,spacedim,number>::set_face(
    const DHFaceIterator& f,
    unsigned int n)
  {
    face = static_cast<typename Triangulation<dim>::face_iterator> (f);
    face_number = n;
    sub_number = deal_II_numbers::invalid_unsigned_int;
  }
  
  
  template<int dim, int spacedim, typename number>
  template <class DHCellIterator, class DHFaceIterator>
  inline void
  DoFInfo<dim,spacedim,number>::reinit(
    const DHCellIterator& c,
    const DHFaceIterator& f,
    unsigned int n)
  {
    if ((cell.state() != IteratorState::valid)
	||  cell != static_cast<typename Triangulation<dim>::cell_iterator> (c))
      get_indices(c);
    cell = static_cast<typename Triangulation<dim>::cell_iterator> (c);
    set_face(f,n);
    
    if (block_info)
      LocalResults<number>::reinit(block_info->local());
    else
      LocalResults<number>::reinit(aux_local_indices);
  }


  template<int dim, int spacedim, typename number>
  template <class DHFaceIterator>
  inline void
  DoFInfo<dim,spacedim,number>::set_subface(
    const DHFaceIterator& f,
    unsigned int n,
    unsigned int s)
  {
    face = static_cast<typename Triangulation<dim>::face_iterator> (f);
    face_number = n;
    sub_number = s;
  }
  
  
  template<int dim, int spacedim, typename number>
  template <class DHCellIterator, class DHFaceIterator>
  inline void
  DoFInfo<dim,spacedim,number>::reinit(
    const DHCellIterator& c,
    const DHFaceIterator& f,
    unsigned int n,
    unsigned int s)
  {
    if (cell.state() != IteratorState::valid
	|| cell != static_cast<typename Triangulation<dim>::cell_iterator> (c))
      get_indices(c);
    cell = static_cast<typename Triangulation<dim>::cell_iterator> (c);
    set_subface(f,n,s);
    
    if (block_info)
      LocalResults<number>::reinit(block_info->local());
    else
      LocalResults<number>::reinit(aux_local_indices);
  }


  template<int dim, int spacedim, typename number>
  inline const BlockIndices&
  DoFInfo<dim,spacedim,number>::local_indices() const
  {
    if (block_info)
      return block_info->local();
    return aux_local_indices;
  }

//----------------------------------------------------------------------//

  template <int dim, class DOFINFO>
  inline
  DoFInfoBox<dim, DOFINFO>::DoFInfoBox(const DOFINFO& seed)
		  :
		  cell(seed)
  {
    for (unsigned int i=0;i<GeometryInfo<dim>::faces_per_cell;++i)
      {
	exterior[i] = seed;
	interior[i] = seed;
	interior_face_available[i] = false;
	exterior_face_available[i] = false;
      }
  }


  template <int dim, class DOFINFO>
  inline
  DoFInfoBox<dim, DOFINFO>::DoFInfoBox(const DoFInfoBox<dim, DOFINFO>& other)
		  :
		  cell(other.cell)
  {
    for (unsigned int i=0;i<GeometryInfo<dim>::faces_per_cell;++i)
      {
	exterior[i] = other.exterior[i];
	interior[i] = other.interior[i];
	interior_face_available[i] = false;
	exterior_face_available[i] = false;
      }
  }


  template <int dim, class DOFINFO>
  inline void
  DoFInfoBox<dim, DOFINFO>::reset ()
  {
    for (unsigned int i=0;i<GeometryInfo<dim>::faces_per_cell;++i)
      {
    	interior_face_available[i] = false;
	exterior_face_available[i] = false;
      }
  }


  template <int dim, class DOFINFO>
  template <class ASSEMBLER>
  inline void
  DoFInfoBox<dim, DOFINFO>::assemble (ASSEMBLER& assembler) const
  {
    assembler.assemble(cell);
    for (unsigned int i=0;i<GeometryInfo<dim>::faces_per_cell;++i)
      {
					 // Only do something if data available
	if (interior_face_available[i])
	  {
					     // If both data
					     // available, it is an
					     // interior face
	    if (exterior_face_available[i])
	      assembler.assemble(interior[i], exterior[i]);
	    else
	      assembler.assemble(interior[i]);
	  }
      }
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
