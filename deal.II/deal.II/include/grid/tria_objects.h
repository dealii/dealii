//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__tria_objects_h
#define __deal2__tria_objects_h

#include <base/config.h>
#include <base/exceptions.h>
#include <base/geometry_info.h>
#include <grid/tria_line.h>
#include <grid/tria_quad.h>
#include <grid/tria_hex.h>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace Triangulation
  {

/**
 * General template for information belonging to the geometrical objects of a
 * triangulation, i.e. lines, quads, hexahedrons...  Apart from the vector of
 * objects additional information is included, namely vectors indicating the
 * children, the used-status, user-flags, material-ids..
 *
 * Objects of these classes are include in the TriaLevel and TriaFaces
 * classes.
 *
 * @ingroup grid
 * @author Tobias Leicht, 2006
 */
    
    template <typename G>
    class TriaObjects
    {
      public:
					 /**
					  *  Vector of the objects belonging to
					  *  this level. The index of the object
					  *  equals the index in this container. 
					  */
	std::vector<G> cells;
					 /**
					  *  Index of the first child of an object.
					  *  Since when objects are refined, all
					  *  children are created at the same
					  *  time, they are appended to the list
					  *  after each other.
					  *  We therefore only store the index
					  *  of the first child, the others
					  *  follow immediately afterwards.
					  *
					  *  If an object has no children, -1 is
					  *  stored in this list. An object is
					  *  called active if it has no
					  *  children. The function
					  *  TriaAccessor::has_children()
					  *  tests for this.
					  */
	std::vector<int>  children;
	
					 /**
					  *  Vector storing whether an object is
					  *  used in the @p cells vector.
					  *
					  *  Since it is difficult to delete
					  *  elements in a @p vector, when an
					  *  element is not needed any more
					  *  (e.g. after derefinement), it is
					  *  not deleted from the list, but
					  *  rather the according @p used flag
					  *  is set to @p false.
					  */
	std::vector<bool> used;
	
					 /**
					  *  Make available a field for user data,
					  *  one bit per object. This field is usually
					  *  used when an operation runs over all
					  *  cells and needs information whether
					  *  another cell (e.g. a neighbor) has
					  *  already been processed.
					  *
					  *  You can clear all used flags using
					  *  Triangulation::clear_user_flags().
					  */
	std::vector<bool> user_flags;
	
					 /**
					  * Store boundary and material data. For
					  * example, in one dimension, this field
					  * stores the material id of a line, which
					  * is a number between 0 and 254. In more
					  * than one dimension, lines have no
					  * material id, but they may be at the
					  * boundary; then, we store the
					  * boundary indicator in this field,
					  * which denotes to which part of the
					  * boundary this line belongs and which
					  * boundary conditions hold on this
					  * part. The boundary indicator also
					  * is a number between zero and 254;
					  * the id 255 is reserved for lines
					  * in the interior and may be used
					  * to check whether a line is at the
					  * boundary or not, which otherwise
					  * is not possible if you don't know
					  * which cell it belongs to.
					  */
	std::vector<unsigned char> material_id;

					 /**
					  * Pointer which is not used by the
					  * library but may be accessed and set
					  * by the user to handle data local to
					  * a line/quad/etc.
					  */
	std::vector<void*> user_pointers;

                                         /**
                                          *  Assert that enough space is
                                          *  allocated to accomodate
                                          *  <tt>new_objs</tt> new objects.
                                          *  This function does not only call
                                          *  <tt>vector::reserve()</tt>, but
                                          *  does really append the needed
                                          *  elements.
                                          */
        void reserve_space (const unsigned int new_objs);

					 /**
					  *  Clear all the data contained in this object.
					  */
	void clear();
	
                                         /**
                                          *  Check the memory consistency of the
                                          *  different containers. Should only be
                                          *  called with the prepro flag @p DEBUG
                                          *  set. The function should be called from
                                          *  the functions of the higher
                                          *  TriaLevel classes.
                                          */
        void monitor_memory (const unsigned int true_dimension) const;

                                         /**
                                          * Determine an estimate for the
                                          * memory consumption (in bytes)
                                          * of this object.
                                          */
        unsigned int memory_consumption () const;

                                         /**
                                          *  Exception
                                          */
        DeclException3 (ExcMemoryWasted,
                        char*, int, int,
                        << "The container " << arg1 << " contains "
                        << arg2 << " elements, but it`s capacity is "
                        << arg3 << ".");
                                         /**
                                          *  Exception
                                          */
        DeclException2 (ExcMemoryInexact,
                        int, int,
                        << "The containers have sizes " << arg1 << " and "
                        << arg2 << ", which is not as expected.");
	
    };

/**
 * For hexahedrons the data of TriaObjects needs to be extended, as we can obtain faces
 * (quads) in non-standard-orientation, therefore we declare a class TriaObjectsHex, which
 * additionaly contains a bool-vector of the face-orientations.
 * @ingroup grid
 */
    
    class TriaObjectsHex: public TriaObjects<Hexahedron>
    {
      public:

					 /**
					  * For edges, we enforce a
					  * standard convention that
					  * opposite edges should be
					  * parallel. Now, that's
					  * enforcable in most cases,
					  * and we have code that
					  * makes sure that if a mesh
					  * allows this to happen,
					  * that we have this
					  * convention. We also know
					  * that it is always possible
					  * to have opposite faces
					  * have parallel normal
					  * vectors. (For both things,
					  * see the Agelek, Anderson,
					  * Bangerth, Barth paper
					  * mentioned in the
					  * publications list.)
					  *
					  * The problem is that we
					  * originally had another
					  * condition, namely that
					  * faces 0, 2 and 6 have
					  * normals that point into
					  * the cell, while the other
					  * faces have normals that
					  * point outward. It turns
					  * out that this is not
					  * always possible. In
					  * effect, we have to store
					  * whether the normal vector
					  * of each face of each cell
					  * follows this convention or
					  * not. If this is so, then
					  * this variable stores a
					  * @p true value, otherwise
					  * a @p false value.
					  *
					  * In effect, this field has
					  * <tt>6*n_cells</tt> elements,
					  * being the number of cells
					  * times the six faces each
					  * has.
					  */
	std::vector<bool> face_orientations;

                                         /**
                                          *  Assert that enough space is
                                          *  allocated to accomodate
                                          *  <tt>new_objs</tt> new objects.
                                          *  This function does not only call
                                          *  <tt>vector::reserve()</tt>, but
                                          *  does really append the needed
                                          *  elements.
                                          */
        void reserve_space (const unsigned int new_objs);

					 /**
					  *  Clear all the data contained in this object.
					  */
	void clear();
	
                                         /**
                                          *  Check the memory consistency of the
                                          *  different containers. Should only be
                                          *  called with the prepro flag @p DEBUG
                                          *  set. The function should be called from
                                          *  the functions of the higher
                                          *  TriaLevel classes.
                                          */
        void monitor_memory (const unsigned int true_dimension) const;

                                         /**
                                          * Determine an estimate for the
                                          * memory consumption (in bytes)
                                          * of this object.
                                          */
        unsigned int memory_consumption () const;	    
    };
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
