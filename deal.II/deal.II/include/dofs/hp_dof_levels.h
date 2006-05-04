//----------------------------  hp_dof_levels.h  ------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  hp_dof_levels.h  ------------------------
#ifndef __deal2__hp_dof_levels_h
#define __deal2__hp_dof_levels_h


#include <base/config.h>
#include <dofs/dof_levels.h>
#include <fe/fe_collection.h>

#include <vector>


namespace hp
{
  template <int dim>
  class FECollection;
}


namespace internal
{
  namespace hp
  {

/**
 * Store the indices of the degrees of freedom which are located on
 * objects of dimension @p N. We declare this general template
 * class, but do not actually use it. Rather, only specializations of
 * this class are used.
 *
 * The things we store here is very similar to what is stored in the
 * internal::DoFHandler::DoFLevel class hierarchy (see there for more
 * information, in particular on the layout of the class hierarchy,
 * and the use of file names and dummy arguments). There are two main
 * differences, discussed in the following subsections. In addition to
 * the data already stored by the internal::DoFHandler::DoFLevel
 * classes, we also have to store which finite element each cell
 * uses. This is done in the DoFLevel<0> class, where for each cell we
 * have an entry within the active_fe_indices field for each cell.
 *
 * 
 * <h4>Offset computations</h4>
 * 
 * For hp methods, not all cells may use the same finite element, and
 * it is consequently more complicated to determine where the DoF
 * indices for a given line, quad, or hex are stored. As described in
 * the documentation of the internal::DoFHandler::DoFLevel class, we
 * can compute the location of the first line DoF, for example, by
 * calculating the offset as <code>line_index *
 * dof_handler.get_fe().dofs_per_line</code>. This of course doesn't
 * work any more if different lines may have different numbers of
 * degrees of freedom associated with them. Consequently, rather than
 * using this simple multiplication, each of the line_dofs, quad_dofs,
 * and hex_dofs arrays has an associated array line_dof_offsets,
 * quad_dof_offsets, and hex_dof_offsets. The data corresponding to a
 * line then starts at index <code>line_dof_offsets[line_index]</code>
 * within the <code>line_dofs</code> array. 
 *
 *
 * <h4>Multiple data sets per object</h4>
 *
 * If an object corresponds to a cell, the global dof indices of this
 * cell are stored at the location indicated above in sequential
 * order.
 * 
 * However, if two adjacent cells use different finite elements, then
 * the face that the share needs to store DoF indices for both
 * involved finite elements. While faces therefore have to have at
 * most two sets of DoF indices, it is easy to see that vertices for
 * example can have as many sets of DoF indices associated with them
 * as there are adjacent cells, and the same holds for lines in 3d.
 *
 * Consequently, for objects that have a lower dimensionality than
 * cells, we have to store a map from the finite element index to the
 * set of DoF indices associated. Since real sets are typically very
 * inefficient to store, and since most of the time we expect the
 * number of individual keys to be small (frequently, adjacent cells
 * will have the same finite element, and only a single entry will
 * exist in the map), what we do is instead to store a linked list. In
 * this format, the first entry starting at position
 * <code>line_dofs[line_dof_offsets[line_index]]</code> will denote
 * the finite element index of the set of DoF indices following; after
 * this set, we will store the finite element index of the second set
 * followed by the corresponding DoF indices; and so on. Finally, when
 * all finite element indices adjacent to this object have been
 * covered, we write a -1 to indicate the end of the list.
 *
 * Access to this kind of data, as well as the distinction between
 * cells and objects of lower dimensionality are encoded in the
 * accessor functions, DoFLevel<1>::set_line_dof_index(),
 * DoFLevel<1>::get_line_dof_index(),
 * DoFLevel<2>::set_quad_dof_index(),
 * DoFLevel<2>::get_quad_dof_index(), and
 * DoFLevel<3>::set_hex_dof_index(),
 * DoFLevel<3>::get_hex_dof_index(). The are able to traverse this
 * list and pick out or set a DoF index given the finite element index
 * and its location within the set of DoFs corresponding to this
 * finite element.
 * 
 * 
 * @ingroup hp
 * @author Wolfgang Bangerth, 1998, 2006, Oliver Kayser-Herold 2003.
 */
    template <int N>
    class DoFLevel
    {
      private:
                                         /**
                                          * Make the constructor private
                                          * to avoid that someone uses
                                          * this class.
                                          */
        DoFLevel ();
    };


/**
 * Storage for degrees of freedom on cells. See the documentation of
 * the DoFLevel class template for more complete information on the
 * purpose and layout of this class.
 * 
 * @ingroup hp
 * @author Wolfgang Bangerth, 1998, 2006, Oliver Kayser-Herold 2003.
 */
    template <>
    class DoFLevel<0>
    {
      public:
                                         /**
                                          *  Indices specifying the finite
                                          *  element of hp::FECollection to use
                                          *  for the different cells. The
                                          *  meaning what a cell is, is
                                          *  dimension specific, therefore also
                                          *  the length of this vector depends
                                          *  on the dimension: in one dimension,
                                          *  the length of this vector equals
                                          *  the length of the @p lines vector,
                                          *  in two dimensions that of the @p
                                          *  quads vector, etc.
                                          */

        std::vector<unsigned int> active_fe_indices;
    };


/**
 * Store the indices of the degrees of freedom which are located on
 * the lines. See the general template DoFLevel for more information.
 *
 * @ingroup hp
 * @author Wolfgang Bangerth, 1998, 2006, Oliver Kayser-Herold 2003.
 */
    template <>
    class DoFLevel<1> : public DoFLevel<0>
    {
      private:

                                         /**
                                          * Store the start index for
                                          * the degrees of freedom of each
                                          * line in the @p line_dofs array.
                                          */
        std::vector<unsigned int> line_dof_offsets;

                                         /**
                                          * Store the global indices of
                                          * the degrees of freedom. See
                                          * DoFLevel() for detailed
                                          * information.
                                          */
        std::vector<unsigned int> line_dofs;

      public:
        
                                         /**
                                          * Set the global index of
                                          * the @p local_index-th
                                          * degree of freedom located
                                          * on the line with number @p
                                          * line_index to the value
                                          * given by the last
                                          * argument. The @p
                                          * dof_handler argument is
                                          * used to access the finite
                                          * element that is to be used
                                          * to compute the location
                                          * where this data is stored.
                                          *
                                          * The third argument, @p
                                          * fe_index, denotes which of
                                          * the finite elements
                                          * associated with this
                                          * object we shall
                                          * access. Refer to the
                                          * general documentation of
                                          * the internal::hp::DoFLevel
                                          * class template for more
                                          * information.
					  *
					  * For the meaning of the
					  * last argument, see the
					  * general documentation of
					  * the
					  * internal::DoFHandler::DoFLevel
					  * class.
                                          */
        template <int dim>
        void
        set_dof_index (const ::hp::DoFHandler<dim> &dof_handler,
		       const unsigned int           line_index,
		       const unsigned int           fe_index,
		       const unsigned int           local_index,
		       const unsigned int           global_index,
		       internal::StructuralDimension<1> dummy);

                                         /**
                                          * Return the global index of
                                          * the @p local_index-th
                                          * degree of freedom located
                                          * on the line with number @p
                                          * line_index. The @p
                                          * dof_handler argument is
                                          * used to access the finite
                                          * element that is to be used
                                          * to compute the location
                                          * where this data is stored.
                                          *
                                          * The third argument, @p
                                          * fe_index, denotes which of
                                          * the finite elements
                                          * associated with this
                                          * object we shall
                                          * access. Refer to the
                                          * general documentation of
                                          * the internal::hp::DoFLevel
                                          * class template for more
                                          * information.
					  *
					  * For the meaning of the
					  * last argument, see the
					  * general documentation of
					  * the
					  * internal::DoFHandler::DoFLevel
					  * class.
                                          */
        template <int dim>
        unsigned int
        get_dof_index (const ::hp::DoFHandler<dim> &dof_handler,
		       const unsigned int           line_index,
		       const unsigned int           fe_index,
		       const unsigned int           local_index,
		       internal::StructuralDimension<1> dummy) const;

                                         /**
                                          * Return the number of
                                          * finite elements that are
                                          * active on a given
                                          * object. If this is a cell,
                                          * the answer is of course
                                          * one. If it is a face, the
                                          * answer may be one or two,
                                          * depending on whether the
                                          * two adjacent cells use the
                                          * same finite element or
                                          * not. If it is an edge in
                                          * 3d, the possible return
                                          * value may be one or any
                                          * other value larger than
                                          * that.
                                          *
                                          * The last argument is used
                                          * in the same way as for the
                                          * other functions, i.e. it
                                          * disambiguates between the
                                          * different functions in the
                                          * class hierarchy, whereas
                                          * the template argument
                                          * denotes the space
                                          * dimension we operate in.
                                          */
        template <int dim>
        unsigned int
        n_active_fe_indices (const ::hp::DoFHandler<dim> &dof_handler,
                             const unsigned int           line_index,
                             const StructuralDimension<1>) const;

                                         /**
                                          * Check whether a given
                                          * finite element index is
                                          * used on the present
                                          * object or not.
                                          *
                                          * The last argument is used
                                          * in the same way as for the
                                          * other functions, i.e. it
                                          * disambiguates between the
                                          * different functions in the
                                          * class hierarchy, whereas
                                          * the template argument
                                          * denotes the space
                                          * dimension we operate in.
                                          */
        template <int dim>
        bool
        fe_index_is_active (const ::hp::DoFHandler<dim> &dof_handler,
                            const unsigned int           line_index,
                            const unsigned int           fe_index,
                            const StructuralDimension<1>) const;

                                         /**
                                          * Determine an estimate for the
                                          * memory consumption (in bytes)
                                          * of this object.
                                          */
        unsigned int memory_consumption () const;

                                         /**
                                          * Make the DoFHandler class
                                          * a friend, so that it can
                                          * resize arrays as
                                          * necessary.
                                          */
        template <int> friend class ::hp::DoFHandler;
    };


/**
 * Store the indices of the degrees of freedom which are located on
 * quads. See the general template DoFLevel for more information.
 *
 * @ingroup hp
 * @author Wolfgang Bangerth, 1998, 2006, Oliver Kayser-Herold 2003.
 */
    template <>
    class DoFLevel<2> : public DoFLevel<1>
    {
      private:

                                         /**
                                          * Store the start index for
                                          * the degrees of freedom of each
                                          * quad in the @p quad_dofs array.
                                          */
        std::vector<unsigned int> quad_dof_offsets;

                                         /**
                                          * Store the global indices of
                                          * the degrees of freedom. See
                                          * DoFLevel() for detailed
                                          * information.
                                          */
        std::vector<unsigned int> quad_dofs;

      public:
        
                                         /**
                                          * Set the global index of
                                          * the @p local_index-th
                                          * degree of freedom located
                                          * on the quad with number @p
                                          * quad_index to the value
                                          * given by the last
                                          * argument. The @p
                                          * dof_handler argument is
                                          * used to access the finite
                                          * element that is to be used
                                          * to compute the location
                                          * where this data is stored.
                                          *
                                          * The third argument, @p
                                          * fe_index, denotes which of
                                          * the finite elements
                                          * associated with this
                                          * object we shall
                                          * access. Refer to the
                                          * general documentation of
                                          * the internal::hp::DoFLevel
                                          * class template for more
                                          * information.
					  *
					  * For the meaning of the
					  * last argument, see the
					  * general documentation of
					  * the
					  * internal::DoFHandler::DoFLevel
					  * class.
                                          */
        template <int dim>
        void
        set_dof_index (const ::hp::DoFHandler<dim> &dof_handler,
		       const unsigned int           quad_index,
		       const unsigned int           fe_index,
		       const unsigned int           local_index,
		       const unsigned int           global_index,
		       internal::StructuralDimension<2> dummy);

                                         /**
                                          * Return the global index of
                                          * the @p local_index-th
                                          * degree of freedom located
                                          * on the quad with number @p
                                          * quad_index. The @p
                                          * dof_handler argument is
                                          * used to access the finite
                                          * element that is to be used
                                          * to compute the location
                                          * where this data is stored.
                                          *
                                          * The third argument, @p
                                          * fe_index, denotes which of
                                          * the finite elements
                                          * associated with this
                                          * object we shall
                                          * access. Refer to the
                                          * general documentation of
                                          * the internal::hp::DoFLevel
                                          * class template for more
                                          * information.
 					  *
					  * For the meaning of the
					  * last argument, see the
					  * general documentation of
					  * the
					  * internal::DoFHandler::DoFLevel
					  * class.
                                         */
        template <int dim>
        unsigned int
        get_dof_index (const ::hp::DoFHandler<dim> &dof_handler,
		       const unsigned int           quad_index,
		       const unsigned int           fe_index,
		       const unsigned int           local_index,
		       internal::StructuralDimension<2> dummy) const;

                                         /**
                                          * Return the number of
                                          * finite elements that are
                                          * active on a given
                                          * object. If this is a cell,
                                          * the answer is of course
                                          * one. If it is a face, the
                                          * answer may be one or two,
                                          * depending on whether the
                                          * two adjacent cells use the
                                          * same finite element or
                                          * not. If it is an edge in
                                          * 3d, the possible return
                                          * value may be one or any
                                          * other value larger than
                                          * that.
                                          *
                                          * The last argument is used
                                          * in the same way as for the
                                          * other functions, i.e. it
                                          * disambiguates between the
                                          * different functions in the
                                          * class hierarchy, whereas
                                          * the template argument
                                          * denotes the space
                                          * dimension we operate in.
                                          */
        template <int dim>
        unsigned int
        n_active_fe_indices (const ::hp::DoFHandler<dim> &dof_handler,
                             const unsigned int           quad_index,
                             const StructuralDimension<2>) const;

                                         /**
                                          * Check whether a given
                                          * finite element index is
                                          * used on the present
                                          * object or not.
                                          *
                                          * The last argument is used
                                          * in the same way as for the
                                          * other functions, i.e. it
                                          * disambiguates between the
                                          * different functions in the
                                          * class hierarchy, whereas
                                          * the template argument
                                          * denotes the space
                                          * dimension we operate in.
                                          */
        template <int dim>
        bool
        fe_index_is_active (const ::hp::DoFHandler<dim> &dof_handler,
                            const unsigned int           quad_index,
                            const unsigned int           fe_index,
                            const StructuralDimension<2>) const;

                                         /**
					  * Import the respective
					  * functions from the base
					  * class.
					  */
	using DoFLevel<1>::set_dof_index;
	using DoFLevel<1>::get_dof_index;
	using DoFLevel<1>::n_active_fe_indices;
	using DoFLevel<1>::fe_index_is_active;
	
                                         /**
                                          * Determine an estimate for the
                                          * memory consumption (in bytes)
                                          * of this object.
                                          */
        unsigned int memory_consumption () const;

                                         /**
                                          * Make the DoFHandler class
                                          * a friend, so that it can
                                          * resize arrays as
                                          * necessary.
                                          */
        template <int> friend class ::hp::DoFHandler;
    };



/**
 * Store the indices of the degrees of freedom which are located on
 * hexhedra. See the general template DoFLevel for more information.
 *
 * @ingroup hp
 * @author Wolfgang Bangerth, 1998, 2006, Oliver Kayser-Herold 2003.
 */
    template <>
    class DoFLevel<3> : public DoFLevel<2>
    {
      private:

                                         /**
                                          * Store the start index for
                                          * the degrees of freedom of each
                                          * hex in the @p hex_dofs array.
                                          */
        std::vector<unsigned int> hex_dof_offsets;

                                         /**
                                          * Store the global indices of
                                          * the degrees of freedom. See
                                          * DoFLevel() for detailed
                                          * information.
                                          */
        std::vector<unsigned int> hex_dofs;

      public:
        
                                         /**
                                          * Set the global index of
                                          * the @p local_index-th
                                          * degree of freedom located
                                          * on the hex with number @p
                                          * hex_index to the value
                                          * given by the last
                                          * argument. The @p
                                          * dof_handler argument is
                                          * used to access the finite
                                          * element that is to be used
                                          * to compute the location
                                          * where this data is stored.
                                          *
                                          * The third argument, @p
                                          * fe_index, denotes which of
                                          * the finite elements
                                          * associated with this
                                          * object we shall
                                          * access. Refer to the
                                          * general documentation of
                                          * the internal::hp::DoFLevel
                                          * class template for more
                                          * information.
					  *
					  * For the meaning of the
					  * last argument, see the
					  * general documentation of
					  * the
					  * internal::DoFHandler::DoFLevel
					  * class.
                                          */
        template <int dim>
        void
        set_dof_index (const ::hp::DoFHandler<dim> &dof_handler,
		       const unsigned int           hex_index,
		       const unsigned int           fe_index,
		       const unsigned int           local_index,
		       const unsigned int           global_index,
		       internal::StructuralDimension<3> dummy);

                                         /**
                                          * Return the global index of
                                          * the @p local_index-th
                                          * degree of freedom located
                                          * on the hex with number @p
                                          * hex_index. The @p
                                          * dof_handler argument is
                                          * used to access the finite
                                          * element that is to be used
                                          * to compute the location
                                          * where this data is stored.
                                          *
                                          * The third argument, @p
                                          * fe_index, denotes which of
                                          * the finite elements
                                          * associated with this
                                          * object we shall
                                          * access. Refer to the
                                          * general documentation of
                                          * the internal::hp::DoFLevel
                                          * class template for more
                                          * information.
 					  *
					  * For the meaning of the
					  * last argument, see the
					  * general documentation of
					  * the
					  * internal::DoFHandler::DoFLevel
					  * class.
					  */
        template <int dim>
        unsigned int
        get_dof_index (const ::hp::DoFHandler<dim> &dof_handler,
		       const unsigned int           hex_index,
		       const unsigned int           fe_index,
		       const unsigned int           local_index,
		       internal::StructuralDimension<3> dummy) const;

                                         /**
                                          * Return the number of
                                          * finite elements that are
                                          * active on a given
                                          * object. If this is a cell,
                                          * the answer is of course
                                          * one. If it is a face, the
                                          * answer may be one or two,
                                          * depending on whether the
                                          * two adjacent cells use the
                                          * same finite element or
                                          * not. If it is an edge in
                                          * 3d, the possible return
                                          * value may be one or any
                                          * other value larger than
                                          * that.
                                          *
                                          * The last argument is used
                                          * in the same way as for the
                                          * other functions, i.e. it
                                          * disambiguates between the
                                          * different functions in the
                                          * class hierarchy, whereas
                                          * the template argument
                                          * denotes the space
                                          * dimension we operate in.
                                          */
        template <int dim>
        unsigned int
        n_active_fe_indices (const ::hp::DoFHandler<dim> &dof_handler,
                             const unsigned int           hex_index,
                             const StructuralDimension<3>) const;

                                         /**
                                          * Check whether a given
                                          * finite element index is
                                          * used on the present
                                          * object or not.
                                          *
                                          * The last argument is used
                                          * in the same way as for the
                                          * other functions, i.e. it
                                          * disambiguates between the
                                          * different functions in the
                                          * class hierarchy, whereas
                                          * the template argument
                                          * denotes the space
                                          * dimension we operate in.
                                          */
        template <int dim>
        bool
        fe_index_is_active (const ::hp::DoFHandler<dim> &dof_handler,
                            const unsigned int           hex_index,
                            const unsigned int           fe_index,
                            const StructuralDimension<3>) const;

                                         /**
					  * Import the respective
					  * functions from the base
					  * classes.
					  */
	using DoFLevel<2>::set_dof_index;
	using DoFLevel<2>::get_dof_index;
	using DoFLevel<2>::n_active_fe_indices;
	using DoFLevel<2>::fe_index_is_active;

                                         /**
                                          * Determine an estimate for the
                                          * memory consumption (in bytes)
                                          * of this object.
                                          */
        unsigned int memory_consumption () const;

                                         /**
                                          * Make the DoFHandler class
                                          * a friend, so that it can
                                          * resize arrays as
                                          * necessary.
                                          */
        template <int> friend class ::hp::DoFHandler;
    };


// ------------------ inline and template functions


    
    
  } // namespace hp
  
} // namespace internal

#endif
