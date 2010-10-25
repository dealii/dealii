//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2008, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__dof_levels_h
#define __deal2__dof_levels_h


#include <base/config.h>
#include <base/exceptions.h>
#include <dofs/dof_objects.h>
#include <vector>

DEAL_II_NAMESPACE_OPEN

template <int, int> class DoFHandler;
template <int, int> class MGDoFHandler;


namespace internal
{  
  namespace DoFHandler
  {
    

/**    
 *
 * <h4>DoFLevel@<0@></h4>
 *
 * This class is the common base class for all the DoFLevel
 * classes. We here store information that is associated with
 * (logical) cells, rather than concrete objects such as lines, quads,
 * or hexes.
 *
 * At present, all we store are cached values for the DoF indices on
 * each cell, since this is a frequently requested operation. The
 * values are set by DoFCellAccessor::update_cell_dof_indices_cache
 * and are used by DoFCellAccessor::get_dof_indices.
 *
 * Note that vertices are separate from, and in fact have nothing to
 * do with cells. The indices of degrees of freedom located on
 * vertices therefore are not stored here, but rather in member
 * variables of the dealii::DoFHandler class.
 *
 * The indices of degrees of freedom located on lower dimensional
 * objects, i.e. on lines for 2D and on quads and lines for 3D are
 * treated similarly than that on cells. However, theses geometrical
 * objects, which are called faces as a generalisation, are not
 * organised in a hierachical structure of levels. Therefore, the
 * degrees of freedom located on these objects are stored in seperate
 * classes, namely the <tt>DoFFaces</tt> classes.
 *
 *
 * <h4>DoFLevel@<1@>, DoFLevel@<2@>, and DoFLevel@<3@> </h4>
 *
 * These classes are used, respectively, to store the indices located
 * on lines, quads, and hexes. The storage format is as laid out
 * above, and data is stored in objects lines, quads, and
 * hexes, which are all instantiations of <tt>DoFObjects</tt>.
 * However, it isn't usually directly accessed. Rather,
 * except for some access from the DoFHandler class, access is usually
 * through the DoFAccessor::set_dof_index() and
 * DoFAccessor::dof_index() functions or similar functions of derived
 * classes that in turn access the member variables using the
 * DoFHandler::get_dof_index() and corresponding setter functions.
 * Knowledge of the actual data format is therefore
 * encapsulated to the present hierarchy of classes as well as the
 * dealii::DoFHandler class.
 *
 * @author Wolfgang Bangerth, 1998, 2006
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
 */
    template <>
    class DoFLevel<0>
    {
      public:
					 /**
					  * Cache for the DoF indices
					  * on cells. The size of this
					  * array equals the number of
					  * cells on a given level
					  * times
					  * selected_fe.dofs_per_cell.
					  */
	std::vector<unsigned int> cell_dof_indices_cache;

        
                                         /**
                                          * Determine an estimate for the
                                          * memory consumption (in bytes)
                                          * of this object.
                                          */
        unsigned int memory_consumption () const;
    };
    
/**
 * Store the indices of the degrees of freedom which are located on
 * lines. See the general template DoFLevel for more information.
 *
 * @author Wolfgang Bangerth, 1998, 2006
 */
    template <>
    class DoFLevel<1> : public DoFLevel<0>
    {
      public:
                                         /**
                                          * The object containing dof-indices
					  * and related access-functions
					  */
	DoFObjects<1> lines;

					 /**
                                          * Determine an estimate for the
                                          * memory consumption (in bytes)
                                          * of this object.
                                          */
        unsigned int memory_consumption () const;
    };


/**
 * Store the indices of the degrees of freedom which are located on
 * quads. See the general template DoFLevel for more information.
 *
 * @author Wolfgang Bangerth, 1998, 2006
 */
    template <>
    class DoFLevel<2> : public DoFLevel<0>
    {
      public:
                                         /**
                                          * The object containing dof-indices
					  * and related access-functions
					  */
	internal::DoFHandler::DoFObjects<2> quads;

					 /**
                                          * Determine an estimate for the
                                          * memory consumption (in bytes)
                                          * of this object.
                                          */
        unsigned int memory_consumption () const;
    };


/**
 * Store the indices of the degrees of freedom which are located on
 * hexhedra. See the general template DoFLevel for more information.
 *
 * @author Wolfgang Bangerth, 1998, 2006
 */
    template <>
    class DoFLevel<3> : public DoFLevel<0>
    {
      public:
                                         /**
                                          * The object containing dof-indices
					  * and related access-functions
					  */
	internal::DoFHandler::DoFObjects<3> hexes;
	
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
