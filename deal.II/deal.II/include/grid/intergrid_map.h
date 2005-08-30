//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__intergrid_map_h
#define __deal2__intergrid_map_h

#include <base/config.h>
#include <base/smartpointer.h>
#include <dofs/dof_accessor.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>


/**
 * This class provides a map between two grids which are derived from
 * the same coarse grid. For each cell iterator of the source map, it provides
 * the respective cell iterator on the destination map, through its
 * <tt>operator []</tt>.
 *
 * Usually, the two grids will be refined differently. Then, the value
 * returned for an iterator on the source grid will be either:
 * <ul>
 *   <li> The same cell on the destination grid, if it exists there;
 *   <li> The most refined cell of the destination grid from which the
 *      pendant of the source cell could be obtained by refinement. This
 *      cell is always active and has a refinement level less than that
 *      of the source cell.
 * </ul>
 * Keys for this map are all cells on the source grid, whether active or
 * not.
 *
 * For example, consider these two one-dimensional grids:
 * @verbatim
 * Grid 1:
 *   x--x--x-----x-----------x
 *    1  2    3        4 
 *
 * Grid 2:
 *   x-----x-----x-----x-----x
 *      1     2     3     4
 * @endverbatim
 * (Cell numbers are only given as an example and will not correspond
 * to real cell iterator's indices.) The mapping from grid 1 to grid 2
 * will then be as follows:
 * @verbatim
 *    Cell on grid 1         Cell on grid 2
 *          1  ------------------>  1
 *          2  ------------------>  1
 *          3  ------------------>  2
 *          4  ------------------>  mother cell of cells 3 and 4
 *                                  (a non-active cell, not shown here)
 * @endverbatim
 * Besides the mappings shown here, the non-active cells on grid 1 are also
 * valid keys. For example, the mapping for the mother cell of cells 1 and 2
 * on the first grid will point to cell 1 on the second grid.
 *
 * The implementation of this class is such that not only cell iterators
 * into triangulations can be mapped, but also iterators into objects of
 * type DoFHandler and MGDoFHandler. The extension to other classes
 * offering iterator functions and some minor additional requirements is
 * simple.
 *
 * Note that this class could in principle be based on the C++ <tt>map<Key,Value></tt>
 * data type. Instead, it uses another data format which is more effective both
 * in terms of computing time for access as well as with regard to memory
 * consumpion.
 *
 *
 * <h3>Usage</h3>
 *
 * In practice, use of this class is as follows:
 * @verbatim
 *                   // have two grids, which are derived from the
 *                   // same coarse grid
 *   Triangulation<dim> tria1, tria2;
 *   DoFHandler<dim> dof_handler_1(tria1), dof_handler_2(tria2);
 *   ...
 *                   // do something with these objects, e.g.
 *                   // refine the triangulations differently,
 *                   // distribute degrees of freedom, etc
 *   ...
 *                   // create the mapping
 *   InterGridMap<DoFHandler,dim> grid_1_to_2_map;
 *   grid_1_to_2_map.make_mapping (dof_handler_1,
 *                                 dof_handler_2);
 *   ...
 *   typename DoFHandler<dim>::cell_iterator cell = dof_handler_1.begin(),
 *                                           endc = dof_handler_1.end();
 *   for (; cell!=endc; ++cell)
 *                    // now do something with the cell of dof_handler_2
 *                    // corresponding to @p cell (which is one of
 *                    // dof_handler_1
 *     f( grid_1_to_2_map[cell]);
 * @endverbatim
 *
 * Note that the template parameters to this class have to be given as
 * <tt>InterGridMap<DoFHandler,2></tt>, i.e. the dimension is given explicitly and
 * no dimension is attributed to the first parameter, which here is
 * DoFHandler (and could equally well be Triangulation or MGDoFHandler).
 *
 * @author Wolfgang Bangerth, 1999
 */
template <template <int> class GridClass, int dim>
class InterGridMap 
{
  public:

#ifdef DEAL_II_TEMPLATE_TEMPLATE_TYPEDEF_BUG
				     // helper class
    struct GridClass_dim : public GridClass<dim>
    {
					 // constructor. will
					 // not be implemented,
					 // but suppresses compiler
					 // warning about non-default
					 // constructor of GridClass
	GridClass_dim ();

					 /**
					  * Declare iterator type, for
					  * access from outside.
					  */
	typedef typename GridClass<dim>::cell_iterator cell_iterator;
    };
    
				     /**
				      * Typedef to the iterator type of
				      * the grid class under consideration.
				      */
    typedef typename GridClass_dim::cell_iterator cell_iterator;

#else

				     /**
				      * Typedef to the iterator type of
				      * the grid class under consideration.
				      */
    typedef typename GridClass<dim>::cell_iterator cell_iterator;
#endif
				     /**
				      * Constructor setting the class
				      * name arguments in the
				      * SmartPointer members.
				      */
    InterGridMap();
    
				     /**
				      * Create the mapping between the two
				      * grids.
				      */
    void make_mapping (const GridClass<dim> &source_grid,
		       const GridClass<dim> &destination_grid);

				     /**
				      * Access operator: give a cell
				      * on the source grid and receive
				      * the respective cell on the
				      * other grid, or if that does not
				      * exist, the most refined cell
				      * of which the source cell would
				      * be created if it were further
				      * refined.
				      */
    cell_iterator operator [] (const cell_iterator &source_cell) const;

				     /**
				      * Delete all data of this class.
				      */
    void clear ();

				     /**
				      * Return a pointer to the source
				      * grid.
				      */
    const GridClass<dim> & get_source_grid () const;
    
				     /**
				      * Return a pointer to the
				      * destination grid.
				      */
    const GridClass<dim> & get_destination_grid () const;

				     /**
				      * Determine an estimate for the
				      * memory consumption (in bytes)
				      * of this object.
				      */
    unsigned int memory_consumption () const;
    
				     /**
				      * Exception
				      */
    DeclException1 (ExcInvalidKey,
		    cell_iterator,
		    << "The iterator " << arg1 << " is not valid as key for "
		    << "this map.");
				     /**
				      * Exception
				      */
    DeclException0 (ExcIncompatibleGrids);
    
  private:
				     /**
				      * The actual data. Hold one iterator
				      * for each cell on each level.
				      */
    std::vector<std::vector<cell_iterator> > mapping;

				     /**
				      * Store a pointer to the source grid.
				      */
    SmartPointer<const GridClass<dim> > source_grid;

				     /**
				      * Likewise for the destination grid.
				      */
    SmartPointer<const GridClass<dim> > destination_grid;

				     /**
				      * Set the mapping for the pair of
				      * cells given. These shall match
				      * in level of refinement and all
				      * other properties.
				      */
    void set_mapping (const cell_iterator &src_cell,
		      const cell_iterator &dst_cell);

				     /**
				      * Set the value of the key @p src_cell
				      * to @p dst_cell. Do so as well for
				      * all the children and their children
				      * of @p src_cell. This function is
				      * used for cells which are more
				      * refined on @p src_grid than on
				      * @p dst_grid; then all values of
				      * the hierarchy of cells and their
				      * children point to one cell on the
				      * @p dst_grid.
				      */
    void set_entries_to_cell (const cell_iterator &src_cell,
			      const cell_iterator &dst_cell);
};


#endif
