/*----------------------------   intergrid_map.h     ---------------------------*/
/*      $Id$                 */
#ifndef __intergrid_map_H
#define __intergrid_map_H
/*----------------------------   intergrid_map.h     ---------------------------*/



/**
 * This class provides a map between two grids which are derived from
 * the same coarse grid. For each cell iterator of the source map, it provides
 * the respective cell iterator on the destination map, through its
 * #operator []#.
 *
 * Usually, the two grids will be refined differently. Then, the value
 * returned for an iterator on the source grid will be either:
 * \begin{itemize}
 *   \item The same cell on the destination grid, if it exists there;
 *   \item The most refined cell of the destination grid from which the
 *      pendant of the source cell could be obtained by refinement. This
 *      cell is always active and has a refinement level less than that
 *      of the source cell.
 * \end{itemize}
 * Keys for this map are all cells on the source grid, whether active or
 * not.
 *
 * For example, consider these two one-dimensional grids:
 * \begin{verbatim}
 * Grid 1:
 *   x--x--x-----x-----------x
 *    1  2    3        4 
 *
 * Grid 2:
 *   x-----x-----x-----x-----x
 *      1     2     3     4
 * \end{verbatim}
 * (Cell numbers are only given as an example and will not correspond
 * to real cell iterator's indices.) The mapping from grid 1 to grid 2
 * will then be as follows:
 * \begin{verbatim}
 *    Cell on grid 1         Cell on grid 2
 *          1  ------------------>  1
 *          2  ------------------>  1
 *          3  ------------------>  2
 *          4  ------------------>  mother cell of cells 3 and 4
 *                                  (a non-active cell, not shown here)
 * \end{verbatim}
 * Besides the mappings shown here, the non-active cells on grid 1 are also
 * valid keys. For example, the mapping for the mother cell of cells 1 and 2
 * on the first grid will point to cell 1 on the second grid.
 *
 * The implementation of this class is such that not only cell iterators
 * into triangulations can be mapped, but also iterators into objects of
 * type #DoFHandler# and #MGDoFHandler#. The extension to other classes
 * offering iterator functions and some minor additional requirements is
 * simple.
 *
 * Note that this class could in principle be based on the C++ #map<Key,Value>#
 * data type. Instead, it uses another data format which is more effective both
 * in terms of computing time for access as well as with regard to memory
 * consumpion.
 *
 *
 * \section{Usage}
 *
 * In practice, use of this class is as follows:
 * \begin{verbatim}
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
 *                    // corresponding to #cell# (which is one of
 *                    // dof_handler_1
 *     f( grid_1_to_2_map[cell]);
 * \end{verbatim}
 *
 * Note that the template parameters to this class have to be given as
 * #InterGridMap<DoFHandler,2>#, i.e. the dimension is given explicitely and
 * no dimension is attributed to the first parameter, which here is
 * #DoFHandler# (and could equally well be #Triangulation# or #MGDoFHandler#).
 *
 * @author Wolfgang Bangerth, 1999
 */
template <template <int> class GridClass, int dim>
class InterGridMap 
{
  public:

#if (__GNUC__==2) && (__GNUC_MINOR__==95)
				     // helper class
    struct GridClass_dim : public GridClass<dim> {
					 // constructor. will
					 // not be implemented,
					 // but suppresses compiler
					 // warning about non-default
					 // constructor of GridClass
	GridClass_dim ();
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
    vector<vector<cell_iterator> > mapping;

				     /**
				      * Set the mapping for the pair of
				      * cells given. These shall match
				      * in level of refinement and all
				      * other properties.
				      */
    void set_mapping (const cell_iterator &src_cell,
		      const cell_iterator &dst_cell);

				     /**
				      * Set the value of the key #src_cell#
				      * to #dst_cell#. Do so as well for
				      * all the children and their children
				      * of #src_cell#. This function is
				      * used for cells which are more
				      * refined on #src_grid# than on
				      * #dst_grid#; then all values of
				      * the hierarchy of cells and their
				      * children point to one cell on the
				      * #dst_grid#.
				      */
    void set_entries_to_cell (const cell_iterator &src_cell,
			      const cell_iterator &dst_cell);
};



/*----------------------------   intergrid_map.h     ---------------------------*/
/* end of #ifndef __intergrid_map_H */
#endif
/*----------------------------   intergrid_map.h     ---------------------------*/
