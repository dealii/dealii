/*----------------------------   persistent_tria.h     ---------------------------*/
/*      $Id$                 */
#ifndef __persistent_tria_H
#define __persistent_tria_H
/*----------------------------   persistent_tria.h     ---------------------------*/


#include <base/smartpointer.h>
#include <grid/tria.h>
#include <vector>


/**
 * This class handles the history of a triangulation and can rebuild it after
 * it was deleted some time before. Its main purpose is support for
 * time-dependent problems where one frequently deletes a triangulation due
 * to memory pressure and later wants to rebuild it; this class has all the
 * information to rebuild it exactly as it was beforem including the mapping
 * of cell numbers to the geometrical cells.
 *
 * Basically, this is a drop-in replacement for the triangulation. Since it
 * is derived from the #Triangulation<dim># class, it shares all the
 * functionality, but it overrides some virtual functions and adds some
 * functions, too. The main change to the base class is that it overrides
 * the #execute_coarsening_and_refinement# function, where the new version
 * first stores all refinement and coarsening flags and only then calls the
 * respective function of the base class. The stored flags may later be
 * used to restore the grid just as it was before. Some other functions
 * have been extended slightly as well, see their documentation for more
 * information.
 *
 * We note that since the triangulation is created in exactly the same state
 * as it was before, other objects working on it should result in the same
 * state as well. This holds in particular for the #DoFHandler# object, which
 * will assign the same degrees of freedom to the original cells and the ones
 * after reconstruction of the triangulation. You can therefore safely use data
 * vectors computed on the original grid on the reconstructed grid as well.
 *
 *
 * \subsection{Usage}
 * You can use objects of this class almost in the same way as objects of the
 * #Triangulation# class. One of the few differences is that you can only
 * construct such an object by giving a coarse grid to the constructor. The
 * coarse grid will be used to base the triangulation on, and therefore the
 * lifetime of the coarse grid has to be longer than the lifetime of the
 * object of this class.
 *
 * Basically, usage looks like this:
 * \begin{verbatim}
 *   Triangulation<dim> coarse_grid;
 *   ...                     // initialize coarse grid
 *
 *   PersistentTriangulation<dim> grid (coarse_grid);
 *
 *   for (...) 
 *     {
 *                           // restore grid from coarse grid
 *                           // and stored refinement flags
 *       grid.restore ();
 *       ...                 // do something with the grid
 *
 *       ...                 // flag some cells for refinement
 *                           // or coarsening
 *       grid.execute_coarsening_and_refinement ();
 *                           // actually refine grid and store
 *                           // the flags
 *
 *       ...                 // so something more with the grid
 *
 *       grid.clear ();      // delete the grid, but keep the
 *                           // refinement flags for later use
 *                           // in grid.restore() above
 *
 *       ...                 // do something where the grid
 *                           // is not needed anymore, e.g.
 *                           // working with another grid
 *     };
 * \end{verbatim}
 *
 * Note that initially, the #PersistentTriangulation# object does not
 * constitute a triangulation; it only becomes one after #restore# is first
 * called. Note also that the #execute_coarsening_and_refinement# stores
 * all necessary flags for later reconstruction using the #restore# function.
 * #Triangulation<dim>::clear ()# resets the underlying triangulation to a
 * virgin state, but does not affect the stored refinement flags needed for
 * later reconstruction and does also not touch the coarse grid which is
 * used withing #restore()#.
 *
 * @author Wolfgang Bangerth, 1999
 */
template <int dim>
class PersistentTriangulation : public Triangulation<dim> 
{
  public:
				     /**
				      * Build up the triangulation from the
				      * coarse grid in future. Copy smoothing
				      * flags, etc from that grid as well.
				      * Note that the initial state of the
				      * triangulation is empty, unless
				      * #restore_grid# is called for the
				      * first time.
				      *
				      * The coarse grid must persist until
				      * the end of this object, since it will
				      * be used upon reconstruction of the
				      * grid.
				      */
    PersistentTriangulation (const Triangulation<dim> &coarse_grid);

				     /**
				      * Copy constructor. This operation
				      * is only allowed, if the triangulation
				      * underlying the object to be copied
				      * is presently empty. Refinement flags
				      * as well as the pointer to the
				      * coarse grid are copied, however.
				      */
    PersistentTriangulation (const PersistentTriangulation<dim> &old_tria);
    
				     /**
				      * Destructor.
				      */
    virtual ~PersistentTriangulation ();
    
				     /**
				      * Overloaded version of the same
				      * function in the base class which
				      * stores the refinement and coarsening
				      * flags for later reconstruction of the
				      * triangulation and after that calls
				      * the respective function of the
				      * base class.
				      */
    virtual void execute_coarsening_and_refinement ();

    				     /**
				      * Restore the grid according to the saved
				      * data. For this, the coarse grid is
				      * copied and the grid is stepwise
				      * rebuilt using the saved flags.
				      *
				      * Note that this function will result in
				      * an error if the underlying triangulation
				      * is not empty, i.e. it will only succeed
				      * if this object is newly created or
				      * #clear()# was called on it before.
				      */
    void restore ();

				     /**
				      * Overload this function to use
				      * #tria# as a new coarse grid. The
				      * present triangulation and all
				      * refinement and coarsening flags
				      * storing its history are deleted,
				      * and the state of the underlying
				      * triangulation is reset to be
				      * empty, until #restore_grid# is
				      * called the next time.
				      *
				      * The coarse grid must persist until
				      * the end of this object, since it will
				      * be used upon reconstruction of the
				      * grid.
				      */
    virtual void copy_triangulation (const Triangulation<dim> &tria);

				     /**
				      * Throw an error, since this function
				      * is not useful in the context of this
				      * class.
				      */
    virtual void block_write (ostream &out) const;

				     /**
				      * Throw an error, since this function
				      * is not useful in the context of this
				      * class.
				      */
    virtual void block_read (istream &in);

				     /**
				      * Throw an error, since this function
				      * is not useful in the context of this
				      * class.
				      */
    virtual void create_triangulation (const vector<Point<dim> >    &vertices,
				       const vector<CellData<dim> > &cells,
				       const SubCellData            &subcelldata);

				     /**
				      * Exception.
				      */
    DeclException0 (ExcFunctionNotUseful);
				     /**
				      * Exception.
				      */
    DeclException0 (ExcTriaNotEmpty);
    
  private:
				     /**
				      * This grid shall be used as coarse
				      * grid.
				      */
    SmartPointer<const Triangulation<dim> > coarse_grid;
    
    				     /**
				      * Vectors holding the refinement and
				      * coarsening flags of the different
				      * sweeps on this time level. The
				      * vectors therefore hold the history
				      * of the grid.
				      */
    vector<vector<bool> >   refine_flags;

				     /**
				      * @see refine_flags
				      */
    vector<vector<bool> >   coarsen_flags;
};



/*----------------------------   persistent_tria.h     ---------------------------*/
/* end of #ifndef __persistent_tria_H */
#endif
/*----------------------------   persistent_tria.h     ---------------------------*/
