//----------------------------  persistent_tria.h  ---------------------------
//    Version: $Name$
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  persistent_tria.h  ---------------------------
#ifndef __deal2__persistent_tria_h
#define __deal2__persistent_tria_h


#include <base/smartpointer.h>
#include <grid/tria.h>
#include <vector>


/**
 * This class handles the history of a triangulation and can rebuild it after
 * it was deleted some time before. Its main purpose is support for
 * time-dependent problems where one frequently deletes a triangulation due
 * to memory pressure and later wants to rebuild it; this class has all the
 * information to rebuild it exactly as it was before including the mapping
 * of cell numbers to the geometrical cells.
 *
 * Basically, this is a drop-in replacement for the triangulation. Since it
 * is derived from the @ref{Triangulation} class, it shares all the
 * functionality, but it overrides some virtual functions and adds some
 * functions, too. The main change to the base class is that it overrides
 * the @p{execute_coarsening_and_refinement} function, where the new version
 * first stores all refinement and coarsening flags and only then calls the
 * respective function of the base class. The stored flags may later be
 * used to restore the grid just as it was before. Some other functions
 * have been extended slightly as well, see their documentation for more
 * information.
 *
 * We note that since the triangulation is created in exactly the same state
 * as it was before, other objects working on it should result in the same
 * state as well. This holds in particular for the @ref{DoFHandler} object, which
 * will assign the same degrees of freedom to the original cells and the ones
 * after reconstruction of the triangulation. You can therefore safely use data
 * vectors computed on the original grid on the reconstructed grid as well.
 *
 *
 * @sect3{Usage}
 * You can use objects of this class almost in the same way as objects of the
 * @ref{Triangulation} class. One of the few differences is that you can only
 * construct such an object by giving a coarse grid to the constructor. The
 * coarse grid will be used to base the triangulation on, and therefore the
 * lifetime of the coarse grid has to be longer than the lifetime of the
 * object of this class.
 *
 * Basically, usage looks like this:
 * @begin{verbatim}
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
 * @end{verbatim}
 *
 * Note that initially, the @ref{PersistentTriangulation} object does not
 * constitute a triangulation; it only becomes one after @p{restore} is first
 * called. Note also that the @p{execute_coarsening_and_refinement} stores
 * all necessary flags for later reconstruction using the @p{restore} function.
 * @ref{Triangulation}@p{<dim>::clear ()} resets the underlying triangulation to a
 * virgin state, but does not affect the stored refinement flags needed for
 * later reconstruction and does also not touch the coarse grid which is
 * used withing @p{restore()}.
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
				      * @p{restore_grid} is called for the
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
				      * Restore the grid according to
				      * the saved data. For this, the
				      * coarse grid is copied and the
				      * grid is stepwise rebuilt using
				      * the saved flags.
				      *
				      * Note that this function will
				      * result in an error if the
				      * underlying triangulation is
				      * not empty, i.e. it will only
				      * succeed if this object is
				      * newly created or @p{clear()}
				      * was called on it before.
				      *
				      * Multiply calles the
				      * @p{restore(unsigned int)}
				      * function in a loop over all
				      * refinement steps.
				      */
    void restore ();

				     /**
				      * Differential restore. Performs
				      * the @p{step_no}th local
				      * refinement and coarsening step.
				      * Step 0 stands for the copying
				      * of the coarse grid.
				      *
				      * This function will only
				      * succeed if the triangulation
				      * is in just the state it were
				      * if restore would have been
				      * called from
				      * @p{step=0...step_no-1} before.
				      */
    void restore (const unsigned int step_no);

				     /**
				      * Returns the number of
				      * refinement and coarsening
				      * steps. This is given by the
				      * size of the @p{refine_flags}
				      * vector.
				      */
    unsigned int n_refinement_steps () const;

				     /**
				      * Overload this function to use
				      * @p{tria} as a new coarse grid. The
				      * present triangulation and all
				      * refinement and coarsening flags
				      * storing its history are deleted,
				      * and the state of the underlying
				      * triangulation is reset to be
				      * empty, until @p{restore_grid} is
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
    virtual void create_triangulation (const typename std::vector<Point<dim> >    &vertices,
				       const typename std::vector<CellData<dim> > &cells,
				       const SubCellData                          &subcelldata);

				     /**
				      * Writes all refine and coarsen
				      * flags to the ostream @p{out}.
				      */
    virtual void write_flags(std::ostream &out) const;

				     /**
				      * Reads all refine and coarsen flags
				      * that previously were written by
				      * @p{write_flags(...)}. This is especially
				      * useful for rebuilding the triangulation
				      * after the end or breakdown of a program
				      * and its restart.
				      */
    virtual void read_flags(std::istream &in);

				     /**
				      * Clears all flags. Retains the
				      * same coarse grid.
				      */
    virtual void clear_flags();

				     /**
				      * Determine an estimate for the
				      * memory consumption (in bytes)
				      * of this object.
				      */
    virtual unsigned int memory_consumption () const;

				     /**
				      * Exception.
				      */
    DeclException0 (ExcFunctionNotUseful);
				     /**
				      * Exception.
				      */
    DeclException0 (ExcTriaNotEmpty);
				     /**
				      * Exception.
				      */
    DeclException0 (ExcFlagsNotCleared);
    
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
    std::vector<std::vector<bool> >   refine_flags;

				     /**
				      * @see refine_flags
				      */
    std::vector<std::vector<bool> >   coarsen_flags;
};


#endif
