/*----------------------------   persistent_tria.h     ---------------------------*/
/*      $Id$                 */
#ifndef __persistent_tria_H
#define __persistent_tria_H
/*----------------------------   persistent_tria.h     ---------------------------*/


#include <base/smartpointer.h>
#include <grid/tria.h>
#include <vector>


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
				      */
    void restore_grid ();

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
