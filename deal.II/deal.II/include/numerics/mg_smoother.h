/*----------------------------   mg_smoother.h     ---------------------------*/
/*      $Id$                 */
#ifndef __mg_smoother_H
#define __mg_smoother_H
/*----------------------------   mg_smoother.h     ---------------------------*/


#include <lac/forward-declarations.h>
#include <basic/forward-declarations.h>
#include <base/subscriptor.h>
#include <vector>



/**
 * Abstract base class for multigrid smoothers. It declares the interface
 * to smoothers and implements some functionality for setting the values
 * of vectors at interior boundaries (i.e. boundaries between differing
 * levels of the triangulation) to zero, by building a list of these degrees
 * of freedom's indices at construction time.
 *
 * @author Wolfgang Bangerth, Guido Kanschat, 1999
 */
class MGSmoother :  public Subscriptor 
{
  private:
				     /**
				      * Default constructor. Made private
				      * to prevent it being called, which
				      * is necessary since this could
				      * circumpass the set-up of the list
				      * if interior boundary dofs.
				      */
    MGSmoother ();
    
  public:

				     /**
				      * Constructor. This one sets clobbers
				      * the indices of the degrees of freedom
				      * on the interior boundaries between
				      * the different levels, which are
				      * needed by the function
				      * #set_zero_interior_boundaries#.
				      */
    template <int dim>
    MGSmoother (const MGDoFHandler<dim> &mg_dof);

				     /**
				      * Destructor, made virtual.
				      */
    virtual ~MGSmoother ();
    
				     /**
				      * Smooth the vector #u# on the given
				      * level. This function should set the
				      * interior boundary values of u to zero
				      * at the end of its work, so you may want
				      * to call #set_zero_interior_boundary#
				      * at the end of your derived function,
				      * or another function doing similar
				      * things.
				      *
				      * This function is responsible for the
				      * presmoothing, postsmoothing is done
				      * by another function, which, however,
				      * calls this function if not overloaded.
				      */
    virtual void pre_smooth (const unsigned int  level,
			     Vector<float>      &u) const = 0;

				     /**
				      * Postsmoothing; the same applies as for
				      * the presmoothing function. If you
				      * want pre- and postsmoothing to do the
				      * same operations, you may want to not
				      * overload this function, since its
				      * default implementation simply calls
				      * #pre_smooth#.
				      */
    virtual void post_smooth (const unsigned int  level,
			      Vector<float>      &u) const;

				     /**
				      * Reset the values of the degrees of
				      * freedom on interior boundaries between
				      * different levels to zero in the given
				      * data vector #u#.
				      */
    void set_zero_interior_boundary (const unsigned int  level,
				     Vector<float>      &u) const;

  private:
				     /**
				      * For each level, we store a list of
				      * degree of freedom indices which are
				      * located on interior boundaries between
				      * differing levels of the triangulation.
				      *
				      * These arrays are set by the constructor.
				      */
    vector<vector<int> > interior_boundary_dofs;
};

    


/*----------------------------   mg_smoother.h     ---------------------------*/
/* end of #ifndef __mg_smoother_H */
#endif
/*----------------------------   mg_smoother.h     ---------------------------*/
			 
