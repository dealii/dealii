//----------------------------  data_out_rotation.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  data_out_rotation.h  ---------------------------
#ifndef __deal2__data_out_rotation_h
#define __deal2__data_out_rotation_h


#include <numerics/data_out.h>

#include <string>
#include <vector>

template <int dim> class DoFHandler;

/**
 *
 * @author Wolfgang Bangerth, 2000
 */
template <int dim>
class DataOutRotation : public DataOut_DoFData<dim,dim+1>
{
  public:
    				     /**
				      * This is the central function of
				      * this class since it builds the list of
				      * patches to be written by the low-level
				      * functions of the base class. See the
				      * general documentation of this class
				      * for further information.
				      *
				      * The function supports
				      * multithreading, if deal.II is
				      * compiled in multithreading
				      * mode. The default number of
				      * threads to be used to build
				      * the patches is set to
				      * @p{multithread_info.n_default_threads}.
				      */
    virtual void build_patches (const unsigned int n_patches_per_circle,
				const unsigned int n_subdivisions = 1,
				const unsigned int n_threads      = multithread_info.n_default_threads);

				     /**
				      * Return the first cell which we
				      * want output for. The default
				      * implementation returns the first
				      * active cell, but you might want
				      * to return other cells in a derived
				      * class.
				      */
    virtual typename DoFHandler<dim>::cell_iterator
    first_cell ();
    
				     /**
				      * Return the next cell after @p{cell} which
				      * we want output for.
				      * If there are no more cells,
				      * @p{dofs->end()} shall be returned.
				      *
				      * The default
				      * implementation returns the next
				      * active cell, but you might want
				      * to return other cells in a derived
				      * class. Note that the default

				      * implementation assumes that
				      * the given @p{cell} is active, which
				      * is guaranteed as long as @p{first_cell}
				      * is also used from the default
				      * implementation. Overloading only one
				      * of the two functions might not be
				      * a good idea.
				      */
    virtual typename DoFHandler<dim>::cell_iterator
    next_cell (const typename DoFHandler<dim>::cell_iterator &cell);

				     /**
				      * Exception
				      */
    DeclException1 (ExcInvalidNumberOfSubdivisions,
		    int,
		    << "The number of subdivisions per patch, " << arg1
		    << ", is not valid.");
    
  private:
				     /**
				      * All data needed in one thread
				      * is gathered in the struct
				      * Data.
				      * The data is handled globally
				      * to avoid allocation of memory
				      * in the threads.
				      */
    struct Data 
    {
	unsigned int n_threads;
	unsigned int this_thread;
	unsigned int n_components;
	unsigned int n_datasets;
	unsigned int n_patches_per_circle;
	unsigned int n_subdivisions;
	vector<double>          patch_values;
	vector<Vector<double> > patch_values_system;
	Data ()
	  {}
    };
				     /**
				      * Builds every @p{n_threads}'s
				      * patch. This function may be
				      * called in parallel.
				      * If multithreading is not
				      * used, the function is called
				      * once and generates all patches.
				      */
    void build_some_patches (Data data);
};


#endif
