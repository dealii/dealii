//----------------------------  data_out_rotation.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
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
 * This class generates output in the full domain of computations that
 * were done using rotational symmetry of domain and solution. In
 * particular, if a computation of a three dimensional problem with
 * rotational symmetry around the @p{z}-axis (i.e. in the
 * @p{r-z}-plane) was done, then this class can be used to generate
 * the output in the original @p{x-y-z} space. In order to do so, it
 * generates from each cell in the computational mesh a cell in the
 * space with dimension one greater than that of the DoFHandler
 * object. The resulting output will then consist of hexahedra forming
 * an object that has rotational symmetry around the z-axis. As most
 * graphical programs can not represent ring-like structures, the
 * angular (rotation) variable is discretized into a finite number of
 * intervals as well; the number of these intervals must be given to
 * the @p{build_patches} function. It is noted, however, that while
 * this function generates nice pictures of the whole domain, it often
 * produces @em{very} large output files.
 *
 *
 * @sect3{Interface}
 *
 * The interface of this class is copied from the @ref{DataOut}
 * class. Furthermore, they share the common parent class
 * @ref{DataOut_DoFData}. See the reference of these two classes for a
 * discussion of the interface and how to extend it by deriving
 * further classes from this class.
 *
 *
 * @sect3{Details for 1d computations}
 *
 * The one coordinate in the triangulation used by the
 * @ref{DoFHandler} object passed to this class is taken as the radial
 * variable, and the output will then be either a circle or a ring
 * domain. It is in the user's responsibility to assure that the
 * radial coordinate only attains non-negative values.
 *
 *
 * @sect3{Details for 2d computations}
 *
 * We consider the computation (represented by the @ref{DoFHandler}
 * object that is attached to this class) to have happened in the
 * @p{r-z}-plane, where @p{r} is the radial variable and @p{z} denotes
 * the axis of revolution around which the solution is symmetric. The
 * output is in @p{x-y-z} space, where the radial dependence is
 * transformed to the @p{x-y} plane. At present, it is not possible to
 * exchange the meaning of the first and second variable of the plane
 * in which the simulation was made, i.e. generate output from a
 * simulation where the first variable denoted the symmetry axis, and
 * the second denoted the radial variable. You have to take that into
 * account when first programming your application.
 *
 * It is in the responsibility of the user to make sure that the
 * radial variable attains only non-negative values.
 *
 * @author Wolfgang Bangerth, 2000
 */
template <int dim>
class DataOutRotation : public DataOut_DoFData<dim,dim+1>
{
  public:
    				     /**
				      * This is the central function
				      * of this class since it builds
				      * the list of patches to be
				      * written by the low-level
				      * functions of the base
				      * class. See the general
				      * documentation of this class
				      * for further information.
				      *
				      * In addition to the same
				      * parameters as found in the
				      * respective function of the
				      * @ref{DataOut} class, the first
				      * parameter denotes into how
				      * many intervals the angular
				      * (rotation) variable is to be
				      * subdivided.
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
				      * implementation returns the
				      * first active cell, but you
				      * might want to return other
				      * cells in a derived class.
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
				     /**
				      * Exception
				      */
    DeclException1 (ExcRadialVariableHasNegativeValues,
		    double,
		    << "The radial variable attains a negative value of " << arg1);
    
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
	std::vector<double>          patch_values;
	std::vector<Vector<double> > patch_values_system;
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


/* -------------- declaration of explicit specializations ------------- */

template <> void DataOutRotation<3>::build_some_patches (Data);


#endif
