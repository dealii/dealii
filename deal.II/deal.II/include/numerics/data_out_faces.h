//----------------------------  data_out_faces.h  ---------------------------
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
//----------------------------  data_out_faces.h  ---------------------------
#ifndef __deal2__data_out_faces_h
#define __deal2__data_out_faces_h


#include <numerics/data_out.h>

#include <string>
#include <vector>

template <int dim> class DoFHandler;


/**
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
 *
 * @author Wolfgang Bangerth, 2000
 */
template <int dim>
class DataOutFaces : public DataOut_DoFData<dim,dim-1,dim>
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
				      * The function supports
				      * multithreading, if deal.II is
				      * compiled in multithreading
				      * mode. The default number of
				      * threads to be used to build
				      * the patches is set to
				      * @p{multithread_info.n_default_threads}.
				      */
    virtual void
    build_patches (const unsigned int n_subdivisions = 1,
		   const unsigned int n_threads      = multithread_info.n_default_threads);

				     /**
				      * Declare a way to describe a
				      * face which we would like to
				      * generate output for. The usual
				      * way would, of course, be to
				      * use an object of type
				      * @ref{DoFHandler}@p{<dim>::fec_iterator},
				      * but since we have to describe
				      * faces to objects of type
				      * @ref{FEValues}, we can only
				      * represent faces by pairs of a
				      * cell and the number of the
				      * face. This pair is here
				      * aliased to a name that is
				      * better to type.
				      */
    typedef pair<typename DoFHandler<dim>::cell_iterator,unsigned int> FaceDescriptor;
    
    
				     /**
				      * Return the first face which we
				      * want output for. The default
				      * implementation returns the
				      * first active face on the
				      * boundary, but you might want
				      * to return another face in a
				      * derived class.
				      */
    virtual FaceDescriptor first_face ();
    
				     /**
				      * Return the next face after
				      * @p{face} which we want output
				      * for.  If there are no more
				      * face, @p{dofs->end()} shall be
				      * returned as the first
				      * component of the return value.
				      *
				      * The default implementation
				      * returns the next active face
				      * on the boundary, but you might
				      * want to return other faces in
				      * a derived class. Note that the
				      * default implementation assumes
				      * that the given @p{face} is
				      * active, which is guaranteed as
				      * long as @p{first_face} is also
				      * used from the default
				      * implementation. Overloading
				      * only one of the two functions
				      * might not be a good idea.
				      */
    virtual FaceDescriptor next_face (const FaceDescriptor &face);

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
