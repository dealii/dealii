//----------------------------  data_out_faces.h  ---------------------------
//    Version: $Name$
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
 * This class generates output from faces of a triangulation rather
 * than from cells, as do for example the @ref{DataOut} and
 * @ref{DataOut_Rotation} classes. It might be used to generate output
 * only for the surface of the triangulation (this is the default of
 * this class), or for another arbitrary set of faces. The output of
 * this class is a set of patches (as defined by the class
 * @ref{DataOutBase::Patch}), one for each face for which output is to
 * be generated. These patches can then be written in several
 * graphical data formats by the functions of the underlying classes.
 *
 * @sect3{Interface}
 *
 * The interface of this class is copied from the @ref{DataOut}
 * class. Furthermore, they share the common parent class
 * @ref{DataOut_DoFData}. See the reference of these two classes for a
 * discussion of the interface.
 *
 *
 * @sect3{Extending this class}
 *
 * The sequence of faces to generate patches from is generated in the
 * same way as in the @ref{DataOut} class, see there for a description
 * of the respective interface. For obvious reasons, the functions
 * generating the sequence of faces which shall be used to generate
 * output, are called @p{first_face} and @p{next_face} in this class,
 * rather than @p{first_cell} and @p{next_cell}.
 *
 * Since we need to initialize objects of type @ref{FEValues} with the
 * faces generated from these functions, it is not sufficient that
 * they only return face iterators. Rather, we need a pair of cell and
 * the number of the face, as the values of finite element fields need
 * not necessarily be unique on a face (think of discontinuous finite
 * elements, where the value of the finite element field depend on the
 * direction from which you approach a face, thus it is necessary to
 * use a pair of cell and face, rather than only a face
 * iterator). Therefore, this class defines a @p{typedef} which
 * creates a type @p{FaceDescriptor} that is an abbreviation for a
 * pair of cell iterator and face number. The functions @p{first_face}
 * and @p{next_face} operate on objects of this type.
 *
 * Extending this class might, for example, be useful if you only want
 * output from certain portions of the boundary, e.g. as indicated by
 * the boundary indicator of the respective faces. However, it is also
 * conceivable that one generates patches not from boundary faces, but
 * from interior faces that are selected due to other criteria; one
 * application might be to use only those faces where one component of
 * the solution attains a certain value, in order to display the
 * values of other solution components on these faces. Other
 * applications certainly exist, for which the author is not
 * imaginative enough.
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
    typedef typename std::pair<typename DoFHandler<dim>::cell_iterator,unsigned int> FaceDescriptor;
    
    
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


#endif
