//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__data_out_faces_h
#define __deal2__data_out_faces_h


#include <deal.II/base/config.h>
#include <deal.II/numerics/data_out.h>

#include <string>
#include <vector>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace DataOutFaces
  {
				     /**
				      * A derived class for use in the
				      * DataOutFaces class. This is
				      * a class for the
				      * AdditionalData kind of data
				      * structure discussed in the
				      * documentation of the
				      * WorkStream context.
				      */
    template <int dim, int spacedim>
    struct ParallelData : public internal::DataOut::ParallelDataBase<dim,spacedim>
    {
	template <class FE>
	ParallelData (const Quadrature<dim-1> &quadrature,
		      const unsigned int n_components,
		      const unsigned int n_datasets,
		      const unsigned int n_subdivisions,
		      const std::vector<unsigned int> &n_postprocessor_outputs,
		      const Mapping<dim,spacedim> &mapping,
		      const FE &finite_elements,
		      const UpdateFlags update_flags);

	const dealii::hp::QCollection<dim-1> q_collection;
	const dealii::hp::MappingCollection<dim,spacedim> mapping_collection;
	dealii::hp::FEFaceValues<dim> x_fe_values;

	std::vector<Point<dim> > patch_normals;
	std::vector<Point<spacedim> > patch_evaluation_points;
    };
  }
}


/**
 * This class generates output from faces of a triangulation rather
 * than from cells, as do for example the DataOut and
 * DataOut_Rotation() classes. It might be used to generate output
 * only for the surface of the triangulation (this is the default of
 * this class), or for another arbitrary set of faces. The output of
 * this class is a set of patches (as defined by the class
 * DataOutBase::Patch()), one for each face for which output is to
 * be generated. These patches can then be written in several
 * graphical data formats by the functions of the underlying classes.
 *
 * <h3>Interface</h3>
 *
 * The interface of this class is copied from the DataOut
 * class. Furthermore, they share the common parent class
 * DataOut_DoFData. See the reference of these two classes for a
 * discussion of the interface.
 *
 *
 * <h3>Extending this class</h3>
 *
 * The sequence of faces to generate patches from is generated in the
 * same way as in the DataOut class, see there for a description
 * of the respective interface. For obvious reasons, the functions
 * generating the sequence of faces which shall be used to generate
 * output, are called @p first_face and @p next_face in this class,
 * rather than @p first_cell and @p next_cell.
 *
 * Since we need to initialize objects of type FEValues with the
 * faces generated from these functions, it is not sufficient that
 * they only return face iterators. Rather, we need a pair of cell and
 * the number of the face, as the values of finite element fields need
 * not necessarily be unique on a face (think of discontinuous finite
 * elements, where the value of the finite element field depend on the
 * direction from which you approach a face, thus it is necessary to
 * use a pair of cell and face, rather than only a face
 * iterator). Therefore, this class defines a @p typedef which
 * creates a type @p FaceDescriptor that is an abbreviation for a
 * pair of cell iterator and face number. The functions @p first_face
 * and @p next_face operate on objects of this type.
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
 * @ingroup output
 * @author Wolfgang Bangerth, 2000
 */
template <int dim, class DH=DoFHandler<dim> >
class DataOutFaces : public DataOut_DoFData<DH,DH::dimension-1,
					    DH::dimension>
{
  public:
				     /**
				      * Typedef to the iterator type
				      * of the dof handler class under
				      * consideration.
				      */
    typedef typename DataOut_DoFData<DH,DH::dimension-1,
				     DH::dimension>::cell_iterator cell_iterator;

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
				      * <tt>multithread_info.n_default_threads</tt>.
				      */
    virtual void
    build_patches (const unsigned int n_subdivisions = 0);

				     /**
				      * Same as above, except that the
				      * additional first parameter
				      * defines a mapping that is to
				      * be used in the generation of
				      * output. If
				      * <tt>n_subdivisions>1</tt>, the
				      * points interior of subdivided
				      * patches which originate from
				      * cells at the boundary of the
				      * domain can be computed using the
				      * mapping, i.e. a higher order
				      * mapping leads to a
				      * representation of a curved
				      * boundary by using more
				      * subdivisions.
				      *
				      * Even for non-curved cells the
				      * mapping argument can be used
				      * for the Eulerian mappings (see
				      * class MappingQ1Eulerian) where
				      * a mapping is used not only to
				      * determine the position of
				      * points interior to a cell, but
				      * also of the vertices.  It
				      * offers an opportunity to watch
				      * the solution on a deformed
				      * triangulation on which the
				      * computation was actually
				      * carried out, even if the mesh
				      * is internally stored in its
				      * undeformed configuration and
				      * the deformation is only
				      * tracked by an additional
				      * vector that holds the
				      * deformation of each vertex.
				      *
				      * @todo The @p mapping argument should be
				      * replaced by a hp::MappingCollection in
				      * case of a hp::DoFHandler.
				      */
    virtual void build_patches (const Mapping<DH::dimension> &mapping,
				const unsigned int n_subdivisions = 0);

				     /**
				      * Declare a way to describe a
				      * face which we would like to
				      * generate output for. The usual
				      * way would, of course, be to
				      * use an object of type
				      * <tt>DoFHandler<dim>::face_iterator</tt>,
				      * but since we have to describe
				      * faces to objects of type
				      * FEValues, we can only
				      * represent faces by pairs of a
				      * cell and the number of the
				      * face. This pair is here
				      * aliased to a name that is
				      * better to type.
				      */
    typedef typename std::pair<cell_iterator,unsigned int> FaceDescriptor;
    
    
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
				      * @p face which we want output
				      * for.  If there are no more
				      * face, <tt>dofs->end()</tt> shall be
				      * returned as the first
				      * component of the return value.
				      *
				      * The default implementation
				      * returns the next active face
				      * on the boundary, but you might
				      * want to return other faces in
				      * a derived class. Note that the
				      * default implementation assumes
				      * that the given @p face is
				      * active, which is guaranteed as
				      * long as @p first_face is also
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

				     /**
				      * Exception
				      */
    DeclException0 (ExcCellNotActiveForCellData);
    
  private:
				     /**
				      * Build one patch. This function
				      * is called in a WorkStream
				      * context.
				      */
    void build_one_patch (const FaceDescriptor *cell_and_face,
			  internal::DataOutFaces::ParallelData<DH::dimension, DH::dimension> &data,
			  DataOutBase::Patch<DH::dimension-1,DH::space_dimension> &patch);
};


DEAL_II_NAMESPACE_CLOSE

#endif
