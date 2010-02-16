//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2006, 2007, 2008, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef __deal2__mesh_worker_info_h
#define __deal2__mesh_worker_info_h

#include <base/config.h>
#include <lac/matrix_block.h>
#include <boost/shared_ptr.hpp>
#include <dofs/block_info.h>
#include <fe/fe_values.h>
#include <numerics/mesh_worker_vector_selector.h>

DEAL_II_NAMESPACE_OPEN

namespace MeshWorker
{
/**
 * The class providing the scrapbook to fill with local integration
 * results. These can be values, local contributions to forms or cell
 * and face matrices.
 *
 * The local matrices initialized by reinit() of the info object and
 * then assembled into the global system by Assembler classes.
 *
 * @ingroup MeshWorker
 * @author Guido Kanschat, 2009
 */
  template <typename number>
  class LocalResults
  {
    private:
				       /**
					* Initialize a single local
					* matrix block. A helper
					* function for initialize()
					*/
      void initialize_local(MatrixBlock<FullMatrix<number> >& M,
			    const unsigned int row,
			    const unsigned int col);

    public:
      void initialize_numbers(const unsigned int n);
      void initialize_vectors(const unsigned int n);
				       /**
					* Allocate @p n local
					* matrices. Additionally,
					* set their block row and
					* column coordinates to
					* zero. The matrices
					* themselves are resized by
					* reinit().
					*
					* The template parameter @p
					* MatrixPtr should point to
					* a MatrixBlock
					* instantiation in order to
					* provide row and column info.
					*/
      void initialize_matrices(unsigned int n, bool both);

				       /**
					* Allocate a local matrix
					* for each of the global
					* ones in @p
					* matrices. Additionally,
					* set their block row and
					* column coordinates. The
					* matrices themselves are
					* resized by reinit().
					*
					* The template parameter @p
					* MatrixPtr should point to
					* a MatrixBlock
					* instantiation in order to
					* provide row and column info.
					*/
      template <class MatrixPtr>
      void initialize_matrices(const std::vector<MatrixPtr>& matrices,
			       bool both);

				       /**
					* Reinitialize matrices for
					* new cell. Resizes the
					* matrices for hp and sets
					* them to zero.
					*/
      void reinit(const BlockIndices& local_sizes);

				       /**
					* The number of scalar values.
					*/
      unsigned int n_values () const;
      
				       /**
					* The number of vectors.
					*/
      unsigned int n_vectors () const;
      
				       /**
					* The number of matrices.
					*/
      unsigned int n_matrices () const;
      
				       /**
					* Access scalar value at index
					* @p i.
					*/
      number& value(unsigned int i);

				       /**
					* Read scalar value at index
					* @p i.
					*/
      number value(unsigned int i) const;

				       /**
					* Access vector at index @p i.
					*/
      BlockVector<number>& vector(unsigned int i);
      
				       /**
					* Read vector at index @p i.
					*/
      const BlockVector<number>& vector(unsigned int i) const;
      
				       /**
					* Access matrix at index @p
					* i. For results on internal
					* faces, a true value for @p
					* external refers to the flux
					* between cells, while false
					* refers to entries coupling
					* inside the cell.
					*/
      MatrixBlock<FullMatrix<number> >& matrix(unsigned int i, bool external = false);
      
				       /**
					* Read matrix at index @p
					* i. For results on internal
					* faces, a true value for @p
					* external refers to the flux
					* between cells, while false
					* refers to entries coupling
					* inside the cell.
					*/
      const MatrixBlock<FullMatrix<number> >& matrix(unsigned int i, bool external = false) const;
      
    private:
				       /**
					* The local numbers,
					* computed on a cell or on a
					* face.
					*/
      std::vector<number> J;

				       /**
					* The local vectors. This
					* field is public, so that
					* local integrators can
					* write to it.
					*/
      std::vector<BlockVector<number> > R;

				       /**
					* The local matrices
					* coupling degrees of
					* freedom in the cell
					* itself or within the
					* first cell on a face.
					*/
      std::vector<MatrixBlock<FullMatrix<number> > > M1;

				       /**
					* The local matrices
					* coupling test functions on
					* the cell with trial
					* functions on the other
					* cell.
					*
					* Only used on interior
					* faces.
					*/
      std::vector<MatrixBlock<FullMatrix<number> > > M2;
  };

  template <int dim, int spacedim> class DoFInfoBox;

/**
 * A class containing information on geometry and degrees of freedom
 * of a mesh object.
 *
 * The information in these objects is usually used by one of the
 * Assembler classes. It is also the kind of information which is
 * needed in mesh based matrices (often referred to as matrix free
 * methods).
 *
 * In addition to the information on degrees of freedom stored in this
 * class, it also provides the local computation space for the worker
 * object operating on it in LocalResults. This base class will automatically
 * reinitialized on each cell, but initial setup is up to the user and
 * should be done when initialize() for this class is called.
 *
 * This class operates in two different modes, corresponding to the
 * data models discussed in the Assembler namespace documentation.
 *
 * The choice of the local data model is triggered by the vector
 * #BlockInfo::local_renumbering, which in turn is usually filled by
 * BlockInfo::initialize_local(). If this function has been used, or
 * the vector has been changed from zero-length, then local dof
 * indices stored in this object will automatically be renumbered to
 * reflect local block structure. This means, the first entries in
 * #indices will refer to the first block of the system, then comes
 * the second block and so on.
 *
 * The BlockInfo object is stored as a pointer. Therefore, if the
 * block structure changes, for instance because of mesh refinement,
 * the DoFInfo class will automatically use the new structures.
 *
 * @ingroup MeshWorker
 * @author Guido Kanschat, 2009
 */
  template<int dim, int spacedim = dim>
  class DoFInfo : public LocalResults<double>
  {
    public:
				       /// The current cell
      typename Triangulation<dim>::cell_iterator cell;

				       /// The current face
      typename Triangulation<dim>::face_iterator face;

				       /**
					* The number of the current
					* face on the current cell.
					*
					* This number is
					* deal_II_numbers::invalid_unsigned_int
					* if the info object was
					* initialized with a cell.
					*/

      unsigned int face_number;
				       /**
					* The number of the current
					* subface on the current
					* face
					*
					* This number is
					* deal_II_numbers::invalid_unsigned_int
					* if the info object was not
					* initialized with a subface.
					*/
      unsigned int sub_number;

				       /// The DoF indices of the current cell
      std::vector<unsigned int> indices;

				       /**
					* Constructor setting the
					* #block_info pointer.
					*/
      DoFInfo(const BlockInfo& block_info);

				       /**
					* Default constructor
					* leaving the #block_info
					* pointer empty, but setting
					* the #aux_local_indices.
					*/
      template <class DH>
      DoFInfo (const DH& dof_handler);
      
				       /**
					* Set the current cell and
					* fill #indices.
					*/
      template <class DHCellIterator>
      void reinit(const DHCellIterator& c);

				       /**
					* Set the current face and
					* fill #indices if the #cell
					* changed.
					*/
      template <class DHCellIterator, class DHFaceIterator>
      void reinit(const DHCellIterator& c,
		  const DHFaceIterator& f,
		  const unsigned int n);

				       /**
					* Set the current subface
					* and fill #indices if the
					* #cell changed.
					*/
      template <class DHCellIterator, class DHFaceIterator>
      void reinit(const DHCellIterator& c,
		  const DHFaceIterator& f,
		  const unsigned int n,
		  const unsigned int s);

      const BlockIndices& local_indices() const;


				       /// The block structure of the system
      SmartPointer<const BlockInfo,DoFInfo<dim,spacedim> > block_info;

    private:
				       /**
					* Standard constructor, not
					* setting any block
					* indices. Use of this
					* constructor is not
					* recommended, but it is
					* needed for the arrays in
					* DoFInfoBox.
					*/
      DoFInfo ();
      
      				       /// Fill index vector with active indices
      void get_indices(const typename DoFHandler<dim, spacedim>::cell_iterator& c);

				       /// Fill index vector with level indices
      void get_indices(const typename MGDoFHandler<dim, spacedim>::cell_iterator& c);

				       /// Auxiliary vector
      std::vector<unsigned int> indices_org;

				       /**
					* An auxiliary local
					* BlockIndices object created
					* if #block_info is not set.
					* It contains just a single
					* block of the size of
					* degrees of freedom per cell.
					*/
      BlockIndices aux_local_indices;
      
      friend class DoFInfoBox<dim, spacedim>;
  };


  /**
 * A class bundling the MeshWorker::DoFInfo objects used on a cell.
 *
 * @todo Currently, we are storing an object for the cells and two for
 * each face. We could gather all face data pertaining to the cell
 * itself in one object, saving a bit of memory and a few operations,
 * but sacrificing some cleanliness.
 *
 * @author Guido Kanschat, 2010
 */
  template < int dim, int spacedim=dim>
  class DoFInfoBox
  {
    public:
				       /**
					* Constructor copying the seed
					* into all other objects.
					*/
      DoFInfoBox(const MeshWorker::DoFInfo<dim, spacedim>& seed);

				       /**
					* Reset all the availability flags.
					*/
      void reset();
      
				       /**
					* The data for the cell.
					*/
      MeshWorker::DoFInfo<dim> cell;
				       /**
					* The data for the faces from inside.
					*/
      MeshWorker::DoFInfo<dim> interior[GeometryInfo<dim>::faces_per_cell];
				       /**
					* The data for the faces from outside.
					*/
      MeshWorker::DoFInfo<dim> exterior[GeometryInfo<dim>::faces_per_cell];

				       /**
					* A set of flags, indicating
					* whether data on an interior
					* face is available.
					*/
      bool interior_face_available[GeometryInfo<dim>::faces_per_cell];
				       /**
					* A set of flags, indicating
					* whether data on an exterior
					* face is available.
					*/
      bool exterior_face_available[GeometryInfo<dim>::faces_per_cell];
  };


  

/**
 * Class for objects handed to local integration functions.
 *
 * Objects of this class contain one or more objects of type FEValues,
 * FEFaceValues or FESubfacevalues to be used in local
 * integration. They are stored in an array of pointers to the base
 * classes FEValuesBase for cells and FEFaceValuesBase for faces and
 * subfaces, respectively. The template parameter VECTOR allows the
 * use of different data types for the global system.
 *
 * The @p FEVALUESBASE template parameter should be either
 * FEValuesBase or FEFaceValuesBase, depending on whether the object
 * is used to integrate over cells or faces. The actual type of @p
 * FEVALUES object is fixed in the constructor and only used to
 * initialize the pointers in #fevalv.
 *
 * Additionally, this function containes space to store the values of
 * finite element functions stored in #global_data in the
 * quadrature points. These vectors are initialized automatically on
 * each cell or face. In order to avoid initializing unused vectors,
 * you can use initialize_selector() to select the vectors by name
 * that you actually want to use.
 *
 * <h3>Integration models</h3>
 *
 * This class supports two local integration models, corresponding to
 * the data models in the documentation of the Assembler namespace.
 * One is the
 * standard model suggested by the use of FESystem. Namely, there is
 * one FEValuseBase object in this class, containing all shape
 * functions of the whole system, and having as many components as the
 * system. Using this model involves loops over all system shape
 * functions. It requires to identify the system components
 * for each shape function and to select the correct bilinear form,
 * usually in an @p if or @p switch statement.
 *
 * The second integration model builds one FEValuesBase object per
 * base element of the system. The degrees of freedom on each cell are
 * renumbered by block, such that they represent the same block
 * structure as the global system. Objects performing the integration
 * can then process each block separately, which improves reusability
 * of code considerably.
 *
 * @note As described in DoFInfo, the use of the local block model is
 * triggered by calling BlockInfo::initialize_local() before
 * using initialize() in this class.
 *
 * @ingroup MeshWorker
 * @author Guido Kanschat, 2009
 */
  template<int dim, class FEVALUESBASE, int spacedim = dim>
  class IntegrationInfo
  {
    private:
				       /// vector of FEValues objects
      std::vector<boost::shared_ptr<FEVALUESBASE> > fevalv;
    public:
				       /**
					* Constructor.
					*/
      IntegrationInfo();

				       /**
					* Build all internal
					* structures, in particular
					* the FEValuesBase objects
					* and allocate space for
					* data vectors.
					*
					* @param el is the finite
					* element of the DoFHandler.
					*
					* @param mapping is the Mapping
					* object used to map the
					* mesh cells.
					*
					* @param quadrature is a
					* Quadrature formula used in
					* the constructor of the
					* FEVALUES objects.
					*
					* @param flags are the
					* UpdateFlags used in
					* the constructor of the
					* FEVALUES objects.
					*
					* @param local_block_info is
					* an optional parameter for
					* systems of PDE. If it is
					* provided with reasonable
					* data, then the degrees of
					* freedom on the cells will be
					* re-ordered to reflect the
					* block structure of the system.
					*/
      template <class FEVALUES>
      void initialize(const FiniteElement<dim,spacedim>& el,
		      const Mapping<dim,spacedim>& mapping,
		      const Quadrature<FEVALUES::integral_dimension>& quadrature,
		      const UpdateFlags flags,
		      const BlockInfo* local_block_info = 0);

				       /**
					* Initialize the data
					* vector and cache the
					* selector.
					*/
      void initialize_data(const boost::shared_ptr<VectorDataBase<dim,spacedim> > data);

				       /**
					* Delete the data created by initialize().
					*/
      void clear();

				       /// This is true if we are assembling for multigrid
      bool multigrid;
				       /// Access to finite element
				       /**
					* This is the access
					* function being used, if
					* the constructor for a
					* single element was
					* used. It throws an
					* exception, if applied to a
					* vector of elements.
					*/
      const FEVALUESBASE& fe_values () const;

				       /// Access to finite elements
				       /**
					* This access function must
					* be used if the constructor
					* for a group of elements
					* was used.
					*
					* @see DGBlockSplitApplication
					*/
      const FEVALUESBASE& fe_values (const unsigned int i) const;

				       /**
					* The vector containing the
					* values of finite element
					* functions in the quadrature
					* points.
					*
					* There is one vector per
					* selected finite element
					* function, containing one
					* vector for each component,
					* containing vectors with
					* values for each quadrature
					* point.
					*/
     std::vector<std::vector<std::vector<double> > > values;

				       /**
					* The vector containing the
					* derivatives of finite
					* element functions in the
					* quadrature points.
					*
					* There is one vector per
					* selected finite element
					* function, containing one
					* vector for each component,
					* containing vectors with
					* values for each quadrature
					* point.
					*/
      std::vector<std::vector<std::vector<Tensor<1,dim> > > > gradients;

				       /**
					* The vector containing the
					* second derivatives of finite
					* element functions in the
					* quadrature points.
					*
					* There is one vector per
					* selected finite element
					* function, containing one
					* vector for each component,
					* containing vectors with
					* values for each quadrature
					* point.
					*/
      std::vector<std::vector<std::vector<Tensor<2,dim> > > > hessians;

				       /**
					* Reinitialize internal data
					* structures for use on a cell.
					*/
      void reinit(const DoFInfo<dim, spacedim>& i);
      
				       /**
					* Use the finite element
					* functions in #global_data
					* and fill the vectors
					* #values, #gradients and
					* #hessians.
					*/
      void fill_local_data(const DoFInfo<dim, spacedim>& info, const bool split_fevalues);

				       /**
					* The global data vector
					* used to compute function
					* values in quadrature
					* points.
					*/
      boost::shared_ptr<VectorDataBase<dim, spacedim> > global_data;

    private:
				       /**
					* Use the finite element
					* functions in #global_data
					* and fill the vectors
					* #values, #gradients and
					* #hessians with values
					* according to the
					* selector.
					*/
      template <typename TYPE>
      void fill_local_data(
	std::vector<std::vector<std::vector<TYPE> > >& data,
	VectorSelector& selector,
	bool split_fevalues) const;
				       /**
					* Cache the number of
					* components of the system element.
					*/
      unsigned int n_components;
  };

/**
 * A simple container collecting the five info objects required by the
 * integration loops. In addition to providing these objects, this
 * class offers two functions to initialize them with reasonable
 * contents all at once.
 *
 * MeshWorker::loop() requires five info objects collected into a
 * single class. These are the data members #cell, #boundary, #face,
 * #subface, and #neighbor. The names of those are expected by
 * MeshWorker::loop().
 *
 * @ingroup MeshWorker
 * @author Guido Kanschat, 2009
 */
  template <int dim, int spacedim=dim>
  class IntegrationInfoBox
  {
    public:

/// The type of the info object for cells
      typedef IntegrationInfo<dim, FEValuesBase<dim, spacedim>, spacedim> CellInfo;

/// The type of the info objects for faces.
      typedef IntegrationInfo<dim, FEFaceValuesBase<dim, spacedim>, spacedim> FaceInfo;
      
      template <class WORKER>
      void initialize(const WORKER&,
		      const FiniteElement<dim, spacedim>& el,
		      const Mapping<dim, spacedim>& mapping);

      template <class WORKER, typename VECTOR>
      void initialize(const WORKER&,
		      const FiniteElement<dim, spacedim>& el,
		      const Mapping<dim, spacedim>& mapping,
		      const NamedData<VECTOR*>& data);

      boost::shared_ptr<MeshWorker::VectorDataBase<dim, spacedim> > cell_data;
      boost::shared_ptr<MeshWorker::VectorDataBase<dim, spacedim> > boundary_data;
      boost::shared_ptr<MeshWorker::VectorDataBase<dim, spacedim> > face_data;

/// The info object for a cell
      CellInfo cell;
/// The info object for a boundary face
      FaceInfo boundary;
/// The info object for a regular interior face, seen from the first cell
      FaceInfo face;
/// The info object for the refined side of an interior face seen from the first cell
      FaceInfo subface;
/// The info object for an interior face, seen from the other cell
      FaceInfo neighbor;
  };


//----------------------------------------------------------------------//

  template <typename number>
  inline void
  LocalResults<number>::initialize_numbers(unsigned int n)
  {
    J.resize(n);
  }


  template <typename number>
  inline void
  LocalResults<number>::initialize_vectors(const unsigned int n)
  {
    R.resize(n);
  }


  template <typename number>
  template <class MatrixPtr>
  inline void
  LocalResults<number>::initialize_matrices(
    const std::vector<MatrixPtr>& matrices,
    bool both)
  {
    M1.resize(matrices.size());
    if (both)
      M2.resize(matrices.size());
    for (unsigned int i=0;i<matrices.size();++i)
      {
	const unsigned int row = matrices[i]->row;
	const unsigned int col = matrices[i]->column;

	M1[i].row = row;
	M1[i].column = col;
	if (both)
	  {
	    M2[i].row = row;
	    M2[i].column = col;
	  }
      }
  }


  template <typename number>
  inline void
  LocalResults<number>::initialize_matrices(const unsigned int n,
					    const bool both)
  {
    M1.resize(n);
    if (both)
      M2.resize(n);
    for (unsigned int i=0;i<n;++i)
      {
	M1[i].row = 0;
	M1[i].column = 0;
	if (both)
	  {
	    M2[i].row = 0;
	    M2[i].column = 0;
	  }
      }
  }


  template <typename number>
  inline void
  LocalResults<number>::reinit(const BlockIndices& bi)
  {
    for (unsigned int i=0;i<J.size();++i)
      J[i] = 0.;
    for (unsigned int i=0;i<R.size();++i)
      R[i].reinit(bi);
    for (unsigned int i=0;i<M1.size();++i)
      M1[i].matrix.reinit(bi.block_size(M1[i].row),
			  bi.block_size(M1[i].column));
    for (unsigned int i=0;i<M2.size();++i)
      M2[i].matrix.reinit(bi.block_size(M2[i].row),
			  bi.block_size(M2[i].column));
  }


  template <typename number>
  inline
  unsigned int
  LocalResults<number>::n_values() const
  {
    return J.size();
  }
  

  template <typename number>
  inline
  unsigned int
  LocalResults<number>::n_vectors() const
  {
    return R.size();
  }
  

  template <typename number>
  inline
  unsigned int
  LocalResults<number>::n_matrices() const
  {
    return M1.size();
  }
  

  template <typename number>
  inline
  number&
  LocalResults<number>::value(unsigned int i)
  {
    AssertIndexRange(i,J.size());
    return J[i];
  }
  

  template <typename number>
  inline
  BlockVector<number>&
  LocalResults<number>::vector(unsigned int i)
  {
    AssertIndexRange(i,R.size());
    return R[i];
  }
  
  
  template <typename number>
  inline
  MatrixBlock<FullMatrix<number> >&
  LocalResults<number>::matrix(unsigned int i, bool external)
  {
    if (external)
      {
	AssertIndexRange(i,M2.size());
	return M2[i];
      }
    AssertIndexRange(i,M1.size());
    return M1[i];
  }
  
  template <typename number>
  inline
  number
  LocalResults<number>::value(unsigned int i) const
  {
    AssertIndexRange(i,J.size());
    return J[i];
  }
  

  template <typename number>
  inline
  const BlockVector<number>&
  LocalResults<number>::vector(unsigned int i) const
  {
    AssertIndexRange(i,R.size());
    return R[i];
  }
  
  
  template <typename number>
  inline
  const MatrixBlock<FullMatrix<number> >&
  LocalResults<number>::matrix(unsigned int i, bool external) const
  {
    if (external)
      {
	AssertIndexRange(i,M2.size());
	return M2[i];
      }
    AssertIndexRange(i,M1.size());
    return M1[i];
  }
  
//----------------------------------------------------------------------//

  template <int dim, int spacedim>
  template <class DH>
  DoFInfo<dim,spacedim>::DoFInfo(const DH& dof_handler)
  {
    std::vector<unsigned int> aux(1);
    aux[0] = dof_handler.get_fe().dofs_per_cell;
    aux_local_indices.reinit(aux);
  }


  template <int dim, int spacedim>
  template <class DHCellIterator>
  inline void
  DoFInfo<dim,spacedim>::reinit(const DHCellIterator& c)
  {
    get_indices(c);
    cell = static_cast<typename Triangulation<dim,spacedim>::cell_iterator> (c);
    face_number = deal_II_numbers::invalid_unsigned_int;
    sub_number = deal_II_numbers::invalid_unsigned_int;
    if (block_info)
      LocalResults<double>::reinit(block_info->local());
    else
      LocalResults<double>::reinit(aux_local_indices);
  }


  template<int dim, int spacedim>
  template <class DHCellIterator, class DHFaceIterator>
  inline void
  DoFInfo<dim, spacedim>::reinit(const DHCellIterator& c,
				 const DHFaceIterator& f,
				 unsigned int n)
  {
    if ((cell.state() != IteratorState::valid)
	||  cell != static_cast<typename Triangulation<dim>::cell_iterator> (c))
      get_indices(c);
    cell = static_cast<typename Triangulation<dim>::cell_iterator> (c);
    face = static_cast<typename Triangulation<dim>::face_iterator> (f);
    face_number = n;
    sub_number = deal_II_numbers::invalid_unsigned_int;
    if (block_info)
      LocalResults<double>::reinit(block_info->local());
    else
      LocalResults<double>::reinit(aux_local_indices);
  }


  template<int dim, int spacedim>
  template <class DHCellIterator, class DHFaceIterator>
  inline void
  DoFInfo<dim, spacedim>::reinit(const DHCellIterator& c,
				 const DHFaceIterator& f,
				 unsigned int n,
				 const unsigned int s)
  {
    if (cell.state() != IteratorState::valid
	|| cell != static_cast<typename Triangulation<dim>::cell_iterator> (c))
      get_indices(c);
    cell = static_cast<typename Triangulation<dim>::cell_iterator> (c);
    face = static_cast<typename Triangulation<dim>::face_iterator> (f);
    face_number = n;
    sub_number = s;
    if (block_info)
      LocalResults<double>::reinit(block_info->local());
    else
      LocalResults<double>::reinit(aux_local_indices);
  }


  template<int dim, int spacedim>
  inline const BlockIndices&
  DoFInfo<dim, spacedim>::local_indices() const
  {
    if (block_info)
      return block_info->local();
    return aux_local_indices;
  }

//----------------------------------------------------------------------//
  
  template < int dim, int spacedim>
  inline
  DoFInfoBox<dim, spacedim>::DoFInfoBox(const DoFInfo<dim, spacedim>& seed)
		  :
		  cell(seed)
  {
    for (unsigned int i=0;i<GeometryInfo<dim>::faces_per_cell;++i)
      {
	exterior[i] = seed;
	interior[i] = seed;
	interior_face_available[i] = false;
	exterior_face_available[i] = false;
      }
  }


  template < int dim, int spacedim>
  inline void
  DoFInfoBox<dim, spacedim>::reset ()
  {
    for (unsigned int i=0;i<GeometryInfo<dim>::faces_per_cell;++i)
      {
    	interior_face_available[i] = false;
	exterior_face_available[i] = false;
      }
  }
  
//----------------------------------------------------------------------//

  template <int dim, class FVB, int spacedim>
  inline const FVB&
  IntegrationInfo<dim,FVB,spacedim>::fe_values() const
  {
    AssertDimension(fevalv.size(), 1);
    return *fevalv[0];
  }


  template <int dim, class FVB, int spacedim>
  inline const FVB&
  IntegrationInfo<dim,FVB,spacedim>::fe_values(unsigned int i) const
  {
    Assert (i<fevalv.size(), ExcIndexRange(i,0,fevalv.size()));
    return *fevalv[i];
  }


  template <int dim, class FVB, int spacedim>
  inline void
  IntegrationInfo<dim,FVB,spacedim>::reinit(const DoFInfo<dim, spacedim>& info)
  {
    for (unsigned int i=0;i<fevalv.size();++i)
      {
	FVB& febase = *fevalv[i];
	if (info.sub_number != deal_II_numbers::invalid_unsigned_int)
	  {
					     // This is a subface
	    FESubfaceValues<dim>& fe = dynamic_cast<FESubfaceValues<dim>&> (febase);
	    fe.reinit(info.cell, info.face_number, info.sub_number);
	  }
	else if (info.face_number != deal_II_numbers::invalid_unsigned_int)
	  {
					     // This is a face
	    FEFaceValues<dim>& fe = dynamic_cast<FEFaceValues<dim>&> (febase);
	    fe.reinit(info.cell, info.face_number);    
	  }
	else
	  {
					     // This is a cell
	    FEValues<dim>& fe = dynamic_cast<FEValues<dim>&> (febase);
	    fe.reinit(info.cell);
	  }
      }

    const bool split_fevalues = info.block_info != 0;
    fill_local_data(info, split_fevalues);
  }


//----------------------------------------------------------------------//

  template <int dim, int sdim>
  template <class WORKER>
  void
  IntegrationInfoBox<dim,sdim>::
  initialize(const WORKER& integrator,
	     const FiniteElement<dim,sdim>& el,
	     const Mapping<dim,sdim>& mapping)
  {
    cell.initialize<FEValues<dim,sdim> >(el, mapping, integrator.cell_quadrature,
			 integrator.cell_flags);
    boundary.initialize<FEFaceValues<dim,sdim> >(el, mapping, integrator.boundary_quadrature,
			 integrator.boundary_flags);
    face.initialize<FEFaceValues<dim,sdim> >(el, mapping, integrator.face_quadrature,
			 integrator.face_flags);
    subface.initialize<FESubfaceValues<dim,sdim> >(el, mapping, integrator.face_quadrature,
			    integrator.face_flags);
    neighbor.initialize<FEFaceValues<dim,sdim> >(el, mapping, integrator.face_quadrature,
						      integrator.neighbor_flags);
  }


  template <int dim, int sdim>
  template <class WORKER, typename VECTOR>
  void
  IntegrationInfoBox<dim,sdim>::initialize(
    const WORKER& integrator,
    const FiniteElement<dim,sdim>& el,
    const Mapping<dim,sdim>& mapping,
    const NamedData<VECTOR*>& data)
  {
    initialize(integrator, el, mapping);
    boost::shared_ptr<VectorData<VECTOR, dim, sdim> > p;

    p = boost::shared_ptr<VectorData<VECTOR, dim, sdim> >(new VectorData<VECTOR, dim, sdim> (integrator.cell_selector));
    p->initialize(data);
    cell_data = p;
    cell.initialize_data(p);

    p = boost::shared_ptr<VectorData<VECTOR, dim, sdim> >(new VectorData<VECTOR, dim, sdim> (integrator.boundary_selector));
    p->initialize(data);
    boundary_data = p;
    boundary.initialize_data(p);

    p = boost::shared_ptr<VectorData<VECTOR, dim, sdim> >(new VectorData<VECTOR, dim, sdim> (integrator.face_selector));
    p->initialize(data);
    face_data = p;
    face.initialize_data(p);
    subface.initialize_data(p);
    neighbor.initialize_data(p);
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
