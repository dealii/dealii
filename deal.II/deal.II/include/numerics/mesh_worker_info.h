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


/**
 * Basic info class only containing information on geometry and
 * degrees of freedom of the mesh object.
 *
 * The information in these objects is usually used by one of the
 * Assembler classes. It is also the kind of information which is
 * needed in mesh based matrices (often referred to as matrix free
 * methods).
 *
 * In addition to the information on degrees of freedom stored in this
 * class, it also provides the local computation space for the worker
 * object operating on it. This space is provided by the base class
 * template DATATYPE. This base class will automatically
 * reinitialized on each cell, but initial setup is up to the user and
 * should be done when initialize() for this class is called. The
 * currently available base classes are
 * <ul>
 * <li> LocalVectors
 * <li> LocalMatrices
 * </ul>
 *
 * This class operates in two different modes, corresponding to the
 * data models discussed in the Assembler namespace documentation.
 *
 * The choice of the local data model is triggered by the vector
 * #BlockInfo::local_renumbering, which in turn is usually filled by
 * BlockInfo::initialize_local(). If this function has been used, or
 * the vector has been changed from zero-length, then local dof
 * indices stored in this object will automatically be renumbered to
 * reflect local block structure.
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
      DoFInfo(const DH& dof_handler);

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
				       /// Fill index vector
      void get_indices(const typename DoFHandler<dim, spacedim>::cell_iterator c);

				       /// Fill index vector with level indices
      void get_indices(const typename MGDoFHandler<dim, spacedim>::cell_iterator c);

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
  class IntegrationInfo : public DoFInfo<dim, spacedim>
  {
    private:
				       /// vector of FEValues objects
      std::vector<boost::shared_ptr<FEVALUESBASE> > fevalv;
    public:
				       /**
					* Constructor forwarding
					* information to DoFInfo.
					*/
      IntegrationInfo(const BlockInfo& block_info);

				       /**
					* Constructor forwarding
					* information to DoFInfo.
					*/
      template <class DH>
      IntegrationInfo(const DH& dof_handler);

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
					*/
      template <class FEVALUES>
      void initialize(const FiniteElement<dim,spacedim>& el,
		      const Mapping<dim,spacedim>& mapping,
		      const Quadrature<FEVALUES::integral_dimension>& quadrature,
		      const UpdateFlags flags);

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
      const FEVALUESBASE& fe() const;

				       /// Access to finite elements
				       /**
					* This access function must
					* be used if the constructor
					* for a group of elements
					* was used.
					*
					* @see DGBlockSplitApplication
					*/
      const FEVALUESBASE& fe(const unsigned int i) const;

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
      template <class DHCellIterator>
      void reinit(const DHCellIterator& c);

				       /**
					* Reinitialize internal data
					* structures for use on a face.
					*/
      template <class DHCellIterator, class DHFaceIterator>
      void reinit(const DHCellIterator& c,
		  const DHFaceIterator& f,
		  const unsigned int fn);

				       /**
					* Reinitialize internal data
					* structures for use on a subface.
					*/
      template <class DHCellIterator, class DHFaceIterator>
      void reinit(const DHCellIterator& c,
		  const DHFaceIterator& f,
		  const unsigned int fn,
		  const unsigned int sn);


				       /**
					* @deprecated This is the
					* old version not using
					* VectorSelector.
					*
					* Use the finite element
					* functions in #global_data
					* and fill the vectors
					* #values, #gradients and
					* #hessians.
					*/
      void fill_local_data(const bool split_fevalues);

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
 * integration loops.
 *
 * @ingroup MeshWorker
 * @author Guido Kanschat, 2009
 */
  template <int dim, int spacedim=dim>
  class IntegrationInfoBox
  {
    public:
      typedef IntegrationInfo<dim, FEValuesBase<dim, spacedim>, spacedim> CellInfo;
      typedef IntegrationInfo<dim, FEFaceValuesBase<dim, spacedim>, spacedim> FaceInfo;

				       /**
					* Initialize all members
					* with the same argument.
					*/
      template <typename T>
      IntegrationInfoBox(const T&);

      template <class WORKER, class ASSEMBLER>
      void initialize(const WORKER&,
		      ASSEMBLER &assembler,
		      const FiniteElement<dim, spacedim>& el,
		      const Mapping<dim, spacedim>& mapping);

      template <class WORKER, typename VECTOR>
      void initialize(const WORKER&,
		      const FiniteElement<dim, spacedim>& el,
		      const Mapping<dim, spacedim>& mapping,
		      const NamedData<VECTOR*>& data);

//      private:

      boost::shared_ptr<MeshWorker::VectorDataBase<dim, spacedim> > cell_data;
      boost::shared_ptr<MeshWorker::VectorDataBase<dim, spacedim> > bdry_data;
      boost::shared_ptr<MeshWorker::VectorDataBase<dim, spacedim> > face_data;

      CellInfo cell_info;
      FaceInfo bdry_info;
      FaceInfo face_info;
      FaceInfo subface_info;
      FaceInfo neighbor_info;
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

  template <int dim, class FVB, int spacedim>
  template <class DH>
  IntegrationInfo<dim,FVB,spacedim>::IntegrationInfo(const DH& dof_handler)
		  :
		  DoFInfo<dim, spacedim>(dof_handler),
		  fevalv(0),
		  multigrid(false),
		  global_data(boost::shared_ptr<VectorDataBase<dim, spacedim> >(new VectorDataBase<dim, spacedim>))
  {}


  template <int dim, class FVB, int spacedim>
  inline const FVB&
  IntegrationInfo<dim,FVB,spacedim>::fe() const
  {
    AssertDimension(fevalv.size(), 1);
    return *fevalv[0];
  }


  template <int dim, class FVB, int spacedim>
  inline const FVB&
  IntegrationInfo<dim,FVB,spacedim>::fe(unsigned int i) const
  {
    Assert (i<fevalv.size(), ExcIndexRange(i,0,fevalv.size()));
    return *fevalv[i];
  }


  template <int dim, class FVB, int spacedim>
  template <class DHCellIterator>
  inline void
  IntegrationInfo<dim,FVB,spacedim>::reinit(const DHCellIterator& c)
  {
    DoFInfo<dim,spacedim>::reinit(c);
    for (unsigned int i=0;i<fevalv.size();++i)
      {
	FVB& febase = *fevalv[i];
	FEValues<dim>& fe = dynamic_cast<FEValues<dim>&> (febase);
	fe.reinit(this->cell);
      }

    const bool split_fevalues = this->block_info != 0;
    fill_local_data(split_fevalues);
  }


  template <int dim, class FVB, int spacedim>
  template <class DHCellIterator, class DHFaceIterator>
  inline void
  IntegrationInfo<dim,FVB,spacedim>::reinit(
    const DHCellIterator& c,
    const DHFaceIterator& f,
    const unsigned int fn)
  {
    DoFInfo<dim,spacedim>::reinit(c, f, fn);
    for (unsigned int i=0;i<fevalv.size();++i)
      {
	FVB& febase = *fevalv[i];
	FEFaceValues<dim>& fe = dynamic_cast<FEFaceValues<dim>&> (febase);
	fe.reinit(this->cell, fn);
      }

    const bool split_fevalues = this->block_info != 0;
    fill_local_data(split_fevalues);
  }


  template <int dim, class FVB, int spacedim>
  template <class DHCellIterator, class DHFaceIterator>
  inline void
  IntegrationInfo<dim,FVB,spacedim>::reinit(
    const DHCellIterator& c,
    const DHFaceIterator& f,
    const unsigned int fn,
    const unsigned int sn)
  {
    DoFInfo<dim,spacedim>::reinit(c, f, fn, sn);
    for (unsigned int i=0;i<fevalv.size();++i)
      {
	FVB& febase = *fevalv[i];
	FESubfaceValues<dim>& fe = dynamic_cast<FESubfaceValues<dim>&> (febase);
	fe.reinit(this->cell, fn, sn);
      }

    const bool split_fevalues = this->block_info != 0;
    fill_local_data(split_fevalues);
  }

//----------------------------------------------------------------------//

  template <int dim, int sdim>
  template <typename T>
  IntegrationInfoBox<dim,sdim>::IntegrationInfoBox(const T& t)
		  :
		  cell_info(t),
		  bdry_info(t),
		  face_info(t),
		  subface_info(t),
		  neighbor_info(t)
  {}


  template <int dim, int sdim>
  template <class WORKER, class ASSEMBLER>
  void
  IntegrationInfoBox<dim,sdim>::
  initialize(const WORKER& integrator,
	     ASSEMBLER &assembler,
	     const FiniteElement<dim,sdim>& el,
	     const Mapping<dim,sdim>& mapping)
  {
    assembler.initialize_info(cell_info, false);
    assembler.initialize_info(bdry_info, false);
    assembler.initialize_info(face_info, true);
    assembler.initialize_info(subface_info, true);
    assembler.initialize_info(neighbor_info, true);

    cell_info.initialize<FEValues<dim,sdim> >(el, mapping, integrator.cell_quadrature,
			 integrator.cell_flags);
    bdry_info.initialize<FEFaceValues<dim,sdim> >(el, mapping, integrator.bdry_quadrature,
			 integrator.bdry_flags);
    face_info.initialize<FEFaceValues<dim,sdim> >(el, mapping, integrator.face_quadrature,
			 integrator.face_flags);
    subface_info.initialize<FESubfaceValues<dim,sdim> >(el, mapping, integrator.face_quadrature,
			    integrator.face_flags);
    neighbor_info.initialize<FEFaceValues<dim,sdim> >(el, mapping, integrator.face_quadrature,
						      integrator.ngbr_flags);
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
    cell_info.initialize_data(p);

    p = boost::shared_ptr<VectorData<VECTOR, dim, sdim> >(new VectorData<VECTOR, dim, sdim> (integrator.bdry_selector));
    p->initialize(data);
    bdry_data = p;
    bdry_info.initialize_data(p);

    p = boost::shared_ptr<VectorData<VECTOR, dim, sdim> >(new VectorData<VECTOR, dim, sdim> (integrator.face_selector));
    p->initialize(data);
    face_data = p;
    face_info.initialize_data(p);
    subface_info.initialize_data(p);
    neighbor_info.initialize_data(p);
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
