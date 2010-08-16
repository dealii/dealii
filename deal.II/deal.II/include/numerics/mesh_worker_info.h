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
#include <base/quadrature_lib.h>
#include <base/std_cxx1x/shared_ptr.h>
#include <dofs/block_info.h>
#include <fe/fe_values.h>
#include <numerics/mesh_worker.h>
#include <numerics/mesh_worker_vector_selector.h>

DEAL_II_NAMESPACE_OPEN

namespace MeshWorker
{
  template <int dim, class DOFINFO> class DoFInfoBox;
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
  template<int dim, int spacedim = dim, typename number = double>
  class DoFInfo : public LocalResults<number>
  {
    public:
				       /// The current cell
      typename Triangulation<dim, spacedim>::cell_iterator cell;

				       /// The current face
      typename Triangulation<dim, spacedim>::face_iterator face;

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
					* Constructor
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

      bool level_cell;
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

      friend class DoFInfoBox<dim, DoFInfo<dim, spacedim, number> >;
  };


  /**
 * A class bundling the MeshWorker::DoFInfo objects used on a cell.
 *
 * @todo Currently, we are storing an object for the cells and two for
 * each face. We could gather all face data pertaining to the cell
 * itself in one object, saving a bit of memory and a few operations,
 * but sacrificing some cleanliness.
 *
 * @ingroup MeshWorker
 * @author Guido Kanschat, 2010
 */
  template <int dim, class DOFINFO>
  class DoFInfoBox
  {
    public:
				       /**
					* Constructor copying the seed
					* into all other objects.
					*/
      DoFInfoBox(const DOFINFO& seed);

				       /**
					* Copy constructor, taking
					* #cell and using it as a seed
					* in the other constructor.
					*/
      DoFInfoBox(const DoFInfoBox<dim, DOFINFO>&);

				       /**
					* Reset all the availability flags.
					*/
      void reset();

				       /**
					* After all info objects have
					* been filled appropriately,
					* use the ASSEMBLER object
					* to assemble them into the
					* global data. See
					* MeshWorker::Assembler for
					* available classes.
					*/
      template <class ASSEMBLER>
      void assemble(ASSEMBLER& ass) const;


				       /**
					* The data for the cell.
					*/
      DOFINFO cell;
				       /**
					* The data for the faces from inside.
					*/
      DOFINFO interior[GeometryInfo<dim>::faces_per_cell];
				       /**
					* The data for the faces from outside.
					*/
      DOFINFO exterior[GeometryInfo<dim>::faces_per_cell];

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
  template<int dim, int spacedim = dim>
  class IntegrationInfo
  {
    private:
				       /// vector of FEValues objects
      std::vector<std_cxx1x::shared_ptr<FEValuesBase<dim, spacedim> > > fevalv;
    public:
      static const unsigned int dimension = dim;
      static const unsigned int space_dimension = spacedim;

				       /**
					* Constructor.
					*/
      IntegrationInfo();

				       /**
					* Copy constructor, creating a
					* clone to be used by
					* WorksTream::run().
					*/
      IntegrationInfo(const IntegrationInfo<dim, spacedim>& other);

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
      void initialize_data(const std_cxx1x::shared_ptr<VectorDataBase<dim,spacedim> > &data);

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
      const FEValuesBase<dim, spacedim>& fe_values () const;

				       /// Access to finite elements
				       /**
					* This access function must
					* be used if the constructor
					* for a group of elements
					* was used.
					*
					* @see DGBlockSplitApplication
					*/
      const FEValuesBase<dim, spacedim>& fe_values (const unsigned int i) const;

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
      template<typename number>
      void fill_local_data(const DoFInfo<dim, spacedim, number>& info, bool split_fevalues);

				       /**
					* The global data vector
					* used to compute function
					* values in quadrature
					* points.
					*/
      std_cxx1x::shared_ptr<VectorDataBase<dim, spacedim> > global_data;

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
 * The object holding the scratch data for integrating over cells and
 * faces. IntegrationInfoBox serves three main purposes:
 *
 * <ol>
 * <li> It provides the interface needed by MeshWorker::loop(), namely
 * the two functions post_cell() and post_faces(), as well as
 * the data members #cell, #boundary, #face,
 * #subface, and #neighbor.
 *
 * <li> It contains all information needed to initialize the FEValues
 * and FEFaceValues objects in the IntegrationInfo data members.
 *
 * <li> It stores information on finite element vectors and whether
 * their data should be used to compute values or derivatives of
 * functions at quadrature points.
 * </ol>
 *
 * In order to allow for sufficient generality, a few steps have to be
 * undertaken to use this class.
 *
 * First, you should consider if you need values from any vectors in a
 * NamedData object. If so, fill the VectorSelector objects
 * #cell_selector, #boundary_selector and #face_selector with their names
 * and the data type (value, gradient, Hessian) to be extracted.
 *
 * Afterwards, you will need to consider UpdateFlags for FEValues
 * objects. A good start is initialize_update_flags(), which looks at
 * the selectors filled before and adds all the flags needed to get
 * the selection. Additional flags can be set with add_update_flags().
 *
 * Finally, we need to choose quadrature formulas. If you choose to
 * use Gauss formulas only, use initialize_gauss_quadrature() with
 * appropriate values. Otherwise, you can fill the variables
 * #cell_quadrature, #boundary_quadrature and #face_quadrature directly.
 *
 * In order to save time, you can set the variables boundary_fluxes
 * and interior_fluxes of the base class to false, thus telling the
 * Meshworker::loop() not to loop over those faces.
 *
 * All the information in here is used to set up IntegrationInfo
 * objects correctly, typically in an IntegrationInfoBox.
 *
 * @ingroup MeshWorker
 * @author Guido Kanschat, 2009
 */
  template <int dim, int spacedim=dim>
  class IntegrationInfoBox
  {
    public:

/// The type of the info object for cells
      typedef IntegrationInfo<dim, spacedim> CellInfo;

      void initialize(const FiniteElement<dim, spacedim>& el,
		      const Mapping<dim, spacedim>& mapping,
		      const BlockInfo* block_info = 0);

      template <typename VECTOR>
      void initialize(const FiniteElement<dim, spacedim>& el,
		      const Mapping<dim, spacedim>& mapping,
		      const NamedData<VECTOR*>& data,
		      const BlockInfo* block_info = 0);

      template <typename VECTOR>
      void initialize(const FiniteElement<dim, spacedim>& el,
		      const Mapping<dim, spacedim>& mapping,
		      const NamedData<MGLevelObject<VECTOR>*>& data,
		      const BlockInfo* block_info = 0);
				       /**
					* @name FEValues setup
					*/
				       /* @{ */
      void initialize_update_flags();

				       /**
					* Add additional values for update.
					*/
      void add_update_flags(const UpdateFlags flags, bool cell = true,
			    bool boundary = true, bool face = true,
			    bool neighbor = true);

				       /** Assign n-point Gauss
					* quadratures to each of the
					* quadrature rules. Here, a
					* size of zero points means
					* that no loop over these grid
					* entities should be
					* performed.
					*/
      void initialize_gauss_quadrature(unsigned int n_cell_points,
				       unsigned int n_boundary_points,
				       unsigned int n_face_points);

				       /**
					* The set of update flags
					* for boundary cell integration.
					*
					* Defaults to
					* #update_JxW_values.
					*/
      UpdateFlags cell_flags;
				       /**
					* The set of update flags
					* for boundary face integration.
					*
					* Defaults to
					* #update_JxW_values and
					* #update_normal_vectors.
					*/
      UpdateFlags boundary_flags;

				       /**
					* The set of update flags
					* for interior face integration.
					*
					* Defaults to
					* #update_JxW_values and
					* #update_normal_vectors.
					*/
      UpdateFlags face_flags;

				       /**
					* The set of update flags
					* for interior face integration.
					*
					* Defaults to
					* #update_default, since
					* quadrature weights are
					* taken from the other cell.
					*/
      UpdateFlags neighbor_flags;

				       /**
					* The quadrature rule used
					* on cells.
					*/
      Quadrature<dim> cell_quadrature;

				       /**
					* The quadrature rule used
					* on boundary faces.
					*/
      Quadrature<dim-1> boundary_quadrature;

				       /**
					* The quadrature rule used
					* on interior faces.
					*/
      Quadrature<dim-1> face_quadrature;
				       /* @} */

				       /**
					* @name Data vectors
					*/
      				       /* @{ */

				       /**
					* Initialize the
					* VectorSelector objects
					* #cell_selector,
					* #boundary_selector and
					* #face_selector in order to
					* save computational
					* eeffort. If no selectors
					* are used, then values for
					* all named vectors in
					* DoFInfo::global_data will be
					* computed in all quadrature
					* points.
					*
					* This function will also
					* add UpdateFlags to the
					* flags stored in this class.
					*/
				       /**
					* Select the vectors from
					* DoFInfo::global_data
					* that should be computed in
					* the quadrature points on cells.
					*/
      MeshWorker::VectorSelector cell_selector;

				       /**
					* Select the vectors from
					* DoFInfo::global_data
					* that should be computed in
					* the quadrature points on
					* boundary faces.
					*/
      MeshWorker::VectorSelector boundary_selector;

				       /**
					* Select the vectors from
					* DoFInfo::global_data
					* that should be computed in
					* the quadrature points on
					* interior faces.
					*/
      MeshWorker::VectorSelector face_selector;

      std_cxx1x::shared_ptr<MeshWorker::VectorDataBase<dim, spacedim> > cell_data;
      std_cxx1x::shared_ptr<MeshWorker::VectorDataBase<dim, spacedim> > boundary_data;
      std_cxx1x::shared_ptr<MeshWorker::VectorDataBase<dim, spacedim> > face_data;
				       /* @} */

				       /**
					* @name Interface for MeshWorker::loop()
					*/
				       /* @{ */
				       /**
					* A callback function which is
					* called in the loop over all
					* cells, after the action on a
					* cell has been performed and
					* before the faces are dealt
					* with.
					*
					* In order for this function
					* to have this effect,
					* at least either of the
					* arguments
					* <tt>boundary_worker</tt> or
					* <tt>face_worker</tt>
					* arguments of loop() should
					* be nonzero. Additionally,
					* <tt>cells_first</tt> should
					* be true. If
					* <tt>cells_first</tt> is
					* false, this function is
					* called before any action on
					* a cell is taken.
					*
					* And empty function in this
					* class, but can be replaced
					* in other classes given to
					* loop() instead.
					*
					* See loop() and cell_action()
					* for more details of how this
					* function can be used.
					*/
      template <class DOFINFO>
      void post_cell(const DoFInfoBox<dim, DOFINFO>&);

				       /**
					* A callback function which is
					* called in the loop over all
					* cells, after the action on
					* the faces of a cell has been
					* performed and before the
					* cell itself is dealt with
					* (assumes
					* <tt>cells_first</tt> is false).
					*
					* In order for this function
					* to have a reasonable effect,
					* at least either of the
					* arguments
					* <tt>boundary_worker</tt> or
					* <tt>face_worker</tt>
					* arguments of loop() should
					* be nonzero. Additionally,
					* <tt>cells_first</tt> should
					* be false.
					*
					* And empty function in this
					* class, but can be replaced
					* in other classes given to
					* loop() instead.
					*
					* See loop() and cell_action()
					* for more details of how this
					* function can be used.
					*/
      template <class DOFINFO>
      void post_faces(const DoFInfoBox<dim, DOFINFO>&);

/// The info object for a cell
      CellInfo cell;
/// The info object for a boundary face
      CellInfo boundary;
/// The info object for a regular interior face, seen from the first cell
      CellInfo face;
/// The info object for the refined side of an interior face seen from the first cell
      CellInfo subface;
/// The info object for an interior face, seen from the other cell
      CellInfo neighbor;

				       /* @} */
  };


//----------------------------------------------------------------------//

  template <int dim, int spacedim, typename number>
  template <class DH>
  DoFInfo<dim,spacedim,number>::DoFInfo(const DH& dof_handler)
  {
    std::vector<unsigned int> aux(1);
    aux[0] = dof_handler.get_fe().dofs_per_cell;
    aux_local_indices.reinit(aux);
  }


  template <int dim, int spacedim, typename number>
  template <class DHCellIterator>
  inline void
  DoFInfo<dim,spacedim,number>::reinit(const DHCellIterator& c)
  {
    get_indices(c);
    cell = static_cast<typename Triangulation<dim,spacedim>::cell_iterator> (c);
    face_number = deal_II_numbers::invalid_unsigned_int;
    sub_number = deal_II_numbers::invalid_unsigned_int;
    if (block_info)
      LocalResults<number>::reinit(block_info->local());
    else
      LocalResults<number>::reinit(aux_local_indices);
  }


  template<int dim, int spacedim, typename number>
  template <class DHCellIterator, class DHFaceIterator>
  inline void
  DoFInfo<dim,spacedim,number>::reinit(const DHCellIterator& c,
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
      LocalResults<number>::reinit(block_info->local());
    else
      LocalResults<number>::reinit(aux_local_indices);
  }


  template<int dim, int spacedim, typename number>
  template <class DHCellIterator, class DHFaceIterator>
  inline void
  DoFInfo<dim,spacedim,number>::reinit(const DHCellIterator& c,
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
      LocalResults<number>::reinit(block_info->local());
    else
      LocalResults<number>::reinit(aux_local_indices);
  }


  template<int dim, int spacedim, typename number>
  inline const BlockIndices&
  DoFInfo<dim,spacedim,number>::local_indices() const
  {
    if (block_info)
      return block_info->local();
    return aux_local_indices;
  }

//----------------------------------------------------------------------//

  template <int dim, class DOFINFO>
  inline
  DoFInfoBox<dim, DOFINFO>::DoFInfoBox(const DOFINFO& seed)
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


  template <int dim, class DOFINFO>
  inline
  DoFInfoBox<dim, DOFINFO>::DoFInfoBox(const DoFInfoBox<dim, DOFINFO>& other)
		  :
		  cell(other.cell)
  {
    for (unsigned int i=0;i<GeometryInfo<dim>::faces_per_cell;++i)
      {
	exterior[i] = other.exterior[i];
	interior[i] = other.interior[i];
	interior_face_available[i] = false;
	exterior_face_available[i] = false;
      }
  }


  template <int dim, class DOFINFO>
  inline void
  DoFInfoBox<dim, DOFINFO>::reset ()
  {
    for (unsigned int i=0;i<GeometryInfo<dim>::faces_per_cell;++i)
      {
    	interior_face_available[i] = false;
	exterior_face_available[i] = false;
      }
  }


  template <int dim, class DOFINFO>
  template <class ASSEMBLER>
  inline void
  DoFInfoBox<dim, DOFINFO>::assemble (ASSEMBLER& assembler) const
  {
    assembler.assemble(cell);
    for (unsigned int i=0;i<GeometryInfo<dim>::faces_per_cell;++i)
      {
					 // Only do something if data available
	if (interior_face_available[i])
	  {
					     // If both data
					     // available, it is an
					     // interior face
	    if (exterior_face_available[i])
	      assembler.assemble(interior[i], exterior[i]);
	    else
	      assembler.assemble(interior[i]);
	  }
      }
  }


//----------------------------------------------------------------------//

  template<int dim, int sdim>
  inline
  IntegrationInfo<dim,sdim>::IntegrationInfo()
		  :
		  fevalv(0),
		  multigrid(false),
		  global_data(std_cxx1x::shared_ptr<VectorDataBase<dim, sdim> >(new VectorDataBase<dim, sdim>))
  {}


  template<int dim, int sdim>
  inline
  IntegrationInfo<dim,sdim>::IntegrationInfo(const IntegrationInfo<dim,sdim>& other)
		  :
		  multigrid(other.multigrid),
		  values(other.values),
		  gradients(other.gradients),
		  hessians(other.hessians),
		  global_data(other.global_data),
		  n_components(other.n_components)
  {
    fevalv.resize(other.fevalv.size());
    for (unsigned int i=0;i<other.fevalv.size();++i)
      {
	const FEValuesBase<dim,sdim>& p = *other.fevalv[i];
	const FEValues<dim,sdim>* pc = dynamic_cast<const FEValues<dim,sdim>*>(&p);
	const FEFaceValues<dim,sdim>* pf = dynamic_cast<const FEFaceValues<dim,sdim>*>(&p);
	const FESubfaceValues<dim,sdim>* ps = dynamic_cast<const FESubfaceValues<dim,sdim>*>(&p);

	if (pc != 0)
	  fevalv[i] = std_cxx1x::shared_ptr<FEValuesBase<dim,sdim> > (
	    reinterpret_cast<FEFaceValuesBase<dim,sdim>*>(
	      new FEValues<dim,sdim> (pc->get_mapping(), pc->get_fe(),
				      pc->get_quadrature(), pc->get_update_flags())));
	else if (pf != 0)
	  fevalv[i] = std_cxx1x::shared_ptr<FEValuesBase<dim,sdim> > (
	    new FEFaceValues<dim,sdim> (pf->get_mapping(), pf->get_fe(), pf->get_quadrature(), pf->get_update_flags()));
	else if (ps != 0)
	  fevalv[i] = std_cxx1x::shared_ptr<FEValuesBase<dim,sdim> > (
	    new FESubfaceValues<dim,sdim> (ps->get_mapping(), ps->get_fe(), ps->get_quadrature(), ps->get_update_flags()));
	else
	  Assert(false, ExcInternalError());
      }
  }


  template <>
  inline
  IntegrationInfo<1,1>::IntegrationInfo(const IntegrationInfo<1,1>& other)
		  :
		  multigrid(other.multigrid),
		  values(other.values),
		  gradients(other.gradients),
		  hessians(other.hessians),
		  global_data(other.global_data),
		  n_components(other.n_components)
  {
    const int dim = 1;
    const int sdim = 1;

    fevalv.resize(other.fevalv.size());
    for (unsigned int i=0;i<other.fevalv.size();++i)
      {
	const FEValuesBase<dim,sdim>& p = *other.fevalv[i];
	const FEValues<dim,sdim>* pc = dynamic_cast<const FEValues<dim,sdim>*>(&p);
	const FEFaceValues<dim,sdim>* pf = dynamic_cast<const FEFaceValues<dim,sdim>*>(&p);
	const FESubfaceValues<dim,sdim>* ps = dynamic_cast<const FESubfaceValues<dim,sdim>*>(&p);

	if (pc != 0)
	  fevalv[i] = std_cxx1x::shared_ptr<FEValuesBase<dim,sdim> > (
	    reinterpret_cast<FEFaceValuesBase<dim,sdim>*>(
	      new FEValues<dim,sdim> (pc->get_mapping(), pc->get_fe(),
				      pc->get_quadrature(), pc->get_update_flags())));
	else if (pf != 0)
	  {
	    Assert (false, ExcImpossibleInDim(1));
	    fevalv[i].reset ();
	  }
	else if (ps != 0)
	  {
	    Assert (false, ExcImpossibleInDim(1));
	    fevalv[i].reset();
	  }
	else
	  Assert(false, ExcInternalError());
      }
  }


  template <>
  inline
  IntegrationInfo<1,2>::IntegrationInfo(const IntegrationInfo<1,2>& other)
		  :
		  multigrid(other.multigrid),
		  values(other.values),
		  gradients(other.gradients),
		  hessians(other.hessians),
		  global_data(other.global_data),
		  n_components(other.n_components)
  {
    const int dim = 1;
    const int sdim = 2;

    fevalv.resize(other.fevalv.size());
    for (unsigned int i=0;i<other.fevalv.size();++i)
      {
	const FEValuesBase<dim,sdim>& p = *other.fevalv[i];
	const FEValues<dim,sdim>* pc = dynamic_cast<const FEValues<dim,sdim>*>(&p);
	const FEFaceValues<dim,sdim>* pf = dynamic_cast<const FEFaceValues<dim,sdim>*>(&p);
	const FESubfaceValues<dim,sdim>* ps = dynamic_cast<const FESubfaceValues<dim,sdim>*>(&p);

	if (pc != 0)
	  fevalv[i] = std_cxx1x::shared_ptr<FEValuesBase<dim,sdim> > (
	    reinterpret_cast<FEFaceValuesBase<dim,sdim>*>(
	      new FEValues<dim,sdim> (pc->get_mapping(), pc->get_fe(),
				      pc->get_quadrature(), pc->get_update_flags())));
	else if (pf != 0)
	  {
	    Assert (false, ExcImpossibleInDim(1));
	    fevalv[i].reset();
	  }
	else if (ps != 0)
	  {
	    Assert (false, ExcImpossibleInDim(1));
	    fevalv[i].reset();
	  }
	else
	  Assert(false, ExcInternalError());
      }
  }



  template<int dim, int sdim>
  template <class FEVALUES>
  inline void
  IntegrationInfo<dim,sdim>::initialize(
    const FiniteElement<dim,sdim>& el,
    const Mapping<dim,sdim>& mapping,
    const Quadrature<FEVALUES::integral_dimension>& quadrature,
    const UpdateFlags flags,
    const BlockInfo* block_info)
  {
    if (block_info == 0 || block_info->local().size() == 0)
      {
	fevalv.resize(1);
	fevalv[0] = std_cxx1x::shared_ptr<FEValuesBase<dim,sdim> > (
	  new FEVALUES (mapping, el, quadrature, flags));
      }
    else
      {
	fevalv.resize(el.n_base_elements());
	for (unsigned int i=0;i<fevalv.size();++i)
	  {
	    fevalv[i] = std_cxx1x::shared_ptr<FEValuesBase<dim,sdim> > (
	      new FEVALUES (mapping, el.base_element(i), quadrature, flags));
	  }
      }
    n_components = el.n_components();
  }


  template <int dim, int spacedim>
  inline const FEValuesBase<dim, spacedim>&
  IntegrationInfo<dim,spacedim>::fe_values() const
  {
    AssertDimension(fevalv.size(), 1);
    return *fevalv[0];
  }


  template <int dim, int spacedim>
  inline const FEValuesBase<dim, spacedim>&
  IntegrationInfo<dim,spacedim>::fe_values(unsigned int i) const
  {
    Assert (i<fevalv.size(), ExcIndexRange(i,0,fevalv.size()));
    return *fevalv[i];
  }


  template <int dim, int spacedim>
  inline void
  IntegrationInfo<dim,spacedim>::reinit(const DoFInfo<dim, spacedim>& info)
  {
    for (unsigned int i=0;i<fevalv.size();++i)
      {
	FEValuesBase<dim, spacedim>& febase = *fevalv[i];
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
	    FEValues<dim,spacedim>& fe = dynamic_cast<FEValues<dim,spacedim>&> (febase);
	    fe.reinit(info.cell);
	  }
      }

    const bool split_fevalues = info.block_info != 0;
    if (!global_data->empty())
      fill_local_data(info, split_fevalues);
  }


  template <>
  inline void
  IntegrationInfo<1,1>::reinit(const DoFInfo<1,1>& info)
  {
    const int dim = 1;
    const int spacedim = 1;

    for (unsigned int i=0;i<fevalv.size();++i)
      {
	FEValuesBase<dim, spacedim>& febase = *fevalv[i];
	if (info.sub_number != deal_II_numbers::invalid_unsigned_int)
	  {
					     // This is a subface
	    Assert (false, ExcImpossibleInDim(1));
	  }
	else if (info.face_number != deal_II_numbers::invalid_unsigned_int)
	  {
					     // This is a face
	    Assert (false, ExcImpossibleInDim(1));
	  }
	else
	  {
					     // This is a cell
	    FEValues<dim,spacedim>& fe = dynamic_cast<FEValues<dim,spacedim>&> (febase);
	    fe.reinit(info.cell);
	  }
      }

    const bool split_fevalues = info.block_info != 0;
    if (!global_data->empty())
      fill_local_data(info, split_fevalues);
  }


  template <>
  inline void
  IntegrationInfo<1,2>::reinit(const DoFInfo<1,2>& info)
  {
    const int dim = 1;
    const int spacedim = 2;

    for (unsigned int i=0;i<fevalv.size();++i)
      {
	FEValuesBase<dim, spacedim>& febase = *fevalv[i];
	if (info.sub_number != deal_II_numbers::invalid_unsigned_int)
	  {
					     // This is a subface
	    Assert (false, ExcImpossibleInDim(1));
	  }
	else if (info.face_number != deal_II_numbers::invalid_unsigned_int)
	  {
					     // This is a face
	    Assert (false, ExcImpossibleInDim(1));
	  }
	else
	  {
					     // This is a cell
	    FEValues<dim,spacedim>& fe = dynamic_cast<FEValues<dim,spacedim>&> (febase);
	    fe.reinit(info.cell);
	  }
      }

    const bool split_fevalues = info.block_info != 0;
    if (!global_data->empty())
      fill_local_data(info, split_fevalues);
  }


  template <int dim, int sdim>
  inline
  void
  IntegrationInfoBox<dim,sdim>::
  initialize(const FiniteElement<dim,sdim>& el,
	     const Mapping<dim,sdim>& mapping,
	     const BlockInfo* block_info)
  {
    cell.template initialize<FEValues<dim,sdim> >(el, mapping, cell_quadrature,
						  cell_flags, block_info);
    boundary.template initialize<FEFaceValues<dim,sdim> >(el, mapping, boundary_quadrature,
							  boundary_flags, block_info);
    face.template initialize<FEFaceValues<dim,sdim> >(el, mapping, face_quadrature,
						      face_flags, block_info);
    subface.template initialize<FESubfaceValues<dim,sdim> >(el, mapping, face_quadrature,
							    face_flags, block_info);
    neighbor.template initialize<FEFaceValues<dim,sdim> >(el, mapping, face_quadrature,
							  neighbor_flags, block_info);
  }


  template <>
  inline
  void
  IntegrationInfoBox<1,1>::
  initialize(const FiniteElement<1,1>& el,
	     const Mapping<1,1>& mapping,
	     const BlockInfo* block_info)
  {
    const int dim = 1;
    const int sdim = 1;

    cell.initialize<FEValues<dim,sdim> >(el, mapping, cell_quadrature,
					 cell_flags, block_info);
  }



  template <>
  inline
  void
  IntegrationInfoBox<1,2>::
  initialize(const FiniteElement<1,2>& el,
	     const Mapping<1,2>& mapping,
	     const BlockInfo* block_info)
  {
    const int dim = 1;
    const int sdim = 2;

    cell.initialize<FEValues<dim,sdim> >(el, mapping, cell_quadrature,
					 cell_flags, block_info);
  }


//----------------------------------------------------------------------//

  template <int dim, int sdim>
  template <typename VECTOR>
  void
  IntegrationInfoBox<dim,sdim>::initialize(
    const FiniteElement<dim,sdim>& el,
    const Mapping<dim,sdim>& mapping,
    const NamedData<VECTOR*>& data,
    const BlockInfo* block_info)
  {
    initialize(el, mapping, block_info);
    std_cxx1x::shared_ptr<VectorData<VECTOR, dim, sdim> > p;

    p = std_cxx1x::shared_ptr<VectorData<VECTOR, dim, sdim> >(new VectorData<VECTOR, dim, sdim> (cell_selector));
    p->initialize(data);
    cell_data = p;
    cell.initialize_data(p);

    p = std_cxx1x::shared_ptr<VectorData<VECTOR, dim, sdim> >(new VectorData<VECTOR, dim, sdim> (boundary_selector));
    p->initialize(data);
    boundary_data = p;
    boundary.initialize_data(p);

    p = std_cxx1x::shared_ptr<VectorData<VECTOR, dim, sdim> >(new VectorData<VECTOR, dim, sdim> (face_selector));
    p->initialize(data);
    face_data = p;
    face.initialize_data(p);
    subface.initialize_data(p);
    neighbor.initialize_data(p);
  }


  template <int dim, int sdim>
  template <typename VECTOR>
  void
  IntegrationInfoBox<dim,sdim>::initialize(
    const FiniteElement<dim,sdim>& el,
    const Mapping<dim,sdim>& mapping,
    const NamedData<MGLevelObject<VECTOR>*>& data,
    const BlockInfo* block_info)
  {
    initialize(el, mapping, block_info);
    std_cxx1x::shared_ptr<MGVectorData<VECTOR, dim, sdim> > p;

    p = std_cxx1x::shared_ptr<MGVectorData<VECTOR, dim, sdim> >(new MGVectorData<VECTOR, dim, sdim> (cell_selector));
    p->initialize(data);
    cell_data = p;
    cell.initialize_data(p);

    p = std_cxx1x::shared_ptr<MGVectorData<VECTOR, dim, sdim> >(new MGVectorData<VECTOR, dim, sdim> (boundary_selector));
    p->initialize(data);
    boundary_data = p;
    boundary.initialize_data(p);

    p = std_cxx1x::shared_ptr<MGVectorData<VECTOR, dim, sdim> >(new MGVectorData<VECTOR, dim, sdim> (face_selector));
    p->initialize(data);
    face_data = p;
    face.initialize_data(p);
    subface.initialize_data(p);
    neighbor.initialize_data(p);
  }


  template <int dim, int sdim>
  template <class DOFINFO>
  void
  IntegrationInfoBox<dim,sdim>::post_cell(const DoFInfoBox<dim, DOFINFO>&)
  {}


  template <int dim, int sdim>
  template <class DOFINFO>
  void
  IntegrationInfoBox<dim,sdim>::post_faces(const DoFInfoBox<dim, DOFINFO>&)
  {}


  template <int dim, int sdim>
  inline
  void
  IntegrationInfoBox<dim,sdim>::initialize_gauss_quadrature(
    unsigned int cp,
    unsigned int bp,
    unsigned int fp)
  {
    cell_quadrature = QGauss<dim>(cp);
    boundary_quadrature = QGauss<dim-1>(bp);
    face_quadrature = QGauss<dim-1>(fp);
  }


  template <>
  inline
  void
  IntegrationInfoBox<1,1>::
  initialize_gauss_quadrature(const unsigned int cp,
			      const unsigned int,
			      const unsigned int)
  {
    cell_quadrature = QGauss<1>(cp);
  }


  template <>
  inline
  void
  IntegrationInfoBox<1,2>::
  initialize_gauss_quadrature(const unsigned int cp,
			      const unsigned int,
			      const unsigned int)
  {
    cell_quadrature = QGauss<1>(cp);
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
