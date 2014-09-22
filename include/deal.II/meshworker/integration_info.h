// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2014 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#ifndef __deal2__mesh_worker_integration_info_h
#define __deal2__mesh_worker_integration_info_h

#include <deal.II/base/config.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/std_cxx11/shared_ptr.h>
#include <deal.II/dofs/block_info.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/meshworker/local_results.h>
#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/vector_selector.h>

DEAL_II_NAMESPACE_OPEN

namespace MeshWorker
{
  /**
   * Class for objects handed to local integration functions.
   *
   * Objects of this class contain one or more objects of type FEValues,
   * FEFaceValues or FESubfaceValues to be used in local
   * integration. They are stored in an array of pointers to the base
   * classes FEValuesBase. The template parameter VECTOR allows the
   * use of different data types for the global system.
   *
   * Additionally, this function contains space to store the values of
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
    std::vector<std_cxx11::shared_ptr<FEValuesBase<dim, spacedim> > > fevalv;
  public:
    static const unsigned int dimension = dim;
    static const unsigned int space_dimension = spacedim;

    /**
     * Constructor.
     */
    IntegrationInfo();

    /**
     * Copy constructor, creating a clone to be used by
     * WorksTream::run().
     */
    IntegrationInfo(const IntegrationInfo<dim, spacedim> &other);

    /**
     * Build all internal structures, in particular the FEValuesBase
     * objects and allocate space for data vectors.
     *
     * @param el is the finite element of the DoFHandler.
     *
     * @param mapping is the Mapping object used to map the mesh
     * cells.
     *
     * @param quadrature is a Quadrature formula used in the
     * constructor of the FEVALUES objects.
     *
     * @param flags are the UpdateFlags used in the constructor of the
     * FEVALUES objects.
     *
     * @param local_block_info is an optional parameter for systems of
     * PDE. If it is provided with reasonable data, then the degrees
     * of freedom on the cells will be re-ordered to reflect the block
     * structure of the system.
     */
    template <class FEVALUES>
    void initialize(const FiniteElement<dim,spacedim> &el,
                    const Mapping<dim,spacedim> &mapping,
                    const Quadrature<FEVALUES::integral_dimension> &quadrature,
                    const UpdateFlags flags,
                    const BlockInfo *local_block_info = 0);

    /**
     * Initialize the data vector and cache the selector.
     */
    void initialize_data(const std_cxx11::shared_ptr<VectorDataBase<dim,spacedim> > &data);

    /**
     * Delete the data created by initialize().
     */
    void clear();

    /**
     * Return a reference to the FiniteElement that was used to
     * initialize this object.
     */
    const FiniteElement<dim, spacedim> &finite_element() const;

    /// This is true if we are assembling for multigrid
    bool multigrid;
    /// Access to finite element
    /**
     * This is the access function being used, if initialize() for a
     * single element was used (without the BlockInfo argument). It
     * throws an exception, if applied to a vector of elements.
     */
    const FEValuesBase<dim, spacedim> &fe_values () const;

    /// Access to finite elements
    /**
     * This access function must be used if the initalize() for a
     * group of elements was used (with a valid BlockInfo object).
     */
    const FEValuesBase<dim, spacedim> &fe_values (const unsigned int i) const;

    /**
     * The vector containing the values of finite element functions in
     * the quadrature points.
     *
     * There is one vector per selected finite element function,
     * containing one vector for each component, containing vectors
     * with values for each quadrature point.
     */
    std::vector<std::vector<std::vector<double> > > values;

    /**
     * The vector containing the derivatives of finite element
     * functions in the quadrature points.
     *
     * There is one vector per selected finite element function,
     * containing one vector for each component, containing vectors
     * with values for each quadrature point.
     */
    std::vector<std::vector<std::vector<Tensor<1,dim> > > > gradients;

    /**
     * The vector containing the second derivatives of finite element
     * functions in the quadrature points.
     *
     * There is one vector per selected finite element function,
     * containing one vector for each component, containing vectors
     * with values for each quadrature point.
     */
    std::vector<std::vector<std::vector<Tensor<2,dim> > > > hessians;

    /**
     * Reinitialize internal data structures for use on a cell.
     */
    template <typename number>
    void reinit(const DoFInfo<dim, spacedim, number> &i);

    /**
     * Use the finite element functions in #global_data and fill the
     * vectors #values, #gradients and #hessians.
     */
    template<typename number>
    void fill_local_data(const DoFInfo<dim, spacedim, number> &info, bool split_fevalues);

    /**
     * The global data vector used to compute function values in
     * quadrature points.
     */
    std_cxx11::shared_ptr<VectorDataBase<dim, spacedim> > global_data;

    /**
     * The memory used by this object.
     */
    std::size_t memory_consumption () const;

  private:
    /**
     * The pointer to the (system) element used for initialization.
     */
    SmartPointer<const FiniteElement<dim, spacedim>, IntegrationInfo<dim, spacedim> > fe_pointer;

    /**
     * Use the finite element functions in #global_data and fill the
     * vectors #values, #gradients and #hessians with values according
     * to the selector.
     */
    template <typename TYPE>
    void fill_local_data(
      std::vector<std::vector<std::vector<TYPE> > > &data,
      VectorSelector &selector,
      bool split_fevalues) const;
    /**
     * Cache the number of components of the system element.
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
   *
   * <li> It makes educated guesses on quadrature rules and update
   * flags, so that minimal code has to be written when default
   * parameters are sufficient.
   * </ol>
   *
   * In order to allow for sufficient generality, a few steps have to be
   * undertaken to use this class.
   *
   * First, you should consider if you need values from any vectors in a
   * AnyData object. If so, fill the VectorSelector objects
   * #cell_selector, #boundary_selector and #face_selector with their names
   * and the data type (value, gradient, Hessian) to be extracted.
   *
   * Afterwards, you will need to consider UpdateFlags for FEValues
   * objects. A good start is initialize_update_flags(), which looks at
   * the selectors filled before and adds all the flags needed to get
   * the selection. Additional flags can be set with add_update_flags().
   *
   * Finally, we need to choose quadrature formulas. In the simplest
   * case, you might be happy with the default settings, which are
   * <i>n</i>-point Gauss formulas. If only derivatives of the shape
   * functions are used (#update_values is not set) <i>n</i> equals the
   * highest polynomial degree in the FiniteElement, if #update_values
   * is set, <i>n</i> is one higher than this degree.  If you choose to
   * use Gauss formulas of other size, use initialize_gauss_quadrature()
   * with appropriate values. Otherwise, you can fill the variables
   * #cell_quadrature, #boundary_quadrature and #face_quadrature
   * directly.
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

    /**
     * The type of the info object
     * for cells.
     */
    typedef IntegrationInfo<dim, spacedim> CellInfo;

    /**
     * Default constructor.
     */
    IntegrationInfoBox ();

    /**
     * Initialize the IntegrationInfo objects contained.
     *
     * Before doing so, add update flags necessary to produce the data
     * needed and also set uninitialized quadrature rules to Gauss
     * formulas, which integrate polynomial bilinear forms exactly.
     */
    void initialize(const FiniteElement<dim, spacedim> &el,
                    const Mapping<dim, spacedim> &mapping,
                    const BlockInfo *block_info = 0);

    /**
     * Initialize the IntegrationInfo objects contained.
     *
     * Before doing so, add update flags necessary to produce the data
     * needed and also set uninitialized quadrature rules to Gauss
     * formulas, which integrate polynomial bilinear forms exactly.
     */
    template <typename VECTOR>
    void initialize(const FiniteElement<dim, spacedim> &el,
                    const Mapping<dim, spacedim> &mapping,
                    const AnyData &data,
                    const VECTOR &dummy,
                    const BlockInfo *block_info = 0);
    /**
     * @deprecated Use AnyData instead of NamedData.
     */
    template <typename VECTOR>
    void initialize(const FiniteElement<dim, spacedim> &el,
                    const Mapping<dim, spacedim> &mapping,
                    const NamedData<VECTOR *> &data,
                    const BlockInfo *block_info = 0);

    /**
     * @deprecated Use AnyData instead of NamedData.
     */
    template <typename VECTOR>
    void initialize(const FiniteElement<dim, spacedim> &el,
                    const Mapping<dim, spacedim> &mapping,
                    const NamedData<MGLevelObject<VECTOR>*> &data,
                    const BlockInfo *block_info = 0);
    /**
     * @name FEValues setup
     */
    /* @{ */

    /**
    * Call this function before initialize() in order to guess the
    * update flags needed, based on the data selected.
    *
    * When computing face fluxes, we normally can use the geometry
    * (integration weights and normal vectors) from the original cell
    * and thus can avoid updating these values on the neighboring
    * cell. Set <tt>neighbor_geometry</tt> to true in order to
    * initialize these values as well.
    */
    void initialize_update_flags(bool neighbor_geometry = false);

    /**
     * Add FEValues UpdateFlags for integration on all objects (cells,
     * boundary faces and all interior faces).
     */
    void add_update_flags_all (const UpdateFlags flags);

    /**
     * Add FEValues UpdateFlags for integration on cells.
     */
    void add_update_flags_cell(const UpdateFlags flags);

    /**
     * Add FEValues UpdateFlags for integration on boundary faces.
     */
    void add_update_flags_boundary(const UpdateFlags flags);

    /**
     * Add FEValues UpdateFlags for integration on interior faces.
     */
    void add_update_flags_face(const UpdateFlags flags);

    /**
     * Add additional update flags to the ones already set in this
     * program. The four boolean flags indicate whether the additional
     * flags should be set for cell, boundary, interelement face for
     * the cell itself or neighbor cell, or any combination thereof.
     */
    void add_update_flags(const UpdateFlags flags,
                          const bool cell = true,
                          const bool boundary = true,
                          const bool face = true,
                          const bool neighbor = true);

    /**
     * Assign n-point Gauss quadratures to each of the quadrature
     * rules. Here, a size of zero points means that no loop over
     * these grid entities should be performed.
     *
     * If the parameter <tt>force</tt> is true, then all quadrature
     * sets are filled with new quadrature rules. If it is false, then
     * only empty rules are changed.
     */
    void initialize_gauss_quadrature(unsigned int n_cell_points,
                                     unsigned int n_boundary_points,
                                     unsigned int n_face_points,
                                     const bool force = true);

    /**
     * The memory used by this object.
     */
    std::size_t memory_consumption () const;

    /**
     * The set of update flags for boundary cell integration.
     *
     * Defaults to #update_JxW_values.
     */
    UpdateFlags cell_flags;
    /**
     * The set of update flags for boundary face integration.
     *
     * Defaults to #update_JxW_values and #update_normal_vectors.
     */
    UpdateFlags boundary_flags;

    /**
     * The set of update flags for interior face integration.
     *
     * Defaults to #update_JxW_values and #update_normal_vectors.
     */
    UpdateFlags face_flags;

    /**
     * The set of update flags for interior face integration.
     *
     * Defaults to #update_default, since quadrature weights are taken
     * from the other cell.
     */
    UpdateFlags neighbor_flags;

    /**
     * The quadrature rule used on cells.
     */
    Quadrature<dim> cell_quadrature;

    /**
     * The quadrature rule used on boundary faces.
     */
    Quadrature<dim-1> boundary_quadrature;

    /**
     * The quadrature rule used on interior faces.
     */
    Quadrature<dim-1> face_quadrature;
    /* @} */

    /**
     * @name Data vectors
     */
    /* @{ */

    /**
     * Initialize the VectorSelector objects #cell_selector,
     * #boundary_selector and #face_selector in order to save
     * computational effort. If no selectors are used, then values for
     * all named vectors in DoFInfo::global_data will be computed in
     * all quadrature points.
     *
     * This function will also add UpdateFlags to the flags stored in
     * this class.
     */
    /**
     * Select the vectors from DoFInfo::global_data that should be
     * computed in the quadrature points on cells.
     */
    MeshWorker::VectorSelector cell_selector;

    /**
     * Select the vectors from DoFInfo::global_data that should be
     * computed in the quadrature points on boundary faces.
     */
    MeshWorker::VectorSelector boundary_selector;

    /**
     * Select the vectors from DoFInfo::global_data that should be
     * computed in the quadrature points on interior faces.
     */
    MeshWorker::VectorSelector face_selector;

    std_cxx11::shared_ptr<MeshWorker::VectorDataBase<dim, spacedim> > cell_data;
    std_cxx11::shared_ptr<MeshWorker::VectorDataBase<dim, spacedim> > boundary_data;
    std_cxx11::shared_ptr<MeshWorker::VectorDataBase<dim, spacedim> > face_data;
    /* @} */

    /**
     * @name Interface for MeshWorker::loop()
     */
    /* @{ */
    /**
     * A callback function which is called in the loop over all cells,
     * after the action on a cell has been performed and before the
     * faces are dealt with.
     *
     * In order for this function to have this effect, at least either
     * of the arguments <tt>boundary_worker</tt> or
     * <tt>face_worker</tt> arguments of loop() should be
     * nonzero. Additionally, <tt>cells_first</tt> should be true. If
     * <tt>cells_first</tt> is false, this function is called before
     * any action on a cell is taken.
     *
     * And empty function in this class, but can be replaced in other
     * classes given to loop() instead.
     *
     * See loop() and cell_action() for more details of how this
     * function can be used.
     */
    template <class DOFINFO>
    void post_cell(const DoFInfoBox<dim, DOFINFO> &);

    /**
     * A callback function which is called in the loop over all cells,
     * after the action on the faces of a cell has been performed and
     * before the cell itself is dealt with (assumes
     * <tt>cells_first</tt> is false).
     *
     * In order for this function to have a reasonable effect, at
     * least either of the arguments <tt>boundary_worker</tt> or
     * <tt>face_worker</tt> arguments of loop() should be
     * nonzero. Additionally, <tt>cells_first</tt> should be false.
     *
     * And empty function in this class, but can be replaced in other
     * classes given to loop() instead.
     *
     * See loop() and cell_action() for more details of how this
     * function can be used.
     */
    template <class DOFINFO>
    void post_faces(const DoFInfoBox<dim, DOFINFO> &);

    /**
     * The info object for a cell.
     */
    CellInfo cell;
    /**
     * The info object for a boundary face.
     */
    CellInfo boundary;
    /**
     * The info object for a regular interior face, seen from the
     * first cell.
     */
    CellInfo face;
    /**
     * The info object for the refined side of an interior face seen
     * from the first cell.
     */
    CellInfo subface;
    /**
     * The info object for an interior face, seen from the other cell.
     */
    CellInfo neighbor;

    /* @} */
  };


//----------------------------------------------------------------------//

  template<int dim, int sdim>
  inline
  IntegrationInfo<dim,sdim>::IntegrationInfo()
    :
    fevalv(0),
    multigrid(false),
    global_data(std_cxx11::shared_ptr<VectorDataBase<dim, sdim> >(new VectorDataBase<dim, sdim>))
  {}


  template<int dim, int sdim>
  inline
  IntegrationInfo<dim,sdim>::IntegrationInfo(const IntegrationInfo<dim,sdim> &other)
    :
    multigrid(other.multigrid),
    values(other.values),
    gradients(other.gradients),
    hessians(other.hessians),
    global_data(other.global_data),
    fe_pointer(other.fe_pointer),
    n_components(other.n_components)
  {
    fevalv.resize(other.fevalv.size());
    for (unsigned int i=0; i<other.fevalv.size(); ++i)
      {
        const FEValuesBase<dim,sdim> &p = *other.fevalv[i];
        const FEValues<dim,sdim> *pc = dynamic_cast<const FEValues<dim,sdim>*>(&p);
        const FEFaceValues<dim,sdim> *pf = dynamic_cast<const FEFaceValues<dim,sdim>*>(&p);
        const FESubfaceValues<dim,sdim> *ps = dynamic_cast<const FESubfaceValues<dim,sdim>*>(&p);

        if (pc != 0)
          fevalv[i] = std_cxx11::shared_ptr<FEValuesBase<dim,sdim> > (
                        new FEValues<dim,sdim> (pc->get_mapping(), pc->get_fe(),
                                                pc->get_quadrature(), pc->get_update_flags()));
        else if (pf != 0)
          fevalv[i] = std_cxx11::shared_ptr<FEValuesBase<dim,sdim> > (
                        new FEFaceValues<dim,sdim> (pf->get_mapping(), pf->get_fe(), pf->get_quadrature(), pf->get_update_flags()));
        else if (ps != 0)
          fevalv[i] = std_cxx11::shared_ptr<FEValuesBase<dim,sdim> > (
                        new FESubfaceValues<dim,sdim> (ps->get_mapping(), ps->get_fe(), ps->get_quadrature(), ps->get_update_flags()));
        else
          Assert(false, ExcInternalError());
      }
  }



  template<int dim, int sdim>
  template <class FEVALUES>
  inline void
  IntegrationInfo<dim,sdim>::initialize(
    const FiniteElement<dim,sdim> &el,
    const Mapping<dim,sdim> &mapping,
    const Quadrature<FEVALUES::integral_dimension> &quadrature,
    const UpdateFlags flags,
    const BlockInfo *block_info)
  {
    fe_pointer = &el;
    if (block_info == 0 || block_info->local().size() == 0)
      {
        fevalv.resize(1);
        fevalv[0] = std_cxx11::shared_ptr<FEValuesBase<dim,sdim> > (
                      new FEVALUES (mapping, el, quadrature, flags));
      }
    else
      {
        fevalv.resize(el.n_base_elements());
        for (unsigned int i=0; i<fevalv.size(); ++i)
          {
            fevalv[i] = std_cxx11::shared_ptr<FEValuesBase<dim,sdim> > (
                          new FEVALUES (mapping, el.base_element(i), quadrature, flags));
          }
      }
    n_components = el.n_components();
  }


  template <int dim, int spacedim>
  inline const FiniteElement<dim, spacedim> &
  IntegrationInfo<dim,spacedim>::finite_element() const
  {
    Assert (fe_pointer !=0, ExcNotInitialized());
    return *fe_pointer;
  }

  template <int dim, int spacedim>
  inline const FEValuesBase<dim, spacedim> &
  IntegrationInfo<dim,spacedim>::fe_values() const
  {
    AssertDimension(fevalv.size(), 1);
    return *fevalv[0];
  }


  template <int dim, int spacedim>
  inline const FEValuesBase<dim, spacedim> &
  IntegrationInfo<dim,spacedim>::fe_values(unsigned int i) const
  {
    Assert (i<fevalv.size(), ExcIndexRange(i,0,fevalv.size()));
    return *fevalv[i];
  }


  template <int dim, int spacedim>
  template <typename number>
  inline void
  IntegrationInfo<dim,spacedim>::reinit(const DoFInfo<dim, spacedim, number> &info)
  {
    for (unsigned int i=0; i<fevalv.size(); ++i)
      {
        FEValuesBase<dim, spacedim> &febase = *fevalv[i];
        if (info.sub_number != deal_II_numbers::invalid_unsigned_int)
          {
            // This is a subface
            FESubfaceValues<dim> &fe = dynamic_cast<FESubfaceValues<dim>&> (febase);
            fe.reinit(info.cell, info.face_number, info.sub_number);
          }
        else if (info.face_number != deal_II_numbers::invalid_unsigned_int)
          {
            // This is a face
            FEFaceValues<dim> &fe = dynamic_cast<FEFaceValues<dim>&> (febase);
            fe.reinit(info.cell, info.face_number);
          }
        else
          {
            // This is a cell
            FEValues<dim,spacedim> &fe = dynamic_cast<FEValues<dim,spacedim>&> (febase);
            fe.reinit(info.cell);
          }
      }

    const bool split_fevalues = info.block_info != 0;
    if (!global_data->empty())
      fill_local_data(info, split_fevalues);
  }




//----------------------------------------------------------------------//

  template <int dim, int sdim>
  inline
  void
  IntegrationInfoBox<dim,sdim>::initialize_gauss_quadrature(
    unsigned int cp,
    unsigned int bp,
    unsigned int fp,
    bool force)
  {
    if (force || cell_quadrature.size() == 0)
      cell_quadrature = QGauss<dim>(cp);
    if (force || boundary_quadrature.size() == 0)
      boundary_quadrature = QGauss<dim-1>(bp);
    if (force || face_quadrature.size() == 0)
      face_quadrature = QGauss<dim-1>(fp);
  }


  template <int dim, int sdim>
  inline
  void
  IntegrationInfoBox<dim,sdim>::add_update_flags_all (const UpdateFlags flags)
  {
    add_update_flags(flags, true, true, true, true);
  }


  template <int dim, int sdim>
  inline
  void
  IntegrationInfoBox<dim,sdim>::add_update_flags_cell (const UpdateFlags flags)
  {
    add_update_flags(flags, true, false, false, false);
  }


  template <int dim, int sdim>
  inline
  void
  IntegrationInfoBox<dim,sdim>::add_update_flags_boundary (const UpdateFlags flags)
  {
    add_update_flags(flags, false, true, false, false);
  }


  template <int dim, int sdim>
  inline
  void
  IntegrationInfoBox<dim,sdim>::add_update_flags_face (const UpdateFlags flags)
  {
    add_update_flags(flags, false, false, true, true);
  }


  template <int dim, int sdim>
  inline
  void
  IntegrationInfoBox<dim,sdim>::initialize(
    const FiniteElement<dim,sdim> &el,
    const Mapping<dim,sdim> &mapping,
    const BlockInfo *block_info)
  {
    initialize_update_flags();
    initialize_gauss_quadrature(
      (cell_flags & update_values) ? (el.tensor_degree()+1) : el.tensor_degree(),
      (boundary_flags & update_values) ? (el.tensor_degree()+1) : el.tensor_degree(),
      (face_flags & update_values) ? (el.tensor_degree()+1) : el.tensor_degree(), false);

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


  template <int dim, int sdim>
  template <typename VECTOR>
  void
  IntegrationInfoBox<dim,sdim>::initialize(
    const FiniteElement<dim,sdim> &el,
    const Mapping<dim,sdim> &mapping,
    const AnyData &data,
    const VECTOR &dummy,
    const BlockInfo *block_info)
  {
    initialize(el, mapping, block_info);
    std_cxx11::shared_ptr<VectorData<VECTOR, dim, sdim> > p;
    VectorDataBase<dim,sdim> *pp;

    p = std_cxx11::shared_ptr<VectorData<VECTOR, dim, sdim> >(new VectorData<VECTOR, dim, sdim> (cell_selector));
    // Public member function of parent class was not found without
    // explicit cast
    pp = &*p;
    pp->initialize(data);
    cell_data = p;
    cell.initialize_data(p);

    p = std_cxx11::shared_ptr<VectorData<VECTOR, dim, sdim> >(new VectorData<VECTOR, dim, sdim> (boundary_selector));
    pp = &*p;
    pp->initialize(data);
    boundary_data = p;
    boundary.initialize_data(p);

    p = std_cxx11::shared_ptr<VectorData<VECTOR, dim, sdim> >(new VectorData<VECTOR, dim, sdim> (face_selector));
    pp = &*p;
    pp->initialize(data);
    face_data = p;
    face.initialize_data(p);
    subface.initialize_data(p);
    neighbor.initialize_data(p);
  }


  template <int dim, int sdim>
  template <typename VECTOR>
  void
  IntegrationInfoBox<dim,sdim>::initialize(
    const FiniteElement<dim,sdim> &el,
    const Mapping<dim,sdim> &mapping,
    const NamedData<VECTOR *> &data,
    const BlockInfo *block_info)
  {
    initialize(el, mapping, block_info);
    std_cxx11::shared_ptr<VectorData<VECTOR, dim, sdim> > p;

    p = std_cxx11::shared_ptr<VectorData<VECTOR, dim, sdim> >(new VectorData<VECTOR, dim, sdim> (cell_selector));
    p->initialize(data);
    cell_data = p;
    cell.initialize_data(p);

    p = std_cxx11::shared_ptr<VectorData<VECTOR, dim, sdim> >(new VectorData<VECTOR, dim, sdim> (boundary_selector));
    p->initialize(data);
    boundary_data = p;
    boundary.initialize_data(p);

    p = std_cxx11::shared_ptr<VectorData<VECTOR, dim, sdim> >(new VectorData<VECTOR, dim, sdim> (face_selector));
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
    const FiniteElement<dim,sdim> &el,
    const Mapping<dim,sdim> &mapping,
    const NamedData<MGLevelObject<VECTOR>*> &data,
    const BlockInfo *block_info)
  {
    initialize(el, mapping, block_info);
    std_cxx11::shared_ptr<MGVectorData<VECTOR, dim, sdim> > p;

    p = std_cxx11::shared_ptr<MGVectorData<VECTOR, dim, sdim> >(new MGVectorData<VECTOR, dim, sdim> (cell_selector));
    p->initialize(data);
    cell_data = p;
    cell.initialize_data(p);

    p = std_cxx11::shared_ptr<MGVectorData<VECTOR, dim, sdim> >(new MGVectorData<VECTOR, dim, sdim> (boundary_selector));
    p->initialize(data);
    boundary_data = p;
    boundary.initialize_data(p);

    p = std_cxx11::shared_ptr<MGVectorData<VECTOR, dim, sdim> >(new MGVectorData<VECTOR, dim, sdim> (face_selector));
    p->initialize(data);
    face_data = p;
    face.initialize_data(p);
    subface.initialize_data(p);
    neighbor.initialize_data(p);
  }


  template <int dim, int sdim>
  template <class DOFINFO>
  void
  IntegrationInfoBox<dim,sdim>::post_cell(const DoFInfoBox<dim, DOFINFO> &)
  {}


  template <int dim, int sdim>
  template <class DOFINFO>
  void
  IntegrationInfoBox<dim,sdim>::post_faces(const DoFInfoBox<dim, DOFINFO> &)
  {}


}

DEAL_II_NAMESPACE_CLOSE

#endif
