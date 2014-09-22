// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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

#ifndef __deal2__mesh_worker_vector_selector_h
#define __deal2__mesh_worker_vector_selector_h

#include <deal.II/base/named_data.h>
#include <deal.II/algorithms/any_data.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/mg_level_object.h>

DEAL_II_NAMESPACE_OPEN

template<int,int> class FEValuesBase;

namespace MeshWorker
{

  /**
   * A class that selects vectors from a list of named vectors.
   *
   * Since the number of vectors in an AnyData object may grow with every
   * nesting of applications or loops, it is important to be able to
   * select those, which are actually used in computing residuals etc.
   * This class organizes the selection.
   *
   * It is used for instance in IntegrationWorker to
   * determine which values, derivatives or second derivatives are
   * actually computed.
   *
   * @ingroup MeshWorker
   * @author Guido Kanschat 2009
   */
  class VectorSelector :
    public Subscriptor
  {
  public:
    /**
     * Add a vector to the selection of finite element functions. The
     * arguments are the name of the vector and indicators, which
     * information is to be extracted from the vector. The name refers
     * to an entry in a AnyData object, which will be identified by
     * initialize().  The three bool parameters indicate, whether
     * values, greadients and Hessians of the finite element function
     * are to be computed on each cell or face.
     */
    void add(const std::string &name,
             const bool values = true,
             const bool gradients = false,
             const bool hessians = false);

    /**
     * Does the same as the function above but it is possible to
     * select a block of the global vector.
     */
//      void add(const std::string& name,
//               const unsigned int selected_block,
//             bool values = true,
//             bool gradients = false,
//             bool hessians = false);

    /**
     * Initialize the selection field with a data vector. While add()
     * only enters the names of vectors, which will be used in the
     * integration loop over cells and faces, this function links the
     * names to actual vectos in a AnyData object.
     *
     * @note This function caches the index associated with a
     * name. Therefore, it must be called every time after the AnyData
     * object has changed.
     */
    void initialize(const AnyData &);

    /**
     * @deprecated Use AnyData instead of NamedData.
     */
    template <class DATA>
    void initialize(const NamedData<DATA> &);

    /**
     * Check whether any vector is selected.
     */
    bool empty () const;

    /**
     * Returns true if values are selected for any vector.
     */
    bool has_values () const;

    /**
     * Returns true if gradients are selected for any vector.
     */
    bool has_gradients () const;

    /**
     * Returns true if hessians are selected for any vector.
     */
    bool has_hessians () const;

    /**
     * Number of vectors for values
     */
    unsigned int n_values () const;

    /**
     * Number of vectors for gradients
     */
    unsigned int n_gradients () const;

    /**
     * Number of vectors for Hessians
     */
    unsigned int n_hessians () const;

    /**
     * The vector index for the ith value
     */
    unsigned int value_index (const unsigned int i) const;

    /**
     * The vector index for the ith gradient
     */
    unsigned int gradient_index (const unsigned int i) const;

    /**
     * The vector index for the ith Hessian
     */
    unsigned int hessian_index (const unsigned int i) const;

    /**
     * Print the contents of the
     * selection to the stream.
     */
    template <class STREAM, typename DATA>
    void print (STREAM &s, const AnyData &v) const;

    /**
     * @deprecated Use AnyData instead of NamedData.
     */
    template <class STREAM, typename DATA>
    void print (STREAM &s, const NamedData<DATA> &v) const;

    /**
     * Print the number of selections to the stream.
     */
    template <class STREAM>
    void print (STREAM &s) const;

    /**
     * The memory used by this object.
     */
    std::size_t memory_consumption () const;

  protected:
    /**
     * Selection of the vectors
     * used to compute values.
     */
    NamedSelection value_selection;

    /**
     * Selection of the vectors used to compute gradients.
     */
    NamedSelection gradient_selection;

    /**
     * Selection of the vectors used to compute hessians.
     */
    NamedSelection hessian_selection;
  };

  /**
   * Based on VectorSelector, this is the class used by IntegrationInfo
   * to compute values of source vectors in quadrature points.
   *
   * @ingroup MeshWorker
   * @author Guido Kanschat, 2009
   */
  template <int dim, int spacedim = dim>
  class VectorDataBase :
    public VectorSelector
  {
  public:
    /**
     * Constructor
     */
    VectorDataBase();

    /**
     * Constructor from a base class object
     */
    VectorDataBase(const VectorSelector &);

    /**
     * Initialize with a AnyData object and cache the indices in the
     * VectorSelector base class.
     *
     * @note Make sure the VectorSelector base class was filled with
     * reasonable data before calling this function.
     */
    void initialize(const AnyData &);

    /**
     * Virtual, but empty destructor.
     */
    virtual ~VectorDataBase();
    /**
     * The only function added to VectorSelector is an abstract
     * virtual function implemented in the derived class template and
     * called by IntegrationInfo.
     *
     * Depending on the selections made in our base class, this fills
     * the first three arguments with the local data of the finite
     * element functions. It is usually called either for the whole
     * FESystem, or for each base element separately.
     *
     * @param values is the vector filled with the values of the
     * finite element function in the quadrature points.
     *
     * @param gradients is the vector filled with the derivatives of
     * the finite element function in the quadrature points.
     *
     * @param hessians is the vector filled with the second
     * derivatives of the finite element function in the quadrature
     * points.
     *
     * @param fe is the FEValuesBase object which is used to compute
     * the function values. Its UpdateFlags must have been set
     * appropriately.
     *
     * @param index is the local index vector. If @p fe refers to base
     * elements of the system, this vector should be sorted by block
     * and the arguments @p start and @p size below specify the subset
     * of @p indices used.
     *
     * @param component is the first index in @p values, @p gradients
     * and @p hessians entered in this function.
     *
     * @param n_comp is the number of components to be filled.
     *
     * @param start is the first index of this block in @p indices, or
     * zero if no base elements are used.
     *
     * @param size is the number of dofs per cell of the current
     * element or base element.
     */
    virtual void fill(
      std::vector<std::vector<std::vector<double> > > &values,
      std::vector<std::vector<std::vector<Tensor<1,dim> > > > &gradients,
      std::vector<std::vector<std::vector<Tensor<2,dim> > > > &hessians,
      const FEValuesBase<dim,spacedim> &fe,
      const std::vector<types::global_dof_index> &index,
      const unsigned int component,
      const unsigned int n_comp,
      const unsigned int start,
      const unsigned int size) const;

    /**
     * Fill the local data vector from level vectors. Performs exactly
     * what the other fill() does, but uses the cell level to access a
     * single level out of a hierarchy of level vectors, instead of a
     * global data vector on the active cells.
     */
    virtual void mg_fill(
      std::vector<std::vector<std::vector<double> > > &values,
      std::vector<std::vector<std::vector<Tensor<1,dim> > > > &gradients,
      std::vector<std::vector<std::vector<Tensor<2,dim> > > > &hessians,
      const FEValuesBase<dim,spacedim> &fe,
      const unsigned int level,
      const std::vector<types::global_dof_index> &index,
      const unsigned int component,
      const unsigned int n_comp,
      const unsigned int start,
      const unsigned int size) const;

    /**
     * The memory used by this object.
     */
    std::size_t memory_consumption () const;
  protected:
    AnyData data;
  };


  /**
   * Based on VectorSelector, this is the class that implements the
   * function VectorDataBase::fill() for a certain type of vector, using
   * AnyData to identify vectors by name.
   *
   * @ingroup MeshWorker
   * @author Guido Kanschat, 2009
   */
  template <class VECTOR, int dim, int spacedim = dim>
  class VectorData :
    public VectorDataBase<dim, spacedim>
  {
  public:
    /**
     * Constructor.
     */
    VectorData();
    /**
     * Constructor using a
     * prefilled VectorSelector
     */
    VectorData(const VectorSelector &);

    /**
     * @deprecated Use AnyData instead of NamedData.
     */
    void initialize(const NamedData<VECTOR *> &);

    /**
     * Initialize with a single vector and cache the indices in the
     * VectorSelector base class.
     *
     * @note Make sure the VectorSelector base class was filled with
     * reasonable data before calling this function.
     */
    void initialize(const VECTOR *, const std::string &name);

    virtual void fill(
      std::vector<std::vector<std::vector<double> > > &values,
      std::vector<std::vector<std::vector<Tensor<1,dim> > > > &gradients,
      std::vector<std::vector<std::vector<Tensor<2,dim> > > > &hessians,
      const FEValuesBase<dim,spacedim> &fe,
      const std::vector<types::global_dof_index> &index,
      const unsigned int component,
      const unsigned int n_comp,
      const unsigned int start,
      const unsigned int size) const;

    virtual void mg_fill(
      std::vector<std::vector<std::vector<double> > > &values,
      std::vector<std::vector<std::vector<Tensor<1,dim> > > > &gradients,
      std::vector<std::vector<std::vector<Tensor<2,dim> > > > &hessians,
      const FEValuesBase<dim,spacedim> &fe,
      const unsigned int level,
      const std::vector<types::global_dof_index> &index,
      const unsigned int component,
      const unsigned int n_comp,
      const unsigned int start,
      const unsigned int size) const;

    /**
    * The memory used by this object.
    */
    std::size_t memory_consumption () const;
  };


  /**
   * Based on VectorSelector, this is the class that implements the
   * function VectorDataBase::fill() for a certain type of multilevel vectors, using
   * AnyData to identify vectors by name.
   *
   * @ingroup MeshWorker
   * @author Guido Kanschat, 2010
   */
  template <class VECTOR, int dim, int spacedim = dim>
  class MGVectorData :
    public VectorData<VECTOR, dim, spacedim>
  {
  public:
    /**
     * Constructor.
     */
    MGVectorData();
    /**
     * Constructor using a prefilled VectorSelector
     */
    MGVectorData(const VectorSelector &);

    /**
     * @deprecated Use AnyData instead of NamedData.
     */
    void initialize(const NamedData<MGLevelObject<VECTOR>*> &);

    /**
     * Initialize with a single vector and cache the indices in the
     * VectorSelector base class.
     *
     * @note Make sure the VectorSelector base class was filled with
     * reasonable data before calling this function.
     */
    void initialize(const MGLevelObject<VECTOR> *, const std::string &name);


    /**
     * The memory used by this object.
     */
    std::size_t memory_consumption () const;
  };


//----------------------------------------------------------------------//

  inline void
  VectorSelector::add(const std::string &name,
                      const bool values,
                      const bool gradients,
                      const bool hessians)
  {
    if (values)
      value_selection.add(name);
    if (gradients)
      gradient_selection.add(name);
    if (hessians)
      hessian_selection.add(name);
  }


  //inline void
  //VectorSelector::add(const std::string& name,
  //   const unsigned int block,
  //   bool values, bool gradients, bool hessians)
  //{
  //  if (values) value_selection.add(name, block);
  //  if (gradients) gradient_selection.add(name, block);
  //  if (hessians) hessian_selection.add(name, block);
  //}


  inline void
  VectorSelector::initialize(const AnyData &src)
  {
    value_selection.initialize(src);
    gradient_selection.initialize(src);
    hessian_selection.initialize(src);
  }

  template <typename DATA>
  inline void
  VectorSelector::initialize(const NamedData<DATA> &src)
  {
    value_selection.initialize(src);
    gradient_selection.initialize(src);
    hessian_selection.initialize(src);
  }

  inline bool
  VectorSelector::empty() const
  {
    return (value_selection.size() == 0 &&
            gradient_selection.size() == 0 &&
            hessian_selection.size() == 0);
  }


  inline bool
  VectorSelector::has_values() const
  {
    return value_selection.size() != 0;
  }


  inline bool
  VectorSelector::has_gradients() const
  {
    return gradient_selection.size() != 0;
  }


  inline bool
  VectorSelector::has_hessians() const
  {
    return hessian_selection.size() != 0;
  }


  inline unsigned int
  VectorSelector::n_values() const
  {
    return value_selection.size();
  }


  inline unsigned int
  VectorSelector::n_gradients() const
  {
    return gradient_selection.size();
  }


  inline unsigned int
  VectorSelector::n_hessians() const
  {
    return hessian_selection.size();
  }


  inline unsigned int
  VectorSelector::value_index(const unsigned int i) const
  {
    return value_selection(i);
  }


  inline unsigned int
  VectorSelector::gradient_index(const unsigned int i) const
  {
    return gradient_selection(i);
  }


  inline unsigned int
  VectorSelector::hessian_index(const unsigned int i) const
  {
    return hessian_selection(i);
  }


  template <class STREAM>
  inline void
  VectorSelector::print(STREAM &s) const
  {
    s << "values: " << n_values()
      << " gradients: " << n_gradients()
      << " hessians: " << n_hessians()
      << std::endl;
  }


  template <class STREAM, typename DATA>
  inline void
  VectorSelector::print(STREAM &s, const AnyData &v) const
  {
    s << "values:   ";
    for (unsigned int i=0; i<n_values(); ++i)
      s << " '" << v.name(value_selection(i)) << '\'';
    s << std::endl << "gradients:";
    for (unsigned int i=0; i<n_gradients(); ++i)
      s << " '" << v.name(gradient_selection(i)) << '\'';
    s << std::endl << "hessians: ";
    for (unsigned int i=0; i<n_hessians(); ++i)
      s << " '" << v.name(hessian_selection(i)) << '\'';
    s << std::endl;
  }


  template <class STREAM, typename DATA>
  inline void
  VectorSelector::print(STREAM &s, const NamedData<DATA> &v) const
  {
    s << "values:   ";
    for (unsigned int i=0; i<n_values(); ++i)
      s << " '" << v.name(value_selection(i)) << '\'';
    s << std::endl << "gradients:";
    for (unsigned int i=0; i<n_gradients(); ++i)
      s << " '" << v.name(gradient_selection(i)) << '\'';
    s << std::endl << "hessians: ";
    for (unsigned int i=0; i<n_hessians(); ++i)
      s << " '" << v.name(hessian_selection(i)) << '\'';
    s << std::endl;
  }


  inline
  std::size_t
  VectorSelector::memory_consumption () const
  {
    return sizeof(*this);
  }
}


DEAL_II_NAMESPACE_CLOSE

#endif
