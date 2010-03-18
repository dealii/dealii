//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__mesh_worker_vector_selector_h
#define __deal2__mesh_worker_vector_selector_h

#include <base/named_data.h>
#include <base/tensor.h>
#include <base/smartpointer.h>
#include <multigrid/mg_level_object.h>

DEAL_II_NAMESPACE_OPEN

template<int,int> class FEValuesBase;

namespace MeshWorker
{
  
/**
 * A class that selects vectors from a list of named vectors.
 *
 * Since the number of vectors in a NamedData object may grow with every
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
					* Add a vector to the
					* selection of finite element
					* functions. The arguments are
					* the name of the vector and
					* indicators, which
					* information is to be
					* extracted from the
					* vector. The name refers to
					* an entry in a NamedData
					* object, which will be
					* identified by initialize().
					* The three bool parameters
					* indicate, whether values,
					* greadients and Hessians of
					* the finite element function
					* are to be computed on each
					* cell or face.
					*/
      void add(const std::string& name,
	       bool values = true,
	       bool gradients = false,
	       bool hessians = false);

				       /**
					* Initialize the selection
					* field with a data
					* vector. While add() only
					* enters the names of vectors,
					* which will be used in the integration
					* loop over cells and faces,
					* this function links the
					* names to actual vectos in a
					* NamedData object.
					*
					* @note This function caches
					* the index associated with a
					* name. Therefore, it must be
					* called everytime after the
					* NamedData object has changed.
					*/
      template <class DATA>
      void initialize(const NamedData<DATA>&);
      
				       /**
					* Check whether any vector is selected.
					*/
      bool empty () const;

				       /**
					* Returns true if values are
					* selected for any vector.
					*/
      bool has_values () const;
      
				       /**
					* Returns true if gradients are
					* selected for any vector.
					*/
      bool has_gradients () const;
      
				       /**
					* Returns true if hessians are
					* selected for any vector.
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
      unsigned int value_index (unsigned int i) const;
      
				       /**
					* The vector index for the ith gradient
					*/
      unsigned int gradient_index (unsigned int i) const;
      
				       /**
					* The vector index for the ith Hessian
					*/
      unsigned int hessian_index (unsigned int i) const;
      
				       /**
					* Print the contents of the
					* selection to the stream.
					*/
      template <class STREAM, typename DATA>
      void print (STREAM& s, const NamedData<DATA>& v) const;
      
				       /**
					* Print the number of selections to the stream.
					*/
      template <class STREAM>
      void print (STREAM& s) const;
      
    protected:
				       /**
					* Selection of the vectors
					* used to compute values.
					*/
      NamedSelection value_selection;
      
				       /**
					* Selection of the vectors
					* used to compute gradients.
					*/
      NamedSelection gradient_selection;
      
				       /**
					* Selection of the vectors
					* used to compute hessians.
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
					* Constructor from a base
					* class object
					*/
      VectorDataBase(const VectorSelector&);
      
				       /**
					* Virtual, but empty
					* destructor.
					*/
      virtual ~VectorDataBase();
				       /**
					* The only function added to
					* VectorSelector is an
					* abstract virtual function
					* implemented in the derived
					* class template and called by
					* IntegrationInfo.
					*
					* Depending on the selections
					* made in our base class, this
					* fills the first three
					* arguments with the local
					* data of the finite element
					* functions. It is usually
					* called either for the whole
					* FESystem, or for each base
					* element separately.
					*
					* @param values is the vector
					* filled with the values of
					* the finite element function
					* in the quadrature points.
					*
					* @param gradients is the vector
					* filled with the derivatives of
					* the finite element function
					* in the quadrature points.
					*
					* @param hessians is the
					* vector filled with the
					* second derivatives of the
					* finite element function in
					* the quadrature points.
					*
					* @param fe is the
					* FEValuesBase object which is
					* used to compute the function
					* values. Its UpdateFlags must
					* have been set appropriately.
					*
					* @param index is the local
					* index vector. If @p fe
					* refers to base elements of
					* the system, this vector
					* should be sorted by block
					* and the arguments @p start
					* and @p size below specify
					* the subset of @p indices
					* used.
					*
					* @param component is the first
					* index in @p values, @p
					* gradients and @p hessians
					* entered in this function.
					*
					* @param n_comp is the number
					* of components to be filled.
					*
					* @param start is the first
					* index of this block in @p
					* indices, or zero if
					* no base elements are used.
					*
					* @param size is the number of
					* dofs per cell of the current
					* element or base element.
					*/
      virtual void fill(
	std::vector<std::vector<std::vector<double> > >& values,
	std::vector<std::vector<std::vector<Tensor<1,dim> > > >& gradients,
	std::vector<std::vector<std::vector<Tensor<2,dim> > > >& hessians,
	const FEValuesBase<dim,spacedim>& fe,
	const std::vector<unsigned int>& index,
	unsigned int component,
	unsigned int n_comp,
	unsigned int start,
	unsigned int size) const;

				       /**
					* Fill the local data vector
					* from level vectors. Performs
					* exactly what the other
					* fill() does, but uses the
					* cell level to access a
					* single level out of a
					* hierarchy of level vectors,
					* instead of a global data
					* vector on the active cells.
					*/
       virtual void mg_fill(
	std::vector<std::vector<std::vector<double> > >& values,
	std::vector<std::vector<std::vector<Tensor<1,dim> > > >& gradients,
	std::vector<std::vector<std::vector<Tensor<2,dim> > > >& hessians,
	const FEValuesBase<dim,spacedim>& fe,
	unsigned int level,
	const std::vector<unsigned int>& index,
	unsigned int component,
	unsigned int n_comp,
	unsigned int start,
	unsigned int size) const;     
  };

  
/**
 * Based on VectorSelector, this is the class that implements the
 * function VectorDataBase::fill() for a certain type of vector, using
 * NamedData to identify vectors by name.
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
      VectorData(const VectorSelector&);

				       /**
					* Initialize with a NamedData
					* object and cache the indices
					* in the VectorSelector base
					* class.
					*
					* @note Make sure the
					* VectorSelector base class
					* was filled with reasonable
					* data before calling this
					* function.
					*/
      void initialize(const NamedData<VECTOR*>&);

				       /**
					* Initialize with a single
					* vector and cache the indices
					* in the VectorSelector base
					* class.
					*
					* @note Make sure the
					* VectorSelector base class
					* was filled with reasonable
					* data before calling this
					* function.
					*/
      void initialize(const VECTOR*, const std::string& name);
      
      virtual void fill(
	std::vector<std::vector<std::vector<double> > >& values,
	std::vector<std::vector<std::vector<Tensor<1,dim> > > >& gradients,
	std::vector<std::vector<std::vector<Tensor<2,dim> > > >& hessians,
	const FEValuesBase<dim,spacedim>& fe,
	const std::vector<unsigned int>& index,
	unsigned int component,
	unsigned int n_comp,
	unsigned int start,
	unsigned int size) const;
    private:
      NamedData<SmartPointer<const VECTOR,VectorData<VECTOR,dim,spacedim> > > data;
  };

  
/**
 * Based on VectorSelector, this is the class that implements the
 * function VectorDataBase::fill() for a certain type of multilevel vectors, using
 * NamedData to identify vectors by name.
 *
 * @ingroup MeshWorker
 * @author Guido Kanschat, 2010
 */
  template <class VECTOR, int dim, int spacedim = dim>
  class MGVectorData :
      public VectorDataBase<dim, spacedim>
  {
    public:
				       /**
					* Constructor.
					*/
      MGVectorData();
				       /**
					* Constructor using a
					* prefilled VectorSelector
					*/
      MGVectorData(const VectorSelector&);

				       /**
					* Initialize with a NamedData
					* object and cache the indices
					* in the VectorSelector base
					* class.
					*
					* @note Make sure the
					* VectorSelector base class
					* was filled with reasonable
					* data before calling this
					* function.
					*/
      void initialize(const NamedData<MGLevelObject<VECTOR>*>&);

				       /**
					* Initialize with a single
					* vector and cache the indices
					* in the VectorSelector base
					* class.
					*
					* @note Make sure the
					* VectorSelector base class
					* was filled with reasonable
					* data before calling this
					* function.
					*/
      void initialize(const MGLevelObject<VECTOR>*, const std::string& name);
      
      virtual void mg_fill(
	std::vector<std::vector<std::vector<double> > >& values,
	std::vector<std::vector<std::vector<Tensor<1,dim> > > >& gradients,
	std::vector<std::vector<std::vector<Tensor<2,dim> > > >& hessians,
	const FEValuesBase<dim,spacedim>& fe,
	unsigned int level,
	const std::vector<unsigned int>& index,
	unsigned int component,
	unsigned int n_comp,
	unsigned int start,
	unsigned int size) const;
    private:
      NamedData<SmartPointer<const MGLevelObject<VECTOR>,MGVectorData<VECTOR,dim,spacedim> > > data;
  };

  
//----------------------------------------------------------------------//

  inline void
  VectorSelector::add(const std::string& name, bool values, bool gradients, bool hessians)
  {
    if (values) value_selection.add(name);
    if (gradients) gradient_selection.add(name);
    if (hessians) hessian_selection.add(name);  
  }


  template <typename DATA>
  inline void
  VectorSelector::initialize(const NamedData<DATA>& src)
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
  VectorSelector::value_index(unsigned int i) const
  {
    return value_selection(i);
  }
  

  inline unsigned int
  VectorSelector::gradient_index(unsigned int i) const
  {
    return gradient_selection(i);
  }
  

  inline unsigned int
  VectorSelector::hessian_index(unsigned int i) const
  {
    return hessian_selection(i);
  }
  

  template <class STREAM>
  inline void
  VectorSelector::print(STREAM& s) const
  {
    s << "values: " << n_values()
      << " gradients: " << n_gradients()
      << " hessians: " << n_hessians()
      << std::endl;
  }


  template <class STREAM, typename DATA>
  inline void
  VectorSelector::print(STREAM& s, const NamedData<DATA>& v) const
  {
    s << "values:   ";
    for (unsigned int i=0;i<n_values();++i)
      s << " '" << v.name(value_selection(i)) << '\'';
    s << std::endl << "gradients:";
    for (unsigned int i=0;i<n_gradients();++i)
      s << " '" << v.name(gradient_selection(i)) << '\'';
    s << std::endl << "hessians: ";
    for (unsigned int i=0;i<n_hessians();++i)
      s << " '" << v.name(hessian_selection(i)) << '\'';
    s << std::endl;
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
