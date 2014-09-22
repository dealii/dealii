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


#include <deal.II/meshworker/vector_selector.h>
#include <deal.II/base/vector_slice.h>
#include <deal.II/fe/fe_values.h>

DEAL_II_NAMESPACE_OPEN

namespace MeshWorker
{
  template <int dim, int spacedim>
  VectorDataBase<dim, spacedim>::~VectorDataBase()
  {}


  template <int dim, int spacedim>
  VectorDataBase<dim, spacedim>::VectorDataBase(const VectorSelector &v)
    :
    VectorSelector(v)
  {}


  template <int dim, int spacedim>
  VectorDataBase<dim, spacedim>::VectorDataBase()
  {}


  template <int dim, int spacedim>
  void
  VectorDataBase<dim, spacedim>::initialize(const AnyData &d)
  {
    this->data = d;
    VectorSelector::initialize(d);
  }


  template <int dim, int spacedim>
  void
  VectorDataBase<dim, spacedim>::fill(
    std::vector<std::vector<std::vector<double> > > &,
    std::vector<std::vector<std::vector<Tensor<1,dim> > > > &,
    std::vector<std::vector<std::vector<Tensor<2,dim> > > > &,
    const FEValuesBase<dim,spacedim> &,
    const std::vector<types::global_dof_index> &,
    const unsigned int,
    const unsigned int,
    const unsigned int,
    const unsigned int) const
  {
    Assert(false, ExcNotImplemented());
  }


  template <int dim, int spacedim>
  void
  VectorDataBase<dim, spacedim>::mg_fill(
    std::vector<std::vector<std::vector<double> > > &,
    std::vector<std::vector<std::vector<Tensor<1,dim> > > > &,
    std::vector<std::vector<std::vector<Tensor<2,dim> > > > &,
    const FEValuesBase<dim,spacedim> &,
    const unsigned int,
    const std::vector<types::global_dof_index> &,
    const unsigned int,
    const unsigned int,
    const unsigned int,
    const unsigned int) const
  {
    Assert(false, ExcNotImplemented());
  }


//----------------------------------------------------------------------//

  template <class VECTOR, int dim, int spacedim>
  VectorData<VECTOR, dim, spacedim>::VectorData()
  {}


  template <class VECTOR, int dim, int spacedim>
  VectorData<VECTOR, dim, spacedim>::VectorData(const VectorSelector &s)
    :
    VectorDataBase<dim, spacedim>(s)
  {}


  template <class VECTOR, int dim, int spacedim>
  void
  VectorData<VECTOR, dim, spacedim>::initialize(const NamedData<VECTOR *> &d)
  {
    this->data = d;
    VectorSelector::initialize(d);
  }


  template <class VECTOR, int dim, int spacedim>
  void
  VectorData<VECTOR, dim, spacedim>::initialize(const VECTOR *v, const std::string &name)
  {
    SmartPointer<const VECTOR,VectorData<VECTOR, dim, spacedim> > p = v;
    this->data.add(p, name);
    VectorSelector::initialize(this->data);
  }


  template <class VECTOR, int dim, int spacedim>
  void
  VectorData<VECTOR, dim, spacedim>::fill(
    std::vector<std::vector<std::vector<double> > > &values,
    std::vector<std::vector<std::vector<Tensor<1,dim> > > > &gradients,
    std::vector<std::vector<std::vector<Tensor<2,dim> > > > &hessians,
    const FEValuesBase<dim,spacedim> &fe,
    const std::vector<types::global_dof_index> &index,
    const unsigned int component,
    const unsigned int n_comp,
    const unsigned int start,
    const unsigned int size) const
  {
    AssertDimension(values.size(), this->n_values());
    AssertDimension(gradients.size(), this->n_gradients());
    AssertDimension(hessians.size(), this->n_hessians());

    const AnyData &data = this->data;
    for (unsigned int i=0; i<this->n_values(); ++i)
      {
        const VECTOR *src = data.read_ptr<VECTOR>(this->value_index(i));
        VectorSlice<std::vector<std::vector<double> > > dst(values[i], component, n_comp);
        fe.get_function_values(*src, make_slice(index, start, size), dst, true);
      }

    for (unsigned int i=0; i<this->n_gradients(); ++i)
      {
        const VECTOR *src = data.read_ptr<VECTOR>(this->gradient_index(i));
        VectorSlice<std::vector<std::vector<Tensor<1,dim> > > > dst(gradients[i], component, n_comp);
        fe.get_function_gradients(*src, make_slice(index, start, size), dst, true);
      }

    for (unsigned int i=0; i<this->n_hessians(); ++i)
      {
        const VECTOR *src = data.read_ptr<VECTOR>(this->hessian_index(i));
        VectorSlice<std::vector<std::vector<Tensor<2,dim> > > > dst(hessians[i], component, n_comp);
        fe.get_function_hessians(*src, make_slice(index, start, size), dst, true);
      }
  }


  template <class VECTOR, int dim, int spacedim>
  std::size_t
  VectorData<VECTOR, dim, spacedim>::memory_consumption () const
  {
    std::size_t mem = VectorSelector::memory_consumption();
    mem += sizeof (this->data);
    return mem;
  }

//----------------------------------------------------------------------//

  template <class VECTOR, int dim, int spacedim>
  MGVectorData<VECTOR, dim, spacedim>::MGVectorData()
  {}


  template <class VECTOR, int dim, int spacedim>
  MGVectorData<VECTOR, dim, spacedim>::MGVectorData(const VectorSelector &s)
    :
    VectorData<VECTOR, dim, spacedim>(s)
  {}


  template <class VECTOR, int dim, int spacedim>
  void
  MGVectorData<VECTOR, dim, spacedim>::initialize(const NamedData<MGLevelObject<VECTOR>*> &d)
  {
    this->data = d;
    VectorSelector::initialize(d);
  }


  template <class VECTOR, int dim, int spacedim>
  void
  MGVectorData<VECTOR, dim, spacedim>::initialize(const MGLevelObject<VECTOR> *v, const std::string &name)
  {
    SmartPointer<const MGLevelObject<VECTOR>, MGVectorData<VECTOR, dim, spacedim> >
    p = v;
    this->data.add(p, name);
    VectorSelector::initialize(this->data);
  }


  template <class VECTOR, int dim, int spacedim>
  void
  VectorData<VECTOR, dim, spacedim>::mg_fill(
    std::vector<std::vector<std::vector<double> > > &values,
    std::vector<std::vector<std::vector<Tensor<1,dim> > > > &gradients,
    std::vector<std::vector<std::vector<Tensor<2,dim> > > > &hessians,
    const FEValuesBase<dim,spacedim> &fe,
    const unsigned int level,
    const std::vector<types::global_dof_index> &index,
    const unsigned int component,
    const unsigned int n_comp,
    const unsigned int start,
    const unsigned int size) const
  {
    AssertDimension(values.size(), this->n_values());
    AssertDimension(gradients.size(), this->n_gradients());
    AssertDimension(hessians.size(), this->n_hessians());

    const AnyData &data = this->data;
    for (unsigned int i=0; i<this->n_values(); ++i)
      {
        const MGLevelObject<VECTOR> *src = data.read_ptr<MGLevelObject<VECTOR> >(this->value_index(i));
        VectorSlice<std::vector<std::vector<double> > > dst(values[i], component, n_comp);
        fe.get_function_values((*src)[level], make_slice(index, start, size), dst, true);
      }

    for (unsigned int i=0; i<this->n_gradients(); ++i)
      {
        const MGLevelObject<VECTOR> *src = data.read_ptr<MGLevelObject<VECTOR> >(this->value_index(i));
        VectorSlice<std::vector<std::vector<Tensor<1,dim> > > > dst(gradients[i], component, n_comp);
        fe.get_function_gradients((*src)[level], make_slice(index, start, size), dst, true);
      }

    for (unsigned int i=0; i<this->n_hessians(); ++i)
      {
        const MGLevelObject<VECTOR> *src = data.read_ptr<MGLevelObject<VECTOR> >(this->value_index(i));
        VectorSlice<std::vector<std::vector<Tensor<2,dim> > > > dst(hessians[i], component, n_comp);
        fe.get_function_hessians((*src)[level], make_slice(index, start, size), dst, true);
      }
  }
}

DEAL_II_NAMESPACE_CLOSE
