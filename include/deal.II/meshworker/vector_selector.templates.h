// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2016 by the deal.II authors
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

#ifndef dealii__vector_selector_templates_h
#define dealii__vector_selector_templates_h


#include <deal.II/meshworker/vector_selector.h>
#include <deal.II/base/vector_slice.h>
#include <deal.II/fe/fe_values.h>

DEAL_II_NAMESPACE_OPEN

namespace MeshWorker
{
  template <int dim, int spacedim, typename Number>
  VectorDataBase<dim, spacedim, Number>::~VectorDataBase()
  {}


  template <int dim, int spacedim, typename Number>
  VectorDataBase<dim, spacedim, Number>::VectorDataBase(const VectorSelector &v)
    :
    VectorSelector(v)
  {}


  template <int dim, int spacedim, typename Number>
  VectorDataBase<dim, spacedim, Number>::VectorDataBase()
  {}


  template <int dim, int spacedim, typename Number>
  void
  VectorDataBase<dim, spacedim, Number>::initialize(const AnyData &d)
  {
    this->data = d;
    VectorSelector::initialize(d);
  }


  template <int dim, int spacedim, typename Number>
  void
  VectorDataBase<dim, spacedim, Number>::fill(
    std::vector<std::vector<std::vector<Number> > > &,
    std::vector<std::vector<std::vector<Tensor<1,dim,Number> > > > &,
    std::vector<std::vector<std::vector<Tensor<2,dim,Number> > > > &,
    const FEValuesBase<dim,spacedim> &,
    const std::vector<types::global_dof_index> &,
    const unsigned int,
    const unsigned int,
    const unsigned int,
    const unsigned int) const
  {
    Assert(false, ExcNotImplemented());
  }


  template <int dim, int spacedim, typename Number>
  void
  VectorDataBase<dim, spacedim, Number>::mg_fill(
    std::vector<std::vector<std::vector<Number> > > &,
    std::vector<std::vector<std::vector<Tensor<1,dim,Number> > > > &,
    std::vector<std::vector<std::vector<Tensor<2,dim,Number> > > > &,
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

  template <typename VectorType, int dim, int spacedim>
  VectorData<VectorType, dim, spacedim>::VectorData()
  {}


  template <typename VectorType, int dim, int spacedim>
  VectorData<VectorType, dim, spacedim>::VectorData(const VectorSelector &s)
    :
    VectorDataBase<dim, spacedim, typename VectorType::value_type>(s)
  {}


  template <typename VectorType, int dim, int spacedim>
  void
  VectorData<VectorType, dim, spacedim>::initialize(const AnyData &d)
  {
    this->data = d;
    VectorSelector::initialize(d);
  }


  template <typename VectorType, int dim, int spacedim>
  void
  VectorData<VectorType, dim, spacedim>::initialize(const VectorType *v, const std::string &name)
  {
    SmartPointer<const VectorType,VectorData<VectorType, dim, spacedim> > p = v;
    this->data.add(p, name);
    VectorSelector::initialize(this->data);
  }


  template <typename VectorType, int dim, int spacedim>
  void
  VectorData<VectorType, dim, spacedim>::fill(
    std::vector<std::vector<std::vector<typename VectorType::value_type> > > &values,
    std::vector<std::vector<std::vector<Tensor<1,dim,typename VectorType::value_type> > > > &gradients,
    std::vector<std::vector<std::vector<Tensor<2,dim,typename VectorType::value_type> > > > &hessians,
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
        const VectorType *src = data.read_ptr<VectorType>(this->value_index(i));
        VectorSlice<std::vector<std::vector<typename VectorType::value_type> > > dst(values[i], component, n_comp);
        fe.get_function_values(*src, make_slice(index, start, size), dst, true);
      }

    for (unsigned int i=0; i<this->n_gradients(); ++i)
      {
        const VectorType *src = data.read_ptr<VectorType>(this->gradient_index(i));
        VectorSlice<std::vector<std::vector<Tensor<1,dim,typename VectorType::value_type> > > > dst(gradients[i], component, n_comp);
        fe.get_function_gradients(*src, make_slice(index, start, size), dst, true);
      }

    for (unsigned int i=0; i<this->n_hessians(); ++i)
      {
        const VectorType *src = data.read_ptr<VectorType>(this->hessian_index(i));
        VectorSlice<std::vector<std::vector<Tensor<2,dim,typename VectorType::value_type> > > > dst(hessians[i], component, n_comp);
        fe.get_function_hessians(*src, make_slice(index, start, size), dst, true);
      }
  }


  template <typename VectorType, int dim, int spacedim>
  std::size_t
  VectorData<VectorType, dim, spacedim>::memory_consumption () const
  {
    std::size_t mem = VectorSelector::memory_consumption();
    mem += sizeof (this->data);
    return mem;
  }

//----------------------------------------------------------------------//

  template <typename VectorType, int dim, int spacedim>
  MGVectorData<VectorType, dim, spacedim>::MGVectorData()
  {}


  template <typename VectorType, int dim, int spacedim>
  MGVectorData<VectorType, dim, spacedim>::MGVectorData(const VectorSelector &s)
    :
    VectorData<VectorType, dim, spacedim>(s)
  {}


  template <typename VectorType, int dim, int spacedim>
  void
  MGVectorData<VectorType, dim, spacedim>::initialize(const AnyData &d)
  {
    this->data = d;
    VectorSelector::initialize(d);
  }


  template <typename VectorType, int dim, int spacedim>
  void
  MGVectorData<VectorType, dim, spacedim>::initialize(const MGLevelObject<VectorType> *v, const std::string &name)
  {
    SmartPointer<const MGLevelObject<VectorType>, MGVectorData<VectorType, dim, spacedim> >
    p = v;
    this->data.add(p, name);
    VectorSelector::initialize(this->data);
  }


  template <typename VectorType, int dim, int spacedim>
  void
  VectorData<VectorType, dim, spacedim>::mg_fill
  (std::vector<std::vector<std::vector<typename VectorType::value_type> > >                &values,
   std::vector<std::vector<std::vector<Tensor<1,dim,typename VectorType::value_type> > > > &gradients,
   std::vector<std::vector<std::vector<Tensor<2,dim,typename VectorType::value_type> > > > &hessians,
   const FEValuesBase<dim,spacedim>           &fe,
   const unsigned int                         level,
   const std::vector<types::global_dof_index> &index,
   const unsigned int                         component,
   const unsigned int                         n_comp,
   const unsigned int                         start,
   const unsigned int                         size)
  const
  {
    AssertDimension(values.size(), this->n_values());
    AssertDimension(gradients.size(), this->n_gradients());
    AssertDimension(hessians.size(), this->n_hessians());

    const AnyData &data = this->data;
    for (unsigned int i=0; i<this->n_values(); ++i)
      {
        const MGLevelObject<VectorType> *src = data.read_ptr<MGLevelObject<VectorType> >(this->value_index(i));
        VectorSlice<std::vector<std::vector<typename VectorType::value_type> > > dst(values[i], component, n_comp);
        fe.get_function_values((*src)[level], make_slice(index, start, size), dst, true);
      }

    for (unsigned int i=0; i<this->n_gradients(); ++i)
      {
        const MGLevelObject<VectorType> *src = data.read_ptr<MGLevelObject<VectorType> >(this->value_index(i));
        VectorSlice<std::vector<std::vector<Tensor<1,dim,typename VectorType::value_type> > > > dst(gradients[i], component, n_comp);
        fe.get_function_gradients((*src)[level], make_slice(index, start, size), dst, true);
      }

    for (unsigned int i=0; i<this->n_hessians(); ++i)
      {
        const MGLevelObject<VectorType> *src = data.read_ptr<MGLevelObject<VectorType> >(this->value_index(i));
        VectorSlice<std::vector<std::vector<Tensor<2,dim,typename VectorType::value_type> > > > dst(hessians[i], component, n_comp);
        fe.get_function_hessians((*src)[level], make_slice(index, start, size), dst, true);
      }
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
