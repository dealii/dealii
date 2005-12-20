//----------------------------  hp_fe_values.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  hp_fe_values.cc  ---------------------------


#include <fe/hp_fe_values.h>
#include <fe/mapping_q1.h>


// -------------------------- FEValuesMap -------------------------

namespace internal 
{
  template <int dim, class FEValues>
  FEValuesMap<dim,FEValues>::~FEValuesMap () 
  {}
  
  
  template <int dim, class FEValues>
  FEValues &
  FEValuesMap<dim,FEValues>::select_fe_values (const FiniteElement<dim> &fe)
  {
                                     // check if the finite element
                                     // does not exist as a key in the
                                     // map
    if (fe_to_fe_values_map.find (&fe) ==
        fe_to_fe_values_map.end())
                                       // a-ha! doesn't yet, so let's
                                       // make it up
      fe_to_fe_values_map[&fe]
        = boost::shared_ptr<FEValues> (create_fe_values (fe));


                                     // now there definitely is one!
    present_fe_values = fe_to_fe_values_map[&fe];

    return *present_fe_values;
  }


// -------------------------- hpFEValuesBase -------------------------

  template <int dim, int q_dim>
  const MappingQ1<dim>
  hpFEValuesBase<dim,q_dim>::default_mapping;

  

  template <int dim, int q_dim>
  hpFEValuesBase<dim,q_dim>::hpFEValuesBase (
      const MappingCollection<dim> &mapping_collection,
      const QCollection<q_dim>     &qcollection,
      const UpdateFlags             update_flags)
      :
                  mapping_collection (mapping_collection),
                  qcollection (qcollection),
                  update_flags (update_flags)
  {}


  template <int dim, int q_dim>
  hpFEValuesBase<dim,q_dim>::hpFEValuesBase (const QCollection<q_dim> &qcollection,
                                             const UpdateFlags         update_flags)
                  :
                  mapping_collection (default_mapping),
                  qcollection (qcollection),
                  update_flags (update_flags)
  {}

}


// -------------------------- hpFEValues -------------------------


template <int dim>
hpFEValues<dim>::hpFEValues (const MappingCollection<dim> &mapping,
                             const FECollection<dim>      &/*fe_collection*/,
                             const QCollection<dim>       &qcollection,
                             const UpdateFlags             update_flags)
                :
                internal::hpFEValuesBase<dim,dim> (mapping,
                                                   qcollection,
                                                   update_flags)
{}


template <int dim>
hpFEValues<dim>::hpFEValues (const FECollection<dim> &/*fe_collection*/,
                             const QCollection<dim>  &qcollection,
                             const UpdateFlags        update_flags)
                :
                internal::hpFEValuesBase<dim,dim> (qcollection,
                                                   update_flags)
{}


template <int dim>
void
hpFEValues<dim>::reinit (const typename hpDoFHandler<dim>::cell_iterator &cell)
{
  this->present_fe_index = cell->active_fe_index ();
  this->select_fe_values (cell->get_fe()).reinit (cell);
}



template <int dim>
FEValues<dim> *
hpFEValues<dim>::create_fe_values (const FiniteElement<dim> &fe) const
{
  return new FEValues<dim> (
      this->mapping_collection.get_mapping (this->present_fe_index), fe,
      this->qcollection.get_quadrature (this->present_fe_index),
      this->update_flags);
}


// -------------------------- hpFEFaceValues -------------------------


template <int dim>
hpFEFaceValues<dim>::hpFEFaceValues (const MappingCollection<dim> &mapping,
                                     const FECollection<dim>  &/*fe_collection*/,
                                     const QCollection<dim-1> &qcollection,
                                     const UpdateFlags         update_flags)
                :
                internal::hpFEValuesBase<dim,dim-1> (mapping,
                                                     qcollection,
                                                     update_flags)
{}


template <int dim>
hpFEFaceValues<dim>::hpFEFaceValues (const FECollection<dim>  &/*fe_collection*/,
                                     const QCollection<dim-1> &qcollection,
                                     const UpdateFlags         update_flags)
                :
                internal::hpFEValuesBase<dim,dim-1> (qcollection,
                                                     update_flags)
{}


template <int dim>
void
hpFEFaceValues<dim>::reinit (const typename hpDoFHandler<dim>::cell_iterator &cell,
                             const unsigned int face_no)
{
  this->present_fe_index = cell->active_fe_index ();
  this->select_fe_values (cell->get_fe()).reinit (cell, face_no);
}



template <int dim>
FEFaceValues<dim> *
hpFEFaceValues<dim>::create_fe_values (const FiniteElement<dim> &fe) const
{
  return new FEFaceValues<dim> (
      this->mapping_collection.get_mapping (this->present_fe_index), fe,
      this->qcollection.get_quadrature (this->present_fe_index),
      this->update_flags);
}


// -------------------------- hpFESubfaceValues -------------------------


template <int dim>
hpFESubfaceValues<dim>::hpFESubfaceValues (const MappingCollection<dim> &mapping,
                                           const FECollection<dim>  &/*fe_collection*/,
                                           const QCollection<dim-1> &qcollection,
                                           const UpdateFlags         update_flags)
                :
                internal::hpFEValuesBase<dim,dim-1> (mapping,
                                                     qcollection,
                                                     update_flags)
{}


template <int dim>
hpFESubfaceValues<dim>::hpFESubfaceValues (const FECollection<dim>  &/*fe_collection*/,
                                           const QCollection<dim-1> &qcollection,
                                           const UpdateFlags         update_flags)
                :
                internal::hpFEValuesBase<dim,dim-1> (qcollection,
                                                     update_flags)
{}


template <int dim>
void
hpFESubfaceValues<dim>::reinit (const typename hpDoFHandler<dim>::cell_iterator &cell,
                                const unsigned int face_no,
                                const unsigned int subface_no)
{
  this->present_fe_index = cell->active_fe_index ();
  this->select_fe_values (cell->get_fe()).reinit (cell, face_no, subface_no);
}



template <int dim>
FESubfaceValues<dim> *
hpFESubfaceValues<dim>::create_fe_values (const FiniteElement<dim> &fe) const
{
  return new FESubfaceValues<dim> (
      this->mapping_collection.get_mapping (this->present_fe_index), fe,
      this->qcollection.get_quadrature (this->present_fe_index),
      this->update_flags);
}




// explicit instantiations
namespace internal
{
  template class FEValuesMap<deal_II_dimension,FEValues<deal_II_dimension> >;
  template class FEValuesMap<deal_II_dimension,FEFaceValues<deal_II_dimension> >;
  template class FEValuesMap<deal_II_dimension,FESubfaceValues<deal_II_dimension> >;

  template class hpFEValuesBase<deal_II_dimension,deal_II_dimension>;
#if deal_II_dimension >= 2
  template class hpFEValuesBase<deal_II_dimension,deal_II_dimension-1>;
#endif
}

template class hpFEValues<deal_II_dimension>;
#if deal_II_dimension >= 2
template class hpFEFaceValues<deal_II_dimension>;
template class hpFESubfaceValues<deal_II_dimension>;
#endif

