//----------------------------  hp_fe_values.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2006 by the deal.II authors
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
  
  namespace hp
  {
    template <int dim, class FEValues>
    FEValuesMap<dim,FEValues>::~FEValuesMap () 
    {}
  
  
    template <int dim, class FEValues>
    FEValues &
    FEValuesMap<dim,FEValues>::select_fe_values (const FiniteElement<dim> &fe,
                                                 const unsigned int active_fe_index)
    {
      std::pair<SmartPointer<const FiniteElement<dim> >, unsigned int> fe_pair=
	std::make_pair(&fe, active_fe_index);
                                       // check if the finite element
                                       // does not exist as a key in the
                                       // map
      if (fe_to_fe_values_map.find (fe_pair) == fe_to_fe_values_map.end())
                                         // a-ha! doesn't yet, so let's
                                         // make it up
        fe_to_fe_values_map[fe_pair]
          = boost::shared_ptr<FEValues> (create_fe_values (fe, active_fe_index));


                                       // now there definitely is one!
      present_fe_values = fe_to_fe_values_map[fe_pair];

      return *present_fe_values;
    }


// -------------------------- FEValuesBase -------------------------

    template <int dim, int q_dim>
    const MappingQ1<dim>
    FEValuesBase<dim,q_dim>::default_mapping;

  

    template <int dim, int q_dim>
    FEValuesBase<dim,q_dim>::FEValuesBase (
      const ::hp::MappingCollection<dim> &mapping_collection,
      const ::hp::QCollection<q_dim>     &q_collection,
      const UpdateFlags             update_flags)
                    :
                    mapping_collection (mapping_collection),
                    q_collection (q_collection),
                    update_flags (update_flags)
    {}


    template <int dim, int q_dim>
    FEValuesBase<dim,q_dim>::FEValuesBase (const ::hp::QCollection<q_dim> &q_collection,
                                           const UpdateFlags         update_flags)
                    :
                    mapping_collection (default_mapping),
                    q_collection (q_collection),
                    update_flags (update_flags)
    {}
    
  }
}



namespace hp
{
  
// -------------------------- FEValues -------------------------


  template <int dim>
  FEValues<dim>::FEValues (const hp::MappingCollection<dim> &mapping,
                           const hp::FECollection<dim>  &/*fe_collection*/,
                           const hp::QCollection<dim>       &q_collection,
                           const UpdateFlags             update_flags)
                  :
                  internal::hp::FEValuesBase<dim,dim> (mapping,
                                                       q_collection,
                                                       update_flags)
  {}


  template <int dim>
  FEValues<dim>::FEValues (const hp::FECollection<dim> &/*fe_collection*/,
                           const hp::QCollection<dim>      &q_collection,
                           const UpdateFlags            update_flags)
                  :
                  internal::hp::FEValuesBase<dim,dim> (q_collection,
                                                       update_flags)
  {}


  template <int dim>
  void
  FEValues<dim>::reinit (const typename hp::DoFHandler<dim>::cell_iterator &cell)
  {
    this->present_fe_index = cell->active_fe_index ();
    this->select_fe_values (cell->get_fe(), this->present_fe_index).reinit (cell);
  }



  template <int dim>
  ::FEValues<dim> *
  FEValues<dim>::create_fe_values (const FiniteElement<dim> &fe,
                                   const unsigned int active_fe_index) const
  {
    return new ::FEValues<dim> (
      this->mapping_collection[active_fe_index], fe,
      this->q_collection[active_fe_index], this->update_flags);
  }


// -------------------------- FEFaceValues -------------------------


  template <int dim>
  FEFaceValues<dim>::FEFaceValues (const hp::MappingCollection<dim> &mapping,
                                   const hp::FECollection<dim>  &/*fe_collection*/,
                                   const hp::QCollection<dim-1> &q_collection,
                                   const UpdateFlags         update_flags)
                  :
                  internal::hp::FEValuesBase<dim,dim-1> (mapping,
                                                         q_collection,
                                                         update_flags)
  {}


  template <int dim>
  FEFaceValues<dim>::FEFaceValues (const hp::FECollection<dim>  &/*fe_collection*/,
                                   const hp::QCollection<dim-1> &q_collection,
                                   const UpdateFlags         update_flags)
                  :
                  internal::hp::FEValuesBase<dim,dim-1> (q_collection,
                                                         update_flags)
  {}


  template <int dim>
  void
  FEFaceValues<dim>::reinit (const typename hp::DoFHandler<dim>::cell_iterator &cell,
                             const unsigned int face_no)
  {
    this->present_fe_index = cell->active_fe_index ();
    this->select_fe_values (cell->get_fe(), this->present_fe_index).reinit (cell, face_no);
  }


  template <int dim>
  void
  FEFaceValues<dim>::reinit (const typename hp::DoFHandler<dim>::cell_iterator &cell,
                             const unsigned int face_no,
                             const unsigned int active_fe_index)
  {
    this->present_fe_index = active_fe_index;
    this->select_fe_values (cell->get_fe(), active_fe_index).reinit (cell, face_no);
  }


  template <int dim>
  ::FEFaceValues<dim> *
  FEFaceValues<dim>::create_fe_values (const FiniteElement<dim> &fe,
                                       const unsigned int active_fe_index) const
  {
    return new ::FEFaceValues<dim> (
      this->mapping_collection[active_fe_index], fe,
      this->q_collection[active_fe_index], this->update_flags);
  }


// -------------------------- FESubfaceValues -------------------------


  template <int dim>
  FESubfaceValues<dim>::FESubfaceValues (const hp::MappingCollection<dim> &mapping,
                                         const hp::FECollection<dim>  &/*fe_collection*/,
                                         const hp::QCollection<dim-1> &q_collection,
                                         const UpdateFlags         update_flags)
                  :
                  internal::hp::FEValuesBase<dim,dim-1> (mapping,
                                                         q_collection,
                                                         update_flags)
  {}


  template <int dim>
  FESubfaceValues<dim>::FESubfaceValues (const hp::FECollection<dim>  &/*fe_collection*/,
                                         const hp::QCollection<dim-1> &q_collection,
                                         const UpdateFlags         update_flags)
                  :
                  internal::hp::FEValuesBase<dim,dim-1> (q_collection,
                                                         update_flags)
  {}


  template <int dim>
  void
  FESubfaceValues<dim>::reinit (const typename hp::DoFHandler<dim>::cell_iterator &cell,
                                const unsigned int face_no,
                                const unsigned int subface_no)
  {
    this->present_fe_index = cell->active_fe_index ();
    this->select_fe_values (cell->get_fe(), this->present_fe_index).reinit (cell, face_no, subface_no);
  }


  template <int dim>
  void
  FESubfaceValues<dim>::reinit (const typename hp::DoFHandler<dim>::cell_iterator &cell,
                                const unsigned int face_no,
                                const unsigned int subface_no,
                                const unsigned int active_fe_index)
  {
    this->present_fe_index = active_fe_index;
    this->select_fe_values (cell->get_fe(), active_fe_index).reinit (cell, face_no, subface_no);
  }


  template <int dim>
  ::FESubfaceValues<dim> *
  FESubfaceValues<dim>::create_fe_values (const FiniteElement<dim> &fe,
                                          const unsigned int active_fe_index) const
  {
    return new ::FESubfaceValues<dim> (
      this->mapping_collection[active_fe_index], fe,
      this->q_collection[active_fe_index], this->update_flags);
  }
}


// explicit instantiations
namespace internal
{
  namespace hp
  {
    template class FEValuesBase<deal_II_dimension,deal_II_dimension>;
#if deal_II_dimension >= 2
    template class FEValuesBase<deal_II_dimension,deal_II_dimension-1>;
#endif
  }
}

namespace hp
{
  template class FEValues<deal_II_dimension>;
#if deal_II_dimension >= 2
  template class FEFaceValues<deal_II_dimension>;
  template class FESubfaceValues<deal_II_dimension>;
#endif
}

// Putting the following explicit instantiations into the brackets 
// of the appropriate namespace somehow causes problems with the 
// Apple gcc3.3. Therefore these are separated.
template class internal::hp::FEValuesMap<deal_II_dimension,FEValues<deal_II_dimension> >;
template class internal::hp::FEValuesMap<deal_II_dimension,FEFaceValues<deal_II_dimension> >;
template class internal::hp::FEValuesMap<deal_II_dimension,FESubfaceValues<deal_II_dimension> >;
