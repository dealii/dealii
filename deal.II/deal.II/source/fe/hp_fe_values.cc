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


namespace internal
{
  
  namespace hp
  {
// -------------------------- FEValuesBase -------------------------

    template <int dim, int q_dim, class FEValues>
    FEValuesBase<dim,q_dim,FEValues>::
    FEValuesBase (const ::hp::MappingCollection<dim> &mapping_collection,
                  const ::hp::FECollection<dim>      &fe_collection,
                  const ::hp::QCollection<q_dim>     &q_collection,
                  const UpdateFlags                   update_flags)
                    :
                    fe_collection (&fe_collection),
                    mapping_collection (&mapping_collection),
                    q_collection (q_collection),
                    fe_values_table (fe_collection.size(),
                                     mapping_collection.size(),
                                     q_collection.size()),
                    present_fe_values_index (deal_II_numbers::invalid_unsigned_int,
                                             deal_II_numbers::invalid_unsigned_int,
                                             deal_II_numbers::invalid_unsigned_int),
                    update_flags (update_flags)
    {}


    template <int dim, int q_dim, class FEValues>
    FEValuesBase<dim,q_dim,FEValues>::
    FEValuesBase (const ::hp::FECollection<dim>      &fe_collection,
                  const ::hp::QCollection<q_dim> &q_collection,
                  const UpdateFlags         update_flags)
                    :
                    fe_collection (&fe_collection),
                    mapping_collection (&::hp::StaticMappingQ1<dim>::mapping_collection),
                    q_collection (q_collection),
                    fe_values_table (fe_collection.size(),
                                     1,
                                     q_collection.size()),
                    present_fe_values_index (deal_II_numbers::invalid_unsigned_int,
                                             deal_II_numbers::invalid_unsigned_int,
                                             deal_II_numbers::invalid_unsigned_int),
                    update_flags (update_flags)
    {}



    template <int dim, int q_dim, class FEValues>
    FEValues &
    FEValuesBase<dim,q_dim,FEValues>::
    select_fe_values (const unsigned int fe_index,
                      const unsigned int mapping_index,
                      const unsigned int q_index)
    {
      Assert (fe_index < fe_collection->size(),
              ExcIndexRange (fe_index, 0, fe_collection->size()));
      Assert (mapping_index < mapping_collection->size(),
              ExcIndexRange (mapping_index, 0, mapping_collection->size()));
      Assert (q_index < q_collection.size(),
              ExcIndexRange (q_index, 0, q_collection.size()));


                                       // set the triple of indices
                                       // that we want to work with
      present_fe_values_index = TableIndices<3> (fe_index,
                                                 mapping_index,
                                                 q_index);
                                       
                                       // first check whether we
                                       // already have an object for
                                       // this particular combination
                                       // of indices
      if (fe_values_table(present_fe_values_index) == 0)
        fe_values_table(present_fe_values_index)
          =
          boost::shared_ptr<FEValues>
          (new FEValues ((*mapping_collection)[mapping_index],
                         (*fe_collection)[fe_index],
                         q_collection[q_index],
                         update_flags));

                                       // now there definitely is one!
      return *fe_values_table(present_fe_values_index);
    }    
  }
}



namespace hp
{
  
// -------------------------- FEValues -------------------------


  template <int dim>
  FEValues<dim>::FEValues (const hp::MappingCollection<dim> &mapping,
                           const hp::FECollection<dim>      &fe_collection,
                           const hp::QCollection<dim>       &q_collection,
                           const UpdateFlags                 update_flags)
                  :
                  internal::hp::FEValuesBase<dim,dim,::FEValues<dim> > (mapping,
                                                                  fe_collection,
                                                                  q_collection,
                                                                  update_flags)
  {}


  template <int dim>
  FEValues<dim>::FEValues (const hp::FECollection<dim> &fe_collection,
                           const hp::QCollection<dim>      &q_collection,
                           const UpdateFlags            update_flags)
                  :
                  internal::hp::FEValuesBase<dim,dim,::FEValues<dim> > (fe_collection,
                                                       q_collection,
                                                       update_flags)
  {}


  template <int dim>
  void
  FEValues<dim>::reinit (const typename hp::DoFHandler<dim>::cell_iterator &cell,
                         const unsigned int q_index,
                         const unsigned int mapping_index,
                         const unsigned int fe_index)
  {
                                     // determine which indices we
                                     // should actually use
    unsigned int real_q_index       = q_index,
                 real_mapping_index = mapping_index,
                 real_fe_index      = fe_index;

    if (real_q_index == deal_II_numbers::invalid_unsigned_int)
      if (this->q_collection.size() > 1)
        real_q_index = cell->active_fe_index();
      else
        real_q_index = 0;
    
    if (real_mapping_index == deal_II_numbers::invalid_unsigned_int)
      if (this->mapping_collection->size() > 1)
        real_mapping_index = cell->active_fe_index();
      else
        real_mapping_index = 0;

    if (real_fe_index == deal_II_numbers::invalid_unsigned_int)
      real_fe_index = cell->active_fe_index();

                                     // some checks
    Assert (real_q_index < this->q_collection.size(),
            ExcIndexRange (real_q_index, 0, this->q_collection.size()));
    Assert (real_mapping_index < this->mapping_collection->size(),
            ExcIndexRange (real_mapping_index, 0, this->mapping_collection->size()));
    Assert (real_fe_index < this->fe_collection->size(),
            ExcIndexRange (real_fe_index, 0, this->fe_collection->size()));
    
                                     // now finally actually get the
                                     // corresponding object and
                                     // initialize it
    this->select_fe_values (real_fe_index,
                            real_mapping_index,
                            real_q_index).reinit (cell);
  }


// -------------------------- FEFaceValues -------------------------


  template <int dim>
  FEFaceValues<dim>::FEFaceValues (const hp::MappingCollection<dim> &mapping,
                                   const hp::FECollection<dim>  &fe_collection,
                                   const hp::QCollection<dim-1> &q_collection,
                                   const UpdateFlags         update_flags)
                  :
                  internal::hp::FEValuesBase<dim,dim-1,::FEFaceValues<dim> > (mapping,
                                                         fe_collection,
                                                         q_collection,
                                                         update_flags)
  {}


  template <int dim>
  FEFaceValues<dim>::FEFaceValues (const hp::FECollection<dim>  &fe_collection,
                                   const hp::QCollection<dim-1> &q_collection,
                                   const UpdateFlags         update_flags)
                  :
                  internal::hp::FEValuesBase<dim,dim-1,::FEFaceValues<dim> > (fe_collection,
                                                         q_collection,
                                                         update_flags)
  {}


  template <int dim>
  void
  FEFaceValues<dim>::reinit (const typename hp::DoFHandler<dim>::cell_iterator &cell,
                             const unsigned int face_no,
                             const unsigned int q_index,
                             const unsigned int mapping_index,
                             const unsigned int fe_index)
  {
                                     // determine which indices we
                                     // should actually use
    unsigned int real_q_index       = q_index,
                 real_mapping_index = mapping_index,
                 real_fe_index      = fe_index;

    if (real_q_index == deal_II_numbers::invalid_unsigned_int)
      if (this->q_collection.size() > 1)
        real_q_index = cell->active_fe_index();
      else
        real_q_index = 0;
    
    if (real_mapping_index == deal_II_numbers::invalid_unsigned_int)
      if (this->mapping_collection->size() > 1)
        real_mapping_index = cell->active_fe_index();
      else
        real_mapping_index = 0;

    if (real_fe_index == deal_II_numbers::invalid_unsigned_int)
      real_fe_index = cell->active_fe_index();

                                     // some checks
    Assert (real_q_index < this->q_collection.size(),
            ExcIndexRange (real_q_index, 0, this->q_collection.size()));
    Assert (real_mapping_index < this->mapping_collection->size(),
            ExcIndexRange (real_mapping_index, 0, this->mapping_collection->size()));
    Assert (real_fe_index < this->fe_collection->size(),
            ExcIndexRange (real_fe_index, 0, this->fe_collection->size()));
    
                                     // now finally actually get the
                                     // corresponding object and
                                     // initialize it
    this->select_fe_values (real_fe_index,
                            real_mapping_index,
                            real_q_index).reinit (cell, face_no);
  }


// -------------------------- FESubfaceValues -------------------------


  template <int dim>
  FESubfaceValues<dim>::FESubfaceValues (const hp::MappingCollection<dim> &mapping,
                                         const hp::FECollection<dim>  &fe_collection,
                                         const hp::QCollection<dim-1> &q_collection,
                                         const UpdateFlags         update_flags)
                  :
                  internal::hp::FEValuesBase<dim,dim-1,::FESubfaceValues<dim> > (mapping,
                                                         fe_collection,
                                                         q_collection,
                                                         update_flags)
  {}


  template <int dim>
  FESubfaceValues<dim>::FESubfaceValues (const hp::FECollection<dim>  &fe_collection,
                                         const hp::QCollection<dim-1> &q_collection,
                                         const UpdateFlags         update_flags)
                  :
                  internal::hp::FEValuesBase<dim,dim-1,::FESubfaceValues<dim> > (fe_collection,
                                                         q_collection,
                                                         update_flags)
  {}


  template <int dim>
  void
  FESubfaceValues<dim>::reinit (const typename hp::DoFHandler<dim>::cell_iterator &cell,
                                const unsigned int face_no,
                                const unsigned int subface_no,
                                const unsigned int q_index,
                                const unsigned int mapping_index,
                                const unsigned int fe_index)
  {
                                     // determine which indices we
                                     // should actually use
    unsigned int real_q_index       = q_index,
                 real_mapping_index = mapping_index,
                 real_fe_index      = fe_index;

    if (real_q_index == deal_II_numbers::invalid_unsigned_int)
      if (this->q_collection.size() > 1)
        real_q_index = cell->active_fe_index();
      else
        real_q_index = 0;
    
    if (real_mapping_index == deal_II_numbers::invalid_unsigned_int)
      if (this->mapping_collection->size() > 1)
        real_mapping_index = cell->active_fe_index();
      else
        real_mapping_index = 0;

    if (real_fe_index == deal_II_numbers::invalid_unsigned_int)
      real_fe_index = cell->active_fe_index();

                                     // some checks
    Assert (real_q_index < this->q_collection.size(),
            ExcIndexRange (real_q_index, 0, this->q_collection.size()));
    Assert (real_mapping_index < this->mapping_collection->size(),
            ExcIndexRange (real_mapping_index, 0, this->mapping_collection->size()));
    Assert (real_fe_index < this->fe_collection->size(),
            ExcIndexRange (real_fe_index, 0, this->fe_collection->size()));
    
                                     // now finally actually get the
                                     // corresponding object and
                                     // initialize it
    this->select_fe_values (real_fe_index,
                            real_mapping_index,
                            real_q_index).reinit (cell, face_no, subface_no);
  }
}


// explicit instantiations
namespace internal
{
  namespace hp
  {
    template class FEValuesBase<deal_II_dimension,deal_II_dimension,
                                ::FEValues<deal_II_dimension> >;
#if deal_II_dimension >= 2
    template class FEValuesBase<deal_II_dimension,deal_II_dimension-1,
                                ::FEFaceValues<deal_II_dimension> >;
    template class FEValuesBase<deal_II_dimension,deal_II_dimension-1,
                                ::FESubfaceValues<deal_II_dimension> >;
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
