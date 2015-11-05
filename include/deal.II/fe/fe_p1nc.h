// ---------------------------------------------------------------------
//
// Copyright (C) 2015 by the deal.II authors
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

#ifndef dealii__fe_p1nc_h
#define dealii__fe_p1nc_h

#include <deal.II/base/config.h>
#include <deal.II/fe/fe.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/qprojector.h>


DEAL_II_NAMESPACE_OPEN


/*!@addtogroup fe */
/*@{*/


// parametric version of P1 Nonconforming FE
class FE_P1NCParametric : public FiniteElement<2,2>
{
private:


  static
  std::vector<ComponentMask>
  get_nonzero_component();


  static
  std::vector<unsigned int>
  get_dpo_vector ();


  class InternalData : public FiniteElement<2,2>::InternalDataBase
  {
  public:

    std::vector<std::vector<double> > shape_values;
    std::vector<std::vector<Tensor<1,2> > > shape_gradients;
  };


  virtual FiniteElement<2,2>::InternalDataBase *
  get_data (const UpdateFlags update_flags,
            const Mapping<2,2> &mapping,
            const Quadrature<2> &quadrature) const ;




  virtual
  void
  fill_fe_values (const Mapping<2,2>                                &mapping,
                  const Triangulation<2,2>::cell_iterator           &cell,
                  const Quadrature<2>                               &quadrature,
                  const Mapping<2,2>::InternalDataBase              &mapping_internal,
                  const FiniteElement<2,2>::InternalDataBase        &fe_internal,
                  const internal::FEValues::MappingRelatedData<2,2> &mapping_data,
                  internal::FEValues::FiniteElementRelatedData<2,2> &output_data,
                  const CellSimilarity::Similarity                   cell_similarity) const;




  virtual
  void
  fill_fe_face_values (const Mapping<2,2>                                &mapping,
                       const Triangulation<2,2>::cell_iterator           &cell,
                       const unsigned int                                 face_no,
                       const Quadrature<1>                               &quadrature,
                       const Mapping<2,2>::InternalDataBase              &mapping_internal,
                       const FiniteElement<2,2>::InternalDataBase        &fe_internal,
                       const internal::FEValues::MappingRelatedData<2,2> &mapping_data,
                       internal::FEValues::FiniteElementRelatedData<2,2> &output_data) const;




  virtual
  void
  fill_fe_subface_values (const Mapping<2,2>                                &mapping,
                          const Triangulation<2,2>::cell_iterator           &cell,
                          const unsigned int                                 face_no,
                          const unsigned int                                 sub_no,
                          const Quadrature<1>                               &quadrature,
                          const Mapping<2,2>::InternalDataBase              &mapping_internal,
                          const FiniteElement<2,2>::InternalDataBase        &fe_internal,
                          const internal::FEValues::MappingRelatedData<2,2> &mapping_data,
                          internal::FEValues::FiniteElementRelatedData<2,2> &output_data) const;



public:
  FE_P1NCParametric() ;



  virtual std::string get_name () const ;



  virtual UpdateFlags     update_once (const UpdateFlags flags) const ;

  virtual UpdateFlags     update_each (const UpdateFlags flags) const ;

  const std::vector<Point<1> > &get_unit_face_support_points () const ;

  Point<1> unit_face_support_point (const unsigned index) const ;

  bool has_face_support_points () const ;

  virtual FiniteElement<2,2> *clone () const ;

  virtual ~FE_P1NCParametric ();





  virtual double shape_value (const unsigned int i,
                              const Point<2>    &p) const ;

  virtual Tensor<1,2> shape_grad (const unsigned int i,
                                  const Point<2>    &p) const ;



protected:

  void initialize_constraints () ;

};





/** @}*/

DEAL_II_NAMESPACE_CLOSE

#endif
